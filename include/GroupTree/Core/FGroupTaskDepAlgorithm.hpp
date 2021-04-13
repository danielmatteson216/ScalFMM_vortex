
// Keep in private GIT
#ifndef FGROUPTASKDEPALGORITHM_HPP
#define FGROUPTASKDEPALGORITHM_HPP


#include "Utils/FGlobal.hpp"
#include "Core/FCoreCommon.hpp"
#include "Utils/FQuickSort.hpp"
#include "Containers/FTreeCoordinate.hpp"
#include "Utils/FLog.hpp"
#include "Utils/FTic.hpp"
#include "Utils/FAlgorithmTimers.hpp"

#include "Utils/FTaskTimer.hpp"

#include "FOutOfBlockInteraction.hpp"

#include <vector>

#include <omp.h>



#undef commute_if_supported
#ifdef OPENMP_SUPPORT_COMMUTE
#define commute_if_supported commute
#else
#define commute_if_supported inout
#endif

#undef priority_if_supported
#ifdef OPENMP_SUPPORT_PRIORITY
#include "GroupTree/StarPUUtils/FOmpPriorities.hpp"
#define priority_if_supported(x) priority(x)
#else
#define priority_if_supported(x)
#endif

#undef taskname_if_supported
#ifdef OPENMP_SUPPORT_TASK_NAME
#define taskname_if_supported(n) taskname(n)
#else
#define taskname_if_supported(n)
#endif


template <class OctreeClass, class CellContainerClass,
          class SymboleCellClass, class PoleCellClass, class LocalCellClass, class KernelClass, class ParticleGroupClass, class ParticleContainerClass>
class FGroupTaskDepAlgorithm : public FAbstractAlgorithm, public FAlgorithmTimers {
protected:
    template <class OtherBlockClass>
    struct BlockInteractions{
        OtherBlockClass* otherBlock;
        std::vector<OutOfBlockInteraction> interactions;
    };

    std::vector< std::vector< std::vector<BlockInteractions<CellContainerClass>>>> externalInteractionsAllLevel;
    std::vector< std::vector<BlockInteractions<ParticleGroupClass>>> externalInteractionsLeafLevel;

    int MaxThreads;         //< The number of threads
    OctreeClass*const tree;       //< The Tree
    KernelClass** kernels;        //< The kernels
    const bool noCommuteAtLastLevel;

#ifdef SCALFMM_TIME_OMPTASKS
    FTaskTimer taskTimeRecorder;
#endif

#ifdef OPENMP_SUPPORT_PRIORITY
    FOmpPriorities priorities;
#endif

public:
    using multipole_t = PoleCellClass;
    using local_expansion_t = LocalCellClass;
    using symbolic_data_t = SymboleCellClass;

    FGroupTaskDepAlgorithm(OctreeClass*const inTree, const KernelClass* const inKernels, const int inMaxThreads = -1)
        : tree(inTree), kernels(nullptr),
          noCommuteAtLastLevel(getenv("SCALFMM_NO_COMMUTE_LAST_L2L") != NULL && getenv("SCALFMM_NO_COMMUTE_LAST_L2L")[0] != '1'?false:true)
#ifdef SCALFMM_TIME_OMPTASKS
            , taskTimeRecorder(MaxThreads)
#endif
#ifdef OPENMP_SUPPORT_PRIORITY
            , priorities(tree->getHeight())
#endif
    {
        FAssertLF(tree, "tree cannot be null");
        FAssertLF(inKernels, "kernels cannot be null");

        FAbstractAlgorithm::setNbLevelsInTree(tree->getHeight());

        MaxThreads = 1;
        #pragma omp parallel
        #pragma omp master
            MaxThreads = omp_get_num_threads();

        this->kernels = new KernelClass*[MaxThreads];
        #pragma omp parallel num_threads(MaxThreads)
        {
          #pragma omp critical (FGroupTaskDepAlgorithm_InitKernels)
            this->kernels[omp_get_thread_num()] = new KernelClass(*inKernels);
        }

        rebuildInteractions();

        FLOG(FLog::Controller << "FGroupTaskDepAlgorithm (Max Thread " << MaxThreads << ")\n");

#ifdef SCALFMM_TIME_OMPTASKS
        #pragma omp parallel num_threads(MaxThreads)
        {
            taskTimeRecorder.init(omp_get_thread_num());
        }
#endif
        FLOG(FLog::Controller << "SCALFMM_NO_COMMUTE_LAST_L2L " << noCommuteAtLastLevel << "\n");
    }

    ~FGroupTaskDepAlgorithm(){
        for(int idxThread = 0 ; idxThread < MaxThreads ; ++idxThread){
            delete this->kernels[idxThread];
        }
        delete[] kernels;
    }

    void rebuildInteractions(){
        #pragma omp parallel num_threads(MaxThreads)
        {
            #pragma omp single nowait
            {
                // For now rebuild all external interaction
                buildExternalInteractionVecs();
            }
        }
    }

protected:
    /**
      * Runs the complete algorithm.
      */
    void executeCore(const unsigned operationsToProceed) override {
        FLOG( FLog::Controller << "\tStart FGroupTaskDepAlgorithm\n" );

        FTIME_TASKS(taskTimeRecorder.start());
        Timers[P2MTimer].tic();
        #pragma omp parallel num_threads(MaxThreads)
        {
            #pragma omp single nowait
            {
                FLOG( FTic timerSoumission; );

                if( operationsToProceed & FFmmP2P ) directPass();

                if(operationsToProceed & FFmmP2M) bottomPass();

                if(operationsToProceed & FFmmM2M) upwardPass();

                if(operationsToProceed & FFmmM2L) transferPass( FAbstractAlgorithm::upperWorkingLevel,FAbstractAlgorithm::lowerWorkingLevel-1);

                if(operationsToProceed & FFmmL2L) downardPass();

                if(operationsToProceed & FFmmM2L) transferPass(FAbstractAlgorithm::lowerWorkingLevel-1, FAbstractAlgorithm::lowerWorkingLevel);

                if( operationsToProceed & FFmmL2P ) mergePass();

                FLOG( FLog::Controller << "\t\t Submitting the tasks took " << timerSoumission.tacAndElapsed() << "s\n" );
                #pragma omp taskwait
            }
        }
        Timers[P2MTimer].tac();
        FTIME_TASKS(taskTimeRecorder.end());
        FTIME_TASKS(taskTimeRecorder.saveToDisk("/tmp/taskstime-FGroupTaskDepAlgorithm.txt"));
    }


    /**
     * This function is creating the interactions vector between blocks.
     * It fills externalInteractionsAllLevel and externalInteractionsLeafLevel.
     * Warning, the omp task for now are using the class attributes!
     *
     */
    void buildExternalInteractionVecs(){
            FLOG( FTic timer; FTic leafTimer; FTic cellTimer; );
            // Reset interactions
            externalInteractionsAllLevel.clear();
            externalInteractionsLeafLevel.clear();
            // One per level + leaf level
            externalInteractionsAllLevel.resize(tree->getHeight());

            // First leaf level
            {
                // We create one big vector per block
                externalInteractionsLeafLevel.resize(tree->getNbParticleGroup());

                for(int idxGroup = 0 ; idxGroup < tree->getNbParticleGroup() ; ++idxGroup){
                    // Create the vector
                    ParticleGroupClass* containers = tree->getParticleGroup(idxGroup);

                    std::vector<BlockInteractions<ParticleGroupClass>>* externalInteractions = &externalInteractionsLeafLevel[idxGroup];

                    #pragma omp task default(none) firstprivate(idxGroup, containers, externalInteractions)
                    { // Can be a task(inout:iterCells)
                        std::vector<OutOfBlockInteraction> outsideInteractions;
                        const MortonIndex blockStartIdx = containers->getStartingIndex();
                        const MortonIndex blockEndIdx   = containers->getEndingIndex();

                        for(int leafIdx = 0 ; leafIdx < containers->getNumberOfLeavesInBlock() ; ++leafIdx){
                            const MortonIndex mindex = containers->getLeafMortonIndex(leafIdx);
                            // ParticleContainerClass particles = containers->template getLeaf<ParticleContainerClass>(leafIdx);

                            MortonIndex interactionsIndexes[26];
                            int interactionsPosition[26];
                            FTreeCoordinate coord(mindex);
                            int counter = coord.getNeighborsIndexes(tree->getHeight(),interactionsIndexes,interactionsPosition);

                            for(int idxInter = 0 ; idxInter < counter ; ++idxInter){
                                if( blockStartIdx <= interactionsIndexes[idxInter] && interactionsIndexes[idxInter] < blockEndIdx ){
                                    // Inside block interaction, do nothing
                                }
                                else if(interactionsIndexes[idxInter] < mindex){
                                    OutOfBlockInteraction property;
                                    property.insideIndex = mindex;
                                    property.outIndex    = interactionsIndexes[idxInter];
                                    property.relativeOutPosition = interactionsPosition[idxInter];
                                    property.insideIdxInBlock = leafIdx;
                                    property.outsideIdxInBlock = -1;
                                    outsideInteractions.push_back(property);
                                }
                            }
                        }

                        // Sort to match external order
                        FQuickSort<OutOfBlockInteraction, int>::QsSequential(outsideInteractions.data(),int(outsideInteractions.size()));

                        int currentOutInteraction = 0;
                        for(int idxLeftGroup = 0 ; idxLeftGroup < idxGroup && currentOutInteraction < int(outsideInteractions.size()) ; ++idxLeftGroup){
                            ParticleGroupClass* leftContainers = tree->getParticleGroup(idxLeftGroup);
                            const MortonIndex blockStartIdxOther    = leftContainers->getStartingIndex();
                            const MortonIndex blockEndIdxOther      = leftContainers->getEndingIndex();

                            while(currentOutInteraction < int(outsideInteractions.size())
                                  && (outsideInteractions[currentOutInteraction].outIndex < blockStartIdxOther
                                      || leftContainers->getLeafIndex(outsideInteractions[currentOutInteraction].outIndex) == -1)
                                  && outsideInteractions[currentOutInteraction].outIndex < blockEndIdxOther){
                                currentOutInteraction += 1;
                            }

                            int lastOutInteraction = currentOutInteraction;
                            int copyExistingInteraction = currentOutInteraction;
                            while(lastOutInteraction < int(outsideInteractions.size()) && outsideInteractions[lastOutInteraction].outIndex < blockEndIdxOther){
                                const int leafPos = leftContainers->getLeafIndex(outsideInteractions[lastOutInteraction].outIndex);
                                if(leafPos != -1){
                                    if(copyExistingInteraction != lastOutInteraction){
                                        outsideInteractions[copyExistingInteraction] = outsideInteractions[lastOutInteraction];
                                    }
                                    outsideInteractions[copyExistingInteraction].outsideIdxInBlock = leafPos;
                                    copyExistingInteraction += 1;
                                }
                                lastOutInteraction += 1;
                            }

                            const int nbInteractionsBetweenBlocks = (copyExistingInteraction-currentOutInteraction);
                            if(nbInteractionsBetweenBlocks){
                                externalInteractions->emplace_back();
                                BlockInteractions<ParticleGroupClass>* interactions = &externalInteractions->back();
                                interactions->otherBlock = leftContainers;
                                interactions->interactions.resize(nbInteractionsBetweenBlocks);
                                std::copy(outsideInteractions.begin() + currentOutInteraction,
                                          outsideInteractions.begin() + copyExistingInteraction,
                                          interactions->interactions.begin());
                            }

                            currentOutInteraction = lastOutInteraction;
                        }
                    }
                }
            }
            FLOG( leafTimer.tac(); );
            FLOG( cellTimer.tic(); );
            {
                for(int idxLevel = tree->getHeight()-1 ; idxLevel >= 2 ; --idxLevel){
                    externalInteractionsAllLevel[idxLevel].resize(tree->getNbCellGroupAtLevel(idxLevel));

                    for(int idxGroup = 0 ; idxGroup < tree->getNbCellGroupAtLevel(idxLevel) ; ++idxGroup){
                        CellContainerClass* currentCells = tree->getCellGroup(idxLevel, idxGroup);

                        std::vector<BlockInteractions<CellContainerClass>>* externalInteractions = &externalInteractionsAllLevel[idxLevel][idxGroup];

                        #pragma omp task default(none) firstprivate(idxGroup, currentCells, idxLevel, externalInteractions)
                        {
                            std::vector<OutOfBlockInteraction> outsideInteractions;
                            const MortonIndex blockStartIdx = currentCells->getStartingIndex();
                            const MortonIndex blockEndIdx   = currentCells->getEndingIndex();

                            for(int cellIdx = 0 ; cellIdx < currentCells->getNumberOfCellsInBlock() ; ++cellIdx){
                                const MortonIndex mindex = currentCells->getCellMortonIndex(cellIdx);

                                MortonIndex interactionsIndexes[189];
                                int interactionsPosition[189];
                                const FTreeCoordinate coord(mindex);
                                int counter = coord.getInteractionNeighbors(idxLevel,interactionsIndexes,interactionsPosition);

                                for(int idxInter = 0 ; idxInter < counter ; ++idxInter){
                                    if( blockStartIdx <= interactionsIndexes[idxInter] && interactionsIndexes[idxInter] < blockEndIdx ){
                                        // Nothing to do
                                    }
                                    else if(interactionsIndexes[idxInter] < mindex){
                                        OutOfBlockInteraction property;
                                        property.insideIndex = mindex;
                                        property.outIndex    = interactionsIndexes[idxInter];
                                        property.relativeOutPosition = interactionsPosition[idxInter];
                                        property.insideIdxInBlock = cellIdx;
                                        property.outsideIdxInBlock = -1;
                                        outsideInteractions.push_back(property);
                                    }
                                }
                            }

                            // Manage outofblock interaction
                            FQuickSort<OutOfBlockInteraction, int>::QsSequential(outsideInteractions.data(),int(outsideInteractions.size()));

                            int currentOutInteraction = 0;
                            for(int idxLeftGroup = 0 ; idxLeftGroup < idxGroup && currentOutInteraction < int(outsideInteractions.size()) ; ++idxLeftGroup){
                                CellContainerClass* leftCells   = tree->getCellGroup(idxLevel, idxLeftGroup);
                                const MortonIndex blockStartIdxOther = leftCells->getStartingIndex();
                                const MortonIndex blockEndIdxOther   = leftCells->getEndingIndex();

                                while(currentOutInteraction < int(outsideInteractions.size())
                                      && (outsideInteractions[currentOutInteraction].outIndex < blockStartIdxOther
                                          || leftCells->getCellIndex(outsideInteractions[currentOutInteraction].outIndex) == -1)
                                      && outsideInteractions[currentOutInteraction].outIndex < blockEndIdxOther){
                                    currentOutInteraction += 1;
                                }

                                int lastOutInteraction = currentOutInteraction;
                                int copyExistingInteraction = currentOutInteraction;
                                while(lastOutInteraction < int(outsideInteractions.size()) && outsideInteractions[lastOutInteraction].outIndex < blockEndIdxOther){
                                    const int cellPos = leftCells->getCellIndex(outsideInteractions[lastOutInteraction].outIndex);
                                    if(cellPos != -1){
                                        if(copyExistingInteraction != lastOutInteraction){
                                            outsideInteractions[copyExistingInteraction] = outsideInteractions[lastOutInteraction];
                                        }
                                        outsideInteractions[copyExistingInteraction].outsideIdxInBlock = cellPos;
                                        copyExistingInteraction += 1;
                                    }
                                    lastOutInteraction += 1;
                                }

                                // Create interactions
                                const int nbInteractionsBetweenBlocks = (copyExistingInteraction-currentOutInteraction);
                                if(nbInteractionsBetweenBlocks){
                                    externalInteractions->emplace_back();
                                    BlockInteractions<CellContainerClass>* interactions = &externalInteractions->back();
                                    interactions->otherBlock = leftCells;
                                    interactions->interactions.resize(nbInteractionsBetweenBlocks);
                                    std::copy(outsideInteractions.begin() + currentOutInteraction,
                                              outsideInteractions.begin() + copyExistingInteraction,
                                              interactions->interactions.begin());
                                }

                                currentOutInteraction = lastOutInteraction;
                            }
                        }
                    }
                }
            }
            FLOG( cellTimer.tac(); );

            #pragma omp taskwait

            FLOG( FLog::Controller << "\t\t Prepare in " << timer.tacAndElapsed() << "s\n" );
            FLOG( FLog::Controller << "\t\t\t Prepare at leaf level in   " << leafTimer.elapsed() << "s\n" );
            FLOG( FLog::Controller << "\t\t\t Prepare at other levels in " << cellTimer.elapsed() << "s\n" );
        }


    void bottomPass(){
        FLOG( FTic timer; );

        for(int idxGroup = 0 ; idxGroup < tree->getNbParticleGroup() ; ++idxGroup){
            CellContainerClass* leafCells  = tree->getCellGroup(tree->getHeight()-1, idxGroup);
            PoleCellClass* cellPoles = leafCells->getRawMultipoleBuffer();

            ParticleGroupClass* containers = tree->getParticleGroup(idxGroup);

            #pragma omp task default(shared) firstprivate(leafCells, cellPoles, containers) depend(inout: cellPoles[0]) priority_if_supported(priorities.getInsertionPosP2M()) taskname_if_supported("P2M")
            {
                FTIME_TASKS(FTaskTimer::ScopeEvent taskTime(omp_get_thread_num(), &taskTimeRecorder, leafCells->getStartingIndex() * 20 * 8, "P2M"));
                KernelClass*const kernel = kernels[omp_get_thread_num()];

                for(int leafIdx = 0 ; leafIdx < leafCells->getNumberOfCellsInBlock() ; ++leafIdx){
                    multipole_t* leaf_multipole = &(leafCells->getMultipole(leafIdx));
                    symbolic_data_t* leaf_symbolic = &(leafCells->getSymbolic(leafIdx));
                    ParticleContainerClass particles
                        = containers->template getLeaf<ParticleContainerClass>(leafIdx);

                    FAssertLF(leafCells->getCellMortonIndex(leafIdx)
                              == containers->getLeafMortonIndex(leafIdx));

                    kernel->P2M(leaf_multipole, leaf_symbolic, &particles);
                }
            }
        }

        FLOG( FLog::Controller << "\t\t bottomPass in " << timer.tacAndElapsed() << "s\n" );
    }

    void upwardPass(){
        FLOG( FTic timer; );
        for(int idxLevel = FMath::Min(tree->getHeight() - 2, FAbstractAlgorithm::lowerWorkingLevel - 1) ; idxLevel >= FAbstractAlgorithm::upperWorkingLevel ; --idxLevel){
            typename OctreeClass::CellGroupIterator iterCells = tree->cellsBegin(idxLevel);
            const typename OctreeClass::CellGroupIterator endCells = tree->cellsEnd(idxLevel);

            typename OctreeClass::CellGroupIterator iterChildCells = tree->cellsBegin(idxLevel+1);
            const typename OctreeClass::CellGroupIterator endChildCells = tree->cellsEnd(idxLevel+1);

            while(iterCells != endCells){
                assert(iterChildCells != endChildCells);
                CellContainerClass*const currentCells = (*iterCells);
                PoleCellClass* cellPoles = currentCells->getRawMultipoleBuffer();

                CellContainerClass* subCellGroup = nullptr;
                PoleCellClass* subCellGroupPoles = nullptr;

                // Skip current group if needed
                if( (*iterChildCells)->getEndingIndex() <= (currentCells->getStartingIndex()<<3) ){
                    ++iterChildCells;
                    FAssertLF( iterChildCells != endChildCells );
                    FAssertLF( ((*iterChildCells)->getStartingIndex()>>3) == currentCells->getStartingIndex() );
                }

                while(true){
                    subCellGroup = (*iterChildCells);
                    subCellGroupPoles = (*iterChildCells)->getRawMultipoleBuffer();

                    #pragma omp task default(none) firstprivate(idxLevel, currentCells, cellPoles, subCellGroup, subCellGroupPoles) depend(commute_if_supported: cellPoles[0]) depend(in: subCellGroupPoles[0]) priority_if_supported(priorities.getInsertionPosM2M(idxLevel)) taskname_if_supported("M2M")
                    {
                        KernelClass*const kernel = kernels[omp_get_thread_num()];
                        const MortonIndex firstParent = FMath::Max(currentCells->getStartingIndex(), subCellGroup->getStartingIndex()>>3);
                        const MortonIndex lastParent = FMath::Min(currentCells->getEndingIndex()-1, (subCellGroup->getEndingIndex()-1)>>3);
                        FTIME_TASKS(FTaskTimer::ScopeEvent taskTime(omp_get_thread_num(), &taskTimeRecorder, ((lastParent * 20) + idxLevel) * 8 + 1, "M2M"));

                        int idxParentCell = currentCells->getCellIndex(firstParent);
                        FAssertLF(idxParentCell != -1);

                        int idxChildCell = subCellGroup->getFistChildIdx(firstParent);
                        FAssertLF(idxChildCell != -1);

                        while(true){
                            multipole_t* parent_multipole
                                = &(currentCells->getMultipole(idxParentCell));
                            const symbolic_data_t* parent_symbolic
                                = &(currentCells->getSymbolic(idxParentCell));
                            std::array<const multipole_t*, 8> child_multipoles {};
                            std::array<const symbolic_data_t*, 8> child_symbolics  {};

                            FAssertLF(parent_symbolic->getMortonIndex()
                                      == currentCells->getCellMortonIndex(idxParentCell));
                            FAssertLF(parent_symbolic->getMortonIndex()
                                      == (subCellGroup->getCellMortonIndex(idxChildCell)>>3));

                            do {
                                const int idxChild = ((subCellGroup->getCellMortonIndex(idxChildCell)) & 7);

                                FAssertLF(child_multipoles[idxChild] == nullptr);
                                FAssertLF(child_symbolics[idxChild] == nullptr);

                                child_multipoles[idxChild] = &(subCellGroup->getMultipole(idxChildCell));
                                child_symbolics[idxChild]  = &(subCellGroup->getSymbolic(idxChildCell));

                                FAssertLF(subCellGroup->getCellMortonIndex(idxChildCell)
                                          == child_symbolics[idxChild]->getMortonIndex());

                                idxChildCell += 1;
                            } while(idxChildCell != subCellGroup->getNumberOfCellsInBlock()
                                    && (parent_symbolic->getMortonIndex()
                                        == (subCellGroup->getCellMortonIndex(idxChildCell)>>3)));

                            kernel->M2M(
                                parent_multipole,
                                parent_symbolic,
                                child_multipoles.data(),
                                child_symbolics.data()
                                );

                            if(currentCells->getCellMortonIndex(idxParentCell) == lastParent){
                                break;
                            }

                            idxParentCell += 1;
                        }
                    }

                    if((*iterChildCells)->getEndingIndex() <= (((currentCells->getEndingIndex()-1)<<3)+7)
                            && (iterChildCells+1) != endChildCells
                            && (*(iterChildCells+1))->getStartingIndex() <= ((currentCells->getEndingIndex()-1)<<3)+7 ){
                        (++iterChildCells);
                    }
                    else{
                        break;
                    }
                }

                ++iterCells;
            }

            FAssertLF(iterCells == endCells);
            FAssertLF((iterChildCells == endChildCells || (++iterChildCells) == endChildCells));
            FAssertLF(iterCells == endCells && (iterChildCells == endChildCells || (++iterChildCells) == endChildCells));
        }
        FLOG( FLog::Controller << "\t\t upwardPass in " << timer.tacAndElapsed() << "s\n" );
    }

    void transferPass(const int startLevel, const int endLevel){
        FLOG( FTic timer; );
        FLOG( FTic timerInBlock; FTic timerOutBlock; );
        for(int idxLevel = startLevel ; idxLevel < endLevel ; ++idxLevel){
            FLOG( timerInBlock.tic() );
            {
                typename OctreeClass::CellGroupIterator iterCells = tree->cellsBegin(idxLevel);
                const typename OctreeClass::CellGroupIterator endCells = tree->cellsEnd(idxLevel);

                while(iterCells != endCells){
                    CellContainerClass* currentCells = (*iterCells);
                    PoleCellClass* cellPoles = currentCells->getRawMultipoleBuffer();
                    LocalCellClass* cellLocals = currentCells->getRawLocalBuffer();

#pragma omp task default(none) firstprivate(currentCells, cellPoles, cellLocals, idxLevel) depend(commute_if_supported: cellLocals[0]) depend(in: cellPoles[0])  priority_if_supported(priorities.getInsertionPosM2L(idxLevel)) taskname_if_supported("M2L")
                    {
                        FTIME_TASKS(FTaskTimer::ScopeEvent taskTime(omp_get_thread_num(), &taskTimeRecorder, ((currentCells->getStartingIndex() *20) + idxLevel ) * 8 + 2, "M2L"));
                        const MortonIndex blockStartIdx = currentCells->getStartingIndex();
                        const MortonIndex blockEndIdx   = currentCells->getEndingIndex();
                        KernelClass * const kernel = kernels[omp_get_thread_num()];
                        std::array<const multipole_t*, 189> inter_multipoles {};
                        std::array<const symbolic_data_t*, 189> inter_symbolics  {};

                        for(int cellIdx = 0 ; cellIdx < currentCells->getNumberOfCellsInBlock() ; ++cellIdx) {
                            local_expansion_t* target_local_expansion = &(currentCells->getLocalExpansion(cellIdx));
                            const symbolic_data_t* target_symbolic = &(currentCells->getSymbolic(cellIdx));

                            FAssertLF(target_symbolic->getMortonIndex()
                                      == currentCells->getCellMortonIndex(cellIdx));

                            MortonIndex interactionsIndexes[189];
                            int interactionsPosition[189];
                            const FTreeCoordinate coord(target_symbolic->getCoordinate());
                            int counter = coord.getInteractionNeighbors(idxLevel,interactionsIndexes,interactionsPosition);

                            int counterExistingCell = 0;

                            for(int idxInter = 0 ; idxInter < counter ; ++idxInter){
                                if( blockStartIdx <= interactionsIndexes[idxInter] && interactionsIndexes[idxInter] < blockEndIdx ){
                                    const int cellPos = currentCells->getCellIndex(interactionsIndexes[idxInter]);
                                    if(cellPos != -1){
                                        const multipole_t* inter_multipole
                                            = &(currentCells->getMultipole(cellPos));
                                        const symbolic_data_t* inter_symbolic
                                            = &(currentCells->getSymbolic(cellPos));

                                        FAssertLF(inter_symbolic->getMortonIndex()
                                                  == interactionsIndexes[idxInter]);

                                        interactionsPosition[counterExistingCell]
                                            = interactionsPosition[idxInter];
                                        inter_multipoles[counterExistingCell] = inter_multipole;
                                        inter_symbolics[counterExistingCell] = inter_symbolic;
                                        counterExistingCell += 1;
                                    }
                                }
                            }

                            kernel->M2L(
                                target_local_expansion,
                                target_symbolic,
                                inter_multipoles.data(),
                                inter_symbolics.data(),
                                interactionsPosition,
                                counterExistingCell
                                );
                        }
                    }
                    ++iterCells;
                }
            }
            FLOG( timerInBlock.tac() );
            FLOG( timerOutBlock.tic() );
            {
                typename OctreeClass::CellGroupIterator iterCells = tree->cellsBegin(idxLevel);
                const typename OctreeClass::CellGroupIterator endCells = tree->cellsEnd(idxLevel);

                typename std::vector<std::vector<BlockInteractions<CellContainerClass>>>::iterator externalInteractionsIter = externalInteractionsAllLevel[idxLevel].begin();

                while(iterCells != endCells){
                    CellContainerClass* currentCells = (*iterCells);
                    PoleCellClass* cellPoles = currentCells->getRawMultipoleBuffer();
                    LocalCellClass* cellLocals = currentCells->getRawLocalBuffer();

                    typename std::vector<BlockInteractions<CellContainerClass>>::iterator currentInteractions = (*externalInteractionsIter).begin();
                    const typename std::vector<BlockInteractions<CellContainerClass>>::iterator currentInteractionsEnd = (*externalInteractionsIter).end();

                    while(currentInteractions != currentInteractionsEnd){
                        CellContainerClass* cellsOther = (*currentInteractions).otherBlock;
                        PoleCellClass* cellOtherPoles = cellsOther->getRawMultipoleBuffer();
                        LocalCellClass* cellOtherLocals = cellsOther->getRawLocalBuffer();
                        const std::vector<OutOfBlockInteraction>* outsideInteractions = &(*currentInteractions).interactions;

                        #pragma omp task default(none) firstprivate(currentCells, cellLocals, outsideInteractions, cellsOther, cellOtherPoles, idxLevel) depend(commute_if_supported: cellLocals[0]) depend(in: cellOtherPoles[0])  priority_if_supported(priorities.getInsertionPosM2LExtern(idxLevel)) taskname_if_supported("M2L-out")
                        {
                            FTIME_TASKS(FTaskTimer::ScopeEvent taskTime(omp_get_thread_num(), &taskTimeRecorder, (((currentCells->getStartingIndex()+1) * (cellsOther->getStartingIndex()+2)) * 20 + idxLevel) * 8 + 3, "M2L-ext"));
                            KernelClass*const kernel = kernels[omp_get_thread_num()];

                            for(int outInterIdx = 0 ; outInterIdx < int(outsideInteractions->size()) ; ++outInterIdx){
                                const auto& inter_data = (*outsideInteractions)[outInterIdx];
                                const auto inter_block_out_idx = inter_data.outsideIdxInBlock;
                                const multipole_t* inter_multipole
                                    = &(cellsOther->getMultipole(inter_block_out_idx));
                                const symbolic_data_t* inter_symbolic
                                    = &(cellsOther->getSymbolic(inter_block_out_idx));

                                FAssertLF(inter_symbolic->getMortonIndex() == inter_data.outIndex);

                                const auto cell_block_in_idx = inter_data.insideIdxInBlock;
                                local_expansion_t* cell_local_exp =
                                    &(currentCells->getLocalExpansion(cell_block_in_idx));
                                const symbolic_data_t* cell_symbolic =
                                    &(currentCells->getSymbolic(cell_block_in_idx));

                                FAssertLF(cell_symbolic->getMortonIndex() == inter_data.insideIndex);

                                kernel->M2L(
                                    cell_local_exp,
                                    cell_symbolic,
                                    &inter_multipole,
                                    &inter_symbolic,
                                    &inter_data.relativeOutPosition,
                                    1
                                    );
                            }
                        }

                        #pragma omp task default(none) firstprivate(currentCells, cellPoles, outsideInteractions, cellsOther, cellOtherLocals, idxLevel) depend(commute_if_supported: cellOtherLocals[0]) depend(in: cellPoles[0])  priority_if_supported(priorities.getInsertionPosM2LExtern(idxLevel)) taskname_if_supported("M2L-out")
                        {
                            FTIME_TASKS(FTaskTimer::ScopeEvent taskTime(omp_get_thread_num(), &taskTimeRecorder, (((currentCells->getStartingIndex()+1) * (cellsOther->getStartingIndex()+1)) * 20 + idxLevel) * 8 + 3, "M2L-ext"));
                            KernelClass*const kernel = kernels[omp_get_thread_num()];

                            for(int outInterIdx = 0 ; outInterIdx < int(outsideInteractions->size()) ; ++outInterIdx){
                                const auto& inter_data = (*outsideInteractions)[outInterIdx];
                                const auto inter_idx = inter_data.outsideIdxInBlock;
                                local_expansion_t* inter_local_expansion
                                    = &(cellsOther->getLocalExpansion(inter_idx));
                                const symbolic_data_t* inter_symbolic
                                    = &(cellsOther->getSymbolic(inter_idx));

                                FAssertLF(inter_symbolic->getMortonIndex() == inter_data.outIndex);

                                const auto cell_idx = inter_data.insideIdxInBlock;
                                const multipole_t* cell_multipole
                                    = &(currentCells->getMultipole(cell_idx));
                                const symbolic_data_t* cell_symbolic
                                    = &(currentCells->getSymbolic(cell_idx));

                                FAssertLF(cell_symbolic->getMortonIndex() == inter_data.insideIndex);

                                const int otherPos = getOppositeInterIndex(inter_data.relativeOutPosition);
                                kernel->M2L(
                                    inter_local_expansion,
                                    inter_symbolic,
                                    &cell_multipole,
                                    &cell_symbolic,
                                    &otherPos,
                                    1);
                            }
                        }

                        ++currentInteractions;
                    }

                    ++iterCells;
                    ++externalInteractionsIter;
                }
            }
            FLOG( timerOutBlock.tac() );
        }
        FLOG( FLog::Controller << "\t\t transferPass in " << timer.tacAndElapsed() << "s\n" );
        FLOG( FLog::Controller << "\t\t\t inblock in  " << timerInBlock.elapsed() << "s\n" );
        FLOG( FLog::Controller << "\t\t\t outblock in " << timerOutBlock.elapsed() << "s\n" );
    }

    void downardPass(){
        FLOG( FTic timer; );
        for(int idxLevel = FAbstractAlgorithm::upperWorkingLevel ; idxLevel < FAbstractAlgorithm::lowerWorkingLevel - 1 ; ++idxLevel){
            typename OctreeClass::CellGroupIterator iterCells = tree->cellsBegin(idxLevel);
            const typename OctreeClass::CellGroupIterator endCells = tree->cellsEnd(idxLevel);

            typename OctreeClass::CellGroupIterator iterChildCells = tree->cellsBegin(idxLevel+1);
            const typename OctreeClass::CellGroupIterator endChildCells = tree->cellsEnd(idxLevel+1);

            while(iterCells != endCells){
                assert(iterChildCells != endChildCells);
                CellContainerClass*const currentCells = (*iterCells);
                LocalCellClass* cellLocals = currentCells->getRawLocalBuffer();

                CellContainerClass* subCellGroup = nullptr;
                LocalCellClass* subCellLocalGroupsLocal = nullptr;

                // Skip current group if needed
                if( (*iterChildCells)->getEndingIndex() <= (currentCells->getStartingIndex()<<3) ){
                    ++iterChildCells;
                    FAssertLF( iterChildCells != endChildCells );
                    FAssertLF( ((*iterChildCells)->getStartingIndex()>>3) == currentCells->getStartingIndex() );
                }

                while(true){
                    subCellGroup = (*iterChildCells);
                    subCellLocalGroupsLocal = (*iterChildCells)->getRawLocalBuffer();

                    if(noCommuteAtLastLevel == false || idxLevel != FAbstractAlgorithm::lowerWorkingLevel - 2){
                        #pragma omp task default(none) firstprivate(idxLevel, currentCells, cellLocals, subCellGroup, subCellLocalGroupsLocal) depend(commute_if_supported: subCellLocalGroupsLocal[0]) depend(in: cellLocals[0])  priority_if_supported(priorities.getInsertionPosL2L(idxLevel)) taskname_if_supported("L2L")
                        {
                            KernelClass*const kernel = kernels[omp_get_thread_num()];

                            const MortonIndex firstParent = FMath::Max(currentCells->getStartingIndex(), subCellGroup->getStartingIndex()>>3);
                            const MortonIndex lastParent = FMath::Min(currentCells->getEndingIndex()-1, (subCellGroup->getEndingIndex()-1)>>3);
                            FTIME_TASKS(FTaskTimer::ScopeEvent taskTime(omp_get_thread_num(), &taskTimeRecorder, ((lastParent * 20) + idxLevel) * 8 + 4, "L2L"));

                            int idxParentCell = currentCells->getCellIndex(firstParent);
                            FAssertLF(idxParentCell != -1);

                            int idxChildCell = subCellGroup->getFistChildIdx(firstParent);
                            FAssertLF(idxChildCell != -1);

                            while(true){
                                const local_expansion_t* parent_local_expansion
                                    = &(currentCells->getLocalExpansion(idxParentCell));
                                const symbolic_data_t* parent_symbolic
                                    = &(currentCells->getSymbolic(idxParentCell));

                                FAssertLF(parent_symbolic->getMortonIndex()
                                          == currentCells->getCellMortonIndex(idxParentCell));

                                std::array<local_expansion_t*, 8> child_local_expansions {};
                                std::array<const symbolic_data_t*, 8> child_symbolics {};

                                FAssertLF(parent_symbolic->getMortonIndex() == (subCellGroup->getCellMortonIndex(idxChildCell)>>3));

                                do{
                                    const int idxChild = ((subCellGroup->getCellMortonIndex(idxChildCell)) & 7);

                                    FAssertLF(child_local_expansions[idxChild] == nullptr);
                                    FAssertLF(child_symbolics[idxChild] == nullptr);

                                    child_local_expansions[idxChild]
                                        = &(subCellGroup->getLocalExpansion(idxChildCell));
                                    child_symbolics[idxChild]
                                        = &(subCellGroup->getSymbolic(idxChildCell));

                                    FAssertLF(subCellGroup->getCellMortonIndex(idxChildCell)
                                              == child_symbolics[idxChild]->getMortonIndex());

                                    idxChildCell += 1;
                                } while(idxChildCell != subCellGroup->getNumberOfCellsInBlock()
                                        && (parent_symbolic->getMortonIndex()
                                            == (subCellGroup->getCellMortonIndex(idxChildCell)>>3)));

                                kernel->L2L(
                                    parent_local_expansion,
                                    parent_symbolic,
                                    child_local_expansions.data(),
                                    child_symbolics.data()
                                    );

                                if(currentCells->getCellMortonIndex(idxParentCell) == lastParent){
                                    break;
                                }

                                idxParentCell += 1;
                            }
                        }
                    }
                    else{
                        #pragma omp task default(none) firstprivate(idxLevel, currentCells, cellLocals, subCellGroup, subCellLocalGroupsLocal) depend(inout: subCellLocalGroupsLocal[0]) depend(in: cellLocals[0])  priority_if_supported(priorities.getInsertionPosL2L(idxLevel)) taskname_if_supported("L2L")
                        {
                            KernelClass*const kernel = kernels[omp_get_thread_num()];

                            const MortonIndex firstParent = FMath::Max(currentCells->getStartingIndex(), subCellGroup->getStartingIndex()>>3);
                            const MortonIndex lastParent = FMath::Min(currentCells->getEndingIndex()-1, (subCellGroup->getEndingIndex()-1)>>3);
                            FTIME_TASKS(FTaskTimer::ScopeEvent taskTime(omp_get_thread_num(), &taskTimeRecorder, ((lastParent * 20) + idxLevel) * 8 + 4, "L2L"));

                            int idxParentCell = currentCells->getCellIndex(firstParent);
                            FAssertLF(idxParentCell != -1);

                            int idxChildCell = subCellGroup->getFistChildIdx(firstParent);
                            FAssertLF(idxChildCell != -1);

                            while(true){
                                const local_expansion_t* parent_local_expansion
                                    = &(currentCells->getLocalExpansion(idxParentCell));
                                const symbolic_data_t* parent_symbolic
                                    = &(currentCells->getSymbolic(idxParentCell));

                                FAssertLF(parent_symbolic->getMortonIndex()
                                          == currentCells->getCellMortonIndex(idxParentCell));

                                std::array<local_expansion_t*, 8> child_local_expansions {};
                                std::array<const symbolic_data_t*, 8> child_symbolics {};

                                FAssertLF(parent_symbolic->getMortonIndex()
                                          == (subCellGroup->getCellMortonIndex(idxChildCell)>>3));

                                do {
                                    const int idxChild = ((subCellGroup->getCellMortonIndex(idxChildCell)) & 7);

                                    FAssertLF(child_symbolics[idxChild] == nullptr);
                                    FAssertLF(child_local_expansions[idxChild] == nullptr);

                                    child_local_expansions[idxChild]
                                        = &(subCellGroup->getLocalExpansion(idxChildCell));
                                    child_symbolics[idxChild]
                                        = &(subCellGroup->getSymbolic(idxChildCell));

                                    FAssertLF(subCellGroup->getCellMortonIndex(idxChildCell)
                                              == child_symbolics[idxChild]->getMortonIndex());

                                    idxChildCell += 1;
                                } while(idxChildCell != subCellGroup->getNumberOfCellsInBlock()
                                        && (parent_symbolic->getMortonIndex()
                                            == (subCellGroup->getCellMortonIndex(idxChildCell)>>3)));

                                kernel->L2L(
                                    parent_local_expansion,
                                    parent_symbolic,
                                    child_local_expansions.data(),
                                    child_symbolics.data()
                                    );

                                if(currentCells->getCellMortonIndex(idxParentCell) == lastParent){
                                    break;
                                }

                                idxParentCell += 1;
                            }
                        }
                    }

                    if((*iterChildCells)->getEndingIndex() <= (((currentCells->getEndingIndex()-1)<<3)+7)
                            && (iterChildCells+1) != endChildCells
                            && (*(iterChildCells+1))->getStartingIndex() <= ((currentCells->getEndingIndex()-1)<<3)+7){
                        (++iterChildCells);
                    }
                    else{
                        break;
                    }
                }



                ++iterCells;
            }

            FAssertLF(iterCells == endCells && (iterChildCells == endChildCells || (++iterChildCells) == endChildCells));
        }
        FLOG( FLog::Controller << "\t\t downardPass in " << timer.tacAndElapsed() << "s\n" );
    }

    void directPass(){
        FLOG( FTic timer; );
        FLOG( FTic timerInBlock; FTic timerOutBlock; );

        FLOG( timerOutBlock.tic() );
        {
            typename OctreeClass::ParticleGroupIterator iterParticles = tree->leavesBegin();
            const typename OctreeClass::ParticleGroupIterator endParticles = tree->leavesEnd();

            typename std::vector<std::vector<BlockInteractions<ParticleGroupClass>>>::iterator externalInteractionsIter = externalInteractionsLeafLevel.begin();

            while(iterParticles != endParticles){
                typename std::vector<BlockInteractions<ParticleGroupClass>>::iterator currentInteractions = (*externalInteractionsIter).begin();
                const typename std::vector<BlockInteractions<ParticleGroupClass>>::iterator currentInteractionsEnd = (*externalInteractionsIter).end();

                ParticleGroupClass* containers = (*iterParticles);
                unsigned char* containersDown = containers->getRawAttributesBuffer();

                while(currentInteractions != currentInteractionsEnd){
                    ParticleGroupClass* containersOther = (*currentInteractions).otherBlock;
                    unsigned char* containersOtherDown = containersOther->getRawAttributesBuffer();
                    const std::vector<OutOfBlockInteraction>* outsideInteractions = &(*currentInteractions).interactions;

#pragma omp task default(none) firstprivate(containers, containersDown, containersOther, containersOtherDown, outsideInteractions) depend(commute_if_supported: containersOtherDown[0], containersDown[0])  priority_if_supported(priorities.getInsertionPosP2PExtern()) taskname_if_supported("P2P-out")
                    {
                        FTIME_TASKS(FTaskTimer::ScopeEvent taskTime(omp_get_thread_num(), &taskTimeRecorder, ((containersOther->getStartingIndex()+1) * (containers->getStartingIndex()+1))*20*8 + 6, "P2P-ext"));
                        KernelClass*const kernel = kernels[omp_get_thread_num()];
                        for(int outInterIdx = 0 ; outInterIdx < int(outsideInteractions->size()) ; ++outInterIdx){
                            ParticleContainerClass interParticles = containersOther->template getLeaf<ParticleContainerClass>((*outsideInteractions)[outInterIdx].outsideIdxInBlock);
                            ParticleContainerClass particles = containers->template getLeaf<ParticleContainerClass>((*outsideInteractions)[outInterIdx].insideIdxInBlock);

                            FAssertLF(containersOther->getLeafMortonIndex((*outsideInteractions)[outInterIdx].outsideIdxInBlock) == (*outsideInteractions)[outInterIdx].outIndex);
                            FAssertLF(containers->getLeafMortonIndex((*outsideInteractions)[outInterIdx].insideIdxInBlock) == (*outsideInteractions)[outInterIdx].insideIndex);

                            ParticleContainerClass* ptrLeaf = &interParticles;
                            kernel->P2POuter( FTreeCoordinate((*outsideInteractions)[outInterIdx].insideIndex),
                                                &particles , &ptrLeaf, &(*outsideInteractions)[outInterIdx].relativeOutPosition, 1);
                            const int otherPosition = getOppositeNeighIndex((*outsideInteractions)[outInterIdx].relativeOutPosition);
                            ptrLeaf = &particles;
                            kernel->P2POuter( FTreeCoordinate((*outsideInteractions)[outInterIdx].outIndex),
                                                &interParticles , &ptrLeaf, &otherPosition, 1);
                        }
                    }

                    ++currentInteractions;
                }

                ++iterParticles;
                ++externalInteractionsIter;
            }
        }
        FLOG( timerOutBlock.tac() );
        FLOG( timerInBlock.tic() );
        {
            typename OctreeClass::ParticleGroupIterator iterParticles = tree->leavesBegin();
            const typename OctreeClass::ParticleGroupIterator endParticles = tree->leavesEnd();

            while(iterParticles != endParticles){
                ParticleGroupClass* containers = (*iterParticles);
                unsigned char* containersDown = containers->getRawAttributesBuffer();

                #pragma omp task default(none) firstprivate(containers, containersDown) depend(commute_if_supported: containersDown[0])  priority_if_supported(priorities.getInsertionPosP2P()) taskname_if_supported("P2P")
                {
                    FTIME_TASKS(FTaskTimer::ScopeEvent taskTime(omp_get_thread_num(), &taskTimeRecorder, containers->getStartingIndex()*20*8 + 5, "P2P"));
                    const MortonIndex blockStartIdx = containers->getStartingIndex();
                    const MortonIndex blockEndIdx = containers->getEndingIndex();
                    KernelClass*const kernel = kernels[omp_get_thread_num()];

                    for(int leafIdx = 0 ; leafIdx < containers->getNumberOfLeavesInBlock() ; ++leafIdx){
                        ParticleContainerClass particles = containers->template getLeaf<ParticleContainerClass>(leafIdx);
                        const MortonIndex mindex = containers->getLeafMortonIndex(leafIdx);

                        MortonIndex interactionsIndexes[26];
                        int interactionsPosition[26];
                        FTreeCoordinate coord(mindex);
                        int counter = coord.getNeighborsIndexes(tree->getHeight(),interactionsIndexes,interactionsPosition);

                        ParticleContainerClass interactionsObjects[26];
                        ParticleContainerClass* interactions[26];
                        int counterExistingCell = 0;

                        for(int idxInter = 0 ; idxInter < counter ; ++idxInter){
                            if( blockStartIdx <= interactionsIndexes[idxInter] && interactionsIndexes[idxInter] < blockEndIdx ){
                                const int leafPos = containers->getLeafIndex(interactionsIndexes[idxInter]);
                                if(leafPos != -1){
                                    interactionsObjects[counterExistingCell] = containers->template getLeaf<ParticleContainerClass>(leafPos);
                                    interactionsPosition[counterExistingCell] = interactionsPosition[idxInter];
                                    interactions[counterExistingCell] = &interactionsObjects[counterExistingCell];
                                    counterExistingCell += 1;
                                }
                            }
                        }

                        kernel->P2P( coord, &particles, &particles , interactions, interactionsPosition, counterExistingCell);
                    }
                }
                ++iterParticles;
            }
        }
        FLOG( timerInBlock.tac() );

        FLOG( FLog::Controller << "\t\t directPass in " << timer.tacAndElapsed() << "s\n" );
        FLOG( FLog::Controller << "\t\t\t inblock  in " << timerInBlock.elapsed() << "s\n" );
        FLOG( FLog::Controller << "\t\t\t outblock in " << timerOutBlock.elapsed() << "s\n" );
    }

    void mergePass(){
        FLOG( FTic timer; );

        for(int idxGroup = 0 ; idxGroup < tree->getNbParticleGroup() ; ++idxGroup){
            CellContainerClass* leafCells  = tree->getCellGroup(tree->getHeight()-1, idxGroup);
            LocalCellClass* cellLocals = leafCells->getRawLocalBuffer();

            ParticleGroupClass* containers = tree->getParticleGroup(idxGroup);
            unsigned char* containersDown = containers->getRawAttributesBuffer();

            #pragma omp task default(shared) firstprivate(leafCells, cellLocals, containers, containersDown) depend(commute_if_supported: containersDown[0]) depend(in: cellLocals[0])  priority_if_supported(priorities.getInsertionPosL2P()) taskname_if_supported("L2P")
            {
                FTIME_TASKS(FTaskTimer::ScopeEvent taskTime(omp_get_thread_num(), &taskTimeRecorder, (leafCells->getStartingIndex()*20*8) + 7, "L2P"));
                KernelClass*const kernel = kernels[omp_get_thread_num()];

                for(int cellIdx = 0 ; cellIdx < leafCells->getNumberOfCellsInBlock() ; ++cellIdx){
                    const local_expansion_t* leaf_local_expansion
                        = &(leafCells->getLocalExpansion(cellIdx));
                    const symbolic_data_t* leaf_symbolic
                        = &(leafCells->getSymbolic(cellIdx));
                    ParticleContainerClass particles
                        = containers->template getLeaf<ParticleContainerClass>(cellIdx);

                    FAssertLF(leaf_symbolic->getMortonIndex()
                              == leafCells->getCellMortonIndex(cellIdx));
                    FAssertLF(leafCells->getCellMortonIndex(cellIdx)
                              == containers->getLeafMortonIndex(cellIdx));

                    kernel->L2P(
                        leaf_local_expansion,
                        leaf_symbolic,
                        &particles
                        );
                }
            }
        }

        FLOG( FLog::Controller << "\t\t L2P in " << timer.tacAndElapsed() << "s\n" );
    }

    int getOppositeNeighIndex(const int index) const {
        // ((idxX+1)*3 + (idxY+1)) * 3 + (idxZ+1)
        return 27-index-1;
    }

    int getOppositeInterIndex(const int index) const {
        // ((( (xdiff+3) * 7) + (ydiff+3))) * 7 + zdiff + 3
        return 343-index-1;
    }
};

#endif // FGROUPTASKDEPALGORITHM_HPP
