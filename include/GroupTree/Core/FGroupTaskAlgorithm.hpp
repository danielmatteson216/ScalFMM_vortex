
// Keep in private GIT
#ifndef FGROUPTASKALGORITHM_HPP
#define FGROUPTASKALGORITHM_HPP

#include "Utils/FGlobal.hpp"
#include "Core/FCoreCommon.hpp"
#include "Utils/FQuickSort.hpp"
#include "Containers/FTreeCoordinate.hpp"
#include "Utils/FLog.hpp"
#include "Utils/FTic.hpp"
#include "Utils/FAlgorithmTimers.hpp"
#include "FOutOfBlockInteraction.hpp"

#include <vector>
#include <vector>

#include <omp.h>

template <class OctreeClass, class CellContainerClass, class KernelClass, class ParticleGroupClass, class ParticleContainerClass>
class FGroupTaskAlgorithm : public FAbstractAlgorithm, public FAlgorithmTimers {
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

public:
    using multipole_t = typename std::remove_pointer<typename OctreeClass::CellGroupIterator::value_type>::type::multipole_t;
    using local_expansion_t = typename std::remove_pointer<typename OctreeClass::CellGroupIterator::value_type>::type::local_expansion_t;
    using symbolic_data_t = typename std::remove_pointer<typename OctreeClass::CellGroupIterator::value_type>::type::symbolic_data_t;

    FGroupTaskAlgorithm(OctreeClass*const inTree, KernelClass* inKernels)
        : tree(inTree), kernels(nullptr)
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
          #pragma omp critical (FGroupTaskAlgorithm_InitKernels)
            this->kernels[omp_get_thread_num()] = new KernelClass(*inKernels);
        }

        rebuildInteractions();

        FLOG(FLog::Controller << "FGroupTaskAlgorithm (Max Thread " << MaxThreads << ")\n");
    }

    ~FGroupTaskAlgorithm(){
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
        FLOG( FLog::Controller << "\tStart FGroupTaskAlgorithm\n" );
        Timers[P2MTimer].tic();
        #pragma omp parallel num_threads(MaxThreads)
        {
            #pragma omp single nowait
            {
                if(operationsToProceed & FFmmP2M) bottomPass();

                if(operationsToProceed & FFmmM2M) upwardPass();

                if(operationsToProceed & FFmmM2L) transferPass();

                if(operationsToProceed & FFmmL2L) downardPass();

            }

            #pragma omp single nowait
            {
                if( operationsToProceed & FFmmP2P ) directPass();
            }

            #pragma omp barrier

            #pragma omp single nowait
            {
                if( operationsToProceed & FFmmL2P ) mergePass();
            }
        }
        Timers[P2MTimer].tac();
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

                    // FIXME: hack around a clang bug
                    // it apparently can't manage a firstprivate const member
                    // such as 'tree', but can manage a local copy...
                    // getting the height first, or making 'tree' shared both
                    // workaround it.
                    OctreeClass * const treeAlias = tree;
                    #pragma omp task default(none) firstprivate(idxGroup, containers, externalInteractions, treeAlias)
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
                            int counter = coord.getNeighborsIndexes(treeAlias->getHeight(),interactionsIndexes,interactionsPosition);

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
                            ParticleGroupClass* leftContainers = treeAlias->getParticleGroup(idxLeftGroup);
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

                        // FIXME: hack around a clang bug
                        // it apparently can't manage a firstprivate const member
                        // such as 'tree', but can manage a local copy...
                        // getting the height first, or making 'tree' shared both
                        // workaround it.
                        OctreeClass * const treeAlias = tree;
                        #pragma omp task default(none) firstprivate(idxGroup, currentCells, idxLevel, externalInteractions, treeAlias)
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
                                CellContainerClass* leftCells   = treeAlias->getCellGroup(idxLevel, idxLeftGroup);
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
            ParticleGroupClass* containers = tree->getParticleGroup(idxGroup);
            #pragma omp task default(shared) firstprivate(leafCells, containers)
            {
                KernelClass*const kernel = kernels[omp_get_thread_num()];

                for(int leafIdx = 0 ; leafIdx < leafCells->getNumberOfCellsInBlock() ; ++leafIdx){
                    multipole_t* leaf_multipole = &leafCells->getMultipole(leafIdx);
                    symbolic_data_t* leaf_symbolic = &leafCells->getSymbolic(leafIdx);
                    ParticleContainerClass particles
                        = containers->template getLeaf<ParticleContainerClass>(leafIdx);

                    FAssertLF(leafCells->getCellMortonIndex(leafIdx) == containers->getLeafMortonIndex(leafIdx));

                    kernel->P2M(leaf_multipole, leaf_symbolic, &particles);
                }
            }
        }
        // Wait for task to complete
        #pragma omp taskwait

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

                CellContainerClass* subCellGroups[9];
                memset(subCellGroups, 0, sizeof(CellContainerClass*) * 9);

                // Skip current group if needed
                if( (*iterChildCells)->getEndingIndex() <= (currentCells->getStartingIndex()<<3) ){
                    ++iterChildCells;
                    FAssertLF( iterChildCells != endChildCells );
                    FAssertLF( ((*iterChildCells)->getStartingIndex()>>3) == currentCells->getStartingIndex() );
                }
                // Copy at max 8 groups
                int nbSubCellGroups = 0;
                subCellGroups[nbSubCellGroups] = (*iterChildCells);
                nbSubCellGroups += 1;
                while((*iterChildCells)->getEndingIndex() <= (((currentCells->getEndingIndex()-1)<<3)+7)
                      && (++iterChildCells) != endChildCells
                      && (*iterChildCells)->getStartingIndex() <= ((currentCells->getEndingIndex()-1)<<3)+7 ){
                    FAssertLF( nbSubCellGroups < 9 );
                    subCellGroups[nbSubCellGroups] = (*iterChildCells);
                    nbSubCellGroups += 1;
                }

                #pragma omp task default(none) firstprivate(idxLevel, currentCells, subCellGroups, nbSubCellGroups, kernels)
                {
                    KernelClass*const kernel = kernels[omp_get_thread_num()];
                    int idxSubCellGroup = 0;
                    int idxChildCell = subCellGroups[0]->getFistChildIdx(currentCells->getCellMortonIndex(0));
                    FAssertLF(idxChildCell != -1);

                    for(int cellIdx = 0 ; cellIdx < currentCells->getNumberOfCellsInBlock() ; ++cellIdx){
                        multipole_t* parent_multipole = &currentCells->getMultipole(cellIdx);
                        symbolic_data_t* parent_symbolic  = &currentCells->getSymbolic(cellIdx);

                        FAssertLF(parent_symbolic->getMortonIndex()
                                  == currentCells->getCellMortonIndex(cellIdx));

                        std::array<const multipole_t*,8> child_multipoles {};
                        std::array<const symbolic_data_t*,8> child_symbolics {};

                        FAssertLF(idxSubCellGroup != nbSubCellGroups);

                        while(idxSubCellGroup != nbSubCellGroups
                              && ((subCellGroups[idxSubCellGroup]->getCellMortonIndex(idxChildCell)>>3)
                                  == parent_symbolic->getMortonIndex()))
                        {
                            const int idxChild = ((subCellGroups[idxSubCellGroup]->getCellMortonIndex(idxChildCell)) & 7);
                            FAssertLF(child_multipoles[idxChild] == nullptr);
                            FAssertLF(child_symbolics[idxChild] == nullptr);

                            child_multipoles[idxChild]
                                = &(subCellGroups[idxSubCellGroup]->getMultipole(idxChildCell));
                            child_symbolics [idxChild]
                                = &(subCellGroups[idxSubCellGroup]->getSymbolic(idxChildCell));

                            FAssertLF(subCellGroups[idxSubCellGroup]->getCellMortonIndex(idxChildCell)
                                      == child_symbolics[idxChild]->getMortonIndex());

                            idxChildCell += 1;
                            if(idxChildCell == subCellGroups[idxSubCellGroup]->getNumberOfCellsInBlock()){
                                idxChildCell = 0;
                                idxSubCellGroup += 1;
                            }
                        }

                        kernel->M2M(
                            parent_multipole,
                            parent_symbolic,
                            child_multipoles.data(),
                            child_symbolics.data()
                            );
                    }
                }

                ++iterCells;
            }
            // Wait this level before the next one
            #pragma omp taskwait

            FAssertLF(iterCells == endCells);
            FAssertLF((iterChildCells == endChildCells || (++iterChildCells) == endChildCells));
            FAssertLF(iterCells == endCells && (iterChildCells == endChildCells || (++iterChildCells) == endChildCells));
        }
        FLOG( FLog::Controller << "\t\t upwardPass in " << timer.tacAndElapsed() << "s\n" );
    }

    void transferPass(){
        FLOG( FTic timer; );
        FLOG( FTic timerInBlock; FTic timerOutBlock; );
        for(int idxLevel = FAbstractAlgorithm::lowerWorkingLevel-1 ; idxLevel >= FAbstractAlgorithm::upperWorkingLevel ; --idxLevel){
            FLOG( timerInBlock.tic() );
            {
                typename OctreeClass::CellGroupIterator iterCells = tree->cellsBegin(idxLevel);
                const typename OctreeClass::CellGroupIterator endCells = tree->cellsEnd(idxLevel);

                while(iterCells != endCells){
                    CellContainerClass* currentCells = (*iterCells);

                    #pragma omp task default(none) firstprivate(currentCells, idxLevel, kernels)
                    {
                        const MortonIndex blockStartIdx = currentCells->getStartingIndex();
                        const MortonIndex blockEndIdx = currentCells->getEndingIndex();
                        KernelClass*const kernel = kernels[omp_get_thread_num()];

                        std::array<const multipole_t*, 189> neighbour_multipoles {};
                        std::array<const symbolic_data_t*, 189> neighbour_symbolics {};

                        for(int cellIdx = 0; cellIdx < currentCells->getNumberOfCellsInBlock(); ++cellIdx) {
                            local_expansion_t* target_local_exp = &currentCells->getLocalExpansion(cellIdx);
                            symbolic_data_t* target_symbolic  = &currentCells->getSymbolic(cellIdx);

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
                                        const multipole_t* neighbour_multipole
                                            = &currentCells->getMultipole(cellPos);
                                        const symbolic_data_t* neighbour_symbolic
                                            = &currentCells->getSymbolic(cellPos);

                                        FAssertLF(neighbour_symbolic->getMortonIndex()
                                                  == interactionsIndexes[idxInter]);

                                        interactionsPosition[counterExistingCell] = interactionsPosition[idxInter];
                                        neighbour_multipoles[counterExistingCell] = neighbour_multipole;
                                        neighbour_symbolics [counterExistingCell] = neighbour_symbolic;
                                        counterExistingCell += 1;
                                    }
                                }
                            }

                            kernel->M2L(
                                target_local_exp,
                                target_symbolic,
                                neighbour_multipoles.data(),
                                neighbour_symbolics.data(),
                                interactionsPosition,
                                counterExistingCell
                                );
                        }
                    }
                    ++iterCells;
                }
                #pragma omp taskwait
            }
            FLOG( timerInBlock.tac() );
            FLOG( timerOutBlock.tic() );
            {
                typename OctreeClass::CellGroupIterator iterCells = tree->cellsBegin(idxLevel);
                const typename OctreeClass::CellGroupIterator endCells = tree->cellsEnd(idxLevel);

                typename std::vector<std::vector<BlockInteractions<CellContainerClass>>>::iterator externalInteractionsIter = externalInteractionsAllLevel[idxLevel].begin();

                while(iterCells != endCells){
                    CellContainerClass* currentCells = (*iterCells);

                    typename std::vector<BlockInteractions<CellContainerClass>>::iterator currentInteractions = (*externalInteractionsIter).begin();
                    const typename std::vector<BlockInteractions<CellContainerClass>>::iterator currentInteractionsEnd = (*externalInteractionsIter).end();

                    while(currentInteractions != currentInteractionsEnd){
                        CellContainerClass* cellsOther = (*currentInteractions).otherBlock;
                        const std::vector<OutOfBlockInteraction>* outsideInteractions = &(*currentInteractions).interactions;

                        #pragma omp task default(none) firstprivate(currentCells, outsideInteractions, cellsOther, idxLevel, kernels)
                        {
                            KernelClass*const kernel = kernels[omp_get_thread_num()];

                            for(int outInterIdx = 0 ; outInterIdx < int(outsideInteractions->size()) ; ++outInterIdx) {
                                const auto& inter_data = (*outsideInteractions)[outInterIdx];
                                const symbolic_data_t* inter_symbolic
                                    = & cellsOther->getSymbolic(inter_data.outsideIdxInBlock);
                                const multipole_t* inter_multipole
                                    = & cellsOther->getMultipole(inter_data.outsideIdxInBlock);
                                local_expansion_t* inter_local_exp
                                    = & cellsOther->getLocalExpansion(inter_data.outsideIdxInBlock);

                                FAssertLF(inter_symbolic->getMortonIndex()
                                          == inter_data.outIndex);

                                const symbolic_data_t* cell_symbolic
                                    = &currentCells->getSymbolic(inter_data.insideIdxInBlock);
                                const multipole_t* cell_multipole
                                    = &currentCells->getMultipole(inter_data.insideIdxInBlock);
                                local_expansion_t* cell_local_exp
                                    = &currentCells->getLocalExpansion(inter_data.insideIdxInBlock);

                                FAssertLF(cell_symbolic->getMortonIndex() == inter_data.insideIndex);

                                const int pos = inter_data.relativeOutPosition;
                                kernel->M2L(
                                    cell_local_exp,
                                    cell_symbolic,
                                    &inter_multipole,
                                    &inter_symbolic,
                                    &pos,
                                    1);

                                const int otherPos = getOppositeInterIndex(pos);
                                kernel->M2L(
                                    inter_local_exp,
                                    inter_symbolic,
                                    &cell_multipole,
                                    &cell_symbolic,
                                    &otherPos,
                                    1);
                            }
                        }

                        #pragma omp taskwait

                        ++currentInteractions;
                    }

                    ++iterCells;
                    ++externalInteractionsIter;
                }
            }
            FLOG( timerOutBlock.tac() );
        }
        FLOG( FLog::Controller << "\t\t transferPass in " << timer.tacAndElapsed() << "s\n" );
        FLOG( FLog::Controller << "\t\t\t inblock in  " << timerInBlock.cumulated() << "s\n" );
        FLOG( FLog::Controller << "\t\t\t outblock in " << timerOutBlock.cumulated() << "s\n" );
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

                CellContainerClass* subCellGroups[9];
                memset(subCellGroups, 0, sizeof(CellContainerClass*) * 9);

                // Skip current group if needed
                if( (*iterChildCells)->getEndingIndex() <= (currentCells->getStartingIndex()<<3) ){
                    ++iterChildCells;
                    FAssertLF( iterChildCells != endChildCells );
                    FAssertLF( ((*iterChildCells)->getStartingIndex()>>3) == currentCells->getStartingIndex() );
                }
                // Copy at max 8 groups
                int nbSubCellGroups = 0;
                subCellGroups[nbSubCellGroups] = (*iterChildCells);
                nbSubCellGroups += 1;
                while((*iterChildCells)->getEndingIndex() <= (((currentCells->getEndingIndex()-1)<<3)+7)
                      && (++iterChildCells) != endChildCells
                      && (*iterChildCells)->getStartingIndex() <= ((currentCells->getEndingIndex()-1)<<3)+7 ){
                    FAssertLF( nbSubCellGroups < 9 );
                    subCellGroups[nbSubCellGroups] = (*iterChildCells);
                    nbSubCellGroups += 1;
                }

                #pragma omp task default(none) firstprivate(idxLevel, currentCells, subCellGroups, nbSubCellGroups, kernels)
                {
                    KernelClass*const kernel = kernels[omp_get_thread_num()];
                    int idxSubCellGroup = 0;
                    int idxChildCell = subCellGroups[0]->getFistChildIdx(currentCells->getCellMortonIndex(0));
                    FAssertLF(idxChildCell != -1);

                    for(int cellIdx = 0 ; cellIdx < currentCells->getNumberOfCellsInBlock() ; ++cellIdx){
                        const local_expansion_t* parent_local_exp = &currentCells->getLocalExpansion(cellIdx);
                        const symbolic_data_t* parent_symbolic  = &currentCells->getSymbolic(cellIdx);

                        FAssertLF(parent_symbolic->getMortonIndex()
                                  == currentCells->getCellMortonIndex(cellIdx));

                        std::array<local_expansion_t*, 8> child_local_expansions {};
                        std::array<symbolic_data_t*, 8> child_symbolics {};

                        while(idxSubCellGroup != nbSubCellGroups
                              && (subCellGroups[idxSubCellGroup]->getCellMortonIndex(idxChildCell)>>3) == parent_symbolic->getMortonIndex()) {
                            const int idxChild = ((subCellGroups[idxSubCellGroup]->getCellMortonIndex(idxChildCell)) & 7);
                            FAssertLF(child_symbolics[idxChild] == nullptr);
                            FAssertLF(child_local_expansions[idxChild] == nullptr);

                            child_symbolics[idxChild]
                                = &subCellGroups[idxSubCellGroup]->getSymbolic(idxChildCell);
                            child_local_expansions[idxChild]
                                = &subCellGroups[idxSubCellGroup]->getLocalExpansion(idxChildCell);

                            FAssertLF(subCellGroups[idxSubCellGroup]->getCellMortonIndex(idxChildCell)
                                      == child_symbolics[idxChild]->getMortonIndex());

                            idxChildCell += 1;
                            if(idxChildCell == subCellGroups[idxSubCellGroup]->getNumberOfCellsInBlock()){
                                idxChildCell = 0;
                                idxSubCellGroup += 1;
                            }
                        }

                        kernel->L2L(
                            parent_local_exp,
                            parent_symbolic,
                            child_local_expansions.data(),
                            child_symbolics.data()
                            );
                    }
                }

                ++iterCells;
            }

            #pragma omp taskwait

            FAssertLF(iterCells == endCells && (iterChildCells == endChildCells || (++iterChildCells) == endChildCells));
        }
        FLOG( FLog::Controller << "\t\t downardPass in " << timer.tacAndElapsed() << "s\n" );
    }

    void directPass(){
        FLOG( FTic timer; );
        FLOG( FTic timerInBlock; FTic timerOutBlock; );

        FLOG( timerInBlock.tic() );
        {
            typename OctreeClass::ParticleGroupIterator iterParticles = tree->leavesBegin();
            const typename OctreeClass::ParticleGroupIterator endParticles = tree->leavesEnd();

            while(iterParticles != endParticles){
                ParticleGroupClass* containers = (*iterParticles);

                // FIXME: hack around a clang bug
                // it apparently can't manage a firstprivate const member
                // such as 'tree', but can manage a local copy...
                // getting the height first, or making 'tree' shared both
                // workaround it.
                OctreeClass * const treeAlias = tree;
                #pragma omp task default(none) firstprivate(containers, kernels, treeAlias)
                {
                    const MortonIndex blockStartIdx = containers->getStartingIndex();
                    const MortonIndex blockEndIdx = containers->getEndingIndex();
                    KernelClass*const kernel = kernels[omp_get_thread_num()];

                    for(int leafIdx = 0 ; leafIdx < containers->getNumberOfLeavesInBlock() ; ++leafIdx){
                        ParticleContainerClass particles = containers->template getLeaf<ParticleContainerClass>(leafIdx);
                        const MortonIndex mindex = containers->getLeafMortonIndex(leafIdx);

                        MortonIndex interactionsIndexes[26];
                        int interactionsPosition[26];
                        FTreeCoordinate coord(mindex);
                        int counter = coord.getNeighborsIndexes(treeAlias->getHeight(),interactionsIndexes,interactionsPosition);

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
            #pragma omp taskwait
        }
        FLOG( timerInBlock.tac() );
        FLOG( timerOutBlock.tic() );
        {
            typename OctreeClass::ParticleGroupIterator iterParticles = tree->leavesBegin();
            const typename OctreeClass::ParticleGroupIterator endParticles = tree->leavesEnd();

            typename std::vector<std::vector<BlockInteractions<ParticleGroupClass>>>::iterator externalInteractionsIter = externalInteractionsLeafLevel.begin();

            while(iterParticles != endParticles){
                typename std::vector<BlockInteractions<ParticleGroupClass>>::iterator currentInteractions = (*externalInteractionsIter).begin();
                const typename std::vector<BlockInteractions<ParticleGroupClass>>::iterator currentInteractionsEnd = (*externalInteractionsIter).end();

                ParticleGroupClass* containers = (*iterParticles);

                while(currentInteractions != currentInteractionsEnd){
                    ParticleGroupClass* containersOther = (*currentInteractions).otherBlock;
                    const std::vector<OutOfBlockInteraction>* outsideInteractions = &(*currentInteractions).interactions;

                    #pragma omp task default(none) firstprivate(containers, containersOther, outsideInteractions, kernels)
                    {
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
                    // only one task but need to wait for it
                    #pragma omp taskwait

                    ++currentInteractions;
                }

                ++iterParticles;
                ++externalInteractionsIter;
            }
        }
        FLOG( timerOutBlock.tac() );

        FLOG( FLog::Controller << "\t\t directPass in " << timer.tacAndElapsed() << "s\n" );
        FLOG( FLog::Controller << "\t\t\t inblock  in " << timerInBlock.cumulated() << "s\n" );
        FLOG( FLog::Controller << "\t\t\t outblock in " << timerOutBlock.cumulated() << "s\n" );
    }

    void mergePass(){
        FLOG( FTic timer; );

        for(int idxGroup = 0 ; idxGroup < tree->getNbParticleGroup() ; ++idxGroup){
            CellContainerClass* leafCells  = tree->getCellGroup(tree->getHeight()-1, idxGroup);
            ParticleGroupClass* containers = tree->getParticleGroup(idxGroup);
            #pragma omp task default(shared) firstprivate(leafCells, containers, kernels)
            {
                KernelClass*const kernel = kernels[omp_get_thread_num()];

                for(int cellIdx = 0 ; cellIdx < leafCells->getNumberOfCellsInBlock() ; ++cellIdx){
                    const local_expansion_t* leaf_local_exp = &leafCells->getLocalExpansion(cellIdx);
                    const symbolic_data_t* leaf_symbolic  = &leafCells->getSymbolic(cellIdx);
                    ParticleContainerClass particles
                        = containers->template getLeaf<ParticleContainerClass>(cellIdx);

                    FAssertLF(leaf_symbolic->getMortonIndex() == leafCells->getCellMortonIndex(cellIdx));
                    FAssertLF(leafCells->getCellMortonIndex(cellIdx) == containers->getLeafMortonIndex(cellIdx));

                    kernel->L2P(leaf_local_exp, leaf_symbolic, &particles);
                }
            }
        }
        // Wait for task to complete
        #pragma omp taskwait

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

#endif // FGROUPTASKALGORITHM_HPP
