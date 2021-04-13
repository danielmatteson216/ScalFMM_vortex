
// Keep in private GIT
#ifndef FGROUPSEQALGORITHM_HPP
#define FGROUPSEQALGORITHM_HPP

#include "Utils/FGlobal.hpp"
#include "Core/FCoreCommon.hpp"
#include "Utils/FQuickSort.hpp"
#include "Containers/FTreeCoordinate.hpp"
#include "Utils/FLog.hpp"
#include "Utils/FTic.hpp"

#include "FOutOfBlockInteraction.hpp"

#include <vector>
#include <vector>

template <class OctreeClass, class CellContainerClass, class CellClass, class KernelClass, class ParticleGroupClass, class ParticleContainerClass>
class FGroupSeqAlgorithm : public FAbstractAlgorithm {
protected:
    const int MaxThreads;         //< The number of threads
    OctreeClass*const tree;       //< The Tree
    KernelClass*const kernels;    //< The kernels

public:
    using multipole_t = typename std::remove_pointer<typename OctreeClass::CellGroupIterator::value_type>::type::multipole_t;
    using local_expansion_t = typename std::remove_pointer<typename OctreeClass::CellGroupIterator::value_type>::type::local_expansion_t;
    using symbolic_data_t = typename std::remove_pointer<typename OctreeClass::CellGroupIterator::value_type>::type::symbolic_data_t;

    FGroupSeqAlgorithm(OctreeClass*const inTree, KernelClass* inKernels) : MaxThreads(1), tree(inTree), kernels(inKernels){
        FAssertLF(tree, "tree cannot be null");
        FAssertLF(kernels, "kernels cannot be null");

        FAbstractAlgorithm::setNbLevelsInTree(tree->getHeight());

        FLOG(FLog::Controller << "FGroupSeqAlgorithm (Max Thread " << MaxThreads << ")\n");
    }

    ~FGroupSeqAlgorithm(){
    }

protected:
    /**
      * Runs the complete algorithm.
      */
    void executeCore(const unsigned operationsToProceed) override {
        FLOG( FLog::Controller << "\tStart FGroupSeqAlgorithm\n" );

        if(operationsToProceed & FFmmP2M) bottomPass();

        if(operationsToProceed & FFmmM2M) upwardPass();

        if(operationsToProceed & FFmmM2L) transferPass();

        if(operationsToProceed & FFmmL2L) downardPass();

        if( (operationsToProceed & FFmmP2P) || (operationsToProceed & FFmmL2P) ){
            directPass((operationsToProceed & FFmmP2P), (operationsToProceed & FFmmL2P));
        }
    }

    void bottomPass(){
        FLOG( FTic timer; );
        typename OctreeClass::ParticleGroupIterator iterParticles = tree->leavesBegin();
        const typename OctreeClass::ParticleGroupIterator endParticles = tree->leavesEnd();

        typename OctreeClass::CellGroupIterator iterCells = tree->cellsBegin(tree->getHeight()-1);
        const typename OctreeClass::CellGroupIterator endCells = tree->cellsEnd(tree->getHeight()-1);

        while(iterParticles != endParticles && iterCells != endCells){
            { // Can be a task(in:iterParticles, out:iterCells)
                FAssertLF((*iterCells)->getNumberOfCellsInBlock() == (*iterParticles)->getNumberOfLeavesInBlock());

                for(int leafIdx = 0 ; leafIdx < (*iterCells)->getNumberOfCellsInBlock() ; ++leafIdx){
                    multipole_t* leaf_multipole = &(*iterCells)->getMultipole(leafIdx);
                    symbolic_data_t* leaf_symbolic = &(*iterCells)->getSymbolic(leafIdx);

                    ParticleContainerClass particles = (*iterParticles)->template getLeaf<ParticleContainerClass>(leafIdx);
                    FAssertLF((*iterCells)->getCellMortonIndex(leafIdx) == (*iterParticles)->getLeafMortonIndex(leafIdx));
                    kernels->P2M(leaf_multipole, leaf_symbolic, &particles);
                }
            }

            ++iterParticles;
            ++iterCells;
        }

        FAssertLF(iterParticles == endParticles && iterCells == endCells);
        FLOG( FLog::Controller << "\t\t bottomPass in " << timer.tacAndElapsed() << "s\n" );
    }

    void upwardPass(){
        FLOG( FTic timer; );
        for(int idxLevel = FMath::Min(tree->getHeight() - 2, FAbstractAlgorithm::lowerWorkingLevel - 1) ; idxLevel >= FAbstractAlgorithm::upperWorkingLevel ; --idxLevel){
            typename OctreeClass::CellGroupIterator iterCells = tree->cellsBegin(idxLevel);
            const typename OctreeClass::CellGroupIterator endCells = tree->cellsEnd(idxLevel);

            typename OctreeClass::CellGroupIterator iterChildCells = tree->cellsBegin(idxLevel+1);
            const typename OctreeClass::CellGroupIterator endChildCells = tree->cellsEnd(idxLevel+1);

            int idxChildCell = 0;

            while(iterCells != endCells && iterChildCells != endChildCells){
                { // Can be a task(in:iterParticles, out:iterChildCells ...)

                    for(int cellIdx = 0 ; cellIdx < (*iterCells)->getNumberOfCellsInBlock() ; ++cellIdx){
                        multipole_t* parent_multipole = &(*iterCells)->getMultipole(cellIdx);
                        symbolic_data_t* parent_symbolic  = &(*iterCells)->getSymbolic(cellIdx);

                        std::array<const multipole_t*,8> child_multipoles {};
                        std::array<const symbolic_data_t*,8> child_symbolics {};

                        FAssertLF(iterChildCells != endChildCells);

                        while(iterChildCells != endChildCells
                              && (((*iterChildCells)->getCellMortonIndex(idxChildCell)>>3)
                                  == parent_symbolic->getMortonIndex())) {
                            const int idxChild = (((*iterChildCells)->getCellMortonIndex(idxChildCell)) & 7);

                            FAssertLF(child_multipoles[idxChild] == nullptr);
                            FAssertLF(child_symbolics [idxChild] == nullptr);

                            child_multipoles[idxChild] = &(*iterChildCells)->getMultipole(idxChildCell);
                            child_symbolics [idxChild] = &(*iterChildCells)->getSymbolic(idxChildCell);

                            FAssertLF((*iterChildCells)->getCellMortonIndex(idxChildCell)
                                      == child_symbolics[idxChild]->getMortonIndex());

                            idxChildCell += 1;
                            if(idxChildCell == (*iterChildCells)->getNumberOfCellsInBlock()){
                                idxChildCell = 0;
                                ++iterChildCells;
                            }
                        }

                        kernels->M2M(
                            parent_multipole,
                            parent_symbolic,
                            child_multipoles.data(),
                            child_symbolics.data()
                            );
                    }
                }

                ++iterCells;
            }

            FAssertLF(iterCells == endCells);
            FAssertLF(iterChildCells == endChildCells);
            FAssertLF(iterCells == endCells && (iterChildCells == endChildCells || (++iterChildCells) == endChildCells));
        }
        FLOG( FLog::Controller << "\t\t upwardPass in " << timer.tacAndElapsed() << "s\n" );
    }

    void transferPass(){
        FLOG( FTic timer; );
        for(int idxLevel = FAbstractAlgorithm::lowerWorkingLevel-1 ; idxLevel >= FAbstractAlgorithm::upperWorkingLevel ; --idxLevel){
            typename OctreeClass::CellGroupIterator iterCells = tree->cellsBegin(idxLevel);
            const typename OctreeClass::CellGroupIterator endCells = tree->cellsEnd(idxLevel);

            while(iterCells != endCells){
                std::vector<OutOfBlockInteraction> outsideInteractions;

                { // Can be a task(inout:iterCells, out:outsideInteractions)
                    std::array<const multipole_t*, 189> neighbour_multipoles {};
                    std::array<const symbolic_data_t*, 189> neighbour_symbolics {};

                    const MortonIndex blockStartIdx = (*iterCells)->getStartingIndex();
                    const MortonIndex blockEndIdx = (*iterCells)->getEndingIndex();

                    for(int cellIdx  = 0 ; cellIdx < (*iterCells)->getNumberOfCellsInBlock() ; ++cellIdx){
                        local_expansion_t* target_local_exp = &(*iterCells)->getLocalExpansion(cellIdx);
                        symbolic_data_t* target_symbolic  = &(*iterCells)->getSymbolic(cellIdx);

                        const MortonIndex mindex = (*iterCells)->getCellMortonIndex(cellIdx);
                        FAssertLF(target_symbolic->getMortonIndex() == mindex);

                        MortonIndex interactionsIndexes[189];
                        int interactionsPosition[189];
                        const FTreeCoordinate coord(target_symbolic->getCoordinate());
                        int counter = coord.getInteractionNeighbors(idxLevel,interactionsIndexes,interactionsPosition);

                        int counterExistingCell = 0;

                        for(int idxInter = 0 ; idxInter < counter ; ++idxInter){
                            if( blockStartIdx <= interactionsIndexes[idxInter] && interactionsIndexes[idxInter] < blockEndIdx ){
                                const int cellPos = (*iterCells)->getCellIndex(interactionsIndexes[idxInter]);
                                if(cellPos != -1){
                                    const multipole_t* neighbour_multipole
                                        = &(*iterCells)->getMultipole(cellPos);
                                    const symbolic_data_t* neighbour_symbolic
                                        = &(*iterCells)->getSymbolic(cellPos);
                                    FAssertLF(neighbour_symbolic->getMortonIndex()
                                              == interactionsIndexes[idxInter]);

                                    interactionsPosition[counterExistingCell] = interactionsPosition[idxInter];
                                    neighbour_multipoles[counterExistingCell] = neighbour_multipole;
                                    neighbour_symbolics [counterExistingCell] = neighbour_symbolic;
                                    counterExistingCell += 1;
                                }
                            }
                            else if(interactionsIndexes[idxInter] < mindex){
                                OutOfBlockInteraction property;
                                property.insideIndex = mindex;
                                property.outIndex    = interactionsIndexes[idxInter];
                                property.relativeOutPosition = interactionsPosition[idxInter];
                                property.insideIdxInBlock = cellIdx;
                                outsideInteractions.push_back(property);
                            }
                        }

                        kernels->M2L(
                            target_local_exp,
                            target_symbolic,
                            neighbour_multipoles.data(),
                            neighbour_symbolics.data(),
                            interactionsPosition,
                            counterExistingCell
                            );
                    }
                }


                // Manage outofblock interaction
                FQuickSort<OutOfBlockInteraction, int>::QsSequential(outsideInteractions.data(),int(outsideInteractions.size()));

                typename OctreeClass::CellGroupIterator iterLeftCells = tree->cellsBegin(idxLevel);
                int currentOutInteraction = 0;
                while(iterLeftCells != iterCells && currentOutInteraction < int(outsideInteractions.size())){
                    const MortonIndex outBlockStartIdx = (*iterLeftCells)->getStartingIndex();
                    const MortonIndex outBlockEndIdx = (*iterLeftCells)->getEndingIndex();

                    while(currentOutInteraction < int(outsideInteractions.size()) && outsideInteractions[currentOutInteraction].outIndex < outBlockStartIdx){
                        currentOutInteraction += 1;
                    }

                    int lastOutInteraction = currentOutInteraction;
                    while(lastOutInteraction < int(outsideInteractions.size()) && outsideInteractions[lastOutInteraction].outIndex < outBlockEndIdx){
                        lastOutInteraction += 1;
                    }

                    { // Can be a task(in:currentOutInteraction, in:outsideInteractions, in:lastOutInteraction, inout:iterLeftCells, inout:iterCells)

                        for(int outInterIdx = currentOutInteraction ; outInterIdx < lastOutInteraction ; ++outInterIdx){
                            const int cellPos = (*iterLeftCells)->getCellIndex(outsideInteractions[outInterIdx].outIndex);
                            if(cellPos != -1){
                                const symbolic_data_t* inter_symbolic
                                    = &(*iterLeftCells)->getSymbolic(cellPos);
                                const multipole_t* inter_multipole
                                    = &(*iterLeftCells)->getMultipole(cellPos);
                                local_expansion_t* inter_local_exp
                                    = &(*iterLeftCells)->getLocalExpansion(cellPos);

                                FAssertLF(inter_symbolic->getMortonIndex()
                                          == outsideInteractions[outInterIdx].outIndex);
                                const symbolic_data_t* cell_symbolic
                                    = &(*iterCells)->getSymbolic(outsideInteractions[outInterIdx].insideIdxInBlock);
                                const multipole_t* cell_multipole
                                    = &(*iterCells)->getMultipole(outsideInteractions[outInterIdx].insideIdxInBlock);
                                local_expansion_t* cell_local_exp
                                    = &(*iterCells)->getLocalExpansion(outsideInteractions[outInterIdx].insideIdxInBlock);

                                FAssertLF(cell_symbolic->getMortonIndex()
                                          == outsideInteractions[outInterIdx].insideIndex);

                                kernels->M2L(
                                    cell_local_exp,
                                    cell_symbolic,
                                    &inter_multipole,
                                    &inter_symbolic,
                                    &outsideInteractions[outInterIdx].relativeOutPosition,
                                    1);
                                const int otherPos = getOppositeInterIndex(outsideInteractions[outInterIdx].relativeOutPosition);
                                kernels->M2L(
                                    inter_local_exp,
                                    inter_symbolic,
                                    &cell_multipole,
                                    &cell_symbolic,
                                    &otherPos,
                                    1);
                            }
                        }
                    }

                    currentOutInteraction = lastOutInteraction;
                    ++iterLeftCells;
                }

                ++iterCells;
            }

        }
        FLOG( FLog::Controller << "\t\t transferPass in " << timer.tacAndElapsed() << "s\n" );
    }

    void downardPass(){
        FLOG( FTic timer; );
        for(int idxLevel = FAbstractAlgorithm::upperWorkingLevel ; idxLevel < FAbstractAlgorithm::lowerWorkingLevel - 1 ; ++idxLevel){
            typename OctreeClass::CellGroupIterator iterCells = tree->cellsBegin(idxLevel);
            const typename OctreeClass::CellGroupIterator endCells = tree->cellsEnd(idxLevel);

            typename OctreeClass::CellGroupIterator iterChildCells = tree->cellsBegin(idxLevel+1);
            const typename OctreeClass::CellGroupIterator endChildCells = tree->cellsEnd(idxLevel+1);

            int idxChildCell = 0;

            while(iterCells != endCells && iterChildCells != endChildCells){
                { // Can be a task(in:iterParticles, inout:iterChildCells ...)

                    for(int cellIdx = 0 ; cellIdx < (*iterCells)->getNumberOfCellsInBlock() ; ++cellIdx){
                        const local_expansion_t* parent_local_exp = &(*iterCells)->getLocalExpansion(cellIdx);
                        const symbolic_data_t* parent_symbolic  = &(*iterCells)->getSymbolic(cellIdx);
                        FAssertLF(parent_symbolic->getMortonIndex()
                                  == (*iterCells)->getCellMortonIndex(cellIdx));
                        std::array<local_expansion_t*, 8> child_local_expansions {};
                        std::array<const symbolic_data_t*, 8> child_symbolics {};

                        FAssertLF(iterChildCells != endChildCells);

                        while(iterChildCells != endChildCells
                              && (((*iterChildCells)->getCellMortonIndex(idxChildCell)>>3)
                                  == parent_symbolic->getMortonIndex())){
                            const int idxChild = (((*iterChildCells)->getCellMortonIndex(idxChildCell)) & 7);

                            FAssertLF(child_symbolics[idxChild] == nullptr);
                            FAssertLF(child_local_expansions[idxChild] == nullptr);

                            child_symbolics[idxChild]
                                = &(*iterChildCells)->getSymbolic(idxChildCell);
                            child_local_expansions[idxChild]
                                = &(*iterChildCells)->getLocalExpansion(idxChildCell);

                            FAssertLF((*iterChildCells)->getCellMortonIndex(idxChildCell)
                                      == child_symbolics[idxChild]->getMortonIndex());

                            idxChildCell += 1;
                            if(idxChildCell == (*iterChildCells)->getNumberOfCellsInBlock()){
                                idxChildCell = 0;
                                ++iterChildCells;
                            }
                        }

                        kernels->L2L(
                            parent_local_exp,
                            parent_symbolic,
                            child_local_expansions.data(),
                            child_symbolics.data()
                            );
                    }
                }

                ++iterCells;
            }

            FAssertLF(iterCells == endCells && iterChildCells == endChildCells);
        }
        FLOG( FLog::Controller << "\t\t downardPass in " << timer.tacAndElapsed() << "s\n" );
    }

    void directPass(const bool p2pEnabled, const bool l2pEnabled){
        FLOG( FTic timer; );
        if(l2pEnabled){
            typename OctreeClass::ParticleGroupIterator iterParticles = tree->leavesBegin();
            const typename OctreeClass::ParticleGroupIterator endParticles = tree->leavesEnd();

            typename OctreeClass::CellGroupIterator iterCells = tree->cellsBegin(tree->getHeight()-1);
            const typename OctreeClass::CellGroupIterator endCells = tree->cellsEnd(tree->getHeight()-1);

            while(iterParticles != endParticles && iterCells != endCells){
                { // Can be a task(in:iterCells, inout:iterParticles)
                    for(int leafIdx = 0 ; leafIdx < (*iterCells)->getNumberOfCellsInBlock() ; ++leafIdx){
                        const local_expansion_t* leaf_local_exp = &(*iterCells)->getLocalExpansion(leafIdx);
                        const symbolic_data_t* leaf_symbolic  = &(*iterCells)->getSymbolic(leafIdx);

                        ParticleContainerClass particles = (*iterParticles)->template getLeaf<ParticleContainerClass>(leafIdx);
                        FAssertLF((*iterCells)->getCellMortonIndex(leafIdx) == (*iterParticles)->getLeafMortonIndex(leafIdx));
                        kernels->L2P(leaf_local_exp, leaf_symbolic, &particles);
                    }
                }

                ++iterParticles;
                ++iterCells;
            }

            FAssertLF(iterParticles == endParticles && iterCells == endCells);
        }
        if(p2pEnabled){
            typename OctreeClass::ParticleGroupIterator iterParticles = tree->leavesBegin();
            const typename OctreeClass::ParticleGroupIterator endParticles = tree->leavesEnd();

            while(iterParticles != endParticles){
                typename std::vector<OutOfBlockInteraction> outsideInteractions;

                { // Can be a task(inout:iterCells, out:outsideInteractions)
                    const MortonIndex blockStartIdx = (*iterParticles)->getStartingIndex();
                    const MortonIndex blockEndIdx = (*iterParticles)->getEndingIndex();

                    for(int leafIdx = 0 ; leafIdx < (*iterParticles)->getNumberOfLeavesInBlock() ; ++leafIdx){
                        ParticleContainerClass particles = (*iterParticles)->template getLeaf<ParticleContainerClass>(leafIdx);

                        const MortonIndex mindex = (*iterParticles)->getLeafMortonIndex(leafIdx);
                        MortonIndex interactionsIndexes[26];
                        int interactionsPosition[26];
                        FTreeCoordinate coord(mindex);
                        int counter = coord.getNeighborsIndexes(tree->getHeight(),interactionsIndexes,interactionsPosition);

                        ParticleContainerClass interactionsObjects[26];
                        ParticleContainerClass* interactions[26];
                        int counterExistingCell = 0;

                        for(int idxInter = 0 ; idxInter < counter ; ++idxInter){
                            if( blockStartIdx <= interactionsIndexes[idxInter] && interactionsIndexes[idxInter] < blockEndIdx ){
                                const int leafPos = (*iterParticles)->getLeafIndex(interactionsIndexes[idxInter]);
                                if(leafPos != -1){
                                    interactionsObjects[counterExistingCell] = (*iterParticles)->template getLeaf<ParticleContainerClass>(leafPos);
                                    interactionsPosition[counterExistingCell] = interactionsPosition[idxInter];
                                    interactions[counterExistingCell] = &interactionsObjects[counterExistingCell];
                                    counterExistingCell += 1;
                                }
                            }
                            else if(interactionsIndexes[idxInter] < mindex){
                                OutOfBlockInteraction property;
                                property.insideIndex = mindex;
                                property.outIndex    = interactionsIndexes[idxInter];
                                property.relativeOutPosition = interactionsPosition[idxInter];
                                property.insideIdxInBlock = leafIdx;
                                outsideInteractions.push_back(property);
                            }
                        }

                        kernels->P2P( coord, &particles, &particles , interactions, interactionsPosition, counterExistingCell);
                    }
                }


                // Manage outofblock interaction
                FQuickSort<OutOfBlockInteraction, int>::QsSequential(outsideInteractions.data(),int(outsideInteractions.size()));

                typename OctreeClass::ParticleGroupIterator iterLeftParticles = tree->leavesBegin();
                int currentOutInteraction = 0;
                while(iterLeftParticles != iterParticles && currentOutInteraction < int(outsideInteractions.size())){
                    const MortonIndex blockStartIdx = (*iterLeftParticles)->getStartingIndex();
                    const MortonIndex blockEndIdx = (*iterLeftParticles)->getEndingIndex();

                    while(currentOutInteraction < int(outsideInteractions.size()) && outsideInteractions[currentOutInteraction].outIndex < blockStartIdx){
                        currentOutInteraction += 1;
                    }

                    int lastOutInteraction = currentOutInteraction;
                    while(lastOutInteraction < int(outsideInteractions.size()) && outsideInteractions[lastOutInteraction].outIndex < blockEndIdx){
                        lastOutInteraction += 1;
                    }

                    { // Can be a task(in:currentOutInteraction, in:outsideInteractions, in:lastOutInteraction, inout:iterLeftParticles, inout:iterParticles)
                        for(int outInterIdx = currentOutInteraction ; outInterIdx < lastOutInteraction ; ++outInterIdx){
                            const int leafPos = (*iterLeftParticles)->getLeafIndex(outsideInteractions[outInterIdx].outIndex);
                            if(leafPos != -1){
                                ParticleContainerClass interParticles = (*iterLeftParticles)->template getLeaf<ParticleContainerClass>(leafPos);
                                ParticleContainerClass particles = (*iterParticles)->template getLeaf<ParticleContainerClass>(outsideInteractions[outInterIdx].insideIdxInBlock);

                                FAssertLF((*iterLeftParticles)->getLeafMortonIndex(leafPos) == outsideInteractions[outInterIdx].outIndex);
                                FAssertLF((*iterParticles)->getLeafMortonIndex(outsideInteractions[outInterIdx].insideIdxInBlock) == outsideInteractions[outInterIdx].insideIndex);

                                ParticleContainerClass* ptrLeaf = &interParticles;
                                kernels->P2POuter( FTreeCoordinate(outsideInteractions[outInterIdx].insideIndex),
                                                    &particles , &ptrLeaf, &outsideInteractions[outInterIdx].relativeOutPosition, 1);
                                const int otherPosition = getOppositeNeighIndex(outsideInteractions[outInterIdx].relativeOutPosition);
                                ptrLeaf = &particles;
                                kernels->P2POuter( FTreeCoordinate(outsideInteractions[outInterIdx].outIndex),
                                                    &interParticles , &ptrLeaf, &otherPosition, 1);
                            }
                        }
                    }

                    currentOutInteraction = lastOutInteraction;
                    ++iterLeftParticles;
                }

                ++iterParticles;
            }
        }
        FLOG( FLog::Controller << "\t\t directPass in " << timer.tacAndElapsed() << "s\n" );
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


#endif // FGROUPSEQALGORITHM_HPP
