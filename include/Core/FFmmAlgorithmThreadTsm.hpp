// See LICENCE file at project root
#ifndef FFMMALGORITHMTHREADTSM_HPP
#define FFMMALGORITHMTHREADTSM_HPP

#include <array>
#include <algorithm>

#include "../Utils/FAssert.hpp"
#include "../Utils/FLog.hpp"

#include "../Utils/FTic.hpp"
#include "../Utils/FGlobal.hpp"
#include "../Utils/FAlgorithmTimers.hpp"
#include "../Containers/FOctree.hpp"
#include "FCoreCommon.hpp"
#include "../Utils/FEnv.hpp"

#include <omp.h>

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FFmmAlgorithmThreadTsm
* @brief
* Please read the license
*
* This class is a threaded FMM algorithm
* It just iterates on a tree and call the kernels with good arguments.
* It used the inspector-executor model :
* iterates on the tree and builds an array to work in parallel on this array
*
* Of course this class does not deallocate pointer given in arguements.
*
* Because this is a Target source model you do not need the P2P to be safe.
* You should not write on sources in the P2P method!
*/
template<class OctreeClass, class CellClass, class ContainerClass, class KernelClass, class LeafClass>
class FFmmAlgorithmThreadTsm : public FAbstractAlgorithm, public FAlgorithmTimers{

    using multipole_t = typename CellClass::multipole_t;
    using local_expansion_t = typename CellClass::local_expansion_t;
    using symbolic_data_t = CellClass;

    OctreeClass* const tree;                  //< The octree to work on
    KernelClass** kernels;                    //< The kernels

    typename OctreeClass::Iterator* iterArray;

    int MaxThreads;

    const int OctreeHeight;

    const int leafLevelSeparationCriteria;

public:
    /** The constructor need the octree and the kernels used for computation
      * @param inTree the octree to work on
      * @param inKernels the kernels to call
      * An assert is launched if one of the arguments is null
      */
    FFmmAlgorithmThreadTsm(OctreeClass* const inTree, const KernelClass* const inKernels, const int inLeafLevelSeperationCriteria = 1)
                      : tree(inTree) , kernels(nullptr), iterArray(nullptr),
                        OctreeHeight(tree->getHeight()), leafLevelSeparationCriteria(inLeafLevelSeperationCriteria) {

        FAssertLF(tree, "tree cannot be null");
        FAssertLF(leafLevelSeparationCriteria < 3, "Separation criteria should be < 3");

        MaxThreads = 1;
        #pragma omp parallel
        #pragma omp master
            MaxThreads = omp_get_num_threads();

        this->kernels = new KernelClass*[MaxThreads];
        #pragma omp parallel num_threads(MaxThreads) 
        {
          #pragma omp critical (InitFFmmAlgorithmTsm)
          {
            this->kernels[omp_get_thread_num()] = new KernelClass(*inKernels);
          }
        }

        FAbstractAlgorithm::setNbLevelsInTree(tree->getHeight());

        FLOG(FLog::Controller << "FFmmAlgorithmThreadTsm\n");
    }

    /** Default destructor */
    virtual ~FFmmAlgorithmThreadTsm(){
        for(int idxThread = 0 ; idxThread < MaxThreads ; ++idxThread){
            delete this->kernels[idxThread];
        }
        delete [] this->kernels;
    }

protected:
    /**
      * To execute the fmm algorithm
      * Call this function to run the complete algorithm
      */
    void executeCore(const unsigned operationsToProceed) override {
        // Count leaf
        int numberOfLeafs = 0;
        typename OctreeClass::Iterator octreeIterator(tree);
        octreeIterator.gotoBottomLeft();
        do{
            ++numberOfLeafs;
        } while(octreeIterator.moveRight());
        iterArray = new typename OctreeClass::Iterator[numberOfLeafs];
        FAssertLF(iterArray, "iterArray bad alloc");

        if(operationsToProceed & FFmmP2M) bottomPass();

        if(operationsToProceed & FFmmM2M) upwardPass();

        if(operationsToProceed & FFmmM2L) transferPass();

        if(operationsToProceed & FFmmL2L) downardPass();

        if((operationsToProceed & FFmmP2P) || (operationsToProceed & FFmmL2P)) directPass((operationsToProceed & FFmmP2P),(operationsToProceed & FFmmL2P));

        delete [] iterArray;
        iterArray = nullptr;


    }

    /** P2M */
    void bottomPass(){
        FLOG( FLog::Controller.write("\tStart Bottom Pass\n").write(FLog::Flush) );
        FLOG( FTic counterTime );

        typename OctreeClass::Iterator octreeIterator(tree);
        int numberOfLeafs = 0;
        // Iterate on leafs
        octreeIterator.gotoBottomLeft();
        do{
            iterArray[numberOfLeafs] = octreeIterator;
            ++numberOfLeafs;
        } while(octreeIterator.moveRight());

        const int chunkSize = FMath::Max(1 , numberOfLeafs/(omp_get_num_threads()*omp_get_num_threads()));

        FLOG(FTic computationCounter);
        #pragma omp parallel num_threads(MaxThreads)
        {
            KernelClass * const myThreadkernels = kernels[omp_get_thread_num()];
            #pragma omp for nowait schedule(dynamic, chunkSize)
            for(int idxLeafs = 0 ; idxLeafs < numberOfLeafs ; ++idxLeafs){
                // We need the current cell that represent the leaf
                // and the list of particles
                ContainerClass* const sources = iterArray[idxLeafs].getCurrentListSrc();
                if(sources->getNbParticles()){
                    iterArray[idxLeafs].getCurrentCell()->setSrcChildTrue();
                    myThreadkernels->P2M(
                        &(iterArray[idxLeafs].getCurrentCell()->getMultipoleData()),
                        iterArray[idxLeafs].getCurrentCell(),
                        sources);
                }
                if(iterArray[idxLeafs].getCurrentListTargets()->getNbParticles()){
                    iterArray[idxLeafs].getCurrentCell()->setTargetsChildTrue();
                }
            }
        }
        FLOG(computationCounter.tac());

        FLOG( counterTime.tac() );
        FLOG( FLog::Controller << "\tFinished (@Bottom Pass (P2M) = "  << counterTime.elapsed() << " s)\n" );
        FLOG( FLog::Controller << "\t\t Computation : " << computationCounter.elapsed() << " s\n" );

    }

    /** M2M */
    void upwardPass(){
        FLOG( FLog::Controller.write("\tStart Upward Pass\n").write(FLog::Flush); );
        FLOG(FTic counterTime);
        FLOG(FTic computationCounter);

        // Start from leal level - 1
        typename OctreeClass::Iterator octreeIterator(tree);
        octreeIterator.gotoBottomLeft();
        octreeIterator.moveUp();

        for(int idxLevel = OctreeHeight - 2 ; idxLevel > FAbstractAlgorithm::lowerWorkingLevel-1 ; --idxLevel){
            octreeIterator.moveUp();
        }

        typename OctreeClass::Iterator avoidGotoLeftIterator(octreeIterator);

        // for each levels
        for(int idxLevel = FMath::Min(OctreeHeight - 2, FAbstractAlgorithm::lowerWorkingLevel - 1) ; idxLevel >= FAbstractAlgorithm::upperWorkingLevel ; --idxLevel ){
            FLOG(FTic counterTimeLevel);
            int numberOfCells = 0;
            // for each cells
            do{
                iterArray[numberOfCells] = octreeIterator;
                ++numberOfCells;
            } while(octreeIterator.moveRight());
            avoidGotoLeftIterator.moveUp();
            octreeIterator = avoidGotoLeftIterator;// equal octreeIterator.moveUp(); octreeIterator.gotoLeft();

            const int chunkSize = FMath::Max(1 , numberOfCells/(omp_get_num_threads()*omp_get_num_threads()));

            FLOG(computationCounter.tic());
            #pragma omp parallel num_threads(MaxThreads)
            {
                KernelClass * const myThreadkernels = kernels[omp_get_thread_num()];
                #pragma omp for nowait schedule(dynamic, chunkSize)
                for(int idxCell = 0 ; idxCell < numberOfCells ; ++idxCell){
                    // We need the current cell and the child
                    // child is an array (of 8 child) that may be null
                    CellClass* source_children[8] = {};
                    CellClass** const children = iterArray[idxCell].getCurrentChild();
                    CellClass* const current_cell = iterArray[idxCell].getCurrentCell();
                    int nbChildWithSrc = 0;
                    for(int idxChild = 0 ; idxChild < 8 ; ++idxChild){
                        if(children[idxChild]){
                            if(children[idxChild]->hasSrcChild()){
                                nbChildWithSrc += 1;
                                source_children[idxChild] = children[idxChild];
                            }
                            if(children[idxChild]->hasTargetsChild()){
                                current_cell->setTargetsChildTrue();
                            }
                        }
                    }
                    if(nbChildWithSrc){
                        current_cell->setSrcChildTrue();

                        multipole_t* const parent_multipole
                            = &(current_cell->getMultipoleData());
                        const symbolic_data_t* const parent_symbolic = current_cell;

                        std::array<const multipole_t*, 8> child_multipoles;
                        std::transform(source_children, source_children+8,
                                       child_multipoles.begin(),
                                       [](CellClass* c) {
                                           return (c == nullptr ? nullptr
                                                   : &(c->getMultipoleData()));
                                       });
                        std::array<const symbolic_data_t*, 8> child_symbolics;
                        std::transform(source_children, source_children+8,
                                       child_symbolics.begin(),
                                       [](CellClass* c) {return c;});

                        myThreadkernels->M2M(parent_multipole,
                                             parent_symbolic,
                                             child_multipoles.data(),
                                             child_symbolics.data());                    }
                }
            }
            FLOG(computationCounter.tac());
            FLOG( FLog::Controller << "\t\t>> Level " << idxLevel << " = "  << counterTimeLevel.tacAndElapsed() << " s\n" );
        }

        FLOG( counterTime.tac() );
        FLOG( FLog::Controller << "\tFinished (@Upward Pass (M2M) = "  << counterTime.elapsed() << " s)\n" );
        FLOG( FLog::Controller << "\t\t Computation : " << computationCounter.cumulated() << " s\n" );

    }

    /** M2L */
    void transferPass(){
            FLOG( FLog::Controller.write("\tStart Downward Pass (M2L)\n").write(FLog::Flush); );
            FLOG(FTic counterTime);
            FLOG(FTic computationCounter);

            typename OctreeClass::Iterator octreeIterator(tree);
            octreeIterator.moveDown();

            for(int idxLevel = 2 ; idxLevel < FAbstractAlgorithm::upperWorkingLevel ; ++idxLevel){
                octreeIterator.moveDown();
            }

            typename OctreeClass::Iterator avoidGotoLeftIterator(octreeIterator);

            // for each levels
            for(int idxLevel = FAbstractAlgorithm::upperWorkingLevel ; idxLevel < FAbstractAlgorithm::lowerWorkingLevel ; ++idxLevel ){
                FLOG(FTic counterTimeLevel);
                const int separationCriteria = (idxLevel != FAbstractAlgorithm::lowerWorkingLevel-1 ? 1 : leafLevelSeparationCriteria);

                int numberOfCells = 0;
                // for each cells
                do{
                    iterArray[numberOfCells] = octreeIterator;
                    ++numberOfCells;
                } while(octreeIterator.moveRight());
                avoidGotoLeftIterator.moveDown();
                octreeIterator = avoidGotoLeftIterator;

                const int chunkSize = FMath::Max(1 , numberOfCells/(omp_get_num_threads()*omp_get_num_threads()));

                FLOG(computationCounter.tic());
                #pragma omp parallel num_threads(MaxThreads)
                {
                    KernelClass * const myThreadkernels = kernels[omp_get_thread_num()];
                    const CellClass* neighbors[342];
                    int neighborPositions[342];

                    #pragma omp for  schedule(dynamic, chunkSize) nowait
                    for(int idxCell = 0 ; idxCell < numberOfCells ; ++idxCell){
                        CellClass* const current_cell = iterArray[idxCell].getCurrentCell();
                        if(current_cell->hasTargetsChild()){
                            const int counter
                                = tree->getInteractionNeighbors(
                                    neighbors,
                                    neighborPositions,
                                    iterArray[idxCell].getCurrentGlobalCoordinate(),
                                    idxLevel,
                                    separationCriteria);
                            if( counter ){
                                int counterWithSrc = 0;
                                for(int idxRealNeighbors = 0 ; idxRealNeighbors < counter ; ++idxRealNeighbors ){
                                    if(neighbors[idxRealNeighbors]->hasSrcChild()){
                                        neighbors[counterWithSrc] = neighbors[idxRealNeighbors];
                                        neighborPositions[counterWithSrc] = neighborPositions[idxRealNeighbors];
                                        ++counterWithSrc;
                                    }
                                }
                                if(counterWithSrc){
                                    local_expansion_t* const target_local_exp
                                        = &(current_cell->getLocalExpansionData());
                                    const symbolic_data_t* const target_symbolic
                                        = current_cell;
                                    std::array<const multipole_t*, 342> neighbor_multipoles;
                                    std::transform(neighbors, neighbors+counterWithSrc,
                                                   neighbor_multipoles.begin(),
                                                   [](const CellClass* c) {
                                                       return (c == nullptr ? nullptr
                                                               : &(c->getMultipoleData()));
                                                   });
                                    std::array<const symbolic_data_t*, 342> neighbor_symbolics;
                                    std::transform(neighbors, neighbors+counterWithSrc,
                                                   neighbor_symbolics.begin(),
                                                   [](const CellClass* c) {return c;});

                                    kernels[omp_get_thread_num()]->M2L(
                                        target_local_exp,
                                        target_symbolic,
                                        neighbor_multipoles.data(),
                                        neighbor_symbolics.data(),
                                        neighborPositions,
                                        counterWithSrc);
                                }
                            }
                        }
                    }

                    FLOG(computationCounter.tic());
                    myThreadkernels->finishedLevelM2L(idxLevel);
                    FLOG(computationCounter.tac());
                }
                FLOG(computationCounter.tac());
                FLOG( FLog::Controller << "\t\t>> Level " << idxLevel << " = "  << counterTimeLevel.tacAndElapsed() << " s\n" );
            }
            FLOG( FLog::Controller << "\tFinished (@Downward Pass (M2L) = "  << counterTime.tacAndElapsed() << " s)\n" );
            FLOG( FLog::Controller << "\t\t Computation : " << computationCounter.cumulated() << " s\n" );
        }

        /* L2L */
        void downardPass(){
            FLOG( FLog::Controller.write("\tStart Downward Pass (L2L)\n").write(FLog::Flush); );
            FLOG(FTic counterTime);
            FLOG(FTic computationCounter);

            typename OctreeClass::Iterator octreeIterator(tree);
            octreeIterator.moveDown();

            for(int idxLevel = 2 ; idxLevel < FAbstractAlgorithm::upperWorkingLevel ; ++idxLevel){
                octreeIterator.moveDown();
            }

            typename OctreeClass::Iterator avoidGotoLeftIterator(octreeIterator);

            const int heightMinusOne = FAbstractAlgorithm::lowerWorkingLevel - 1;
            // for each levels exepted leaf level
            for(int idxLevel = FAbstractAlgorithm::upperWorkingLevel ; idxLevel < heightMinusOne ; ++idxLevel ){
                FLOG(FTic counterTimeLevel);
                int numberOfCells = 0;
                // for each cells
                do{
                    iterArray[numberOfCells] = octreeIterator;
                    ++numberOfCells;
                } while(octreeIterator.moveRight());
                avoidGotoLeftIterator.moveDown();
                octreeIterator = avoidGotoLeftIterator;

                const int chunkSize = FMath::Max(1 , numberOfCells/(omp_get_num_threads()*omp_get_num_threads()));

                FLOG(computationCounter.tic());
                #pragma omp parallel num_threads(MaxThreads)
                {
                    KernelClass * const myThreadkernels = kernels[omp_get_thread_num()];
                    #pragma omp for nowait schedule(dynamic, chunkSize)
                    for(int idxCell = 0 ; idxCell < numberOfCells ; ++idxCell){
                        if( iterArray[idxCell].getCurrentCell()->hasTargetsChild() ){
                            CellClass* target_children[8] = {};
                            CellClass** const realChild = iterArray[idxCell].getCurrentChild();
                            CellClass* const current_cell = iterArray[idxCell].getCurrentCell();
                            for(int idxChild = 0 ; idxChild < 8 ; ++idxChild){
                                if(realChild[idxChild] && realChild[idxChild]->hasTargetsChild()){
                                    target_children[idxChild] = realChild[idxChild];
                                }
                            }


                            local_expansion_t* const parent_expansion
                                = &(current_cell->getLocalExpansionData());
                            const symbolic_data_t* const parent_symbolic = current_cell;

                            std::array<local_expansion_t*, 8> child_expansions;
                            std::transform(target_children, target_children+8,
                                           child_expansions.begin(),
                                           [](CellClass* c) {
                                               return (c == nullptr ? nullptr
                                                       : &(c->getLocalExpansionData()));
                                           });
                            std::array<const symbolic_data_t*, 8> child_symbolics;
                            std::transform(target_children, target_children+8,
                                           child_symbolics.begin(),
                                           [](CellClass* c) {return c;});

                            myThreadkernels->L2L(parent_expansion,
                                                 parent_symbolic,
                                                 child_expansions.data(),
                                                 child_symbolics.data());
                        }
                    }
                }
                FLOG(computationCounter.tac());
                FLOG( FLog::Controller << "\t\t>> Level " << idxLevel << " = "  << counterTimeLevel.tacAndElapsed() << " s\n" );
            }
            FLOG( FLog::Controller << "\tFinished (@Downward Pass (L2L) = "  << counterTime.tacAndElapsed() << " s)\n" );
            FLOG( FLog::Controller << "\t\t Computation : " << computationCounter.cumulated() << " s\n" );
        }


    /** P2P */
    void directPass(const bool p2pEnabled, const bool l2pEnabled){
        FLOG( FLog::Controller.write("\tStart Direct Pass\n").write(FLog::Flush); );
        FLOG(FTic counterTime);

        int numberOfLeafs = 0;
        {
            typename OctreeClass::Iterator octreeIterator(tree);
            octreeIterator.gotoBottomLeft();
            // for each leaf
            do{
                iterArray[numberOfLeafs] = octreeIterator;
                ++numberOfLeafs;
            } while(octreeIterator.moveRight());
        }

        const int chunkSize = FMath::Max(1 , numberOfLeafs/(omp_get_num_threads()*omp_get_num_threads()));

        const int heightMinusOne = OctreeHeight - 1;
        FLOG(FTic computationCounter);
        #pragma omp parallel num_threads(MaxThreads)
        {
            KernelClass * const myThreadkernels = kernels[omp_get_thread_num()];
            // There is a maximum of 26 neighbors
            ContainerClass* neighbors[26];
            int neighborPositions[26];

            #pragma omp for schedule(dynamic, chunkSize) nowait
            for(int idxLeafs = 0 ; idxLeafs < numberOfLeafs ; ++idxLeafs){
                if( iterArray[idxLeafs].getCurrentCell()->hasTargetsChild() ){
                    if(l2pEnabled){
                        myThreadkernels->L2P(
                            &(iterArray[idxLeafs].getCurrentCell()->getLocalExpansionData()),
                            iterArray[idxLeafs].getCurrentCell(),
                            iterArray[idxLeafs].getCurrentListTargets());
                    }
                    if(p2pEnabled){
                        // need the current particles and neighbors particles
                        if(iterArray[idxLeafs].getCurrentCell()->hasSrcChild()){
                            myThreadkernels->P2P( iterArray[idxLeafs].getCurrentGlobalCoordinate(), iterArray[idxLeafs].getCurrentListTargets(),
                                          iterArray[idxLeafs].getCurrentListSrc() , neighbors, neighborPositions, 0);
                        }
                        const int counter = tree->getLeafsNeighbors(neighbors, neighborPositions, iterArray[idxLeafs].getCurrentGlobalCoordinate(),heightMinusOne);
                        myThreadkernels->P2PRemote( iterArray[idxLeafs].getCurrentGlobalCoordinate(), iterArray[idxLeafs].getCurrentListTargets(),
                                      iterArray[idxLeafs].getCurrentListSrc() , neighbors, neighborPositions, counter);
                    }
                }
            }
        }
        FLOG(computationCounter.tac());

        FLOG( counterTime.tac() );
        FLOG( FLog::Controller << "\tFinished (@Direct Pass (L2P + P2P) = "  << counterTime.tacAndElapsed() << " s)\n" );
        FLOG( FLog::Controller << "\t\t Computation L2P + P2P : " << computationCounter.elapsed() << " s\n" );

    }

};


#endif //FFMMALGORITHMTHREADTSM_HPP
