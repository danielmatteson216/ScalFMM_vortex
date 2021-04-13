// See LICENCE file at project root
#ifndef FFMMALGORITHMPERIODIC_HPP
#define FFMMALGORITHMPERIODIC_HPP

#include <array>
#include <algorithm>

#include "../Utils/FGlobal.hpp"
#include "../Utils/FGlobalPeriodic.hpp"
#include "../Utils/FAssert.hpp"
#include "../Utils/FLog.hpp"
#include "../Utils/FAlgorithmTimers.hpp"

#include "../Utils/FTic.hpp"
#include "../Utils/FMemUtils.hpp"

#include "../Containers/FOctree.hpp"
#include "../Containers/FVector.hpp"

#include "FCoreCommon.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FFmmAlgorithmPeriodic
* @brief
* Please read the license
*
* This class is a basic FMM algorithm with periodic behavior
* It just iterates on a tree and call the kernels with good arguments.
*
* Of course this class does not deallocate pointer given in arguments.
*/
template<class FReal, class OctreeClass, class CellClass, class ContainerClass, class KernelClass, class LeafClass>
class FFmmAlgorithmPeriodic : public FAbstractAlgorithm, public FAlgorithmTimers {

    using multipole_t       = typename CellClass::multipole_t;
    using local_expansion_t = typename CellClass::local_expansion_t;
    using symbolic_data_t   = CellClass;

    OctreeClass* const tree;        //< The octree to work on
    KernelClass* kernels;           //< The kernels

    const int OctreeHeight;         //< The height of the octree (real height)
    const int nbLevelsAboveRoot;    //< The nb of level the user ask to go above the tree (>= -1)
    const int offsetRealTree;       //< nbLevelsAboveRoot GetFakeLevel

    const int leafLevelSeperationCriteria;

public:
    /** The constructor need the octree and the kernels used for computation
      * @param inTree the octree to work on
      * @param inKernels the kernels to call
      * An assert is launched if one of the arguments is null
      * @param inUpperLevel this parameter defines the behavior of the periodicity refer to the main doc
      *
      */
    FFmmAlgorithmPeriodic(OctreeClass* const inTree, const int inUpperLevel = 0, const int inLeafLevelSeperationCriteria = 1)
        : tree(inTree) , kernels(nullptr), OctreeHeight(tree->getHeight()),
          nbLevelsAboveRoot(inUpperLevel), offsetRealTree(inUpperLevel + 2), leafLevelSeperationCriteria(inLeafLevelSeperationCriteria) {

        FAssertLF(tree, "tree cannot be null");
        FAssertLF(-1 <= inUpperLevel, "inUpperLevel cannot be < -1");
        FAssertLF(leafLevelSeperationCriteria < 3, "Separation criteria should be < 3");

        FAbstractAlgorithm::setNbLevelsInTree(extendedTreeHeight());

        FLOG(FLog::Controller << "FFmmAlgorithmPeriodic\n");
    }

    /** Default destructor */
    virtual ~FFmmAlgorithmPeriodic(){
    }

    void setKernel(KernelClass*const inKernel){
        kernels = inKernel;
    }




    long long int theoricalRepetition() const {
        if( nbLevelsAboveRoot == -1 ){
            // we know it is 3 (-1;+1)
            return 3;
        }
        return 6 * (1 << nbLevelsAboveRoot);
    }


    void repetitionsIntervals(FTreeCoordinate*const min, FTreeCoordinate*const max) const {
        if( nbLevelsAboveRoot == -1 ){
            // We know it is (-1;1)
            min->setPosition(-1,-1,-1);
            max->setPosition(1,1,1);
        }
        else{
            const int halfRepeated = int(theoricalRepetition()/2);
            min->setPosition(-halfRepeated,-halfRepeated,-halfRepeated);
            // if we repeat the box 8 times, we go from [-4 to 3]
            max->setPosition(halfRepeated-1,halfRepeated-1,halfRepeated-1);
        }
    }


    FReal extendedBoxWidth() const {
        if( nbLevelsAboveRoot == -1 ){
            return tree->getBoxWidth()*2;
        }
        else{
            return tree->getBoxWidth() * FReal(4<<(nbLevelsAboveRoot));
        }
    }

    FReal extendedBoxWidthBoundary() const {
        if( nbLevelsAboveRoot == -1 ){
            return tree->getBoxWidth()*4;
        }
        else{
            return tree->getBoxWidth() * FReal(8<<(nbLevelsAboveRoot));
        }
    }

    /** This function has to be used to init the kernel with correct args
      * it return the box cneter seen from a kernel point of view from the periodicity the user ask for
      * this is computed using the originalBoxWidth and originalBoxCenter given in parameter
      * @param originalBoxCenter the real system center
      * @param originalBoxWidth the real system size
      * @return the center the kernel should use
      */
    FPoint<FReal> extendedBoxCenter() const {
        if( nbLevelsAboveRoot == -1 ){
            const FReal originalBoxWidth            = tree->getBoxWidth();
            const FPoint<FReal> originalBoxCenter   = tree->getBoxCenter();
            const FReal originalBoxWidthDiv2        = originalBoxWidth/2.0;
            return FPoint<FReal>( originalBoxCenter.getX() + originalBoxWidthDiv2,
                                         originalBoxCenter.getY() + originalBoxWidthDiv2,
                                         originalBoxCenter.getZ() + originalBoxWidthDiv2);
        }
        else{
            const FReal originalBoxWidth     = tree->getBoxWidth();
            const FReal originalBoxWidthDiv2 = originalBoxWidth/2.0;
            const FPoint<FReal> originalBoxCenter   = tree->getBoxCenter();

            const FReal offset = extendedBoxWidth()/FReal(2.0);
            return FPoint<FReal>( originalBoxCenter.getX() - originalBoxWidthDiv2 + offset,
                       originalBoxCenter.getY() - originalBoxWidthDiv2 + offset,
                       originalBoxCenter.getZ() - originalBoxWidthDiv2 + offset);
        }
    }

    FPoint<FReal> extendedBoxCenterBoundary() const {
        if( nbLevelsAboveRoot == -1 ){
            const FReal originalBoxWidth            = tree->getBoxWidth();
            const FPoint<FReal> originalBoxCenter   = tree->getBoxCenter();
            const FReal originalBoxWidthDiv2        = originalBoxWidth/2.0;
            return FPoint<FReal>( originalBoxCenter.getX() + originalBoxWidthDiv2,
                                         originalBoxCenter.getY() + originalBoxWidthDiv2,
                                         originalBoxCenter.getZ() + originalBoxWidthDiv2);
        }
        else{
            const FReal originalBoxWidth     = tree->getBoxWidth();
            const FReal originalBoxWidthDiv2 = originalBoxWidth/2.0;
            const FPoint<FReal> originalBoxCenter   = tree->getBoxCenter();

            return FPoint<FReal>( originalBoxCenter.getX() + originalBoxWidthDiv2,
                       originalBoxCenter.getY() + originalBoxWidthDiv2,
                       originalBoxCenter.getZ() + originalBoxWidthDiv2);
        }
    }

    /** This function has to be used to init the kernel with correct args
      * it return the tree heigh seen from a kernel point of view from the periodicity the user ask for
      * this is computed using the originalTreeHeight given in parameter
      * @param originalTreeHeight the real tree heigh
      * @return the heigh the kernel should use
      */
    int extendedTreeHeight() const {
        // The real height
        return OctreeHeight + offsetRealTree;
    }

    int extendedTreeHeightBoundary() const {
        // The real height
        return OctreeHeight + offsetRealTree + 1;
    }

protected:
    /**
      * To execute the fmm algorithm
      * Call this function to run the complete algorithm
      */
    void executeCore(const unsigned operationsToProceed) override {
        FAssertLF(kernels, "kernels cannot be null");

        tree->forEachCell([this](CellClass* node){
                node->setLevel(node->getLevel() + offsetRealTree);
            });

        if(operationsToProceed & FFmmP2M) bottomPass();

        if(operationsToProceed & FFmmM2M) {
            upwardPass();
            // before downward pass we have to perform the periodicity
            this->processPeriodicLevels();
          }

        if(operationsToProceed & FFmmM2L){
            transferPass();
        }

        if(operationsToProceed & FFmmL2L) downardPass();

        if((operationsToProceed & FFmmP2P) || (operationsToProceed & FFmmL2P)) directPass((operationsToProceed & FFmmP2P),(operationsToProceed & FFmmL2P));

        tree->forEachCell([this](CellClass* node){
                node->setLevel(node->getLevel() - offsetRealTree);
            });

    }


    /////////////////////////////////////////////////////////////////////////////
    // P2M
    /////////////////////////////////////////////////////////////////////////////

    /** P2M */
    void bottomPass(){
        FLOG( FLog::Controller.write("\tStart Bottom Pass\n").write(FLog::Flush) );
        FLOG(FTic counterTime);
        FLOG(FTic computationCounter);

        typename OctreeClass::Iterator octreeIterator(tree);

        // Iterate on leafs
        octreeIterator.gotoBottomLeft();
        do{
            // We need the current cell that represent the leaf
            // and the list of particles
            multipole_t* const leaf_multipole
                = &(octreeIterator.getCurrentCell()->getMultipoleData());
            const symbolic_data_t* const leaf_symbolic
                = octreeIterator.getCurrentCell();
            FLOG(computationCounter.tic());
            kernels->P2M(leaf_multipole,
                         leaf_symbolic,
                         octreeIterator.getCurrentListSrc());
            FLOG(computationCounter.tac());
        } while(octreeIterator.moveRight());

        FLOG( FLog::Controller << "\tFinished (@Bottom Pass (P2M) = "  << counterTime.tacAndElapsed() << " s)\n" );
        FLOG( FLog::Controller << "\t\t Computation : " << computationCounter.cumulated() << " s\n" );
    }

    /////////////////////////////////////////////////////////////////////////////
    // Upward
    /////////////////////////////////////////////////////////////////////////////

    /** M2M */
    void upwardPass(){
        FLOG( FLog::Controller.write("\tStart Upward Pass\n").write(FLog::Flush); );
        FLOG(FTic counterTime);
        FLOG(FTic computationCounter);

        // Start from leal level - 1
        typename OctreeClass::Iterator octreeIterator(tree);
        octreeIterator.gotoBottomLeft();
        octreeIterator.moveUp();

        typename OctreeClass::Iterator avoidGotoLeftIterator(octreeIterator);

        // for each levels
        for(int idxLevel = OctreeHeight - 2 ; idxLevel > 0 ; --idxLevel ){
            FLOG(FTic counterTimeLevel);
            // for each cells
            do{
                // We need the current cell and the child
                // child is an array (of 8 child) that may be null
                multipole_t* const parent_multipole
                    = &(octreeIterator.getCurrentCell()->getMultipoleData());
                const symbolic_data_t* const parent_symbolic
                    = octreeIterator.getCurrentCell();

                CellClass** children = octreeIterator.getCurrentChildren();
                std::array<const multipole_t*, 8> child_multipoles;
                std::transform(children, children+8, child_multipoles.begin(),
                               [](CellClass* c) {
                                   return (c == nullptr ? nullptr
                                           : &(c->getMultipoleData()));
                               });
                std::array<const symbolic_data_t*, 8> child_symbolics;
                std::transform(children, children+8, child_symbolics.begin(),
                               [](CellClass* c) {return c;});
                FLOG(computationCounter.tic());
                kernels->M2M(parent_multipole,
                             parent_symbolic,
                             child_multipoles.data(),
                             child_symbolics.data());
                FLOG(computationCounter.tac());
            } while(octreeIterator.moveRight());

            FLOG( FLog::Controller << "\t\t>> Level " << idxLevel << "(" << (idxLevel + offsetRealTree) << ", " << octreeIterator.getCurrentCell()->getLevel() << ") = "  << counterTimeLevel.tacAndElapsed() << " s\n" );
            avoidGotoLeftIterator.moveUp();
            octreeIterator = avoidGotoLeftIterator;// equal octreeIterator.moveUp(); octreeIterator.gotoLeft();
        }


        FLOG( FLog::Controller << "\tFinished (@Upward Pass (M2M) = "  << counterTime.tacAndElapsed() << " s)\n" );
        FLOG( FLog::Controller << "\t\t Computation : " << computationCounter.cumulated() << " s\n" );
    }

    /////////////////////////////////////////////////////////////////////////////
    // Transfer
    /////////////////////////////////////////////////////////////////////////////

    /** M2L L2L */
    void transferPass(){
        FLOG( FLog::Controller.write("\tStart Downward Pass (M2L)\n").write(FLog::Flush); );
        FLOG(FTic counterTime);
        FLOG(FTic computationCounter);

        typename OctreeClass::Iterator octreeIterator(tree);
        typename OctreeClass::Iterator avoidGotoLeftIterator(octreeIterator);

        const CellClass* neighbors[342];
        int neighborPositions[342];

        // for each levels
        for(int idxLevel = 1 ; idxLevel < OctreeHeight ; ++idxLevel ){
            FLOG(FTic counterTimeLevel);
            const int fakeLevel = idxLevel + offsetRealTree;
            const int separationCriteria = (idxLevel != OctreeHeight-1 ? 1 : leafLevelSeperationCriteria);
            // for each cells
            do{
                const int counter = tree->getPeriodicInteractionNeighbors(neighbors, neighborPositions, octreeIterator.getCurrentGlobalCoordinate(), idxLevel, AllDirs, separationCriteria);

                if(counter == 0) {
                    continue;
                }

                local_expansion_t* const target_local_exp
                    = &(octreeIterator.getCurrentCell()->getLocalExpansionData());
                const symbolic_data_t* const target_symbolic
                    = octreeIterator.getCurrentCell();
                std::array<const multipole_t*, 342> neighbor_multipoles;
                std::array<const symbolic_data_t*, 342> neighbor_symbolics;
                std::transform(neighbors, neighbors+counter, neighbor_multipoles.begin(),
                               [](const CellClass* c) {
                                   return (c == nullptr ? nullptr
                                           : &(c->getMultipoleData()));
                               });
                std::transform(neighbors, neighbors+counter, neighbor_symbolics.begin(),
                               [](const CellClass* c) {return c;});


                FLOG(computationCounter.tic());
                kernels->M2L(
                    target_local_exp,
                    target_symbolic,
                    neighbor_multipoles.data(),
                    neighbor_symbolics.data(),
                    neighborPositions,
                    counter);
                FLOG(computationCounter.tac());
            } while(octreeIterator.moveRight());
            FLOG( FLog::Controller << "\t\t>> Level " << idxLevel << "(" << fakeLevel << ", " << octreeIterator.getCurrentCell()->getLevel() << ") = "  << counterTimeLevel.tacAndElapsed() << " s\n" );
            avoidGotoLeftIterator.moveDown();
            octreeIterator = avoidGotoLeftIterator;

            FLOG(computationCounter.tic());
            kernels->finishedLevelM2L(fakeLevel);
            FLOG(computationCounter.tac());
        }
        FLOG( FLog::Controller << "\tFinished (@Downward Pass (M2L) = "  << counterTime.tacAndElapsed() << " s)\n" );
        FLOG( FLog::Controller << "\t\t Computation : " << computationCounter.cumulated() << " s\n" );
    }

    /////////////////////////////////////////////////////////////////////////////
    // Downward
    /////////////////////////////////////////////////////////////////////////////


    void downardPass(){ // second L2L
        FLOG( FLog::Controller.write("\tStart Downward Pass (L2L)\n").write(FLog::Flush); );
        FLOG(FTic counterTime);
        FLOG(FTic computationCounter );

        typename OctreeClass::Iterator octreeIterator(tree);
        typename OctreeClass::Iterator avoidGotoLeftIterator(octreeIterator);

        const int heightMinusOne = OctreeHeight - 1;
        // for each levels exepted leaf level
        for(int idxLevel = 1 ; idxLevel < heightMinusOne ; ++idxLevel ){
            FLOG(FTic counterTimeLevel);

            // for each cells
            do{
                local_expansion_t* const parent_local_exp
                    = &(octreeIterator.getCurrentCell()->getLocalExpansionData());
                const symbolic_data_t* const parent_symbolic
                    = octreeIterator.getCurrentCell();
                CellClass** children = octreeIterator.getCurrentChildren();
                std::array<local_expansion_t*, 8> child_local_expansions;
                std::transform(children, children+8, child_local_expansions.begin(),
                               [](CellClass* c) {return (c == nullptr ? nullptr
                                                         : &(c->getLocalExpansionData()));
                               });
                std::array<symbolic_data_t*, 8> child_symbolics;
                std::transform(children, children+8, child_symbolics.begin(),
                               [](CellClass* c) {return c;});
                FLOG(computationCounter.tic());
                kernels->L2L(
                    parent_local_exp,
                    parent_symbolic,
                    child_local_expansions.data(),
                    child_symbolics.data()
                    );
                FLOG(computationCounter.tac());
            } while(octreeIterator.moveRight());

            FLOG( FLog::Controller << "\t\t>> Level " << idxLevel << "(" << (idxLevel + offsetRealTree) << ", " << octreeIterator.getCurrentCell()->getLevel() << ") = "  << counterTimeLevel.tacAndElapsed() << " s\n" );
            avoidGotoLeftIterator.moveDown();
            octreeIterator = avoidGotoLeftIterator;
        }

        FLOG( FLog::Controller << "\tFinished (@Downward Pass (L2L) = "  << counterTime.tacAndElapsed() << " s)\n" );
        FLOG( FLog::Controller << "\t\t Computation : " << computationCounter.cumulated() << " s\n" );


    }

    /////////////////////////////////////////////////////////////////////////////
    // Direct
    /////////////////////////////////////////////////////////////////////////////

    /** P2P */
    void directPass(const bool p2pEnabled, const bool l2pEnabled){
        FLOG( FLog::Controller.write("\tStart Direct Pass\n").write(FLog::Flush); );
        FLOG(FTic counterTime);
        FLOG(FTic computationCounterL2P);
        FLOG(FTic computationCounterP2P);

        const int heightMinusOne = OctreeHeight - 1;
        const FReal boxWidth = tree->getBoxWidth();

        typename OctreeClass::Iterator octreeIterator(tree);
        octreeIterator.gotoBottomLeft();
        // There is a maximum of 26 neighbors
        ContainerClass* neighbors[26];
        FTreeCoordinate offsets[26];
        int neighborPositions[26];
        bool hasPeriodicLeaves;
        // for each leafs
        do{
            if(l2pEnabled){
                FLOG(computationCounterL2P.tic());
                kernels->L2P(
                    &(octreeIterator.getCurrentCell()->getLocalExpansionData()),
                    octreeIterator.getCurrentCell(),
                    octreeIterator.getCurrentListTargets()
                    );
                FLOG(computationCounterL2P.tac());
            }
            if(p2pEnabled){
                // need the current particles and neighbors particles
                const FTreeCoordinate centerOfLeaf = octreeIterator.getCurrentGlobalCoordinate();
                const int counter = tree->getPeriodicLeafsNeighbors( neighbors, neighborPositions, offsets,
                                                                     hasPeriodicLeaves, centerOfLeaf,
                                                                     heightMinusOne, AllDirs);
                int periodicNeighborsCounter = 0;

                if(hasPeriodicLeaves){
                    ContainerClass* periodicNeighbors[26];
                    int periodicNeighborPositions[26];

                    for(int idxNeig = 0 ; idxNeig < counter ; ++idxNeig){
                        if( !offsets[idxNeig].equals(0,0,0) ){
                            FAssertLF(octreeIterator.getCurrentListTargets() != neighbors[idxNeig]);
                            // Put periodic neighbors into other array
                            FReal*const positionsX = neighbors[idxNeig]->getPositions()[0];
                            FReal*const positionsY = neighbors[idxNeig]->getPositions()[1];
                            FReal*const positionsZ = neighbors[idxNeig]->getPositions()[2];

                            for(FSize idxPart = 0; idxPart < neighbors[idxNeig]->getNbParticles() ; ++idxPart){
                                positionsX[idxPart] += boxWidth * FReal(offsets[idxNeig].getX());
                                positionsY[idxPart] += boxWidth * FReal(offsets[idxNeig].getY());
                                positionsZ[idxPart] += boxWidth * FReal(offsets[idxNeig].getZ());
                            }

                            offsets[periodicNeighborsCounter] = offsets[idxNeig];
                            periodicNeighbors[periodicNeighborsCounter] = neighbors[idxNeig];
                            periodicNeighborPositions[periodicNeighborsCounter] = neighborPositions[idxNeig];
                            ++periodicNeighborsCounter;
                        }
                        else{
                            neighbors[idxNeig-periodicNeighborsCounter] = neighbors[idxNeig];
                            neighborPositions[idxNeig-periodicNeighborsCounter] = neighborPositions[idxNeig];
                        }
                    }

                    FLOG(computationCounterP2P.tic());
                    kernels->P2PRemote(octreeIterator.getCurrentGlobalCoordinate(),octreeIterator.getCurrentListTargets(),
                                 octreeIterator.getCurrentListSrc(), periodicNeighbors, periodicNeighborPositions, periodicNeighborsCounter);
                    FLOG(computationCounterP2P.tac());

                    for(int idxNeig = 0 ; idxNeig < periodicNeighborsCounter ; ++idxNeig){
                        FReal*const positionsX = periodicNeighbors[idxNeig]->getPositions()[0];
                        FReal*const positionsY = periodicNeighbors[idxNeig]->getPositions()[1];
                        FReal*const positionsZ = periodicNeighbors[idxNeig]->getPositions()[2];

                        for(FSize idxPart = 0; idxPart < periodicNeighbors[idxNeig]->getNbParticles() ; ++idxPart){
                            positionsX[idxPart] -= boxWidth * FReal(offsets[idxNeig].getX());
                            positionsY[idxPart] -= boxWidth * FReal(offsets[idxNeig].getY());
                            positionsZ[idxPart] -= boxWidth * FReal(offsets[idxNeig].getZ());
                        }
                    }
                }

                FLOG(computationCounterP2P.tic());
                kernels->P2P(octreeIterator.getCurrentGlobalCoordinate(),octreeIterator.getCurrentListTargets(),
                             octreeIterator.getCurrentListSrc(), neighbors, neighborPositions, counter - periodicNeighborsCounter);
                FLOG(computationCounterP2P.tac());
            }
        } while(octreeIterator.moveRight());


        FLOG( FLog::Controller << "\tFinished (@Direct Pass (L2P + P2P) = "  << counterTime.tacAndElapsed() << " s)\n" );
        FLOG( FLog::Controller << "\t\t Computation L2P : " << computationCounterL2P.cumulated() << " s\n" );
        FLOG( FLog::Controller << "\t\t Computation P2P : " << computationCounterP2P.cumulated() << " s\n" );

    }

    /////////////////////////////////////////////////////////////////////////////
    // Periodic levels = levels <= 0
    /////////////////////////////////////////////////////////////////////////////

    /** Get the index of a interaction neighbors (for M2L)
      * @param x the x position in the interactions (from -3 to +3)
      * @param y the y position in the interactions (from -3 to +3)
      * @param z the z position in the interactions (from -3 to +3)
      * @return the index (from 0 to 342)
      */
    int neighIndex(const int x, const int y, const int z) const {
        return (((x+3)*7) + (y+3))*7 + (z + 3);
    }

    /** Periodicity Core
      * This function is split in several part:
      * 1 - special case managment
      * There is nothing to do if nbLevelsAboveRoot == -1 and only
      * a M2L if nbLevelsAboveRoot == 0
      * 2 - if nbLevelsAboveRoot > 0
      * First we compute M2M and special M2M if needed for the border
      * Then the M2L by taking into account the periodicity directions
      * Then the border by using the precomputed M2M
      * Finally the L2L
      */
    void processPeriodicLevels(){
        FLOG( FLog::Controller.write("\tStart Periodic Pass\n").write(FLog::Flush); );
        FLOG(FTic counterTime);

        if( nbLevelsAboveRoot != -1 ){
            // we will use offsetRealTree-1 cells but for simplicity allocate offsetRealTree
            // upperCells[offsetRealTree-1] is root cell
            CellClass*const upperCells = new CellClass[offsetRealTree];
            for(int i = 0; i < offsetRealTree; ++i) {
                upperCells[i].setLevel(i+1);
            }

            {  // Build Expansion at level 0
                typename OctreeClass::Iterator octreeIterator(tree);
                octreeIterator.gotoLeft();

                multipole_t* const real_tree_root_multipole
                    = &((&upperCells[offsetRealTree-1])->getMultipoleData());
                const symbolic_data_t* const real_tree_root_symbolic
                    = &(upperCells[offsetRealTree-1]);

                CellClass** children = octreeIterator.getCurrentBox();

                std::array<const multipole_t*, 8> level_1_multipoles;
                std::transform(children, children+8, level_1_multipoles.begin(),
                               [](CellClass* c) {
                                   return (c == nullptr ? nullptr
                                           : &(c->getMultipoleData()));
                               });
                std::array<const symbolic_data_t*, 8> level_1_symbolics;
                std::transform(children, children+8, level_1_symbolics.begin(),
                               [](CellClass* c) {return c;});
                //
                kernels->M2M(
                    real_tree_root_multipole,
                    real_tree_root_symbolic,
                    level_1_multipoles.data(),
                    level_1_symbolics.data()
                    );
            }
            {
                CellClass* children[8] = {};
                for(int idxLevel = offsetRealTree-1 ; idxLevel > 1  ; --idxLevel){
                    // Copy virtual cell at given level into virtual children
                    FMemUtils::setall(children,&upperCells[idxLevel],8);


                    multipole_t* const virtual_parent_multipole
                        = &((&upperCells[idxLevel-1])->getMultipoleData());
                    const symbolic_data_t* const virtual_parent_symbolic
                        = &(upperCells[idxLevel-1]);

                    std::array<const multipole_t*, 8> virtual_child_multipoles;
                    std::transform(children, children+8, virtual_child_multipoles.begin(),
                                   [](CellClass* c) {return &(c->getMultipoleData());});
                    std::array<const symbolic_data_t*, 8> virtual_child_symbolics;
                    std::transform(children, children+8, virtual_child_symbolics.begin(),
                                   [](CellClass* c) {return c;});

                    kernels->M2M(
                        virtual_parent_multipole,
                        virtual_parent_symbolic,
                        virtual_child_multipoles.data(),
                        virtual_child_symbolics.data()
                        );
                }
            }
            //
            //        M2L PASS
            //
            CellClass*const downerCells = new CellClass[offsetRealTree];
            for(int i = 0; i < offsetRealTree; ++i) {
                downerCells[i].setLevel(i+1); // TODO: check that the level is correct (see other todos under)
            }
            // Build a virtual environment for the M2L at the topmost level
            {
                const int idxUpperLevel = 2;

                const CellClass* neighbors[342]{};
                int neighborPositions[342]{};
                int counter = 0;

                // The cells are all the same, duplicate the central cell and
                // create fake neighbour positions
                for(int idxX = -3 ; idxX <= 2 ; ++idxX){
                    for(int idxY = -3 ; idxY <= 2 ; ++idxY){
                        for(int idxZ = -3 ; idxZ <= 2 ; ++idxZ){
                            if( FMath::Abs(idxX) > 1 || FMath::Abs(idxY) > 1 || FMath::Abs(idxZ) > 1){
                                neighbors[counter]         = &upperCells[idxUpperLevel-1];
                                neighborPositions[counter] = neighIndex(idxX,idxY,idxZ);
                                ++counter;
                            }
                        }
                    }
                }

                local_expansion_t* const target_local_exp
                    = &(downerCells[idxUpperLevel-1].getLocalExpansionData());
                const symbolic_data_t* const target_symbolic
                    = &downerCells[idxUpperLevel-1];

                std::array<const multipole_t*, 342>     neighbor_multipoles{};
                std::array<const symbolic_data_t*, 342> neighbor_symbolics{};
                std::transform(neighbors, neighbors+counter, neighbor_multipoles.begin(),
                               [](const CellClass* c) {return &(c->getMultipoleData());});
                std::transform(neighbors, neighbors+counter, neighbor_symbolics.begin(),
                               [](const CellClass* c) {return c;});
                // compute M2L
                kernels->M2L(
                    target_local_exp,
                    target_symbolic,
                    neighbor_multipoles.data(),
                    neighbor_symbolics.data(),
                    neighborPositions,
                    counter);
            }


            // note: the triple loop bounds are not the same than for the previous piece of code
            // which handles the topmost virtual cell

            for(int idxUpperLevel = 3 ; idxUpperLevel <= offsetRealTree ; ++idxUpperLevel){
                const CellClass* neighbors[342];
                int neighborPositions[342];
                int counter = 0;
                for(int idxX = -2 ; idxX <= 3 ; ++idxX){
                    for(int idxY = -2 ; idxY <= 3 ; ++idxY){
                        for(int idxZ = -2 ; idxZ <= 3 ; ++idxZ){
                            if( FMath::Abs(idxX) > 1 || FMath::Abs(idxY) > 1 || FMath::Abs(idxZ) > 1){
                                neighbors[counter]         = &upperCells[idxUpperLevel-1];
                                neighborPositions[counter] = neighIndex(idxX,idxY,idxZ);
                                ++counter;
                            }
                        }
                    }
                }

                local_expansion_t* const target_local_exp
                    = &(downerCells[idxUpperLevel-1].getLocalExpansionData());
                const symbolic_data_t* const target_symbolic
                    = &downerCells[idxUpperLevel-1];

                std::array<const multipole_t*, 342> neighbor_multipoles;
                std::array<const symbolic_data_t*, 342> neighbor_symbolics;
                std::transform(neighbors, neighbors+counter, neighbor_multipoles.begin(),
                               [](const CellClass* c) {
                                   return (c == nullptr ? nullptr
                                           : &(c->getMultipoleData()));
                               });
                std::transform(neighbors, neighbors+counter, neighbor_symbolics.begin(),
                               [](const CellClass* c) {return c;});

                // compute M2L
                kernels->M2L(
                    target_local_exp,
                    target_symbolic,
                    neighbor_multipoles.data(),
                    neighbor_symbolics.data(),
                    neighborPositions,
                    counter);
            }
            // Run the L2L for all but the lowest virtual levels
            //  2 is the highest level
            //  offsetRealTree-1 is the root level
            {

                std::array<local_expansion_t*, 8> virtual_child_local_exps = {};
                std::array<symbolic_data_t*, 8> virtual_child_symbolics = {};
                for(int idxLevel = 2 ; idxLevel < offsetRealTree-1  ; ++idxLevel){
                    local_expansion_t* const virtual_parent_local_exp
                        = &(downerCells[idxLevel-1].getLocalExpansionData());
                    const symbolic_data_t* const virtual_parent_symbolic
                        = &downerCells[idxLevel-1];

                    virtual_child_local_exps[0]
                        = &(downerCells[idxLevel].getLocalExpansionData());
                    virtual_child_symbolics[0]
                        = &(downerCells[idxLevel]);

                    kernels->L2L(
                        virtual_parent_local_exp,
                        virtual_parent_symbolic,
                        virtual_child_local_exps.data(),
                        virtual_child_symbolics.data()
                        );
                }
            }
            // Run the L2L for the lowest virtual level
            if(offsetRealTree >2){
                std::array<local_expansion_t*, 8> virtual_child_local_exps  {};
                std::array<symbolic_data_t*, 8> virtual_child_symbolics  {};
                const int idxLevel = offsetRealTree-1;
                local_expansion_t* const virtual_parent_local_exp
                    = &(downerCells[idxLevel-1].getLocalExpansionData());
                const symbolic_data_t* const virtual_parent_symbolic
                    = &downerCells[idxLevel-1];

                virtual_child_local_exps[7] = &(downerCells[idxLevel].getLocalExpansionData());
                virtual_child_symbolics[7]  = &downerCells[idxLevel];
                kernels->L2L(
                    virtual_parent_local_exp,
                    virtual_parent_symbolic,
                    virtual_child_local_exps.data(),
                    virtual_child_symbolics.data()
                    );
            }

            // Run L2L from the lowest virtual level to the highest real tree level
            {
                typename OctreeClass::Iterator octreeIterator(tree);
                octreeIterator.gotoLeft();

                local_expansion_t* const virtual_parent_local_exp
                    = &(downerCells[offsetRealTree-1].getLocalExpansionData());
                const symbolic_data_t* const virtual_parent_symbolic
                    = &downerCells[offsetRealTree-1];

                CellClass** children = octreeIterator.getCurrentBox();
                std::array<local_expansion_t*, 8> child_local_expansions;
                std::transform(children, children+8, child_local_expansions.begin(),
                               [](CellClass* c) {return (c == nullptr ? nullptr
                                                         : &(c->getLocalExpansionData()));
                               });
                std::array<symbolic_data_t*, 8> child_symbolics;
                std::transform(children, children+8, child_symbolics.begin(),
                               [](CellClass* c) {return c;});

                kernels->L2L(
                    virtual_parent_local_exp,
                    virtual_parent_symbolic,
                    child_local_expansions.data(),
                    child_symbolics.data()
                    );
            }
            delete[] upperCells;
            delete[] downerCells;
        }

        FLOG( FLog::Controller << "\tFinished (@Periodic = "  << counterTime.tacAndElapsed() << " s)\n" );
    }


};


#endif // FFMMALGORITHMPERIODIC_HPP
