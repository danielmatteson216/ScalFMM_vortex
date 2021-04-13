// See LICENCE file at project root

#include <iostream>

#include <cstdio>
#include <cstdlib>
#include <time.h>


#include "Utils/FParameters.hpp"

#include "Containers/FOctree.hpp"
#include "Containers/FVector.hpp"
#include "Components/FSimpleLeaf.hpp"

#include "Utils/FAssert.hpp"
#include "Utils/FPoint.hpp"

#include "Components/FBasicParticleContainer.hpp"
#include "Components/FBasicCell.hpp"

#include "Utils/FTic.hpp"

#include "Utils/FParameterNames.hpp"

/**
* In this file we show how to use octree with iteration
* This is a good example to understand FOctree::Iterator.
*/

int main(int argc, char ** argv){
    FHelpDescribeAndExit(argc, argv,
                         "Show how to iterate on an octree (only the code is interesting)",
                         FParameterDefinitions::NbParticles, FParameterDefinitions::OctreeHeight,
                         FParameterDefinitions::OctreeSubHeight);

    typedef double FReal;
    typedef FBasicParticleContainer<FReal,0,FReal>     ContainerClass;
    typedef FSimpleLeaf<FReal, ContainerClass >                     LeafClass;
    typedef FOctree<FReal, FBasicCell, ContainerClass , LeafClass >  OctreeClass;
    ///////////////////////What we do/////////////////////////////
    std::cout << ">> This executable is useless to execute.\n";
    std::cout << ">> It is only interesting to wath the code to understand\n";
    std::cout << ">> how to use octree iterator.\n";
    //////////////////////////////////////////////////////////////

    const int NbLevels = FParameters::getValue(argc,argv,FParameterDefinitions::OctreeHeight.options, 9);
    const int NbSubLevels = FParameters::getValue(argc,argv,FParameterDefinitions::OctreeSubHeight.options, 3);
    const FSize NbPart = FParameters::getValue(argc,argv,FParameterDefinitions::NbParticles.options, FSize(2000000));

    FTic counterTime;

    srand48 ( 1 ); // volontary set seed to constant
    // -----------------------------------------------------

    OctreeClass tree(NbLevels, NbSubLevels, 1.0, FPoint<FReal>(0.5,0.5,0.5));

    // -----------------------------------------------------
    std::cout << "Creating and inserting " << NbPart << " particles ..." << std::endl;
    counterTime.tic();
    {
        FPoint<FReal> particle;
        for(long idxPart = 0 ; idxPart < NbPart ; ++idxPart){
            particle.setPosition(FReal(drand48()),FReal(drand48()),FReal(drand48()));
            tree.insert(particle);
        }
    }
    counterTime.tac();
    std::cout << "Done  " << "(" << counterTime.elapsed() << "s)." << std::endl;

    // -----------------------------------------------------
    {
        std::cout << "Itering on Cells ..." << std::endl;
        counterTime.tic();

        OctreeClass::Iterator octreeIterator(&tree);
        octreeIterator.gotoBottomLeft();
        for(int idxLevel = NbLevels - 1 ; idxLevel >= 1 ; --idxLevel ){
            int counter = 0;
            do{
                ++counter;
                //counter += octreeIterator.getCurrentList()->getSize();
            } while(octreeIterator.moveRight());           
            octreeIterator.moveUp();
            octreeIterator.gotoLeft();
            std::cout << "Cells at this level " << counter << " ...\n";
        }
        counterTime.tac();
        std::cout << "Done  " << "(" << counterTime.elapsed() << "s)." << std::endl;
    }
    // -----------------------------------------------------
    {
        std::cout << "Itering on particles fast ..." << std::endl;
        counterTime.tic();

        OctreeClass::Iterator octreeIterator(&tree);
        octreeIterator.gotoBottomLeft();

        OctreeClass::Iterator avoidGoLeft(octreeIterator);

        for(int idx = 0 ; idx < NbLevels - 1; ++idx ){
            do{
            } while(octreeIterator.moveRight());
            avoidGoLeft.moveUp();
            octreeIterator = avoidGoLeft;
        }
        counterTime.tac();
        std::cout << "Done  " << "(" << counterTime.elapsed() << "s)." << std::endl;
    }
    // -----------------------------------------------------

    return 0;
}



