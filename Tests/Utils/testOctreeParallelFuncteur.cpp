// See LICENCE file at project root

#include <iostream>

#include <cstdio>
#include <cstdlib>
#include <time.h>

#include "Utils/FTic.hpp"
#include "Utils/FParameters.hpp"

#include "Containers/FOctree.hpp"
#include "Containers/FParForEachOctree.hpp"
#include "Containers/FVector.hpp"

#include "Utils/FAssert.hpp"
#include "Utils/FPoint.hpp"
#include "Utils/FParObject.hpp"

#include "Components/FBasicParticleContainer.hpp"
#include "Components/FBasicCell.hpp"
#include "Components/FSimpleLeaf.hpp"

#include "Files/FRandomLoader.hpp"

#include "Utils/FParameterNames.hpp"

/**
* In this file we show how to use octree's functeur
*/

int main(int argc, char ** argv){
    FHelpDescribeAndExit(argc, argv,
                         "Show how to use an octree functeur parallelized (only the code is interesting)",
                         FParameterDefinitions::NbParticles);

    typedef double FReal;
    typedef FBasicCell CellClass;
    typedef FBasicParticleContainer<FReal,0,FReal>      ContainerClass;
    typedef FSimpleLeaf<FReal, ContainerClass >                     LeafClass;
    typedef FOctree<FReal, CellClass, ContainerClass , LeafClass >  OctreeClass;
    ///////////////////////What we do/////////////////////////////
    std::cout << ">> This executable is useless to execute.\n";
    std::cout << ">> It is only interesting to wath the code to understand\n";
    std::cout << ">> how to use the Octree\n";
    //////////////////////////////////////////////////////////////
    const FSize NbPart = FParameters::getValue(argc,argv,FParameterDefinitions::NbParticles.options, FSize(2000));
    FTic counter;

    FRandomLoader<FReal> loader(NbPart, 1, FPoint<FReal>(0.5,0.5,0.5), 1);
    OctreeClass tree(10, 3, loader.getBoxWidth(), loader.getCenterOfBox());

    // -----------------------------------------------------
    std::cout << "Creating and inserting " << NbPart << " particles ..." << std::endl;
    counter.tic();

    {
        FPoint<FReal> particlePosition;
        for(FSize idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
            loader.fillParticle(&particlePosition);
            tree.insert(particlePosition);
        }
    }

    counter.tac();
    std::cout << "Done  " << "(" << counter.elapsed() << ")." << std::endl;
    // -----------------------------------------------------

    // Call a function on each leaf
    FParObject<long> nbParticles(0);
    FParForEachOctree::forEachLeaf(&tree, [&](LeafClass* leaf){
        nbParticles.getMine() += leaf->getSrc()->getNbParticles();
    });
    const long totalNbParticles = nbParticles.reduce([](long v1, long v2) -> long {return v1 + v2;});
    std::cout << "There are " << totalNbParticles << " particles " << std::endl;

    // Call a function on each cell
    FParObject<long> nbCells;
    nbCells = 0;
    FParForEachOctree::forEachCell(&tree, [&nbCells](CellClass* /*cell*/){
        nbCells.getMine() += 1;
    });
    std::cout << "There are " << nbCells.reduce([](long v1, long v2) -> long {return v1 + v2;}) << " cells " << std::endl;

    // To get cell and particles at leaf level
    FParForEachOctree::forEachCellLeaf(&tree, [&](CellClass* cell, LeafClass* /*leaf*/){
        cell->resetToInitialState();
    });

    // If we need an array of long
    FParArray<long> arrayExample(10);
    arrayExample.resize(50);

    FParForEachOctree::forEachCellLeaf(&tree, [&](CellClass* /*cell*/, LeafClass* /*leaf*/){
        arrayExample.getMine()[0] += 10;
    });

    return 0;
}




