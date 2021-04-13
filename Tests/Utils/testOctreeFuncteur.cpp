// See LICENCE file at project root

#include <iostream>

#include <cstdio>
#include <cstdlib>
#include <time.h>

#include "Utils/FTic.hpp"
#include "Utils/FParameters.hpp"

#include "Containers/FOctree.hpp"
#include "Containers/FVector.hpp"

#include "Utils/FAssert.hpp"
#include "Utils/FPoint.hpp"

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
                         "Show how to use an octree functeur (only the code is interesting)",
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
    long nbParticles = 0;
    tree.forEachLeaf([&](LeafClass* leaf){
        nbParticles += leaf->getSrc()->getNbParticles();
    });
    std::cout << "There are " << nbParticles << " particles " << std::endl;

    // Call a function on each cell
    long nbCells = 0;
    tree.forEachCell([&nbCells](CellClass* /*cell*/){
        nbCells += 1;
    });
    std::cout << "There are " << nbCells << " cells " << std::endl;

    // To get cell and particles at leaf level
    tree.forEachCellLeaf([&](CellClass* /*cell*/, LeafClass* /*leaf*/){
    });

    return 0;
}



