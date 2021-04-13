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

#include "Components/FBasicCell.hpp"

#include "Components/FSimpleLeaf.hpp"

#include "Files/FBasicLoader.hpp"

#include "Components/FBasicParticleContainer.hpp"

#include "Utils/FParameterNames.hpp"

/**
  * In this file we show an example of FBasicLoader use
* Inserting 2000000 particles ...
* Done  (5.77996).
* Deleting particles ...
* Done  (0.171918).
  */

int main(int argc, char ** argv){
    FHelpDescribeAndExit(argc, argv, "Load a file and put the particles in a tree",
                         FParameterDefinitions::InputFile, FParameterDefinitions::OctreeHeight);

    typedef double FReal;
    typedef FBasicParticleContainer<FReal,0,FReal>      ContainerClass;
    typedef FSimpleLeaf<FReal, ContainerClass >                     LeafClass;
    typedef FOctree< FReal, FBasicCell, ContainerClass , LeafClass >  OctreeClass;
    ///////////////////////What we do/////////////////////////////
    std::cout << ">> This executable is useless to execute.\n";
    std::cout << ">> It is only interesting to wath the code to understand\n";
    std::cout << ">> how to use the basic loader.\n";
    //////////////////////////////////////////////////////////////

    // Use testLoaderCreate.exe to create this file
    FTic counter;

    const char* const filename = FParameters::getStr(argc,argv,FParameterDefinitions::InputFile.options, "../Data/test20k.basic");
    std::cout << "Opening : " << filename << "\n";

    // open basic particles loader
    FBasicLoader<FReal> loader(filename);
    if(!loader.isOpen()){
        std::cout << "Loader Error, " << filename << "is missing\n";
        return 1;
    }

    {
        // otree
        OctreeClass tree(FParameters::getValue(argc,argv,FParameterDefinitions::OctreeHeight.options, 5), FParameters::getValue(argc,argv,FParameterDefinitions::OctreeSubHeight.options, 3),
                         loader.getBoxWidth(), loader.getCenterOfBox());

        // -----------------------------------------------------
        std::cout << "Inserting " << loader.getNumberOfParticles() << " particles ..." << std::endl;
        counter.tic();

        FPoint<FReal> particlePosition;
        for(FSize idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
            loader.fillParticle(&particlePosition);
            tree.insert(particlePosition);
        }

        counter.tac();
        std::cout << "Done  " << "(" << counter.elapsed() << ")." << std::endl;

        // -----------------------------------------------------
        std::cout << "Deleting particles ..." << std::endl;
        counter.tic();
    }
    counter.tac();
    std::cout << "Done  " << "(" << counter.elapsed() << ")." << std::endl;
    // -----------------------------------------------------

    return 0;
}



