// See LICENCE file at project root

#include <iostream>

#include <cstdio>
#include <cstdlib>
#include <time.h>

#include "Utils/FTic.hpp"

#include "Containers/FOctree.hpp"
#include "Containers/FVector.hpp"

#include "Utils/FAssert.hpp"
#include "Utils/FPoint.hpp"

#include "Components/FBasicCell.hpp"

#include "Components/FTypedLeaf.hpp"

#include "Files/FFmaTsmLoader.hpp"

#include "Components/FBasicParticleContainer.hpp"

#include "Utils/FParameters.hpp"

#include "Utils/FParameterNames.hpp"


int main(int argc, char ** argv ){
    FHelpDescribeAndExit(argc, argv, "Load a file and put the particles in a tree with the TSM (target source model)",
                         FParameterDefinitions::InputFile, FParameterDefinitions::OctreeHeight);

    typedef double FReal;
    typedef FBasicParticleContainer<FReal,1,FReal>     ContainerClass;
    typedef FTypedLeaf< FReal, ContainerClass >                     LeafClass;
    typedef FOctree< FReal, FBasicCell, ContainerClass , LeafClass >  OctreeClass;
    ///////////////////////What we do/////////////////////////////
    std::cout << ">> This executable is useless to execute.\n";
    std::cout << ">> It is only interesting to wath the code to understand\n";
    std::cout << ">> how to use the Tsm loader\n";
    //////////////////////////////////////////////////////////////

    // Use testLoaderCreate.exe to create this file
    FTic counter;
    const char* const filename = FParameters::getStr(argc,argv,FParameterDefinitions::InputFile.options, "../Data/test20k.tsm.fma");
    std::cout << "Opening : " << filename << "\n";

    // open basic particles loader
    FFmaTsmLoader<FReal> loader(filename);
    if(!loader.isOpen()){
        std::cout << "Loader Error, " << filename << "is missing\n";
        return 1;
    }
    {
        // otree
        OctreeClass tree(FParameters::getValue(argc,argv,FParameterDefinitions::OctreeHeight.options, 5), FParameters::getValue(argc,argv,FParameterDefinitions::OctreeSubHeight.options, 3),
                         loader.getBoxWidth(),loader.getCenterOfBox());

        // -----------------------------------------------------
        std::cout << "Inserting " << loader.getNumberOfParticles() << " particles ..." << std::endl;
        counter.tic();

        FPoint<FReal> particlePosition;
        FReal physicalValue = 0.0;
        FParticleType particleType;
        for(FSize idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
            loader.fillParticle(&particlePosition,&physicalValue, &particleType);
            tree.insert(particlePosition, particleType, physicalValue);
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



