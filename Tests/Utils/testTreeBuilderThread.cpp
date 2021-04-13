// See LICENCE file at project root

#include <limits>
#include <iostream>


#include "Components/FSimpleLeaf.hpp"
#include "Containers/FVector.hpp"

#include "Containers/FOctree.hpp"

#include "Kernels/P2P/FP2PParticleContainer.hpp"

#include "Kernels/Rotation/FRotationKernel.hpp"
#include "Kernels/Rotation/FRotationCell.hpp"

#include "Utils/FMath.hpp"
#include "Utils/FMemUtils.hpp"
#include "Utils/FParameters.hpp"

#include "Files/FFmaGenericLoader.hpp"

#include "Utils/FParameterNames.hpp"

#include "Files/FTreeBuilder.hpp"

int main(int argc, char** argv){
    FHelpDescribeAndExit(argc, argv,
                         "Test the insertion of particle using TreeBuilder.",
                         FParameterDefinitions::InputFile, FParameterDefinitions::OctreeHeight,
                         FParameterDefinitions::OctreeSubHeight);

    typedef double FReal;

    static const int P = 9;

    typedef FRotationCell<FReal,P>               CellClass;
    typedef FP2PParticleContainer<FReal>          ContainerClass;

    typedef FSimpleLeaf<FReal, ContainerClass >                     LeafClass;
    typedef FOctree<FReal, CellClass, ContainerClass , LeafClass >  OctreeClass;

    const int NbLevels      = FParameters::getValue(argc,argv,FParameterDefinitions::OctreeHeight.options, 5);
    const int SizeSubLevels = FParameters::getValue(argc,argv,FParameterDefinitions::OctreeSubHeight.options, 3);
    const char* const filename = FParameters::getStr(argc,argv,FParameterDefinitions::InputFile.options, "../Data/test20k.fma");

    // -----------------------------------------------------

    FFmaGenericLoader<FReal> loader(filename);
    if(!loader.isOpen()){
        std::cout << "Loader Error, " << filename << " is missing\n";
        return 1;
    }

    // -----------------------------------------------------

    {
        ContainerClass particles;
        particles.reserve(loader.getNumberOfParticles());

        FTic timer;

        for(FSize idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
            FPoint<FReal> particlePosition;
            FReal physicalValue;
            loader.fillParticle(&particlePosition,&physicalValue);
            particles.push(particlePosition, physicalValue );
        }
        std::cout << "Load the file in " << timer.tacAndElapsed() << "s\n";

        timer.tic();
        {
            const FReal*const partX = particles.getPositions()[0];
            const FReal*const partY = particles.getPositions()[1];
            const FReal*const partZ = particles.getPositions()[2];
            const FReal*const physicalValues = particles.getPhysicalValues();

            OctreeClass tree(NbLevels, SizeSubLevels, loader.getBoxWidth(), loader.getCenterOfBox());

            for(FSize idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
                const FPoint<FReal> particlePosition(
                            partX[idxPart], partY[idxPart], partZ[idxPart]
                            );
                tree.insert(particlePosition, physicalValues[idxPart] );
            }
        }
        std::cout << "Create the tree in sequential " << timer.tacAndElapsed() << "s\n";

        timer.tic();
        {
            OctreeClass tree(NbLevels, SizeSubLevels, loader.getBoxWidth(), loader.getCenterOfBox());
            // This will modify particles!
            FTreeBuilder<FReal,OctreeClass, LeafClass>::BuildTreeFromArray(&tree, particles);
        }
        std::cout << "Create the tree in parallel " << timer.tacAndElapsed() << "s\n";

    }

    // -----------------------------------------------------

    return 0;
}

