// See LICENCE file at project root

#include <iostream>

#include <cstdio>
#include <cstdlib>

#include "Utils/FTic.hpp"
#include "Utils/FParameters.hpp"

#include "Components/FTypedLeaf.hpp"

#include "Containers/FOctree.hpp"
#include "Containers/FVector.hpp"

#include "Core/FFmmAlgorithmTsm.hpp"

#include "Kernels/Spherical/FSphericalKernel.hpp"
#include "Kernels/Spherical/FSphericalCell.hpp"

#include "Kernels/Rotation/FRotationKernel.hpp"
#include "Kernels/Rotation/FRotationCell.hpp"

#include "Files/FFmaTsmLoader.hpp"

#include "Kernels/P2P/FP2PParticleContainer.hpp"

#include "Utils/FParameterNames.hpp"

/** This program show an example of use of
  * the fmm basic algo
  * it also check that eachh particles is little or longer
  * related that each other
  */


// Simply create particles and try the kernels
template <class FReal, class CellClass, class ContainerClass, class LeafClass, class OctreeClass,
          class KernelClass, class FmmClass, typename... Args>
int testFunction(int argc, char ** argv, Args ... kernelPreArgs){
    FTic counter;
    // Retrieve parameters
    const int NbLevels = FParameters::getValue(argc,argv,FParameterDefinitions::OctreeHeight.options, 5);
    const int SizeSubLevels = FParameters::getValue(argc,argv,FParameterDefinitions::OctreeSubHeight.options, 3);
    // Get working file
    const char* const filename = FParameters::getStr(argc,argv,FParameterDefinitions::InputFile.options, "../Data/test20k.tsm.fma");
    std::cout << "Opening : " << filename << "\n";
    // Create particles loader
    FFmaTsmLoader<FReal> loader(filename);
    if(!loader.isOpen()){
        std::cout << "Loader Error, " << filename << " is missing\n";
        return 1;
    }

    // -----------------------------------------------------
    // Build the tree
    OctreeClass tree(NbLevels, SizeSubLevels,loader.getBoxWidth(),loader.getCenterOfBox());

    // -----------------------------------------------------

    std::cout << "Creating & Inserting " << loader.getNumberOfParticles() << " particles ..." << std::endl;
    std::cout << "\tHeight : " << NbLevels << " \t sub-height : " << SizeSubLevels << std::endl;
    counter.tic();

    for(FSize idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
        FPoint<FReal> particlePosition;
        FReal physicalValue = 0.0;
        FParticleType particleType;
        loader.fillParticle(&particlePosition,&physicalValue,&particleType);
        tree.insert(particlePosition, particleType, physicalValue );
    }

    counter.tac();
    std::cout << "Done  " << "(@Creating and Inserting Particles = " << counter.elapsed() << "s)." << std::endl;

    // -----------------------------------------------------

    std::cout << "Create kernel ..." << std::endl;
    counter.tic();
    //    KernelClass kernels( kernelPreArgs... , NbLevels, loader.getBoxWidth(), loader.getCenterOfBox());
    KernelClass kernels( kernelPreArgs... , loader.getBoxWidth(), loader.getCenterOfBox());

    counter.tac();
    std::cout << "Done  " << " in " << counter.elapsed() << "s)." << std::endl;

    // -----------------------------------------------------

    std::cout << "Working on particles ..." << std::endl;

    FmmClass algo(&tree,&kernels);

    counter.tic();
    algo.execute();
    counter.tac();

    std::cout << "Done  " << "(@Algorithm = " << counter.elapsed() << "s)." << std::endl;

    { // get sum forces&potential
        FReal potential = 0;
        FReal fx = 0.0, fy = 0.0, fz = 0.0;

        tree.forEachLeaf([&](LeafClass* leaf){
            const FReal*const potentials = leaf->getTargets()->getPotentials();
            const FReal*const forcesX = leaf->getTargets()->getForcesX();
            const FReal*const forcesY = leaf->getTargets()->getForcesY();
            const FReal*const forcesZ = leaf->getTargets()->getForcesZ();
            const FSize nbParticlesInLeaf = leaf->getTargets()->getNbParticles();

            for(FSize idxPart = 0 ; idxPart < nbParticlesInLeaf ; ++idxPart){
                potential += potentials[idxPart];
                fx += forcesX[idxPart];
                fy += forcesY[idxPart];
                fz += forcesZ[idxPart];
            }
        });

        std::cout << "Foces Sum  x = " << fx << " y = " << fy << " z = " << fz << std::endl;
        std::cout << "Potential = " << potential << std::endl;
    }


    // -----------------------------------------------------

    return 0;
}

// This is the real main!
int main(int argc, char ** argv){
    FHelpDescribeAndExit(argc, argv,
                         "Test the TSM (target source model) using the Rotation or the Spherical Harmonic Old implementations.",
                         FParameterDefinitions::OctreeHeight,FParameterDefinitions::SHDevelopment,
                         FParameterDefinitions::OctreeSubHeight, FParameterDefinitions::InputFile,
                         FParameterDefinitions::SphericalKernel, FParameterDefinitions::RotationKernel);

    std::cout << "[PARAM] Use Parameters -spherical -rotation -chebyshev\n";
    const int NbLevels = FParameters::getValue(argc,argv,FParameterDefinitions::OctreeHeight.options, 5);

    if( FParameters::existParameter(argc,argv,FParameterDefinitions::SphericalKernel.options) ){
        std::cout << "[INFO] -spherical is used\n";
        // Create template
        typedef double FReal;
        typedef FTypedSphericalCell< FReal>            CellClass;
        typedef FP2PParticleContainer<FReal>         ContainerClass;

        typedef FTypedLeaf< FReal, ContainerClass >                      LeafClass;
        typedef FOctree<FReal, CellClass, ContainerClass , LeafClass >  OctreeClass;
        typedef FSphericalKernel< FReal, CellClass, ContainerClass >          KernelClass;

        typedef FFmmAlgorithmTsm<OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass > FmmClass;

        const int DevP = FParameters::getValue(argc,argv,FParameterDefinitions::SHDevelopment.options, 8);
        CellClass::Init(DevP);

        // Call Main function
        testFunction< FReal, CellClass, ContainerClass, LeafClass, OctreeClass, KernelClass, FmmClass>(argc, argv, DevP,NbLevels);
    }

    if( FParameters::existParameter(argc,argv,FParameterDefinitions::RotationKernel.options) ){
        std::cout << "[INFO] -rotation is used\n";
        // Create template
        typedef double FReal;
        static const int P = 9;
        typedef FTypedRotationCell<FReal,P>            CellClass;
        typedef FP2PParticleContainer<FReal>         ContainerClass;

        typedef FTypedLeaf< FReal, ContainerClass >                      LeafClass;
        typedef FOctree<FReal, CellClass, ContainerClass , LeafClass >  OctreeClass;
        typedef FRotationKernel< FReal, CellClass, ContainerClass, P >          KernelClass;

        typedef FFmmAlgorithmTsm<OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass > FmmClass;

        // Call Main function
        testFunction< FReal, CellClass, ContainerClass, LeafClass, OctreeClass, KernelClass, FmmClass>(argc, argv,NbLevels);
    }

    return 0;
}
