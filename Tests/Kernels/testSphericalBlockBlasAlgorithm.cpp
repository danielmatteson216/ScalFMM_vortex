// See LICENCE file at project root

// ==== CMAKE =====
// @FUSE_BLAS
// ================

#include <iostream>

#include <cstdio>
#include <cstdlib>

#include "Utils/FTic.hpp"
#include "Utils/FParameters.hpp"

#include "Containers/FOctree.hpp"
#include "Containers/FVector.hpp"

#include "Core/FFmmAlgorithm.hpp"
#include "Core/FFmmAlgorithmThread.hpp"
#include "Core/FFmmAlgorithmTask.hpp"

#include "Components/FSimpleLeaf.hpp"
#include "Components/FBasicCell.hpp"

#include "Kernels/Spherical/FSphericalBlockBlasKernel.hpp"
#include "Kernels/Spherical/FSphericalCell.hpp"

#include "Files/FFmaScanfLoader.hpp"

#include "Kernels/P2P/FP2PParticleContainer.hpp"

#include "Utils/FParameterNames.hpp"

/** This program show an example of use of
  * the fmm basic algo
  * it also check that eachh particles is little or longer
  * related that each other
  */



// Simply create particles and try the kernels
int main(int argc, char ** argv){
    FHelpDescribeAndExit(argc, argv,
                         "Run a Spherical Harmonic (Block BLAS Implementation) FMM kernel and compare the accuracy with a direct computation.",
                         FParameterDefinitions::InputFile, FParameterDefinitions::OctreeHeight,
                         FParameterDefinitions::OctreeSubHeight, FParameterDefinitions::SequentialFmm,
                         FParameterDefinitions::TaskFmm, FParameterDefinitions::SHDevelopment);

    typedef double FReal;
    typedef FSphericalCell<FReal>                 CellClass;
    typedef FP2PParticleContainer<FReal>        ContainerClass;

    typedef FSimpleLeaf<FReal, ContainerClass >                     LeafClass;
    typedef FOctree< FReal,CellClass, ContainerClass , LeafClass >  OctreeClass;
    typedef FSphericalBlockBlasKernel<FReal, CellClass, ContainerClass > KernelClass;

    typedef FFmmAlgorithm<OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass > FmmClass;
    typedef FFmmAlgorithmThread<OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass > FmmClassThread;
    typedef FFmmAlgorithmTask<OctreeClass,  CellClass, ContainerClass, KernelClass, LeafClass > FmmClassTask;
    ///////////////////////What we do/////////////////////////////
    std::cout << ">> This executable has to be used to test Spherical Block Blas algorithm.\n";
    std::cout << ">> You can pass -sequential or -task (thread by default).\n";
    //////////////////////////////////////////////////////////////
    const int DevP = FParameters::getValue(argc,argv,FParameterDefinitions::SHDevelopment.options, 8);
    const int NbLevels = FParameters::getValue(argc,argv,FParameterDefinitions::OctreeHeight.options, 5);
    const int SizeSubLevels = FParameters::getValue(argc,argv,FParameterDefinitions::OctreeSubHeight.options, 3);
    FTic counter;

    const char* const filename = FParameters::getStr(argc,argv,FParameterDefinitions::InputFile.options, "../Data/test20k.fma");
    std::cout << "Opening : " << filename << "\n";

    FFmaScanfLoader<FReal> loader(filename);
    if(!loader.isOpen()){
        std::cout << "Loader Error, " << filename << " is missing\n";
        return 1;
    }

    // -----------------------------------------------------
    CellClass::Init(DevP, true);
    OctreeClass tree(NbLevels, SizeSubLevels,loader.getBoxWidth(),loader.getCenterOfBox());

    // -----------------------------------------------------

    std::cout << "Creating & Inserting " << loader.getNumberOfParticles() << " particles ..." << std::endl;
    std::cout << "\tHeight : " << NbLevels << " \t sub-height : " << SizeSubLevels << std::endl;
    counter.tic();

    for(FSize idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
        FPoint<FReal> particlePosition;
        FReal physicalValue = 0.0;
        loader.fillParticle(&particlePosition,&physicalValue);
        tree.insert(particlePosition, physicalValue );
    }

    counter.tac();
    std::cout << "Done  " << "(@Creating and Inserting Particles = " << counter.elapsed() << "s)." << std::endl;

    // -----------------------------------------------------

    std::cout << "Create kernel ..." << std::endl;
    counter.tic();

    KernelClass kernels(DevP, NbLevels, loader.getBoxWidth(), loader.getCenterOfBox());

    counter.tac();
    std::cout << "Done  " << " in " << counter.elapsed() << "s)." << std::endl;

    // -----------------------------------------------------

    std::cout << "Working on particles ..." << std::endl;

    if( FParameters::findParameter(argc,argv,FParameterDefinitions::SequentialFmm.options) != FParameters::NotFound){
        FmmClass algo(&tree,&kernels);
        counter.tic();
        algo.execute();
    }
    else if( FParameters::findParameter(argc,argv,FParameterDefinitions::TaskFmm.options) != FParameters::NotFound){
        FmmClassTask algo(&tree,&kernels);
        counter.tic();
        algo.execute();
    }
    else {
        FmmClassThread algo(&tree,&kernels);
        counter.tic();
        algo.execute();
    }

    counter.tac();
    std::cout << "Done  " << "(@Algorithm = " << counter.elapsed() << "s)." << std::endl;

    { // get sum forces&potential
        FReal potential = 0;
        FReal fx = 0.0, fy = 0.0, fz = 0.0;

        tree.forEachLeaf([&](LeafClass* leaf){
            const FReal*const physicalValues = leaf->getTargets()->getPhysicalValues();
            const FReal*const potentials = leaf->getTargets()->getPotentials();
            const FReal*const forcesX = leaf->getTargets()->getForcesX();
            const FReal*const forcesY = leaf->getTargets()->getForcesY();
            const FReal*const forcesZ = leaf->getTargets()->getForcesZ();
            const FSize nbParticlesInLeaf = leaf->getTargets()->getNbParticles();

            for(FSize idxPart = 0 ; idxPart < nbParticlesInLeaf ; ++idxPart){
                potential += potentials[idxPart] * physicalValues[idxPart];
                fx += forcesX[idxPart];
                fy += forcesY[idxPart];
                fz += forcesZ[idxPart];
            }
        });

        std::cout << "Foces Sum  x = " << fx << " y = " << fy << " z = " << fz << std::endl;
        std::cout << "Potential = " << potential << std::endl;
    }

    return 0;
}



