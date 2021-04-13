// See LICENCE file at project root

// ==== CMAKE =====
// @FUSE_FFT
// @FUSE_BLAS
//  ==== Git =====
// ================

/**
 * \brief ScalFMM P2P-only algorithm for the 1/r kernel example.
 * \authors O. Coulaud
 *
 * This program runs the only the P2P Algorithm with the cutoff kernel (1/r) and
 * compares the results with a direct computation.
 */

#include <iostream>

#include <cstdio>
#include <cstdlib>
#include <string>

#include "ScalFmmConfig.h"
#include "Utils/FGlobal.hpp"

#include "Utils/FParameters.hpp"
#include "Utils/FParameterNames.hpp"
//  Our particles loader
#include "Files/FFmaGenericLoader.hpp"
//
// Cells without local and multipole values
#include "Kernels/Generic/FGenericData.hpp"
// Leaves
#include "Components/FSimpleLeaf.hpp"
#include "Kernels/P2P/FP2PParticleContainerIndexed.hpp"
#include "Kernels/Interpolation/FInterpMatrixKernel.hpp"  // For the predefined Matrix kernels
#include "Kernels/Interpolation/FCutOffKernel.hpp"        // Generic CutOff kernel defined with Matrix Kernel
// The classical Octree
#include "Containers/FOctree.hpp"

#ifdef _OPENMP
#include "Core/FFmmAlgorithmThread.hpp"
#else
#include "Core/FFmmAlgorithm.hpp"
#endif

#include <memory>

int main(int argc, char* argv[])
{
    FHelpDescribeAndExit(argc, argv,
                         "Driver for CutOff  kernel  (1/r kernel).",
                         FParameterDefinitions::InputFile, FParameterDefinitions::OctreeHeight,
                         FParameterDefinitions::OctreeSubHeight, FParameterDefinitions::OutputFile,
                         FParameterDefinitions::NbThreads);


    const std::string defaultFile(SCALFMMDataPath+"unitCubeXYZQ100.bfma" );
    const std::string filename       = FParameters::getStr(argc,argv,FParameterDefinitions::InputFile.options, defaultFile.c_str());
    const unsigned int TreeHeight    = FParameters::getValue(argc, argv, FParameterDefinitions::OctreeHeight.options, 2);
    const unsigned int SubTreeHeight = FParameters::getValue(argc, argv, FParameterDefinitions::OctreeSubHeight.options, 1);
    const unsigned int NbThreads     = FParameters::getValue(argc, argv, FParameterDefinitions::NbThreads.options, 4);


    #ifdef _OPENMP
    omp_set_num_threads(NbThreads);
    std::cout << "\n>> Using " << omp_get_max_threads() << " threads.\n" << std::endl;
    #else
    std::cout << "\n>> Sequential version.\n" << std::endl;
    #endif
    //
    std::cout << "Parameters  "<< std::endl
              << "      Octree Depth      " << TreeHeight    << std::endl
              << "      SubOctree depth   " << SubTreeHeight << std::endl
              << "      Input file  name: " << filename      << std::endl
              << "      Thread number:    " << NbThreads     << std::endl
              << std::endl;
    //
    // init timer
    FTic time;

    // open particle file
    ////////////////////////////////////////////////////////////////////
    //
    typedef double FReal;
    FFmaGenericLoader<FReal> loader(filename);
    //
    ////////////////////////////////////////////////////////////////////
    // begin CutOff kernel
    // accuracy
    // typedefs
    typedef FP2PParticleContainerIndexed<FReal>                 ContainerClass;
    typedef FSimpleLeaf<FReal,  ContainerClass >                LeafClass;
    typedef FGenericData<char,char>                             CellClass;
    typedef FOctree<FReal, CellClass,ContainerClass,LeafClass>  OctreeClass;
    //
    typedef FInterpMatrixKernelR<FReal>                                          MatrixKernelClass;
    const MatrixKernelClass MatrixKernel;
    typedef FCutOffKernel<FReal,CellClass,ContainerClass,MatrixKernelClass>  KernelClass;
    //
    #ifdef _OPENMP
    typedef FFmmAlgorithmThread<OctreeClass,CellClass,ContainerClass,KernelClass,LeafClass> FmmClass;
    #else
    typedef FFmmAlgorithm<OctreeClass,CellClass,ContainerClass,KernelClass,LeafClass> FmmClass;
    #endif
    // init oct-tree
    OctreeClass tree(TreeHeight, SubTreeHeight, loader.getBoxWidth(), loader.getCenterOfBox());


    { // -----------------------------------------------------
        std::cout << "Creating & Inserting " << loader.getNumberOfParticles()
                  << " particles ..." << std::endl;
        std::cout << "\tHeight : " << TreeHeight << " \t sub-height : " << SubTreeHeight << std::endl;
        time.tic();
        //
        FPoint<FReal> position;
        FReal physicalValue = 0.0;
        //
        for(FSize idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
            //
            // Read particle per particle from file
            loader.fillParticle(&position,&physicalValue);
            //
            // put particle in octree
            tree.insert(position, idxPart, physicalValue);
        }

        time.tac();
        std::cout << "Done  " << "(@Creating and Inserting Particles = "
                  << time.elapsed() << " s) ." << std::endl;
    } // -----------------------------------------------------

    { // -----------------------------------------------------
        std::cout << " CutOff Algorithm  ... " << std::endl;

        time.tic();
        //
        std::unique_ptr<KernelClass> kernels(new KernelClass(TreeHeight, loader.getBoxWidth(), loader.getCenterOfBox(),&MatrixKernel));
        //
        FmmClass algo(&tree, kernels.get());
        //
        algo.execute( FFmmP2P);   // Here the call of the FMM algorithm
        //
        time.tac();
        std::cout << "Timers Far Field \n"
                  << "P2M " << algo.getTime(FAlgorithmTimers::P2MTimer) << " seconds\n"
                  << "M2M " << algo.getTime(FAlgorithmTimers::M2MTimer) << " seconds\n"
                  << "M2L " << algo.getTime(FAlgorithmTimers::M2LTimer) << " seconds\n"
                  << "L2L " << algo.getTime(FAlgorithmTimers::L2LTimer) << " seconds\n"
                  << "P2P and L2P " << algo.getTime(FAlgorithmTimers::NearTimer) << " seconds\n"
                  << std::endl;


        std::cout << "Done  " << "(@Algorithm = " << time.elapsed() << " s) ." << std::endl;
    }
    // -----------------------------------------------------
    //
    // Some output
    //
    //
    { // -----------------------------------------------------
        FSize N1=0, N2= loader.getNumberOfParticles()/2, N3= loader.getNumberOfParticles() -1; ;
        FReal energy =0.0 ;
        //
        //   Loop over all leaves
        //
        std::cout <<std::endl<<" &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& "<<std::endl;
        std::cout << std::scientific;
        std::cout.precision(10) ;

        tree.forEachLeaf([&](LeafClass* leaf){
                const FReal*const posX = leaf->getTargets()->getPositions()[0];
                const FReal*const posY = leaf->getTargets()->getPositions()[1];
                const FReal*const posZ = leaf->getTargets()->getPositions()[2];

                const FReal*const potentials = leaf->getTargets()->getPotentials();
                const FReal*const forcesX = leaf->getTargets()->getForcesX();
                const FReal*const forcesY = leaf->getTargets()->getForcesY();
                const FReal*const forcesZ = leaf->getTargets()->getForcesZ();
                const FSize nbParticlesInLeaf = leaf->getTargets()->getNbParticles();
                const FReal*const physicalValues = leaf->getTargets()->getPhysicalValues();

                const FVector<FSize>& indexes = leaf->getTargets()->getIndexes();

                for(FSize idxPart = 0 ; idxPart < nbParticlesInLeaf ; ++idxPart){
                    const FSize indexPartOrig = indexes[idxPart];
                    if ((indexPartOrig == N1) || (indexPartOrig == N2) || (indexPartOrig == N3)  ) {
                        std::cout << "Index "<< indexPartOrig <<"  potential  " << potentials[idxPart]
                                  << " Pos "<<posX[idxPart]<<" "<<posY[idxPart]<<" "<<posZ[idxPart]
                                  << "   Forces: " << forcesX[idxPart] << " " << forcesY[idxPart] << " "<< forcesZ[idxPart] <<std::endl;
                    }
                    energy += potentials[idxPart]*physicalValues[idxPart] ;
                }
            });
        std::cout <<std::endl<<"Energy: "<< energy<<std::endl;
        std::cout <<std::endl<<" &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& "<<std::endl<<std::endl;

    }
    // -----------------------------------------------------
    if(FParameters::existParameter(argc, argv, FParameterDefinitions::OutputFile.options)){
        std::string name(FParameters::getStr(argc,argv,FParameterDefinitions::OutputFile.options,   "output.fma"));
        FFmaGenericWriter<FReal> writer(name) ;
        //
        FSize NbPoints = loader.getNumberOfParticles();
        FReal * particles ;
        particles = new FReal[8*NbPoints] ;
        memset(particles,0,8*NbPoints*sizeof(FReal));
        FSize j = 0 ;
        tree.forEachLeaf([&](LeafClass* leaf){
                //
                // Input
                const FReal*const posX = leaf->getTargets()->getPositions()[0];
                const FReal*const posY = leaf->getTargets()->getPositions()[1];
                const FReal*const posZ = leaf->getTargets()->getPositions()[2];
                const FReal*const physicalValues = leaf->getTargets()->getPhysicalValues();
                const FVector<FSize>& indexes = leaf->getTargets()->getIndexes();
                //
                // Computed data
                const FReal*const potentials = leaf->getTargets()->getPotentials();
                const FReal*const forcesX = leaf->getTargets()->getForcesX();
                const FReal*const forcesY = leaf->getTargets()->getForcesY();
                const FReal*const forcesZ = leaf->getTargets()->getForcesZ();
                //
                const FSize nbParticlesInLeaf = leaf->getTargets()->getNbParticles();
                for(FSize idxPart = 0 ; idxPart < nbParticlesInLeaf ; ++idxPart){
                    j = 8*indexes[idxPart];
                    particles[j]      = posX[idxPart] ;
                    particles[j+1]  = posY[idxPart] ;
                    particles[j+2]  = posZ[idxPart] ;
                    particles[j+3]  = physicalValues[idxPart] ;
                    particles[j+4]  = potentials[idxPart] ;
                    particles[j+5]  =  forcesX[idxPart] ;
                    particles[j+6]  =  forcesY[idxPart] ;
                    particles[j+7]  =  forcesZ[idxPart] ;
                }
            });

        writer.writeHeader( loader.getCenterOfBox(), loader.getBoxWidth() ,  NbPoints, sizeof(FReal), 8) ;
        writer.writeArrayOfReal(particles,  8 , NbPoints);

        delete[] particles;

        //
        std::string name1( "output.fma");
        //
        FFmaGenericWriter<FReal> writer1(name1) ;
        writer1.writeDistributionOfParticlesFromOctree(&tree,NbPoints) ;
    }


    return 0;
}
