// ===================================================================================
// Copyright ScalFmm 2015 INRIA, Olivier Coulaud, Berenger Bramas
// olivier.coulaud@inria.fr, berenger.bramas@inria.fr
// This software is a computer program whose purpose is to compute the FMM.
//
// This software is governed by the CeCILL-C and LGPL licenses and
// abiding by the rules of distribution of free software.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public and CeCILL-C Licenses for more details.
// "http://www.cecill.info".
// "http://www.gnu.org/licenses".
// ===================================================================================

// ==== CMAKE =====
// @FUSE_FFT
// @FUSE_BLAS
//  ==== Git =====
//
// ================

/** \brief Uniform FMM example
 *
 * \file
 * \authors B. Bramas, O. Coulaud
 *
 * This program runs the FMM Algorithm with the interpolation kernel based on
 * uniform (grid points) interpolation (1/r kernel). It then compares the
 * results with a direct computation.
 */


#include <iostream>
#include <iomanip>

#include <cstdio>
#include <cstdlib>
#include <memory>
#include <string>

#ifdef _OPENMP
#include <omp.h>
#endif
#include "ScalFmmConfig.h"

#include "Files/FFmaGenericLoader.hpp"

#include "Kernels/Interpolation/FInterpMatrixKernel.hpp"

#include "Components/FSimpleLeaf.hpp"
#include "Kernels/P2P/FP2PParticleContainerIndexed.hpp"

#include "Utils/FParameters.hpp"
#include "Utils/FParameterNames.hpp"

#ifdef _OPENMP
#include "Adaptive/FAdaptiveTask.hpp"
#endif
#include "Adaptive/FAdaptiveSequential.hpp"

#ifdef SCALFMM_USE_STARPU
#include "Adaptive/FAdaptiveStarPU.hpp"
#endif

#include "Adaptive/FTree.hpp"


// Types definition

// accuracy
typedef double FReal;
constexpr unsigned int ORDER = 7;

//
// Specification
using MatrixKernelClass = FInterpMatrixKernelR<FReal>;
//
//
/// definition of the common tree structure
//
using ContainerClass = FP2PParticleContainerIndexed<FReal>;
using CellClass      = FInterpolationCell<FReal, ORDER>;
using OctreeClass    = FTree<ContainerClass,CellClass>;
using LeafClass      = typename OctreeClass::node_t;

using KernelClass    = FInterpolationAdaptiveKernel<FReal,CellClass,ContainerClass,MatrixKernelClass,ORDER> ;

using FmmClass       = FAdaptiveTask<OctreeClass, KernelClass>;
using FmmClassSeq    = FAdaptiveSequential<OctreeClass, KernelClass>;

#ifdef SCALFMM_USE_STARPU
using FmmClassStarPU = FAdaptiveStarPU<OctreeClass, KernelClass>;
#endif

namespace FPD = FParameterDefinitions;
namespace FParameterDefinitions {
    FParameterNames Density = {
        {"--max-density", "-d"}, "Maximun particle count per leaf."
    };
    FParameterNames PrintParticles = {
        {"--print-particles"}, "Print particles in tree."
    };
    FParameterNames UseSequential = {
        {"--seq"}, "Use sequential algorithm"
    };
    FParameterNames UseStarPU = {
        {"--starpu"}, "Use StarPU algorithm"
    };
}

int main(int argc, char* argv[]) {
    FHelpDescribeAndExit(
        argc, argv,
        "Driver for Lagrange interpolation kernel  (1/r kernel).",
        FPD::InputFile,
        FPD::OutputFile,
        FPD::NbThreads,
        FPD::Density,
        FPD::UseSequential,
        FPD::UseStarPU,
        );


    const std::string defaultFile(SCALFMMDataPath+"unitCubeXYZQ100.bfma" );

    const std::string filename =
        FParameters::getStr(argc, argv, FPD::InputFile.options, defaultFile.c_str());

    const unsigned int NbThreads =
        FParameters::getValue(argc, argv, FPD::NbThreads.options, omp_get_max_threads());
    omp_set_num_threads(NbThreads);

    const unsigned int maxDensity =
        FParameters::getValue(argc, argv, FPD::Density.options, 1U);

    {
        std::string indent("    ");
        auto w = std::setw(18);
        std::cout << "Parameters" << std::endl << std::left
                  << indent << w << "Input file   : " << filename      << std::endl
                  << indent << w << "Thread number: "    << NbThreads     << std::endl
                  << std::endl;
    }
    //
    // init timer
    FTic time;

    // open particle file
    FFmaGenericLoader<FReal> loader(filename);

    const MatrixKernelClass MatrixKernel;

    // init oct-tree
    OctreeClass tree({loader.getBoxWidth(), loader.getCenterOfBox()});
    tree.leaf_max_particle_count(maxDensity);
    if(FParameters::existParameter(argc, argv, FPD::PrintParticles.options)) {
        tree.print_particles = true;
    }

    { // -----------------------------------------------------
        std::cout << "Creating & Inserting " << loader.getNumberOfParticles()
                  << " particles ..." << std::endl;
        time.tic();
        //
        FPoint<FReal> position;
        FReal physicalValue = 0.0;
        //
        for(FSize idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
            // Read particle per particle from file
            loader.fillParticle(&position,&physicalValue);
            // put particle in octree
            tree.insert(position, idxPart, physicalValue);
        }

        time.tac();
        std::cout << "Done  " << "(@Creating and Inserting Particles = "
                  << time.elapsed() << " s)." << std::endl;
        std::cout << "Tree height: " << tree.height() << std::endl;

        std::fstream("tree.txt", std::ios::out) << tree;
    } // -----------------------------------------------------

    { // -----------------------------------------------------
        std::cout << "\n" << interpolationType << "  FMM (ORDER= "<< ORDER << ") ... " << std::endl;
        time.tic();
        //
        //

        KernelClass kernel(
            static_cast<int>(tree.height()),
            loader.getBoxWidth(),
            loader.getCenterOfBox(),
            &MatrixKernel);


        unsigned int operations = FFmmNearAndFarFields;

        if(FParameters::existParameter(argc, argv, FPD::UseSequential.options)) {
            FmmClassSeq algo(&tree, &kernel);
            algo.execute(operations);   // FMM algorithm call
#ifdef SCALFMM_USE_STARPU
        } else if(FParameters::existParameter(argc, argv, FPD::UseStarPU.options)) {
            FmmClassStarPU algo(&tree, &kernel);
            algo.execute(operations);   // FMM algorithm call
#endif
        } else {
#ifdef SCALFMM_USE_OMP4
            FmmClass algo(&tree, &kernel);
            algo.execute(operations);   // FMM algorithm call
#endif
        }


        //
        //
        time.tac();
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
    if(FParameters::existParameter(argc, argv, FPD::OutputFile.options)){
        std::string name(FParameters::getStr(argc,argv,FPD::OutputFile.options, "output.fma"));
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
                particles[j]    = posX[idxPart] ;
                particles[j+1]  = posY[idxPart] ;
                particles[j+2]  = posZ[idxPart] ;
                particles[j+3]  = physicalValues[idxPart] ;
                particles[j+4]  = potentials[idxPart] ;
                particles[j+5]  = forcesX[idxPart] ;
                particles[j+6]  = forcesY[idxPart] ;
                particles[j+7]  = forcesZ[idxPart] ;
            }
        });

        writer.writeHeader( loader.getCenterOfBox(), loader.getBoxWidth() ,  NbPoints, sizeof(FReal), 8) ;
        writer.writeArrayOfReal(particles,  8 , NbPoints);

        delete[] particles;

        //
//        std::string name1( "output.fma");
//
        // FFmaGenericWriter<FReal> writer1(name1) ;
        // writer1.writeDistributionOfParticlesFromOctree(&tree,NbPoints) ;
    }


    return 0;
}
