// See LICENCE file at project root
//
// ==== CMAKE =====
//
// ================

// Keep in private GIT
//


#include <iostream>
#include <stdexcept>
#include <string>
#include <cstdlib>
#include <cstdio>
//
#include "ScalFmmConfig.h"
#include "Utils/FTic.hpp"
#include "Utils/FParameters.hpp"

#include "Files/FFmaGenericLoader.hpp"

#include "Containers/FOctree.hpp"
#include "Containers/FVector.hpp"

#include "Core/FFmmAlgorithmThread.hpp"

#ifdef SCALFMM_USE_BLAS
// chebyshev kernel

#include "Kernels/Chebyshev/FChebCell.hpp"
#include "Kernels/Interpolation/FInterpMatrixKernel.hpp"
#include "Kernels/Chebyshev/FChebKernel.hpp"
#include "Kernels/Chebyshev/FChebSymKernel.hpp"
#endif
//
// spherical kernel
#include "Kernels/Spherical/FSphericalCell.hpp"
#ifdef SCALFMM_USE_BLAS
#include "Kernels/Spherical/FSphericalBlasKernel.hpp"
#include "Kernels/Spherical/FSphericalBlockBlasKernel.hpp"
#endif
//
// taylor kernel
#include "Kernels/Taylor/FTaylorCell.hpp"
#include "Kernels/Taylor/FTaylorKernel.hpp"
//
#include "Components/FSimpleLeaf.hpp"
#include "Kernels/P2P/FP2PParticleContainerIndexed.hpp"

//Classical Spherical kernel
#include "Kernels/Spherical/FSphericalCell.hpp"
#include "Kernels/Spherical/FSphericalKernel.hpp"

//Rotation kernel
#include "Kernels/Rotation/FRotationKernel.hpp"
#include "Kernels/Rotation/FRotationCell.hpp"

#ifdef SCALFMM_USE_FFT
// Uniform grid kernel
#include "Kernels/Uniform/FUnifCell.hpp"
#include "Kernels/Interpolation/FInterpMatrixKernel.hpp"

#include "Kernels/Uniform/FUnifKernel.hpp"
#endif

#include "Utils/FParameterNames.hpp"

/**
 * This program compares two different kernels, eg., the Chebyshev kernel with
 * the SphericalBlas kernel.
 */
//
/// \file  compareAllPoissonKernels.cpp
//!
//! \brief compareAllPoissonKernels: Driver to compare all different implementations for 1/r kernel.
//!
//! \details  This driver  compare all different implementations provided by our library for the classical Poisson kernel (1/r)<br>
//!     We check the following kernels <br>
//!      - Spherical expansion; Spherical with rotation optimization
//!      - Taylor expansion; Spherical with rotation optimization
//!      - if BLAS is activated
//!               - Blas And Block Blas optiization of the M2L operator
//!               - Chebychev and symetric Chebychev interpolation
//!      - if FFT is activated:  interpolation on uniform grid
//!<br>


// Simply create particles and try the kernels
int main(int argc, char* argv[])
{
    FHelpDescribeAndExit(argc, argv,
                         "Driver for testing different approximations  for the  1/r kernel.",
                         FParameterDefinitions::InputFile, FParameterDefinitions::OctreeHeight,
                         FParameterDefinitions::OctreeSubHeight, FParameterDefinitions::InputFile,
                         FParameterDefinitions::NbThreads, FParameterDefinitions::SHDevelopment);

    typedef double FReal;
        // get info from commande line
    const std::string  filename(FParameters::getStr(argc,argv,FParameterDefinitions::InputFile.options, "../Data/UTest/unitCubeRef20kDouble.bfma"));
    const unsigned int TreeHeight    = FParameters::getValue(argc, argv, FParameterDefinitions::OctreeHeight.options, 5);
    const unsigned int SubTreeHeight = FParameters::getValue(argc, argv, FParameterDefinitions::OctreeSubHeight.options, 2);
    const unsigned int NbThreads     = FParameters::getValue(argc, argv, FParameterDefinitions::NbThreads.options, omp_get_max_threads());
    const int DevP                   = FParameters::getValue(argc, argv, FParameterDefinitions::SHDevelopment.options, 11);

        //
#ifdef _OPENMP
        omp_set_num_threads(NbThreads);
#else
    std::cout << "\n>> Sequential version.\n" << std::endl;
#endif

        std::cout <<     "Parameters  "<< std::endl
                        <<     "      Octree Depth      \t"<< TreeHeight <<std::endl
                        <<        "      SubOctree depth \t"<< SubTreeHeight <<std::endl
                        <<     "      Input file  name: \t" <<filename <<std::endl
                        <<     "      Thread number:  \t" << NbThreads <<std::endl
                        <<std::endl;

        // init timer
        FTic time;

    FFmaGenericLoader<FReal> loader(filename);
        //  if(!loader.isOpen()) throw std::runtime_error("Particle file couldn't be opened!");
        //

        FSize nbParticles = loader.getNumberOfParticles() ;
    FmaRWParticle<FReal, 8,8>* const particles = new FmaRWParticle<FReal, 8,8>[nbParticles];
        //
        loader.fillParticle(particles,nbParticles);
        //
        ////////////////////////////////////////////////////////////////////
        //  Compute direct energy
        FReal energyD =0.0, totPhysicalValue =0.0;

#pragma omp parallel for reduction(+:energyD,totPhysicalValue)
    for(FSize idx = 0 ; idx < loader.getNumberOfParticles()  ; ++idx){
                energyD             +=  particles[idx].getPotential()*particles[idx].getPhysicalValue() ;
                totPhysicalValue += particles[idx].getPhysicalValue() ;
        }
        std::cout << " Total Physical values: "<< totPhysicalValue <<std::endl;
        std::cout << " Energy of the system: "<< energyD <<std::endl;
        ////////////////////////////////////////////////////////////////////

#ifdef  SCALFMM_USE_BLAS
        {	// begin Chebyshev kernel

                // accuracy
                const unsigned int ORDER = 7;
                std::cout << "\nFChebKernel FMM ... ORDER: " << ORDER <<std::endl;

                // typedefs
        typedef FP2PParticleContainerIndexed<FReal> ContainerClass;
        typedef FSimpleLeaf<FReal, ContainerClass> LeafClass;
        typedef FInterpMatrixKernelR<FReal> MatrixKernelClass;
        typedef FChebCell<FReal, ORDER> CellClass;
        typedef FOctree<FReal, CellClass,ContainerClass,LeafClass> OctreeClass;

        typedef FChebKernel<FReal, CellClass,ContainerClass,MatrixKernelClass,ORDER> KernelClass;
                typedef FFmmAlgorithmThread<OctreeClass,CellClass,ContainerClass,KernelClass,LeafClass> FmmClass;


                // init oct-tree
                OctreeClass tree(TreeHeight, SubTreeHeight, loader.getBoxWidth(), loader.getCenterOfBox());

    // Create Matrix Kernel
    const MatrixKernelClass MatrixKernel;

                { // -----------------------------------------------------
                        time.tic();
            for(FSize idxPart = 0 ; idxPart < nbParticles; ++idxPart){
                                // put in tree
                                tree.insert(particles[idxPart].getPosition(), idxPart, particles[idxPart].getPhysicalValue());
                        }
                        time.tac();
                        std::cout << "(FChebKernel @ Inserting Particles = "<< time.elapsed() << " s)." << std::endl;
                } // -----------------------------------------------------

                { // -----------------------------------------------------
                        time.tic();
                        KernelClass kernels(TreeHeight, loader.getBoxWidth(), loader.getCenterOfBox(),&MatrixKernel);
                        FmmClass algorithm(&tree, &kernels);
                        algorithm.execute();
                        time.tac();
                        std::cout <<"(FChebKernel @Algorithm = " << time.elapsed() << " s)." << std::endl;
                } // -----------------------------------------------------
                FReal energy = 0.0;
        FMath::FAccurater<FReal> potentialDiff;
        FMath::FAccurater<FReal> fx, fy, fz;
                { // Check that each particle has been summed with all other

                        tree.forEachLeaf([&](LeafClass* leaf){
                                const FReal*const physicalValues = leaf->getTargets()->getPhysicalValues();
                                const FReal*const potentials        = leaf->getTargets()->getPotentials();
                                const FReal*const forcesX            = leaf->getTargets()->getForcesX();
                                const FReal*const forcesY            = leaf->getTargets()->getForcesY();
                                const FReal*const forcesZ            = leaf->getTargets()->getForcesZ();
                const FSize nbParticlesInLeaf           = leaf->getTargets()->getNbParticles();
                const FVector<FSize>& indexes       = leaf->getTargets()->getIndexes();

                for(FSize idxPart = 0 ; idxPart < nbParticlesInLeaf ; ++idxPart){
                    const FSize indexPartOrig = indexes[idxPart];
                                        potentialDiff.add(particles[indexPartOrig].getPotential(),potentials[idxPart]);
                                        fx.add(particles[indexPartOrig].getForces()[0],forcesX[idxPart]);
                                        fy.add(particles[indexPartOrig].getForces()[1],forcesY[idxPart]);
                                        fz.add(particles[indexPartOrig].getForces()[2],forcesZ[idxPart]);
                                        energy += potentials[idxPart]*physicalValues[idxPart] ;
                                }
                        });
                }

                // Print for information
                std::cout << "FChebKernel Energy "  << FMath::Abs(energy-energyD) /energyD << std::endl;
                std::cout << "FChebKernel Potential " << potentialDiff << std::endl;
                std::cout << "FChebKernel Fx " << fx << std::endl;
                std::cout << "FChebKernel Fy " << fy << std::endl;
                std::cout << "FChebKernel Fz " << fz << std::endl;

        } // end Chebyshev kernel
        {	// begin ChebSymKernel kernel

                // accuracy
                const unsigned int ORDER = 7;
                std::cout << "\nFChebSymKernel FMM ... ORDER: " << ORDER <<std::endl;

                // typedefs
        typedef FP2PParticleContainerIndexed<FReal> ContainerClass;
        typedef FSimpleLeaf<FReal, ContainerClass> LeafClass;
        typedef FInterpMatrixKernelR<FReal> MatrixKernelClass;
        typedef FChebCell<FReal, ORDER> CellClass;
        typedef FOctree<FReal, CellClass,ContainerClass,LeafClass> OctreeClass;

        typedef FChebSymKernel<FReal, CellClass,ContainerClass,MatrixKernelClass,ORDER> KernelClass;
                typedef FFmmAlgorithmThread<OctreeClass,CellClass,ContainerClass,KernelClass,LeafClass> FmmClass;


                // init oct-tree
                OctreeClass tree(TreeHeight, SubTreeHeight, loader.getBoxWidth(), loader.getCenterOfBox());

    // Create Matrix Kernel
    const MatrixKernelClass MatrixKernel;

                { // -----------------------------------------------------
                        time.tic();

            for(FSize idxPart = 0 ; idxPart < nbParticles; ++idxPart){
                                // put in tree
                                tree.insert(particles[idxPart].getPosition(), idxPart, particles[idxPart].getPhysicalValue());
                        }

                        time.tac();
                        std::cout <<  "(FChebSymKernel @Inserting Particles = "<< time.elapsed() << " s)." << std::endl;
                } // -----------------------------------------------------

                { // -----------------------------------------------------
                        time.tic();
                        KernelClass kernels(TreeHeight, loader.getBoxWidth(), loader.getCenterOfBox(),&MatrixKernel);
                        FmmClass algorithm(&tree, &kernels);
                        algorithm.execute();
                        time.tac();
                        std::cout << "(FChebSymKernel @Algorithm = " << time.elapsed() << " s)." << std::endl;
                } // -----------------------------------------------------
                FReal energy = 0.0;
        FMath::FAccurater<FReal> potentialDiff;
        FMath::FAccurater<FReal> fx, fy, fz;
                { // Check that each particle has been summed with all other

                        tree.forEachLeaf([&](LeafClass* leaf){
                                const FReal*const physicalValues = leaf->getTargets()->getPhysicalValues();
                                const FReal*const potentials        = leaf->getTargets()->getPotentials();
                                const FReal*const forcesX            = leaf->getTargets()->getForcesX();
                                const FReal*const forcesY            = leaf->getTargets()->getForcesY();
                                const FReal*const forcesZ            = leaf->getTargets()->getForcesZ();
                const FSize nbParticlesInLeaf           = leaf->getTargets()->getNbParticles();
                const FVector<FSize>& indexes       = leaf->getTargets()->getIndexes();

                for(FSize idxPart = 0 ; idxPart < nbParticlesInLeaf ; ++idxPart){
                    const FSize indexPartOrig = indexes[idxPart];
                                        potentialDiff.add(particles[indexPartOrig].getPotential(),potentials[idxPart]);
                                        fx.add(particles[indexPartOrig].getForces()[0],forcesX[idxPart]);
                                        fy.add(particles[indexPartOrig].getForces()[1],forcesY[idxPart]);
                                        fz.add(particles[indexPartOrig].getForces()[2],forcesZ[idxPart]);
                                        energy += potentials[idxPart]*physicalValues[idxPart] ;
                                }
                        });
                }

                // Print for information
                std::cout << "FChebSymKernel Energy "  << FMath::Abs(energy-energyD) /energyD << std::endl;
                std::cout << "FChebSymKernel Potential " << potentialDiff << std::endl;
                std::cout << "FChebSymKernel Fx " << fx << std::endl;
                std::cout << "FChebSymKernel Fy " << fy << std::endl;
                std::cout << "FChebSymKernel Fz " << fz << std::endl;

        } // end Chebyshev kernel
        //
        ////////////////////////////////////////////////////////////////////
        //
        {	// begin FFmaBlas kernel FSphericalBlockBlasKernel
                        std::cout << "\nFFmaBlas FMM ... P: " <<DevP << std::endl;

                        // typedefs
            typedef FSphericalCell<FReal>                 CellClass;
            typedef FP2PParticleContainerIndexed<FReal>         ContainerClass;
            typedef FSimpleLeaf<FReal,  ContainerClass >                     LeafClass;
            typedef FOctree<FReal,  CellClass, ContainerClass , LeafClass >  OctreeClass;
            typedef FSphericalBlasKernel<FReal,  CellClass, ContainerClass > KernelClass;
                        typedef FFmmAlgorithmThread<OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass > FmmClass;

                        // init cell class and oct-tree
                        CellClass::Init(DevP, true); // only for blas
                        OctreeClass tree(TreeHeight, SubTreeHeight, loader.getBoxWidth(), loader.getCenterOfBox());

                        { // -----------------------------------------------------
                                time.tic();

                for(FSize idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
                                        // put in tree
                                        tree.insert(particles[idxPart].getPosition(), idxPart, particles[idxPart].getPhysicalValue());
                                }

                                time.tac();
                                std::cout << "( FFmaBlas@Inserting Particles = " << time.elapsed() << " s)." << std::endl;
                        } // -----------------------------------------------------

                        // -----------------------------------------------------
                        time.tic();
                        KernelClass kernels(DevP, TreeHeight, loader.getBoxWidth(), loader.getCenterOfBox());
                        FmmClass algorithm(&tree, &kernels);
                        algorithm.execute();
                        time.tac();
                        std::cout << "(FFmaBlas @Algorithm = " << time.elapsed() << " s)." << std::endl;
                        // -----------------------------------------------------

                        FReal energy = 0.0;
            FMath::FAccurater<FReal> potentialDiff;
            FMath::FAccurater<FReal> fx, fy, fz;
                        { // Check that each particle has been summed with all other

                                tree.forEachLeaf([&](LeafClass* leaf){
                                        const FReal*const physicalValues = leaf->getTargets()->getPhysicalValues();
                                        const FReal*const potentials        = leaf->getTargets()->getPotentials();
                                        const FReal*const forcesX            = leaf->getTargets()->getForcesX();
                                        const FReal*const forcesY            = leaf->getTargets()->getForcesY();
                                        const FReal*const forcesZ            = leaf->getTargets()->getForcesZ();
                    const FSize nbParticlesInLeaf           = leaf->getTargets()->getNbParticles();
                    const FVector<FSize>& indexes       = leaf->getTargets()->getIndexes();

                    for(FSize idxPart = 0 ; idxPart < nbParticlesInLeaf ; ++idxPart){
                        const FSize indexPartOrig = indexes[idxPart];
                                                potentialDiff.add(particles[indexPartOrig].getPotential(),potentials[idxPart]);
                                                fx.add(particles[indexPartOrig].getForces()[0],forcesX[idxPart]);
                                                fy.add(particles[indexPartOrig].getForces()[1],forcesY[idxPart]);
                                                fz.add(particles[indexPartOrig].getForces()[2],forcesZ[idxPart]);
                                                energy += potentials[idxPart]*physicalValues[idxPart] ;
                                        }
                                });
                        }

                        // Print for information
                        std::cout << "FFmaBlas Energy "  << FMath::Abs(energy-energyD) /energyD << std::endl;
                        std::cout << "FFmaBlas Potential " << potentialDiff << std::endl;
                        std::cout << "FFmaBlas Fx " << fx << std::endl;
                        std::cout << "FFmaBlas Fy " << fy << std::endl;
                        std::cout << "FFmaBlas Fz " << fz << std::endl;
                } // end FFmaBlas kernel
        //
        ////////////////////////////////////////////////////////////////////
        //
        {	// begin FFmaBlockBlas kernel
                        std::cout << "\nFFmaBlockBlas FMM ... P: " <<DevP << std::endl;

                        // typedefs
            typedef FSphericalCell<FReal>                 CellClass;
            typedef FP2PParticleContainerIndexed<FReal>         ContainerClass;
            typedef FSimpleLeaf<FReal,  ContainerClass >                     LeafClass;
            typedef FOctree<FReal,  CellClass, ContainerClass , LeafClass >  OctreeClass;
            typedef FSphericalBlockBlasKernel<FReal,  CellClass, ContainerClass > KernelClass;
                        typedef FFmmAlgorithmThread<OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass > FmmClass;

                        // init cell class and oct-tree
                        CellClass::Init(DevP, true); // only for blas
                        OctreeClass tree(TreeHeight, SubTreeHeight, loader.getBoxWidth(), loader.getCenterOfBox());

                        { // -----------------------------------------------------
                                time.tic();

                for(FSize idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
                                        // put in tree
                                        tree.insert(particles[idxPart].getPosition(), idxPart, particles[idxPart].getPhysicalValue());
                                }

                                time.tac();
                                std::cout << "( FFmaBlockBlas@Inserting Particles = " << time.elapsed() << " s)." << std::endl;
                        } // -----------------------------------------------------

                        // -----------------------------------------------------
                        time.tic();
                        KernelClass kernels(DevP, TreeHeight, loader.getBoxWidth(), loader.getCenterOfBox());
                        FmmClass algorithm(&tree, &kernels);
                        algorithm.execute();
                        time.tac();
                        std::cout << "(FFmaBlockBlas @Algorithm = " << time.elapsed() << " s)." << std::endl;
                        // -----------------------------------------------------

                        FReal energy = 0.0;
            FMath::FAccurater<FReal> potentialDiff;
            FMath::FAccurater<FReal> fx, fy, fz;
                        { // Check that each particle has been summed with all other

                                tree.forEachLeaf([&](LeafClass* leaf){
                                        const FReal*const physicalValues = leaf->getTargets()->getPhysicalValues();
                                        const FReal*const potentials        = leaf->getTargets()->getPotentials();
                                        const FReal*const forcesX            = leaf->getTargets()->getForcesX();
                                        const FReal*const forcesY            = leaf->getTargets()->getForcesY();
                                        const FReal*const forcesZ            = leaf->getTargets()->getForcesZ();
                    const FSize nbParticlesInLeaf           = leaf->getTargets()->getNbParticles();
                    const FVector<FSize>& indexes       = leaf->getTargets()->getIndexes();

                    for(FSize idxPart = 0 ; idxPart < nbParticlesInLeaf ; ++idxPart){
                        const FSize indexPartOrig = indexes[idxPart];
                                                potentialDiff.add(particles[indexPartOrig].getPotential(),potentials[idxPart]);
                                                fx.add(particles[indexPartOrig].getForces()[0],forcesX[idxPart]);
                                                fy.add(particles[indexPartOrig].getForces()[1],forcesY[idxPart]);
                                                fz.add(particles[indexPartOrig].getForces()[2],forcesZ[idxPart]);
                                                energy += potentials[idxPart]*physicalValues[idxPart] ;
                                        }
                                });
                        }

                        // Print for information
                        std::cout << "FFmaBlockBlas Energy "  << FMath::Abs(energy-energyD) /energyD << std::endl;
                        std::cout << "FFmaBlockBlas Potential " << potentialDiff << std::endl;
                        std::cout << "FFmaBlockBlas Fx " << fx << std::endl;
                        std::cout << "FFmaBlockBlas Fy " << fy << std::endl;
                        std::cout << "FFmaBlockBlas Fz " << fz << std::endl;
                } // end FFmaBlas kernel

#endif

#ifdef  SCALFMM_USE_FFT
        //
        ////////////////////////////////////////////////////////////////////
        //
        {	// begin Lagrange/Uniform Grid kernel

                // TODO

                // accuracy
                const unsigned int ORDER = 7;
                std::cout << "\nLagrange FMM ... ORDER " << ORDER <<std::endl;

                // typedefs

        typedef FP2PParticleContainerIndexed<FReal> ContainerClass;
        typedef FSimpleLeaf<FReal,  ContainerClass >  LeafClass;
        typedef FInterpMatrixKernelR<FReal> MatrixKernelClass;
        typedef FUnifCell<FReal, ORDER> CellClass;
        typedef FOctree<FReal, CellClass,ContainerClass,LeafClass> OctreeClass;
        typedef FUnifKernel<FReal, CellClass,ContainerClass,MatrixKernelClass,ORDER> KernelClass;
                typedef FFmmAlgorithmThread<OctreeClass,CellClass,ContainerClass,KernelClass,LeafClass> FmmClass;


                // init oct-tree
                OctreeClass tree(TreeHeight, SubTreeHeight, loader.getBoxWidth(), loader.getCenterOfBox());

    // Create Matrix Kernel
    const MatrixKernelClass MatrixKernel;

                { // -----------------------------------------------------
                                time.tic();

            for(FSize idxPart = 0 ; idxPart <nbParticles ; ++idxPart){
                                // put in tree
                                tree.insert(particles[idxPart].getPosition(), idxPart, particles[idxPart].getPhysicalValue());
                        }

                        time.tac();
                        std::cout << "(Lagrange @Inserting Particles = " << time.elapsed() << " s)." << std::endl;
                } // -----------------------------------------------------

                { // -----------------------------------------------------
                        time.tic();
                        KernelClass kernels(TreeHeight, loader.getBoxWidth(), loader.getCenterOfBox(),&MatrixKernel);
                        FmmClass algorithm(&tree, &kernels);
                        algorithm.execute();
                        time.tac();
                        std::cout <<  "(Lagrange @Algorithm = " << time.elapsed() << " s)." << std::endl;
                } // -----------------------------------------------------

                FReal energy = 0.0;
        FMath::FAccurater<FReal> potentialDiff;
        FMath::FAccurater<FReal> fx, fy, fz;
                { // Check that each particle has been summed with all other

                        tree.forEachLeaf([&](LeafClass* leaf){
                                const FReal*const physicalValues = leaf->getTargets()->getPhysicalValues();
                                const FReal*const potentials        = leaf->getTargets()->getPotentials();
                                const FReal*const forcesX            = leaf->getTargets()->getForcesX();
                                const FReal*const forcesY            = leaf->getTargets()->getForcesY();
                                const FReal*const forcesZ            = leaf->getTargets()->getForcesZ();
                const FSize nbParticlesInLeaf           = leaf->getTargets()->getNbParticles();
                const FVector<FSize>& indexes       = leaf->getTargets()->getIndexes();

                for(FSize idxPart = 0 ; idxPart < nbParticlesInLeaf ; ++idxPart){
                    const FSize indexPartOrig = indexes[idxPart];
                                        potentialDiff.add(particles[indexPartOrig].getPotential(),potentials[idxPart]);
                                        fx.add(particles[indexPartOrig].getForces()[0],forcesX[idxPart]);
                                        fy.add(particles[indexPartOrig].getForces()[1],forcesY[idxPart]);
                                        fz.add(particles[indexPartOrig].getForces()[2],forcesZ[idxPart]);
                                        energy += potentials[idxPart]*physicalValues[idxPart] ;
                                }
                        });
                }

                // Print for information
                std::cout << "Lagrange Energy "  << FMath::Abs(energy-energyD) /energyD << std::endl;
                std::cout << "Lagrange Potential " << potentialDiff << std::endl;
                std::cout << "Lagrange Fx " << fx << std::endl;
                std::cout << "Lagrange Fy " << fy << std::endl;
                std::cout << "Lagrange Fz " << fz << std::endl;

        } // end Lagrange/Uniform Grid kernel
#endif
        //
        //         Spherical approximation
        //
        {
                //const static int P = 10;
        typedef FSphericalCell<FReal>               CellClass;
        typedef FP2PParticleContainerIndexed<FReal>          ContainerClass;
        typedef FSimpleLeaf<FReal,  ContainerClass >                     LeafClass;
        typedef FOctree<FReal,  CellClass, ContainerClass , LeafClass >  OctreeClass;
        typedef FSphericalKernel<FReal,  CellClass, ContainerClass >     KernelClass;
                typedef FFmmAlgorithmThread<OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass > FmmClass;

#ifndef  SCALFMM_USE_BLAS
                CellClass::Init(DevP, true);
#endif
                OctreeClass tree(TreeHeight, SubTreeHeight, loader.getBoxWidth(), loader.getCenterOfBox());
                std::cout << "\nFFmaSpherical FMM ... P: " << DevP<< std::endl;

                { // -----------------------------------------------------
                        time.tic();

            for(FSize idxPart = 0 ; idxPart < nbParticles; ++idxPart){
                                // put in tree
                                tree.insert(particles[idxPart].getPosition(), idxPart, particles[idxPart].getPhysicalValue());
                        }

                        time.tac();
                        std::cout << "(FFmaSpherical @Inserting Particles = "<< time.elapsed() << " s)." << std::endl;
                } // -----------------------------------------------------

                // -----------------------------------------------------
                time.tic();
            CellClass::Init(DevP);
                KernelClass *kernels = new KernelClass(DevP, TreeHeight, loader.getBoxWidth(), loader.getCenterOfBox());
                FmmClass algorithm(&tree, kernels);


                algorithm.execute();
                time.tac();
                std::cout << "(FFmaSpherical @Algorithm = " << time.elapsed() << " s)." << std::endl;
                // -----------------------------------------------------

                FReal energy = 0.0;
        FMath::FAccurater<FReal> potentialDiff;
        FMath::FAccurater<FReal> fx, fy, fz;
                { // Check that each particle has been summed with all other

                        tree.forEachLeaf([&](LeafClass* leaf){
                                const FReal*const physicalValues = leaf->getTargets()->getPhysicalValues();
                                const FReal*const potentials        = leaf->getTargets()->getPotentials();
                                const FReal*const forcesX            = leaf->getTargets()->getForcesX();
                                const FReal*const forcesY            = leaf->getTargets()->getForcesY();
                                const FReal*const forcesZ            = leaf->getTargets()->getForcesZ();
                const FSize nbParticlesInLeaf           = leaf->getTargets()->getNbParticles();
                const FVector<FSize>& indexes       = leaf->getTargets()->getIndexes();

                for(FSize idxPart = 0 ; idxPart < nbParticlesInLeaf ; ++idxPart){
                    const FSize indexPartOrig = indexes[idxPart];
                                        potentialDiff.add(particles[indexPartOrig].getPotential(),potentials[idxPart]);
                                        fx.add(particles[indexPartOrig].getForces()[0],forcesX[idxPart]);
                                        fy.add(particles[indexPartOrig].getForces()[1],forcesY[idxPart]);
                                        fz.add(particles[indexPartOrig].getForces()[2],forcesZ[idxPart]);
                                        energy += potentials[idxPart]*physicalValues[idxPart] ;
                                }
                        });
                }

                // Print for information
                std::cout << "FFmaSpherical Energy "  << FMath::Abs(energy-energyD) /energyD << std::endl;
                std::cout << "FFmaSpherical Potential " << potentialDiff << std::endl;
                std::cout << "FFmaSpherical Fx " << fx << std::endl;
                std::cout << "FFmaSpherical Fy " << fy << std::endl;
                std::cout << "FFmaSpherical Fz " << fz << std::endl;

        }
        //
        //         Spherical Rotation
        //
        {
                const static int P = 11;
        typedef FRotationCell<FReal, P>               CellClass;
        typedef FP2PParticleContainerIndexed<FReal>          ContainerClass;
        typedef FSimpleLeaf<FReal,  ContainerClass >                     LeafClass;
        typedef FOctree<FReal,  CellClass, ContainerClass , LeafClass >  OctreeClass;
        typedef FRotationKernel<FReal,  CellClass, ContainerClass , P>   KernelClass;
                typedef FFmmAlgorithmThread<OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass > FmmClass;

                OctreeClass tree(TreeHeight, SubTreeHeight, loader.getBoxWidth(), loader.getCenterOfBox());
                std::cout << "\nFFmaRotation FMM ... P: " << P<< std::endl;

                { // -----------------------------------------------------
                        time.tic();

            for(FSize idxPart = 0 ; idxPart < nbParticles; ++idxPart){
                                // put in tree
                                tree.insert(particles[idxPart].getPosition(), idxPart, particles[idxPart].getPhysicalValue());
                        }

                        time.tac();
                        std::cout << "(FFmaRotation @Inserting Particles = "<< time.elapsed() << " s)." << std::endl;
                } // -----------------------------------------------------

                // -----------------------------------------------------
                time.tic();
                KernelClass *kernels = new KernelClass(TreeHeight, loader.getBoxWidth(), loader.getCenterOfBox());
           FmmClass algorithm(&tree, kernels);


                algorithm.execute();
                time.tac();
                std::cout << "(FFmaRotation @Algorithm = " << time.elapsed() << " s)." << std::endl;
                // -----------------------------------------------------

                FReal energy = 0.0;
        FMath::FAccurater<FReal> potentialDiff;
        FMath::FAccurater<FReal> fx, fy, fz;
                { // Check that each particle has been summed with all other

                        tree.forEachLeaf([&](LeafClass* leaf){
                                const FReal*const physicalValues = leaf->getTargets()->getPhysicalValues();
                                const FReal*const potentials        = leaf->getTargets()->getPotentials();
                                const FReal*const forcesX            = leaf->getTargets()->getForcesX();
                                const FReal*const forcesY            = leaf->getTargets()->getForcesY();
                                const FReal*const forcesZ            = leaf->getTargets()->getForcesZ();
                const FSize nbParticlesInLeaf           = leaf->getTargets()->getNbParticles();
                const FVector<FSize>& indexes       = leaf->getTargets()->getIndexes();

                for(FSize idxPart = 0 ; idxPart < nbParticlesInLeaf ; ++idxPart){
                    const FSize indexPartOrig = indexes[idxPart];
                                        potentialDiff.add(particles[indexPartOrig].getPotential(),potentials[idxPart]);
                                        fx.add(particles[indexPartOrig].getForces()[0],forcesX[idxPart]);
                                        fy.add(particles[indexPartOrig].getForces()[1],forcesY[idxPart]);
                                        fz.add(particles[indexPartOrig].getForces()[2],forcesZ[idxPart]);
                                        energy += potentials[idxPart]*physicalValues[idxPart] ;
                                }
                        });
                }

                // Print for information
                std::cout << "FFmaRotation Energy "  << FMath::Abs(energy-energyD) /energyD << std::endl;
                std::cout << "FFmaRotation Potential " << potentialDiff << std::endl;
                std::cout << "FFmaRotation Fx " << fx << std::endl;
                std::cout << "FFmaRotation Fy " << fy << std::endl;
                std::cout << "FFmaRotation Fz " << fz << std::endl;

        }

        ////////////////////////////////////////////////////////////////////
        {	// begin Taylor kernel

                // accuracy
                const unsigned int ORDER = 10;

                // typedefs
        typedef FTaylorCell<FReal, ORDER,1>                                 CellClass;
                std::cout << "\nFFmaTaylor FMM ... ORDER: " << ORDER << std::endl;

        typedef FP2PParticleContainerIndexed<FReal>                          ContainerClass;
        typedef FSimpleLeaf<FReal,  ContainerClass >                         LeafClass;
        typedef FOctree<FReal,  CellClass, ContainerClass , LeafClass >      OctreeClass;
        typedef FTaylorKernel<FReal, CellClass,ContainerClass,ORDER,1>       KernelClass;
                typedef FFmmAlgorithmThread<OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass > FmmClass;

                // init cell class and oct-tree
                OctreeClass tree(TreeHeight, SubTreeHeight, loader.getBoxWidth(), loader.getCenterOfBox());

                { // -----------------------------------------------------
                        time.tic();

            for(FSize idxPart = 0 ; idxPart <nbParticles ; ++idxPart){
                                // put in tree
                                tree.insert(particles[idxPart].getPosition(), idxPart, particles[idxPart].getPhysicalValue());
                        }

                        time.tac();
                        std::cout <<"(FFmaTaylor @Inserting Particles = " << time.elapsed() << " s)." << std::endl;
                } // -----------------------------------------------------

                // -----------------------------------------------------
                time.tic();
                KernelClass kernels(TreeHeight, loader.getBoxWidth(), loader.getCenterOfBox());
                FmmClass algorithm(&tree, &kernels);
                algorithm.execute();
                time.tac();
                std::cout << "(FFmaTaylor @Algorithm = " << time.elapsed() << " s)." << std::endl;
                // -----------------------------------------------------

                FReal energy = 0.0;
        FMath::FAccurater<FReal> potentialDiff;
        FMath::FAccurater<FReal> fx, fy, fz;
                { // Check that each particle has been summed with all other

                        tree.forEachLeaf([&](LeafClass* leaf){
                                const FReal*const physicalValues = leaf->getTargets()->getPhysicalValues();
                                const FReal*const potentials        = leaf->getTargets()->getPotentials();
                                const FReal*const forcesX            = leaf->getTargets()->getForcesX();
                                const FReal*const forcesY            = leaf->getTargets()->getForcesY();
                                const FReal*const forcesZ            = leaf->getTargets()->getForcesZ();
                const FSize nbParticlesInLeaf           = leaf->getTargets()->getNbParticles();
                const FVector<FSize>& indexes       = leaf->getTargets()->getIndexes();

                for(FSize idxPart = 0 ; idxPart < nbParticlesInLeaf ; ++idxPart){
                    const FSize indexPartOrig = indexes[idxPart];
                                        potentialDiff.add(particles[indexPartOrig].getPotential(),potentials[idxPart]);
                                        fx.add(particles[indexPartOrig].getForces()[0],forcesX[idxPart]);
                                        fy.add(particles[indexPartOrig].getForces()[1],forcesY[idxPart]);
                                        fz.add(particles[indexPartOrig].getForces()[2],forcesZ[idxPart]);
                                        energy += potentials[idxPart]*physicalValues[idxPart] ;
                                }
                        });
                }

                // Print for information
                std::cout << "FFmaTaylor Energy "  << FMath::Abs(energy-energyD) /energyD << std::endl;
                std::cout << "FFmaTaylor Potential " << potentialDiff << std::endl;
                std::cout << "FFmaTaylor Fx " << fx << std::endl;
                std::cout << "FFmaTaylor Fy " << fy << std::endl;
                std::cout << "FFmaTaylor Fz " << fz << std::endl;
        } // end FFTaylor kernel
        delete[] particles;

        return 0;

}
