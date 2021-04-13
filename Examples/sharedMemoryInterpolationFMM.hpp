// -*-c++-*-
// See LICENCE file at project root

// ==== CMAKE =====
// @FUSE_FFT
// @FUSE_BLAS
//  ==== Git =====

// ================

/** \brief Uniform FMM example
 *
 * \file
 * \authors B. Bramas, O. Coulaud
 *
 * This program runs the FMM Algorithm with the interpolation kernel based on
 * either uniform (grid points) or Chebychev interpolation (1/r kernel). 
 * It then compares the results with a direct computation.
 */

#include <iostream>
#include <iomanip>
#include <memory>

#include <cstdio>
#include <cstdlib>
#include <string>


#include "ScalFmmConfig.h"
#include "Utils/FGlobal.hpp"

#include "Utils/FParameters.hpp"
#include "Utils/FParameterNames.hpp"

#include "Files/FFmaGenericLoader.hpp"
// Leaves
#include "Components/FSimpleLeaf.hpp"
#include "Kernels/P2P/FP2PParticleContainerIndexed.hpp"
// Octree
#include "Containers/FOctree.hpp"
//
#include "Kernels/Interpolation/FInterpMatrixKernel.hpp"
//
//
// Order of the Interpolation approximation
static constexpr unsigned ORDER = 6 ;
using FReal                 = double;
//   1/r kernel
//
using MatrixKernelClass     = FInterpMatrixKernelR<FReal> ;

//
/// definition of the common tree structure
using CellClass      = FInterpolationCell<FReal, ORDER>;
//using CellUpClass    = typename CellClass::multipole_t;
//using CellDownClass  = typename CellClass::local_expansion_t;
//using CellSymbClass  = FSymbolicData;  
using ContainerClass = FP2PParticleContainerIndexed<FReal>;
using LeafClass      = FSimpleLeaf<FReal,  ContainerClass >   ;
using OctreeClass    = FOctree<FReal, CellClass,ContainerClass,LeafClass>;
using KernelClass    = FInterpolationKernel<FReal,CellClass,ContainerClass,MatrixKernelClass,ORDER> ;


#ifdef _OPENMP
//#include "Core/FFmmAlgorithmThread.hpp"
#include "Core/FFmmAlgorithmPeriodic.hpp"
#include "Core/FFmmAlgorithmSectionTask.hpp"
#else
#include "Core/FFmmAlgorithm.hpp"
#endif


#ifdef _OPENMP
//using FmmClass = FFmmAlgorithmThread<OctreeClass,CellClass,ContainerClass,KernelClass,LeafClass> FmmClass;
using FmmClass    = FFmmAlgorithmSectionTask<OctreeClass,CellClass,ContainerClass,KernelClass,LeafClass> ;
using FmmClassPer = FFmmAlgorithmPeriodic<FReal,OctreeClass,CellClass,ContainerClass,KernelClass,LeafClass>;
#else
using FmmClass  = FFmmAlgorithm<OctreeClass,CellClass,ContainerClass,KernelClass,LeafClass>;
#endif



// Simply create particles and try the kernels
int main(int argc, char* argv[])
{
  FHelpDescribeAndExit(
		       argc, argv,
		       "Driver for Lagrange interpolation kernel  (1/r kernel).",
		       FParameterDefinitions::OctreeHeight,
		       FParameterDefinitions::OctreeSubHeight,
		       FParameterDefinitions::InputFile,
		       FParameterDefinitions::OutputFile,
		       FParameterDefinitions::NbThreads,
		       FParameterDefinitions::PeriodicityNbLevels
		       );


  const std::string defaultFile(SCALFMMDataPath+"unitCubeXYZQ100.bfma" );
  const std::string filename       = FParameters::getStr(argc,argv,FParameterDefinitions::InputFile.options, defaultFile.c_str());
  const int TreeHeight    = FParameters::getValue(argc, argv, FParameterDefinitions::OctreeHeight.options, 5);
  const int SubTreeHeight = FParameters::getValue(argc, argv, FParameterDefinitions::OctreeSubHeight.options, 2);
  bool periodicCondition = false ;
  if(FParameters::existParameter(argc, argv, FParameterDefinitions::PeriodicityNbLevels.options)){
      periodicCondition = true;
    }
  const unsigned int aboveTree = FParameters::getValue(argc, argv, FParameterDefinitions::PeriodicityNbLevels.options, -1);

#ifdef _OPENMP
  const unsigned int NbThreads = periodicCondition ? 1 : FParameters::getValue(argc, argv, FParameterDefinitions::NbThreads.options, omp_get_max_threads());

   omp_set_num_threads(NbThreads);

  std::cout << "\n>> Using " << NbThreads << " threads.\n" << std::endl;
#else
  const int NbThreads =  1;
  std::cout << "\n>> Sequential version.\n" << std::endl;
#endif
  //
  {
    std::string indent("    ");
    auto w = std::setw(18);
    std::cout << "Parameters" << std::endl << std::left
              << indent << w << "Octree Depth: "     << TreeHeight    << std::endl
              << indent << w << "SubOctree depth: "  << SubTreeHeight << std::endl;
    if(periodicCondition){
        std::cout << indent << w << "AboveTree    "<< aboveTree <<std::endl;

      }
    std::cout << indent << w << "Input file  name: " << filename      << std::endl
              << indent << w << "Thread number: "    << NbThreads     << std::endl
              << std::endl;
  }
  //
  // init timer
  FTic time;

  // open particle file
  ////////////////////////////////////////////////////////////////////
  //
  FFmaGenericLoader<FReal> loader(filename);
  //
  // init oct-tree
  OctreeClass tree(TreeHeight, SubTreeHeight, loader.getBoxWidth(), loader.getCenterOfBox());
  //
  // Read particles and insert them in octree
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
  }
  ////////////////////////////////////////////////////////////////////
  //
  //    Execute FMM Algorithm
  //
  ////////////////////////////////////////////////////////////////////

 // { // -----------------------------------------------------
    std::cout << "\n" << interpolationType << "  FMM (ORDER= "<< ORDER << ") ... " << std::endl;

    const MatrixKernelClass    MatrixKernel;
    time.tic();
    //
    std::unique_ptr<KernelClass> kernelsNoPer(new KernelClass(TreeHeight, loader.getBoxWidth(), loader.getCenterOfBox(),&MatrixKernel));
    //
    FmmClass algoNoPer(&tree, kernelsNoPer.get());
    // periodic FMM algorithm
    FmmClassPer algoPer(&tree, aboveTree);
    KernelClass kernelsPer(algoPer.extendedTreeHeight(), algoPer.extendedBoxWidth(),
                           algoPer.extendedBoxCenter(),&MatrixKernel);
    algoPer.setKernel(&kernelsPer);
    //
    FAbstractAlgorithm * algorithm  = nullptr;
    FAlgorithmTimers   * timer      = nullptr;
    if(! periodicCondition) { // Non periodic case
        algorithm  = &algoNoPer ;
        timer      = &algoNoPer ;
      }
    else {                    // Periodic case
        algorithm  = &algoPer ;
        timer      = &algoPer ;

      }

    //
    algorithm->execute();   // Here the call of the FMM algorithm
    //
    time.tac();
    std::cout << "Timers Far Field \n"
              << "P2M " << timer->getTime(FAlgorithmTimers::P2MTimer) << " seconds\n"
              << "M2M " << timer->getTime(FAlgorithmTimers::M2MTimer) << " seconds\n"
              << "M2L " << timer->getTime(FAlgorithmTimers::M2LTimer) << " seconds\n"
              << "L2L " << timer->getTime(FAlgorithmTimers::L2LTimer) << " seconds\n"
              << "P2P and L2P " << timer->getTime(FAlgorithmTimers::NearTimer) << " seconds\n"
	      << std::endl;


    std::cout << "Done  " << "(@Algorithm = " << time.elapsed() << " s) ." << std::endl;
  //}
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
    std::cout.precision(15) ;
//     std::size_t count =0 , numLeaf =0 ;
    FReal TotalPhysicalValue=0.0 ;
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
	  TotalPhysicalValue += physicalValues[idxPart]  ;
	//  ++count ;

	}
//	std::cout << "numLeaf "<< numLeaf << "  " << energy <<"  " << count<<std::endl<<std::endl ;
//	++numLeaf ;
      });
    std::cout <<std::endl<<"aboveRoot: " << aboveTree << "  Energy: "<< energy<<"  TotalPhysicalValue: " << TotalPhysicalValue<< std::endl;
    std::cout <<std::endl<<" &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& "<<std::endl<<std::endl;

  }
  //
  // -----------------------------------------------------
  //
  if(FParameters::existParameter(argc, argv, FParameterDefinitions::OutputFile.options)){
    std::string name(FParameters::getStr(argc,argv,FParameterDefinitions::OutputFile.options,   "output.fma"));
    
    FFmaGenericWriter<FReal> writer(name) ;
    writer.writeDataFromOctree(&tree,loader.getNumberOfParticles());

  }


  return 0;
}
