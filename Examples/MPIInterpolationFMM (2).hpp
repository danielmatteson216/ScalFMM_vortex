// See LICENCE file at project root

// ==== CMAKE =====
// @FUSE_MPI
// @FUSE_BLAS
// ================

#include <iostream>
#include <stdexcept>
#include <cstdio>
#include <cstdlib>


#include "ScalFmmConfig.h"
#include "Containers/FOctree.hpp"
#include "Utils/FMpi.hpp"
#include "Core/FFmmAlgorithmThreadProc.hpp" // FMM code --> runs each stage of the algo
#include "Core/FFmmAlgorithmThreadProcPeriodic.hpp"

#include "Files/FFmaGenericLoader.hpp"      // particle loader
#include "Files/FMpiFmaGenericLoader.hpp"   // particle loader
#include "Files/FMpiTreeBuilder.hpp"        // tree builder

#include "Utils/FLeafBalance.hpp"

#include "Kernels/Interpolation/FInterpMatrixKernel.hpp"

#include "Components/FSimpleLeaf.hpp"

#include "Kernels/P2P/FP2PParticleContainerVortexIndexed.hpp"

#include "Utils/FParameters.hpp"
#include "Utils/FParameterNames.hpp"



//
// Order of the Interpolation approximation
static constexpr unsigned ORDER = 7 ;
using FReal                 = double;
// 
// MATRIX KERNEL CLASS
//using MatrixKernelClass     = FInterpMatrixKernelR<FReal> ;  // OLD KERNEL                                                                           //updated
using MatrixKernelClass     = FInterpMatrixKernelVORTEX<FReal>;  // VORTEX KERNEL                                                                      //updated

// CONTAINER CLASS
using ContainerClass = FP2PParticleContainerVortexIndexed<FReal>; // VORTEX PARTICAL CONTAINER                                 		                   //update with correct container
//using ContainerClass = FP2PParticleContainerIndexed<FReal>; // OLD PARTICAL CONTAINER                                        		                   //update

// LEAF CLASS
using LeafClass      = FSimpleLeaf<FReal, ContainerClass>;

// CELL CLASS
using CellClass      = FInterpolationCell<FReal, ORDER>;

// OCTREE CLASS
using OctreeClass    = FOctree<FReal,CellClass,ContainerClass,LeafClass>;

// MATRIX KERNEL CLASS
//using MatrixKernelClass = FInterpMatrixKernelR<FReal>;  // OLD KERNEL                                                       		                   //updated
using MatrixKernelClass = FInterpMatrixKernelVORTEX<FReal>; // VORTEX KERNEL                                                   		                   //updated
const MatrixKernelClass MatrixKernel;

// KERNEL CLASS
using KernelClass    = FInterpolationKernel<FReal,CellClass,ContainerClass,MatrixKernelClass,ORDER> ;

// FMM CLASS
using FmmClassProc     = FFmmAlgorithmThreadProc<OctreeClass,CellClass,ContainerClass,KernelClass,LeafClass>;
using FmmClassProcPER  = FFmmAlgorithmThreadProcPeriodic<FReal,OctreeClass,CellClass,ContainerClass,KernelClass,LeafClass>;





/// \file
//!
//! \brief This program runs the MPI FMM with Chebyshev interpolation of vortex kernel [ (1/r^2) * mollifier]
//!  \authors B. Bramas, O. Coulaud -> modified by Daniel Matteson for Vortex Thesis Project 06/03/2020
//!
//!  This code is a short example to use the FMM Algorithm Proc with Chebyshev or equispaced grid points Interpolation for the vortex kernel


// Create Particles then run the Kernel
int main(int argc, char* argv[])
{
	
		
// ---------------------- set MPI handling -----------------------------------	
  ///////// PARAMETERS HANDLING //////////////////////////////////////
  const FParameterNames  localIncreaseBox = { {"ratio","-L"}, "Increase the Box size by a factor L:= ratio"};
  FHelpDescribeAndExit(argc, argv,
                       "Driver for Chebyshev Interpolation kernel using MPI  (1/r kernel).\n "
                       "Usully run using : mpirun -np nb_proc_needed ./ChebyshevInterpolationAlgorithm [params].",
                       FParameterDefinitions::OctreeHeight,
                       FParameterDefinitions::OctreeSubHeight,
                       FParameterDefinitions::InputFile,
                       FParameterDefinitions::OutputFile,
                       FParameterDefinitions::NbThreads,
                       FParameterDefinitions::PeriodicityNbLevels,
                       localIncreaseBox
                       ) ;

  // Initialize values for MPI
  FMpi app(argc,argv);
  const bool masterIO = ( app.global().processId() == 0 );
  //
//const std::string defaultFile(SCALFMMDataPath+"unitCubeXYZF121.bfma");  // create string to load file with particle data
//const std::string defaultFile(SCALFMMDataPath+"unitCubeXYZF289.bfma");  // create string to load file with particle data
//const std::string defaultFile(SCALFMMDataPath+"unitCubeXYZF529.bfma");  // create string to load file with particle data   
//const std::string defaultFile(SCALFMMDataPath+"unitCubeXYZF1024.bfma");  // create string to load file with particle data   
//const std::string defaultFile(SCALFMMDataPath+"unitCubeXYZF2025.bfma");  // create string to load file with particle data   
//const std::string defaultFile(SCALFMMDataPath+"unitCubeXYZF3969.bfma");  // create string to load file with particle data   
//const std::string defaultFile(SCALFMMDataPath+"unitCubeXYZF7921.bfma");  // create string to load file with particle data   
//const std::string defaultFile(SCALFMMDataPath+"unitCubeXYZF15625.bfma");  // create string to load file with particle data   
const std::string defaultFile(SCALFMMDataPath+"unitCubeXYZF31329.bfma");  // create string to load file with particle data    
//const std::string defaultFile(SCALFMMDataPath+"unitCubeXYZF62001.bfma");  // create string to load file with particle data 
//const std::string defaultFile(SCALFMMDataPath+"unitCubeXYZF123904.bfma");  // create string to load file with particle data 
//const std::string defaultFile(SCALFMMDataPath+"unitCubeXYZF248004.bfma");  // create string to load file with particle data 
//const std::string defaultFile(SCALFMMDataPath+"unitCubeXYZF301401.bfma");  // create string to load file with particle data

  const std::string  filename      = FParameters::getStr(argc,argv,    FParameterDefinitions::InputFile.options, defaultFile.c_str());
  const unsigned int TreeHeight    = FParameters::getValue(argc, argv, FParameterDefinitions::OctreeHeight.options, 10);    //set tree depth
  const unsigned int SubTreeHeight = FParameters::getValue(argc, argv, FParameterDefinitions::OctreeSubHeight.options, 2);
  const unsigned int NbThreads     = FParameters::getValue(argc, argv, FParameterDefinitions::NbThreads.options, 1);
  bool periodicCondition = false ;
  
  if(FParameters::existParameter(argc, argv, FParameterDefinitions::PeriodicityNbLevels.options)){
      periodicCondition = true;
    }
  const unsigned int aboveTree = FParameters::getValue(argc, argv, FParameterDefinitions::PeriodicityNbLevels.options, 5);

  omp_set_num_threads(NbThreads);
  if(masterIO){
    std::cout << "\n>> Using " << omp_get_max_threads() << " threads.\n" << std::endl;
    
    //
    std::cout << "Parameters"<< std::endl
	      << "      Octree Depth      " << TreeHeight    << std::endl
	      << "      SubOctree depth   " << SubTreeHeight << std::endl;
    if(periodicCondition){
      std::cout << "      AboveTree    "<< aboveTree <<std::endl;
      
    }
    std::cout    << "      Input file  name: " << filename      << std::endl
		 << "      Thread count :    " << NbThreads     << std::endl
		 << std::endl;
  }




  // Initialize timer
  FTic time;

	  ///////// VAR INIT /////////////////////////////////////////////////
  // Creation of the particle loader
  FMpiFmaGenericLoader<FReal> loader(filename,app.global());
  if(!loader.isOpen()) {
      throw std::runtime_error("Particle file couldn't be opened!") ;
    }
  auto boxWidth = loader.getBoxWidth() ;
  //
  if(FParameters::existParameter(argc, argv, localIncreaseBox.options)){
      FReal ratio=  FParameters::getValue(argc, argv, localIncreaseBox.options, 1.0);
      boxWidth *= ratio;
    }

  // Initialize empty oct-tree
  OctreeClass tree(TreeHeight, SubTreeHeight, boxWidth, loader.getCenterOfBox());

  FSize localParticlesNumber = 0 ;

  // -----------------------------------------------------
  if(masterIO){
      std::cout << "Loading & Inserting " << loader.getNumberOfParticles()
                << " particles ..." << std::endl
                <<" Box: "<< std::endl
               << "    width  " << boxWidth << std::endl
               << "    Centre " << loader.getCenterOfBox()<< std::endl;
      std::cout << "\tHeight : " << TreeHeight << " \t sub-height : " << SubTreeHeight << std::endl;
    }
  time.tic();

  /* Mock particle structure to balance the tree over the processes. */
  struct TestParticle{
    FSize index;             // Index of the particle in the original file.
    FPoint<FReal> position;  // Spatial position of the particle.
    FReal physicalValue;     // Physical value of the particle.
    /* Returns the particle position. */
    const FPoint<FReal>& getPosition(){
      return position;
    }
  };

  // Temporary array of particles read by this process.
  TestParticle* particles = new TestParticle[loader.getMyNumberOfParticles()];
  memset(particles, 0, (sizeof(TestParticle) * loader.getMyNumberOfParticles()));

  // Index (in file) of the first particle that will be read by this process.
  FSize idxStart = loader.getStart();
  std::cout << "Proc:" << app.global().processId() << " start-index: " << idxStart << std::endl;

  // Read particles from parts.
  for(FSize idxPart = 0 ; idxPart < loader.getMyNumberOfParticles() ; ++idxPart){
      // Store the index (in the original file) the particle.
      particles[idxPart].index = idxPart + idxStart;
      // Read particle from file
      loader.fillParticle(&particles[idxPart].position,
                          &particles[idxPart].physicalValue);
    }

  // Final vector of particles
  FVector<TestParticle> finalParticles;
  FLeafBalance balancer;
  // Redistribute particules between processes
  FMpiTreeBuilder< FReal, TestParticle >::
      DistributeArrayToContainer(app.global(),
                                 particles,
                                 loader.getMyNumberOfParticles(),
                                 tree.getBoxCenter(),
                                 tree.getBoxWidth(),
                                 tree.getHeight(),
                                 &finalParticles,
                                 &balancer);

  // Free temporary array memory.
  delete[] particles;

  // Insert final particles into tree.

  for(FSize idx = 0 ; idx < finalParticles.getSize(); ++idx){
      tree.insert(finalParticles[idx].position,
                  finalParticles[idx].index,
                  finalParticles[idx].physicalValue);
    }

  time.tac();
// ---------------------- particles inserted into tree -----------------------------------



  localParticlesNumber = finalParticles.getSize() ;

  double timeUsed = time.elapsed();
  double minTime,maxTime;
  std::cout << "Proc:" << app.global().processId()
            << " "     << finalParticles.getSize()
            << " particles have been inserted in the tree. (@Reading and Inserting Particles = "
            << time.elapsed() << " s)."
            << std::endl;

  MPI_Reduce(&timeUsed,&minTime,1,MPI_DOUBLE,MPI_MIN,0,app.global().getComm());
  MPI_Reduce(&timeUsed,&maxTime,1,MPI_DOUBLE,MPI_MAX,0,app.global().getComm());
  if(masterIO){
      std::cout << "readinsert-time-min:" << minTime
                << " readinsert-time-max:" << maxTime
                << std::endl;
    }
  
  // -----------------------------------------------------
  FAbstractAlgorithm * algorithm  = nullptr;
  FAlgorithmTimers   * timer      = nullptr;
   // -----------------------------------------------------
    if(masterIO) {
        std::cout << "\n"<<interpolationType<<" FMM Proc (ORDER="<< ORDER << ") ... " << std::endl;
    }

    // Kernels to use (pointer because of the limited size of the stack)

    // non periodic FMM algorithm ---> to build constructor?
    std::unique_ptr<KernelClass> kernelsNoPer(new KernelClass(TreeHeight, boxWidth,
                                                              loader.getCenterOfBox(),
							      &MatrixKernel));

    FmmClassProc    algoNoPer(app.global(),&tree, kernelsNoPer.get());  

    // periodic FMM algorithm
    FmmClassProcPER algoPer(app.global(),&tree, aboveTree);

	// to build constructor?
	KernelClass kernelsPer(algoPer.extendedTreeHeight(), algoPer.extendedBoxWidth(),
                           algoPer.extendedBoxCenter(),&MatrixKernel);

	algoPer.setKernel(&kernelsPer);  //copy constructor here


	
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    if(! periodicCondition) {// Non periodic case
        algorithm  = &algoNoPer ;
        timer      = &algoNoPer ;
      }
    else {  // Periodic case
        algorithm  = &algoPer ;
        timer      = &algoPer ;
      }
    ///////////////////////////////////////////////////////////////////////////////////////////////////
	
  //
    // FMM exectution  FFmmFarField FFmmNearField  FFmmP2M|FFmmM2M|FFmmM2L|FFmmL2L
    time.tic();
    algorithm->execute();
    time.tac();
    //



   // if(masterIO)
   //  {
   //    std::cout << app.global().processId()  <<"  Morton distribution "<< std::endl ;
   //    for(auto v : mortonLeafDistribution)
   //      std::cout << v << " ";

   //    std::cout <<  std::endl;
   //  }
    //    app.global().barrier();
    timeUsed = time.elapsed();
    MPI_Reduce(&timeUsed,&minTime,1,MPI_DOUBLE,MPI_MIN,0,app.global().getComm());
    MPI_Reduce(&timeUsed,&maxTime,1,MPI_DOUBLE,MPI_MAX,0,app.global().getComm());
	
    if(masterIO){
        std::cout << "Done  " << "(@Algorithm = " << time.elapsed() << "   s)." << std::endl;
        std::cout << "exec-time-min:   " << minTime
                  << " exec-time-max:   " << maxTime
                  << std::endl;
    std::cout << "Timers Far Field \n"
              << "P2M " << timer->getTime(FAlgorithmTimers::P2MTimer) << " seconds\n"
              << "M2M " << timer->getTime(FAlgorithmTimers::M2MTimer) << " seconds\n"
              << "M2L " << timer->getTime(FAlgorithmTimers::M2LTimer) << " seconds\n"
              << "L2L " << timer->getTime(FAlgorithmTimers::L2LTimer) << " seconds\n"
              << "P2P and L2P " << timer->getTime(FAlgorithmTimers::NearTimer) << " seconds\n"
	      << std::endl;

      }
  

  // -----------------------------------------------------
  //
  // Some output
  //
  //
  { // -----------------------------------------------------
    FSize N1=0, N2= loader.getNumberOfParticles()/2, N3= (loader.getNumberOfParticles()-1); ;
    FReal energy =0.0 ;
    //
    //   Loop over all leaves
    //
    std::cout <<std::endl<<" &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& "<<std::endl;
    std::cout << std::scientific;
    std::cout.precision(15) ;

    FReal locTotalPhysicalValue=0.0 ;

    tree.forEachLeaf([&](LeafClass* leaf){
      const FReal*const posX = leaf->getTargets()->getPositions()[0];
      const FReal*const posY = leaf->getTargets()->getPositions()[1];
      const FReal*const posZ = leaf->getTargets()->getPositions()[2];

      //const FReal*const potentials = leaf->getTargets()->getPotentials();
      const FReal*const potentials = leaf->getTargets()->getPotentials_real();	
      const FReal*const potentials_i = leaf->getTargets()->getPotentials_imag();	
	  //needs updated
      const FReal*const forcesX = leaf->getTargets()->getForcesX_real();																									
      const FReal*const forcesY = leaf->getTargets()->getForcesY_real();																									
      const FReal*const forcesZ = leaf->getTargets()->getForcesZ_real();																									
      const FReal*const forcesX_i = leaf->getTargets()->getForcesX_imag();																									
      const FReal*const forcesY_i = leaf->getTargets()->getForcesY_imag();																									
      const FReal*const forcesZ_i = leaf->getTargets()->getForcesZ_imag();																										  
      const FSize nbParticlesInLeaf = leaf->getTargets()->getNbParticles();
      const FReal*const physicalValues = leaf->getTargets()->getPhysicalValues();

      const FVector<FSize>& indexes = leaf->getTargets()->getIndexes();



      for(FSize idxPart = 0 ; idxPart < nbParticlesInLeaf ; ++idxPart){
	const FSize indexPartOrig = indexes[idxPart];
//	if ((indexPartOrig == N1) || (indexPartOrig == N2) || (indexPartOrig == N3)  )
//	  {
	//*
	     std::cout << "Proc "<< app.global().processId() << " Index "<< indexPartOrig <<"  potential  " << potentials[idxPart]										
		      << " Pos "<<posX[idxPart]<<" "<<posY[idxPart]<<" "<<posZ[idxPart]
		      << "   ForcesReal: " << forcesX[idxPart] << " " << forcesY[idxPart] << " "<< forcesZ[idxPart]
			  << "   ForcesImag: " << forcesX_i[idxPart] << " " << forcesY_i[idxPart] << " "<< forcesZ_i[idxPart] << std::endl;												
	//*/
//	  }
	energy += potentials[idxPart]*physicalValues[idxPart] ;																											

	//energy += potentials_i[idxPart]*physicalValues[idxPart] ;																											
	
	locTotalPhysicalValue += physicalValues[idxPart]  ;
      }
	  
	  
      });
    FReal gloEnergy          = app.global().reduceSum(energy);
    FReal TotalPhysicalValue = app.global().reduceSum(locTotalPhysicalValue);
    if(masterIO){
      std::cout <<std::endl<<"aboveRoot: " << aboveTree << " Energy: "<< gloEnergy <<"  TotalPhysicalValue: " << TotalPhysicalValue<< std::endl; 
      std::cout <<std::endl <<" &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& "<<std::endl<<std::endl;
    }
  }


  // -----------------------------------------------------
  if(FParameters::existParameter(argc, argv, FParameterDefinitions::OutputFile.options)){
    std::vector<MortonIndex> mortonLeafDistribution(2*app.global().processCount());
    algorithm->getMortonLeafDistribution(mortonLeafDistribution);
    std::string name(FParameters::getStr(argc,argv,FParameterDefinitions::OutputFile.options, "output.fma"));
    FMpiFmaGenericWriter<FReal> paraWriter(name,app);
    paraWriter.writeDistributionOfParticlesFromOctree(tree,loader.getNumberOfParticles(),localParticlesNumber,
						      mortonLeafDistribution);

  }
  
  return 0;
}
