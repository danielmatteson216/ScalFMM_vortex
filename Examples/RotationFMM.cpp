// See LICENCE file at project root

// ==== CMAKE =====
//
// ================

#include <iostream>
#include <stdexcept>
#include <string>
#include <cstdlib>
#include <cstdio>

#include "ScalFmmConfig.h"
#include "Utils/FParameters.hpp"
#include "Files/FFmaGenericLoader.hpp"

#include "Kernels/Rotation/FRotationKernel.hpp"
#include "Kernels/Rotation/FRotationCell.hpp"

#include "Components/FSimpleLeaf.hpp"
#include "Kernels/P2P/FP2PParticleContainerIndexed.hpp"


#include "Containers/FOctree.hpp"

#ifdef _OPENMP
#include "Core/FFmmAlgorithmThread.hpp"
#else
#include "Core/FFmmAlgorithm.hpp"
#endif

#include "Utils/FParameterNames.hpp"

/// \file  RotationFMM.cpp
//!
//! \brief This program runs the FMM Algorithm with harmonic spherical approximation of 1/r kernel
//!  \authors B. Bramas, O. Coulaud
//!
//!  This code is a short example to use the  rotation harmonic spherical approximation  for the 1/r kernel
//!
//!
//!  <b> General arguments:</b>
//!     \param   -help(-h)      to see the parameters available in this driver
//!     \param   -depth          The depth of the octree
//!     \param   -subdepth     Specifies the size of the sub octree
//!     \param   -t                   The number of threads
//!
//!     \param   -f name          Name of the particles file with extension (.fma or .bfma). The data in  file have to be in our FMA format
//!
//


// Simply create particles and try the kernels
int main(int argc, char* argv[])
{
    FHelpDescribeAndExit(argc, argv,
                         "Driver for HArmonic Spherical + Rotation  --  kernel  (1/r kernel).",
                         FParameterDefinitions::InputFile, FParameterDefinitions::OctreeHeight,
                         FParameterDefinitions::OctreeSubHeight, FParameterDefinitions::OutputFile,
                         FParameterDefinitions::NbThreads);

    const std::string defaultFile(/*SCALFMMDataPath+*/"../Data/test20k.fma");
    const std::string filename(FParameters::getStr(argc,argv,FParameterDefinitions::InputFile.options, defaultFile.c_str()));
    const unsigned int TreeHeight       = FParameters::getValue(argc, argv, FParameterDefinitions::OctreeHeight.options, 5);
    const unsigned int SubTreeHeight  = FParameters::getValue(argc, argv, FParameterDefinitions::OctreeSubHeight.options, 2);
    const unsigned int NbThreads        = FParameters::getValue(argc, argv, FParameterDefinitions::NbThreads.options, 1);

#ifdef _OPENMP
	omp_set_num_threads(NbThreads) ;
	std::cout << "\n>> Using " << omp_get_max_threads() << " threads.\n" << std::endl;
#else
	std::cout << "\n>> Sequential version.\n" << std::endl;
#endif
	//
	std::cout <<	 "Parameters  "<< std::endl
			<<     "      Octree Depth      "<< TreeHeight <<std::endl
			<<	  "      SubOctree depth "<< SubTreeHeight <<std::endl
			<<     "      Input file  name: " <<filename <<std::endl
			<<     "      Thread number:  " << NbThreads <<std::endl
			<<std::endl;
	//
	// init timer
	FTic time;
	//
	// open particle file
	//
	////////////////////////////////////////////////////////////////////
	//
    typedef double FReal;
    FFmaGenericLoader<FReal> loader(filename);
	//
	////////////////////////////////////////////////////////////////////
	//
    // begin spherical kernel
	// accuracy
	const unsigned int P = 22;
	// typedefs
    typedef FP2PParticleContainerIndexed<FReal>                     ContainerClass;
    typedef FSimpleLeaf< FReal, ContainerClass >                       LeafClass;
    typedef FRotationCell<FReal, P>                                             CellClass;
    typedef FOctree<FReal, CellClass,ContainerClass,LeafClass>  OctreeClass;
	//
    typedef FRotationKernel< FReal, CellClass, ContainerClass , P>   KernelClass;
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
			// Read particles from file
			loader.fillParticle(&position,&physicalValue);

			// put particles in octree
			tree.insert(position, idxPart, physicalValue);
		}

		time.tac();
		std::cout << "Done  " << "(@Creating and Inserting Particles = "
				<< time.elapsed() << "  s) ." << std::endl;
	} // -----------------------------------------------------

	{ // -----------------------------------------------------
		std::cout << "\nRotation harmonic Spherical FMM (P="<< P << ") ... " << std::endl;

		time.tic();
		//
		// Here we use a pointer due to the limited size of the stack
		//
		KernelClass *kernels = new KernelClass(TreeHeight, loader.getBoxWidth(), loader.getCenterOfBox());
		//
		FmmClass algorithm(&tree, kernels) ;
		//
		algorithm.execute();   // Here the call of the FMM algorithm
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
					                                                                      << "   Forces: " << forcesX[idxPart] << " " << forcesY[idxPart] << " "<< forcesZ[idxPart] <<std::endl;
				}
				energy += potentials[idxPart]*physicalValues[idxPart] ;
			}
		});
		std::cout <<std::endl<<"Energy: "<< energy<<std::endl;
		std::cout <<std::endl<<" &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& "<<std::endl<<std::endl;

	}
	// -----------------------------------------------------
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
