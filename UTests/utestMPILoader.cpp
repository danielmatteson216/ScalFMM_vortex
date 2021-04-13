// See LICENCE file at project root

// ==== CMAKE =====
// @FUSE_MPI
// ================

#include <iostream>
#include <stdexcept>
#include <cstdio>
#include <cstdlib>


#include "ScalFmmConfig.h"
#include "Utils/FMpi.hpp"

#include "Files/FFmaGenericLoader.hpp"
#include "Files/FMpiFmaGenericLoader.hpp"
#include "FUTester.hpp"


#include "Utils/FParameters.hpp"

#include "Utils/FParameterNames.hpp"

/// \file  utestMPILoader.cpp
//!
//! \brief This program check the MPI Loader
//!  \authors B. Bramas, O. Coulaud
//!
//!
//
/** the test class
 *
 */
class TestMpiLoader : public FUTesterMpi<TestMpiLoader> {

  ///////////////////////////////////////////////////////////
  // The tests!
  ///////////////////////////////////////////////////////////

  template <class FReal>
  void RunTest()	{
    FReal eps = 1.0e-10;
    //
    if(sizeof(FReal) == sizeof(float) ) {
        std::cerr << "No input data available for Float "<< std::endl;
        std::exit(EXIT_FAILURE);
      }
    const std::string parFile("test20k.fma");
    //
    std::string filename(SCALFMMDataPath+parFile);
    //
    FMpiFmaGenericLoader<FReal> mpiLoader(filename,app.global());
    FFmaGenericLoader<FReal>    loader(filename);            // sequential loader

    if(!mpiLoader.isOpen()) throw std::runtime_error("Particle file couldn't be opened!") ;
    ////////////////////////////////////////////////////////////////////

    { // -----------------------------------------------------
      if(app.global().processId() == 0){
          std::cout << "Creating & Inserting " << mpiLoader.getNumberOfParticles()
                    << " particles ..." << std::endl;
        }
      //

      struct TestParticle{
        FSize indexInFile;
        FPoint<FReal> position;
        FReal physicalValue;
        const FPoint<FReal>& getPosition(){
          return position;
        }
        const FReal diffPosition(TestParticle& otherParticles){
          return (position - otherParticles.getPosition()).norm();
        }
      };
      FSize localNumberOfParticles  = mpiLoader.getMyNumberOfParticles() ;
      FSize globalNumberOfParticles = loader.getNumberOfParticles() ;
      TestParticle* particles = new TestParticle[localNumberOfParticles];
      memset(particles, 0, (sizeof(TestParticle) * localNumberOfParticles));
      TestParticle* allParticles = new TestParticle[globalNumberOfParticles];
      memset(allParticles, 0, (sizeof(TestParticle) * globalNumberOfParticles));
      //  READ IN PARALLEL the particles
      //idx (in file) of the first part that will be used by this proc.
      std::cout << "Load particles in parallel" <<std::endl;
      FSize idxStart = mpiLoader.getStart();
      for(FSize idxPart = 0 ; idxPart < localNumberOfParticles ; ++idxPart){
          //Storage of the index (in the original file) of each part.
          particles[idxPart].indexInFile = idxPart + idxStart;
          // Read particles from file
          mpiLoader.fillParticle(&particles[idxPart].position,&particles[idxPart].physicalValue);
        }
      //
      // Read All Particles
      std::cout << "Load all particles in sequential" <<std::endl;

      for(FSize idxPart = 0 ; idxPart < globalNumberOfParticles ; ++idxPart){
          //Storage of the index (in the original file) of each part.
          allParticles[idxPart].indexInFile = idxStart;
          // Read particles from file
          loader.fillParticle(&allParticles[idxPart].position,&allParticles[idxPart].physicalValue);
        }
      //
      // Check if the parallel read is OK
      //
      std::cout << " Compare the txo set of particles" <<std::endl;

      bool check ;
      for(FSize idxPart = 0 ; idxPart < localNumberOfParticles ; ++idxPart){
          //Storage of the index (in the original file) of each part.
          check = (particles[idxPart].diffPosition( allParticles[idxStart+idxPart])) < eps ;
          if(!check) {
              std::cout << "Proc "<< app.global().processId()
                           << " seqPos [ "<<allParticles[idxStart+idxPart].getPosition()<<" ]  "
                           << " parPos: "<<particles[idxPart].getPosition()<<" ]  " <<std::endl;
              break ;
            }
          // Read particles from file
        }
      MPI_Reduce(MPI_IN_PLACE,&check,1,MPI_LOGICAL,MPI_LAND,0,app.global().getComm());
//
      Print("Test1 - Read in parallel  ");
      uassert(check);
      //
    } // -----------------------------------------------------
  }
  /** If memstas is running print the memory used */
  void PostTest() {
    if( FMemStats::controler.isUsed() ){
        std::cout << app.global().processId() << "-> Memory used at the end " << FMemStats::controler.getCurrentAllocated()
                  << " Bytes (" << FMemStats::controler.getCurrentAllocatedMB() << "MB)\n";
        std::cout << app.global().processId() << "-> Max memory used " << FMemStats::controler.getMaxAllocated()
                  << " Bytes (" << FMemStats::controler.getMaxAllocatedMB() << "MB)\n";
        std::cout << app.global().processId() << "-> Total memory used " << FMemStats::controler.getTotalAllocated()
                  << " Bytes (" << FMemStats::controler.getTotalAllocatedMB() << "MB)\n";
      }
  }
  ///////////////////////////////////////////////////////////
  // Set the tests!
  ///////////////////////////////////////////////////////////
  /** TestUnifKernel */
  void TestLoader(){
      typedef double FReal;
      // run test
      RunTest<FReal>();
  }

  /** set test */
  void SetTests(){
    AddTest(&TestMpiLoader::TestLoader,"Compare MpiLoader with sequential Loader");
  }
public:
  TestMpiLoader(int argc,char ** argv) : FUTesterMpi(argc,argv){
  }

};

TestClassMpi(TestMpiLoader);
