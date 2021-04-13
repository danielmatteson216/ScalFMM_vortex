/*
 * genarateDistributions.cpp
 *
 *  Created on: 23 mars 2014
 *      Author: Olivier Coulaud
 */


#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
//
#include "Utils/FGlobal.hpp"
#include "Utils/FPoint.hpp"
#include "Files/FFmaGenericLoader.hpp"
#include "Utils/FParameters.hpp"
#include "Utils/FCompareResults.hpp"
#include "Utils/FParameterNames.hpp"
#include "Adaptive/FBox.hpp"
//
/// \file  compare2files.cpp
//!
//! \brief compare2files: Gives the error between
//!
//!  Driver to transform a FMA format and/or to build a visualization file<br>
//! For a description of the FMA format see FFmaGenericLoader<br>
//!  <b> General arguments:</b>
//!     \param   -help (-h)      to see the parameters available in this driver
//!     \param   -fmmfile1 name1:  first file name to compare (with extension .fma (ascii) or bfma (binary)
//!     \param   -fmmfile2 name2: second file name to compare (with extension .fma (ascii) or bfma (binary)
//!     \param   -ewaldfile2 name2 if name2 contains the result done by the ewald method for 1/r kernel.
//!
//! Error code (return of the main)
//!      -1   Number of points is different in the two files.
//!
//! \b examples
//!
//!
//!   compare2files -file1 unitCubeXYZQ100.fma  -file2 unitCubeXYZQ100

using FReal = double;

int main(int argc, char ** argv){
      const FParameterNames LocalSortParameter { {"-sort"}, "sort files according to the maximal morton index available"};
  FHelpDescribeAndExit(argc, argv,
                       "Driver to change the format of the input file.",
                       FParameterDefinitions::InputFileOne, FParameterDefinitions::InputFileTwow,
                       LocalSortParameter);

  const std::string filename1(FParameters::getStr(argc,argv,FParameterDefinitions::InputFileOne.options,   "data.fma"));
  const std::string filename2(FParameters::getStr(argc,argv,FParameterDefinitions::InputFileTwow.options,  "data.fma"));

  FFmaGenericLoader<FReal> loader1(filename1);
  FFmaGenericLoader<FReal> loader2(filename2);
  //
  // Allocation
  //
  FSize nbParticles          = loader1.getNumberOfParticles();
  const unsigned int nbData   = loader1.getNbRecordPerline() ;
  if(nbParticles !=  loader2.getNumberOfParticles()){
      std::cerr << "Number of points is different in the two files."<<std::endl ;
      return -1 ;
    }
  if( (nbData !=  8) && (loader2.getNbRecordPerline()) ){
      std::cerr << "Wrong files only " << std::min(loader2.getNbRecordPerline(),nbData)<<  " to read."<<std::endl ;
      return -2 ;
    }
  FmaRWParticle<FReal,8,8>*  particles1 = new FmaRWParticle<FReal,8,8>[nbParticles];
  FmaRWParticle<FReal,8,8>*  particles2 = new FmaRWParticle<FReal,8,8>[nbParticles];


  //
  loader1.fillParticle(particles1,nbParticles);
  loader2.fillParticle(particles2,nbParticles);
  if(FParameters::existParameter(argc, argv, LocalSortParameter.options)) {
      // define a box, used in the sort
      const FBox<FPoint<FReal>> box{loader1.getBoxWidth(),loader1.getCenterOfBox()};
      std::cout << "Sort needed !! " << std::endl;
      sortArrayWithMortonIndex(box, nbParticles, particles1) ;
      sortArrayWithMortonIndex(box, nbParticles, particles2) ;
    }

  const FSize error = compareTwoArrays<FReal, FmaRWParticle<FReal,8,8>* >("TAG", nbParticles, particles1, particles2);

  //
  delete[] particles1 ;
  delete[] particles2 ;

  //
  return int(error);
}
