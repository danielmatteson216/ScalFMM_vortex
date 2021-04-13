#ifndef _FGROUPTREE_CHECK_RESULTS_
#define _FGROUPTREE_CHECK_RESULTS_

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

#include "Utils/FGlobal.hpp"
#include "Utils/FAssert.hpp"
#include "Utils/FMath.hpp"
#include "Files/FFmaGenericLoader.hpp"
#include "Utils/FPoint.hpp"
// FBox
#include "Adaptive/FBox.hpp"
// Group linear tree
#include "GroupTree/Core/FGroupLinearTree.hpp"
//
#include "GroupTree/Core/FGroupTools.hpp"
//
//
// param[in] FMpiComm
// param[in] seqLoader
// param[in]  box
// param[in]  TreeHeight
// param[inout]  myParticles
//
template < class LOADER_T, typename PARTICLE_T, typename BOX_T >
void readAndSortAllParticle(LOADER_T & seqLoader, const  BOX_T & box,
                            std::vector<PARTICLE_T> &myParticles, const int TreeHeight ){

  using REAL= typename LOADER_T::dataType ;
  FAssertLF(seqLoader.isOpen());
  const FSize NbParticles   = seqLoader.getNumberOfParticles();
  //
  // Read File
  myParticles.clear() ;

  myParticles.resize(NbParticles) ;
  const std::size_t max_level = sizeof(PARTICLE_T::morton_index) * 8 / 3;

  for(FSize idxPart = 0 ; idxPart < NbParticles; ++idxPart){
      FPoint<REAL> pos ;
      REAL physicalValue ;
      seqLoader.fillParticle(&pos, &physicalValue);//Same with file or not
      //
       MortonIndex morton = inria::linear_tree::get_morton_index( pos, box, max_level);
      myParticles[idxPart].fill(pos,physicalValue,morton) ;
    }
  std::sort(myParticles.begin(), myParticles.end(), [&](const PARTICLE_T& a, const PARTICLE_T& b) {
      return (a.getMorton() < b.getMorton()  ) ;
    }
  );
  // Set the right MortonIndex
  for(FSize idxPart = 0 ; idxPart < NbParticles ; ++idxPart){
      myParticles[idxPart].morton_index = inria::linear_tree::get_morton_index(  myParticles[idxPart] .pos, box,
                                                                                 TreeHeight-1);
    }
}
//
// param[in] FMpiComm
// param[in] elapsedTime time on each processor
// param[out]  minTime  the minimum time on each processor
// param[out]  maxTime  the maximal time on each processor
// param[out]  meanTime  the mean time on each processor
//
template <typename PARTICLE, class REAl, typename OCTREECLASS1,
          typename OCTREECLASS2,class FmmClass1, class FmmClass2>
void checkWithDuplicatedTree( const int& myrank, const PARTICLE &arrayParticles,
                              OCTREECLASS1    & treeCheck,
                              FmmClass1       & algorithm,
                              OCTREECLASS2    & grouptree,
                              FmmClass2       & groupalgo,
                              const int &operationsToProceed,
                              const REAl& epsilon ) {
  //
  std::cout << "checkWithDuplicatedTree - nb part " <<  arrayParticles.size()  <<std::endl;

  //  Compute a sequential FMM
  algorithm.execute(operationsToProceed);
  //
//  std::string fileName("output-Let-") ;
//  fileName += std::to_string(myrank) + ".fma" ;
//  groupTree::saveSolutionInFile(fileName, arrayParticles.size() ,treeCheck) ;

  groupTree::checkCellTree(grouptree, groupalgo, treeCheck, epsilon) ;
  groupTree::checkLeaves(grouptree, groupalgo, treeCheck, epsilon) ;

  std::cout << "Comparing is over" << std::endl;
}


#endif
