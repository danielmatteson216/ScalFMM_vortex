#ifndef FGROUPTOOLS_HPP
#define FGROUPTOOLS_HPP

#include <vector>

#include "ScalFmmConfig.h"
#include "Utils/FGlobal.hpp"
#include "Utils/FPoint.hpp"
#ifdef SCALFMM_USE_MPI
#include "Utils/FMpi.hpp"
#endif


namespace groupTree {
  // Structure for 1 particle
  template<typename FReal>
  struct particle_t {
    using position_t = FPoint<FReal>;
    position_t pos;
    FReal      phi;
    MortonIndex  morton_index;
    const auto& position() const {
      return pos;
    }
    const FPoint<FReal>& getPosition(){
      return pos;
    }
    const FPoint<FReal>& getPosition() const{
      return pos;
    }
    const FReal& physicalValue() const{
      return phi;
    }

    void fill(const position_t &inPos, const FReal &inPhyVal, MortonIndex & inMortonIndex){
      pos = inPos ; phi = inPhyVal ; morton_index = inMortonIndex ;
    }

    int weight() const { return 1;}
    MortonIndex getMorton() const{
      return morton_index;

    }
    friend constexpr MortonIndex morton_index(const particle_t& p) {
      return p.morton_index;
    }
  };


  //
  // param[in] FMpiComm
  // param[in] elapsedTime time on each processor
  // param[out]  minTime  the minimum time on each processor
  // param[out]  maxTime  the maximal time on each processor
  // param[out]  meanTime  the mean time on each processor
  //
  void timeAverage(const FMpi &FMpiComm, double &elapsedTime , double &minTime,
                   double &maxTime, double &meanTime)
  {
    double * allTimes = nullptr ;
    int myrank = FMpiComm.global().processId() , nprocs=FMpiComm.global().processCount() ;
    if(myrank == 0)
      {
        allTimes  = new double[nprocs] ;
      }
#ifdef SCALFMM_USE_MPI
    MPI_Gather(&elapsedTime,1, MPI_DOUBLE, allTimes, 1, MPI_DOUBLE,0 /* root*/,FMpiComm.global().getComm()) ;
#endif
    if(myrank == 0)
      {
        minTime = allTimes[0],  maxTime = allTimes[0], meanTime = allTimes[0] ;

        for (int i = 1 ; i < nprocs ; ++i)  {
            minTime   = std::min(minTime, allTimes[i]) ;
            maxTime   = std::max(maxTime, allTimes[i]) ;
            meanTime += allTimes[i] ;
          }
        meanTime /= nprocs;

      }
    FMpiComm.global().barrier() ;
  }
  template <class OCTREECLASS>
  void saveSolutionInFile(const std::string &fileName, const std::size_t&  NbPoints,
                           OCTREECLASS &tree) {
    using REALTYPE =  typename OCTREECLASS::FRealType ;
    FFmaGenericWriter<REALTYPE> writer(fileName) ;
    //
    REALTYPE * particles = new REALTYPE[8*NbPoints] ;
    memset(particles,0,8*NbPoints*sizeof(REALTYPE));
    FSize j = 0 ;
  #ifdef _VERBOSE_LEAF
    int countLeaf = 0, coutPart=0;
  #endif
    tree.forEachLeaf([&](typename OCTREECLASS::LeafClass* leaf){
      //
      // Input
      const REALTYPE*const posX = leaf->getTargets()->getPositions()[0];
      const REALTYPE*const posY = leaf->getTargets()->getPositions()[1];
      const REALTYPE*const posZ = leaf->getTargets()->getPositions()[2];
      const REALTYPE*const physicalValues = leaf->getTargets()->getPhysicalValues();
    //  const FVector<FSize>& indexes = leaf->getTargets()->getIndexes();
      //
      // Computed data
      const REALTYPE*const potentials = leaf->getTargets()->getPotentials();
      const REALTYPE*const forcesX = leaf->getTargets()->getForcesX();
      const REALTYPE*const forcesY = leaf->getTargets()->getForcesY();
      const REALTYPE*const forcesZ = leaf->getTargets()->getForcesZ();
      //
      //
      const FSize nbParticlesInLeaf = leaf->getTargets()->getNbParticles();
  #ifdef _VERBOSE_LEAF
      std::cout << "Leaf " << countLeaf << " Particles : [ " << coutPart << ", " <<coutPart+nbParticlesInLeaf -1 <<  " ] "  << nbParticlesInLeaf << std::endl;
      coutPart += nbParticlesInLeaf ; ++countLeaf;
  #endif
      for(FSize idxPart = 0 ; idxPart < nbParticlesInLeaf ; ++idxPart,j+=8){
       //   j = 8*indexes[idxPart];
         // j = 8*idxPart;
          particles[j]    = posX[idxPart] ;
          particles[j+1]  = posY[idxPart] ;
          particles[j+2]  = posZ[idxPart] ;
          particles[j+3]  = physicalValues[idxPart] ;
          particles[j+4]  = potentials[idxPart] ;
          particles[j+5]  =  forcesX[idxPart] ;
          particles[j+6]  =  forcesY[idxPart] ;
          particles[j+7]  =  forcesZ[idxPart] ;
        }
    });

    writer.writeHeader( tree.getBoxCenter(), tree.getBoxWidth() ,  NbPoints, sizeof(REALTYPE), 8) ;
    writer.writeArrayOfReal(particles,  8 , NbPoints);

    delete[] particles;
  }

  template< typename FReal, class GROUPTREE_T,class GROUPALGO_T, class OCTTREE_T>
  void checkCellTree(GROUPTREE_T &groupedTree, GROUPALGO_T & groupalgo, OCTTREE_T &treeCheck, const FReal &epsilon){
    //

    std::vector<bool>  OK(groupedTree.getHeight(),true) ;
    groupedTree.forEachMyCellWithLevel(
          [&](typename GROUPTREE_T::GroupSymbolCellClass_T*         gsymb ,
                       typename GROUPTREE_T::GroupCellUpClass_T*    gmul ,
                       typename GROUPTREE_T::GroupCellDownClass_T*  gloc ,
                       const int level)
    {
        if(groupalgo.isDataOwnedBerenger(gsymb->getMortonIndex(), level))
          {
            const auto * cell = treeCheck.getCell(gsymb->getMortonIndex(), level);
            if(cell == nullptr){
                std::cout << "[Empty] Error cell should exist " << gsymb->getMortonIndex() << "\n";
                OK[level] = false ;
              }
            else {
                FMath::FAccurater<FReal> diffUp;
                diffUp.add(cell->getMultipoleData().get(0), gmul->get(0), gmul->getVectorSize());
                if(diffUp.getRelativeInfNorm() > epsilon || diffUp.getRelativeL2Norm() > epsilon){
                    auto * data1 = gmul->get(0);
                    auto * data2 = cell->getMultipoleData().get(0);
                    for (int i = 0; i < gmul->getVectorSize(); ++i) {
                        std::cout << i << "  "<< data1[i] << "  seq "<< data2[i] <<std::endl;
                      }
                    std::cout << "[Up] Up is different at index " << gsymb->getMortonIndex() << " level " << level << " is " << diffUp << "\n";
                    OK[level] = false ;

                  }
                FMath::FAccurater<FReal> diffDown;
                diffDown.add(cell->getLocalExpansionData().get(0), gloc->get(0), gloc->getVectorSize());
                if(diffDown.getRelativeInfNorm() > epsilon || diffDown.getRelativeL2Norm() > epsilon){
                    std::cout << "[Down] Down is different at index " << gsymb->getMortonIndex() << " level " << level << " is " << diffDown << "\n";
                    OK[level] = false ;

                  }
              }
          }
      });
    for (std::size_t l = 0 ; l < OK.size(); ++l){
        std:: cout << "      Level ( " << l << " ) --> " << (OK[l] ? " Ok" : "Error " ) <<std::endl;
      }
    std:: cout << " checkCellTree --> done" <<std::endl;

  }
  template< typename FReal, class GROUPTREE_T,class GROUPALGO_T, class OCTTREE_T>
  void checkLeaves(GROUPTREE_T &groupedTree, GROUPALGO_T & groupalgo, OCTTREE_T &treeCheck, const FReal &epsilon){
    //
    FMath::FAccurater<FReal> potentialGlobalDiff;
    const int leafLevel  =  groupedTree.getHeight() - 1;
    bool  OK = true ;
    groupedTree.template forEachCellLeaf<typename GROUPTREE_T::LeafClass_T >(
          [&](typename GROUPTREE_T::GroupSymbolCellClass_T* gsymb ,
                       typename GROUPTREE_T::GroupCellUpClass_T*   /* gmul */,
                       typename GROUPTREE_T::GroupCellDownClass_T* /* gloc */,
                       typename GROUPTREE_T::LeafClass_T * leafTarget
                       )
    {

        if(groupalgo.isDataOwnedBerenger(gsymb->getMortonIndex(), leafLevel) )// just needed if we have  duplicated tree
          {
            const auto * targets = treeCheck.getLeafSrc(gsymb->getMortonIndex());
            if(targets == nullptr){
                std::cout << "[Empty] Error leaf should exist " << gsymb->getMortonIndex() << "\n";
                OK = false ;

              }
            else{
                const FReal*const gposX = leafTarget->getPositions()[0];
                const FReal*const gposY = leafTarget->getPositions()[1];
                const FReal*const gposZ = leafTarget->getPositions()[2];
                const FSize gnbPartsInLeafTarget = leafTarget->getNbParticles();
                const FReal*const gforceX = leafTarget->getForcesX();
                const FReal*const gforceY = leafTarget->getForcesY();
                const FReal*const gforceZ = leafTarget->getForcesZ();
                const FReal*const gpotential = leafTarget->getPotentials();

                const FReal*const posX = targets->getPositions()[0];
                const FReal*const posY = targets->getPositions()[1];
                const FReal*const posZ = targets->getPositions()[2];
                const FSize nbPartsInLeafTarget = targets->getNbParticles();
                const FReal*const forceX = targets->getForcesX();
                const FReal*const forceY = targets->getForcesY();
                const FReal*const forceZ = targets->getForcesZ();
                const FReal*const potential = targets->getPotentials();

                if(gnbPartsInLeafTarget != nbPartsInLeafTarget){
                    std::cout << "[Empty] Not the same number of particles at " << gsymb->getMortonIndex()
                              << " gnbPartsInLeafTarget " << gnbPartsInLeafTarget << " nbPartsInLeafTarget " << nbPartsInLeafTarget << "\n";
                    OK = false ;
                  }else{
                    FMath::FAccurater<FReal> potentialDiff;
                    FMath::FAccurater<FReal> fx, fy, fz;
                    for(FSize idxPart = 0 ; idxPart < nbPartsInLeafTarget ; ++idxPart){
                        if(gposX[idxPart] != posX[idxPart] || gposY[idxPart] != posY[idxPart] || gposZ[idxPart] != posZ[idxPart]){
                            std::cout << "[Empty] Not the same particle at " << gsymb->getMortonIndex() << " idx " << idxPart << " "
                                      << gposX[idxPart] << " " << posX[idxPart] << " " << gposY[idxPart] << " " << posY[idxPart]
                                         << " " << gposZ[idxPart] << " " << posZ[idxPart] << "\n";
                            OK = false ;
                          }
                        else{
                            potentialGlobalDiff.add(potential[idxPart], gpotential[idxPart]);
                            potentialDiff.add(potential[idxPart], gpotential[idxPart]);
                            fx.add(forceX[idxPart], gforceX[idxPart]);
                            fy.add(forceY[idxPart], gforceY[idxPart]);
                            fz.add(forceZ[idxPart], gforceZ[idxPart]);
                          }
                      }
                    if(potentialDiff.getRelativeInfNorm() > epsilon || potentialDiff.getRelativeL2Norm() > epsilon){
                        std::cout << " potentialDiff is different at index " << gsymb->getMortonIndex() << " is " << potentialDiff << "\n";
                        OK = false ;
                      }
                    if(fx.getRelativeInfNorm() > epsilon || fx.getRelativeL2Norm() > epsilon){
                        std::cout << " fx is different at index " << gsymb->getMortonIndex() << " is " << fx << "\n";
                        OK = false ;
                      }
                    if(fy.getRelativeInfNorm() > epsilon || fy.getRelativeL2Norm() > epsilon){
                        std::cout << " fy is different at index " << gsymb->getMortonIndex() << " is " << fy << "\n";
                        OK = false ;
                      }
                    if(fz.getRelativeInfNorm() > epsilon || fz.getRelativeL2Norm() > epsilon){
                        OK = false ;
                        std::cout << " fz is different at index " << gsymb->getMortonIndex() << " is " << fz << "\n";
                      }
                  }
              }
          }
      });
    std::cout << " potentialDiff  is " << potentialGlobalDiff << "\n";
    std:: cout << " checkLeaves --> " << (OK ? " Ok" : "Error " ) <<std::endl;

  }
}
#endif // FGROUPTOOLS_HPP
