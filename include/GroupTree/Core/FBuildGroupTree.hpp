// ==== CMAKE =====
// @FUSE_MPI
// @FUSE_STARPU
//
#ifndef FBuildGroupTree
#define FBuildGroupTree

#include <vector>
#include <string>
#include "Utils/FGlobal.hpp"

// include algo for linear tree
#include "inria/algorithm/distributed/mpi.hpp"
#include "inria/linear_tree/node.hpp"
#include "inria/linear_tree/balance_tree.hpp"
// tree class
#include "GroupTree/Core/FGroupTree.hpp"
// symbolic data
#include "Components/FSymbolicData.hpp"
// GroupParticleContainer
#include "GroupTree/Core/FP2PGroupParticleContainer.hpp"
// file loader
#include "Files/FMpiFmaGenericLoader.hpp"
// FBox
#include "Adaptive/FBox.hpp"
// Group linear tree
#include "GroupTree/Core/FGroupLinearTree.hpp"
// Function for GroupLinearTree
#include "GroupTree/Core/FDistributedGroupTreeBuilder.hpp"
//
#include "GroupTree/Core/FGroupTools.hpp"

// To construct either the duplicated Octree or the LET
#include "Utils/FLeafBalance.hpp"

namespace groupTree {
  //
  // @param[in]    mpi_comm   the MPI communicator
  // @param[inout] myParticleslocal array of particles on my node. On output the array is sorted
  // @param[in]    total number of particles in the simulation
  // @param[in]    box  size of the simulation box
  // @param[in]    TreeHeight    Height of the tree
  // @param[inout]    localGroupTree  the LET of the octree
  // @param[out]    m_idx_distribution  Distribution of the leaves on the processors
  // @param[out]   nb_blocks
  template <class LOADER, class particleType , class OCTREEGRPOUPCLASS, class  FBox>
  void buildLetTree( inria::mpi::communicator & mpi_comm, LOADER& loader, std::vector<particleType> &myParticles,
                     const  FBox /*<FPoint<FReal>> */box,
                     const int TreeHeight, const int groupSize,
                     OCTREEGRPOUPCLASS * &localGroupTree,
                     std::vector<MortonIndex> &m_idx_distribution, int & nb_blocks
                     ){
    //
    const std::size_t max_level = sizeof(MortonIndex) * 8 / 3;
    const FSize localNumberOfParticles = loader.getMyNumberOfParticles() ;

    myParticles.resize(localNumberOfParticles) ;

    // iterate on all of my particles
    for(FSize idxPart = 0; idxPart < static_cast<FSize>(localNumberOfParticles );++idxPart){
        particleType   tmp;
        // get the current particles
        loader.fillParticle(&tmp.pos,&tmp.phi);
        // set the morton index of the current particle at the max_level
        tmp.morton_index = inria::linear_tree::get_morton_index(tmp.pos, box, max_level);
        // set the weight of the particle
      //  tmp.phi = 0.1;
        // add the particle to my vector of particle
        myParticles[idxPart].fill(tmp.pos, tmp.phi,tmp.morton_index);
      }
    // Now i have all of my particles in a vector, they all have a morton index
    // now we will sort them
    inria::sort(mpi_comm,myParticles, [](const auto& p1, const auto& p2) {
        return p1.morton_index < p2.morton_index;
      });
#ifdef FMM_VERBOSE_DATA_DISTRIUTION
    const FSize totalNumberOfParticles = loader.getNumberOfParticles() ;
    std::cout << " I have "         << myParticles.size() << " particles ..." << std::endl;
    std::cout << "For a total of "
              << totalNumberOfParticles << " particles ..." << std::endl;
#endif
    // create the linear tree
    // a linear tree is a tree, with only the leaf
    int level = TreeHeight -1 ;
    auto linear_tree = inria::linear_tree::create_balanced_linear_tree_at_level(
          mpi_comm,
          level,
          box,
          myParticles);

    // create GroupLinearTree
    FGroupLinearTree<typename decltype(linear_tree)::value_type>  group_linear_tree{mpi_comm};
    group_linear_tree.create_local_group_linear_tree(  &linear_tree, groupSize );

    // group_linear_tree.print_info_tree() ;

    // Redistribute the particle according to the linear tree
    // Redistribution of particles
    inria::linear_tree::redistribute_particles(mpi_comm,
                                               linear_tree,
                                               myParticles);

    // Now we need to modify the morton index of of all particle to
    // have the morton index at TreeHeight-1
#pragma omp parallel for shared(myParticles)
    for(unsigned i = 0 ; i < myParticles.size(); ++i){
        myParticles[i].morton_index = inria::linear_tree::get_morton_index(myParticles[i].pos, box, level);
      }

    // Now we need to share the particle distribution to build the GroupTree
    group_linear_tree.set_index_particle_distribution(myParticles);

    // Now i can declare my groupTree
    // it's a empty instance of the FGroupTree
    FReal width = std::max(box.width(0) , std::max(box.width(1) ,box.width(2) )) ;
    //   using test = typename std::remove_pointer<typename std::remove_reference<decltype(localGroupTree)>::type >::type;
    // //   std::cout << "&&&&&"<<typeid(test).name() <<std::endl;
    //
    localGroupTree = new OCTREEGRPOUPCLASS (TreeHeight,groupSize, box.center(), box.c1() /* corner*/,
                                            width, width/FReal(1<<(TreeHeight-1)));
    // Now i can fill the localGroupTree
    localGroupTree->create_tree(group_linear_tree,myParticles);
#ifdef FMM_VERBOSE_DATA_DISTRIUTION
    localGroupTree->printInfoBlocks();
#endif
    // get the index particle distribution (needed by the algorithm)

    m_idx_distribution = group_linear_tree.get_index_particle_distribution_implicit();
    //  for(int i = 0 ; i < mpi_comm.size() ;++i)
    //    m_idx_distribution[2*i] += 1;
    nb_blocks = dstr_grp_tree_builder::set_cell_group_global_index(*localGroupTree,mpi_comm);
    // now we create the LET
    localGroupTree->create_LET(group_linear_tree);

 //   std::cout << " End buildLetTree function " << std::endl;
  }
  // BuilddMortonDistributionForCGroupCellInTree
  //
  // @param[in]  parallelManager            The Height of the octree
  // @param[in]  mortonLeaves            The Height of the octree
  // @param[in]  TreeHeight            The Height of the octree
  // @param[in]  groupSize
  // @param[in]  MortonIndexDistribution The Morton distribution at the leaf level
  // @param[in]  nodeRepartition
  // @param[out] sizeForEachGroup         For each level give the size of all group cell in the process
  void BuilddMortonDistributionForCGroupCellInTree(const FMpi &parallelManager,  std::vector<MortonIndex> &mortonLeaves,
                                                   const int & TreeHeight, const int groupSize,
                                                   const std::vector<MortonIndex> &MortonIndexDistribution,
                                                   std::vector<std::vector<std::vector<MortonIndex>>>& nodeRepartition,
                                                   std::vector< std::vector<int>> &sizeForEachGroup ){
    //
    const int nproc = parallelManager.global().processCount() ;
    //
    // Build the groupe size of all groups in the Tree (sequential one)
    //
#ifdef FMM_VERBOSE_DATA_DISTRIUTION
    std::cout << "Morton distribution inside BuilddMortonDistributionForCGroupCellInTree " <<std::endl;
    for (auto v : MortonIndexDistribution)
      std::cout << "  " << v ;
    std::cout << std::endl;
#endif
    int processId ;
    for( processId = 0; processId < nproc; ++processId)
      {
        FSize size_last, countGroup;
        // pas de +1 si on ne commence pas à 0
        FSize leafOnProcess = MortonIndexDistribution[2*processId+1] - MortonIndexDistribution[2*processId]  ;
        size_last = leafOnProcess%groupSize;
        countGroup = (leafOnProcess - size_last)/groupSize;
        for(int i = 0; i < countGroup; ++i)
          sizeForEachGroup[TreeHeight-1].push_back(groupSize);
        if(size_last > 0)
          sizeForEachGroup[TreeHeight-1].push_back((int)size_last);
      }
    //
    //Pour chaque niveau calcul de la taille des groupe
    for(int idxLevel = TreeHeight - 2; idxLevel >= 0; --idxLevel)
      {
        processId = 0;
        int countCellsInTheGroup = 0;
        MortonIndex previousMortonCell = -1;

//        std::cout << "Compute Level " << idxLevel << std::endl;
        for(std::size_t idxLeaf = 0; idxLeaf < mortonLeaves.size(); ++idxLeaf)
          {
            MortonIndex mortonCell = (mortonLeaves[idxLeaf]) >> (3*(TreeHeight - 1 - idxLevel));
            if(mortonCell <= nodeRepartition[idxLevel][processId][1]) //Si l'indice est dans le working interval
              {
                if(mortonCell != previousMortonCell) //Si c'est un nouvelle indice
                  {
                    ++countCellsInTheGroup; //On le compte dans le groupe
                    previousMortonCell = mortonCell;
                    if(countCellsInTheGroup == groupSize) //Si le groupe est plein on ajoute le compte
                      {
                        sizeForEachGroup[idxLevel].push_back(groupSize);
                        countCellsInTheGroup = 0;
                      }
                  }
              }
            else //Si l'on change d'interval de process on ajoute ce que l'on a compté
              {
                if(countCellsInTheGroup > 0)
                  sizeForEachGroup[idxLevel].push_back(countCellsInTheGroup);
                countCellsInTheGroup = 1;
                previousMortonCell = mortonCell;
                ++processId;
              }
          }
        if(countCellsInTheGroup > 0)
          sizeForEachGroup[idxLevel].push_back(countCellsInTheGroup);
        //
#ifdef FMM_VERBOSE_DATA_DISTRIUTION
        // Print sizeForEachGroup at the current level
        for( auto v : sizeForEachGroup[idxLevel])
          std::cout << "  "<< v ;
        std::cout << std::endl;
#endif
      }


  }
  // Build Node distribution for all LEvel starting with Leaf Distribution
  //
  void createNodeRepartition(std::vector<MortonIndex> distributedMortonIndex,
                             std::vector<std::vector<std::vector<MortonIndex>>>& nodeRepartition,
                             int nproc, int treeHeight) {
    //
    nodeRepartition.resize(treeHeight, std::vector<std::vector<MortonIndex>>(nproc, std::vector<MortonIndex>(2)));
    for(int node_id = 0; node_id < nproc; ++node_id){
        nodeRepartition[treeHeight-1][node_id][0] = distributedMortonIndex[node_id*2];
        nodeRepartition[treeHeight-1][node_id][1] = distributedMortonIndex[node_id*2+1];
      }
    for(int idxLevel = treeHeight - 2; idxLevel >= 0  ; --idxLevel){
        nodeRepartition[idxLevel][0][0] = nodeRepartition[idxLevel+1][0][0] >> 3;
        nodeRepartition[idxLevel][0][1] = nodeRepartition[idxLevel+1][0][1] >> 3;
        for(int node_id = 1; node_id < nproc; ++node_id){
            nodeRepartition[idxLevel][node_id][0] = FMath::Max(nodeRepartition[idxLevel+1][node_id][0] >> 3, nodeRepartition[idxLevel][node_id-1][0]+1); //Berenger phd :)
            nodeRepartition[idxLevel][node_id][1] = nodeRepartition[idxLevel+1][node_id][1] >> 3;
          }
      }
  }

  //
  //
  // @param[in]    mpi_comm   the MPI communicator
  // @param[in]    filename  Particles name file
  // @param[in]    option   option to build the groupTree
  //                     1 we use a given Morton distribution
  //                     2 we build the Morton distribution
  // @param[out] myParticleslocal array of particles on my node. On output the array is sorted
  // @param[in]    box  size of the simulation box
  // @param[in]    TreeHeight    Height of the tree
  // @param[in]    localGroupTree  the LET of the octree
  // @param[inout]    localGroupTree  the LET of the octree
  // @param[inout]    m_idx_distribution  Distribution of the leaves on the processors
  // @param[out]    nb_blocks
  template <class PARTICLE_T , class OCTREEGRPOUPCLASS>
  void buildDuplicatedTree( const FMpi &parallelManager, const int option, const std::string &filename,
                            std::vector<PARTICLE_T> &myParticles, const  FBox<FPoint<FReal>>& box,
                            const int TreeHeight, const int groupSize, OCTREEGRPOUPCLASS * &GroupTree,
                            std::vector<MortonIndex> &MortonIndexDistribution ,int & nb_block)
  {

    //
    //loader
//    std::cout << "Opening : " << filename << " ...";
    FFmaGenericLoader<FReal>  loader(filename);
    FAssertLF(loader.isOpen());
 //   std::cout << " done." << std::endl;
    const FSize totalNbParticles   = loader.getNumberOfParticles();
    //
    const std::size_t max_level = sizeof(PARTICLE_T::morton_index) * 8 / 3;
    //
    myParticles.resize(totalNbParticles);
    for(FSize idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
        FPoint<FReal> pos ;
        FReal physicalValue ;
        loader.fillParticle(&pos, &physicalValue);//Same with file or not
        //
   //     physicalValue = 0.1 ;
        MortonIndex morton = inria::linear_tree::get_morton_index( pos, box, max_level);
        myParticles[idxPart].fill(pos,physicalValue,morton) ;
      }
    std::sort(myParticles.begin(), myParticles.end(), [&](const PARTICLE_T& a, const PARTICLE_T& b) {
        return (a.getMorton() < b.getMorton()  ) ;
      }
    );
    //
    FP2PParticleContainer<FReal> allParticles;

    // Set the right MortonIndex
  //  MortonIndex mm = 0 ;
    for(FSize idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
        myParticles[idxPart].morton_index = inria::linear_tree::get_morton_index(  myParticles[idxPart] .pos, box,
                                                                                   TreeHeight-1);
  //      mm = std::max(mm,myParticles[idxPart].morton_index );
        allParticles.push(myParticles[idxPart].getPosition(), myParticles[idxPart].physicalValue() );
      }
    // Create the linear tree
    // a linear tree is a tree, with only the leaf
    // Build a vector of MortonIndex at Leaf level from particles
    //
    std::size_t nbLeaves = 1 , pos=0;
    MortonIndex previousMorton = myParticles[0].morton_index;
    for(std::size_t idxPart = 1 ; idxPart < myParticles.size(); ++idxPart){
        if(previousMorton != myParticles[idxPart].morton_index){
            previousMorton = myParticles[idxPart].morton_index ;
            ++nbLeaves ;
          }
      }
 //   std::cout<< "Number of leaves" << nbLeaves <<std::endl ;
    std::vector<MortonIndex> mortonLeaves(nbLeaves,-1) ;

    previousMorton    = myParticles[0].morton_index;
    mortonLeaves[pos] = myParticles[0].morton_index;

    for(std::size_t idxPart = 1 ; idxPart < myParticles.size(); ++idxPart){
        if(previousMorton != myParticles[idxPart].getMorton() ){
            ++pos ;
            previousMorton    = myParticles[idxPart].morton_index ;
            mortonLeaves[pos] = myParticles[idxPart].morton_index ;
          }
      }
    const int nproc = parallelManager.global().processCount() ; //
    if(option >1 ) {
        std::cout << " Construct the distribution used in Beregnger's thesis "<< std::endl;
        FLeafBalance balancer;
        MortonIndexDistribution.clear() ;
        //
        // Build the Morton index as in Berenger's thesis
        //Calcul du working interval au niveau des feuilles
        previousMorton = -1;
        int countLeaf  = 0;
        int processId  = 0;
        FSize leafOnProcess = balancer.getRight(nbLeaves, nproc, 0) - balancer.getLeft(nbLeaves, nproc, 0);
        std::cout << "  leafOnProcess  " << leafOnProcess <<  " empty? " << MortonIndexDistribution.empty() << "  "  << MortonIndexDistribution.size() <<std::endl;
        MortonIndexDistribution.push_back(previousMorton);
        for(FSize idxPart = 0 ; idxPart < loader.getNumberOfParticles()  ; ++idxPart)
          {
            if(myParticles[idxPart].morton_index != previousMorton)
              {
                previousMorton = myParticles[idxPart].morton_index ;
                ++countLeaf;
                if(countLeaf == leafOnProcess)
                  {
                    ++processId;
                    if (processId < nproc){
                      leafOnProcess = balancer.getRight(nbLeaves, nproc, processId) - balancer.getLeft(nbLeaves, nproc, processId);
                      MortonIndexDistribution.push_back(previousMorton);
                      MortonIndexDistribution.push_back(previousMorton);
                      countLeaf = 0;
                      }
                  }
              }
          }
        MortonIndexDistribution.push_back(myParticles[loader.getNumberOfParticles() - 1].morton_index) ;
        //
      }
    // otherwise we use the given Morton distribution
#ifdef FMM_VERBOSE_DATA_DISTRIUTION
    std::cout << "    Morton distribution to build the duplicated tree " <<MortonIndexDistribution.size() << " "<<std::endl<<std::flush;
    for (auto v : MortonIndexDistribution)
      std::cout << "  " << v ;
    std::cout << std::endl;
#endif
    //////////////////////////////////////////////////////////////////////////
    std::vector< std::vector<std::vector<MortonIndex>>> nodeRepartition;
    std::vector< std::vector<int>>                          sizeForEachGroup(TreeHeight);
    createNodeRepartition(MortonIndexDistribution, nodeRepartition, nproc, TreeHeight) ;
#ifdef FMM_VERBOSE_DATA_DISTRIUTION
    for ( std::size_t idLevel=0;  idLevel< nodeRepartition.size() ; ++idLevel){
        std::cout << "  nodeRepartition at level " << idLevel << std::endl ;
        for ( std::size_t procID=0 ;  procID<  nodeRepartition[idLevel].size();  ++procID){
            std::cout << "  n  proc( " << procID << "  ) " <<
                         " [ " << nodeRepartition[idLevel][procID][0] << ", "
                      << nodeRepartition[idLevel][procID][1] <<" ]" <<std::endl ;
          }
      }
#endif

    BuilddMortonDistributionForCGroupCellInTree(parallelManager,mortonLeaves,TreeHeight,groupSize,
                                                MortonIndexDistribution,nodeRepartition,sizeForEachGroup ) ;
#ifdef FMM_VERBOSE_DATA_DISTRIUTION
    //
    // Print group size per level
    std::cout << std::endl<< "  Group size at the leaf level " << std::endl ;
    int totalLeaves = 0 ;
    for ( std::size_t idLevel=2;  idLevel< sizeForEachGroup.size() ; ++idLevel){
        std::cout << "  Group size at level " << idLevel << std::endl ;
        totalLeaves = 0 ;
        for ( auto v : sizeForEachGroup[idLevel]){
            totalLeaves += v;
            std::cout << "   " << v ;
          }
        std::cout << std::endl ;std::cout << " Total number of leaves: " <<totalLeaves << std::endl;
      }
#endif
    //
    GroupTree  = new OCTREEGRPOUPCLASS (TreeHeight, loader.getBoxWidth(), loader.getCenterOfBox(),
                                        groupSize, &allParticles, sizeForEachGroup, true);
    //
    //
  }
}

#endif
