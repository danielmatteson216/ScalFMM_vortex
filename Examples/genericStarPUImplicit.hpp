// -*-c++-*-
// ==== CMAKE =====
// @FUSE_BLAS
// @FUSE_MPI
// @FUSE_STARPU
// ================
//
#include "Utils/FGlobal.hpp"
// parameters
#include "Utils/FParameters.hpp"
#include "Utils/FParameterNames.hpp"

// include algo for linear tree
#include "inria/algorithm/distributed/mpi.hpp"
#include "inria/linear_tree/balance_tree.hpp"
// tree class
#include "GroupTree/Core/FGroupTree.hpp"
// symbolic data
#include "Components/FSymbolicData.hpp"
//

// GroupParticleContainer
#include "GroupTree/Core/FP2PGroupParticleContainer.hpp"
// file loader
#include "Files/FMpiFmaGenericLoader.hpp"
#include "Files/FTreeMpiCsvSaver.hpp"
// FBox
#include "Adaptive/FBox.hpp"
// Group linear tree
#include "GroupTree/Core/FGroupLinearTree.hpp"
// Function for GroupLinearTree
#include "GroupTree/Core/FDistributedGroupTreeBuilder.hpp"
//
// Algorithm include
#include "GroupTree/StarPUUtils/FStarPUKernelCapacities.hpp"
#include "GroupTree/StarPUUtils/FStarPUCpuWrapper.hpp"
#include "GroupTree/Core/FGroupTaskStarpuImplicitAlgorithm.hpp"
//
// To construct either the duplicated Octree or the LET
//
static const int ORDER  = 6;
using FReal             = double;

#include "GroupTree/Core/FBuildGroupTree.hpp"
//For validation
#include "GroupTree/Core/FGroupTools.hpp"
#include "GroupTree/Core/FCheckResults.hpp"
#include "Components/FSimpleLeaf.hpp"
#include "Core/FFmmAlgorithm.hpp"
// Four output
#include "json.hpp"


//
//   1/r kernel
using MatrixKernelClass = FInterpMatrixKernelR<FReal> ;
//
// definition of the common tree structure
using CellClass           = FInterpolationCell<FReal, ORDER>;
using GroupCellUpClass    = typename CellClass::multipole_t;
using GroupCellDownClass  = typename CellClass::local_expansion_t;
using GroupCellSymbClass  = FSymbolicData;
using GroupContainerClass = FP2PGroupParticleContainer<FReal>;
using GroupOctreeClass    = FGroupTree<FReal,GroupCellSymbClass,GroupCellUpClass,
                                                                   GroupCellDownClass, GroupContainerClass, 1, 4, FReal>;

// definition of algorithm structure
using GroupKernelClass    = FStarPUAllCpuCapacities<FInterpolationKernel<FReal, CellClass,
                                                                         GroupContainerClass,MatrixKernelClass,ORDER>>;

using GroupCpuWrapper     = FStarPUCpuWrapper<
typename GroupOctreeClass::CellGroupClass, CellClass, GroupKernelClass,
typename GroupOctreeClass::ParticleGroupClass, GroupContainerClass> ;
//
using GroupAlgorithm      = FGroupTaskStarPUImplicitAlgorithm<GroupOctreeClass,
typename GroupOctreeClass::CellGroupClass,GroupKernelClass,
typename GroupOctreeClass::ParticleGroupClass, GroupCpuWrapper>;
//////////////////////////////////////////////////////////////////


int main(int argc, char *argv[]) {
  // Parameter definition
  const FParameterNames LocalOptionBlocSize { {"-bs"},"The size of the block of the blocked tree"};
  const FParameterNames LocalOptionValidate { {"-check-result"}, "To compare with direct computation"};
  const FParameterNames LocalOptionBuildTree { {"-tree"}, "To compare with direct computation 0 let, 1 duplicate tree (let distribution) 2 duplicate tree "};
  const std::string TreeBuilderOption[3]={"Let tree", "Duplicated tree with Let distribution", "Duplicated tree"};
  // Parameter help
  FHelpDescribeAndExit(argc, argv,
                       "Test the blocked tree created with linear tree." ,
                       FParameterDefinitions::OctreeHeight,
                       FParameterDefinitions::InputFile,
                       FParameterDefinitions::OutputFile,
                       LocalOptionBlocSize,
                       LocalOptionValidate, LocalOptionBuildTree);

  // Get parameters
  // Get the groupSize
  const int groupSize =            FParameters::getValue(argc,argv,LocalOptionBlocSize.options, 250);
  // Get the file input
  const char* const filename       =
      FParameters::getStr(argc,argv,FParameterDefinitions::InputFile.options, "../Data/test20k.fma");
  // Get the treeHeight
  const unsigned int TreeHeight    =
      FParameters::getValue(argc, argv, FParameterDefinitions::OctreeHeight.options, 5);
  const int optionBuildTree = FParameters::getValue(argc, argv, LocalOptionBuildTree.options,0) ;


  // Init MPI communicator
  // Initialisation MPI Berenger
  FMpi FMpiComm;
  // Initialisation MPI Quentin
  inria::mpi::communicator mpi_comm(FMpiComm.global().getComm());
  // Show job information
  std::cout << "JOB INFORMATION " << std::endl;
  std::cout << "File name:    " << filename              << std::endl;
  std::cout <<  "TreeHeight:  "   <<   TreeHeight    << std::endl;
  std::cout <<  "Block size:  "   <<    groupSize   << std::endl;
  std::cout <<  "Tree type:   "    <<  TreeBuilderOption[optionBuildTree] << std::endl;
  std::cout << "------------------------------------------" << std::endl;
  FTic time;
  // Use FMpiFmaGenericLoader to read the box simulation size
  std::cout << "Opening : " << filename << " ...";
  FMpiFmaGenericLoader<FReal> loader(filename, FMpiComm.global());
  std::cout << " done." << std::endl;

  // define a box, used in the sort
  const FBox<FPoint<FReal>> box{loader.getBoxWidth(),loader.getCenterOfBox()};
  FReal width = std::max(box.width(0) , std::max(box.width(1) ,box.width(2) )) ;
  //
  // The group tree used for the computation
  GroupOctreeClass * computeOnGroupTree = nullptr ;
  //
  ///////////////////////////////////////////////////////////////////////////////
  //          Build Let or duplicated tree
  ///////////////////////////////////////////////////////////////////////////////
  //
  std::string title   ;
  int nb_block ;
  std::vector<MortonIndex> mortonCellDistribution ;
  // vector to stock all particles
  std::vector<groupTree::particle_t<FReal> > myParticles ;
  //
  // define the max level to sort particle

  std::string octreeType;
  if(optionBuildTree <= 1 ){
      title = "Distribution LETGroupTree in ";
      octreeType = "Let group"   ;
      //
      GroupOctreeClass *localGroupTree = nullptr;
      //
      groupTree::buildLetTree(mpi_comm, loader, myParticles,
                   box, TreeHeight, groupSize, localGroupTree, mortonCellDistribution ,nb_block);
      computeOnGroupTree   = localGroupTree ;
    }
  //
  if(optionBuildTree > 0 ){
      title ="duplicate GroupTree in ";
      octreeType = "duplicate group"   ;
      //
      GroupOctreeClass *fullGroupTree = nullptr;
      //
      groupTree::buildDuplicatedTree( FMpiComm, optionBuildTree, filename, myParticles, box, TreeHeight,
                           groupSize, fullGroupTree,mortonCellDistribution,nb_block);
      computeOnGroupTree   = fullGroupTree ;
      nb_block =0;
    }
  time.tac();
#ifdef FMM_VERBOSE_DATA_DISTRIUTION
  computeOnGroupTree->printInfoBlocks();
#endif
  std::cout << title <<  mortonCellDistribution.size() << std::endl;
  for ( auto v :  mortonCellDistribution)
    std::cout << "  " << v     ;
  std::cout << std::endl;
  std::cout << " nb_block: " << nb_block <<std::endl;
  std::cout << " Creating GroupTree in " << time.elapsed() << "s)." << std::endl;
  //
  ///////////////////////////////////////////////////////////////////////////////
  //          Computation part
  ///////////////////////////////////////////////////////////////////////////////
  //
  // define the operation to proceed
  //  FFmmNearField only Near field
  //  FFmmFarField  only Far field
  //  FFmmNearAndFarFields  full FMM
  //  By operator FFmmP2P| FFmmP2M | |  FFmmM2M  FFmmM2L | FFmmL2L | FFmmL2P
  const unsigned int operationsToProceed =   FFmmP2P |  FFmmP2M |  FFmmM2M |  FFmmM2L | FFmmL2L | FFmmL2P ;

  const MatrixKernelClass MatrixKernel;
  GroupKernelClass groupkernel(TreeHeight, width, box.center() , &MatrixKernel);
  GroupAlgorithm groupalgo(computeOnGroupTree,&groupkernel, mortonCellDistribution,nb_block);
  // wait all proc
  FTic timerExecute;
  FMpiComm.global().barrier();    // Synchronization for timer
  // start new timer
  timerExecute.tic();
  //  starpu_fxt_start_profiling();
  groupalgo.execute(operationsToProceed);
  computeOnGroupTree->printInfoBlocks();

  FMpiComm.global().barrier();
  //   starpu_fxt_stop_profiling();
  timerExecute.tac();
  auto timeElapsed = timerExecute.elapsed();
  // print times
  double minTime,maxTime,meanTime ;
  groupTree::timeAverage(FMpiComm, timeElapsed, minTime, maxTime, meanTime) ;
  std::cout <<  " FMM time (in sec.)  on node: " << timeElapsed
             << " min " << minTime << " max " << maxTime
             << " mean " << meanTime << std::endl;
  //
  ///////////////////////////////////////////////////////////////////////////////
  //          Extraction des resultats
  ///////////////////////////////////////////////////////////////////////////////
  //
  nlohmann::json result;
  std::string name = std::to_string(TreeHeight);
  name += "_" + std::to_string(groupSize)+"_"+std::to_string(loader.getNumberOfParticles()) + ".json";
  result["TreeHeight"] = TreeHeight;
  result["GroupSize"]  = groupSize;
  result["Filename"]   = filename;
  result["NbParticle"] = loader.getNumberOfParticles();
  result["Octree"]     = octreeType;
  result["Algorithm"]["time"] = timeElapsed;
  result["Algorithm"]["mean"] = meanTime;
  result["Algorithm"]["min"]  = minTime;
  result["Algorithm"]["max"]  = maxTime;
  result["kernel"]            = interpolationKernel ;
  std::ofstream out(name);
  out << result << std::endl;
  //
  ///////////////////////////////////////////////////////////////////////////////
  //          Validation
  ///////////////////////////////////////////////////////////////////////////////
  //
  // Validate the result
  if(FParameters::existParameter(argc, argv, LocalOptionValidate.options) == true){
      // Check the result with a previous computation
      // The resuls are stored in the files

      typedef FP2PParticleContainer<FReal>         ContainerClass;
      typedef FSimpleLeaf<FReal, ContainerClass >  LeafClass;
      using OctreeClass = FOctree<FReal, CellClass,ContainerClass,LeafClass>       ;
      using KernelClass = FInterpolationKernel<FReal,CellClass,ContainerClass,MatrixKernelClass,ORDER>    ;
      using FmmClass    = FFmmAlgorithm<OctreeClass,CellClass,ContainerClass,KernelClass,LeafClass> ;
      const int SubTreeHeight=3;
      OctreeClass treeCheck(TreeHeight, SubTreeHeight,width,box.center());
      const FReal epsilon = 1E-10;
      KernelClass kernels(TreeHeight, width, box.center(), &MatrixKernel);
      //
      if(optionBuildTree > 0 ){
          for(std::size_t idxPart = 0 ; idxPart < myParticles.size() ; ++idxPart){
              // put in tree
              treeCheck.insert(myParticles[idxPart].getPosition(),
                               myParticles[idxPart].physicalValue());
            }
        }
      else {
          //  We introduce a sequetial loader to read the particles from the file
          // the particles are sort as in the Let tree construction
          // then we insert the particles in an FOctree (not a grouped one)
          FFmaGenericLoader<FReal>  seqLoader(filename);
          readAndSortAllParticle(seqLoader, box,  myParticles, TreeHeight ) ;
          for(std::size_t idxPart = 0 ; idxPart < myParticles.size() ; ++idxPart){
              // put in tree
              treeCheck.insert(myParticles[idxPart].getPosition(),
                               myParticles[idxPart].physicalValue());
            }
        }
      FmmClass algorithm(&treeCheck, &kernels);

      checkWithDuplicatedTree(FMpiComm.global().processId(), myParticles,treeCheck, algorithm,
                                 *computeOnGroupTree,groupalgo,operationsToProceed,epsilon );

  }
  //
  //  Store the output position,rho, Potential, Force in file
  //  Only binary write is available
  if(FParameters::existParameter(argc, argv, FParameterDefinitions::OutputFile.options)){
      std::string outputFile(FParameters::getStr(argc,argv,FParameterDefinitions::OutputFile.options,   "outputMPI.bfma"));

      FMpiFmaGenericWriter<FReal> paraWriter(outputFile,FMpiComm);
      paraWriter.writeDistributionOfParticlesFromGroupedOctree(*computeOnGroupTree,loader.getNumberOfParticles(),mortonCellDistribution);
    }
  return 0;
}
