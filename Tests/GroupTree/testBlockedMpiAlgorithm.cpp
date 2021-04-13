
// Keep in private GIT
// @FUSE_MPI
// @FUSE_STARPU

#include "Utils/FGlobal.hpp"
#include "Utils/FMpi.hpp"

#include "GroupTree/Core/FGroupTree.hpp"

#include "Components/FSimpleLeaf.hpp"
#include "Components/FSymbolicData.hpp"
#include "Containers/FVector.hpp"


#include "Utils/FMath.hpp"
#include "Utils/FMemUtils.hpp"
#include "Utils/FParameters.hpp"

#include "Files/FRandomLoader.hpp"

#include "GroupTree/Core/FGroupTaskStarpuMpiAlgorithm.hpp"

#include "GroupTree/Core/FP2PGroupParticleContainer.hpp"
#include "GroupTree/Core/FGroupTaskAlgorithm.hpp"

#include "Utils/FParameterNames.hpp"

#include "Components/FTestParticleContainer.hpp"
#include "Components/FTestKernels.hpp"
#include "Components/FTestCell.hpp"
#include "GroupTree/TestKernel/FGroupTestParticleContainer.hpp"

#include "Utils/FLeafBalance.hpp"
#include "Files/FMpiTreeBuilder.hpp"

#include "Core/FFmmAlgorithm.hpp"
#include "Containers/FCoordinateComputer.hpp"

#include "GroupTree/StarPUUtils/FStarPUKernelCapacities.hpp"
#include "GroupTree/StarPUUtils/FStarPUCpuWrapper.hpp"




int main(int argc, char* argv[]){
    const FParameterNames LocalOptionBlocSize {
        {"-bs"},
        "The size of the block of the blocked tree"
    };
    FHelpDescribeAndExit(argc, argv, "Test the blocked tree by counting the particles.",
                         FParameterDefinitions::OctreeHeight,
                         FParameterDefinitions::NbParticles,
                         LocalOptionBlocSize);
    typedef double FReal;
    // Initialize the types
    using GroupCellClass     = FTestCell;
    using GroupCellUpClass   = typename FTestCell::multipole_t;
    using GroupCellDownClass = typename FTestCell::local_expansion_t;
    using GroupCellSymbClass = FSymbolicData;

    typedef FGroupTestParticleContainer<FReal>                                     GroupContainerClass;
    typedef FGroupTree< FReal, GroupCellSymbClass, GroupCellUpClass, GroupCellDownClass,
            GroupContainerClass, 0, 1, long long int>  GroupOctreeClass;
    typedef FStarPUAllCpuCapacities<FTestKernels< GroupCellClass, GroupContainerClass >>  GroupKernelClass;
    typedef FStarPUCpuWrapper<typename GroupOctreeClass::CellGroupClass, GroupCellClass, GroupKernelClass, typename GroupOctreeClass::ParticleGroupClass, GroupContainerClass> GroupCpuWrapper;
    typedef FGroupTaskStarPUMpiAlgorithm<GroupOctreeClass, typename GroupOctreeClass::CellGroupClass, GroupKernelClass, typename GroupOctreeClass::ParticleGroupClass, GroupCpuWrapper> GroupAlgorithm;


    FMpi mpiComm(argc, argv);
    // Get params
    const int NbLevels      = FParameters::getValue(argc,argv,FParameterDefinitions::OctreeHeight.options, 5);
    const FSize NbParticles   = FParameters::getValue(argc,argv,FParameterDefinitions::NbParticles.options, FSize(20));
    const int groupSize      = FParameters::getValue(argc,argv,LocalOptionBlocSize.options, 250);
    const FSize totalNbParticles = (NbParticles*mpiComm.global().processCount());

    // Load the particles
    FRandomLoader<FReal> loader(NbParticles, 1.0, FPoint<FReal>(0,0,0), mpiComm.global().processId());
    FAssertLF(loader.isOpen());

    // Fill the particles
    struct TestParticle{
        FPoint<FReal> position;
        const FPoint<FReal>& getPosition(){
            return position;
        }
    };

    std::unique_ptr<TestParticle[]> particles(new TestParticle[loader.getNumberOfParticles()]);
    memset(particles.get(), 0, sizeof(TestParticle) * loader.getNumberOfParticles());
    for(FSize idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
        loader.fillParticle(&particles[idxPart].position);
    }
    // Sort in parallel
    FVector<TestParticle> myParticles;
    FLeafBalance balancer;
    FMpiTreeBuilder<FReal, TestParticle >::DistributeArrayToContainer(mpiComm.global(),
                                                                particles.get(),
                                                                loader.getNumberOfParticles(),
                                                                loader.getCenterOfBox(),
                                                                loader.getBoxWidth(),
                                                                NbLevels,
                                                                &myParticles,
                                                                &balancer);

    FTestParticleContainer<FReal> allParticles;
    for(FSize idxPart = 0 ; idxPart < myParticles.getSize() ; ++idxPart){
        allParticles.push(myParticles[idxPart].position);
    }

    // Each proc need to know the righest morton index
    const FTreeCoordinate host = FCoordinateComputer::GetCoordinateFromPosition<FReal>(
                loader.getCenterOfBox(),
                loader.getBoxWidth(),
                NbLevels,
                myParticles[myParticles.getSize()-1].position );
    const MortonIndex myLeftLimite = host.getMortonIndex();
    MortonIndex leftLimite = -1;
    if(mpiComm.global().processId() != 0){
        FMpi::Assert(MPI_Recv(&leftLimite, sizeof(leftLimite), MPI_BYTE,
                              mpiComm.global().processId()-1, 0,
                              mpiComm.global().getComm(), MPI_STATUS_IGNORE), __LINE__);
    }
    if(mpiComm.global().processId() != mpiComm.global().processCount()-1){
        FMpi::Assert(MPI_Send(const_cast<MortonIndex*>(&myLeftLimite), sizeof(myLeftLimite), MPI_BYTE,
                              mpiComm.global().processId()+1, 0,
                              mpiComm.global().getComm()), __LINE__);
    }
    FLOG(std::cout << "My last index is " << leftLimite << "\n");
    FLOG(std::cout << "My left limite is " << myLeftLimite << "\n");


    // Put the data into the tree
    GroupOctreeClass groupedTree(NbLevels, loader.getBoxWidth(), loader.getCenterOfBox(), groupSize,
                                 &allParticles, true, leftLimite);
    groupedTree.printInfoBlocks();

    // Run the algorithm
    GroupKernelClass groupkernel;
    GroupAlgorithm groupalgo(mpiComm.global(), &groupedTree,&groupkernel);
    groupalgo.execute();

    std::cout << "Wait Others... " << std::endl;
    mpiComm.global().barrier();

    groupedTree.forEachCellLeaf<GroupContainerClass>(
        [&](GroupCellSymbClass*  gsymb,
            GroupCellUpClass*    /* gmul */,
            GroupCellDownClass*  /* gloc */,
            GroupContainerClass* leaf)
        {
            const FSize nbPartsInLeaf = leaf->getNbParticles();
            const long long int* dataDown = leaf->getDataDown();
            for(FSize idxPart = 0 ; idxPart < nbPartsInLeaf ; ++idxPart){
                if(dataDown[idxPart] != totalNbParticles-1){
                    std::cout << "[Full] Error a particle has " << dataDown[idxPart]
                              << " (it should be " << (totalNbParticles-1)
                              << ") at index " << gsymb->getMortonIndex()
                              << "\n";
                }
            }
        });



    typedef FTestCell                   CellClass;
    typedef FTestParticleContainer<FReal>      ContainerClass;
    typedef FSimpleLeaf<FReal, ContainerClass >                     LeafClass;
    typedef FOctree<FReal, CellClass, ContainerClass , LeafClass >  OctreeClass;
    typedef FTestKernels< CellClass, ContainerClass >         KernelClass;
    typedef FFmmAlgorithm<OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass >     FmmClass;

    {
        // Usual octree
        OctreeClass tree(NbLevels, 2, loader.getBoxWidth(), loader.getCenterOfBox());
        for(int idxProc = 0 ; idxProc < mpiComm.global().processCount() ; ++idxProc){
            FRandomLoader<FReal> loaderAll(NbParticles, 1.0, FPoint<FReal>(0,0,0), idxProc);
            for(FSize idxPart = 0 ; idxPart < loaderAll.getNumberOfParticles() ; ++idxPart){
                FPoint<FReal> pos;
                loaderAll.fillParticle(&pos);
                tree.insert(pos);
            }
        }
        // Usual algorithm
        KernelClass kernels;            // FTestKernels FBasicKernels
        FmmClass algo(&tree,&kernels);  //FFmmAlgorithm FFmmAlgorithmThread
        algo.execute();

        // Compare the results
        groupedTree.forEachCellWithLevel(
            [&](GroupCellSymbClass* gsymb,
                GroupCellUpClass*   gmul,
                GroupCellDownClass* gloc,
                const int level)
            {
                const CellClass* cell = tree.getCell(gsymb->getMortonIndex(), level);
                if(cell == nullptr){
                    std::cout << "[Empty] Error cell should not exist " << gsymb->getMortonIndex() << "\n";
                }
                else {
                    if(*gmul != cell->getDataUp()){
                        std::cout << "[Up] Up is different at index " << gsymb->getMortonIndex() << " level " << level << " is " << *gmul << " should be " << cell->getDataUp() << "\n";
                    }
                    if(*gloc != cell->getDataDown()){
                        std::cout << "[Down] Down is different at index " << gsymb->getMortonIndex() << " level " << level << " is " << *gloc << " should be " << cell->getDataDown() << "\n";
                    }
                }
            });
    }

    return 0;
}
