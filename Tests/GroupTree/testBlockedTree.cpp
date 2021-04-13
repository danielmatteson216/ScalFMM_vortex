
// Keep in private GIT

#include "Utils/FGlobal.hpp"

#include "GroupTree/Core/FGroupTree.hpp"

#include "Components/FSimpleLeaf.hpp"
#include "Components/FSymbolicData.hpp"
#include "Containers/FVector.hpp"

#include "Containers/FOctree.hpp"

#include "Core/FFmmAlgorithm.hpp"

#include "Kernels/P2P/FP2PParticleContainer.hpp"

#include "Kernels/Rotation/FRotationKernel.hpp"
#include "Kernels/Rotation/FRotationCell.hpp"

#include "Utils/FMath.hpp"
#include "Utils/FMemUtils.hpp"
#include "Utils/FParameters.hpp"

#include "Core/FFmmAlgorithm.hpp"
#include "Core/FFmmAlgorithmThread.hpp"
#include "Core/FFmmAlgorithmTask.hpp"

#include "Files/FFmaGenericLoader.hpp"

#include "GroupTree/Core/FGroupSeqAlgorithm.hpp"
#include "GroupTree/Core/FGroupTaskAlgorithm.hpp"
#ifdef SCALFMM_USE_OMP4
#include "GroupTree/Core/FGroupTaskDepAlgorithm.hpp"
#endif
#ifdef SCALFMM_USE_STARPU
#include "GroupTree/Core/FGroupTaskStarpuAlgorithm.hpp"
#include "GroupTree/StarPUUtils/FStarPUKernelCapacities.hpp"
#endif
#include "GroupTree/Core/FP2PGroupParticleContainer.hpp"

#include "Utils/FParameterNames.hpp"


int main(int argc, char* argv[]){
    const FParameterNames LocalOptionBlocSize {
        {"-bs"},
        "The size of the block of the blocked tree"
    };
    FHelpDescribeAndExit(argc, argv,
                         "Test the blocked tree.",
                         FParameterDefinitions::OctreeHeight, FParameterDefinitions::OctreeSubHeight,
                         FParameterDefinitions::InputFile, LocalOptionBlocSize);

    typedef double FReal;
    static const int P = 3;
    typedef FRotationCell<FReal,P>               CellClass;
    typedef FP2PParticleContainer<FReal>          ContainerClass;
    typedef FSimpleLeaf<FReal, ContainerClass >                     LeafClass;
    typedef FOctree<FReal, CellClass, ContainerClass , LeafClass >  OctreeClass;

    using GroupCellClass     = FRotationCell<FReal, P>;
    using GroupCellUpClass   = typename GroupCellClass::multipole_t;
    using GroupCellDownClass = typename GroupCellClass::local_expansion_t;
    using GroupCellSymbClass = FSymbolicData;


    typedef FP2PGroupParticleContainer<FReal>          GroupContainerClass;
    typedef FGroupTree< FReal, GroupCellSymbClass, GroupCellUpClass, GroupCellDownClass, GroupContainerClass, 1, 4, FReal>  GroupOctreeClass;


    FTic counter;
    const int NbLevels      = FParameters::getValue(argc,argv,FParameterDefinitions::OctreeHeight.options, 5);
    const int SizeSubLevels = FParameters::getValue(argc,argv,FParameterDefinitions::OctreeSubHeight.options, 3);
    const char* const filename = FParameters::getStr(argc,argv,FParameterDefinitions::InputFile.options, "../Data/test20k.bin.fma.double");

    FFmaGenericLoader<FReal> loader(filename);
    FAssertLF(loader.isOpen());

    OctreeClass tree(NbLevels, SizeSubLevels, loader.getBoxWidth(), loader.getCenterOfBox());

    FP2PParticleContainer<FReal> allParticles;

    for(FSize idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
        FPoint<FReal> particlePosition;
        FReal physicalValue;
        loader.fillParticle(&particlePosition,&physicalValue);
        tree.insert(particlePosition, physicalValue );
        allParticles.push(particlePosition, physicalValue);
    }

    std::cout << "Done  " << "(@Creating and Inserting Particles = " << counter.tacAndElapsed() << "s)." << std::endl;

    const int groupSize      = FParameters::getValue(argc,argv,LocalOptionBlocSize.options, 250);

    counter.tic();
    GroupOctreeClass groupedTree2(NbLevels, groupSize, &tree);
    std::cout << "Done  " << "(@Converting the tree with all Octree = " << counter.tacAndElapsed() << "s). Group size is " << groupSize << "." << std::endl;

    counter.tic();
    GroupOctreeClass groupedTree3(NbLevels, loader.getBoxWidth(), loader.getCenterOfBox(), groupSize, &allParticles);
    std::cout << "Done  " << "(@Converting the tree with all Octree = " << counter.tacAndElapsed() << "s). Group size is " << groupSize << "." << std::endl;

    groupedTree2.printInfoBlocks();
    groupedTree3.printInfoBlocks();


#ifdef SCALFMM_USE_STARPU
    typedef FStarPUAllCpuCapacities<FRotationKernel< FReal, GroupCellClass, GroupContainerClass , P>>   GroupKernelClass;
    typedef FStarPUCpuWrapper<typename GroupOctreeClass::CellGroupClass, GroupCellClass, GroupKernelClass, typename GroupOctreeClass::ParticleGroupClass, GroupContainerClass> GroupCpuWrapper;
    typedef FGroupTaskStarPUAlgorithm<GroupOctreeClass, typename GroupOctreeClass::CellGroupClass, GroupKernelClass, typename GroupOctreeClass::ParticleGroupClass, GroupCpuWrapper, GroupContainerClass > GroupAlgorithm;
#elif defined(SCALFMM_USE_OMP4)
    typedef FRotationKernel< FReal, GroupCellClass, GroupContainerClass , P>  GroupKernelClass;
    typedef FGroupTaskDepAlgorithm<GroupOctreeClass, typename GroupOctreeClass::CellGroupClass,
            GroupCellSymbClass, GroupCellUpClass, GroupCellDownClass, GroupKernelClass, typename GroupOctreeClass::ParticleGroupClass, GroupContainerClass > GroupAlgorithm;
#else
    typedef FRotationKernel< FReal, GroupCellClass, GroupContainerClass , P>  GroupKernelClass;
    //typedef FGroupSeqAlgorithm<GroupOctreeClass, typename GroupOctreeClass::CellGroupClass, GroupCellClass, GroupKernelClass, typename GroupOctreeClass::ParticleGroupClass, GroupContainerClass > GroupAlgorithm;
    typedef FGroupTaskAlgorithm<GroupOctreeClass, typename GroupOctreeClass::CellGroupClass, GroupKernelClass, typename GroupOctreeClass::ParticleGroupClass, GroupContainerClass > GroupAlgorithm;
#endif

    GroupKernelClass kernel(NbLevels, loader.getBoxWidth(), loader.getCenterOfBox());
    GroupAlgorithm algo(&groupedTree2,&kernel);

    algo.execute();

    return 0;
}
