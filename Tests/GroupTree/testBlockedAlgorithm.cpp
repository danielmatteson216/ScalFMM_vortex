
// Keep in private GIT


#include "Utils/FGlobal.hpp"

//#undef  SCALFMM_USE_STARPU
//#undef  SCALFMM_USE_OMP4

#include "GroupTree/Core/FGroupTree.hpp"

#include "Components/FSimpleLeaf.hpp"
#include "Containers/FVector.hpp"


#include "Utils/FMath.hpp"
#include "Utils/FMemUtils.hpp"
#include "Utils/FParameters.hpp"

#include "Files/FRandomLoader.hpp"

#include "GroupTree/Core/FGroupSeqAlgorithm.hpp"
#ifdef SCALFMM_USE_OMP4
#include "GroupTree/Core/FGroupTaskDepAlgorithm.hpp"
#endif
#ifdef SCALFMM_USE_STARPU
#include "GroupTree/Core/FGroupTaskStarpuAlgorithm.hpp"
#include "GroupTree/StarPUUtils/FStarPUKernelCapacities.hpp"

#include "GroupTree/StarPUUtils/FStarPUCpuWrapper.hpp"
#endif
#include "GroupTree/Core/FP2PGroupParticleContainer.hpp"
#include "GroupTree/Core/FGroupTaskAlgorithm.hpp"

#include "Utils/FParameterNames.hpp"

#include "Components/FTestParticleContainer.hpp"
#include "Components/FTestCell.hpp"
#include "Components/FSymbolicData.hpp"
#include "Components/FTestKernels.hpp"
#include "GroupTree/TestKernel/FGroupTestParticleContainer.hpp"

#include "Files/FFmaGenericLoader.hpp"
#include "Core/FFmmAlgorithm.hpp"

int main(int argc, char* argv[]){
    setenv("STARPU_NCPU","1",1);
    const FParameterNames LocalOptionBlocSize {
        {"-bs"},
        "The size of the block of the blocked tree"
    };
    FHelpDescribeAndExit(argc, argv, "Test the blocked tree by counting the particles.",
                         FParameterDefinitions::OctreeHeight, FParameterDefinitions::NbParticles,
                         FParameterDefinitions::OctreeSubHeight, LocalOptionBlocSize);

    typedef double FReal;

    // Initialize the types
    using GroupCellClass     = FTestCell;
    using GroupCellUpClass   = typename FTestCell::multipole_t;
    using GroupCellDownClass = typename FTestCell::local_expansion_t;
    using GroupCellSymbClass = FSymbolicData;

    typedef FGroupTestParticleContainer<FReal>                                GroupContainerClass;
    typedef FGroupTree< FReal, GroupCellSymbClass, GroupCellUpClass, GroupCellDownClass,
            GroupContainerClass, 0, 1, long long int>  GroupOctreeClass;
#ifdef SCALFMM_USE_STARPU
    typedef FStarPUAllCpuCapacities<FTestKernels< GroupCellClass, GroupContainerClass >>  GroupKernelClass;
    typedef FStarPUCpuWrapper<typename GroupOctreeClass::CellGroupClass, GroupCellClass, GroupKernelClass, typename GroupOctreeClass::ParticleGroupClass, GroupContainerClass> GroupCpuWrapper;
    typedef FGroupTaskStarPUAlgorithm<GroupOctreeClass, typename GroupOctreeClass::CellGroupClass, GroupKernelClass, typename GroupOctreeClass::ParticleGroupClass, GroupCpuWrapper, GroupContainerClass > GroupAlgorithm;
#elif defined(SCALFMM_USE_OMP4)
    typedef FTestKernels< GroupCellClass, GroupContainerClass >  GroupKernelClass;
    typedef FGroupTaskDepAlgorithm<GroupOctreeClass, typename GroupOctreeClass::CellGroupClass,
            GroupCellSymbClass, GroupCellUpClass, GroupCellDownClass, GroupKernelClass, typename GroupOctreeClass::ParticleGroupClass, GroupContainerClass > GroupAlgorithm;
#else
    typedef FTestKernels< GroupCellClass, GroupContainerClass >  GroupKernelClass;
    //typedef FGroupSeqAlgorithm<GroupOctreeClass, typename GroupOctreeClass::CellGroupClass, GroupCellClass, GroupKernelClass, typename GroupOctreeClass::ParticleGroupClass, GroupContainerClass > GroupAlgorithm;
    typedef FGroupTaskAlgorithm<GroupOctreeClass, typename GroupOctreeClass::CellGroupClass, GroupKernelClass, typename GroupOctreeClass::ParticleGroupClass, GroupContainerClass > GroupAlgorithm;
#endif

    typedef FTestCell                   CellClass;
    typedef FTestParticleContainer<FReal>      ContainerClass;
    typedef FSimpleLeaf<FReal, ContainerClass >                     LeafClass;
    typedef FOctree<FReal, CellClass, ContainerClass , LeafClass >  OctreeClass;
    typedef FTestKernels< CellClass, ContainerClass >         KernelClass;

    // FFmmAlgorithmTask FFmmAlgorithmThread
    typedef FFmmAlgorithm<OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass >     FmmClass;

    // Get params
    const int NbLevels      = FParameters::getValue(argc,argv,FParameterDefinitions::OctreeHeight.options, 5);
    const int groupSize      = FParameters::getValue(argc,argv,LocalOptionBlocSize.options, 250);

//#define LOAD_FILE
#ifndef LOAD_FILE
    const FSize NbParticles   = FParameters::getValue(argc,argv,FParameterDefinitions::NbParticles.options, FSize(20));
    FRandomLoader<FReal> loader(NbParticles, 1.0, FPoint<FReal>(0,0,0), 0);
#else
    // Load the particles
    const char* const filename = FParameters::getStr(argc,argv,FParameterDefinitions::InputFile.options, "../Data/test20k.fma");
    FFmaGenericLoader<FReal> loader(filename);
#endif
    FAssertLF(loader.isOpen());

    // Usual octree
    OctreeClass tree(NbLevels, FParameters::getValue(argc,argv,FParameterDefinitions::OctreeSubHeight.options, 2),
                     loader.getBoxWidth(), loader.getCenterOfBox());

    FTestParticleContainer<FReal> allParticles;
    for(FSize idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
        FPoint<FReal> particlePosition;
#ifndef LOAD_FILE
        loader.fillParticle(&particlePosition);
#else
        FReal ph;
        loader.fillParticle(&particlePosition, &ph);
#endif
        allParticles.push(particlePosition);
        tree.insert(particlePosition);
    }

    // Put the data into the tree
    //GroupOctreeClass groupedTree(NbLevels, groupSize, &tree);
    //GroupOctreeClass groupedTree(NbLevels, loader.getBoxWidth(), loader.getCenterOfBox(), groupSize, &allParticles);
    //GroupOctreeClass groupedTree(NbLevels, loader.getBoxWidth(), loader.getCenterOfBox(), groupSize, &allParticles, false, true);
    GroupOctreeClass groupedTree(NbLevels, loader.getBoxWidth(), loader.getCenterOfBox(), groupSize, &allParticles, false, true, 0.2);
    groupedTree.printInfoBlocks();

    // Check tree structure at leaf level
    // groupedTree.forEachCellLeaf<GroupContainerClass>([&](GroupCellClass gcell, GroupContainerClass* gleaf){
    groupedTree.forEachCellLeaf<GroupContainerClass>(
        [&](GroupCellSymbClass*  gsymb,
            GroupCellUpClass*    /* gmul */,
            GroupCellDownClass*  /* gloc */,
            GroupContainerClass* gleaf)
        {
            const ContainerClass* src = tree.getLeafSrc(gsymb->getMortonIndex());
            if(src == nullptr){
                std::cout << "[PartEmpty] Error cell should not exist "
                          << gsymb->getMortonIndex() << "\n";
            } else {
                if(src->getNbParticles() != gleaf->getNbParticles()) {
                    std::cout << "[Part] Nb particles is different at index "
                              << gsymb->getMortonIndex()
                              << " is " << gleaf->getNbParticles()
                              << " should be " << src->getNbParticles() << "\n";
                }
            }
        });

    // Run the algorithm
    GroupKernelClass groupkernel;
    GroupAlgorithm groupalgo(&groupedTree,&groupkernel);
    groupalgo.execute();

    // Usual algorithm
    KernelClass kernels;            // FTestKernels FBasicKernels
    FmmClass algo(&tree,&kernels);  //FFmmAlgorithm FFmmAlgorithmThread
    algo.execute();

    // Validate the result
    groupedTree.forEachCellLeaf<GroupContainerClass>(
        [&](GroupCellSymbClass* /*gsymb*/,
            GroupCellUpClass* gmul,
            GroupCellDownClass* /*gloc*/,
            GroupContainerClass* leaf)
        {
            const FSize nbPartsInLeaf = leaf->getNbParticles();
            if(gmul->get() != nbPartsInLeaf){
                std::cout << "[P2M] Error a Cell has " << gmul->get()
                          << " (it should be " << nbPartsInLeaf << ")\n";
            }
        });
    groupedTree.forEachCellLeaf<GroupContainerClass>(
        [&](GroupCellSymbClass* gsymb,
            GroupCellUpClass* /*gmul*/,
            GroupCellDownClass* /*gloc*/,
            GroupContainerClass* leaf)
        {
            const FSize nbPartsInLeaf = leaf->getNbParticles();
            const long long int* dataDown = leaf->getDataDown();
            for(FSize idxPart = 0 ; idxPart < nbPartsInLeaf ; ++idxPart){
                if(dataDown[idxPart] != loader.getNumberOfParticles()-1){
                    std::cout << "[Full] Error a particle has " << dataDown[idxPart] << " (it should be " << (loader.getNumberOfParticles()-1) << ") at index " << gsymb->getMortonIndex() << "\n";
                }
            }
        });
    // Compare the results
    groupedTree.forEachCellWithLevel(
        [&](GroupCellSymbClass* gsymb,
            GroupCellUpClass* gmul,
            GroupCellDownClass* gloc,
            const int level)
        {
            const CellClass* cell = tree.getCell(gsymb->getMortonIndex(), level);
            if(cell == nullptr){
                std::cout << "[Empty] Error cell should not exist " << gsymb->getMortonIndex() << "\n";
            }
            else {
                if(gmul->get() != cell->getDataUp().get()){
                    std::cout << "[Up] Up is different at index " << gsymb->getMortonIndex() << " level " << level << " is " << *gsymb << " should be " << cell->getDataUp() << "\n";
                }
                if(gloc->get() != cell->getDataDown().get()){
                    std::cout << "[Down] Down is different at index " << gsymb->getMortonIndex() << " level " << level << " is " << *gloc << " should be " << cell->getDataDown() << "\n";
                }
            }
        });

    return 0;
}
