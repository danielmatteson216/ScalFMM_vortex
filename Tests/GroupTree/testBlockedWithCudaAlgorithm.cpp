
// Keep in private GIT

// @FUSE_STARPU
// @FUSE_CUDA

#include "Utils/FGlobal.hpp"

#include "GroupTree/Core/FGroupTree.hpp"

#include "Components/FSimpleLeaf.hpp"
#include "Containers/FVector.hpp"


#include "Utils/FMath.hpp"
#include "Utils/FMemUtils.hpp"
#include "Utils/FParameters.hpp"

#include "Files/FRandomLoader.hpp"

#include "GroupTree/Core/FGroupSeqAlgorithm.hpp"

#include "GroupTree/Core/FGroupTaskStarpuAlgorithm.hpp"

#include "GroupTree/Core/FP2PGroupParticleContainer.hpp"
#include "GroupTree/Core/FGroupTaskAlgorithm.hpp"

#include "Utils/FParameterNames.hpp"

#include "GroupTree/StarPUUtils/FStarPUCpuWrapper.hpp"
#include "GroupTree/StarPUUtils/FStarPUCudaWrapper.hpp"

#include "Components/FTestParticleContainer.hpp"
#include "Components/FTestCell.hpp"
#include "Components/FTestKernels.hpp"
#include "GroupTree/TestKernel/FGroupTestParticleContainer.hpp"

#include "Files/FFmaGenericLoader.hpp"
#include "Core/FFmmAlgorithm.hpp"

#include "GroupTree/StarPUUtils/FStarPUKernelCapacities.hpp"

#include "GroupTree/TestKernel/FTestCellPOD.hpp"

//#include "GroupTree/Cuda/FCudaTestKernels.hpp"
//#include "GroupTree/Cuda/FCudaGroupOfParticles.hpp"
//#include "GroupTree/Cuda/FCudaGroupAttachedLeaf.hpp"
//#include "GroupTree/Cuda/FCudaGroupOfCells.hpp"
//#include "GroupTree/Cuda/FCudaDeviceWrapper.hpp"


template <class FReal>
class FTestCudaKernels;

template <class FReal, unsigned NbSymbAttributes, unsigned NbAttributesPerParticle, class AttributeClass>
class FCudaGroupAttachedLeaf;

template <class FReal, unsigned NbSymbAttributes, unsigned NbAttributesPerParticle, class AttributeClass>
class FCudaGroupOfParticles;

template <class SymboleCellClass, class PoleCellClass, class LocalCellClass>
class FCudaGroupOfCells;

int main(int argc, char* argv[]){
    const FParameterNames LocalOptionBlocSize {
        {"-bs"},
        "The size of the block of the blocked tree"
    };
    FHelpDescribeAndExit(argc, argv, "Test the blocked tree by counting the particles.",
                         FParameterDefinitions::OctreeHeight, FParameterDefinitions::NbThreads,
                         FParameterDefinitions::NbParticles, LocalOptionBlocSize);
    // Initialize the types
    typedef double FReal;

    typedef FTestCellPODCore  GroupCellSymbClass;
    typedef FTestCellPODData  GroupCellUpClass;
    typedef FTestCellPODData  GroupCellDownClass;
    typedef FTestCellPOD      GroupCellClass;

    typedef FGroupTestParticleContainer<FReal>                                     GroupContainerClass;
    typedef FGroupTree< FReal, GroupCellClass, GroupCellSymbClass, GroupCellUpClass, GroupCellDownClass,
            GroupContainerClass, 0, 1, long long int>  GroupOctreeClass;
    typedef FStarPUAllCudaCapacities<FTestKernels< GroupCellClass, GroupContainerClass >>  GroupKernelClass;

    typedef FStarPUCpuWrapper<typename GroupOctreeClass::CellGroupClass, GroupCellClass, GroupKernelClass, typename GroupOctreeClass::ParticleGroupClass, GroupContainerClass> GroupCpuWrapper;
    typedef FStarPUCudaWrapper<GroupKernelClass, GroupCellSymbClass, GroupCellUpClass, GroupCellDownClass,
            FCudaGroupOfCells<GroupCellSymbClass, GroupCellUpClass, GroupCellDownClass>,
            FCudaGroupOfParticles<FReal, 0, 1, long long int>, FCudaGroupAttachedLeaf<FReal, 0, 1, long long int>, FTestCudaKernels<FReal> > GroupCudaWrapper;

    typedef FGroupTaskStarPUAlgorithm<GroupOctreeClass, typename GroupOctreeClass::CellGroupClass, GroupKernelClass, typename GroupOctreeClass::ParticleGroupClass,
            GroupCpuWrapper, GroupCudaWrapper> GroupAlgorithm;

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
    OctreeClass tree(NbLevels, 2, loader.getBoxWidth(), loader.getCenterOfBox());

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
    GroupOctreeClass groupedTree(NbLevels, loader.getBoxWidth(), loader.getCenterOfBox(), groupSize, &allParticles);
    groupedTree.printInfoBlocks();

    // Check tree structure at leaf level
    groupedTree.forEachCellLeaf<GroupContainerClass>([&](GroupCellClass gcell, GroupContainerClass* gleaf){
        const ContainerClass* src = tree.getLeafSrc(gcell.getMortonIndex());
        if(src == nullptr){
            std::cout << "[PartEmpty] Error cell should not exist " << gcell.getMortonIndex() << "\n";
        }
        else {
            if(src->getNbParticles() != gleaf->getNbParticles()){
                std::cout << "[Part] Nb particles is different at index " << gcell.getMortonIndex() << " is " << gleaf->getNbParticles() << " should be " << src->getNbParticles() << "\n";
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
    groupedTree.forEachCellLeaf<GroupContainerClass>([&](GroupCellClass cell, GroupContainerClass* leaf){
        const FSize nbPartsInLeaf = leaf->getNbParticles();
        if(cell.getDataUp() != nbPartsInLeaf){
            std::cout << "[P2M] Error a Cell has " << cell.getDataUp() << " (it should be " << nbPartsInLeaf << ")\n";
        }
    });
    groupedTree.forEachCellLeaf<GroupContainerClass>([&](GroupCellClass cell, GroupContainerClass* leaf){
        const FSize nbPartsInLeaf = leaf->getNbParticles();
        const long long int* dataDown = leaf->getDataDown();
        for(FSize idxPart = 0 ; idxPart < nbPartsInLeaf ; ++idxPart){
            if(dataDown[idxPart] != loader.getNumberOfParticles()-1){
                std::cout << "[Full] Error a particle has " << dataDown[idxPart] << " (it should be " << (loader.getNumberOfParticles()-1) << ") at index " << cell.getMortonIndex() << "\n";
            }
        }
    });
    // Compare the results
    groupedTree.forEachCellWithLevel([&](GroupCellClass gcell, const int level){
        const CellClass* cell = tree.getCell(gcell.getMortonIndex(), level);
        if(cell == nullptr){
            std::cout << "[Empty] Error cell should not exist " << gcell.getMortonIndex() << "\n";
        }
        else {
            if(gcell.getDataUp() != cell->getDataUp()){
                std::cout << "[Up] Up is different at index " << gcell.getMortonIndex() << " level " << level << " is " << gcell.getDataUp() << " should be " << cell->getDataUp() << "\n";
            }
            if(gcell.getDataDown() != cell->getDataDown()){
                std::cout << "[Down] Down is different at index " << gcell.getMortonIndex() << " level " << level << " is " << gcell.getDataDown() << " should be " << cell->getDataDown() << "\n";
            }
        }
    });

    return 0;
}
