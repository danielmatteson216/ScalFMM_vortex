
// Keep in private GIT


#include "Utils/FGlobal.hpp"

#undef SCALFMM_USE_STARPU
#undef SCALFMM_USE_OMP4

#include "GroupTree/Core/FGroupTreeDyn.hpp"

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
#include "GroupTree/Core/FGroupTaskAlgorithm.hpp"

#include "Utils/FParameterNames.hpp"

#include "Components/FTestParticleContainer.hpp"
#include "Components/FTestCell.hpp"
#include "Components/FTestKernels.hpp"
#include "Components/FSymbolicData.hpp"

#include "Files/FFmaGenericLoader.hpp"
#include "Core/FFmmAlgorithm.hpp"


template <class FReal>
class FTestAttachedLeafDyn : public FGroupAttachedLeafDyn<FReal> {
    typedef FGroupAttachedLeafDyn<FReal> Parent;
public:
    using Parent::Parent;

    void init(const MortonIndex inIndex, const UnknownDescriptor<FReal> /*inParticles*/[],
                  const FSize inNbParticles, const size_t inSymbSize, const size_t inDownSize) override {
        memset(Parent::symbPart, 0, inSymbSize);
        memset(Parent::downPart, 0, inDownSize);
        *(FSize*)Parent::symbPart = inNbParticles;
        *(MortonIndex*)(Parent::symbPart+sizeof(FSize))= inIndex;
    }

    static void GetSizeFunc(const MortonIndex /*inIndex*/, const UnknownDescriptor<FReal> /*inParticles*/[],
                            const FSize inNbParticles, size_t* inSymbSize, size_t* inDownSize){
        *inSymbSize = sizeof(FSize) + sizeof(MortonIndex);
        *inDownSize = (inNbParticles*sizeof(long long int));
    }

    template <class ParticleClassContainer>
    static void GetSizeContainerFunc(const MortonIndex /*inIndex*/, const void* container,
                            size_t* inSymbSize, size_t* inDownSize){
        *inSymbSize = sizeof(FSize) + sizeof(MortonIndex);
        *inDownSize = (((const ParticleClassContainer*)container)->getNbParticles()*sizeof(long long int));
    }

    template<class ParticleClassContainer>
    void copyFromContainer(const MortonIndex inMindex, const ParticleClassContainer* particles){
        FAssertLF(Parent::isAttachedToSomething());
        *(FSize*)Parent::symbPart = particles->getNbParticles();
        *(MortonIndex*)(Parent::symbPart+sizeof(FSize))= inMindex;
        FAssertLF(getNbParticles() == particles->getNbParticles());
        memcpy(Parent::downPart, particles->getDataDown(), particles->getNbParticles()*sizeof(long long int));
    }

    FSize getNbParticles() const{
        return *(FSize*)Parent::symbPart;
    }

    MortonIndex getMortonIndex() const{
        return *(MortonIndex*)(Parent::symbPart+sizeof(FSize));
    }

    long long int* getDataDown(){
        return (long long int*) Parent::downPart;
    }

    const long long int* getDataDown()const{
        return (const long long int*) Parent::downPart;
    }
};


class FTestCellDyn :public FTestCell {
public:

    using FTestCell::FTestCell;

    void release(){
        // nothing
    }
};



int main(int argc, char* argv[]){
    setenv("STARPU_NCPU","1",1);
    const FParameterNames LocalOptionBlocSize {
        {"-bs"},
        "The size of the block of the blocked tree"
    };
    FHelpDescribeAndExit(argc, argv, "Test the blocked tree by counting the particles.",
                         FParameterDefinitions::OctreeHeight, FParameterDefinitions::NbParticles, LocalOptionBlocSize);

    typedef double FReal;

    using GroupCellClass     = FTestCell;
    using GroupCellUpClass   = typename FTestCell::multipole_t;
    using GroupCellDownClass = typename FTestCell::local_expansion_t;
    using GroupCellSymbClass = FSymbolicData;

    typedef FTestAttachedLeafDyn<FReal> GroupContainerClass;

    typedef FGroupTreeDyn<FReal, GroupCellSymbClass, GroupCellUpClass, GroupCellDownClass, GroupContainerClass>  GroupOctreeClass;
#ifdef SCALFMM_USE_STARPU
    typedef FStarPUAllCpuCapacities<FTestKernels< GroupCellClass, GroupContainerClass >>  GroupKernelClass;
    typedef FStarPUCpuWrapper<typename GroupOctreeClass::CellGroupClass, GroupCellClass, GroupKernelClass, typename GroupOctreeClass::ParticleGroupClass, GroupContainerClass> GroupCpuWrapper;
    typedef FGroupTaskStarPUAlgorithm<GroupOctreeClass, typename GroupOctreeClass::CellGroupClass, GroupKernelClass, typename GroupOctreeClass::ParticleGroupClass, GroupCpuWrapper, GroupContainerClass > GroupAlgorithm;
#elif defined(SCALFMM_USE_OMP4)
    typedef FTestKernels< GroupCellClass, GroupContainerClass >  GroupKernelClass;
    typedef FGroupTaskDepAlgorithm<GroupOctreeClass, typename GroupOctreeClass::CellGroupClass,
            unsigned char, unsigned char, unsigned char, GroupKernelClass, typename GroupOctreeClass::ParticleGroupClass, GroupContainerClass > GroupAlgorithm;
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
    OctreeClass tree(NbLevels, 2, loader.getBoxWidth(), loader.getCenterOfBox());

    std::unique_ptr<UnknownDescriptor<FReal>[]> allParticles(new UnknownDescriptor<FReal>[loader.getNumberOfParticles()]);
    for(FSize idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
        FPoint<FReal> particlePosition;
#ifndef LOAD_FILE
        loader.fillParticle(&particlePosition);
#else
        FReal ph;
        loader.fillParticle(&particlePosition, &ph);
#endif
        allParticles[idxPart].pos = (particlePosition);
        tree.insert(particlePosition);
    }

    std::unique_ptr<size_t[]> cellSymbSizePerLevel(new size_t[NbLevels]);
    std::unique_ptr<size_t[]> cellUpSizePerLevel(new size_t[NbLevels]);
    std::unique_ptr<size_t[]> cellDownSizePerLevel(new size_t[NbLevels]);
     for(int idx = 0 ; idx < NbLevels ; ++idx){
         cellSymbSizePerLevel[idx] = sizeof(GroupCellSymbClass);
         cellUpSizePerLevel[idx] = sizeof(GroupCellUpClass);
         cellDownSizePerLevel[idx] = sizeof(GroupCellDownClass);
     }

    // Put the data into the tree
//    GroupOctreeClass groupedTree(NbLevels, groupSize, &tree,
//              cellSymbSizePerLevel.get(), cellUpSizePerLevel.get(), cellDownSizePerLevel.get(),
//                                 [](const MortonIndex inIndex, const void* inParticles,
//                                    size_t* inSymbSize, size_t* inDownSize) {
//                                        GroupContainerClass::GetSizeContainerFunc<ContainerClass>(
//                                                    inIndex, inParticles, inSymbSize, inDownSize);
//                                 },
//     [](const MortonIndex /*mindex*/,
//                        unsigned char* symbBuff, const size_t /*symbSize*/,
//                        unsigned char* upBuff, const size_t /*upSize*/,
//                        unsigned char* downBuff, const size_t /*downSize*/,
     //                         const int /*inLevel*/){
//         GroupCellClass cell(symbBuff, upBuff, downBuff);
//     });

     GroupOctreeClass groupedTree(
         NbLevels, loader.getBoxWidth(), loader.getCenterOfBox(), groupSize,
         cellSymbSizePerLevel.get(), cellUpSizePerLevel.get(), cellDownSizePerLevel.get(),
         allParticles.get(), loader.getNumberOfParticles(),
         [](const MortonIndex inIndex, const UnknownDescriptor<FReal> inParticles[],
            const FSize inNbParticles, size_t* inSymbSize, size_t* inDownSize){
             GroupContainerClass::GetSizeFunc(inIndex, inParticles, inNbParticles, inSymbSize, inDownSize);
         },
         [](const MortonIndex inIndex, const UnknownDescriptor<FReal> inParticles[],
            const FSize inNbParticles, unsigned char* symbBuffer, const size_t inSymbSize,
            unsigned char* downBuffer, const size_t inDownSize){
             GroupContainerClass leaf(symbBuffer, downBuffer);
             leaf.init(inIndex, inParticles, inNbParticles, inSymbSize, inDownSize);
         },
         [](const MortonIndex /*mindex*/,
            unsigned char* /*symbBuff*/, const size_t /*symbSize*/,
            unsigned char* /*upBuff*/, const size_t /*upSize*/,
            unsigned char* /*downBuff*/, const size_t /*downSize*/,
            const int /*inLevel*/){
             //GroupCellClass cell(symbBuff, upBuff, downBuff);
         });

//    GroupOctreeClass groupedTree(NbLevels, loader.getBoxWidth(), loader.getCenterOfBox(), groupSize,
//                                 cellSymbSizePerLevel.get(), cellUpSizePerLevel.get(), cellDownSizePerLevel.get(),
//                                  allParticles.get(), loader.getNumberOfParticles(),
//                                  [](const MortonIndex inIndex, const UnknownDescriptor<FReal> inParticles[],
//                                     const FSize inNbParticles, size_t* inSymbSize, size_t* inDownSize){
//                                         GroupContainerClass::GetSizeFunc(inIndex, inParticles, inNbParticles, inSymbSize, inDownSize);
//                                  },
//    [](const MortonIndex inIndex, const UnknownDescriptor<FReal> inParticles[],
//                  const FSize inNbParticles, unsigned char* symbBuffer, const size_t inSymbSize,
//    unsigned char* downBuffer, const size_t inDownSize){
//        GroupContainerClass leaf(symbBuffer, downBuffer);
//        leaf.init(inIndex, inParticles, inNbParticles, inSymbSize, inDownSize);
//    },
//    [](const MortonIndex /*mindex*/,
//                       unsigned char* symbBuff, const size_t /*symbSize*/,
//                       unsigned char* upBuff, const size_t /*upSize*/,
//                       unsigned char* downBuff, const size_t /*downSize*/,
     //                         const int /*inLevel*/){
//        GroupCellClass cell(symbBuff, upBuff, downBuff);
//    }
//                                 false, true);
//    GroupOctreeClass groupedTree(NbLevels, loader.getBoxWidth(), loader.getCenterOfBox(), groupSize,
//                                cellSymbSizePerLevel.get(), cellUpSizePerLevel.get(), cellDownSizePerLevel.get(),
//                                 allParticles.get(), loader.getNumberOfParticles(),
//                                 [](const MortonIndex inIndex, const UnknownDescriptor<FReal> inParticles[],
//                                    const FSize inNbParticles, size_t* inSymbSize, size_t* inDownSize){
//                                        GroupContainerClass::GetSizeFunc(inIndex, inParticles, inNbParticles, inSymbSize, inDownSize);
//                                 },
//    [](const MortonIndex inIndex, const UnknownDescriptor<FReal> inParticles[],
//                  const FSize inNbParticles, unsigned char* symbBuffer, const size_t inSymbSize,
//    unsigned char* downBuffer, const size_t inDownSize){
//        GroupContainerClass leaf(symbBuffer, downBuffer);
//        leaf.init(inIndex, inParticles, inNbParticles, inSymbSize, inDownSize);
//    },
//    [](const MortonIndex /*mindex*/,
//                       unsigned char* symbBuff, const size_t /*symbSize*/,
//                       unsigned char* upBuff, const size_t /*upSize*/,
//                       unsigned char* downBuff, const size_t /*downSize*/,
//                         const int /*inLevel*/){
//        GroupCellClass cell(symbBuff, upBuff, downBuff);
//    }
//                                false, true, 0.2);
     groupedTree.printInfoBlocks();

    // Check tree structure at leaf level
     groupedTree.forEachCellLeaf<GroupContainerClass>(
         [&](GroupCellSymbClass* gsymb,
             GroupCellUpClass* /*gmul*/,
             GroupCellDownClass* /*gloc*/,
             GroupContainerClass* gleaf){
             const ContainerClass* src = tree.getLeafSrc(gsymb->getMortonIndex());
             if(src == nullptr){
                 std::cout << "[PartEmpty] Error cell should not exist " << gsymb->getMortonIndex() << "\n";
             }
             else {
                 if(src->getNbParticles() != gleaf->getNbParticles()){
                     std::cout << "[Part] Nb particles is different at index "
                               << gsymb->getMortonIndex()
                               << " is " << gleaf->getNbParticles()
                               << " should be " << src->getNbParticles()
                               << "\n";
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
                          << " (it should be " << nbPartsInLeaf
                          << ")\n";
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
                    std::cout << "[Full] Error a particle has " << dataDown[idxPart]
                              << " (it should be " << (loader.getNumberOfParticles()-1)
                              << ") at index " << gsymb->getMortonIndex()
                              << "\n";
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
                if(*gmul != cell->getDataUp()){
                    std::cout << "[Up] Up is different at index " << gsymb->getMortonIndex() << " level " << level << " is " << *gmul << " should be " << cell->getDataUp() << "\n";
                }
                if(*gloc != cell->getDataDown()){
                    std::cout << "[Down] Down is different at index " << gsymb->getMortonIndex() << " level " << level << " is " << *gloc << " should be " << cell->getDataDown() << "\n";
                }
            }
        });

    return 0;
}
