// ==== CMAKE =====
// @FUSE_BLAS
// @FUSE_FFT
// @FUSE_STARPU
// @FUSE_CUDA
// ================
// Keep in private GIT


#include "Utils/FGlobal.hpp"

#include "GroupTree/Core/FGroupTree.hpp"

#include "Components/FSimpleLeaf.hpp"
#include "Containers/FVector.hpp"

#include "Kernels/P2P/FP2PParticleContainer.hpp"

#include "Kernels/Interpolation/FInterpMatrixKernel.hpp"
#include "Kernels/Uniform/FUnifKernel.hpp"

#include "GroupTree/Uniform/FUnifCellPOD.hpp"
#include "GroupTree/Uniform/FUnifCudaCellPOD.hpp"

#include "Utils/FMath.hpp"
#include "Utils/FMemUtils.hpp"
#include "Utils/FParameters.hpp"

#include "Files/FRandomLoader.hpp"
#include "Files/FFmaGenericLoader.hpp"

#include "GroupTree/Core/FGroupSeqAlgorithm.hpp"
#include "GroupTree/Core/FGroupTaskAlgorithm.hpp"

#include "GroupTree/Core/FGroupTaskStarpuAlgorithm.hpp"
#include "GroupTree/StarPUUtils/FStarPUKernelCapacities.hpp"

#include "GroupTree/Core/FP2PGroupParticleContainer.hpp"

#include "GroupTree/Cuda/FCudaDeviceWrapper.hpp"
#include "GroupTree/Cuda/FCudaEmptyCellSymb.hpp"
#include "GroupTree/Cuda/FCudaGroupOfParticles.hpp"
#include "GroupTree/Cuda/FCudaGroupOfCells.hpp"

#include "GroupTree/Uniform/FUnifCudaSharedData.hpp"

#include "Utils/FParameterNames.hpp"

#include <memory>

template <class FReal, int ORDER>
class FUnifCuda;

#define RANDOM_PARTICLES

template <class BaseClass>
class FStarPUP2PM2LCudaOnlyCapacities : public BaseClass, public FStarPUAbstractCapacities {
    bool check(const FStarPUTypes inPu) const override{
        return inPu == FSTARPU_CPU_IDX;
    }
public:
    using BaseClass::BaseClass;

    bool supportP2P(const FStarPUTypes inPu) const override {
        return inPu == FSTARPU_CUDA_IDX;
    }
    bool supportP2PExtern(const FStarPUTypes inPu) const override {
        return inPu == FSTARPU_CUDA_IDX;
    }
    bool supportP2PMpi(const FStarPUTypes inPu) const override {
        return inPu == FSTARPU_CUDA_IDX;
    }

    bool supportM2L(const FStarPUTypes inPu) const override {
        return inPu == FSTARPU_CUDA_IDX;
    }
    bool supportM2LExtern(const FStarPUTypes inPu) const override {
        return inPu == FSTARPU_CUDA_IDX;
    }
    bool supportM2LMpi(const FStarPUTypes inPu) const override {
        return inPu == FSTARPU_CUDA_IDX;
    }
};


const FParameterNames LocalOptionBlocSize { {"-bs"}, "The size of the block of the blocked tree"};
const FParameterNames LocalOptionValidate { {"-validation"}, "To compare with direct computation"};

template <template <class> class PUCapacities>
int mainCore(int argc, char* argv[]){
    // Initialize the types
    typedef double FReal;
    static const int ORDER = 5;
    typedef FInterpMatrixKernelR<FReal> MatrixKernelClass;

    typedef FUnifCellPODCore         GroupCellSymbClass;
    typedef FUnifCellPODPole<FReal,ORDER>  GroupCellUpClass;
    typedef FUnifCellPODLocal<FReal,ORDER> GroupCellDownClass;
    typedef FUnifCellPOD<FReal,ORDER>      GroupCellClass;

    typedef FP2PGroupParticleContainer<FReal>          GroupContainerClass;
    typedef FGroupTree< FReal, GroupCellSymbClass, GroupCellUpClass, GroupCellDownClass, GroupContainerClass, 1, 4, FReal>  GroupOctreeClass;

    typedef PUCapacities<FUnifKernel<FReal,GroupCellClass,GroupContainerClass,MatrixKernelClass,ORDER>> GroupKernelClass;
    typedef FStarPUCpuWrapper<typename GroupOctreeClass::CellGroupClass, GroupCellClass, GroupKernelClass, typename GroupOctreeClass::ParticleGroupClass, GroupContainerClass> GroupCpuWrapper;

    typedef FStarPUCudaWrapper<GroupKernelClass,
            FBasicCellPOD, FCudaUnifCellPODPole<FReal,ORDER>,FCudaUnifCellPODLocal<FReal,ORDER>,
            FCudaGroupOfCells<FBasicCellPOD, FCudaUnifCellPODPole<FReal,ORDER>,FCudaUnifCellPODLocal<FReal,ORDER> >,
            FCudaGroupOfParticles<FReal, 1, 4, FReal>, FCudaGroupAttachedLeaf<FReal, 1, 4, FReal>, FUnifCuda<FReal,ORDER> > GroupCudaWrapper;

    typedef FGroupTaskStarPUAlgorithm<GroupOctreeClass, typename GroupOctreeClass::CellGroupClass, GroupKernelClass, typename GroupOctreeClass::ParticleGroupClass,
            GroupCpuWrapper, GroupCudaWrapper > GroupAlgorithm;

    // Get params
    const int NbLevels      = FParameters::getValue(argc,argv,FParameterDefinitions::OctreeHeight.options, 5);
    const int groupSize     = FParameters::getValue(argc,argv,LocalOptionBlocSize.options, 250);

    // Load the particles
#ifdef RANDOM_PARTICLES
    FRandomLoader<FReal> loader(FParameters::getValue(argc,argv,FParameterDefinitions::NbParticles.options, 2000), 1.0, FPoint<FReal>(0,0,0), 0);
#else
    const char* const filename = FParameters::getStr(argc,argv,FParameterDefinitions::InputFile.options, "../Data/test20k.fma");
    FFmaGenericLoader<FReal> loader(filename);
#endif
    FAssertLF(loader.isOpen());
    FTic timer;

    FP2PParticleContainer<FReal> allParticles;
    for(FSize idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
        FPoint<FReal> particlePosition;
        FReal physicalValue;
#ifdef RANDOM_PARTICLES
        physicalValue = 0.10;
        loader.fillParticle(&particlePosition);
#else
        loader.fillParticle(&particlePosition, &physicalValue);
#endif
        allParticles.push(particlePosition, physicalValue);
    }
    std::cout << "Particles loaded in " << timer.tacAndElapsed() << "s\n";

    // Put the data into the tree
    timer.tic();
    GroupOctreeClass groupedTree(NbLevels, loader.getBoxWidth(), loader.getCenterOfBox(), groupSize, &allParticles);
    groupedTree.printInfoBlocks();
    std::cout << "Tree created in " << timer.tacAndElapsed() << "s\n";

    // Run the algorithm
    const MatrixKernelClass MatrixKernel;
    GroupKernelClass groupkernel(NbLevels, loader.getBoxWidth(), loader.getCenterOfBox(), &MatrixKernel);
    GroupAlgorithm groupalgo(&groupedTree,&groupkernel);

    {
        typedef FUnifM2LHandler<FReal, ORDER,MatrixKernelClass::Type> M2LHandlerClass;
        const M2LHandlerClass M2LHandler(&MatrixKernel,
                                         NbLevels,
                                         loader.getBoxWidth(),
                                         1);
        // Copy precomputed matrix to cuda
        FUnifCudaSharedData<FReal, ORDER> hostData;
        hostData.BoxWidth = loader.getBoxWidth();

        std::cout << "Copy precomputed tab of size " << hostData.ninteractions << " " << hostData.opt_rc << std::endl;

        for(int idxInter = 0 ; idxInter < hostData.ninteractions ; ++idxInter){
            for(int idx = 0 ; idx < hostData.opt_rc ; ++idx){
                hostData.FC[hostData.opt_rc*idxInter + idx].complex[0] = M2LHandler.getFc(idxInter, idx).getReal();
                hostData.FC[hostData.opt_rc*idxInter + idx].complex[1] = M2LHandler.getFc(idxInter, idx).getImag();
            }
        }
        std::cout << "Copy to GPU" << std::endl;

        groupalgo.forEachCudaWorker([&](void* cudaKernel){
            FUnifCudaFillObject(cudaKernel,hostData);
        });

        std::cout << "Done" << std::endl;
    }

    timer.tic();
    groupalgo.execute();
    std::cout << "@EXEC TIME = " << timer.tacAndElapsed() << "s\n";

    // Validate the result
    if(FParameters::existParameter(argc, argv, LocalOptionValidate.options)){
        FSize offsetParticles = 0;
        FReal*const allPhysicalValues = allParticles.getPhysicalValues();
        FReal*const allPosX = const_cast<FReal*>( allParticles.getPositions()[0]);
        FReal*const allPosY = const_cast<FReal*>( allParticles.getPositions()[1]);
        FReal*const allPosZ = const_cast<FReal*>( allParticles.getPositions()[2]);

        groupedTree.forEachCellLeaf<FP2PGroupParticleContainer<FReal> >([&](GroupCellClass cellTarget, FP2PGroupParticleContainer<FReal> * leafTarget){
            const FReal*const physicalValues = leafTarget->getPhysicalValues();
            const FReal*const posX = leafTarget->getPositions()[0];
            const FReal*const posY = leafTarget->getPositions()[1];
            const FReal*const posZ = leafTarget->getPositions()[2];
            const FSize nbPartsInLeafTarget = leafTarget->getNbParticles();

            for(FSize idxPart = 0 ; idxPart < nbPartsInLeafTarget ; ++idxPart){
                allPhysicalValues[offsetParticles + idxPart] = physicalValues[idxPart];
                allPosX[offsetParticles + idxPart] = posX[idxPart];
                allPosY[offsetParticles + idxPart] = posY[idxPart];
                allPosZ[offsetParticles + idxPart] = posZ[idxPart];
            }

            offsetParticles += nbPartsInLeafTarget;
        });

        FAssertLF(offsetParticles == loader.getNumberOfParticles());

        FReal*const allDirectPotentials = allParticles.getPotentials();
        FReal*const allDirectforcesX = allParticles.getForcesX();
        FReal*const allDirectforcesY = allParticles.getForcesY();
        FReal*const allDirectforcesZ = allParticles.getForcesZ();

        for(int idxTgt = 0 ; idxTgt < offsetParticles ; ++idxTgt){
            for(int idxMutual = idxTgt + 1 ; idxMutual < offsetParticles ; ++idxMutual){
                FP2PR::MutualParticles(
                    allPosX[idxTgt],allPosY[idxTgt],allPosZ[idxTgt], allPhysicalValues[idxTgt],
                    &allDirectforcesX[idxTgt], &allDirectforcesY[idxTgt], &allDirectforcesZ[idxTgt], &allDirectPotentials[idxTgt],
                    allPosX[idxMutual],allPosY[idxMutual],allPosZ[idxMutual], allPhysicalValues[idxMutual],
                    &allDirectforcesX[idxMutual], &allDirectforcesY[idxMutual], &allDirectforcesZ[idxMutual], &allDirectPotentials[idxMutual]
                );
            }
        }

        FMath::FAccurater<FReal> potentialDiff;
        FMath::FAccurater<FReal> fx, fy, fz;
        offsetParticles = 0;
        groupedTree.forEachCellLeaf<FP2PGroupParticleContainer<FReal> >([&](GroupCellClass cellTarget, FP2PGroupParticleContainer<FReal> * leafTarget){
            const FReal*const potentials = leafTarget->getPotentials();
            const FReal*const forcesX = leafTarget->getForcesX();
            const FReal*const forcesY = leafTarget->getForcesY();
            const FReal*const forcesZ = leafTarget->getForcesZ();
            const FSize nbPartsInLeafTarget = leafTarget->getNbParticles();

            for(int idxTgt = 0 ; idxTgt < nbPartsInLeafTarget ; ++idxTgt){
                potentialDiff.add(allDirectPotentials[idxTgt + offsetParticles], potentials[idxTgt]);
                fx.add(allDirectforcesX[idxTgt + offsetParticles], forcesX[idxTgt]);
                fy.add(allDirectforcesY[idxTgt + offsetParticles], forcesY[idxTgt]);
                fz.add(allDirectforcesZ[idxTgt + offsetParticles], forcesZ[idxTgt]);
            }

            offsetParticles += nbPartsInLeafTarget;
        });

        std::cout << "Error : Potential " << potentialDiff << "\n";
        std::cout << "Error : fx " << fx << "\n";
        std::cout << "Error : fy " << fy << "\n";
        std::cout << "Error : fz " << fz << "\n";
    }

    return 0;
}

int main(int argc, char* argv[]){
    const FParameterNames LocalOptionP2PM2LCudaOnly { {"-p2p-m2l-cuda-only"},
        "To compute the P2P and the M2L on CUDA only (the rest is still "
        "on the CPU."};

    FHelpDescribeAndExit(argc, argv, "Perform Lagrange Kernel based simulation with CUDA+StarPU.",
                         FParameterDefinitions::OctreeHeight,
#ifdef RANDOM_PARTICLES
                         FParameterDefinitions::NbParticles,
#else
                         FParameterDefinitions::InputFile,
#endif
                         FParameterDefinitions::NbThreads,
                         LocalOptionBlocSize, LocalOptionValidate,
                         LocalOptionP2PM2LCudaOnly);

    if(FParameters::existParameter(argc, argv, LocalOptionP2PM2LCudaOnly.options)){
        std::cout << "/!\\ P2P and M2L only on GPU\n";
        return mainCore<FStarPUP2PM2LCudaOnlyCapacities>(argc, argv);
    }
    else{
        return mainCore<FStarPUCudaP2PM2LCapacities>(argc, argv);
    }
}
