

// @FUSE_STARPU
// @FUSE_BLAS
// @FUSE_FFT

// @SCALFMM_PRIVATE

#include "Utils/FGlobal.hpp"

#include "GroupTree/Core/FGroupTree.hpp"

#include "Components/FSimpleLeaf.hpp"
#include "Components/FSymbolicData.hpp"
#include "Containers/FVector.hpp"

#include "Kernels/P2P/FP2PParticleContainer.hpp"

#include "Kernels/Interpolation/FInterpMatrixKernel.hpp"
#include "Kernels/Uniform/FUnifKernel.hpp"
#include "Kernels/Uniform/FUnifCell.hpp"

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


#ifdef SCALFMM_USE_CUDA
#include "GroupTree/Cuda/FCudaDeviceWrapper.hpp"
#include "GroupTree/Cuda/FCudaEmptyCellSymb.hpp"
#include "GroupTree/Cuda/FCudaGroupOfParticles.hpp"
#include "GroupTree/Cuda/FCudaGroupOfCells.hpp"
#include "GroupTree/Uniform/FUnifCudaCellPOD.hpp"
#include "GroupTree/Uniform/FUnifCudaSharedData.hpp"
#include "GroupTree/Cuda/FCudaData.hpp"
#include "GroupTree/Cuda/FCudaTic.hpp"
template <class FReal, int ORDER>
class FUnifCuda;
#endif // SCALFMM_USE_CUDA

#include "Utils/FParameterNames.hpp"

#include "GroupTree/StarPUUtils/FStarPUCpuWrapper.hpp"
#include "GroupTree/StarPUUtils/FStarPUCptInteractionsWrapper.hpp"

#include <memory>

enum {
    TIME_P2M = 0,
    TIME_M2M,
    TIME_M2L,
    TIME_L2L,
    TIME_L2P,
    TIME_P2P,
    TIME_NB
};

const char* TIME_STR[TIME_NB] = {
    "TIME_P2M",
    "TIME_M2M",
    "TIME_M2L",
    "TIME_L2L",
    "TIME_L2P",
    "TIME_P2P"
};


template <int ORDER, class FReal>
void CheckCpu(FReal timePerOperators[TIME_NB], const FSize nbPartsPerLeaf, const int groupSize){
    typedef FInterpMatrixKernelR<FReal> MatrixKernelClass;

    using GroupCellClass     = FUnifCell<FReal, ORDER>;
    using GroupCellUpClass   = typename GroupCellClass::multipole_t;
    using GroupCellDownClass = typename GroupCellClass::local_expansion_t;
    using GroupCellSymbClass = FSymbolicData;

    typedef FP2PGroupParticleContainer<FReal>          GroupContainerClass;
    typedef FGroupTree< FReal, GroupCellSymbClass, GroupCellUpClass, GroupCellDownClass, GroupContainerClass, 1, 4, FReal>  GroupOctreeClass;

    typedef FUnifKernel<FReal,GroupCellClass,GroupContainerClass,MatrixKernelClass,ORDER> GroupKernelClass;
    typedef FStarPUCpuWrapper<typename GroupOctreeClass::CellGroupClass, GroupCellClass, GroupKernelClass, typename GroupOctreeClass::ParticleGroupClass, GroupContainerClass> GroupCpuWrapper;
    typedef FStarPUCptInteractionsWrapper<typename GroupOctreeClass::CellGroupClass, GroupCellClass, GroupKernelClass, typename GroupOctreeClass::ParticleGroupClass, GroupContainerClass> GroupCptInterWrapper;

    // Build the groups, find the minimum h to have a group of groupSize
    int NbLevels = 2;
    while( (1 << ((NbLevels-1-1)*3)) < groupSize ){
        NbLevels += 1;
    }

    const int nbLeavesInTree = ((NbLevels-1)*3);
    const FSize totalParts = nbLeavesInTree * nbPartsPerLeaf;

    FRandomLoader<FReal> loader(totalParts, 1.0, FPoint<FReal>(0,0,0), 0);

    FP2PParticleContainer<FReal> allParticles;
    for(FSize idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
        FPoint<FReal> particlePosition;
        FReal physicalValue;
        physicalValue = 0.10;
        loader.fillParticle(&particlePosition);
        allParticles.push(particlePosition, physicalValue);
    }

    GroupOctreeClass groupedTree(NbLevels, loader.getBoxWidth(), loader.getCenterOfBox(), groupSize, &allParticles);
    // groupedTree.printInfoBlocks();

    const MatrixKernelClass MatrixKernel;
    GroupKernelClass groupkernel(NbLevels, loader.getBoxWidth(), loader.getCenterOfBox(), &MatrixKernel);

    GroupCpuWrapper cpuWorker(NbLevels);
    cpuWorker.initKernel(0, &groupkernel);

    GroupCptInterWrapper cptWorker(NbLevels);
    cptWorker.initKernel(0, &groupkernel);
    cptWorker.compute(false);

    FTic timer;

    // We test the P2M
    timer.tic();
    cpuWorker.bottomPassPerform(groupedTree.getCellGroup(NbLevels-1,0), groupedTree.getParticleGroup(0));
    timer.tac();
    cptWorker.bottomPassPerform(groupedTree.getCellGroup(NbLevels-1,0), groupedTree.getParticleGroup(0));
    timePerOperators[TIME_P2M] = timer.elapsed()/double(cptWorker.getNbInteractions(0, GroupCptInterWrapper::INTER_P2M));

    // We test the M2M
    timer.tic();
    cpuWorker.upwardPassPerform(groupedTree.getCellGroup(NbLevels-2,0), groupedTree.getCellGroup(NbLevels-1,0), NbLevels-2);
    timer.tac();
    cptWorker.upwardPassPerform(groupedTree.getCellGroup(NbLevels-2,0), groupedTree.getCellGroup(NbLevels-1,0), NbLevels-2);
    timePerOperators[TIME_M2M] = timer.elapsed()/double(cptWorker.getNbInteractions(0, GroupCptInterWrapper::INTER_M2M));

    // We test the M2L
    timer.tic();
    cpuWorker.transferInPassPerform(groupedTree.getCellGroup(NbLevels-1,0), NbLevels-1);
    timer.tac();
    cptWorker.transferInPassPerform(groupedTree.getCellGroup(NbLevels-1,0), NbLevels-1);
    timePerOperators[TIME_M2L] = timer.elapsed()/double(cptWorker.getNbInteractions(0, GroupCptInterWrapper::INTER_M2L));

    // We test the L2L
    timer.tic();
    cpuWorker.downardPassPerform(groupedTree.getCellGroup(NbLevels-2,0), groupedTree.getCellGroup(NbLevels-1,0), NbLevels-2);
    timer.tac();
    cptWorker.downardPassPerform(groupedTree.getCellGroup(NbLevels-2,0), groupedTree.getCellGroup(NbLevels-1,0), NbLevels-2);
    timePerOperators[TIME_L2L] = timer.elapsed()/double(cptWorker.getNbInteractions(0, GroupCptInterWrapper::INTER_L2L));

    // We test the L2P
    timer.tic();
    cpuWorker.mergePassPerform(groupedTree.getCellGroup(NbLevels-1,0), groupedTree.getParticleGroup(0));
    timer.tac();
    cptWorker.mergePassPerform(groupedTree.getCellGroup(NbLevels-1,0), groupedTree.getParticleGroup(0));
    timePerOperators[TIME_L2P] = timer.elapsed()/double(cptWorker.getNbInteractions(0, GroupCptInterWrapper::INTER_L2P));

    // We test the P2P
    timer.tic();
    cpuWorker.directInPassPerform(groupedTree.getParticleGroup(0));
    timer.tac();
    cptWorker.directInPassPerform(groupedTree.getParticleGroup(0));
    timePerOperators[TIME_P2P] = timer.elapsed()/double(cptWorker.getNbInteractions(0, GroupCptInterWrapper::INTER_P2P));

    cptWorker.releaseKernel(0);
    cpuWorker.releaseKernel(0);
}

#ifdef SCALFMM_USE_CUDA

template <int ORDER, class FReal>
void CheckGpu(FReal timePerOperators[TIME_NB], const FSize nbPartsPerLeaf, const int groupSize){
    typedef FInterpMatrixKernelR<FReal> MatrixKernelClass;

    typedef FUnifCellPODCore         GroupCellSymbClass;
    typedef FUnifCellPODPole<FReal,ORDER>  GroupCellUpClass;
    typedef FUnifCellPODLocal<FReal,ORDER> GroupCellDownClass;
    typedef FUnifCellPOD<FReal,ORDER>      GroupCellClass;

    typedef FP2PGroupParticleContainer<FReal>          GroupContainerClass;
    typedef FGroupTree< FReal, GroupCellSymbClass, GroupCellUpClass, GroupCellDownClass, GroupContainerClass, 1, 4, FReal>  GroupOctreeClass;

    typedef FUnifKernel<FReal,GroupCellClass,GroupContainerClass,MatrixKernelClass,ORDER> GroupKernelClass;

    typedef FStarPUCudaWrapper<GroupKernelClass,
            FBasicCellPOD, FCudaUnifCellPODPole<FReal,ORDER>,FCudaUnifCellPODLocal<FReal,ORDER>,
            FCudaGroupOfCells<FBasicCellPOD, FCudaUnifCellPODPole<FReal,ORDER>,FCudaUnifCellPODLocal<FReal,ORDER> >,
            FCudaGroupOfParticles<FReal, 1, 4, FReal>, FCudaGroupAttachedLeaf<FReal, 1, 4, FReal>, FUnifCuda<FReal,ORDER> > GroupCudaWrapper;

    typedef FStarPUCptInteractionsWrapper<typename GroupOctreeClass::CellGroupClass, GroupCellClass, GroupKernelClass, typename GroupOctreeClass::ParticleGroupClass, GroupContainerClass> GroupCptInterWrapper;

    // Build the groups, find the minimum h to have a group of groupSize
    int NbLevels = 2;
    while( (1 << ((NbLevels-1-1)*3)) < groupSize ){
        NbLevels += 1;
    }

    const int nbLeavesInTree = ((NbLevels-1)*3);
    const FSize totalParts = nbLeavesInTree * nbPartsPerLeaf;

    FRandomLoader<FReal> loader(totalParts, 1.0, FPoint<FReal>(0,0,0), 0);

    FP2PParticleContainer<FReal> allParticles;
    for(FSize idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
        FPoint<FReal> particlePosition;
        FReal physicalValue;
        physicalValue = 0.10;
        loader.fillParticle(&particlePosition);
        allParticles.push(particlePosition, physicalValue);
    }

    GroupOctreeClass groupedTree(NbLevels, loader.getBoxWidth(), loader.getCenterOfBox(), groupSize, &allParticles);
    // groupedTree.printInfoBlocks();

    const MatrixKernelClass MatrixKernel;
    GroupKernelClass groupkernel(NbLevels, loader.getBoxWidth(), loader.getCenterOfBox(), &MatrixKernel);

    GroupCudaWrapper gpuWorker(NbLevels);
    gpuWorker.initKernel(0, &groupkernel);

    GroupCptInterWrapper cptWorker(NbLevels);
    cptWorker.initKernel(0, &groupkernel);
    cptWorker.compute(false);

    FCudaCheck(cudaSetDevice(0));
    cudaStream_t currentStream;
    FCudaCheck(cudaStreamCreate(&currentStream));

    // We test the M2L
    {
        typename GroupOctreeClass::CellGroupClass* cells = groupedTree.getCellGroup(NbLevels-1,0);
        // gpuWorker.transferInPassPerform(groupedTree.getCellGroup(NbLevels-1,0), NbLevels-1);
        FCudaData<unsigned char> rawBuffer((const unsigned char*)cells->getRawBuffer(),
                                           cells->getBufferSizeInByte());
        FCudaData<unsigned char> rawBufferPole((const unsigned char*)cells->getRawMultipoleBuffer(),
                                           cells->getMultipoleBufferSizeInByte());
        FCudaData<unsigned char> rawBufferLocal((const unsigned char*)cells->getRawLocalBuffer(),
                                           cells->getLocalBufferSizeInByte());


        FCudaTic timer(currentStream);

        FCuda__transferInPassCallback< FBasicCellPOD, FCudaUnifCellPODPole<FReal,ORDER>,FCudaUnifCellPODLocal<FReal,ORDER>,
                FCudaGroupOfCells<FBasicCellPOD, FCudaUnifCellPODPole<FReal,ORDER>,FCudaUnifCellPODLocal<FReal,ORDER> >,
                FCudaGroupOfParticles<FReal, 1, 4, FReal>,
                FCudaGroupAttachedLeaf<FReal, 1, 4, FReal>,
                FUnifCuda<FReal,ORDER> >(
                    rawBuffer.get(), rawBuffer.getSize(),
                    rawBufferPole.get(),
                    rawBufferLocal.get(),
                    NbLevels-1, gpuWorker.getKernel(0), currentStream,
                    FCuda__GetGridSize(gpuWorker.getKernel(0),cells->getNumberOfCellsInBlock()),
                    FCuda__GetBlockSize(gpuWorker.getKernel(0)));

        timer.tac();
        cptWorker.transferInPassPerform(cells, NbLevels-1);
        timePerOperators[TIME_M2L] = timer.elapsed()/double(cptWorker.getNbInteractions(0, GroupCptInterWrapper::INTER_M2L));
    }

    // We test the P2P
    {
        typename GroupOctreeClass::ParticleGroupClass* particles = groupedTree.getParticleGroup(0);
        //gpuWorker.directInPassPerform(groupedTree.getParticleGroup(0));

        FCudaData<unsigned char> rawBuffer((const unsigned char*)particles->getRawBuffer(),
                                           particles->getBufferSizeInByte());
        FCudaData<unsigned char> rawAttributesBuffer((const unsigned char*)particles->getRawAttributesBuffer(),
                                           particles->getAttributesBufferSizeInByte());

        FCudaTic timer(currentStream);

        FCuda__directInPassCallback< FBasicCellPOD, FCudaUnifCellPODPole<FReal,ORDER>,FCudaUnifCellPODLocal<FReal,ORDER>,
                FCudaGroupOfCells<FBasicCellPOD, FCudaUnifCellPODPole<FReal,ORDER>,FCudaUnifCellPODLocal<FReal,ORDER> >,
                FCudaGroupOfParticles<FReal, 1, 4, FReal>,
                FCudaGroupAttachedLeaf<FReal, 1, 4, FReal>,
                FUnifCuda<FReal,ORDER>>(
                    rawBuffer.get(), rawBuffer.getSize(),
                    rawAttributesBuffer.get(),
                    NbLevels, gpuWorker.getKernel(0), currentStream,
                    FCuda__GetGridSize(gpuWorker.getKernel(0),particles->getNumberOfLeavesInBlock()),
                    FCuda__GetBlockSize(gpuWorker.getKernel(0)));

        timer.tac();
        cptWorker.directInPassPerform(particles);
        timePerOperators[TIME_P2P] = timer.elapsed()/double(cptWorker.getNbInteractions(0, GroupCptInterWrapper::INTER_P2P));
    }

    cptWorker.releaseKernel(0);
    gpuWorker.releaseKernel(0);
}

#endif // SCALFMM_USE_CUDA


int main(int argc, char** argv){
    typedef double FReal;

    {
        const int startGroupSize = 100;
        const int endGroupSize   = 10000;
        const int nbStepsSize    = 10;

        std::cout << "Operations" << "\t";
        for(int idx = 0 ; idx < TIME_NB ; ++idx){
            std::cout << TIME_STR[idx] << "\t";
        }
        std::cout << std::endl;

        const int startNbParts = 50;
        const int endNbParts   = 1000;
        const int nbStepsNbParts = 5;

        for(int idxNbPart = startNbParts ; idxNbPart <= endNbParts ; idxNbPart += (endNbParts-startNbParts)/nbStepsNbParts){
            for(int idxGroupSize = startGroupSize ; idxGroupSize <= endGroupSize ; idxGroupSize += (endGroupSize-startGroupSize)/nbStepsSize){
                FReal timePerOperators[TIME_NB] = {0};
                CheckCpu<5,FReal>(timePerOperators, idxNbPart, idxGroupSize);

                std::cout << "(" << idxNbPart << "P/" << idxGroupSize << "G)\t";
                for(int idx = 0 ; idx < TIME_NB ; ++idx){
                    std::cout << timePerOperators[idx] << "\t";
                }
                std::cout << std::endl;
            }
        }
    }

#ifdef SCALFMM_USE_CUDA

    {
        const int startGroupSize = 100;
        const int endGroupSize   = 4000;
        const int nbStepsSize    = 10;

        std::cout << "Operations" << "\t";
        for(int idx = 0 ; idx < TIME_NB ; ++idx){
            std::cout << TIME_STR[idx] << "\t";
        }
        std::cout << std::endl;

        const int startNbParts = 100;
        const int endNbParts   = 10000;
        const int nbStepsNbParts = 10;

        for(int idxNbPart = startNbParts ; idxNbPart <= endNbParts ; idxNbPart += (endNbParts-startNbParts)/nbStepsNbParts){
            for(int idxGroupSize = startGroupSize ; idxGroupSize <= endGroupSize ; idxGroupSize += (endGroupSize-startGroupSize)/nbStepsSize){
                FReal timePerOperators[TIME_NB] = {0};
                CheckGpu<5,FReal>(timePerOperators, idxNbPart, idxGroupSize);

                std::cout << "(" << idxNbPart << "P/" << idxGroupSize << "G)\t";
                for(int idx = 0 ; idx < TIME_NB ; ++idx){
                    std::cout << timePerOperators[idx] << "\t";
                }
                std::cout << std::endl;
            }
        }
    }

#endif // SCALFMM_USE_CUDA

    return 0;
}
