#ifndef _F_BLOCKED_MPI_INTERPOLATION_HPP_
#define _F_BLOCKED_MPI_INTERPOLATION_HPP_


#include "Utils/FGlobal.hpp"

#include "GroupTree/Core/FGroupTree.hpp"

#include "Components/FSimpleLeaf.hpp"
#include "Components/FSymbolicData.hpp"
#include "Containers/FVector.hpp"

#include "Kernels/P2P/FP2PParticleContainer.hpp"

#include "Utils/FMath.hpp"
#include "Utils/FMemUtils.hpp"
#include "Utils/FParameters.hpp"

#include "Files/FRandomLoader.hpp"
#include "Files/FFmaGenericLoader.hpp"

#include "GroupTree/Core/FGroupSeqAlgorithm.hpp"
#include "GroupTree/Core/FGroupTaskAlgorithm.hpp"
#include "GroupTree/Core/FGroupTaskStarpuAlgorithm.hpp"
#include "GroupTree/Core/FP2PGroupParticleContainer.hpp"

#include "Utils/FParameterNames.hpp"

#include "Components/FTestParticleContainer.hpp"
#include "Components/FTestCell.hpp"
#include "Components/FTestKernels.hpp"

#include "Core/FFmmAlgorithmThreadProc.hpp"
#include "Files/FMpiTreeBuilder.hpp"
#include "GroupTree/Core/FGroupTaskStarpuMpiAlgorithm.hpp"
#include "Core/FFmmAlgorithm.hpp" //For validation

#include "Files/FMpiFmaGenericLoader.hpp"
#include "Containers/FCoordinateComputer.hpp"

#include "GroupTree/StarPUUtils/FStarPUKernelCapacities.hpp"

#include <memory>
using namespace std;

namespace blockedMpiInterpolation{

//Function header
void timeAverage(int mpi_rank, int nproc, double elapsedTime);
FSize getNbParticlesPerNode(FSize mpi_count, FSize mpi_rank, FSize total);

template<
    class GroupCellClass,
    class GroupCellUpClass,
    class GroupCellDownClass,
    class GroupCellSymbClass,
    class KernelClass,
    class MatrixKernelClass
    >
auto execute_algorithm(int argc, char* argv[]){
    //Define parameters
    const FParameterNames LocalOptionBlocSize { {"-bs"}, "The size of the block of the blocked tree"};
    const FParameterNames LocalOptionEllipsoid = {{"-ellipsoid"} , " non uniform distribution on  an ellipsoid of aspect ratio given by a=0.5 b=0.25 c=0.125"};
    const FParameterNames LocalOptionCube = {{"-cube", "-uniform"} , " uniform distribution on cube (default)"};
    // Define types
    using FReal = double;
    using GroupContainerClass =
        FP2PGroupParticleContainer<FReal>;
    using GroupOctreeClass =
        FGroupTree< FReal, GroupCellSymbClass, GroupCellUpClass, GroupCellDownClass, GroupContainerClass, 1, 4, FReal>;
    using GroupKernelClass =
        FStarPUAllCpuCapacities<KernelClass>;
    using GroupCpuWrapper =
        FStarPUCpuWrapper<typename GroupOctreeClass::CellGroupClass, GroupCellClass, GroupKernelClass, typename GroupOctreeClass::ParticleGroupClass, GroupContainerClass> ;
    using GroupAlgorithm =
        FGroupTaskStarPUMpiAlgorithm<GroupOctreeClass, typename GroupOctreeClass::CellGroupClass, GroupKernelClass, typename GroupOctreeClass::ParticleGroupClass, GroupCpuWrapper> ;
    // Init MPI_COM
    FMpi mpiComm(argc,argv);

    // Init timer
    FTic timer;

    // Getting parameters
    const int groupSize =
        FParameters::getValue(argc,argv,LocalOptionBlocSize.options, 250);
    const unsigned int TreeHeight    =
        FParameters::getValue(argc, argv, FParameterDefinitions::OctreeHeight.options, 5);
        
        const FSize totalNbParticles =
        FParameters::getValue(argc,argv,FParameterDefinitions::NbParticles.options, FSize(20));

    const FSize NbParticles   =
        getNbParticlesPerNode(mpiComm.global().processCount(), mpiComm.global().processId(), totalNbParticles);

    // init particles position and physical value
    struct TestParticle{
        FPoint<FReal> position;
        FReal physicalValue;
        const FPoint<FReal>& getPosition(){
            return position;
        }
		const unsigned int getWriteDataSize(void) const {
			return sizeof(FReal);
		}
		const unsigned int getWriteDataNumber(void) const {
			return 3;
		}
		const FReal* getPtrFirstData(void) const {
			return position.data();
		}
    };

    // LOADING PARTICLE
    #ifndef LOAD_FILE
        FRandomLoader<FReal> loader(NbParticles, 1.0, FPoint<FReal>(0,0,0), mpiComm.global().processId());
        FAssertLF(loader.isOpen());
        TestParticle* allParticles = new TestParticle[loader.getNumberOfParticles()];
        memset(allParticles,0,(unsigned int) (sizeof(TestParticle)* loader.getNumberOfParticles()));
        for(FSize idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
            loader.fillParticle(&allParticles[idxPart].position);
            allParticles[idxPart].physicalValue = 0.1;
        }
    #else
        const char* const filename = FParameters::getStr(argc,argv,FParameterDefinitions::InputFile.options, "../Data/test20k.fma");
        FMpiFmaGenericLoader<FReal> loader(filename,mpiComm.global());
        FAssertLF(loader.isOpen());
        TestParticle* allParticles = new TestParticle[loader.getMyNumberOfParticles()];
        memset(allParticles,0,(unsigned int) (sizeof(TestParticle)* loader.getMyNumberOfParticles()));
        for(FSize idxPart = 0 ; idxPart < loader.getMyNumberOfParticles() ; ++idxPart){
            loader.fillParticle(&allParticles[idxPart].position,&allParticles[idxPart].physicalValue);
        }
    #endif

    FVector<TestParticle> myParticles;
    FLeafBalance balancer;
    FMpiTreeBuilder< FReal,
                    TestParticle >::DistributeArrayToContainer(
                                        mpiComm.global(),
                                        allParticles,
                                        loader.getNumberOfParticles(),
                                        loader.getCenterOfBox(),
                                        loader.getBoxWidth(),
                                        TreeHeight,
                                        &myParticles,
                                        &balancer);

    // Each proc need to know the righest morton index
    const FTreeCoordinate host = FCoordinateComputer::GetCoordinateFromPosition<FReal>(
                loader.getCenterOfBox(),
                loader.getBoxWidth(),
                TreeHeight,
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
    FP2PParticleContainer<FReal> myParticlesInContainer;
    for(FSize idxPart = 0 ; idxPart < myParticles.getSize() ; ++idxPart){
        myParticlesInContainer.push(myParticles[idxPart].position,
                                    myParticles[idxPart].physicalValue);
    }
    GroupOctreeClass groupedTree(TreeHeight, loader.getBoxWidth(), loader.getCenterOfBox(), groupSize,
                                 &myParticlesInContainer, true, leftLimite);
    timer.tac();
	std::cerr << "Done  " << "(@Creating and Inserting Particles = " << timer.elapsed() << "s)." << std::endl;

    int operationsToProceed =  FFmmP2P | FFmmP2M | FFmmM2M | FFmmM2L | FFmmL2L | FFmmL2P;
    { // -----------------------------------------------------


        const MatrixKernelClass MatrixKernel;
        // Create Matrix Kernel
        GroupKernelClass groupkernel(TreeHeight, loader.getBoxWidth(), loader.getCenterOfBox(), &MatrixKernel);
        // Run the algorithm
        GroupAlgorithm groupalgo(mpiComm.global(), &groupedTree,&groupkernel);
		mpiComm.global().barrier();
        timer.tic();
		starpu_fxt_start_profiling();
        groupalgo.execute(operationsToProceed);
		mpiComm.global().barrier();
		starpu_fxt_stop_profiling();
        timer.tac();
		timeAverage(mpiComm.global().processId(), mpiComm.global().processCount(), timer.elapsed());
    }
    return &groupedTree;
}

void timeAverage(int mpi_rank, int nproc, double elapsedTime){
    if(mpi_rank == 0){
                double sumElapsedTimeMin = elapsedTime;
                double sumElapsedTimeMax = elapsedTime;
        for(int i = 1; i < nproc; ++i){
            double tmp;
            MPI_Recv(&tmp, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            if(tmp < sumElapsedTimeMin)
                sumElapsedTimeMin = tmp;
            if(tmp > sumElapsedTimeMax)
                sumElapsedTimeMax = tmp;
        }
        std::cout << "Min time per node (MPI)  : " << sumElapsedTimeMin << "s" << std::endl;
        std::cout << "Max time per node (MPI)  : " << sumElapsedTimeMax << "s" << std::endl;
    } else {
        MPI_Send(&elapsedTime, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

FSize getNbParticlesPerNode(FSize mpi_count, FSize mpi_rank, FSize total){
	if(mpi_rank < (total%mpi_count))
		return ((total - (total%mpi_count))/mpi_count)+1;
	return ((total - (total%mpi_count))/mpi_count);
}

}

#endif
