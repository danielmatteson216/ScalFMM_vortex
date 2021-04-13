// ==== CMAKE =====
// @FUSE_BLAS
// ================
// Keep in private GIT
// @FUSE_MPI
// @FUSE_STARPU

// A supprimer ou merger Unif et Cheb
#include "Files/FBlockedMpiInterpolation.hpp"
#include "Utils/FGlobal.hpp"


#include "Kernels/P2P/FP2PParticleContainer.hpp"

#include "Kernels/Chebyshev/FChebSymKernel.hpp"
#include "Kernels/Chebyshev/FChebCell.hpp"
#include "Kernels/Interpolation/FInterpMatrixKernel.hpp"

#include "GroupTree/Core/FGroupSeqAlgorithm.hpp"
#include "GroupTree/Core/FGroupTaskAlgorithm.hpp"
#include "GroupTree/Core/FGroupTaskStarpuAlgorithm.hpp"
#include "GroupTree/Core/FP2PGroupParticleContainer.hpp"

#include "Components/FTestParticleContainer.hpp"
#include "Components/FTestCell.hpp"
#include "Components/FTestKernels.hpp"

#include "Core/FFmmAlgorithmThreadProc.hpp"
#include "Files/FMpiTreeBuilder.hpp"
#include "GroupTree/Core/FGroupTaskStarpuMpiAlgorithm.hpp"

#include "Files/FMpiFmaGenericLoader.hpp"
#include "Containers/FCoordinateComputer.hpp"

#include "GroupTree/StarPUUtils/FStarPUKernelCapacities.hpp"

#include <memory>


int main(int argc, char *argv[]){
    const FParameterNames LocalOptionBlocSize { {"-bs"}, "The size of the block of the blocked tree"};
    const FParameterNames LocalOptionNoValidate { {"-no-validation"}, "To avoid comparing with direct computation"};
    FHelpDescribeAndExit(argc, argv, "Test the blocked tree by counting the particles.",
                         FParameterDefinitions::OctreeHeight,FParameterDefinitions::InputFile,
                         FParameterDefinitions::NbParticles,
                         LocalOptionBlocSize,LocalOptionNoValidate);

    // using FReal = double;
    // static const int ORDER = 6;
    // using GroupContainerClass = FP2PGroupParticleContainer<FReal>;
    // using MatrixKernelClass   =  FInterpMatrixKernelR<FReal>;
    // using GroupCellClass      = FChebCell<FReal, ORDER>;
    // using GroupCellUpClass    = typename GroupCellClass::multipole_t;
    // using GroupCellDownClass  = typename GroupCellClass::local_expansion_t;
    // using GroupCellSymbClass = FSymbolicData;
    // using kernelClass = FChebSymKernel<FReal,GroupCellClass,GroupContainerClass,MatrixKernelClass,ORDER>;

    // auto groupedTree = blockedMpiInterpolation::execute_algorithm<
    //   GroupCellClass,
    //   GroupCellUpClass,
    //   GroupCellDownClass,
    //   GroupCellSymbClass,
    //   kernelClass,
    //   MatrixKernelClass
    //   >(argc,argv);

    // Validation

}
