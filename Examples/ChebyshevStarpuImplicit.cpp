// ==== CMAKE =====
// @FUSE_BLAS
// @FUSE_MPI
// @FUSE_STARPU
//
//
#include <string>
//
// Chebychev cell class
#include "Kernels/Chebyshev/FChebCell.hpp"
#include "Kernels/Chebyshev/FChebSymKernel.hpp"
//
template<typename FReal, int ORDER> 
using FInterpolationCell =  FChebCell<FReal, ORDER>;

template<typename FReal, typename GroupCellClass,
	 typename GroupContainerClass,
	 typename MatrixKernelClass, int ORDER>  

using FInterpolationKernel = FChebSymKernel<FReal,
					    GroupCellClass,
					    GroupContainerClass,
					    MatrixKernelClass,
					    ORDER> ;

static std::string interpolationKernel("Chebyshev");
#include "genericStarPUImplicit.hpp"

