// ==== CMAKE =====
// @FUSE_BLAS
// @FUSE_MPI
// @FUSE_STARPU
//
#include <string>
//
// Uniform Grid points cell class
#include "Kernels/Uniform/FUnifCell.hpp"
#include "Kernels/Uniform/FUnifKernel.hpp"
//
template<typename FReal, int ORDER> 
using FInterpolationCell =  FUnifCell<FReal, ORDER>;

template<typename FReal, typename GroupCellClass,
	 typename GroupContainerClass,
	 typename MatrixKernelClass, int ORDER>  
using FInterpolationKernel = FUnifKernel<FReal,
					    GroupCellClass,
					    GroupContainerClass,
					    MatrixKernelClass,
					    ORDER> ;

static std::string interpolationKernel("Uniform");

#include "genericStarPUImplicit.hpp"

