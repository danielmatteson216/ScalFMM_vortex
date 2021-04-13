// ==== CMAKE =====
// @FUSE_BLAS
// @FUSE_FFT
// ================
// Keep in private GIT
//
//
/** \brief Uniform FMM example
 *
 * \file
 * \authors B. Bramas, O. Coulaud
 *
 * This program runs the FMM Algorithm with the interpolation kernel based on
 * uniform (grid points) interpolation (1/r kernel). It then compares the
 * results with a direct computation.
 */
#include <string> 
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
const std::string interpolationType("Uniform interpolation");

#include "sharedMemoryInterpolationFMM.hpp" 
