// ==== CMAKE =====
// @FUSE_BLAS
// ================
// Keep in private GIT
//
//
/** \brief Chebyshev FMM example
 *
 * \file
 * \authors B. Bramas, O. Coulaud
 *
 * This program runs the FMM Algorithm with the interpolation kernel based on
 * Chebyshev (grid points) interpolation (1/r kernel). It then compares the
 * results with a direct computation.
 */
#include <string> 
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
const std::string interpolationType("Chebyshev interpolation");

#include "sharedMemoryInterpolationFMM.hpp" 
