// ==== CMAKE =====
// @FUSE_BLAS
// @FUSE_MPI
// ================
//
// ChebyshevHybridFMM.cpp
//
/** \brief Chebyshev FMM example
 *
 * \file
 * \authors O. Coulaud 
 *
 * This program runs the FMM Algorithm with the interpolation kernel based on
 * Chebyshev (grid points) interpolation (1/r kernel). It then compares the
 * results with a direct computation.
 */
#include <string> 
#include "Kernels/Chebyshev/FChebCell.hpp"

#include "Kernels/Chebyshev/FChebSymKernel_i.hpp"

template<typename FReal, int ORDER> 
using FInterpolationCell =  FChebCell<FReal, ORDER>;

template<typename FReal, typename GroupCellClass,
	 typename GroupContainerClass,
	 typename MatrixKernelClass, int ORDER>  
					
using FInterpolationKernel = FChebSymKernel_i<FReal,
					    GroupCellClass,
					    GroupContainerClass,
					    MatrixKernelClass,
					    ORDER> ;
						
const std::string interpolationType("Chebyshev interpolation");

#include "MPIInterpolationFMM.hpp" 
