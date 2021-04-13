// ==== CMAKE =====
// @FUSE_FFT
// @FUSE_BLAS
//  ==== Git =====
//
// ================

/** \brief equispaced point (Uniform or Lagrange) Adaptive FMM example
 *
 * \file
 * \authors  O. Coulaud
 *
 * This program runs the Adaptive FMM Algorithm with the interpolation kernel based on
 * uniform (grid points) interpolation (1/r kernel). It then compares the
 * results with a direct computation.
 */

#include <string>
#include "Kernels/Uniform/FUnifCell.hpp"
#include "Adaptive/FAdaptUnifKernel.hpp"
//
template<typename FReal, int ORDER>
using FInterpolationCell =  FUnifCell<FReal, ORDER>;

template<typename FReal, typename GroupCellClass,
	 typename GroupContainerClass,
	 typename MatrixKernelClass, int ORDER>
using FInterpolationAdaptiveKernel = FAdaptUnifKernel<FReal,
					    GroupCellClass,
					    GroupContainerClass,
					    MatrixKernelClass,
					    ORDER> ;

const std::string interpolationType("Lagrange equispaced points interpolation");

#include "sharedMemoryInterpolationAdaptiveFMM.hpp"
