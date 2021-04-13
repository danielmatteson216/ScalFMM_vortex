// ===================================================================================
// Copyright ScalFmm 2011 INRIA, Olivier Coulaud, Berenger Bramas
// olivier.coulaud@inria.fr, berenger.bramas@inria.fr
// This software is a computer program whose purpose is to compute the FMM.
//
// This software is governed by the CeCILL-C and LGPL licenses and
// abiding by the rules of distribution of free software.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public and CeCILL-C Licenses for more details.
// "http://www.cecill.info".
// "http://www.gnu.org/licenses".
// ===================================================================================

// ==== CMAKE =====
// @FUSE_FFT
// @FUSE_BLAS
//
//  ==== Git =====
//
// ================

/** \brief Chebyshev adaptive FMM example
 *
 * \file
 * \authors  O. Coulaud
 *
 * This program runs the Adaptive FMM Algorithm with the interpolation kernel based on
 * uniform (grid points) interpolation (1/r kernel). It then compares the
 * results with a direct computation.
 */

#include <string>
#include "Kernels/Chebyshev/FChebCell.hpp"
#include "Adaptive/FAdaptChebKernel.hpp"
//
template<typename FReal, int ORDER>
using FInterpolationCell =  FChebCell<FReal, ORDER>;

template<typename FReal, typename GroupCellClass,
	 typename GroupContainerClass,
	 typename MatrixKernelClass, int ORDER>
using FInterpolationAdaptiveKernel = FAdaptChebKernel<FReal,
					    GroupCellClass,
					    GroupContainerClass,
					    MatrixKernelClass,
					    ORDER> ;

const std::string interpolationType("Chebyshev interpolation");

#include "sharedMemoryInterpolationAdaptiveFMM.hpp"
