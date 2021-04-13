#ifndef FADAPT_UNIF_KERNEL_HPP
#define FADAPT_UNIF_KERNEL_HPP

#include <cassert>

#include "Kernels/Uniform/FUnifKernel.hpp"
#include "InastempCompileConfig.h"

#include "Utils/FMath.hpp"

#include <fstream>


/**
 * \brief Adaptive FMM Lagrange kernel adaptor.
 *
 * Inherits from the original Lagrange kernel and adds the needed operators to
 * make it adaptive.
 */
template<
    class FReal,
    class CellClass,
    class ContainerClass,
    class MatrixKernelClass,
    int ORDER,
    int NVALS = 1
    >
class FAdaptUnifKernel :
    public FUnifKernel<FReal, CellClass, ContainerClass, MatrixKernelClass, ORDER, NVALS>
{
    using FBase =  FUnifKernel<FReal, CellClass, ContainerClass, MatrixKernelClass, ORDER, NVALS>;

public:

    using FBase::FBase;

    using FBase::M2M;
    using FBase::M2L;
    using FBase::L2L;
    using FBase::P2P;
    using FBase::P2M;
    using FBase::L2P;

    // /**
    //  * \brief Adaptive P2M
    //  *
    //  * Computes the leaf box width for every leaf. The uniform one relies on the
    //  * fact that all leaves are at the same level to pre-compute this width.
    //  */
    // template<class SymbolicData>
    // void P2M(typename CellClass::multipole_t* const LeafMultipole,
    //          const SymbolicData* const LeafSymbData,
    //          const ContainerClass* const SourceParticles)
    // {
    //     FReal leafBoxWidth = FBase::BoxWidth / FReal(1 <<  LeafSymbData->depth);

    //     const FPoint<FReal> LeafCellCenter(
    //         FBase::getCellCenter(LeafSymbData->getCoordinate(),
    //                              static_cast<int>(LeafSymbData->depth)));

    //     FBase::Interpolator->applyP2M(LeafCellCenter, leafBoxWidth,
    //                                   LeafMultipole->get(0), SourceParticles);

    //     for(int idxRhs = 0 ; idxRhs < NVALS ; ++idxRhs){
    //         // 2) apply Discrete Fourier Transform
    //         this->M2LHandler.applyZeroPaddingAndDFT(LeafMultipole->get(idxRhs),
    //                                                 LeafMultipole->getTransformed(idxRhs));
    //     }
    // }

    // /**
    //  * \brief Adaptive L2P
    //  *
    //  * Computes the leaf box width for every leaf. The uniform one relies on the
    //  * fact that all leaves are at the same level to pre-compute this width.
    //  */
    // template<class SymbolicData>
    // void L2P(const typename CellClass::local_expansion_t* const LeafCell,
    //          const SymbolicData* const LeafSymbData,
    //          ContainerClass* const TargetParticles)
    // {
    //     FReal leafBoxWidth = FBase::BoxWidth / FReal(1 <<  LeafSymbData->depth);

    //     const FPoint<FReal> LeafCellCenter(
    //         FBase::getCellCenter(LeafSymbData->getCoordinate(),
    //                              static_cast<int>(LeafSymbData->depth)));

    //     FReal localExp[NVALS*FBase::nnodes];

    //     for(int idxRhs = 0 ; idxRhs < NVALS ; ++idxRhs){
    //         // 1)  Apply Inverse Discete Fourier Transform
    //         FBase::M2LHandler.unapplyZeroPaddingAndDFT(LeafCell->getTransformed(idxRhs),
    //                                             localExp + idxRhs*FBase::nnodes);
    //         FBlas::add(FBase::nnodes,const_cast<FReal*>(LeafCell->get(idxRhs)),localExp + idxRhs*FBase::nnodes);
    //     }
    //     // 2.a) apply Sx
    //     FBase::Interpolator->applyL2P(LeafCellCenter, leafBoxWidth,
    //                                   localExp, TargetParticles);
    //     // 2.b) apply Px (grad Sx)
    //     FBase::Interpolator->applyL2PGradient(LeafCellCenter, leafBoxWidth,
    //                                           localExp, TargetParticles);
    // }

    /**
     * \brief P2L operator
     *
     * Computes interactions from a leaf particles to a neighbouring local
     * expansion.
     */
    template<class SymbolicData>
    void P2L(typename CellClass::local_expansion_t* const local,
             const SymbolicData * const symb,
             const ContainerClass* const particles,
             const SymbolicData * const /*source_symb*/)
    {
        using ComputeClass = InaVecBestType<FReal>;
        constexpr int FRealCount = ComputeClass::VecLength;

        // Target cell: local
        const FReal localCellWidth(FBase::BoxWidth / FReal(1 << symb->getLevel()));
        const FPoint<FReal> localCellCenter =
            FBase::getCellCenter(symb->getCoordinate(),static_cast<int>(symb->getLevel()));

        // interpolation points of target (X) cell
        FPoint<FReal> X[FBase::nnodes];
        FUnifTensor<FReal,ORDER>::setRoots(localCellCenter, localCellWidth, X);

        // Particles attributes
        const FReal * const posX = particles->getPositions()[0];
        const FReal * const posY = particles->getPositions()[1];
        const FReal * const posZ = particles->getPositions()[2];
        const FReal * const physicalValues = particles->getPhysicalValues();

//        const FReal* pX = particles->getPositions()[0];
//        const FReal* pY = particles->getPositions()[1];
//        const FReal* pZ = particles->getPositions()[2];
//        const FReal* pV = particles->getPhysicalValues();

        // apply P2L
        for(int idxRhs = 0 ; idxRhs < NVALS ; ++idxRhs){
            for (unsigned int m = 0; m < FBase::nnodes; ++m) {
                ComputeClass XX = ComputeClass(X[m].getX());
                ComputeClass XY = ComputeClass(X[m].getY());
                ComputeClass XZ = ComputeClass(X[m].getZ());

                ComputeClass tmpLocalExp = ComputeClass::GetZero();
                // Compute using vectorization for all but the last array elements

                for (std::size_t idxPart = 0 ; idxPart < particles->getNbParticles() ; idxPart += FRealCount)
                {
                    tmpLocalExp +=
                        FBase::MatrixKernel->evaluate(
                            XX, XY, XZ,
                            ComputeClass(&posX[idxPart]), ComputeClass(&posY[idxPart]), ComputeClass(&posZ[idxPart]))
                        * physicalValues[idxPart];
                }

                local->get(idxRhs)[m] += (tmpLocalExp.horizontalSum());
            }
        }// NVALS
    }

    /**
     * \brief M2P operator
     *
     * Computes interactions from a multipole to a leaf.
     */
    template<class SymbolicData>
    void M2P(const typename CellClass::multipole_t* const pole,
             const SymbolicData* const symb,
             ContainerClass* const particles,
             const SymbolicData * const /*target_symb*/)
    {
        using ComputeClass = InaVecBestType<FReal>;
        constexpr int FRealCount = ComputeClass::VecLength;

        // Source cell: pole
        const FReal poleCellWidth(FBase::BoxWidth / FReal(1 << symb->getLevel()));
        const FPoint<FReal> poleCellCenter(
            FBase::getCellCenter(symb->getCoordinate(),static_cast<int>(symb->getLevel())));

        // interpolation points of source (Y) cell
        FPoint<FReal> Y[FBase::nnodes];
        FUnifTensor<FReal,ORDER>::setRoots(poleCellCenter, poleCellWidth, Y);

        // read positions
        const FReal* const posX = (particles->getPositions()[0]);
        const FReal* const posY = (particles->getPositions()[1]);
        const FReal* const posZ = (particles->getPositions()[2]);

        // get potential
        FReal* const physVal = (particles->getPhysicalValues());
        FReal* const potentials = (particles->getPotentials());
        FReal* const fx = (particles->getForcesX());
        FReal* const fy = (particles->getForcesY());
        FReal* const fz = (particles->getForcesZ());

        for(int idxRhs = 0 ; idxRhs < NVALS ; ++idxRhs){

            // apply M2P
            for (unsigned int n=0; n<FBase::nnodes; ++n){

                ComputeClass  MultipoleExpansion =
                    ComputeClass(pole->get(idxRhs)[n]);

                ComputeClass YX = ComputeClass(Y[n].getX());
                ComputeClass YY = ComputeClass(Y[n].getY());
                ComputeClass YZ = ComputeClass(Y[n].getZ());

                for(std::size_t idxPart = 0;
                    idxPart < particles->getNbParticles();
                    idxPart += FRealCount)
                {

                    ComputeClass Kxy[1];
                    ComputeClass dKxy[3];
                    FBase::MatrixKernel->evaluateBlockAndDerivative(
                        ComputeClass(&posX[idxPart]),
                        ComputeClass(&posY[idxPart]),
                        ComputeClass(&posZ[idxPart]),
                        YX, YY, YZ,
                        Kxy,dKxy);

                    (ComputeClass(&potentials[idxPart]) + Kxy[0] * MultipoleExpansion).storeInArray(&potentials[idxPart]);
                    (ComputeClass(&fx[idxPart]) + dKxy[0] * physVal[idxPart] * MultipoleExpansion).storeInArray(&fx[idxPart]);
                    (ComputeClass(&fy[idxPart]) + dKxy[1] * physVal[idxPart] * MultipoleExpansion).storeInArray(&fy[idxPart]);
                    (ComputeClass(&fz[idxPart]) + dKxy[2] * physVal[idxPart] * MultipoleExpansion).storeInArray(&fz[idxPart]);

                }

            }// Particles

        }// NVALS
    }
};



#endif /* FADAPT_UNIF_KERNEL_HPP */
