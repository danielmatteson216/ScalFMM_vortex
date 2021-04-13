// See LICENCE file at project root
#ifndef FCHEBTENSORIALKERNEL_HPP
#define FCHEBTENSORIALKERNEL_HPP

#include "Utils/FGlobal.hpp"

#include "Utils/FSmartPointer.hpp"

#include "Kernels/Chebyshev/FAbstractChebKernel.hpp"
#include "Kernels/Chebyshev/FChebTensorialM2LHandler.hpp" //PB: temporary version

class FTreeCoordinate;

/**
 * @author Matthias Messner(matthias.messner@inria.fr)
 * @class FChebTensorialKernel
 * @brief
 * Please read the license
 *
 * This kernels implement the Chebyshev interpolation based FMM operators. It
 * implements all interfaces (P2P, P2M, M2M, M2L, L2L, L2P) which are required by
 * the FFmmAlgorithm and FFmmAlgorithmThread.
 *
 * @tparam CellClass Type of cell
 * @tparam ContainerClass Type of container to store particles
 * @tparam MatrixKernelClass Type of matrix kernel function
 * @tparam ORDER Chebyshev interpolation order
 */
template < class FReal, class CellClass, class ContainerClass,   class MatrixKernelClass, int ORDER, int NVALS = 1>
class FChebTensorialKernel
    : public FAbstractChebKernel< FReal, CellClass, ContainerClass, MatrixKernelClass, ORDER, NVALS>
{
    enum {nRhs = MatrixKernelClass::NRHS,
          nLhs = MatrixKernelClass::NLHS,
          nPot = MatrixKernelClass::NPOT,
          nPv = MatrixKernelClass::NPV};

protected://PB: for OptiDis

    // private types
    typedef FChebTensorialM2LHandler<FReal, ORDER,MatrixKernelClass,MatrixKernelClass::Type> M2LHandlerClass;

    // using from
    typedef FAbstractChebKernel< FReal, CellClass, ContainerClass, MatrixKernelClass, ORDER, NVALS>
    AbstractBaseClass;

    /// Needed for P2P and M2L operators
    const MatrixKernelClass *const MatrixKernel;

    /// Needed for M2L operator
    FSmartPointer<  M2LHandlerClass,FSmartPointerMemory> M2LHandler;

public:
    /**
     * The constructor initializes all constant attributes and it reads the
     * precomputed and compressed M2L operators from a binary file (an
     * runtime_error is thrown if the required file is not valid).
     */
    FChebTensorialKernel(const int inTreeHeight,
                         const FReal inBoxWidth,
                         const FPoint<FReal>& inBoxCenter,
                         const MatrixKernelClass *const inMatrixKernel,
                         const FReal inBoxWidthExtension)
    : FAbstractChebKernel< FReal, CellClass, ContainerClass, MatrixKernelClass, ORDER, NVALS>(inTreeHeight,inBoxWidth,inBoxCenter,inBoxWidthExtension),
      MatrixKernel(inMatrixKernel),
      M2LHandler(new M2LHandlerClass(MatrixKernel,
                                     inTreeHeight,
                                     inBoxWidth,
                                     inBoxWidthExtension))
    { }


    template<class SymbolicData>
    void P2M(typename CellClass::multipole_t* const LeafMultipole,
             const SymbolicData* const LeafSymbData,
             const ContainerClass* const SourceParticles)
    {
        const FPoint<FReal> LeafCellCenter(AbstractBaseClass::getLeafCellCenter(LeafSymbData->getCoordinate()));
        const FReal ExtendedLeafCellWidth(AbstractBaseClass::BoxWidthLeaf
                                          + AbstractBaseClass::BoxWidthExtension);

        for(int idxV = 0 ; idxV < NVALS ; ++idxV){
            // 1) apply Sy
            AbstractBaseClass::Interpolator->applyP2M(
                LeafCellCenter,
                ExtendedLeafCellWidth,
                LeafMultipole->get(idxV*nRhs),
                SourceParticles);
        }
    }


    template<class SymbolicData>
    void M2M(typename CellClass::multipole_t* const FRestrict ParentMultipole,
             const SymbolicData* const ParentSymb,
             const typename CellClass::multipole_t*const FRestrict *const FRestrict ChildMultipoles,
             const SymbolicData* const /*ChildSymbs*/[])
    {
        int TreeLevel = static_cast<int>(ParentSymb->getLevel());
        for(int idxV = 0 ; idxV < NVALS ; ++idxV){
            for(int idxRhs = 0 ; idxRhs < nRhs ; ++idxRhs){
                // update multipole index
                int idxMul = idxV*nRhs + idxRhs;
                FBlas::scal(AbstractBaseClass::nnodes*2, FReal(0.), ParentMultipole->get(idxMul));
                for (unsigned int ChildIndex=0; ChildIndex < 8; ++ChildIndex){
                    if (ChildMultipoles[ChildIndex]){
                        AbstractBaseClass::Interpolator->applyM2M(
                            ChildIndex,
                            ChildMultipoles[ChildIndex]->get(idxMul),
                            ParentMultipole->get(idxMul),
                            TreeLevel/*Cell width extension specific*/);
                    }
                }
            }
        }
    }

    template<class SymbolicData>
    void M2L(typename CellClass::local_expansion_t * const FRestrict TargetExpansion,
             const SymbolicData* const TargetSymb,
             const typename CellClass::multipole_t * const FRestrict SourceMultipoles[],
             const SymbolicData* const FRestrict /*SourceSymbs*/[],
             const int neighborPositions[],
             const int inSize)
    {
        int TreeLevel = static_cast<int>(TargetSymb->getLevel());
        // scale factor for homogeneous case
        const FReal CellWidth(AbstractBaseClass::BoxWidth / FReal(FMath::pow(2, TreeLevel)));
        const FReal ExtendedCellWidth(CellWidth + AbstractBaseClass::BoxWidthExtension);
        const FReal scale(MatrixKernel->getScaleFactor(ExtendedCellWidth));

        for(int idxV = 0 ; idxV < NVALS ; ++idxV){
            for (int idxLhs=0; idxLhs < nLhs; ++idxLhs){
                // update local index
                const int idxLoc = idxV*nLhs + idxLhs;

                FReal *const localExpansion = TargetExpansion->get(idxLoc);

                // update idxRhs
                const int idxRhs = idxLhs % nPv;
                // update multipole index
                const int idxMul = idxV*nRhs + idxRhs;

                // get index in matrix kernel
                const unsigned int d = MatrixKernel->getPosition(idxLhs);

                for(int idxExistingNeigh = 0 ; idxExistingNeigh < inSize ; ++idxExistingNeigh){
                    const int idx = neighborPositions[idxExistingNeigh];
                    M2LHandler->applyC(
                        idx, TreeLevel, scale, d,
                        SourceMultipoles[idxExistingNeigh]->get(idxMul),
                        localExpansion);
                }
            }// NLHS=NPOT*NPV
        }// NVALS

    }

    template<class SymbolicData>
    void L2L(const typename CellClass::local_expansion_t * const FRestrict ParentExpansion,
             const SymbolicData* const ParentSymb,
             typename CellClass::local_expansion_t * FRestrict * const FRestrict ChildExpansions,
             const SymbolicData* const /*ChildSymbs*/[])
    {
        int TreeLevel = static_cast<int>(ParentSymb->getLevel());
        for(int idxV = 0 ; idxV < NVALS ; ++idxV){
            for(int idxLhs = 0 ; idxLhs < nLhs ; ++idxLhs){
                int idxLoc = idxV*nLhs + idxLhs;
                for (unsigned int ChildIndex=0; ChildIndex < 8; ++ChildIndex){
                    if (ChildExpansions[ChildIndex]){
                        AbstractBaseClass::Interpolator->applyL2L(ChildIndex,
                                                                  ParentExpansion->get(idxLoc),
                                                                  ChildExpansions[ChildIndex]->get(idxLoc),
                                                                  TreeLevel/*Cell width extension specific*/);
                    }
                }
            }//NLHS
        }// NVALS
    }


    template<class SymbolicData>
    void L2P(const typename CellClass::local_expansion_t* const LeafExpansion,
             const SymbolicData* const LeafSymbData,
             ContainerClass* const TargetParticles)
    {
        const FPoint<FReal> LeafCellCenter(AbstractBaseClass::getLeafCellCenter(LeafSymbData->getCoordinate()));
        const FReal ExtendedLeafCellWidth(AbstractBaseClass::BoxWidthLeaf
                                          + AbstractBaseClass::BoxWidthExtension);

        for(int idxV = 0 ; idxV < NVALS ; ++idxV){
            AbstractBaseClass::Interpolator->applyL2PTotal(LeafCellCenter, ExtendedLeafCellWidth,
                                                           LeafExpansion->get(idxV*nLhs), TargetParticles);
        }
    }

    void P2P(const FTreeCoordinate& inPosition,
             ContainerClass* const /*FFRestrict */ inTargets,
             const ContainerClass* const /*FFRestrict */ inSources,
             ContainerClass* const inNeighbors[], const int neighborPositions[],
             const int inSize) override {
        if(inTargets == inSources){
            P2POuter(inPosition, inTargets, inNeighbors, neighborPositions, inSize);
            DirectInteractionComputer<FReal, MatrixKernelClass::NCMP, NVALS>::P2PInner(inTargets,MatrixKernel);
        }
        else{
            const ContainerClass* const srcPtr[1] = {inSources};
            DirectInteractionComputer<FReal, MatrixKernelClass::NCMP, NVALS>::P2PRemote(inTargets,srcPtr,1,MatrixKernel);
            DirectInteractionComputer<FReal, MatrixKernelClass::NCMP, NVALS>::P2PRemote(inTargets,inNeighbors,inSize,MatrixKernel);
        }
    }

    void P2POuter(const FTreeCoordinate& /*inLeafPosition*/,
             ContainerClass* const FRestrict inTargets,
             ContainerClass* const inNeighbors[], const int neighborPositions[],
             const int inSize) override {
        int nbNeighborsToCompute = 0;
        while(nbNeighborsToCompute < inSize
              && neighborPositions[nbNeighborsToCompute] < 14){
            nbNeighborsToCompute += 1;
        }
        DirectInteractionComputer<FReal, MatrixKernelClass::NCMP, NVALS>::P2P(inTargets,inNeighbors,nbNeighborsToCompute,MatrixKernel);
    }


    void P2PRemote(const FTreeCoordinate& /*inPosition*/,
                   ContainerClass* const FRestrict inTargets, const ContainerClass* const FRestrict /*inSources*/,
                   const ContainerClass* const inNeighbors[], const int /*neighborPositions*/[],
                   const int inSize) override {
        DirectInteractionComputer<FReal, MatrixKernelClass::NCMP, NVALS>::P2PRemote(inTargets,inNeighbors,inSize,MatrixKernel);
    }

};


#endif //FCHEBTENSORIALKERNELS_HPP

// [--END--]
