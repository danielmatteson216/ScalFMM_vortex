// See LICENCE file at project root
#ifndef FCHEBKERNEL_HPP
#define FCHEBKERNEL_HPP

#include "Utils/FGlobal.hpp"

#include "Utils/FSmartPointer.hpp"

#include "Kernels/Chebyshev/FAbstractChebKernel.hpp"

#include "Kernels/Chebyshev/FChebM2LHandler.hpp"

class FTreeCoordinate;

/**
 * @author Matthias Messner(matthias.messner@inria.fr)
 * @class FChebKernel
 * @brief  Chebyshev interpolation based FMM operators for general non oscillatory kernels..
 *
 * This class implements the Chebyshev interpolation based FMM operators. It
 * implements all interfaces (P2P, P2M, M2M, M2L, L2L, L2P) which are required by
 * the FFmmAlgorithm, FFmmAlgorithmThread ...
 *
 * @tparam CellClass Type of cell
 * @tparam ContainerClass Type of container to store particles
 * @tparam MatrixKernelClass Type of matrix kernel function
 * @tparam ORDER Chebyshev interpolation order
 *
 * The ORDER sets the accuracy of the Chebyshev FMM while the EPSILON parameter introduces extra error but optimize the M2L step.
 *  In fact, in the Chebyshev FMM we perform compression on the M2L operators using various low rank approximation techniques
 *  (see https://arxiv.org/abs/1210.7292 for further details). Therefore we use a second accuracy criterion, namely EPSILON,
 *  in order to set the accuracy of these methods. For most kernels that we tested and in particular for 1/r, setting EPSILON=10^-ORDER d
 *  oes not introduce extra error in the FMM and captures the rank efficiently. If you think that for your kernel you need a better
 *  approximation of the M2L operators, then you can try to set EPSILON to 10 ^- (ORDER+{1,2,...}).
 */
template < class FReal, class CellClass, class ContainerClass,   class MatrixKernelClass, int ORDER, int NVALS = 1>
class FChebKernel
    : public FAbstractChebKernel< FReal, CellClass, ContainerClass, MatrixKernelClass, ORDER, NVALS>
{
    // private types
    typedef FChebM2LHandler<FReal, ORDER,MatrixKernelClass> M2LHandlerClass;

    // using from
    typedef FAbstractChebKernel<FReal, CellClass, ContainerClass, MatrixKernelClass, ORDER, NVALS>
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
     *
     * @param[in] epsilon  The compression parameter for M2L operator.
     *
     * The M2L optimized Chebyshev FMM implemented in ScalFMM are kernel dependent, but keeping EPSILON=10^-ORDER is usually fine.
     *  On the other hand you can short-circuit this feature by setting EPSILON to the machine accuracy,
     *  but this will significantly slow down the computations.
     *
     */
    FChebKernel(const int inTreeHeight,  const FReal inBoxWidth, const FPoint<FReal>& inBoxCenter, const MatrixKernelClass *const inMatrixKernel,
                const FReal Epsilon)
    : FAbstractChebKernel< FReal, CellClass, ContainerClass, MatrixKernelClass, ORDER, NVALS>(inTreeHeight,
                                                                                       inBoxWidth,
                                                                                       inBoxCenter),
      MatrixKernel(inMatrixKernel),
      M2LHandler(new M2LHandlerClass(MatrixKernel,Epsilon))
    {
        // read precomputed compressed m2l operators from binary file
        //M2LHandler->ReadFromBinaryFileAndSet();
        M2LHandler->ComputeAndCompressAndSet();
    }


    /**
     * The constructor initializes all constant attributes and it reads the
     * precomputed and compressed M2L operators from a binary file (an
     * runtime_error is thrown if the required file is not valid).
     * Same as \see above constructor  but the epsilon is automatically set to EPSILON=10^-ORDER
     */
    FChebKernel(const int inTreeHeight, const FReal inBoxWidth, const FPoint<FReal>& inBoxCenter, const MatrixKernelClass *const inMatrixKernel)
        :   FChebKernel(inTreeHeight, inBoxWidth,inBoxCenter,inMatrixKernel,FMath::pow(10.0,static_cast<FReal>(-ORDER)))
    {

    }


    template<class SymbolicData>
    void P2M(typename CellClass::multipole_t* const LeafMultipole,
             const SymbolicData* const LeafSymbData,
             const ContainerClass* const SourceParticles)
    {
        int leafLevel = static_cast<int>(LeafSymbData->getLevel());
        FReal leafBoxWidth = AbstractBaseClass::BoxWidth / FReal(1 << leafLevel);

        const FPoint<FReal> LeafCellCenter(
            AbstractBaseClass::getCellCenter(LeafSymbData->getCoordinate(),
                                             leafLevel));

        // 1) apply Sy
        AbstractBaseClass::Interpolator->applyP2M(LeafCellCenter, leafBoxWidth,
                                                  LeafMultipole->get(0), SourceParticles);

        for(int idxRhs = 0 ; idxRhs < NVALS ; ++idxRhs){
            // 2) apply B
            M2LHandler->applyB(LeafMultipole->get(idxRhs),
                               LeafMultipole->get(idxRhs) + AbstractBaseClass::nnodes);
        }
    }


    template<class SymbolicData>
    void M2M(typename CellClass::multipole_t * const FRestrict ParentMultipole,
             const SymbolicData* const /*ParentSymb*/,
             const typename CellClass::multipole_t * const FRestrict * const FRestrict ChildMultipoles,
             const SymbolicData* const /*ChildSymbs*/[])
    {
        for(int idxRhs = 0 ; idxRhs < NVALS ; ++idxRhs){
            // 1) apply Sy
            for (unsigned int ChildIndex=0; ChildIndex < 8; ++ChildIndex){
                if (ChildMultipoles[ChildIndex]){
                    AbstractBaseClass::Interpolator->applyM2M(
                        ChildIndex, ChildMultipoles[ChildIndex]->get(idxRhs),
                        ParentMultipole->get(idxRhs));
                }
            }
            // 2) apply B
            M2LHandler->applyB(ParentMultipole->get(idxRhs),
                               ParentMultipole->get(idxRhs) + AbstractBaseClass::nnodes);
        }
    }


//  void M2L(CellClass* const FRestrict TargetCell,
//                   const CellClass* SourceCells[],
//                   const int NumSourceCells,
//                   const int TreeLevel) const
//  {
//      const FReal CellWidth(BoxWidth / FReal(FMath::pow(2, TreeLevel)));
//      const FTreeCoordinate& cx = TargetCell->getCoordinate();
//      for (int idx=0; idx<NumSourceCells; ++idx) {
//          const FTreeCoordinate& cy = SourceCells[idx]->getCoordinate();
//          const int transfer[3] = {cy.getX()-cx.getX(),
//                                                           cy.getY()-cx.getY(),
//                                                           cy.getZ()-cx.getZ()};
//          M2LHandler->applyC(transfer, CellWidth,
//                                              SourceCells[idx]->get() + AbstractBaseClass::nnodes,
//                                              TargetCell->get() + AbstractBaseClass::nnodes);
//      }
//  }

    template<class SymbolicData>
    void M2L(typename CellClass::local_expansion_t * const FRestrict TargetExpansion,
             const SymbolicData* const TargetSymb,
             const typename CellClass::multipole_t * const FRestrict SourceMultipoles[],
             const SymbolicData* const FRestrict /*SourceSymbs*/[],
             const int neighborPositions[],
             const int inSize)
    {
        for(int idxRhs = 0 ; idxRhs < NVALS ; ++idxRhs){
            FReal *const CompressedLocalExpansion = TargetExpansion->get(idxRhs) + AbstractBaseClass::nnodes;
            const FReal CellWidth(AbstractBaseClass::BoxWidth / FReal(1 << TargetSymb->getLevel()));
            for(int idxExistingNeigh = 0 ; idxExistingNeigh < inSize ; ++idxExistingNeigh){
                const int idx = neighborPositions[idxExistingNeigh];
                M2LHandler->applyC(idx, CellWidth,
                                   SourceMultipoles[idxExistingNeigh]->get(idxRhs)
                                   + AbstractBaseClass::nnodes,
                                   CompressedLocalExpansion);
            }
        }
    }

//  void M2L(CellClass* const FRestrict TargetCell, const CellClass* SourceCells[],
//    const int neighborPositions[], const int inSize, const int TreeLevel)  override {
//      const unsigned int rank = M2LHandler.getRank();
//      FBlas::scal(343*rank, FReal(0.), MultipoleExpansion);
//      const FReal CellWidth(BoxWidth / FReal(FMath::pow(2, TreeLevel)));
//      for(int idxExistingNeigh = 0 ; idxExistingNeigh < inSize ; ++idxExistingNeigh){
//          const int idx = neighborPositions[idxExistingNeigh];
//              FBlas::copy(rank, const_cast<FReal *const>(SourceCells[idxExistingNeigh]->get())+AbstractBaseClass::nnodes,
//                                      MultipoleExpansion+idx*rank);
//
//      M2LHandler->applyC(CellWidth, MultipoleExpansion, TargetCell->get() + AbstractBaseClass::nnodes);
//  }


    template<class SymbolicData>
    void L2L(const typename CellClass::local_expansion_t * const FRestrict ParentExpansion,
             const SymbolicData* const /*ParentSymb*/,
             typename CellClass::local_expansion_t * FRestrict *const FRestrict ChildExpansions,
             const SymbolicData* const /*ChildSymbs*/[])
    {
        for(int idxRhs = 0 ; idxRhs < NVALS ; ++idxRhs){
            // 1) apply U
            M2LHandler->applyU(ParentExpansion->get(idxRhs) + AbstractBaseClass::nnodes,
                               const_cast<typename CellClass::local_expansion_t*>(ParentExpansion)->get(idxRhs));
            // 2) apply Sx
            for (unsigned int ChildIndex=0; ChildIndex < 8; ++ChildIndex){
                if (ChildExpansions[ChildIndex]){
                    AbstractBaseClass::Interpolator->applyL2L(
                        ChildIndex,
                        ParentExpansion->get(idxRhs),
                        ChildExpansions[ChildIndex]->get(idxRhs));
                }
            }
        }
    }

    template<class SymbolicData>
    void L2P(const typename CellClass::local_expansion_t* const LeafExpansion,
             const SymbolicData* const LeafSymb,
             ContainerClass* const TargetParticles)
    {
        int leafLevel = static_cast<int>(LeafSymb->getLevel());
        FReal leafBoxWidth = AbstractBaseClass::BoxWidth / FReal(1 << leafLevel);

        const FPoint<FReal> LeafCellCenter(
            AbstractBaseClass::getCellCenter(LeafSymb->getCoordinate(),
                                             leafLevel));

        for(int idxRhs = 0 ; idxRhs < NVALS ; ++idxRhs){
            // 1) apply U
            M2LHandler->applyU(LeafExpansion->get(idxRhs) + AbstractBaseClass::nnodes,
                               const_cast<typename CellClass::local_expansion_t*>(LeafExpansion)->get(idxRhs));
        }

        //// 2.a) apply Sx
        //AbstractBaseClass::Interpolator->applyL2P(LeafCellCenter,
        //                                          AbstractBaseClass::BoxWidthLeaf,
        //                                          LeafCell->get(0),
        //                                          TargetParticles);
        //// 2.b) apply Px (grad Sx)
        //AbstractBaseClass::Interpolator->applyL2PGradient(LeafCellCenter,
        //                                                  AbstractBaseClass::BoxWidthLeaf,
        //                                                  LeafCell->get(0),
        //                                                  TargetParticles);

        // 2.c) apply Sx and Px (grad Sx)
        AbstractBaseClass::Interpolator->applyL2PTotal(LeafCellCenter, leafBoxWidth,
                                                       LeafExpansion->get(0), TargetParticles);

    }

    void P2P(const FTreeCoordinate& inPosition,
             ContainerClass* const /*FFRestrict */ inTargets,
             const ContainerClass* const /*FFRestrict */ inSources,
             ContainerClass* const inNeighbors[], const int neighborPositions[],
             const int inSize)
        override
    {
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


#endif //FCHEBKERNELS_HPP

// [--END--]
