// See LICENCE file at project root
// Keep in private GIT

#ifndef FUNIFKERNEL_HPP
#define FUNIFKERNEL_HPP

#include <algorithm>
#include <vector>

#include "Utils/FGlobal.hpp"

#include "Utils/FSmartPointer.hpp"

#include "FAbstractUnifKernel.hpp"
#include "FUnifM2LHandler.hpp"

class FTreeCoordinate;

/**
 * @author Pierre Blanchard (pierre.blanchard@inria.fr)
 * @brief
 *
 * Please read the license
 *
 * This kernels implement the Lagrange interpolation based FMM operators. It
 * implements all interfaces (P2P,P2M,M2M,M2L,L2L,L2P) which are required by
 * the FFmmAlgorithm and FFmmAlgorithmThread.
 *
 * @tparam CellClass Type of cell
 * @tparam ContainerClass Type of container to store particles
 * @tparam MatrixKernelClass Type of matrix kernel function
 * @tparam ORDER Lagrange interpolation order
 */
template < class FReal, class CellClass, class ContainerClass,   class MatrixKernelClass, int ORDER, int NVALS = 1>
class FUnifKernel
  : public FAbstractUnifKernel<FReal, CellClass, ContainerClass, MatrixKernelClass, ORDER, NVALS>
{
protected:
    // private types
    typedef FUnifM2LHandler<FReal, ORDER,MatrixKernelClass::Type> M2LHandlerClass;

    // using from
    typedef FAbstractUnifKernel< FReal, CellClass, ContainerClass, MatrixKernelClass, ORDER, NVALS>
    AbstractBaseClass;

    /// Needed for P2P and M2L operators
    const MatrixKernelClass *const MatrixKernel;

    /// Needed for M2L operator
    const M2LHandlerClass M2LHandler;

    /// Leaf level separation criterion
    const int LeafLevelSeparationCriterion;

public:
    /**
    * The constructor initializes all constant attributes and it reads the
    * precomputed and compressed M2L operators from a binary file (an
    * runtime_error is thrown if the required file is not valid).
    */
    FUnifKernel(const int inTreeHeight,
                const FReal inBoxWidth,
                const FPoint<FReal>& inBoxCenter,
                const MatrixKernelClass *const inMatrixKernel,
                const int inLeafLevelSeparationCriterion = 1)
    : FAbstractUnifKernel< FReal, CellClass, ContainerClass, MatrixKernelClass, ORDER, NVALS>(inTreeHeight,inBoxWidth,inBoxCenter),
      MatrixKernel(inMatrixKernel),
      M2LHandler(MatrixKernel,
                 inTreeHeight,
                 inBoxWidth,
                 inLeafLevelSeparationCriterion),
      LeafLevelSeparationCriterion(inLeafLevelSeparationCriterion)
    { }


    template<class SymbolicData>
    void P2M(typename CellClass::multipole_t* const LeafMultipole,
             const SymbolicData* const LeafSymbData,
             const ContainerClass* const SourceParticles)
    {
        FReal leafBoxWidth = AbstractBaseClass::BoxWidth / FReal(1 <<  LeafSymbData->getLevel());

        const FPoint<FReal> LeafCellCenter(
            AbstractBaseClass::getCellCenter(LeafSymbData->getCoordinate(),
                                             static_cast<int>(LeafSymbData->getLevel())));

        AbstractBaseClass::Interpolator->applyP2M(LeafCellCenter, leafBoxWidth,
                                                  LeafMultipole->get(0), SourceParticles);

        for(int idxRhs = 0 ; idxRhs < NVALS ; ++idxRhs){
            // 2) apply Discrete Fourier Transform
            this->M2LHandler.applyZeroPaddingAndDFT(LeafMultipole->get(idxRhs),
                                                    LeafMultipole->getTransformed(idxRhs));
        }
    }


    template<class SymbolicData>
    void M2M(typename CellClass::multipole_t* const FRestrict ParentMultipole,
             const SymbolicData* const /*ParentSymb*/,
             const typename CellClass::multipole_t*const FRestrict *const FRestrict ChildMultipoles,
             const SymbolicData* const /*ChildSymbs*/[])
    {
        for(int idxRhs = 0 ; idxRhs < NVALS ; ++idxRhs){
            // 1) apply Sy
            //FBlas::scal(AbstractBaseClass::nnodes, FReal(0.), ParentMultipole->get(idxRhs));
            for (unsigned int ChildIndex=0; ChildIndex < 8; ++ChildIndex){
                if (ChildMultipoles[ChildIndex]){
                    AbstractBaseClass::Interpolator->applyM2M(
                        ChildIndex, ChildMultipoles[ChildIndex]->get(idxRhs),
                        ParentMultipole->get(idxRhs));
                }
            }
            // 2) Apply Discete Fourier Transform
            M2LHandler.applyZeroPaddingAndDFT(ParentMultipole->get(idxRhs),
                                              ParentMultipole->getTransformed(idxRhs));
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
        const FReal CellWidth(AbstractBaseClass::BoxWidth / FReal(1 << TargetSymb->getLevel()));
        const FReal scale(MatrixKernel->getScaleFactor(CellWidth));

        for(int idxRhs = 0 ; idxRhs < NVALS ; ++idxRhs){
            stdComplex<FReal> *const TransformedLocalExpansion = TargetExpansion->getTransformed(idxRhs);

            for(int idxExistingNeigh = 0 ; idxExistingNeigh < inSize ; ++idxExistingNeigh){
                const int idxNeigh = neighborPositions[idxExistingNeigh];
                M2LHandler.applyFC(idxNeigh, static_cast<int>(TargetSymb->getLevel()), scale,
                                   SourceMultipoles[idxExistingNeigh]->getTransformed(idxRhs),
                                   TransformedLocalExpansion);
            }
        }
    }


    template<class SymbolicData>
    void L2L(const typename CellClass::local_expansion_t * const FRestrict ParentExpansion,
             const SymbolicData* const /*ParentSymb*/,
             typename CellClass::local_expansion_t * FRestrict * const FRestrict ChildExpansions,
             const SymbolicData* const /*ChildSymbs*/[])
    {
        for(int idxRhs = 0 ; idxRhs < NVALS ; ++idxRhs){
            // 1) Apply Inverse Discete Fourier Transform
            FReal localExp[AbstractBaseClass::nnodes];
            M2LHandler.unapplyZeroPaddingAndDFT(ParentExpansion->getTransformed(idxRhs),
                                                localExp);
            FBlas::add(AbstractBaseClass::nnodes,
                       const_cast<FReal*>(ParentExpansion->get(idxRhs)),localExp);
            // 2) apply Sx
            for (unsigned int ChildIndex=0; ChildIndex < 8; ++ChildIndex){
                if (ChildExpansions[ChildIndex]){
                    AbstractBaseClass::Interpolator->applyL2L(
                        ChildIndex, localExp, ChildExpansions[ChildIndex]->get(idxRhs));
                }
            }
        }
    }


    template<class SymbolicData>
    void L2P(const typename CellClass::local_expansion_t* const LeafExpansion,
             const SymbolicData* const LeafSymbData,
             ContainerClass* const TargetParticles)
    {
        FReal leafBoxWidth = AbstractBaseClass::BoxWidth / FReal(1 <<  LeafSymbData->getLevel());

        const FPoint<FReal> LeafCellCenter(
            AbstractBaseClass::getCellCenter(LeafSymbData->getCoordinate(),
                                             static_cast<int>(LeafSymbData->getLevel())));

        FReal localExp[NVALS*AbstractBaseClass::nnodes];

        for(int idxRhs = 0 ; idxRhs < NVALS ; ++idxRhs){
            // 1)  Apply Inverse Discete Fourier Transform
            this->M2LHandler.unapplyZeroPaddingAndDFT(LeafExpansion->getTransformed(idxRhs),
                                                localExp + idxRhs*AbstractBaseClass::nnodes);
            FBlas::add(AbstractBaseClass::nnodes,const_cast<FReal*>(LeafExpansion->get(idxRhs)),localExp + idxRhs*AbstractBaseClass::nnodes);
        }
        // 2.a) apply Sx
        AbstractBaseClass::Interpolator->applyL2P(LeafCellCenter, leafBoxWidth,
                                      localExp, TargetParticles);
        // 2.b) apply Px (grad Sx)
        AbstractBaseClass::Interpolator->applyL2PGradient(LeafCellCenter, leafBoxWidth,
                                              localExp, TargetParticles);
    }
    ///
    /// \brief P2P  Particle to particle operator
    /// \param inPosition
    /// \param inTargets
    /// \param inSources
    /// \param inNeighbors
    /// \param neighborPositions
    /// \param inSize
    ///
    void P2P(const FTreeCoordinate& inPosition,
             ContainerClass* const /*FFRestrict */ inTargets,
             const ContainerClass* const /*FRestrict*/ inSources, // no restrict inTargets == inSources)
             ContainerClass* const inNeighbors[],
             const int neighborPositions[],
             const int inSize)
        override
    {
        this->P2P(inPosition, inTargets, inSources, inNeighbors, neighborPositions, inSize, true);
    }

    void P2P(const FTreeCoordinate& inPosition,
             ContainerClass* const FRestrict inTargets,
             const ContainerClass* const FRestrict inSources,   // no restrict inTargets == inSources)
             ContainerClass* const inNeighbors[], const int neighborPositions[],
             const int inSize, bool do_inner)
    {
        // Standard FMM separation criterion, i.e. max 27 neighbor clusters per leaf
        if(LeafLevelSeparationCriterion==1) {
            if(inTargets == inSources) {
                P2POuter(inPosition, inTargets, inNeighbors, neighborPositions, inSize);

                if(do_inner) {
                    DirectInteractionComputer<FReal, MatrixKernelClass::NCMP, NVALS>::P2PInner(inTargets,MatrixKernel);
                }
            } else {
                const ContainerClass* const srcPtr[1] = {inSources};
                if(do_inner) {
                    DirectInteractionComputer<FReal, MatrixKernelClass::NCMP, NVALS>::P2PRemote(inTargets,srcPtr,1,MatrixKernel);
                }
                DirectInteractionComputer<FReal, MatrixKernelClass::NCMP, NVALS>::P2PRemote(inTargets,inNeighbors,inSize,MatrixKernel);
            }
        } // Nearfield interactions are only computed within the target leaf
        else if(LeafLevelSeparationCriterion==0){
            DirectInteractionComputer<FReal,MatrixKernelClass::NCMP, NVALS>::P2PRemote(inTargets,inNeighbors,inSize,MatrixKernel);
        }
        // If criterion equals -1 then no P2P need to be performed.
    }

    void P2POuter(const FTreeCoordinate& /*inLeafPosition*/,
             ContainerClass* const FRestrict inTargets,
             ContainerClass* const inNeighbors[], const int neighborPositions[],
             const int inSize) override {
        std::vector<ContainerClass*> neighbours;

        for(int i = 0; i < inSize; ++i) {
            if(neighborPositions[i] < 14) {
                neighbours.push_back(inNeighbors[i]);
            }
        }

        DirectInteractionComputer<FReal, MatrixKernelClass::NCMP, NVALS>::
            P2P(inTargets, neighbours.data(), static_cast<int>(neighbours.size()), MatrixKernel);
    }

    void P2PRemote(const FTreeCoordinate& /*inPosition*/,
                   ContainerClass* const FRestrict inTargets, const ContainerClass* const FRestrict /*inSources*/,
                   const ContainerClass* const inNeighbors[], const int /*neighborPositions*/[],
                   const int inSize) override {
        // Standard FMM separation criterion, i.e. max 27 neighbor clusters per leaf
        if(LeafLevelSeparationCriterion==1)
            DirectInteractionComputer<FReal,MatrixKernelClass::NCMP, NVALS>::P2PRemote(inTargets,inNeighbors,inSize,MatrixKernel);
        // Nearfield interactions are only computed within the target leaf
        if(LeafLevelSeparationCriterion==0)
            DirectInteractionComputer<FReal,MatrixKernelClass::NCMP, NVALS>::P2PRemote(inTargets,inNeighbors,0,MatrixKernel);
        // If criterion equals -1 then no P2P need to be performed.
    }

};


#endif //FUNIFKERNEL_HPP

// [--END--]
