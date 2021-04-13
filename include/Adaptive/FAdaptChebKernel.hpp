#ifndef FADAPTCHEBKERNEL_HPP
#define FADAPTCHEBKERNEL_HPP

#include "Kernels/Chebyshev/FChebSymKernel_i.hpp"
#include "InastempCompileConfig.h"

#ifdef _OPENMP
#include <omp.h>
#endif

template<
    class FReal,
    class CellClass,
    class ContainerClass,
    class MatrixKernelClass,
    int ORDER,
    int NVALS = 1
    >
class FAdaptChebKernel : public FChebSymKernel<FReal, CellClass, ContainerClass, MatrixKernelClass, ORDER, NVALS>{
    using FBase =  FChebSymKernel_i<FReal, CellClass, ContainerClass, MatrixKernelClass, ORDER, NVALS>;

public:
    using FBase::FBase;

    using FBase::M2M;
    using FBase::M2L;
    using FBase::L2L;
    using FBase::P2P;

    void setup() {

        #ifdef _OPENMP
        if(omp_get_thread_num() == 0)
        #endif
        {
            std::cout << "Symetric Chebyshev adaptive kernel" << '\n';
        }
    }

    template<class SymbolicData>
    void P2M(typename CellClass::multipole_t* const LeafCell,
             const SymbolicData * const LeafSymbData,
             const ContainerClass* const SourceParticles)
    {
        FReal leafBoxWidth = FBase::BoxWidth / FReal(1 <<  LeafSymbData->getLevel());
        const FPoint<FReal> LeafCellCenter =
            FBase::getCellCenter(LeafSymbData->getCoordinate(), static_cast<int>(LeafSymbData->getLevel()));
        FBase::Interpolator->applyP2M(LeafCellCenter, leafBoxWidth,
                                      LeafCell->get(0), SourceParticles);
    }

    template<class SymbolicData>
    void L2P(const typename CellClass::local_expansion_t* const LeafCell,
             const SymbolicData * const LeafSymbData,
             ContainerClass* const TargetParticles)
    {
        FReal leafBoxWidth = FBase::BoxWidth / FReal(1 <<  LeafSymbData->getLevel());
        const FPoint<FReal> LeafCellCenter =
            FBase::getCellCenter(LeafSymbData->getCoordinate(), static_cast<int>(LeafSymbData->getLevel()));
        // apply Sx and Px (grad Sx)
        FBase::Interpolator->applyL2PTotal(LeafCellCenter, leafBoxWidth,
                                           LeafCell->get(0), TargetParticles);
    }

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
        const FPoint<FReal> localCellCenter(
            FBase::getCellCenter(symb->getCoordinate(),static_cast<int>(symb->getLevel())));

        // interpolation points of target (X) cell
        FPoint<FReal> X[FBase::nnodes];
        FChebTensor<FReal,ORDER>::setRoots(localCellCenter, localCellWidth, X);

        // Particles attributes
        const FReal * const posX = particles->getPositions()[0];
        const FReal * const posY = particles->getPositions()[1];
        const FReal * const posZ = particles->getPositions()[2];
        const FReal * const physicalValues = particles->getPhysicalValues();

        // apply P2L
        for(int idxRhs = 0 ; idxRhs < NVALS ; ++idxRhs){
            for (unsigned int m = 0; m<FBase::nnodes; ++m) {
                ComputeClass XX = ComputeClass(X[m].getX());
                ComputeClass XY = ComputeClass(X[m].getY());
                ComputeClass XZ = ComputeClass(X[m].getZ());

                // Compute using vectorization for all but the last array elements
                ComputeClass tmpLocalExp = ComputeClass::GetZero();
                for (std::size_t idxPart = 0 ; idxPart < particles->getNbParticles() ; idxPart += FRealCount)
                {
                    tmpLocalExp +=
                        FBase::MatrixKernel->evaluate(																				// gotta be updated due to FInterpMatrixKernel change.
                            XX, XY, XZ,
                                ComputeClass(&posX[idxPart]), ComputeClass(&posY[idxPart]), ComputeClass(&posZ[idxPart]))
                        * physicalValues[idxPart];
                }

                local->get(idxRhs)[m] += (tmpLocalExp.horizontalSum());
            }
        }// NVALS
    }

    template<class SymbolicData>
    void M2P(const typename CellClass::multipole_t* const pole,
             const SymbolicData * const symb,
             ContainerClass* const particles,
             const SymbolicData * const /*source_symb*/)
    {
			std::cout << "\n"<< "\n"<<" ----------------------------------------------------------------------------------------------"  << std::endl; 		
			std::cout <<" --------------------------------------------   THIS DOES NOT PRINT  (DNP) --------------------------------------------------"  << std::endl; 			
			std::cout <<" ----------------------------------------------------------------------------------------------"   << "\n"<< "\n" << std::endl; 					
		
        using ComputeClass = InaVecBestType<FReal>;
        constexpr int FRealCount = ComputeClass::VecLength;

        // Source cell: pole
        const FReal poleCellWidth(FBase::BoxWidth / FReal(1 << symb->getLevel()));
        const FPoint<FReal> poleCellCenter(
            FBase::getCellCenter(symb->getCoordinate(),static_cast<int>(symb->getLevel())));

        // interpolation points of source (Y) cell
        FPoint<FReal> Y[FBase::nnodes];
        FChebTensor<FReal,ORDER>::setRoots(poleCellCenter, poleCellWidth, Y);

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
                    FBase::MatrixKernel->evaluateBlockAndDerivative(																	// gotta be updated due to FInterpMatrixKernel change.
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

#endif /* FADAPTCHEBKERNEL_HPP */
