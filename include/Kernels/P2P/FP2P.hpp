// See LICENCE file at project root
#ifndef FP2P_HPP
#define FP2P_HPP

#include "Utils/FPoint.hpp"


namespace FP2P {

/**
   * @brief MutualParticles (generic version)
   * P2P mutual interaction,
   * this function computes the interaction for 2 particles.
   *
   * Formulas are:
   * \f[ F = - q_1 * q_2 * grad(K_{12}) \f]
   * \f[ P_1 = q_2 * K_{12} \f]
   * \f[ P_2 = q_1 * K_{12} \f]
   * In details for \f$\displaystyle K(x,y)=\frac{1}{|x-y|}=\frac{1}{r}\f$ :
   * \f[\displaystyle F(x) = \frac{ \Delta_x * q_1 * q_2 }{ r^2 } \f]
   * \f[\displaystyle P_1 = \frac{ q_2 }{ r } \f]
   * \f[\displaystyle P_2 = \frac{ q_1 }{ r } \f]
   *
   * @param sourceX
   * @param sourceY
   * @param sourceZ
   * @param sourcePhysicalValue
   * @param targetX
   * @param targetY
   * @param targetZ
   * @param targetPhysicalValue
   * @param targetForceX
   * @param targetForceY
   * @param targetForceZ
   * @param targetPotential
   * @param MatrixKernel pointer to an interaction kernel evaluator
   */
template <class FReal, typename MatrixKernelClass>
inline void MutualParticles(const FReal targetX,const FReal targetY,const FReal targetZ, const FReal targetPhysicalValue,
                            FReal* targetForceX, FReal* targetForceY, FReal* targetForceZ, FReal* targetPotential,
                            const FReal sourceX,const FReal sourceY,const FReal sourceZ, const FReal sourcePhysicalValue,
                            FReal* sourceForceX, FReal* sourceForceY, FReal* sourceForceZ, FReal* sourcePotential,
                            const MatrixKernelClass *const MatrixKernel){

    // Compute kernel of interaction...
    const FPoint<FReal> sourcePoint(sourceX,sourceY,sourceZ);
    const FPoint<FReal> targetPoint(targetX,targetY,targetZ);
    FReal Kxy[1];
    FReal dKxy[3];
    MatrixKernel->evaluateBlockAndDerivative(targetPoint,sourcePoint,Kxy,dKxy);
    const FReal mutual_coeff = MatrixKernel->getMutualCoefficient(); // 1 if symmetric; -1 if antisymmetric

    FReal coef = (targetPhysicalValue * sourcePhysicalValue);

    (*targetForceX) += dKxy[0] * coef;
    (*targetForceY) += dKxy[1] * coef;
    (*targetForceZ) += dKxy[2] * coef;
    (*targetPotential) += ( Kxy[0] * sourcePhysicalValue );

    (*sourceForceX) -= dKxy[0] * coef;
    (*sourceForceY) -= dKxy[1] * coef;
    (*sourceForceZ) -= dKxy[2] * coef;
    (*sourcePotential) += ( mutual_coeff * Kxy[0] * targetPhysicalValue );
}

/**
   * @brief NonMutualParticles (generic version)
   * P2P mutual interaction,
   * this function computes the interaction for 2 particles.
   *
   * Formulas are:
   * \f[
   * F = - q_1 * q_2 * grad K{12}
   * P_1 = q_2 * K{12} ; P_2 = q_1 * K_{12}
   * \f]
   * In details for \f$K(x,y)=1/|x-y|=1/r\f$ :
   * \f$ P_1 = \frac{ q_2 }{ r } \f$
   * \f$ P_2 = \frac{ q_1 }{ r } \f$
   * \f$ F(x) = \frac{ \Delta_x * q_1 * q_2 }{ r^2 } \f$
   */
template <class FReal, typename MatrixKernelClass>
inline void NonMutualParticles(const FReal targetX,const FReal targetY,const FReal targetZ, const FReal targetPhysicalValue,
                               FReal* targetForceX, FReal* targetForceY, FReal* targetForceZ, FReal* targetPotential,
                               const FReal sourceX,const FReal sourceY,const FReal sourceZ, const FReal sourcePhysicalValue,
                               const MatrixKernelClass *const MatrixKernel){

    // Compute kernel of interaction...
    const FPoint<FReal> sourcePoint(sourceX,sourceY,sourceZ);
    const FPoint<FReal> targetPoint(targetX,targetY,targetZ);
    FReal Kxy[1];
    FReal dKxy[3];
    MatrixKernel->evaluateBlockAndDerivative(targetPoint,sourcePoint,Kxy,dKxy);

    FReal coef = (targetPhysicalValue * sourcePhysicalValue);

    (*targetForceX) += dKxy[0] * coef;
    (*targetForceY) += dKxy[1] * coef;
    (*targetForceZ) += dKxy[2] * coef;
    (*targetPotential) += ( Kxy[0] * sourcePhysicalValue );
}




template <class FReal, class ContainerClass, class MatrixKernelClass, class ComputeClass, int NbFRealInComputeClass>
static void GenericFullMutual(ContainerClass* const FRestrict inTargets,
                              ContainerClass* const inNeighbors[],
                              const int limiteNeighbors,
                              const MatrixKernelClass *const MatrixKernel){

    const FSize nbParticlesTargets = inTargets->getNbParticles();
    const FReal*const targetsPhysicalValues = inTargets->getPhysicalValues();
    const FReal*const targetsX = inTargets->getPositions()[0];
    const FReal*const targetsY = inTargets->getPositions()[1];
    const FReal*const targetsZ = inTargets->getPositions()[2];
    FReal*const targetsForcesX = inTargets->getForcesX();
    FReal*const targetsForcesY = inTargets->getForcesY();
    FReal*const targetsForcesZ = inTargets->getForcesZ();
    FReal*const targetsPotentials = inTargets->getPotentials();

    for(FSize idxNeighbors = 0 ; idxNeighbors < limiteNeighbors ; ++idxNeighbors){
        if( inNeighbors[idxNeighbors] ){
            const FSize nbParticlesSources = inNeighbors[idxNeighbors]->getNbParticles();
            const FReal*const sourcesPhysicalValues = inNeighbors[idxNeighbors]->getPhysicalValues();
            const FReal*const sourcesX = inNeighbors[idxNeighbors]->getPositions()[0];
            const FReal*const sourcesY = inNeighbors[idxNeighbors]->getPositions()[1];
            const FReal*const sourcesZ = inNeighbors[idxNeighbors]->getPositions()[2];
            FReal*const sourcesForcesX = inNeighbors[idxNeighbors]->getForcesX();
            FReal*const sourcesForcesY = inNeighbors[idxNeighbors]->getForcesY();
            FReal*const sourcesForcesZ = inNeighbors[idxNeighbors]->getForcesZ();
            FReal*const sourcesPotentials = inNeighbors[idxNeighbors]->getPotentials();

            for(FSize idxTarget = 0 ; idxTarget < nbParticlesTargets ; ++idxTarget){
                FSize idxSource = 0;
                {
                    const FSize nbVectorizedInteractions = (nbParticlesSources/NbFRealInComputeClass)*NbFRealInComputeClass;
                    const ComputeClass tx = ComputeClass(targetsX[idxTarget]);
                    const ComputeClass ty = ComputeClass(targetsY[idxTarget]);
                    const ComputeClass tz = ComputeClass(targetsZ[idxTarget]);
                    const ComputeClass tv = ComputeClass(targetsPhysicalValues[idxTarget]);
                    ComputeClass  tfx = ComputeClass::GetZero();
                    ComputeClass  tfy = ComputeClass::GetZero();
                    ComputeClass  tfz = ComputeClass::GetZero();
                    ComputeClass  tpo = ComputeClass::GetZero();

                    for( ; idxSource < nbVectorizedInteractions ; idxSource += NbFRealInComputeClass){
                        ComputeClass Kxy[1];
                        ComputeClass dKxy[3];
                        MatrixKernel->evaluateBlockAndDerivative(tx,ty,tz,
                                                                 ComputeClass(&sourcesX[idxSource]),
                                                                 ComputeClass(&sourcesY[idxSource]),
                                                                 ComputeClass(&sourcesZ[idxSource]),
                                                                 Kxy,dKxy);
                        const ComputeClass mutual_coeff = ComputeClass(MatrixKernel->getMutualCoefficient());; // 1 if symmetric; -1 if antisymmetric

                        const ComputeClass coef = (tv * ComputeClass(&sourcesPhysicalValues[idxSource]));

                        dKxy[0] *= coef;
                        dKxy[1] *= coef;
                        dKxy[2] *= coef;

                        tfx += dKxy[0];
                        tfy += dKxy[1];
                        tfz += dKxy[2];
                        tpo += Kxy[0] * ComputeClass(&sourcesPhysicalValues[idxSource]);

                        (ComputeClass(&sourcesForcesX[idxSource]) - dKxy[0]).storeInArray(&sourcesForcesX[idxSource]);
                        (ComputeClass(&sourcesForcesY[idxSource]) - dKxy[1]).storeInArray(&sourcesForcesY[idxSource]);
                        (ComputeClass(&sourcesForcesZ[idxSource]) - dKxy[2]).storeInArray(&sourcesForcesZ[idxSource]);
                        (ComputeClass(&sourcesPotentials[idxSource]) + mutual_coeff * Kxy[0] * tv).storeInArray(&sourcesPotentials[idxSource]);
                    }

                    targetsForcesX[idxTarget] += tfx.horizontalSum();
                    targetsForcesY[idxTarget] += tfy.horizontalSum();
                    targetsForcesZ[idxTarget] += tfz.horizontalSum();
                    targetsPotentials[idxTarget] += tpo.horizontalSum();
                }
                {
                    const FReal tx = FReal(targetsX[idxTarget]);
                    const FReal ty = FReal(targetsY[idxTarget]);
                    const FReal tz = FReal(targetsZ[idxTarget]);
                    const FReal tv = FReal(targetsPhysicalValues[idxTarget]);
                    FReal  tfx = FReal(0.);
                    FReal  tfy = FReal(0.);
                    FReal  tfz = FReal(0.);
                    FReal  tpo = FReal(0.);

                    for( ; idxSource < nbParticlesSources ; idxSource += 1){
                        FReal Kxy[1];
                        FReal dKxy[3];
                        MatrixKernel->evaluateBlockAndDerivative(tx,ty,tz,
                                                                 FReal(sourcesX[idxSource]),
                                                                 FReal(sourcesY[idxSource]),
                                                                 FReal(sourcesZ[idxSource]),
                                                                 Kxy,dKxy);
                        const FReal mutual_coeff = FReal(MatrixKernel->getMutualCoefficient());; // 1 if symmetric; -1 if antisymmetric

                        const FReal coef = (tv * FReal(sourcesPhysicalValues[idxSource]));

                        dKxy[0] *= coef;
                        dKxy[1] *= coef;
                        dKxy[2] *= coef;

                        tfx += dKxy[0];
                        tfy += dKxy[1];
                        tfz += dKxy[2];
                        tpo += Kxy[0] * FReal(sourcesPhysicalValues[idxSource]);

                        sourcesForcesX[idxSource] -= dKxy[0];
                        sourcesForcesY[idxSource] -= dKxy[1];
                        sourcesForcesZ[idxSource] -= dKxy[2];
                        sourcesPotentials[idxSource] += mutual_coeff * Kxy[0] * tv;
                    }

                    targetsForcesX[idxTarget] += tfx;
                    targetsForcesY[idxTarget] += tfy;
                    targetsForcesZ[idxTarget] += tfz;
                    targetsPotentials[idxTarget] += tpo;
                }
            }
        }
    }
}

template <class FReal, class ContainerClass, class MatrixKernelClass, class ComputeClass, int NbFRealInComputeClass>
static void GenericInner(ContainerClass* const FRestrict inTargets, const MatrixKernelClass *const MatrixKernel){

    const FSize nbParticlesTargets = inTargets->getNbParticles();
    const FReal*const targetsPhysicalValues = inTargets->getPhysicalValues();
    const FReal*const targetsX = inTargets->getPositions()[0];
    const FReal*const targetsY = inTargets->getPositions()[1];
    const FReal*const targetsZ = inTargets->getPositions()[2];
    FReal*const targetsForcesX = inTargets->getForcesX();
    FReal*const targetsForcesY = inTargets->getForcesY();
    FReal*const targetsForcesZ = inTargets->getForcesZ();
    FReal*const targetsPotentials = inTargets->getPotentials();

    {//In this part, we compute (vectorially) the interaction
        //within the target leaf.

        const FSize nbParticlesSources = nbParticlesTargets;
        const FReal*const sourcesPhysicalValues = targetsPhysicalValues;
        const FReal*const sourcesX = targetsX;
        const FReal*const sourcesY = targetsY;
        const FReal*const sourcesZ = targetsZ;
        FReal*const sourcesForcesX = targetsForcesX;
        FReal*const sourcesForcesY = targetsForcesY;
        FReal*const sourcesForcesZ = targetsForcesZ;
        FReal*const sourcesPotentials = targetsPotentials;

        for(FSize idxTarget = 0 ; idxTarget < nbParticlesTargets ; ++idxTarget){
            FSize idxSource = idxTarget+1;
            {
                const FSize nbVectorizedInteractions = ((nbParticlesSources-idxSource)/NbFRealInComputeClass)*NbFRealInComputeClass + idxSource;
                const ComputeClass tx = ComputeClass(targetsX[idxTarget]);
                const ComputeClass ty = ComputeClass(targetsY[idxTarget]);
                const ComputeClass tz = ComputeClass(targetsZ[idxTarget]);
                const ComputeClass tv = ComputeClass(targetsPhysicalValues[idxTarget]);
                ComputeClass  tfx = ComputeClass::GetZero();
                ComputeClass  tfy = ComputeClass::GetZero();
                ComputeClass  tfz = ComputeClass::GetZero();
                ComputeClass  tpo = ComputeClass::GetZero();

                for( ; idxSource < nbVectorizedInteractions ; idxSource += NbFRealInComputeClass){
                    ComputeClass Kxy[1];
                    ComputeClass dKxy[3];
                    MatrixKernel->evaluateBlockAndDerivative(tx,ty,tz,
                                                             ComputeClass(&sourcesX[idxSource]),
                                                             ComputeClass(&sourcesY[idxSource]),
                                                             ComputeClass(&sourcesZ[idxSource]),
                                                             Kxy,dKxy);
                    const ComputeClass mutual_coeff = ComputeClass(MatrixKernel->getMutualCoefficient()); // 1 if symmetric; -1 if antisymmetric

                    const ComputeClass coef = (tv * ComputeClass(&sourcesPhysicalValues[idxSource]));

                    dKxy[0] *= coef;
                    dKxy[1] *= coef;
                    dKxy[2] *= coef;

                    tfx += dKxy[0];
                    tfy += dKxy[1];
                    tfz += dKxy[2];
                    tpo += Kxy[0]*ComputeClass(&sourcesPhysicalValues[idxSource]);

                    (ComputeClass(&sourcesForcesX[idxSource]) - dKxy[0]).storeInArray(&sourcesForcesX[idxSource]);
                    (ComputeClass(&sourcesForcesY[idxSource]) - dKxy[1]).storeInArray(&sourcesForcesY[idxSource]);
                    (ComputeClass(&sourcesForcesZ[idxSource]) - dKxy[2]).storeInArray(&sourcesForcesZ[idxSource]);
                    (ComputeClass(&sourcesPotentials[idxSource]) + mutual_coeff * Kxy[0] * tv).storeInArray(&sourcesPotentials[idxSource]);
                }

                targetsForcesX[idxTarget] += tfx.horizontalSum();
                targetsForcesY[idxTarget] += tfy.horizontalSum();
                targetsForcesZ[idxTarget] += tfz.horizontalSum();
                targetsPotentials[idxTarget] += tpo.horizontalSum();
            }
            {
                const FReal tx = FReal(targetsX[idxTarget]);
                const FReal ty = FReal(targetsY[idxTarget]);
                const FReal tz = FReal(targetsZ[idxTarget]);
                const FReal tv = FReal(targetsPhysicalValues[idxTarget]);
                FReal  tfx = FReal(0.);
                FReal  tfy = FReal(0.);
                FReal  tfz = FReal(0.);
                FReal  tpo = FReal(0.);

                for( ; idxSource < nbParticlesSources ; idxSource += 1){
                    FReal Kxy[1];
                    FReal dKxy[3];
                    MatrixKernel->evaluateBlockAndDerivative(tx,ty,tz,
                                                             FReal(sourcesX[idxSource]),
                                                             FReal(sourcesY[idxSource]),
                                                             FReal(sourcesZ[idxSource]),
                                                             Kxy,dKxy);
                    const FReal mutual_coeff = FReal(MatrixKernel->getMutualCoefficient()); // 1 if symmetric; -1 if antisymmetric

                    const FReal coef = (tv * FReal(sourcesPhysicalValues[idxSource]));

                    dKxy[0] *= coef;
                    dKxy[1] *= coef;
                    dKxy[2] *= coef;

                    tfx += dKxy[0];
                    tfy += dKxy[1];
                    tfz += dKxy[2];
                    tpo += Kxy[0]*sourcesPhysicalValues[idxSource];

                    sourcesForcesX[idxSource] -= dKxy[0];
                    sourcesForcesY[idxSource] -= dKxy[1];
                    sourcesForcesZ[idxSource] -= dKxy[2];
                    sourcesPotentials[idxSource] += mutual_coeff * Kxy[0] * tv;
                }

                targetsForcesX[idxTarget] += tfx;
                targetsForcesY[idxTarget] += tfy;
                targetsForcesZ[idxTarget] += tfz;
                targetsPotentials[idxTarget] += tpo;
            }
        }
    }
}

template <class FReal, class ContainerClass, class MatrixKernelClass, class ComputeClass, int NbFRealInComputeClass>
static void GenericFullRemote(ContainerClass* const FRestrict inTargets, const ContainerClass* const inNeighbors[],
                              const int limiteNeighbors, const MatrixKernelClass *const MatrixKernel){

    const FSize nbParticlesTargets = inTargets->getNbParticles();
    const FReal*const targetsPhysicalValues = inTargets->getPhysicalValues();
    const FReal*const targetsX = inTargets->getPositions()[0];
    const FReal*const targetsY = inTargets->getPositions()[1];
    const FReal*const targetsZ = inTargets->getPositions()[2];
    FReal*const targetsForcesX = inTargets->getForcesX();
    FReal*const targetsForcesY = inTargets->getForcesY();
    FReal*const targetsForcesZ = inTargets->getForcesZ();
    FReal*const targetsPotentials = inTargets->getPotentials();

    for(FSize idxNeighbors = 0 ; idxNeighbors < limiteNeighbors ; ++idxNeighbors){
        if( inNeighbors[idxNeighbors] ){
            const FSize nbParticlesSources = inNeighbors[idxNeighbors]->getNbParticles();
            const FReal*const sourcesPhysicalValues = inNeighbors[idxNeighbors]->getPhysicalValues();
            const FReal*const sourcesX = inNeighbors[idxNeighbors]->getPositions()[0];
            const FReal*const sourcesY = inNeighbors[idxNeighbors]->getPositions()[1];
            const FReal*const sourcesZ = inNeighbors[idxNeighbors]->getPositions()[2];

            for(FSize idxTarget = 0 ; idxTarget < nbParticlesTargets ; ++idxTarget){
                FSize idxSource = 0;
                {
                    const FSize nbVectorizedInteractions = (nbParticlesSources/NbFRealInComputeClass)*NbFRealInComputeClass;
                    const ComputeClass tx = ComputeClass(targetsX[idxTarget]);
                    const ComputeClass ty = ComputeClass(targetsY[idxTarget]);
                    const ComputeClass tz = ComputeClass(targetsZ[idxTarget]);
                    const ComputeClass tv = ComputeClass(targetsPhysicalValues[idxTarget]);
                    ComputeClass  tfx = ComputeClass::GetZero();
                    ComputeClass  tfy = ComputeClass::GetZero();
                    ComputeClass  tfz = ComputeClass::GetZero();
                    ComputeClass  tpo = ComputeClass::GetZero();

                    for( ; idxSource < nbVectorizedInteractions ; idxSource += NbFRealInComputeClass){
                        ComputeClass Kxy[1];
                        ComputeClass dKxy[3];
                        MatrixKernel->evaluateBlockAndDerivative(tx,ty,tz,
                                                                 ComputeClass(&sourcesX[idxSource]),
                                                                 ComputeClass(&sourcesY[idxSource]),
                                                                 ComputeClass(&sourcesZ[idxSource]),
                                                                 Kxy,dKxy);
                        const ComputeClass coef = (tv * ComputeClass(&sourcesPhysicalValues[idxSource]));

                        dKxy[0] *= coef;
                        dKxy[1] *= coef;
                        dKxy[2] *= coef;

                        tfx += dKxy[0];
                        tfy += dKxy[1];
                        tfz += dKxy[2];
                        tpo += Kxy[0] * ComputeClass(&sourcesPhysicalValues[idxSource]);
                    }

                    targetsForcesX[idxTarget] += tfx.horizontalSum();
                    targetsForcesY[idxTarget] += tfy.horizontalSum();
                    targetsForcesZ[idxTarget] += tfz.horizontalSum();
                    targetsPotentials[idxTarget] += tpo.horizontalSum();
                }
                {
                    const FReal tx = FReal(targetsX[idxTarget]);
                    const FReal ty = FReal(targetsY[idxTarget]);
                    const FReal tz = FReal(targetsZ[idxTarget]);
                    const FReal tv = FReal(targetsPhysicalValues[idxTarget]);
                    FReal  tfx = FReal(0.);
                    FReal  tfy = FReal(0.);
                    FReal  tfz = FReal(0.);
                    FReal  tpo = FReal(0.);

                    for( ; idxSource < nbParticlesSources ; idxSource += 1){
                        FReal Kxy[1];
                        FReal dKxy[3];
                        MatrixKernel->evaluateBlockAndDerivative(tx,ty,tz,
                                                                 FReal(sourcesX[idxSource]),
                                                                 FReal(sourcesY[idxSource]),
                                                                 FReal(sourcesZ[idxSource]),
                                                                 Kxy,dKxy);
                        const FReal coef = (tv * sourcesPhysicalValues[idxSource]);

                        dKxy[0] *= coef;
                        dKxy[1] *= coef;
                        dKxy[2] *= coef;

                        tfx += dKxy[0];
                        tfy += dKxy[1];
                        tfz += dKxy[2];
                        tpo += Kxy[0] * sourcesPhysicalValues[idxSource];
                    }

                    targetsForcesX[idxTarget] += tfx;
                    targetsForcesY[idxTarget] += tfy;
                    targetsForcesZ[idxTarget] += tfz;
                    targetsPotentials[idxTarget] += tpo;
                }
            }
        }
    }
}

} // End namespace

template <class FReal>
struct FP2PT{
};

#include "InastempCompileConfig.h"

template <>
struct FP2PT<double>{
    template <class ContainerClass, class MatrixKernelClass>
    static void FullMutual(ContainerClass* const FRestrict inTargets, ContainerClass* const inNeighbors[],
                           const int limiteNeighbors, const MatrixKernelClass *const MatrixKernel){
        FP2P::GenericFullMutual<double, ContainerClass, MatrixKernelClass, InaVecBestTypeDouble, InaVecBestTypeDouble::VecLength>(inTargets, inNeighbors, limiteNeighbors, MatrixKernel);
    }


    template <class ContainerClass, class MatrixKernelClass>
    static void Inner(ContainerClass* const FRestrict inTargets, const MatrixKernelClass *const MatrixKernel){
        FP2P::GenericInner<double, ContainerClass, MatrixKernelClass, InaVecBestTypeDouble, InaVecBestTypeDouble::VecLength>(inTargets, MatrixKernel);
    }

    template <class ContainerClass, class MatrixKernelClass>
    static void FullRemote(ContainerClass* const FRestrict inTargets, const ContainerClass* const inNeighbors[],
                           const int limiteNeighbors, const MatrixKernelClass *const MatrixKernel){
        FP2P::GenericFullRemote<double, ContainerClass, MatrixKernelClass, InaVecBestTypeDouble, InaVecBestTypeDouble::VecLength>(inTargets, inNeighbors, limiteNeighbors, MatrixKernel);
    }
};

template <>
struct FP2PT<float>{
    template <class ContainerClass, class MatrixKernelClass>
    static void FullMutual(ContainerClass* const FRestrict inTargets, ContainerClass* const inNeighbors[],
                           const int limiteNeighbors, const MatrixKernelClass *const MatrixKernel){
        FP2P::GenericFullMutual<float, ContainerClass, MatrixKernelClass, InaVecBestTypeFloat, InaVecBestTypeFloat::VecLength>(inTargets, inNeighbors, limiteNeighbors, MatrixKernel);
    }

    template <class ContainerClass, class MatrixKernelClass>
    static void Inner(ContainerClass* const FRestrict inTargets, const MatrixKernelClass *const MatrixKernel){
        FP2P::GenericInner<float, ContainerClass, MatrixKernelClass, InaVecBestTypeFloat, InaVecBestTypeFloat::VecLength>(inTargets, MatrixKernel);
    }

    template <class ContainerClass, class MatrixKernelClass>
    static void FullRemote(ContainerClass* const FRestrict inTargets, const ContainerClass* const inNeighbors[],
                           const int limiteNeighbors, const MatrixKernelClass *const MatrixKernel){
        FP2P::GenericFullRemote<float, ContainerClass, MatrixKernelClass, InaVecBestTypeFloat, InaVecBestTypeFloat::VecLength>(inTargets, inNeighbors, limiteNeighbors, MatrixKernel);
    }
};


#include "FP2PTensorialKij.hpp"

#include "FP2PMultiRhs.hpp"

#endif // FP2P_HPP
