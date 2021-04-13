// See LICENCE file at project root
#ifndef FP2P_i_HPP
#define FP2P_i_HPP

#include "Utils/FPoint.hpp"
#include <math.h>

namespace FP2P_i {

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
   */

template <class FReal, typename MatrixKernelClass>
inline void MutualParticles_i(const FReal targetX,const FReal targetY,const FReal targetZ, const FReal targetPhysicalValue,
                            FReal* targetForceX_real, FReal* targetForceY_real, FReal* targetForceZ_real, FReal* targetPotential_real, FReal* targetForceX_imag, FReal* targetForceY_imag, FReal* targetForceZ_imag, FReal* targetPotential_imag,
                            const FReal sourceX,const FReal sourceY,const FReal sourceZ, const FReal sourcePhysicalValue,
                            FReal* sourceForceX_real, FReal* sourceForceY_real, FReal* sourceForceZ_real, FReal* sourcePotential_real, FReal* sourceForceX_imag, FReal* sourceForceY_imag, FReal* sourceForceZ_imag, FReal* sourcePotential_imag,
                            const MatrixKernelClass *const MatrixKernel)
{

				
    // Compute kernel of interaction...
    const FPoint<FReal> sourcePoint(sourceX,sourceY,sourceZ);
    const FPoint<FReal> targetPoint(targetX,targetY,targetZ);
    FReal Kxy[2];
    FReal dKxy[6];
    MatrixKernel->evaluateBlockAndDerivative(targetPoint,sourcePoint,Kxy,dKxy);
    const FReal mutual_coeff = MatrixKernel->getMutualCoefficient(); // 1 if symmetric; -1 if antisymmetric

    FReal coef = (targetPhysicalValue * sourcePhysicalValue);

    (*targetForceX_real) += dKxy[0] * coef;
    (*targetForceY_real) += dKxy[1] * coef;
    (*targetForceZ_real) += dKxy[2] * coef;

    (*targetForceX_imag) += dKxy[3] * coef;
    (*targetForceY_imag) += dKxy[4] * coef;
    (*targetForceZ_imag) += dKxy[5] * coef;

	
    (*targetPotential_real) += ( Kxy[0] * sourcePhysicalValue );
    (*targetPotential_imag) += ( Kxy[1] * sourcePhysicalValue );
	

    (*sourceForceX_real) -= dKxy[0] * coef;
    (*sourceForceY_real) -= dKxy[1] * coef;
    (*sourceForceZ_real) -= dKxy[2] * coef;

    (*sourceForceX_imag) -= dKxy[3] * coef;
    (*sourceForceY_imag) -= dKxy[4] * coef;
    (*sourceForceZ_imag) -= dKxy[5] * coef;


    (*sourcePotential_real) += ( mutual_coeff * Kxy[0] * targetPhysicalValue );	
    (*sourcePotential_imag) += ( mutual_coeff * Kxy[1] * targetPhysicalValue );
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
inline void NonMutualParticles_i(const FReal targetX,const FReal targetY,const FReal targetZ, const FReal targetPhysicalValue,
                               FReal* targetForceX_real, FReal* targetForceY_real, FReal* targetForceZ_real, FReal* targetPotential_real, FReal* targetForceX_imag, FReal* targetForceY_imag, FReal* targetForceZ_imag, FReal* targetPotential_imag,
                               const FReal sourceX,const FReal sourceY,const FReal sourceZ, const FReal sourcePhysicalValue,
                               const MatrixKernelClass *const MatrixKernel)
{


    // Compute kernel of interaction...
    const FPoint<FReal> sourcePoint(sourceX,sourceY,sourceZ);
    const FPoint<FReal> targetPoint(targetX,targetY,targetZ);
    FReal Kxy[2];
    FReal dKxy[6];
    MatrixKernel->evaluateBlockAndDerivative(targetPoint,sourcePoint,Kxy,dKxy);

    FReal coef = (targetPhysicalValue * sourcePhysicalValue);

    (*targetForceX_real) += dKxy[0] * coef;
    (*targetForceY_real) += dKxy[1] * coef;
    (*targetForceZ_real) += dKxy[2] * coef;

    (*targetForceX_imag) += dKxy[3] * coef;
    (*targetForceY_imag) += dKxy[4] * coef;
    (*targetForceZ_imag) += dKxy[5] * coef;

    (*targetPotential_real) += ( Kxy[0] * sourcePhysicalValue );	
    (*targetPotential_imag) += ( Kxy[1] * sourcePhysicalValue );
}

//------------------------------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------------------------------

// 0 -- getPhysicalValues
// 1 -- getPotentials real
// 2 -- getForcesX real
// 3 -- getForcesY real 
// 4 -- getForcesZ real
// 5 -- getPotentials imag
// 6 -- getForcesX imag 
// 7 -- getForcesY imag 
// 8 -- getForcesZ imag 

template <class FReal, class ContainerClass, class MatrixKernelClass, class ComputeClass, int NbFRealInComputeClass>
static void GenericFullMutual_i(ContainerClass* const FRestrict inTargets,
                              ContainerClass* const inNeighbors[],
                              const int limiteNeighbors,
                              const MatrixKernelClass *const MatrixKernel)
{


    const FSize nbParticlesTargets = inTargets->getNbParticles();
    const FReal*const targetsPhysicalValues = inTargets->getPhysicalValues();
    const FReal*const targetsX = inTargets->getPositions()[0];
    const FReal*const targetsY = inTargets->getPositions()[1];
    const FReal*const targetsZ = inTargets->getPositions()[2];
	
    FReal* targetsPotentials_real = inTargets->getPotentials_real();  			//1	
    FReal* targetsForcesX_real = inTargets->getForcesX_real();  				//2	
    FReal* targetsForcesY_real = inTargets->getForcesY_real();  				//3	
    FReal* targetsForcesZ_real = inTargets->getForcesZ_real();  				//4
	

    FReal* targetsPotentials_imag = inTargets->getPotentials_imag();  			//5		
    FReal* targetsForcesX_imag = inTargets->getForcesX_imag();  				//6
    FReal* targetsForcesY_imag = inTargets->getForcesY_imag();  				//7
    FReal* targetsForcesZ_imag = inTargets->getForcesZ_imag();  				//8


    for(FSize idxNeighbors = 0 ; idxNeighbors < limiteNeighbors ; ++idxNeighbors){
		int count = 0;
        if( inNeighbors[idxNeighbors] ){
            const FSize nbParticlesSources = inNeighbors[idxNeighbors]->getNbParticles();
            const FReal*const sourcesPhysicalValues = inNeighbors[idxNeighbors]->getPhysicalValues();
            const FReal*const sourcesX = inNeighbors[idxNeighbors]->getPositions()[0];
            const FReal*const sourcesY = inNeighbors[idxNeighbors]->getPositions()[1];
            const FReal*const sourcesZ = inNeighbors[idxNeighbors]->getPositions()[2];

			FReal* sourcesPotentials_real = inNeighbors[idxNeighbors]->getPotentials_real();  		        //1			
			FReal* sourcesForcesX_real = inNeighbors[idxNeighbors]->getForcesX_real();  				    //2
			FReal* sourcesForcesY_real = inNeighbors[idxNeighbors]->getForcesY_real();  				    //3		
			FReal* sourcesForcesZ_real = inNeighbors[idxNeighbors]->getForcesZ_real();  				    //4
			

			FReal* sourcesPotentials_imag = inNeighbors[idxNeighbors]->getPotentials_imag();  		        //5			
			FReal* sourcesForcesX_imag = inNeighbors[idxNeighbors]->getForcesX_imag();  				    //6
			FReal* sourcesForcesY_imag = inNeighbors[idxNeighbors]->getForcesY_imag();  				    //7
			FReal* sourcesForcesZ_imag = inNeighbors[idxNeighbors]->getForcesZ_imag();  				    //8
			

            for(FSize idxTarget = 0 ; idxTarget < nbParticlesTargets ; ++idxTarget){
                FSize idxSource = 0;
                {
						//	std::cout <<"DOES PRINT !!!!!!  DOES PRINT !!!!!! DOES PRINT !!!!!! DOES PRINT !!!!!! 0" << std::endl; // 						
                    const FSize nbVectorizedInteractions = (nbParticlesSources/NbFRealInComputeClass)*NbFRealInComputeClass;
                    const ComputeClass tx = ComputeClass(targetsX[idxTarget]);
                    const ComputeClass ty = ComputeClass(targetsY[idxTarget]);
                    const ComputeClass tz = ComputeClass(targetsZ[idxTarget]);
                    const ComputeClass tv = ComputeClass(targetsPhysicalValues[idxTarget]);

					ComputeClass  tpo_real = ComputeClass::GetZero();	 		        //1			
					ComputeClass  tfx_real = ComputeClass::GetZero();	 		        //2
					ComputeClass  tfy_real = ComputeClass::GetZero();	 		        //3
					ComputeClass  tfz_real = ComputeClass::GetZero();	 		        //4

					ComputeClass  tpo_imag = ComputeClass::GetZero();	 		        //5					
					ComputeClass  tfx_imag = ComputeClass::GetZero();	 		        //6	
					ComputeClass  tfy_imag = ComputeClass::GetZero();	 		        //7	
					ComputeClass  tfz_imag = ComputeClass::GetZero();	 		        //8	
					

                    for( ; idxSource < nbVectorizedInteractions ; idxSource += NbFRealInComputeClass){
                        ComputeClass Kxy[2];
                        ComputeClass dKxy[6];

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
		
                        dKxy[3] *= coef;
                        dKxy[4] *= coef;
                        dKxy[5] *= coef;
						
						
                        tfx_real += dKxy[0];
                        tfy_real += dKxy[1];
                        tfz_real += dKxy[2];
						
                        tfx_imag += dKxy[3];
                        tfy_imag += dKxy[4];
                        tfz_imag += dKxy[5];


                        tpo_real += Kxy[0] * ComputeClass(&sourcesPhysicalValues[idxSource]);						
                        tpo_imag += Kxy[1] * ComputeClass(&sourcesPhysicalValues[idxSource]);


                        (ComputeClass(&sourcesForcesX_real[idxSource]) - dKxy[0]).storeInArray(&sourcesForcesX_real[idxSource]);
                        (ComputeClass(&sourcesForcesY_real[idxSource]) - dKxy[1]).storeInArray(&sourcesForcesY_real[idxSource]);
                        (ComputeClass(&sourcesForcesZ_real[idxSource]) - dKxy[2]).storeInArray(&sourcesForcesZ_real[idxSource]);

                        (ComputeClass(&sourcesForcesX_imag[idxSource]) - dKxy[3]).storeInArray(&sourcesForcesX_imag[idxSource]);
                        (ComputeClass(&sourcesForcesY_imag[idxSource]) - dKxy[4]).storeInArray(&sourcesForcesY_imag[idxSource]);
                        (ComputeClass(&sourcesForcesZ_imag[idxSource]) - dKxy[5]).storeInArray(&sourcesForcesZ_imag[idxSource]);

                        (ComputeClass(&sourcesPotentials_real[idxSource]) + mutual_coeff * Kxy[0] * tv).storeInArray(&sourcesPotentials_real[idxSource]);						
                        (ComputeClass(&sourcesPotentials_imag[idxSource]) + mutual_coeff * Kxy[1] * tv).storeInArray(&sourcesPotentials_imag[idxSource]);				

						
                    }


					targetsForcesX_real[idxTarget] += tfx_real.horizontalSum();					
					targetsForcesY_real[idxTarget] += tfy_real.horizontalSum();		
					targetsForcesZ_real[idxTarget] += tfz_real.horizontalSum();		
					
                    targetsForcesX_imag[idxTarget] += tfx_imag.horizontalSum();
                    targetsForcesY_imag[idxTarget] += tfy_imag.horizontalSum();
                    targetsForcesZ_imag[idxTarget] += tfz_imag.horizontalSum();		

					
                    targetsPotentials_real[idxTarget] += tpo_real.horizontalSum();					
                    targetsPotentials_imag[idxTarget] += tpo_imag.horizontalSum();
                }
                {
                    const FReal tx = FReal(targetsX[idxTarget]);
                    const FReal ty = FReal(targetsY[idxTarget]);
                    const FReal tz = FReal(targetsZ[idxTarget]);
                    const FReal tv = FReal(targetsPhysicalValues[idxTarget]);
					
                    FReal  tfx_real = FReal(0.);
                    FReal  tfy_real = FReal(0.);
                    FReal  tfz_real = FReal(0.);
					
                    FReal  tfx_imag = FReal(0.);
                    FReal  tfy_imag = FReal(0.);
                    FReal  tfz_imag = FReal(0.);

                    FReal  tpo_real = FReal(0.);					
                    FReal  tpo_imag = FReal(0.);

                    for( ; idxSource < nbParticlesSources ; idxSource += 1){
                        FReal Kxy[2];
                        FReal dKxy[6];					
                        
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

                        dKxy[3] *= coef;
                        dKxy[4] *= coef;
                        dKxy[5] *= coef;						

                        tfx_real += dKxy[0];
                        tfy_real += dKxy[1];
                        tfz_real += dKxy[2];

                        tfx_imag += dKxy[3];
                        tfy_imag += dKxy[4];
                        tfz_imag += dKxy[5];


                        tpo_real += Kxy[0] * FReal(sourcesPhysicalValues[idxSource]);						
                        tpo_imag += Kxy[1] * FReal(sourcesPhysicalValues[idxSource]);
						
                        sourcesForcesX_real[idxSource] -= dKxy[0];
                        sourcesForcesY_real[idxSource] -= dKxy[1];
                        sourcesForcesZ_real[idxSource] -= dKxy[2];
					
                        sourcesForcesX_imag[idxSource] -= dKxy[3];
                        sourcesForcesY_imag[idxSource] -= dKxy[4];
                        sourcesForcesZ_imag[idxSource] -= dKxy[5];			

                        sourcesPotentials_real[idxSource] += mutual_coeff * Kxy[0] * tv;						
                        sourcesPotentials_imag[idxSource] += mutual_coeff * Kxy[1] * tv;
                    }

                    targetsForcesX_real[idxTarget] += tfx_real;					
                    targetsForcesY_real[idxTarget] += tfy_real;
                    targetsForcesZ_real[idxTarget] += tfz_real;

                    targetsForcesX_imag[idxTarget] += tfx_imag;
                    targetsForcesY_imag[idxTarget] += tfy_imag;
                    targetsForcesZ_imag[idxTarget] += tfz_imag;
				
                    targetsPotentials_real[idxTarget] += tpo_real;					
                    targetsPotentials_imag[idxTarget] += tpo_imag;
                }
            }
        count = count + 1;
		}
    }

}



// 0 -- getPhysicalValues
// 1 -- getPotentials real
// 2 -- getForcesX real
// 3 -- getForcesY real 
// 4 -- getForcesZ real
// 5 -- getPotentials imag
// 6 -- getForcesX imag 
// 7 -- getForcesY imag 
// 8 -- getForcesZ imag 


template <class FReal, class ContainerClass, class MatrixKernelClass, class ComputeClass, int NbFRealInComputeClass>
static void GenericInner_i(ContainerClass* const FRestrict inTargets, const MatrixKernelClass *const MatrixKernel)
{

    const FSize nbParticlesTargets = inTargets->getNbParticles();
    const FReal*const targetsPhysicalValues = inTargets->getPhysicalValues();
    const FReal*const targetsX = inTargets->getPositions()[0];
    const FReal*const targetsY = inTargets->getPositions()[1];
    const FReal*const targetsZ = inTargets->getPositions()[2];

    FReal*const targetsPotentials_real = inTargets->getPotentials_real();  			//1	
    FReal*const targetsForcesX_real = inTargets->getForcesX_real();  				//2
    FReal*const targetsForcesY_real = inTargets->getForcesY_real();  				//3
    FReal*const targetsForcesZ_real = inTargets->getForcesZ_real();  				//4

    FReal*const targetsPotentials_imag = inTargets->getPotentials_imag();  			//5		
    FReal*const targetsForcesX_imag = inTargets->getForcesX_imag();  				//6
    FReal*const targetsForcesY_imag = inTargets->getForcesY_imag();  				//7
    FReal*const targetsForcesZ_imag = inTargets->getForcesZ_imag();  				//8
	

    {//In this part, we compute (vectorially) the interaction
        //within the target leaf.

        const FSize nbParticlesSources = nbParticlesTargets;
        const FReal*const sourcesPhysicalValues = targetsPhysicalValues;
        const FReal*const sourcesX = targetsX;
        const FReal*const sourcesY = targetsY;
        const FReal*const sourcesZ = targetsZ;

        FReal*const sourcesPotentials_real = targetsPotentials_real;  		        //1			
        FReal*const sourcesForcesX_real = targetsForcesX_real;  				    //2
        FReal*const sourcesForcesY_real = targetsForcesY_real;  				    //3
        FReal*const sourcesForcesZ_real = targetsForcesZ_real;  				    //4

        FReal*const sourcesPotentials_imag = targetsPotentials_imag;  		        //5			
        FReal*const sourcesForcesX_imag = targetsForcesX_imag;  				    //6
        FReal*const sourcesForcesY_imag = targetsForcesY_imag;  				    //7
        FReal*const sourcesForcesZ_imag = targetsForcesZ_imag;  				    //8
	

        for(FSize idxTarget = 0 ; idxTarget < nbParticlesTargets ; ++idxTarget){
            FSize idxSource = idxTarget+1;
            {
                const FSize nbVectorizedInteractions = ((nbParticlesSources-idxSource)/NbFRealInComputeClass)*NbFRealInComputeClass + idxSource;
                const ComputeClass tx = ComputeClass(targetsX[idxTarget]);
                const ComputeClass ty = ComputeClass(targetsY[idxTarget]);
                const ComputeClass tz = ComputeClass(targetsZ[idxTarget]);
                const ComputeClass tv = ComputeClass(targetsPhysicalValues[idxTarget]);

                ComputeClass  tpo_real = ComputeClass::GetZero();	 		        //1			
                ComputeClass  tfx_real = ComputeClass::GetZero();	 		        //2
                ComputeClass  tfy_real = ComputeClass::GetZero();	 		        //3
                ComputeClass  tfz_real = ComputeClass::GetZero();	 		        //4

                ComputeClass  tpo_imag = ComputeClass::GetZero();	 		        //5					
                ComputeClass  tfx_imag = ComputeClass::GetZero();	 		        //6	
                ComputeClass  tfy_imag = ComputeClass::GetZero();	 		        //7	
                ComputeClass  tfz_imag = ComputeClass::GetZero();	 		        //8	


                for( ; idxSource < nbVectorizedInteractions ; idxSource += NbFRealInComputeClass){
                    ComputeClass Kxy[2];
                    ComputeClass dKxy[6];
					
					
                    MatrixKernel->evaluateBlockAndDerivative(tx,ty,tz,
                                                             ComputeClass(&sourcesX[idxSource]),
                                                             ComputeClass(&sourcesY[idxSource]),
                                                             ComputeClass(&sourcesZ[idxSource]),
                                                             Kxy,dKxy);
                    const ComputeClass mutual_coeff = ComputeClass(MatrixKernel->getMutualCoefficient()); // 1 if symmetric; -1 if antisymmetric

                    const ComputeClass coef = (tv * ComputeClass(&sourcesPhysicalValues[idxSource]));

                    dKxy[0] *= coef;	 		        //0   		Ptot1_real		// Xreal
                    dKxy[1] *= coef;	 		        //1    						// Yreal
                    dKxy[2] *= coef;	 		        //2			Ptot2_real		// Zreal
					
                    dKxy[3] *= coef;	 		        //3			Ptot1_imag		// Ximag
                    dKxy[4] *= coef;	 		        //4   						// Yimag
                    dKxy[5] *= coef;	 		        //5   		Ptot2_imag		// Zimag
					

                    tfx_real += dKxy[0];
                    tfy_real += dKxy[1];
                    tfz_real += dKxy[2];

                    tfx_imag += dKxy[3];
                    tfy_imag += dKxy[4];
                    tfz_imag += dKxy[5];					
					
					
                    tpo_real += Kxy[0]*ComputeClass(&sourcesPhysicalValues[idxSource]);
                    tpo_imag += Kxy[1]*ComputeClass(&sourcesPhysicalValues[idxSource]);
										
                    (ComputeClass(&sourcesForcesX_real[idxSource]) - dKxy[0]).storeInArray(&sourcesForcesX_real[idxSource]);		//IS THE PROBLEM HERE???
                    (ComputeClass(&sourcesForcesY_real[idxSource]) - dKxy[1]).storeInArray(&sourcesForcesY_real[idxSource]);
                    (ComputeClass(&sourcesForcesZ_real[idxSource]) - dKxy[2]).storeInArray(&sourcesForcesZ_real[idxSource]);

	
                    (ComputeClass(&sourcesForcesX_imag[idxSource]) - dKxy[3]).storeInArray(&sourcesForcesX_imag[idxSource]);
                    (ComputeClass(&sourcesForcesY_imag[idxSource]) - dKxy[4]).storeInArray(&sourcesForcesY_imag[idxSource]);
                    (ComputeClass(&sourcesForcesZ_imag[idxSource]) - dKxy[5]).storeInArray(&sourcesForcesZ_imag[idxSource]);	
					
                    (ComputeClass(&sourcesPotentials_real[idxSource]) + mutual_coeff * Kxy[0] * tv).storeInArray(&sourcesPotentials_real[idxSource]);
                    (ComputeClass(&sourcesPotentials_imag[idxSource]) + mutual_coeff * Kxy[1] * tv).storeInArray(&sourcesPotentials_imag[idxSource]);					
 

				}


                targetsForcesX_real[idxTarget] += tfx_real.horizontalSum();
                targetsForcesY_real[idxTarget] += tfy_real.horizontalSum();
                targetsForcesZ_real[idxTarget] += tfz_real.horizontalSum();
				
                targetsForcesX_imag[idxTarget] += tfx_imag.horizontalSum();
                targetsForcesY_imag[idxTarget] += tfy_imag.horizontalSum();
                targetsForcesZ_imag[idxTarget] += tfz_imag.horizontalSum();
				
                targetsPotentials_real[idxTarget] += tpo_real.horizontalSum();
                targetsPotentials_imag[idxTarget] += tpo_imag.horizontalSum();
            }
            {
                const FReal tx = FReal(targetsX[idxTarget]);
                const FReal ty = FReal(targetsY[idxTarget]);
                const FReal tz = FReal(targetsZ[idxTarget]);
                const FReal tv = FReal(targetsPhysicalValues[idxTarget]);
				
                FReal  tpo_real = FReal(0.);		//1
                FReal  tfx_real = FReal(0.);		//2
                FReal  tfy_real = FReal(0.);		//3
                FReal  tfz_real = FReal(0.);		//4
				
                FReal  tpo_imag = FReal(0.);		//5
                FReal  tfx_imag = FReal(0.);		//6
                FReal  tfy_imag = FReal(0.);		//7
                FReal  tfz_imag = FReal(0.);		//8
				

                for( ; idxSource < nbParticlesSources ; idxSource += 1){
                    FReal Kxy[2];
                    FReal dKxy[6];
			
							
                    MatrixKernel->evaluateBlockAndDerivative(tx,ty,tz,
                                                             FReal(sourcesX[idxSource]),
                                                             FReal(sourcesY[idxSource]),
                                                             FReal(sourcesZ[idxSource]),
                                                             Kxy,dKxy);
															 
                    const FReal mutual_coeff = FReal(MatrixKernel->getMutualCoefficient()); // 1 if symmetric; -1 if antisymmetric

                    const FReal coef = (tv * FReal(sourcesPhysicalValues[idxSource]));
					
                    dKxy[0] *= coef;	 		        //0   		Ptot1_real		// Xreal
                    dKxy[1] *= coef;	 		        //1    						// Yreal
                    dKxy[2] *= coef;	 		        //2			Ptot2_real		// Zreal
					
                    dKxy[3] *= coef;	 		        //3			Ptot1_imag		// Ximag
                    dKxy[4] *= coef;	 		        //4   						// Yimag
                    dKxy[5] *= coef;	 		        //5   		Ptot2_imag		// Zimag					


                    tfx_real += dKxy[0];		
                    tfy_real += dKxy[1];
                    tfz_real += dKxy[2];			

                    tfx_imag += dKxy[3];
                    tfy_imag += dKxy[4];
                    tfz_imag += dKxy[5];


                    tpo_real += Kxy[0]*sourcesPhysicalValues[idxSource];					
                    tpo_imag += Kxy[1]*sourcesPhysicalValues[idxSource];
									
                    sourcesForcesX_real[idxSource] -= dKxy[0];		
                    sourcesForcesY_real[idxSource] -= dKxy[1];
                    sourcesForcesZ_real[idxSource] -= dKxy[2];		
				
                    sourcesForcesX_imag[idxSource] -= dKxy[3];
                    sourcesForcesY_imag[idxSource] -= dKxy[4];
                    sourcesForcesZ_imag[idxSource] -= dKxy[5];		
					
                    sourcesPotentials_real[idxSource] += mutual_coeff * Kxy[0] * tv;					
                    sourcesPotentials_imag[idxSource] += mutual_coeff * Kxy[1] * tv;
					
				
                }

                targetsForcesX_real[idxTarget] += tfx_real;				
                targetsForcesY_real[idxTarget] += tfy_real;	
                targetsForcesZ_real[idxTarget] += tfz_real;
				
                targetsForcesX_imag[idxTarget] += tfx_imag;
                targetsForcesY_imag[idxTarget] += tfy_imag;
                targetsForcesZ_imag[idxTarget] += tfz_imag;	
				
                targetsPotentials_real[idxTarget] += tpo_real;
                targetsPotentials_imag[idxTarget] += tpo_imag;				
            }
        }
    }
	
}

/// NOTE: GENERICFULLREMOTE IS NOT A PART OF THE NAMESPACE!!!

// 0 -- getPhysicalValues
// 1 -- getPotentials real
// 2 -- getForcesX real
// 3 -- getForcesY real 
// 4 -- getForcesZ real
// 5 -- getPotentials imag
// 6 -- getForcesX imag 
// 7 -- getForcesY imag 
// 8 -- getForcesZ imag 




template <class FReal, class ContainerClass, class MatrixKernelClass, class ComputeClass, int NbFRealInComputeClass>
static void GenericFullRemote_i(ContainerClass* const FRestrict inTargets, const ContainerClass* const inNeighbors[],
                              const int limiteNeighbors, const MatrixKernelClass *const MatrixKernel)
{
		
    const FSize nbParticlesTargets = inTargets->getNbParticles();
    const FReal*const targetsPhysicalValues = inTargets->getPhysicalValues();
    const FReal*const targetsX = inTargets->getPositions()[0];
    const FReal*const targetsY = inTargets->getPositions()[1];
    const FReal*const targetsZ = inTargets->getPositions()[2];

    FReal*const targetsPotentials_real = inTargets->getPotentials_real();  			//1	
    FReal*const targetsForcesX_real = inTargets->getForcesX_real();  				//2
    FReal*const targetsForcesY_real = inTargets->getForcesY_real();  				//3
    FReal*const targetsForcesZ_real = inTargets->getForcesZ_real();  				//4

    FReal*const targetsPotentials_imag = inTargets->getPotentials_imag();  			//5		
    FReal*const targetsForcesX_imag = inTargets->getForcesX_imag();  				//6
    FReal*const targetsForcesY_imag = inTargets->getForcesY_imag();  				//7
    FReal*const targetsForcesZ_imag = inTargets->getForcesZ_imag();  				//8
	

	
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
					
					ComputeClass  tpo_real = ComputeClass::GetZero();	 		        //1			
					ComputeClass  tfx_real = ComputeClass::GetZero();	 		        //2
					ComputeClass  tfy_real = ComputeClass::GetZero();	 		        //3
					ComputeClass  tfz_real = ComputeClass::GetZero();	 		        //4

					ComputeClass  tpo_imag = ComputeClass::GetZero();	 		        //5					
					ComputeClass  tfx_imag = ComputeClass::GetZero();	 		        //6	
					ComputeClass  tfy_imag = ComputeClass::GetZero();	 		        //7	
					ComputeClass  tfz_imag = ComputeClass::GetZero();	 		        //8	

                    for( ; idxSource < nbVectorizedInteractions ; idxSource += NbFRealInComputeClass){
                        ComputeClass Kxy[2];
                        ComputeClass dKxy[6];
                        MatrixKernel->evaluateBlockAndDerivative(tx,ty,tz,
                                                                 ComputeClass(&sourcesX[idxSource]),
                                                                 ComputeClass(&sourcesY[idxSource]),
                                                                 ComputeClass(&sourcesZ[idxSource]),
                                                                 Kxy,dKxy);
                        const ComputeClass coef = (tv * ComputeClass(&sourcesPhysicalValues[idxSource]));

                        dKxy[0] *= coef;	 		        //0   		Ptot1_real		// Xreal
                        dKxy[1] *= coef;	 		        //1    						// Yreal
                        dKxy[2] *= coef;	 		        //2			Ptot2_real		// Zreal
						
                        dKxy[3] *= coef;	 		        //3			Ptot1_imag		// Ximag
                        dKxy[4] *= coef;	 		        //4   						// Yimag
                        dKxy[5] *= coef;	 		        //5   		Ptot2_imag		// Zimag						

                        tfx_real += dKxy[0];
                        tfy_real += dKxy[1];
                        tfz_real += dKxy[2];
						
                        tfx_imag += dKxy[3];
                        tfy_imag += dKxy[4];
                        tfz_imag += dKxy[5];


                        tpo_real += Kxy[0] * ComputeClass(&sourcesPhysicalValues[idxSource]);						
                        tpo_imag += Kxy[1] * ComputeClass(&sourcesPhysicalValues[idxSource]);
                    }

                    targetsForcesX_real[idxTarget] += tfx_real.horizontalSum();
                    targetsForcesY_real[idxTarget] += tfy_real.horizontalSum();
                    targetsForcesZ_real[idxTarget] += tfz_real.horizontalSum();

                    targetsForcesX_imag[idxTarget] += tfx_imag.horizontalSum();
                    targetsForcesY_imag[idxTarget] += tfy_imag.horizontalSum();
                    targetsForcesZ_imag[idxTarget] += tfz_imag.horizontalSum();
					
                    targetsPotentials_real[idxTarget] += tpo_real.horizontalSum();					
                    targetsPotentials_imag[idxTarget] += tpo_imag.horizontalSum();					
                }
                {
                    const FReal tx = FReal(targetsX[idxTarget]);
                    const FReal ty = FReal(targetsY[idxTarget]);
                    const FReal tz = FReal(targetsZ[idxTarget]);
                    const FReal tv = FReal(targetsPhysicalValues[idxTarget]);
                    FReal  tfx_real = FReal(0.);
                    FReal  tfy_real = FReal(0.);
                    FReal  tfz_real = FReal(0.);					
                    FReal  tpo_real = FReal(0.);
					
                    FReal  tfx_imag = FReal(0.);
                    FReal  tfy_imag = FReal(0.);
                    FReal  tfz_imag = FReal(0.);					
                    FReal  tpo_imag = FReal(0.);					

                    for( ; idxSource < nbParticlesSources ; idxSource += 1){
                        FReal Kxy[2];
                        FReal dKxy[6];
                        MatrixKernel->evaluateBlockAndDerivative(tx,ty,tz,
                                                                 FReal(sourcesX[idxSource]),
                                                                 FReal(sourcesY[idxSource]),
                                                                 FReal(sourcesZ[idxSource]),
                                                                 Kxy,dKxy);
                        const FReal coef = (tv * sourcesPhysicalValues[idxSource]);

                        dKxy[0] *= coef;
                        dKxy[1] *= coef;
                        dKxy[2] *= coef;
						
                        dKxy[3] *= coef;
                        dKxy[4] *= coef;
                        dKxy[5] *= coef;
						
                        tfx_real += dKxy[0];
                        tfy_real += dKxy[1];
                        tfz_real += dKxy[2];
						
                        tfx_imag += dKxy[3];
                        tfy_imag += dKxy[4];
                        tfz_imag += dKxy[5];

                        tpo_real += Kxy[0] * sourcesPhysicalValues[idxSource];						
                        tpo_imag += Kxy[1] * sourcesPhysicalValues[idxSource];
                    }

                    targetsForcesX_real[idxTarget] += tfx_real;
                    targetsForcesY_real[idxTarget] += tfy_real;
                    targetsForcesZ_real[idxTarget] += tfz_real;
					
                    targetsForcesX_imag[idxTarget] += tfx_imag;
                    targetsForcesY_imag[idxTarget] += tfy_imag;
                    targetsForcesZ_imag[idxTarget] += tfz_imag;

                    targetsPotentials_real[idxTarget] += tpo_real;					
                    targetsPotentials_imag[idxTarget] += tpo_imag;			

                }
            }
        }
    }

}

} // End namespace


template <class FReal>
struct FP2PT_i{

};

#include "InastempCompileConfig.h"

template <>
struct FP2PT_i<double>{
    template <class ContainerClass, class MatrixKernelClass>
    static void FullMutual_i(ContainerClass* const FRestrict inTargets, ContainerClass* const inNeighbors[],
                           const int limiteNeighbors, const MatrixKernelClass *const MatrixKernel)
	{	

        FP2P_i::GenericFullMutual_i<double, ContainerClass, MatrixKernelClass, InaVecBestTypeDouble, InaVecBestTypeDouble::VecLength>(inTargets, inNeighbors, limiteNeighbors, MatrixKernel);

	}
	


    template <class ContainerClass, class MatrixKernelClass>
    static void Inner_i(ContainerClass* const FRestrict inTargets, const MatrixKernelClass *const MatrixKernel)
	{

        FP2P_i::GenericInner_i<double, ContainerClass, MatrixKernelClass, InaVecBestTypeDouble, InaVecBestTypeDouble::VecLength>(inTargets, MatrixKernel);	
    }

    template <class ContainerClass, class MatrixKernelClass>
    static void FullRemote_i(ContainerClass* const FRestrict inTargets, const ContainerClass* const inNeighbors[],
                           const int limiteNeighbors, const MatrixKernelClass *const MatrixKernel)
	{

	        FP2P_i::GenericFullRemote_i<double, ContainerClass, MatrixKernelClass, InaVecBestTypeDouble, InaVecBestTypeDouble::VecLength>(inTargets, inNeighbors, limiteNeighbors, MatrixKernel);
    }
};

template <>
struct FP2PT_i<float>{
    template <class ContainerClass, class MatrixKernelClass>
    static void FullMutual_i(ContainerClass* const FRestrict inTargets, ContainerClass* const inNeighbors[],
                           const int limiteNeighbors, const MatrixKernelClass *const MatrixKernel){

        FP2P_i::GenericFullMutual_i<float, ContainerClass, MatrixKernelClass, InaVecBestTypeFloat, InaVecBestTypeFloat::VecLength>(inTargets, inNeighbors, limiteNeighbors, MatrixKernel);
    }

    template <class ContainerClass, class MatrixKernelClass>
    static void Inner_i(ContainerClass* const FRestrict inTargets, const MatrixKernelClass *const MatrixKernel){

      FP2P_i::GenericInner_i<float, ContainerClass, MatrixKernelClass, InaVecBestTypeFloat, InaVecBestTypeFloat::VecLength>(inTargets, MatrixKernel);
    }

    template <class ContainerClass, class MatrixKernelClass>
    static void FullRemote_i(ContainerClass* const FRestrict inTargets, const ContainerClass* const inNeighbors[],
                           const int limiteNeighbors, const MatrixKernelClass *const MatrixKernel){
	
    FP2P_i::GenericFullRemote_i<float, ContainerClass, MatrixKernelClass, InaVecBestTypeFloat, InaVecBestTypeFloat::VecLength>(inTargets, inNeighbors, limiteNeighbors, MatrixKernel);
    }
};


#include "FP2PTensorialKij.hpp"

#include "FP2PMultiRhs.hpp"

#endif // FP2P_HPP
