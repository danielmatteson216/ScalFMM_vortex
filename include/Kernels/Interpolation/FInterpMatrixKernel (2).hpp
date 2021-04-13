// See LICENCE file at project root
#ifndef FINTERPMATRIXKERNEL_HPP
#define FINTERPMATRIXKERNEL_HPP

#include <iostream>
#include <stdexcept>
#include <math.h>

#include "Utils/FPoint.hpp"
#include "Utils/FNoCopyable.hpp"
#include "Utils/FMath.hpp"
#include "Utils/FGlobal.hpp"

#include <sstream>
#include <fstream>
#include <stdlib.h>


// probably not extendable :)
enum KERNEL_FUNCTION_TYPE {HOMOGENEOUS, NON_HOMOGENEOUS};


/**
 * @author Matthias Messner (matthias.messner@inria.fr)
 * @author Pierre Blanchard (pierre.blanchard@inria.fr)
 * @class FInterpMatrixKernels
 * Please read the license
 *
 * This class provides the evaluators and scaling functions of the matrix
 * kernels. A matrix kernel should be understood in the sense of a kernel
 * of interaction (or the fundamental solution of a given equation).
 * It can either be scalar (NCMP=1) or tensorial (NCMP>1) depending on the
 * dimension of the equation considered. NCMP denotes the number of components
 * that are actually stored (e.g. 6 for a \f$3\times3\f$ symmetric tensor).
 *
 * Notes on the application scheme:
 * Let there be a kernel \f$K\f$ such that \f$Potential_i(X)X=K_{ij}(X,Y)PhysicalValue_j(Y)\f$
 * with \f$Potential\f$ the lhs of size NLHS and \f$PhysicalValues\f$ the rhs of size NRHS.
 * The table applyTab provides the indices in the reduced storage table
 * corresponding to the application scheme depicted earlier.
 *
 * PB: BEWARE! Homogeneous matrix kernels do not support cell width extension
 * yet. Is it possible to find a reference width and a scale factor such that
 * only 1 set of M2L opt can be used for all levels??
 *
 */
template <class FReal>
struct FInterpAbstractMatrixKernel : FNoCopyable
{ 
    virtual ~FInterpAbstractMatrixKernel(){} // to remove warning
    //virtual FReal evaluate(const FPoint<FReal>&, const FPoint<FReal>&) const = 0;
    // I need both functions because required arguments are not always given
    virtual FReal getScaleFactor(const FReal, const int) const = 0;
    virtual FReal getScaleFactor(const FReal) const = 0;
};

/// VORTEX KERNEL
template <class FReal>
struct FInterpMatrixKernelVORTEX : FInterpAbstractMatrixKernel<FReal>
{
    static const KERNEL_FUNCTION_TYPE Type = HOMOGENEOUS;
    static const unsigned int NCMP = 1; //< number of components
    static const unsigned int NPV  = 1; //< dim of physical values
    static const unsigned int NPOT = 1; //< dim of potentials
    static const unsigned int NRHS = 1; //< dim of mult exp
    static const unsigned int NLHS = 1; //< dim of loc exp

	
	//const FReal numofparticles = 10;     //11*11 = 121
	//const FReal numofparticles = 16;   //17*17 = 289
	//const FReal numofparticles = 22;   //23*23 = 529
	//const FReal numofparticles = 31;   //32*32 = 	1024
	//const FReal numofparticles = 44;   //45*45 = 2025
	//const FReal numofparticles = 62;   //63*63 = 3969
	//const FReal numofparticles = 88;   //89*89 = 7921
	//const FReal numofparticles = 124;	 //125*125 = 15625
	const FReal numofparticles = 176;  //177*177 = 31329
	//const FReal numofparticles = 248;  //249 *249 = 62001
	//const FReal numofparticles = 351;  //352 *352 = 123904
	//const FReal numofparticles = 497;  //498 *498 = 248004
	//const FReal numofparticles = 548;  //549 *549 = 301401
	
	
	const FReal rvalsq = (2 / (numofparticles * numofparticles)); //this needs to be dynamically driven by the input file...
	const double pi = FReal(M_PI);
	const FReal P2M = pi / 10;		

    FInterpMatrixKernelVORTEX() {}

    // copy ctor
    FInterpMatrixKernelVORTEX(const FInterpMatrixKernelVORTEX& /*other*/) {}

    static const char* getID() { return "ONE_OVER_R_SQUARED"; }

    static void printInfo() { std::cout << "K(x,y)=1/r^2 with r=|x-y|" << std::endl; }

    // returns position in reduced storage
    int getPosition(const unsigned int) const
    {return 0;}

    // returns coefficient of mutual interaction
    // 1 for symmetric kernels
    // -1 for antisymmetric kernels
    // Something else if other property of symmetry
    FReal getMutualCoefficient() const{ return FReal(1.); }





// ===========================================            EVALUATE           ===========================================




    // evaluate interaction
    template <class ValueClass>
    void evaluate(const ValueClass& xt, const ValueClass& yt, const ValueClass& zt, 
                        const ValueClass& xs, const ValueClass& ys, const ValueClass& zs, FReal& Ptot_real, FReal& Ptot_img) const
    {					

// difference in locations of source and target points			
		const FReal dx = (double)(xt-xs);
        const FReal dy = (double)(yt-ys);
        const FReal dz = (double)(zt-zs);		
		
        const FReal dzp = (double)(zt+zs);			
	
// RESULT --->  ptot = ( (P1 + P2) - (P3 + P4) )

//					        	(SCV) = P2
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------		
        const FReal diff = ((dx * dx) + (dz * dz));

		FReal D;		
		FReal P;		
		
		if ( rvalsq > diff / 100){
			D = FMath::Exp((-1*diff)/(rvalsq));			
			P = (D*D);
		} else{
			D = 0;
			P = 0;
		}			
		
		const FReal SCV_denom = ((P2M*dx*dx)+(P2M*dz*dz));

		FReal scv_real = (((dx*D)+(dx*(-2*P)))/SCV_denom);										//real part of SCV
		FReal scv_img = (((-1*dz*D)+(dz*(2*P)))/SCV_denom);										// imag part of SCV		

		if (SCV_denom < .000000001) {
			scv_real = 0;
			scv_img = 0;
		}
		
//					        	(1 / tan(P2M*dzz) = P1
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------					
		
		//  dzz																												
		const double cos_real_dzz = cos(P2M*dx);		
		const double sin_real_dzz = sin(P2M*dx);		

		const double sinh_img_dzz = sinh(P2M*dz);			
		const double cosh_img_dzz =	cosh(P2M*dz);		

		// dzz denom																					
		const double denom_dzz_p1 = FMath::pow((sin_real_dzz*cosh_img_dzz),2.0);		
		const double denom_dzz_p2 = FMath::pow((cos_real_dzz*sinh_img_dzz),2.0);		

		FReal T_real = 0;
		FReal T_img	= 0;	

		T_real = ((sin_real_dzz*cos_real_dzz)/(denom_dzz_p1 + denom_dzz_p2)); 								//real part of P1
		T_img =  ((-1*sinh_img_dzz*cosh_img_dzz)/(denom_dzz_p1 + denom_dzz_p2));								//imag part of P1
		
		
		if ((denom_dzz_p1+denom_dzz_p2) < .000000001) {
			 T_real = 0;
			 T_img	= 0;
		}		
		
        const FReal diffp = ((dx * dx) + (dzp * dzp));
        const FReal dzzp = dx + dzp;		

		FReal Dp;		
		FReal Pp;		
		
		if ( rvalsq > diffp / 100){
			Dp = FMath::Exp((-1*diff)/(rvalsq));			
			Pp = (Dp*Dp);
		} else{
			Dp = 0;
			Pp = 0;
		}	

		const FReal SCVP_denom = ((P2M*dx*dx)+(P2M*dzp*dzp));

		FReal scvp_real = (((dx*Dp)+(dx*(-2*Pp)))/SCVP_denom);									//real part of SCVP
		FReal scvp_img = (((-1*dzp*Dp)+(dzp*(2*Pp)))/SCVP_denom);									// imag part of SCVP	

			
		if (SCVP_denom < .000000001) {
			scvp_real = 0;
			scvp_img = 0;
		}	
			
		double P2M_Realdzzp = (P2M*dx);
		double P2M_Imgdzzp = (P2M*dzp);	
		
		//  dzzp																												
		const double cos_real_dzzp = cos(P2M_Realdzzp);		
		const double sin_real_dzzp = sin(P2M_Realdzzp);		

		const double sinh_img_dzzp = sinh(P2M_Imgdzzp);			
		const double cosh_img_dzzp = cosh(P2M_Imgdzzp);		
		
		// dzzp denom																					
		const double denom_dzzp_p1 = FMath::pow((sin_real_dzzp*cosh_img_dzzp),2.0);		
		const double denom_dzzp_p2 = FMath::pow((cos_real_dzzp*sinh_img_dzzp),2.0);	
	
	
		FReal Tp_real = 0;
		FReal Tp_img =  0;
	
		Tp_real = ((sin_real_dzzp*cos_real_dzzp)/(denom_dzzp_p1 + denom_dzzp_p2));								//real part of P3
		Tp_img =  ((-1*cosh_img_dzzp*sinh_img_dzzp)/(denom_dzzp_p1 + denom_dzzp_p2));	 						//imag part of P3
		
		
		if ((denom_dzzp_p1+denom_dzzp_p2) < .000000001) {
			 Tp_real = 0;
			 Tp_img	= 0;
		}	

		Ptot_real = ((T_real + scv_real) - (Tp_real + scvp_real));		//ptot_real = (P1 + P2)_real  -  (P3 + P4)_real
		Ptot_img =  ((T_img + scv_img) - (Tp_img + scvp_img));		//ptot_imag = (P1 + P2)_imag  -  (P3 + P4)_imag

		 return;
    }



// ===========================================            BLOCK            ===========================================



    // evaluate interaction (blockwise)
    template <class ValueClass>
    void evaluateBlock(const ValueClass& xt, const ValueClass& yt, const ValueClass& zt, 
                       const ValueClass& xs, const ValueClass& ys, const ValueClass& zs,
                       ValueClass block[2]) const
    {


		FReal blockReal;
		FReal blockImg;
		evaluate(xt,yt,zt,xs,ys,zs,blockReal,blockImg);
        block[0] = blockReal;
		block[1] = blockImg;

    }



// ===============================================  BLOCK AND DERIVATIVE ===========================================


    // evaluate interaction and derivative (blockwise)
    template <class ValueClass>
    void evaluateBlockAndDerivative(const ValueClass& xt, const ValueClass& yt, const ValueClass& zt,
                                    const ValueClass& xs, const ValueClass& ys, const ValueClass& zs,
                                    ValueClass block[2], ValueClass blockDerivative[6]) const
    {
		std::cout <<" ------ derivative call  ---------"<<std::endl;
		const ValueClass dx = (xt-xs);
        const ValueClass dy = (yt-ys);
        const ValueClass dz = (zt-zs);	
		
        const ValueClass dzp = (zt+zs);
	

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------		
		
        const ValueClass diff = ((dx * dx) + (dz * dz));
        const ValueClass diffp = ((dx * dx) + (dzp * dzp));	
		
		
		ValueClass E1;		
		ValueClass E2;	
		ValueClass EP1;
		ValueClass EP2;		

		if ( rvalsq > (double) diff / 100){
			E1 = FMath::Exp((-1*diff)/(rvalsq));			
			E2 = FMath::Exp((-2*diff)/(rvalsq)); 			
		} else{
			E1 = 0;
			E2 = 0;
		}		
		
		if ( rvalsq > (double) diffp / 100){
			EP1 = FMath::Exp((-1*diffp)/(rvalsq));
			EP2 = FMath::Exp((-2*diffp)/(rvalsq));	
		} else{
			EP1 = 0;
			EP2 = 0;
		}		
		
		const ValueClass dx2 = dx*dx;
		const ValueClass dx4 = dx2*dx2;		
        const ValueClass dz2 = dz*dz;				
        const ValueClass dz4 = dz2*dz2;	
        const ValueClass dzp2 = dzp*dzp;				
        const ValueClass dzp4 = dzp2*dzp2;			


//					A Vector	(SCV)
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------
		const ValueClass A_denom = (P2M * rvalsq * ((dx2 + dz2)*(dx2 + dz2)) );

		const ValueClass A1_real_p1 = E2 *((2*dx2*rvalsq) + (-2*dz2*rvalsq) + (8*dx4) + (8*dx2*dz2));
		const ValueClass A1_real_p2 = E1 *((-2*dx4) +(-1*dx2*rvalsq) + (dz2*rvalsq) + (-2*dx2*dz2) );
	
		const ValueClass A1_real = (( (A1_real_p1) + (A1_real_p2) ) / (A_denom) ); 

		const ValueClass A1_img_p1 = E2 *((-4*dx*dz*rvalsq) + (-8*dx*dz*dx2) + (-8*dx*dz*dz2));
		const ValueClass A1_img_p2 = E1 *((2*dx*dz*rvalsq) + (2*dx*dz*dx2) + (2*dx*dz*dz2) );
		
		const ValueClass A1_img  = (( A1_img_p1 + A1_img_p2 ) / (A_denom) ); 	

		const ValueClass A2_real =  (-1*A1_img);
	
		const ValueClass A2_img_p1 = E2 *((2*dx2*rvalsq) + (-2*dz2*rvalsq) + (-8*dz4) + (-8*dx2*dz2) );
		const ValueClass A2_img_p2 = E1 *((2*dz4) + (-1*dx2*rvalsq) + (dz2*rvalsq) + (2*dx2*dz2));

		const ValueClass A2_img  = (( A2_img_p1 + A2_img_p2 ) / (A_denom) ); 
	

//					B Vector	(SCVP)
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------		
		const ValueClass B_denom = (P2M * rvalsq * ((dx2 + dzp2)*(dx2 + dzp2)) );

		const ValueClass B1_real_p1 = EP2 *((2*dx2*rvalsq) + (-2*dzp2*rvalsq) + (8*dx4) + (8*dx2*dzp2));
		const ValueClass B1_real_p2 = EP1 *((-2*dx4) +(-1*dx2*rvalsq) + (dzp2*rvalsq) + (-2*dx2*dzp2) );           //this blows up. floating point underflow?
	
		const ValueClass B1_real = (( (B1_real_p1) + (B1_real_p2) ) / (B_denom) ); 		

		const ValueClass B1_img_p1 = EP2 *((-4*dx*dzp*rvalsq) + (-8*dx*dzp*dx2) + (-8*dx*dzp*dzp2));
		const ValueClass B1_img_p2 = EP1 *((2*dx*dzp*rvalsq) + (2*dx*dzp*dx2) + (2*dx*dzp*dzp2) );		
		
		const ValueClass B1_img  = (( B1_img_p1 + B1_img_p2 ) / (B_denom) ); 		

		const ValueClass B2_real =  (-1*B1_img);
	
		const ValueClass B2_img_p1 = EP2 *((2*dx2*rvalsq) + (-2*dzp2*rvalsq) + (-8*dzp4) + (-8*dx2*dzp2) );	
		const ValueClass B2_img_p2 = EP1 *((2*dzp4) + (-1*dx2*rvalsq) + (dzp2*rvalsq) + (2*dx2*dzp2));

		const ValueClass B2_img  = (( B2_img_p1 + B2_img_p2 ) / (B_denom) ); 		

		
		const ValueClass X_A = (cosh(P2M*(double)dz)*cos(P2M*(double)dx));  
		const ValueClass X_B = (sinh(P2M*(double)dz)*sin(P2M*(double)dx));
		const ValueClass X_C = (cosh(P2M*(double)dz)*sin(P2M*(double)dx));
		const ValueClass X_D = (sinh(P2M*(double)dz)*cos(P2M*(double)dx));

		const ValueClass X_C2 = (X_C * X_C);
		const ValueClass X_D2 = (X_D * X_D);

		const ValueClass X_denom = ( ( (X_C2 - X_D2)*(X_C2 - X_D2) )  +  (4*X_C2*X_D2) );
		
		const ValueClass X1_real_p1 = (-1*P2M);
		const ValueClass X1_real_t1 = ( (X_A*(X_C2 - X_D2)) - (2*X_B*X_C*X_D) );												
		const ValueClass X1_real = ( (((X1_real_p1)*(X1_real_t1))/X_denom) + X1_real_p1 );													

		const ValueClass X1_img_t1 = ( (X_B*(X_C2 - X_D2)) + (2*X_A*X_C*X_D) );												
		const ValueClass X1_img = ( ((-1*X1_real_p1)*(X1_img_t1))/X_denom ); 												

		const ValueClass X2_real_t1 = ( (X_B*(X_C2 - X_D2)) + (2*X_A*X_C*X_D) );										
		const ValueClass X2_real = ( (((X1_real_p1)*(X1_img_t1))/X_denom));												

		const ValueClass X2_img = X1_real;
		

//					Y Vector	(1 / tan(P2M*(dx+idz))	
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------		
	
		const ValueClass Y_A = (cosh((P2M*(double)dzp))*cos((P2M*(double)dx)));
		const ValueClass Y_B = (sinh((P2M*(double)dzp))*sin((P2M*(double)dx)));
		const ValueClass Y_C = (cosh((P2M*(double)dzp))*sin((P2M*(double)dx)));
		const ValueClass Y_D = (sinh((P2M*(double)dzp))*cos((P2M*(double)dx)));

		const ValueClass Y_C2 = (Y_C * Y_C);
		const ValueClass Y_D2 = (Y_D * Y_D);
		
		const ValueClass Y_denom = (( ( (Y_C2 - Y_D2)*(Y_C2 - Y_D2)) )  +  (4*Y_C2*Y_D2) );
		
		const ValueClass Y1_real_t1 = ( (Y_A*((Y_C2 - Y_D2)) - (2*Y_B*Y_C*Y_D) ));											
		const ValueClass Y1_real = ( (((X1_real_p1)*(Y1_real_t1))/Y_denom) + X1_real_p1 );												

		const ValueClass Y1_img_t1 = ( (Y_B*(Y_C2 - Y_D2)) + (2*Y_A*Y_C*Y_D) );												
		const ValueClass Y1_img = ( ((-1*X1_real_p1)*(Y1_img_t1))/Y_denom ); 																

		const ValueClass Y2_real_t1 = ( (Y_B*(Y_C2 - Y_D2)) + (2*Y_A*Y_C*Y_D) );								
		const ValueClass Y2_real = ( (((X1_real_p1)*(Y1_img_t1))/Y_denom));								

		const ValueClass Y2_img = Y1_real;


//----------------------------------------------------------------------------------------------------------------------------------------------------------------------		
//						FINAL CALCULATIONS
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------		

		const ValueClass Ptot1_real = ( (A1_real + X1_real) - (B1_real + Y1_real) );
		const ValueClass Ptot1_img = ( (A1_img + X1_img) - (B1_img + Y1_img) );		
		const ValueClass Ptot2_real = ( (A2_real + X2_real) - (B2_real + Y2_real) );
		const ValueClass Ptot2_img = ( (A2_img + X2_img) - (B2_img + Y2_img) );	

		FReal blockReal;
		FReal blockImg;		
		evaluate(xt,yt,zt,xs,ys,zs,blockReal,blockImg);
        block[0] = blockReal;
		block[1] = blockImg;

//========================  set the outputs ==================== 
 blockDerivative[0] =  Ptot1_real;
 blockDerivative[1] = 0;
 blockDerivative[2] =  Ptot2_real; 
 blockDerivative[3] =  Ptot1_img;
 blockDerivative[4] = 0;
 blockDerivative[5] =  Ptot2_img;


    }
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------		




    FReal getScaleFactor(const FReal RootCellWidth, const int TreeLevel) const
    {
        const FReal CellWidth(RootCellWidth / FReal(FMath::pow(2, TreeLevel)));
        return getScaleFactor(CellWidth);
    }

    FReal getScaleFactor(const FReal CellWidth) const
    {
        return FReal(4.) / (CellWidth*CellWidth);
    }

    void evaluate(const FPoint<FReal>& pt, const FPoint<FReal>& ps, FReal& blockReal, FReal& blockImg) const{    //updated
        evaluate<FReal>(pt.getX(), pt.getY(), pt.getZ(), ps.getX(), ps.getY(), ps.getZ(), blockReal, blockImg);  //updated
    }
    void evaluateBlock(const FPoint<FReal>& pt, const FPoint<FReal>& ps, FReal* blockReal, FReal* blockImg) const{   //updated
        evaluateBlock<FReal>(pt.getX(), pt.getY(), pt.getZ(), ps.getX(), ps.getY(), ps.getZ(), blockReal, blockImg);  //updated
    }
    void evaluateBlockAndDerivative(const FPoint<FReal>& pt, const FPoint<FReal>& ps,
                                    FReal block[2], FReal blockDerivative[6]) const {  //updated
        evaluateBlockAndDerivative<FReal>(pt.getX(), pt.getY(), pt.getZ(), ps.getX(), ps.getY(), ps.getZ(), block, blockDerivative); //updated
    }
};

/// One over r
template <class FReal>
struct FInterpMatrixKernelR : FInterpAbstractMatrixKernel<FReal>
{
    static const KERNEL_FUNCTION_TYPE Type = HOMOGENEOUS;
    static const unsigned int NCMP = 1; //< number of components
    static const unsigned int NPV  = 1; //< dim of physical values
    static const unsigned int NPOT = 1; //< dim of potentials
    static const unsigned int NRHS = 1; //< dim of mult exp
    static const unsigned int NLHS = 1; //< dim of loc exp

    FInterpMatrixKernelR() {}

    // copy ctor
    FInterpMatrixKernelR(const FInterpMatrixKernelR& /*other*/) {}

    static const char* getID() { return "ONE_OVER_R"; }

    static void printInfo() { std::cout << "K(x,y)=1/r with r=|x-y|" << std::endl; }

    // returns position in reduced storage
    int getPosition(const unsigned int) const
    {return 0;}

    // returns coefficient of mutual interaction
    // 1 for symmetric kernels
    // -1 for antisymmetric kernels
    // Something else if other property of symmetry
    FReal getMutualCoefficient() const{ return FReal(1.); }

    // evaluate interaction
    template <class ValueClass>
    ValueClass evaluate(const ValueClass& xt, const ValueClass& yt, const ValueClass& zt, 
                        const ValueClass& xs, const ValueClass& ys, const ValueClass& zs) const
    {

        const ValueClass diffx = (xt-xs);
        const ValueClass diffy = (yt-ys);
        const ValueClass diffz = (zt-zs);
        return ValueClass(1) / FMath::Sqrt(diffx*diffx + diffy*diffy + diffz*diffz);
    }

    // evaluate interaction (blockwise)
    template <class ValueClass>
    void evaluateBlock(const ValueClass& xt, const ValueClass& yt, const ValueClass& zt, 
                       const ValueClass& xs, const ValueClass& ys, const ValueClass& zs,
                       ValueClass* block) const
    {
        block[0] = this->evaluate(xt,yt,zt,xs,ys,zs);
    }

    // evaluate interaction and derivative (blockwise)
    template <class ValueClass>
    void evaluateBlockAndDerivative(const ValueClass& xt, const ValueClass& yt, const ValueClass& zt,
                                    const ValueClass& xs, const ValueClass& ys, const ValueClass& zs,
                                    ValueClass block[1], ValueClass blockDerivative[3]) const
    {

        const ValueClass diffx = (xt-xs);
        const ValueClass diffy = (yt-ys);
        const ValueClass diffz = (zt-zs);
        const ValueClass one_over_r = ValueClass(1) / FMath::Sqrt(diffx*diffx + diffy*diffy + diffz*diffz);

        const ValueClass one_over_r3 = one_over_r*one_over_r*one_over_r;

        block[0] = one_over_r;

        blockDerivative[0] = - one_over_r3 * diffx;
        blockDerivative[1] = - one_over_r3 * diffy;
        blockDerivative[2] = - one_over_r3 * diffz;
    }

    FReal getScaleFactor(const FReal RootCellWidth, const int TreeLevel) const
    {
        const FReal CellWidth(RootCellWidth / FReal(FMath::pow(2, TreeLevel)));
        return getScaleFactor(CellWidth);
    }

    FReal getScaleFactor(const FReal CellWidth) const
    {
        return FReal(2.) / CellWidth;
    }

    FReal evaluate(const FPoint<FReal>& pt, const FPoint<FReal>& ps) const {
        return evaluate<FReal>(pt.getX(), pt.getY(), pt.getZ(), ps.getX(), ps.getY(), ps.getZ());
    }
    void evaluateBlock(const FPoint<FReal>& pt, const FPoint<FReal>& ps, FReal* block) const{
        evaluateBlock<FReal>(pt.getX(), pt.getY(), pt.getZ(), ps.getX(), ps.getY(), ps.getZ(), block);
    }
    void evaluateBlockAndDerivative(const FPoint<FReal>& pt, const FPoint<FReal>& ps,
                                    FReal block[1], FReal blockDerivative[3]) const {
        evaluateBlockAndDerivative<FReal>(pt.getX(), pt.getY(), pt.getZ(), ps.getX(), ps.getY(), ps.getZ(), block, blockDerivative);
    }
};

/// One over r when the box size is rescaled to 1
template <class FReal>
struct FInterpMatrixKernelRH :FInterpMatrixKernelR<FReal>{
    static const KERNEL_FUNCTION_TYPE Type = HOMOGENEOUS;
    static const unsigned int NCMP = 1; //< number of components
    static const unsigned int NPV  = 1; //< dim of physical values
    static const unsigned int NPOT = 1; //< dim of potentials
    static const unsigned int NRHS = 1; //< dim of mult exp
    static const unsigned int NLHS = 1; //< dim of loc exp
    FReal LX,LY,LZ ;

    FInterpMatrixKernelRH() 
    : LX(1.0),LY(1.0),LZ(1.0)
    { }

    // copy ctor
    FInterpMatrixKernelRH(const FInterpMatrixKernelRH& other)
    : FInterpMatrixKernelR<FReal>(other), LX(other.LX), LY(other.LY), LZ(other.LZ)
    {}

    static const char* getID() { return "ONE_OVER_RH"; }

    static void printInfo() { std::cout << "K(x,y)=1/rh with rh=sqrt(L_i*(x_i-y_i)^2)" << std::endl; }

    // evaluate interaction
    template <class ValueClass>
    ValueClass evaluate(const ValueClass& xt, const ValueClass& yt, const ValueClass& zt, 
                        const ValueClass& xs, const ValueClass& ys, const ValueClass& zs) const
    {
        const ValueClass diffx = (xt-xs);
        const ValueClass diffy = (yt-ys);
        const ValueClass diffz = (zt-zs);
        return ValueClass(1) / FMath::Sqrt(ValueClass(LX)*diffx*diffx +
                                       ValueClass(LY)*diffy*diffy +
                                       ValueClass(LZ)*diffz*diffz);
    }
    void setCoeff(const FReal& a,  const FReal& b, const FReal& c)
    {LX= a*a ; LY = b*b ; LZ = c *c;}
    // returns position in reduced storage
    int getPosition(const unsigned int) const
    {return 0;}
    // returns coefficient of mutual interaction
    // 1 for symmetric kernels
    // -1 for antisymmetric kernels
    // Something else if other property of symmetry
    FReal getMutualCoefficient() const{ return FReal(1.); }

    template <class ValueClass>
    void evaluateBlock(const ValueClass& xt, const ValueClass& yt, const ValueClass& zt, 
                       const ValueClass& xs, const ValueClass& ys, const ValueClass& zs,
                       ValueClass* block) const
    {
        block[0]=this->evaluate(xt,yt,zt,xs,ys,zs);
    }

    // evaluate interaction and derivative (blockwise)
    template <class ValueClass>
    void evaluateBlockAndDerivative(const ValueClass& xt, const ValueClass& yt, const ValueClass& zt,
                                    const ValueClass& xs, const ValueClass& ys, const ValueClass& zs,
                                    ValueClass block[1], ValueClass blockDerivative[3]) const
    {
        const ValueClass diffx = (xt-xs);
        const ValueClass diffy = (yt-ys);
        const ValueClass diffz = (zt-zs);
        const ValueClass one_over_rL = ValueClass(1) / (ValueClass(LX)*diffx*diffx +
                                                          ValueClass(LY)*diffy*diffy +
                                                          ValueClass(LZ)*diffz*diffz);
        const ValueClass one_over_rL3 = one_over_rL*one_over_rL*one_over_rL;

        block[0] = one_over_rL;

        blockDerivative[0] = ValueClass(LX) * one_over_rL3 * diffx;
        blockDerivative[1] = ValueClass(LY)* one_over_rL3 * diffy;
        blockDerivative[2] = ValueClass(LZ)* one_over_rL3 * diffz;

    }

    FReal getScaleFactor(const FReal RootCellWidth, const int TreeLevel) const
    {
        const FReal CellWidth(RootCellWidth / FReal(FMath::pow(2, TreeLevel)));
        return getScaleFactor(CellWidth);
    }

    FReal getScaleFactor(const FReal CellWidth) const
    {
        return FReal(2.) / CellWidth;
    }

    FReal evaluate(const FPoint<FReal>& pt, const FPoint<FReal>& ps) const{
        return evaluate<FReal>(pt.getX(), pt.getY(), pt.getZ(), ps.getX(), ps.getY(), ps.getZ());
    }
    void evaluateBlock(const FPoint<FReal>& pt, const FPoint<FReal>& ps, FReal* block) const{
        evaluateBlock<FReal>(pt.getX(), pt.getY(), pt.getZ(), ps.getX(), ps.getY(), ps.getZ(), block);
    }
    void evaluateBlockAndDerivative(const FPoint<FReal>& pt, const FPoint<FReal>& ps,
                                    FReal block[1], FReal blockDerivative[3]) const {
        evaluateBlockAndDerivative<FReal>(pt.getX(), pt.getY(), pt.getZ(), ps.getX(), ps.getY(), ps.getZ(), block, blockDerivative);
    }
};


/// One over r^2
template <class FReal>
struct FInterpMatrixKernelRR : FInterpAbstractMatrixKernel<FReal>
{
    static const KERNEL_FUNCTION_TYPE Type = HOMOGENEOUS;
    static const unsigned int NCMP = 1; //< number of components
    static const unsigned int NPV  = 1; //< dim of physical values
    static const unsigned int NPOT = 1; //< dim of potentials
    static const unsigned int NRHS = 1; //< dim of mult exp
    static const unsigned int NLHS = 1; //< dim of loc exp

    FInterpMatrixKernelRR() {}

    // copy ctor
    FInterpMatrixKernelRR(const FInterpMatrixKernelRR& /*other*/) {}

    static const char* getID() { return "ONE_OVER_R_SQUARED"; }

    static void printInfo() { std::cout << "K(x,y)=1/r^2 with r=|x-y|" << std::endl; }

    // returns position in reduced storage
    int getPosition(const unsigned int) const
    {return 0;}

    // returns coefficient of mutual interaction
    // 1 for symmetric kernels
    // -1 for antisymmetric kernels
    // Something else if other property of symmetry
    FReal getMutualCoefficient() const{ return FReal(1.); }

    // evaluate interaction
    template <class ValueClass>
    ValueClass evaluate(const ValueClass& xt, const ValueClass& yt, const ValueClass& zt, 
                        const ValueClass& xs, const ValueClass& ys, const ValueClass& zs) const
    {
        const ValueClass diffx = (xt-xs);
        const ValueClass diffy = (yt-ys);
        const ValueClass diffz = (zt-zs);
        return ValueClass(1) / FReal(diffx*diffx+diffy*diffy+diffz*diffz);
    }

    // evaluate interaction (blockwise)
    template <class ValueClass>
    void evaluateBlock(const ValueClass& xt, const ValueClass& yt, const ValueClass& zt, 
                       const ValueClass& xs, const ValueClass& ys, const ValueClass& zs,
                       ValueClass* block) const
    {
        block[0]=this->evaluate(xt,yt,zt,xs,ys,zs);
    }

    // evaluate interaction and derivative (blockwise)
    template <class ValueClass>
    void evaluateBlockAndDerivative(const ValueClass& xt, const ValueClass& yt, const ValueClass& zt,
                                    const ValueClass& xs, const ValueClass& ys, const ValueClass& zs,
                                    ValueClass block[1], ValueClass blockDerivative[3]) const
    {
        const ValueClass diffx = (xt-xs);
        const ValueClass diffy = (yt-ys);
        const ValueClass diffz = (zt-zs);
        const ValueClass r2 = (diffx*diffx+diffy*diffy+diffz*diffz);
        const ValueClass one_over_r2 = ValueClass(1) / (r2);
        const ValueClass one_over_r4 = one_over_r2*one_over_r2;

        block[0] = one_over_r2;

        const ValueClass coef = ValueClass(-2.) * one_over_r4;
        blockDerivative[0] = coef * diffx;
        blockDerivative[1] = coef * diffy;
        blockDerivative[2] = coef * diffz;

    }

    FReal getScaleFactor(const FReal RootCellWidth, const int TreeLevel) const
    {
        const FReal CellWidth(RootCellWidth / FReal(FMath::pow(2, TreeLevel)));
        return getScaleFactor(CellWidth);
    }

    FReal getScaleFactor(const FReal CellWidth) const
    {
        return FReal(4.) / (CellWidth*CellWidth);
    }

    FReal evaluate(const FPoint<FReal>& pt, const FPoint<FReal>& ps) const{
        return evaluate<FReal>(pt.getX(), pt.getY(), pt.getZ(), ps.getX(), ps.getY(), ps.getZ());
    }
    void evaluateBlock(const FPoint<FReal>& pt, const FPoint<FReal>& ps, FReal* block) const{
        evaluateBlock<FReal>(pt.getX(), pt.getY(), pt.getZ(), ps.getX(), ps.getY(), ps.getZ(), block);
    }
    void evaluateBlockAndDerivative(const FPoint<FReal>& pt, const FPoint<FReal>& ps,
                                    FReal block[1], FReal blockDerivative[3]) const {
        evaluateBlockAndDerivative<FReal>(pt.getX(), pt.getY(), pt.getZ(), ps.getX(), ps.getY(), ps.getZ(), block, blockDerivative);
    }
};



/// One over r^12 - One over r^6
template <class FReal>
struct FInterpMatrixKernelLJ : FInterpAbstractMatrixKernel<FReal>
{
    static const KERNEL_FUNCTION_TYPE Type = NON_HOMOGENEOUS;
    static const unsigned int NCMP = 1; //< number of components
    static const unsigned int NPV  = 1; //< dim of physical values
    static const unsigned int NPOT = 1; //< dim of potentials
    static const unsigned int NRHS = 1; //< dim of mult exp
    static const unsigned int NLHS = 1; //< dim of loc exp

    FInterpMatrixKernelLJ() {}

    // copy ctor
    FInterpMatrixKernelLJ(const FInterpMatrixKernelLJ& /*other*/) {}

    static const char* getID() { return "LENNARD_JONES_POTENTIAL"; }

    static void printInfo() { std::cout << "K(x,y)=1/r with r=|x-y|" << std::endl; }

    // returns position in reduced storage
    int getPosition(const unsigned int) const
    {return 0;}

    // returns coefficient of mutual interaction
    // 1 for symmetric kernels
    // -1 for antisymmetric kernels
    // somethings else if other property of symmetry
    FReal getMutualCoefficient() const{ return FReal(1.); }

    // evaluate interaction
    template <class ValueClass>
    ValueClass evaluate(const ValueClass& xt, const ValueClass& yt, const ValueClass& zt, 
                        const ValueClass& xs, const ValueClass& ys, const ValueClass& zs) const
    {
        const ValueClass diffx = (xt-xs);
        const ValueClass diffy = (yt-ys);
        const ValueClass diffz = (zt-zs);
        const ValueClass r = FMath::Sqrt(diffx*diffx+diffy*diffy+diffz*diffz);
        const ValueClass r3 = r*r*r;
        const ValueClass one_over_r6 = ValueClass(1) / (r3*r3);
        //return one_over_r6 * one_over_r6;
        //return one_over_r6;
        return one_over_r6 * one_over_r6 - one_over_r6;
    }

    // evaluate interaction (blockwise)
    template <class ValueClass>
    void evaluateBlock(const ValueClass& xt, const ValueClass& yt, const ValueClass& zt, 
                       const ValueClass& xs, const ValueClass& ys, const ValueClass& zs,
                       ValueClass* block) const
    {
        block[0]=this->evaluate(xt,yt,zt,xs,ys,zs);
    }

    // evaluate interaction and derivative (blockwise)
    template <class ValueClass>
    void evaluateBlockAndDerivative(const ValueClass& xt, const ValueClass& yt, const ValueClass& zt,
                                    const ValueClass& xs, const ValueClass& ys, const ValueClass& zs,
                                    ValueClass block[1], ValueClass blockDerivative[3]) const
    {
        const ValueClass diffx = (xt-xs);
        const ValueClass diffy = (yt-ys);
        const ValueClass diffz = (zt-zs);
        const ValueClass r = FMath::Sqrt(diffx*diffx+diffy*diffy+diffz*diffz);
        const ValueClass r2 = r*r;
        const ValueClass r3 = r2*r;
        const ValueClass one_over_r6 = ValueClass(1) / (r3*r3);
        const ValueClass one_over_r8 = one_over_r6 / (r2);

        block[0] = one_over_r6 * one_over_r6 - one_over_r6;

        const FReal coef = ValueClass(12.0)*one_over_r6*one_over_r8 - ValueClass(6.0)*one_over_r8;
        blockDerivative[0]= coef * diffx;
        blockDerivative[1]= coef * diffy;
        blockDerivative[2]= coef * diffz;

    }

    FReal getScaleFactor(const FReal, const int) const
    {
        // return 1 because non homogeneous kernel functions cannot be scaled!!!
        return FReal(1.0);
    }

    FReal getScaleFactor(const FReal) const
    {
        // return 1 because non homogeneous kernel functions cannot be scaled!!!
        return FReal(1.0);
    }



    FReal evaluate(const FPoint<FReal>& pt, const FPoint<FReal>& ps) const{
        return evaluate<FReal>(pt.getX(), pt.getY(), pt.getZ(), ps.getX(), ps.getY(), ps.getZ());
    }
    void evaluateBlock(const FPoint<FReal>& pt, const FPoint<FReal>& ps, FReal* block) const{
        evaluateBlock<FReal>(pt.getX(), pt.getY(), pt.getZ(), ps.getX(), ps.getY(), ps.getZ(), block);
    }
    void evaluateBlockAndDerivative(const FPoint<FReal>& pt, const FPoint<FReal>& ps,
                                    FReal block[1], FReal blockDerivative[3]) const {
        evaluateBlockAndDerivative<FReal>(pt.getX(), pt.getY(), pt.getZ(), ps.getX(), ps.getY(), ps.getZ(), block, blockDerivative);
    }
};


/// One over (a+r^2)
template <class FReal>
struct FInterpMatrixKernelAPLUSRR : FInterpAbstractMatrixKernel<FReal>
{
    static const KERNEL_FUNCTION_TYPE Type = NON_HOMOGENEOUS;
    static const unsigned int NCMP = 1; //< number of components
    static const unsigned int NPV  = 1; //< dim of physical values
    static const unsigned int NPOT = 1; //< dim of potentials
    static const unsigned int NRHS = 1; //< dim of mult exp
    static const unsigned int NLHS = 1; //< dim of loc exp

    const FReal CoreWidth;

    FInterpMatrixKernelAPLUSRR(const FReal inCoreWidth = .25)
    : CoreWidth(inCoreWidth)
    {}

    // copy ctor
    FInterpMatrixKernelAPLUSRR(const FInterpMatrixKernelAPLUSRR& other)
    : CoreWidth(other.CoreWidth)
    {}

    static const char* getID() { return "ONE_OVER_A_PLUS_RR"; }

    static void printInfo() { std::cout << "K(x,y)=1/r with r=|x-y|" << std::endl; }

    // returns position in reduced storage
    int getPosition(const unsigned int) const
    {return 0;}

    // returns coefficient of mutual interaction
    // 1 for symmetric kernels
    // -1 for antisymmetric kernels
    // something else if other property of symmetry
    FReal getMutualCoefficient() const{ return FReal(1.); }

    // evaluate interaction
    template <class ValueClass>
    ValueClass evaluate(const ValueClass& xt, const ValueClass& yt, const ValueClass& zt,
                        const ValueClass& xs, const ValueClass& ys, const ValueClass& zs) const
    {
        const ValueClass diffx = (xt-xs);
        const ValueClass diffy = (yt-ys);
        const ValueClass diffz = (zt-zs);
        const ValueClass r2 = (diffx*diffx+diffy*diffy+diffz*diffz);
        return ValueClass(1) / (r2 + ValueClass(CoreWidth));
    }

    // evaluate interaction (blockwise)
    template <class ValueClass>
    void evaluateBlock(const ValueClass& xt, const ValueClass& yt, const ValueClass& zt, 
                       const ValueClass& xs, const ValueClass& ys, const ValueClass& zs,
                       ValueClass* block) const
    {
        block[0]=this->evaluate(xt,yt,zt,xs,ys,zs);
    }

    // evaluate interaction and derivative (blockwise)
    template <class ValueClass>
    void evaluateBlockAndDerivative(const ValueClass& xt, const ValueClass& yt, const ValueClass& zt,
                                    const ValueClass& xs, const ValueClass& ys, const ValueClass& zs,
                                    ValueClass block[1], ValueClass blockDerivative[3]) const
    {
        const ValueClass diffx = (xt-xs);
        const ValueClass diffy = (yt-ys);
        const ValueClass diffz = (zt-zs);
        const ValueClass r2 = (diffx*diffx+diffy*diffy+diffz*diffz);
        const ValueClass one_over_a_plus_r2 = ValueClass(1) / (r2 + ValueClass(CoreWidth));
        const ValueClass one_over_a_plus_r2_squared = one_over_a_plus_r2*one_over_a_plus_r2;

        block[0] = one_over_a_plus_r2;

        // TODO Fix derivative
        const ValueClass coef = ValueClass(-2.) * one_over_a_plus_r2_squared;
        blockDerivative[0] = coef * diffx;
        blockDerivative[1] = coef * diffy;
        blockDerivative[2] = coef * diffz;

    }

    FReal getScaleFactor(const FReal, const int) const
    {
        // return 1 because non homogeneous kernel functions cannot be scaled!!!
        return FReal(1.0);
    }

    FReal getScaleFactor(const FReal) const
    {
        // return 1 because non homogeneous kernel functions cannot be scaled!!!
        return FReal(1.0);    
    }
	
	    FReal evaluate(const FPoint<FReal>& pt, const FPoint<FReal>& ps) const{
        return evaluate<FReal>(pt.getX(), pt.getY(), pt.getZ(), ps.getX(), ps.getY(), ps.getZ());
    }
    void evaluateBlock(const FPoint<FReal>& pt, const FPoint<FReal>& ps, FReal* block) const{
        evaluateBlock<FReal>(pt.getX(), pt.getY(), pt.getZ(), ps.getX(), ps.getY(), ps.getZ(), block);
    }
    void evaluateBlockAndDerivative(const FPoint<FReal>& pt, const FPoint<FReal>& ps,
                                    FReal block[1], FReal blockDerivative[3]) const {
        evaluateBlockAndDerivative<FReal>(pt.getX(), pt.getY(), pt.getZ(), ps.getX(), ps.getY(), ps.getZ(), block, blockDerivative);
    }
	
};



/*!  Functor which provides the interface to assemble a matrix based on the
  number of rows and cols and on the coordinates s and t and the type of the
  generating matrix-kernel function.
*/
template <class FReal, typename MatrixKernelClass>
class EntryComputer
{
    const MatrixKernelClass *const MatrixKernel;

    const unsigned int nt, ns;
    const FPoint<FReal> *const pt, *const ps;
    const FReal *const weights;

public:
    explicit EntryComputer(const MatrixKernelClass *const inMatrixKernel,
                           const unsigned int _nt, const FPoint<FReal> *const _ps,
                           const unsigned int _ns, const FPoint<FReal> *const _pt,
                           const FReal *const _weights = NULL)
        : MatrixKernel(inMatrixKernel),	nt(_nt), ns(_ns), pt(_pt), ps(_ps), weights(_weights) {}
	
	
	 void operator()(const unsigned int tbeg, const unsigned int tend,
                    const unsigned int sbeg, const unsigned int send,
                    FReal *const data, FReal *const data_i) const
    {	

		FReal pt_x, pt_y, pt_z, ps_x, ps_y, ps_z;

		FReal Ptot_real;     //updated
		FReal Ptot_img;		//updated
        unsigned int idx = 0;

		
		if (weights) {

            for (unsigned int j=tbeg; j<tend; ++j)
                for (unsigned int i=sbeg; i<send; ++i)
					{	

					MatrixKernel->evaluate(pt[i], ps[j], Ptot_real, Ptot_img);				//updated	

					
					data[idx] = (weights[i] * weights[j] * Ptot_real);					//updated		
					data_i[idx++] = (weights[i] * weights[j] * Ptot_img);	//updated
				
					}
        } else {
		
            for (unsigned int j=tbeg; j<tend; ++j)
                for (unsigned int i=sbeg; i<send; ++i)
					{	
					MatrixKernel->evaluate(pt[i], ps[j], Ptot_real, Ptot_img);  //updated

                    data[idx] = Ptot_real;  //updated
					data_i[idx++] = Ptot_img;	//updated

					}
        }
	
	}
};





#endif // FINTERPMATRIXKERNEL_HPP

// [--END--]
