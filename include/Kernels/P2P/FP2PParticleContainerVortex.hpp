// See LICENCE file at project root
#ifndef FP2PPARTICLECONTAINERVORTEX_HPP
#define FP2PPARTICLECONTAINERVORTEX_HPP

#include "Components/FBasicParticleContainer_i.hpp"
#include "Utils/FComplex.hpp"

template<class FReal, int NRHS = 1, int NLHS = 1, int NVALS = 1, class... OtherAttrs>
class FP2PParticleContainerVortex : public FBasicParticleContainer_i<FReal, NVALS*(NRHS+8*NLHS), FReal, 3, FAlignedAllocator<FP2PDefaultAlignement,char>, OtherAttrs...> {
    using Parent = FBasicParticleContainer_i<FReal, NVALS*(NRHS+8*NLHS), FReal,
                                           3, FAlignedAllocator<FP2PDefaultAlignement,char>, OtherAttrs...>;

public:
    static const int NbAttributes = NVALS*(NRHS+8*NLHS);
    typedef FReal AttributesClass;

// Total size = NVALS*(NRHS+8*NLHS) = 9 bytes
// 0 -- getPhysicalValues
// 1 -- getPotentials real
// 2 -- getForcesX real
// 3 -- getForcesY real 
// 4 -- getForcesZ real
// 5 -- getPotentials imag
// 6 -- getForcesX imag 
// 7 -- getForcesY imag 
// 8 -- getForcesZ imag 


// NRHS = 1
// NLHS = 1
// NVALS = 1
// idxRhs = 0
// idxVals = 0
// (0+idxRhs)*NVALS+idxVals = ((0+0)*1) + 0 = --> 0
// -------------------------------------------------------------------------------------------------------
    FReal* getPhysicalValues(const int idxVals = 0, const int idxRhs = 0){
      return Parent::getAttribute((0+idxRhs)*NVALS+idxVals);//0
    }

    const FReal* getPhysicalValues(const int idxVals = 0, const int idxRhs = 0) const {
        return Parent::getAttribute((0+idxRhs)*NVALS+idxVals);//0
    }

    FReal* getPhysicalValuesArray(const int idxVals = 0, const int idxRhs = 0){
        return Parent::getRawData() + ((0+idxRhs)*NVALS+idxVals)*Parent::getLeadingRawData();//0
    }

    const FReal* getPhysicalValuesArray(const int idxVals = 0, const int idxRhs = 0) const {
        return Parent::getRawData() + ((0+idxRhs)*NVALS+idxVals)*Parent::getLeadingRawData();//0
    }


// -------------------------------------------------------------------------------------------------------
    FSize getLeadingDimension() const {
        return Parent::getLeadingRawData();
    }

// NRHS = 1
// NLHS = 1
// NVALS = 1
// idxLhs = 0
// idxVals = 0
// (NRHS+idxLhs)*NVALS+idxVals = ((1+0)*1) + 0 = --> 1
// (NRHS+4*NLHS+idxLhs)*NVALS+idxVals = ((1+4+0)*1) + 0 = --> 5
// -------------------------------------------------------------------------------------------------------
    FReal* getPotentials(const int idxVals = 0, const int idxLhs = 0){
        return Parent::getAttribute((NRHS+idxLhs)*NVALS+idxVals);//1
    }

    const FReal* getPotentials(const int idxVals = 0, const int idxLhs = 0) const {	
        return Parent::getAttribute((NRHS+idxLhs)*NVALS+idxVals);//1
    }

    FReal* getPotentials_real(const int idxVals = 0, const int idxLhs = 0)  {										// getPotentials_real -- 1 
        return Parent::getAttribute((NRHS+idxLhs)*NVALS+idxVals);//1
    }

    FReal* getPotentials_imag(const int idxVals = 0, const int idxLhs = 0)  {										// getPotentials_imag -- 5 
        return Parent::getAttribute((NRHS+4*NLHS+idxLhs)*NVALS+idxVals);//5
    }	

    FReal* getPotentialsArray(const int idxVals = 0, const int idxLhs = 0){
        return Parent::getRawData() + ((NRHS+idxLhs)*NVALS+idxVals)*Parent::getLeadingRawData();//1
    }

    const FReal* getPotentialsArray(const int idxVals = 0, const int idxLhs = 0) const {
        return Parent::getRawData() + ((NRHS+idxLhs)*NVALS+idxVals)*Parent::getLeadingRawData();//1
    }

// NRHS = 1
// NLHS = 1
// NVALS = 1
// idxLhs = 0
// idxVals = 0
// ((NRHS+NLHS+idxLhs)*NVALS+idxVals) = ((1+(1*1)+0)*1) + 0 = --> 2    
// ((NRHS+5*NLHS+idxLhs)*NVALS+idxVals) = ((1+(5*1)+0)*1) + 0 = --> 6    (NRHS+3*NLHS+idxLhs)*NVALS+idxVals = ((1+(3*1)+0)*1) + 0 = --> 4   
// -------------------------------------------------------------------------------------------------------
    FReal* getForcesX(const int idxVals = 0, const int idxLhs = 0){
			std::cout <<" --------------------------------------------   THIS DOES NOT PRINT  (DNP) --------------------------------------------------"  << std::endl; 				
        return Parent::getAttribute((NRHS+NLHS+idxLhs)*NVALS+idxVals);//2    ((NRHS+2*NLHS+idxLhs)*NVALS+idxVals) --> 3
    }

    FReal* getForcesX_real(const int idxVals = 0, const int idxLhs = 0)  {										// getForcesX_real -- 2 
		
        return Parent::getAttribute((NRHS+NLHS+idxLhs)*NVALS+idxVals);//2
    }

    FReal* getForcesX_imag(const int idxVals = 0, const int idxLhs = 0)  {										// getForcesX_imag -- 6	
        return Parent::getAttribute((NRHS+5*NLHS+idxLhs)*NVALS+idxVals);//6  
    }


    const FReal* getForcesX(const int idxVals = 0, const int idxLhs = 0) const {
			std::cout <<" --------------------------------------------   THIS DOES NOT PRINT  (DNP) --------------------------------------------------"  << std::endl; 				
        return Parent::getAttribute((NRHS+NLHS+idxLhs)*NVALS+idxVals);//2
    }

    FReal* getForcesXArray(const int idxVals = 0, const int idxLhs = 0){
			std::cout <<" --------------------------------------------   THIS DOES NOT PRINT  (DNP) --------------------------------------------------"  << std::endl; 				
        return Parent::getRawData() + ((NRHS+NLHS+idxLhs)*NVALS+idxVals)*Parent::getLeadingRawData();//2
    }

    const FReal* getForcesXArray(const int idxVals = 0, const int idxLhs = 0) const {
			std::cout <<" --------------------------------------------   THIS DOES NOT PRINT  (DNP) --------------------------------------------------"  << std::endl; 				
        return Parent::getRawData() + ((NRHS+NLHS+idxLhs)*NVALS+idxVals)*Parent::getLeadingRawData();//2
    }

// NRHS = 1
// NLHS = 1
// NVALS = 1
// idxLhs = 0
// idxVals = 0
// (NRHS+2*NLHS+idxLhs)*NVALS+idxVals = ((1+(2*1)+0)*1) + 0 = --> 3 
// (NRHS+6*NLHS+idxLhs)*NVALS+idxVals = ((1+(6*1)+0)*1) + 0 = --> 7
// -------------------------------------------------------------------------------------------------------
    FReal* getForcesY(const int idxVals = 0, const int idxLhs = 0){
        return Parent::getAttribute((NRHS+2*NLHS+idxLhs)*NVALS+idxVals);//3
    }

    FReal* getForcesY_real(const int idxVals = 0, const int idxLhs = 0)  {										// getForcesY_real -- 3 
		
        return Parent::getAttribute((NRHS+2*NLHS+idxLhs)*NVALS+idxVals);//3
    }
	
    FReal* getForcesY_imag(const int idxVals = 0, const int idxLhs = 0)  {										// getForcesY_imag -- 7 
        return Parent::getAttribute((NRHS+6*NLHS+idxLhs)*NVALS+idxVals);//7
    }

    const FReal* getForcesY(const int idxVals = 0, const int idxLhs = 0) const {
        return Parent::getAttribute((NRHS+2*NLHS+idxLhs)*NVALS+idxVals);//3
    }
	
    FReal* getForcesYArray(const int idxVals = 0, const int idxLhs = 0){
        return Parent::getRawData() + ((NRHS+2*NLHS+idxLhs)*NVALS+idxVals)*Parent::getLeadingRawData();//3
    }

    const FReal* getForcesYArray(const int idxVals = 0, const int idxLhs = 0) const {
        return Parent::getRawData() + ((NRHS+2*NLHS+idxLhs)*NVALS+idxVals)*Parent::getLeadingRawData();//3
    }

// NRHS = 1
// NLHS = 1
// NVALS = 1
// idxLhs = 0
// idxVals = 0
// (NRHS+4*NLHS+idxLhs)*NVALS+idxVals = ((1+(3*1)+0)*1) + 0 = --> 4
// (NRHS+3*NLHS+idxLhs)*NVALS+idxVals = ((1+(7*1)+0)*1) + 0 = --> 8
// -------------------------------------------------------------------------------------------------------
    FReal* getForcesZ(const int idxVals = 0, const int idxLhs = 0){
			std::cout <<" --------------------------------------------   THIS DOES NOT PRINT  (DNP) --------------------------------------------------"  << std::endl; 				
        return Parent::getAttribute((NRHS+3*NLHS+idxLhs)*NVALS+idxVals);//4
    }

    const FReal* getForcesZ(const int idxVals = 0, const int idxLhs = 0) const {	
			std::cout <<" --------------------------------------------   THIS DOES NOT PRINT  (DNP) --------------------------------------------------"  << std::endl; 			
        return Parent::getAttribute((NRHS+3*NLHS+idxLhs)*NVALS+idxVals);//4
    }

    FReal* getForcesZ_real(const int idxVals = 0, const int idxLhs = 0)  {										// getForcesZ_real -- 4
		
        return Parent::getAttribute((NRHS+3*NLHS+idxLhs)*NVALS+idxVals);//4
    }

    FReal* getForcesZ_imag(const int idxVals = 0, const int idxLhs = 0)  {										// getForcesZ_imag -- 8 
        return Parent::getAttribute((NRHS+7*NLHS+idxLhs)*NVALS+idxVals);//8
    }
	
    FReal* getForcesZArray(const int idxVals = 0, const int idxLhs = 0){
        return Parent::getRawData() + ((NRHS+3*NLHS+idxLhs)*NVALS+idxVals)*Parent::getLeadingRawData();//4
    }

    const FReal* getForcesZArray(const int idxVals = 0, const int idxLhs = 0) const {
        return Parent::getRawData() + ((NRHS+3*NLHS+idxLhs)*NVALS+idxVals)*Parent::getLeadingRawData();//4
    }


// -------------------------------------------------------------------------------------------------------
    void resetForcesAndPotential(){
        for(int idx = 0 ; idx < 8*NLHS*NVALS ; ++idx){
            Parent::resetToInitialState(idx + NRHS*NVALS);
        }
    }

    int getNVALS() const {
        return NVALS;
    }

};

#endif // FP2PPARTICLECONTAINER_HPP
