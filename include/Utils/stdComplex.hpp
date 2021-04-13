// ===================================================================================
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
#ifndef STDCOMPLEXE_HPP
#define STDCOMPLEXE_HPP

#include <complex>
template<typename T>
//using stdComplex = std::complex<T> ;
class stdComplex : public std::complex<T>  {
public:
    using Base = std::complex<T> ; 
    using Base::Base ; 

    /** Mul other and another and add the result to current complexe 
    
    
    */
    /*!
     * \brief addMul perform  z += other*another without using the  function __muldc3 from libgcc. ;
     * 
     * "Without -ffast-math or other, complex multiplication yields a call to the function __muldc3 from libgcc. "
     * @see https://stackoverflow.com/questions/42659668/stdcomplex-multiplication-is-extremely-slow 
     * 
     * \param other 
     * \param another
     */
    void addMul(const stdComplex<T>& other, const stdComplex<T>& another){
     //   this->complex[0] += (other.complex[0] * another.complex[0]) - (other.complex[1] * another.complex[1]);
     //     this->complex[1] += (other.complex[0] * another.complex[1]) + (other.complex[1] * another.complex[0]);
        
        T realPart = this->real() + (other.real() * another.real()) - (other.imag() * another.imag());
        T imagPart = this->imag() + (other.real() * another.imag()) + (other.imag() * another.real());
        this->real(realPart) ; 
        this->imag(imagPart) ; 
    }
} ; 
#endif //STDCOMPLEXE_HPP


