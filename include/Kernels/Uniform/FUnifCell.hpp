// See LICENCE file at project root
// Keep in private GIT

#ifndef FUNIFCELL_HPP
#define FUNIFCELL_HPP

#include "Utils/stdComplex.hpp"

#include "FUnifTensor.hpp"
#include "Components/FBasicCell.hpp"
#include "Extensions/FExtendCellType.hpp"


/**
 * @author Pierre Blanchard (pierre.blanchard@inria.fr)
 *
 * Please read the license
 *
 * This class defines a cell used in the Lagrange based FMM.
 *
 * PB: This class also contains the storage and accessors for the transformed
 * expansion (in Fourier space, i.e. complex valued).
 *
 * @tparam NVALS is the number of right hand side.
 */
template < class FReal, int ORDER, int NRHS = 1, int NLHS = 1, int NVALS = 1>
class FUnifCell : public FBasicCell, public FAbstractSendable
{
  static const int VectorSize            = TensorTraits<ORDER>::nnodes;
  static const int TransformedVectorSize = (2*ORDER-1)*(2*ORDER-1)*(2*ORDER-1);

public:

  template<class Tag, std::size_t N>
  struct exp_impl {
    /// Multipole expansion in Real space
    FReal exp[N * NVALS * VectorSize]; //< Multipole expansion
    /// Multipole expansion in Fourier space
    stdComplex<FReal> transformed_exp[N * NVALS * TransformedVectorSize];

    
    const FReal* get(const int inRhs) const
    { return this->exp + inRhs*VectorSize; }
      
    FReal* get(const int inRhs)
    { return this->exp + inRhs*VectorSize; }
      
    const stdComplex<FReal>* getTransformed(const int inRhs) const
    { return this->transformed_exp + inRhs*TransformedVectorSize;
    }
    stdComplex<FReal>* getTransformed(const int inRhs)
    { return this->transformed_exp + inRhs*TransformedVectorSize; }

    constexpr int getVectorSize() const {
      return VectorSize;
    }

    void reset(){
      std::memset(this->exp, 0, sizeof(FReal) * N * NVALS * VectorSize);
      std::memset(this->transformed_exp, 0,
              sizeof(stdComplex<FReal>) * N * NVALS * TransformedVectorSize);
    }
    // to extend FAbstractSendable
    template <class BufferWriterClass>
    void serialize(BufferWriterClass& buffer) const{
      buffer.write(this->exp, VectorSize*NVALS*N);
      buffer.write(this->transformed_exp, TransformedVectorSize*NVALS*N);
    }

    template <class BufferReaderClass>
    void deserialize(BufferReaderClass& buffer){
      buffer.fillArray(this->exp, VectorSize*NVALS*N);
      buffer.fillArray(this->transformed_exp, TransformedVectorSize*NVALS*N);
    }

    ///////////////////////////////////////////////////////
    // to extend Serializable
    ///////////////////////////////////////////////////////
    template <class BufferWriterClass>
    void save(BufferWriterClass& buffer) const{
      buffer.write(this->exp, VectorSize*NVALS*NRHS);
      buffer.write(this->transformed_exp, TransformedVectorSize*NVALS*NRHS);
    }

    template <class BufferReaderClass>
    void restore(BufferReaderClass& buffer){
      buffer.fillArray(this->exp, VectorSize*NVALS*NRHS);
      buffer.fillArray(this->transformed_exp, TransformedVectorSize*NVALS*NRHS);
    }

    FSize getSavedSize() const {
      return N * NVALS * VectorSize * (FSize) sizeof(FReal)
	+ N * NVALS * TransformedVectorSize * (FSize) sizeof(stdComplex<FReal>);
    }
  };

  using multipole_t       = exp_impl<class multipole_tag, NRHS>;
  using local_expansion_t = exp_impl<class local_expansion_tag, NLHS>;

  multipole_t m_data {};
  local_expansion_t l_data {};

  bool hasMultipoleData() const noexcept {
    return true;
  }
  bool hasLocalExpansionData() const noexcept {
    return true;
  }


  multipole_t& getMultipoleData() noexcept {
    return m_data;
  }
  const multipole_t& getMultipoleData() const noexcept {
    return m_data;
  }

  local_expansion_t& getLocalExpansionData() noexcept {
    return l_data;
  }
  const local_expansion_t& getLocalExpansionData() const noexcept {
    return l_data;
  }


  ///////////////////////////////////////////////////////
  // to extend FAbstractSendable
  ///////////////////////////////////////////////////////
  template <class BufferWriterClass>
  void serializeUp(BufferWriterClass& buffer) const{
    m_data.serialize(buffer);
  }

  template <class BufferReaderClass>
  void deserializeUp(BufferReaderClass& buffer){
    m_data.deserialize(buffer);
  }

  template <class BufferWriterClass>
  void serializeDown(BufferWriterClass& buffer) const{
    l_data.serialize(buffer);
  }

  template <class BufferReaderClass>
  void deserializeDown(BufferReaderClass& buffer){
    l_data.deserialize(buffer);
  }

  ///////////////////////////////////////////////////////
  // to extend Serializable
  ///////////////////////////////////////////////////////
  template <class BufferWriterClass>
  void save(BufferWriterClass& buffer) const{
    FBasicCell::save(buffer);
    m_data.save(buffer);
    l_data.save(buffer);
  }

  template <class BufferReaderClass>
  void restore(BufferReaderClass& buffer){
    FBasicCell::restore(buffer);
    m_data.restore(buffer);
    l_data.restore(buffer);
  }

  FSize getSavedSize() const {
    return m_data.getSavedSize() + l_data.getSavedSize()
      + FBasicCell::getSavedSize();
  }

  FSize getSavedSizeUp() const {
    return m_data.getSavedSize();
  }

  FSize getSavedSizeDown() const {
    return l_data.getSavedSize();
  }
  ///
  /// \brief reset To zero all multipole and local arrays
  ///
  void resetToInitialState(){
    l_data.reset() ;
    m_data.reset();
  }
  template <class StreamClass>
  friend StreamClass& operator<<(StreamClass& output, const FUnifCell<FReal,ORDER, NRHS, NLHS, NVALS>&  cell){
    output << "Multipole exp NRHS " << NRHS
	   << " NVALS "  << NVALS
	   << " VectorSize "  << cell.getVectorSize()
	   << '\n';
    for (int rhs= 0 ; rhs < NRHS ; ++rhs) {
      const FReal* pole = cell.getMultipole(rhs);
      for (int val= 0 ; val < NVALS ; ++val) {
	output<< "      val : " << val << " exp: " ;
	for (int i= 0 ; i < cell.getVectorSize()  ; ++i) {
	  output<< pole[i] << " ";
	}
	output << std::endl;
      }
    }
    return output;
  }

};

template <class FReal, int ORDER, int NRHS = 1, int NLHS = 1, int NVALS = 1>
class FTypedUnifCell : public FUnifCell<FReal,ORDER,NRHS,NLHS,NVALS>, public FExtendCellType {
public:
  template <class BufferWriterClass>
  void save(BufferWriterClass& buffer) const{
    FUnifCell<FReal,ORDER,NRHS,NLHS,NVALS>::save(buffer);
    FExtendCellType::save(buffer);
  }
  template <class BufferReaderClass>
  void restore(BufferReaderClass& buffer){
    FUnifCell<FReal,ORDER,NRHS,NLHS,NVALS>::restore(buffer);
    FExtendCellType::restore(buffer);
  }
  void resetToInitialState(){
    FUnifCell<FReal,ORDER,NRHS,NLHS,NVALS>::resetToInitialState();
    FExtendCellType::resetToInitialState();
  }
  FSize getSavedSize() const {
    return FExtendCellType::getSavedSize() + FUnifCell<FReal, ORDER,NRHS,NLHS,NVALS>::getSavedSize();
  }
};

#endif //FUNIFCELL_HPP
