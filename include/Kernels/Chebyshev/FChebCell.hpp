// See LICENCE file at project root

#ifndef FCHEBCELL_HPP
#define FCHEBCELL_HPP
#include <iostream>

#include "Components/FBasicCell.hpp"

#include "FChebTensor.hpp"
#include "Extensions/FExtendCellType.hpp"

/**
 * @author Matthias Messner (matthias.messner@inria.fr)
 * @class FChebCell
 * Please read the license
 *
 * This class defines a cell used in the Chebyshev based FMM.
 * @tparam NVALS is the number of right hand side.
 */
template <class FReal, int ORDER, int NRHS = 1, int NLHS = 1, int NVALS = 1>
class FChebCell : public FBasicCell, public FAbstractSendable
{
    // nnodes = ORDER^3
    // we multiply by 2 because we store the  Multipole expansion end the compressed one.
    static constexpr int VectorSize = TensorTraits<ORDER>::nnodes * 2;

public:

    template<class Tag, std::size_t N>
    struct exp_impl {
        FReal exp[N * NVALS * VectorSize];

        const FReal* get(const int inRhs) const
        { return this->exp + inRhs*VectorSize; }
        FReal* get(const int inRhs)
        { return this->exp + inRhs*VectorSize; }

        constexpr int getVectorSize() const {
            return VectorSize;
        }

        // to extend FAbstractSendable
        template <class BufferWriterClass>
        void serialize(BufferWriterClass& buffer) const{
            buffer.write(this->exp, VectorSize*NVALS*NRHS);
        }
        template <class BufferReaderClass>
        void deserialize(BufferReaderClass& buffer){
            buffer.fillArray(this->exp, VectorSize*NVALS*NRHS);
        }

        void reset() {
            memset(this->exp, 0, sizeof(FReal) * N * NVALS * VectorSize);
        }

        FSize getSavedSize() const {
            return FSize(sizeof(FReal)) * VectorSize * N * NVALS;
        }


    };

    using multipole_t       = exp_impl<class multipole_tag, NRHS>;
    using local_expansion_t = exp_impl<class local_expansion_tag, NLHS>;

    multipole_t       m_data {};
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



    /** To get the leading dim of a vec */
    int getVectorSize() const{
        return VectorSize;
    }

    ///
    /// Make it like the begining
    ///
    void resetToInitialState(){
        m_data.reset();
        l_data.reset();
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
        m_data.serialize(buffer);
        l_data.serialize(buffer);
    }
    template <class BufferReaderClass>
    void restore(BufferReaderClass& buffer){
        FBasicCell::restore(buffer);
        m_data.deserialize(buffer);
        l_data.deserialize(buffer);
    }

    FSize getSavedSize() const {
        return FSize(sizeof(FReal)) * VectorSize*(NRHS+NLHS)*NVALS + FBasicCell::getSavedSize();
    }

    FSize getSavedSizeUp() const {
        return FSize(sizeof(FReal)) * VectorSize*(NRHS)*NVALS;
    }

    FSize getSavedSizeDown() const {
        return FSize(sizeof(FReal)) * VectorSize*(NLHS)*NVALS;
    }

    //	template <class StreamClass>
    //	const void print(StreamClass& output) const{
    template <class StreamClass>
    friend StreamClass& operator<<(StreamClass& output, const FChebCell<FReal, ORDER, NRHS, NLHS, NVALS>&  cell){
        //	const void print() const{
        output <<"  Multipole exp NRHS " <<NRHS <<" NVALS "  <<NVALS << " VectorSize/2 "  << cell.getVectorSize() *0.5<< std::endl;
        for (int rhs= 0 ; rhs < NRHS ; ++rhs) {
            const FReal* pole = cell.get(rhs);
            for (int val= 0 ; val < NVALS ; ++val) {
                output<< "      val : " << val << " exp: " ;
                for (int i= 0 ; i < cell.getVectorSize()/2  ; ++i) {
                    output<< pole[i] << " ";
                }
                output << std::endl;
            }
        }
        return output;
    }

};

template <class FReal, int ORDER, int NRHS = 1, int NLHS = 1, int NVALS = 1>
class FTypedChebCell : public FChebCell<FReal, ORDER,NRHS,NLHS,NVALS>, public FExtendCellType {
public:
    template <class BufferWriterClass>
    void save(BufferWriterClass& buffer) const{
        FChebCell<FReal,ORDER,NRHS,NLHS,NVALS>::save(buffer);
        FExtendCellType::save(buffer);
    }
    template <class BufferReaderClass>
    void restore(BufferReaderClass& buffer){
        FChebCell<FReal,ORDER,NRHS,NLHS,NVALS>::restore(buffer);
        FExtendCellType::restore(buffer);
    }
    void resetToInitialState(){
        FChebCell<FReal,ORDER,NRHS,NLHS,NVALS>::resetToInitialState();
        FExtendCellType::resetToInitialState();
    }


    FSize getSavedSize() const {
        return FExtendCellType::getSavedSize() + FChebCell<FReal, ORDER,NRHS,NLHS,NVALS>::getSavedSize();
    }

};
#endif //FCHEBCELL_HPP
