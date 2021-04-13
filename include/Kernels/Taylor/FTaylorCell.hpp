// See LICENCE file at project root
#ifndef FTAYLORCELL_HPP
#define FTAYLORCELL_HPP

#include "Components/FBasicCell.hpp"
#include "Containers/FVector.hpp"
#include "Utils/FMemUtils.hpp"
#include "Extensions/FExtendCellType.hpp"

/**
 *@author Cyrille Piacibello
 *@class FTaylorCell
 *
 *This class is a cell used for the Taylor Expansion Kernel.
 *
 *
 */
template < class FReal, int P, int order>
class FTaylorCell : public FBasicCell, public FAbstractSendable {
protected:

    template<std::size_t S, class Tag>
    struct expansion_impl {
        enum {Size = S};
        FReal exp[Size];

        FReal* get() noexcept {
            return exp;
        }
        const FReal* get() const noexcept {
            return exp;
        }

        void reset() noexcept {
            for(int idx = 0; idx < Size; ++idx) {
                exp[idx].setRealImag(FReal(0.0), FReal(0.0));
            }
        }

        int getSize() const noexcept {
            return Size;
        }

        template<class BufferWriterClass>
        void serialize(BufferWriterClass& buffer) const {
            buffer.write(exp, Size);
        }
        template<class BufferReaderClass>
        void deserialize(BufferReaderClass& buffer) {
            buffer.fillArray(exp, Size);
        }
    };

public:

    using multipole_t = expansion_impl<((P+1)*(P+2)*(P+3))*order/6, class multipole_tag>;
    using local_expansion_t = expansion_impl<((P+1)*(P+2)*(P+3))*order/6, class local_expansion_tag>;

protected:

    multipole_t m_data;
    local_expansion_t l_data;

public:

    const multipole_t& getMultipoleData() const noexcept {
        return m_data;
    }
    multipole_t& getMultipoleData() {
        return m_data;
    }
    const local_expansion_t& getLocalExpansionData() const noexcept {
        return l_data;
    }
    local_expansion_t& getLocalExpansionData() {
        return l_data;
    }

    /** Make it like the begining */
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

    FSize getSavedSizeUp() const {
        return ((FSize) sizeof(FReal)) * (multipole_t::Size);
    }

    FSize getSavedSizeDown() const {
        return ((FSize) sizeof(FReal)) * (local_expansion_t::Size);
    }

    ///////////////////////////////////////////////////////
    // to extend Serializable
    ///////////////////////////////////////////////////////
    template <class BufferWriterClass>
    void save(BufferWriterClass& buffer) const{
        FBasicCell::save(buffer);
        m_data.serialize();
        l_data.serialize();
    }
    template <class BufferReaderClass>
    void restore(BufferReaderClass& buffer){
        FBasicCell::restore(buffer);
        m_data.deserialize();
        l_data.deserialize();
    }

    FSize getSavedSize() const {
        return FSize((sizeof(FReal) * (multipole_t::Size + local_expansion_t::Size)
                     + FBasicCell::getSavedSize()));
    }
};

template <class FReal, int P, int order>
class FTypedTaylorCell : public FTaylorCell<FReal, P,order>, public FExtendCellType {
public:
    template <class BufferWriterClass>
    void save(BufferWriterClass& buffer) const{
        FTaylorCell<FReal,P,order>::save(buffer);
        FExtendCellType::save(buffer);
    }
    template <class BufferReaderClass>
    void restore(BufferReaderClass& buffer){
        FTaylorCell<FReal,P,order>::restore(buffer);
        FExtendCellType::restore(buffer);
    }
    void resetToInitialState(){
        FTaylorCell<FReal,P,order>::resetToInitialState();
        FExtendCellType::resetToInitialState();
    }

    FSize getSavedSize() const {
        return FExtendCellType::getSavedSize() + FTaylorCell<FReal, P,order>::getSavedSize();
    }
};

#endif
