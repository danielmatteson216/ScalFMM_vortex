// See LICENCE file at project root
#ifndef FROTATIONCELL_HPP
#define FROTATIONCELL_HPP

#include "Utils/FComplex.hpp"
#include "Utils/FMemUtils.hpp"

#include "Extensions/FExtendCellType.hpp"

#include "Components/FBasicCell.hpp"


/** This class is a cell used for the rotation based kernel
  * The size of the multipole and local vector are based on a template
  * User should choose this parameter P carrefuly to match with the
  * P of the kernel.
  *
  * Multipole/Local vectors contain value as:
  * {0,0}{1,0}{1,1}...{P,P-1}{P,P}
  * So the size of such vector can be obtained by a suite:
  * (n+1)*n/2 => (P+2)*(P+1)/2
  */
template <class FReal, int P>
class FRotationCell : public FBasicCell, public FAbstractSendable {
protected:

    template<std::size_t S, class Tag>
    struct expansion_impl {
        enum {Size = S};
        FComplex<FReal> exp[Size];

        FComplex<FReal>* get() noexcept {
            return exp;
        }
        const FComplex<FReal>* get() const noexcept {
            return exp;
        }

        FSize getSize() noexcept {
            return Size;
        }

        FSize getSavedSize() const noexcept {
            return ((FSize) sizeof(exp[0])) * Size;
        }

        void reset() noexcept {
            for(int idx = 0; idx < Size; ++idx) {
                exp[idx].setRealImag(FReal(0.0), FReal(0.0));
            }
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

    using multipole_t = expansion_impl<((P+2)*(P+1))/2, class multipole_tag>;
    using local_expansion_t = expansion_impl<((P+2)*(P+1))/2, class local_expansion_tag>;

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
        return ((FSize) sizeof(FComplex<FReal>)) * (multipole_t::Size);
    }

    FSize getSavedSizeDown() const {
        return ((FSize) sizeof(FComplex<FReal>)) * (local_expansion_t::Size);
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
        return FSize(((int) sizeof(FComplex<FReal>)) * (multipole_t::Size + local_expansion_t::Size)
                + FBasicCell::getSavedSize());
    }
};

template <class FReal, int P>
class FTypedRotationCell : public FRotationCell<FReal, P>, public FExtendCellType {
public:
    template <class BufferWriterClass>
    void save(BufferWriterClass& buffer) const{
        FRotationCell<FReal, P>::save(buffer);
        FExtendCellType::save(buffer);
    }
    template <class BufferReaderClass>
    void restore(BufferReaderClass& buffer){
        FRotationCell<FReal, P>::restore(buffer);
        FExtendCellType::restore(buffer);
    }
    void resetToInitialState(){
        FRotationCell<FReal, P>::resetToInitialState();
        FExtendCellType::resetToInitialState();
    }

    FSize getSavedSize() const {
        return FExtendCellType::getSavedSize() + FRotationCell<FReal, P>::getSavedSize();
    }
};

#endif // FROTATIONCELL_HPP
