// See LICENCE file at project root
#ifndef FSPHERICALCELL_HPP
#define FSPHERICALCELL_HPP

#include <algorithm>

#include "Utils/FComplex.hpp"
#include "Utils/FMemUtils.hpp"

#include "Extensions/FExtendCellType.hpp"

#include "Components/FBasicCell.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
*/
template <class FReal>
class FSphericalCell : public FBasicCell, public FAbstractSendable {
protected:
    static int DevP;
    static bool UseBlas;

    template<class PoleTag>
    struct exp_impl {
        using cell_t = FSphericalCell;
        static int Size;

        static int getSize() noexcept {
            return Size;
        }

        FComplex<FReal>* exp;

        exp_impl() : exp{nullptr} {
            this->exp = new FComplex<FReal>[Size];
        }

        exp_impl(const exp_impl& other) : exp_impl() {
            *this = other;
        }

        exp_impl(exp_impl&& other) : exp(nullptr) {
            *this = std::move(other);
        }

        exp_impl& operator=(const exp_impl& other) {
            FMemUtils::copyall(exp, other.exp, Size);
            return *this;
        }

        exp_impl& operator=(exp_impl&& other) {
            using std::swap;
            swap(this->exp, other.exp);
            return *this;
        }

        ~exp_impl() {
            delete[] this->exp;
        }

        FComplex<FReal>* get() noexcept {
            return exp;
        }
        const FComplex<FReal>* get() const noexcept {
            return exp;
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

        FSize getSavedSize() const {
            return ((FSize) sizeof(FComplex<FReal>)) * Size;
        }

        friend bool operator==(const exp_impl& lhs, const exp_impl& rhs) {
            return std::equal(lhs, lhs + Size, rhs.exp);
        }

        friend bool operator!=(const exp_impl& lhs, const exp_impl& rhs) {
            return ! (lhs == rhs);
        }
    };


public:

    using multipole_t = exp_impl<struct multipole_tag>;
    using local_expansion_t = exp_impl<struct local_expansion_tag>;



protected:

    multipole_t m_data {};
    local_expansion_t l_data {};

public:
    static void Init(const int inDevP, const bool inUseBlas = false){
        DevP  = inDevP;
        const int ExpP  = int((inDevP+1) * (inDevP+2) * 0.5);
        const int NExpP = (inDevP+1) * (inDevP+1);

        local_expansion_t::Size = ExpP;
        if(inUseBlas) {
            multipole_t::Size = NExpP;
        }
        else{
            multipole_t::Size = ExpP;
        }
    }


    /** Default constructor */
    FSphericalCell() : m_data(), l_data() {}
    /** Copy constructor */
    FSphericalCell(const FSphericalCell& other) = default;
    /** Move constructor */
    FSphericalCell(FSphericalCell&& other) = default;
    /** Copy operator */
    FSphericalCell& operator=(const FSphericalCell&) = default;
    /** Move operator */
    FSphericalCell& operator=(FSphericalCell&&) = default;
    /** Default destructor */
    virtual ~FSphericalCell(){}

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


    /** Make it like the begining */
    void resetToInitialState(){
        m_data.reset();
        l_data.reset();
    }

    ///////////////////////////////////////////////////////
    // to extend FAbstractSendable
    ///////////////////////////////////////////////////////
    template <class BufferWriterClass>
    void serializeUp(BufferWriterClass& buffer) const {
        m_data.serialize(buffer);
    }
    template <class BufferReaderClass>
    void deserializeUp(BufferReaderClass& buffer){
        m_data.deserialize(buffer);
    }

    template <class BufferWriterClass>
    void serializeDown(BufferWriterClass& buffer) const {
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
    void save(BufferWriterClass& buffer) const {
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

    FSize getSavedSize() const noexcept {
        return m_data.getSavedSize() + l_data.getSavedSize()
            + FBasicCell::getSavedSize();
    }

    FSize getSavedSizeUp() const noexcept {
        return m_data.getSavedSize();
    }

    FSize getSavedSizeDown() const noexcept {
        return l_data.getSavedSize();
    }
};

template <class FReal>
int FSphericalCell<FReal>::DevP(-1);

template <class FReal> template<class Tag>
int FSphericalCell<FReal>::exp_impl<Tag>::Size(-1);


/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
*/
template <class FReal>
class FTypedSphericalCell : public FSphericalCell<FReal>, public FExtendCellType {
public:
    template <class BufferWriterClass>
    void save(BufferWriterClass& buffer) const{
        FSphericalCell<FReal>::save(buffer);
        FExtendCellType::save(buffer);
    }
    template <class BufferReaderClass>
    void restore(BufferReaderClass& buffer){
        FSphericalCell<FReal>::restore(buffer);
        FExtendCellType::restore(buffer);
    }
    void resetToInitialState(){
        FSphericalCell<FReal>::resetToInitialState();
        FExtendCellType::resetToInitialState();
    }

    FSize getSavedSize() const noexcept {
        return FExtendCellType::getSavedSize() + FSphericalCell<FReal>::getSavedSize();
    }
};



#endif //FSPHERICALCELL_HPP
