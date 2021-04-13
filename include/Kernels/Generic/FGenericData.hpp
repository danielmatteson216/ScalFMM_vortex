#ifndef _FGENERICCELL_HPP_
#define _FGENERICCELL_HPP_

#include "Components/FBasicCell.hpp"


template<class Multipole, class LocalExpansion>
class FGenericData : public FBasicCell {
public:
    using multipole_t = Multipole;
    using local_expansion_t = LocalExpansion;
private:

    multipole_t m_data {};
    local_expansion_t l_data {};
public:
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

    void setMultipoleData(const multipole_t& value) noexcept(noexcept(m_data = value)) {
        m_data = value;
    }
    void setMultipoleData(multipole_t&& value) noexcept(noexcept(m_data = std::move(value))) {
        m_data = std::move(value);
    }

    local_expansion_t& getLocalExpansionData() noexcept {
        return l_data;
    }
    const local_expansion_t& getLocalExpansionData() const noexcept {
        return l_data;
    }

    void setLocalExpansionData(const local_expansion_t& value) noexcept(noexcept(l_data = value)) {
        l_data = value;
    }
    void setLocalExpansionData(local_expansion_t&& value) noexcept(noexcept(l_data = std::move(value))) {
        l_data = std::move(value);
    }

    multipole_t& getDataUp() noexcept {
        return m_data;
    }
    const multipole_t& getDataUp() const noexcept {
        return m_data;
    }

    void setDataUp(const multipole_t& value) noexcept(noexcept(m_data = value)) {
        m_data = value;
    }
    void setDataUp(multipole_t&& value) noexcept(noexcept(m_data = std::move(value))) {
        m_data = std::move(value);
    }

    local_expansion_t& getDataDown() noexcept {
        return l_data;
    }
    const local_expansion_t& getDataDown() const noexcept {
        return l_data;
    }

    void setDataDown(const local_expansion_t& value) noexcept(noexcept(l_data = value)) {
        l_data = value;
    }
    void setDataDown(local_expansion_t&& value) noexcept(noexcept(l_data = std::move(value))) {
        l_data = std::move(value);
    }

    /* TODO: Add getSavedSize if multipole and local_exp have it */

    template<class T>
    struct has_getSavedSize {
        template<class U, class V = decltype(std::declval<U>().getSavedSize())>
        static constexpr bool check(U*) {return true;}
        static constexpr bool check(...) {return false;}
        enum {value = check(static_cast<T*>(nullptr))};
    };


    template<class T,
             typename std::enable_if<has_getSavedSize<T>::value, int>::type = 0
             >
    auto call_getSavedSize(const T& obj) const noexcept -> decltype(obj.getSavedSize()) {
        return obj.getSavedSize();
    }

    template<class T,
             typename std::enable_if<! has_getSavedSize<T>::value, int>::type = 0,
             typename std::enable_if<std::is_trivially_copyable<T>::value, int>::type = 0
             >
    FSize call_getSavedSize(const T&) const noexcept {
        return sizeof(T);
    }



    FSize getSavedSize() const noexcept {
        return call_getSavedSize(m_data) + call_getSavedSize(l_data)
            + FBasicCell::getSavedSize();
    }

    FSize getSavedSizeUp() const noexcept {
        return call_getSavedSize(m_data);
    }

    FSize getSavedSizeDown() const noexcept {
        return call_getSavedSize(l_data);
    }

    /* TODO: Add save/restore if multipole and local_exp have it */


    ///////////////////////////////////////////////////////
    // to extend FAbstractSendable
    ///////////////////////////////////////////////////////
    template <class BufferWriterClass>
    void serializeUp(BufferWriterClass& buffer) const {
        m_data.serialize(buffer);
    }

    template <class BufferReaderClass>
    void deserializeUp(BufferReaderClass& buffer) {
        m_data.deserialize(buffer);
    }

    template <class BufferWriterClass>
    void serializeDown(BufferWriterClass& buffer) const {
        l_data.serialize(buffer);
    }

    template <class BufferReaderClass>
    void deserializeDown(BufferReaderClass& buffer) {
        l_data.deserialize(buffer);
    }

};


#endif /* _FGENERICCELL_HPP_ */
