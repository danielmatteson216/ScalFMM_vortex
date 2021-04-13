#ifndef _FSYMBOLICDATA_HPP_
#define _FSYMBOLICDATA_HPP_

#include <cstdlib>

#include "Utils/FGlobal.hpp"

#include "inria/logic.hpp"
#include "Containers/FTreeCoordinate.hpp"

/** \brief Node structural data
 * Keeps data about the node that may be read by kernels or algorithms.
 */
struct FSymbolicData {
    /// Node depth in its tree
    int depth;
    /// Node index in parent child array
    std::size_t m_idx;

    void setLevel(int value) noexcept {
        this->depth = value;
    }
    int getLevel() const noexcept {
        return this->depth;
    }
    void setMortonIndex(MortonIndex value) noexcept {
        this->m_idx = value;
    }
    MortonIndex getMortonIndex() const noexcept {
        return this->m_idx;
    }
    void setCoordinate(const FTreeCoordinate& coord) noexcept {
        this->m_idx = coord.getMortonIndex();
    }
    FTreeCoordinate getCoordinate() const noexcept {
        return FTreeCoordinate(this->m_idx);
    }
    friend std::ostream& operator<<(std::ostream& os, const FSymbolicData& d) {
        return (os << '{' << "\"depth\":" <<d.depth << ',' << "\"index\":" << d.m_idx << '}');
    }
};

namespace inria {
    struct any_t {};
}

#ifdef MAKE_HAS_TRAIT
#error "MAKE_HAS_METHOD_TRAIT macro is already set..."
#else
#define MAKE_HAS_METHOD_TRAIT(method_name_to_be_checked)                \
                                                                        \
    template<class T, class Signature>                                  \
    struct has_method_##method_name_to_be_checked;                      \
                                                                        \
    template<class T, class Ret, class... Args>                         \
    struct has_method_##method_name_to_be_checked<T, Ret(Args...)> {    \
        template<class U,                                               \
                 class R = decltype(std::declval<U>().method_name_to_be_checked(std::declval<Args>()...)), \
                 class = typename std::enable_if<                       \
                     std::conditional<std::is_same<inria::any_t, Ret>::value, \
                                      std::true_type,                   \
                                      std::is_convertible<R, Ret>>::type::value>::type \
                 >                                                      \
        static constexpr bool check(U*){return true;}                   \
        static constexpr bool check(...){return false;}                 \
        enum {value = check(static_cast<T*>(nullptr))};                 \
    };

#endif

namespace scalfmm {
    namespace meta {

        MAKE_HAS_METHOD_TRAIT(getLevel);
        MAKE_HAS_METHOD_TRAIT(getMortonIndex);
        MAKE_HAS_METHOD_TRAIT(getCoordinate);

    }
}

#undef MAKE_HAS_METHOD_TRAIT


template<class T>
struct has_type_symbolic_data_t {
    template<class U, class = typename U::symbolic_data_t>
    static constexpr bool check(U*) {return true;}
    static constexpr bool check(...) {return false;}
    enum {value = check(static_cast<T*>(nullptr))};
};

template<class Default, class T, bool = has_type_symbolic_data_t<T>::value>
struct extract_symbolic_data_t_or_fallback_to_default {
    using type = Default;
};
template<class Default, class T>
struct extract_symbolic_data_t_or_fallback_to_default<Default, T, true> {
    using type = typename T::symbolic_data_t;
};


template<class T>
struct models_symbolic_data {
    template<class U,
             class = inria::require<
                 scalfmm::meta::has_method_getLevel<const T, inria::any_t()>,
                 scalfmm::meta::has_method_getCoordinate<const T, FTreeCoordinate()>,
                 scalfmm::meta::has_method_getMortonIndex<const T, inria::any_t()>
                 >
             >
    static constexpr bool check(U*) {return true;}
    static constexpr bool check(...) {return false;}
    enum {value = check(static_cast<T*>(nullptr))};
};




#endif /* _FSYMBOLICDATA_HPP_ */
