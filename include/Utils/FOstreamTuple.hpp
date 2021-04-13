#ifndef OSTREAM_TUPLE_HPP
#define OSTREAM_TUPLE_HPP

#include "inria/ostream_joiner.hpp"
#include "inria/integer_sequence.hpp"

#include <ostream>
#include <tuple>

namespace scalfmm {
namespace details {
namespace tuple_helper {
/** \brief Helper for tuple formatted output
 *
 * \param os Output stream
 * \param t  Printed tuple
 * \param index_sequence unnamed parameter for automatic deduction of
 * Indices template parameter
 *
 * \tparam Types Types contained in the tuple, automatically deduced
 * \tparam Indices Type indices, automatically deduced
 */
template<typename...Types, std::size_t... Indices>
inline void formated_output(std::ostream& os, const std::tuple<Types...>& t,
                            inria::index_sequence<Indices...>)
{
    os << "(" ;     // opening brace
    // initializer_list members are evaluated from left to right; this
    // evaluates to an equivalent to:
    //
    //     os << std::get<0>(t) << ", ";
    //     os << std::get<1>(t) << ", "
    //     ...
    //     std::get<N>(t) << "";
    auto l = { (os << std::get<Indices>(t)
                << (Indices != sizeof...(Types)-1 ? ", " : ""), 0) ... };
    (void)l;        // ignore fact that initializer_list is not used
    os << ")" ;     // closing brace
}

template<class... Ts>
struct tuple_wrapper {
    const std::tuple<Ts...>& t;
    friend std::ostream& operator<<(std::ostream& os, const tuple_wrapper& w) {
        using details::tuple_helper::formated_output;
        formated_output(os, w.t, inria::index_sequence_for<Ts...>{});
        return os;
    }
};

} } // close namespace [scalfmm]::details::tuple_helper


namespace {

constexpr struct {
    template<class... Ts>
    details::tuple_helper::tuple_wrapper<Ts...> operator()(const std::tuple<Ts...>& t) const {
        return {t};
    }
} tuple_out {};


} // close namsespace scalfmm::details

namespace details {
/**
 * \brief Silence the `unused variable` warnings about tuple_out
 */
template<class T>
void silence_tuple_out_warning() {
    tuple_out(std::tuple<>{});
}

} // close namespace scalfmm::details

} // close namespace scalfmm

#endif
