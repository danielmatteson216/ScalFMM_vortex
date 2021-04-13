#ifndef FTYPEPACK_HPP
#define FTYPEPACK_HPP

#include <tuple>

namespace scalfmm {

    template<std::size_t I, typename T>
    struct pack {
        enum {size = I};
        using type = T;
    };

    namespace details {

        template<
            template<class...> class Final,
            typename ExpandedTuple,
            typename... ToExpand
            >
        struct pack_expand_impl;

        template<
            template<class...> class Final,
            typename... Expanded
            >
        struct pack_expand_impl<Final, std::tuple<Expanded...> > {
            using type = Final<Expanded...>;
        };

        template<
            template<class...> class Final,
            typename T,
            typename... Args,
            typename... Expanded
            >
        struct pack_expand_impl<Final, std::tuple<Expanded...>, pack<0, T>, Args...> {
            using type = typename pack_expand_impl<Final, std::tuple<Expanded...>, Args...>::type;
        };


        template<
            template<class...> class Final,
            std::size_t count,
            typename T,
            typename... Args,
            typename... Expanded
            >
        struct pack_expand_impl<Final, std::tuple<Expanded...>, pack<count, T>, Args...> {
            using type = typename pack_expand_impl<Final, std::tuple<Expanded..., T>, pack<count-1, T>, Args...>::type;
        };

        template<
            template<class...> class Final,
            typename T,
            typename... Args,
            typename... Expanded
            >
        struct pack_expand_impl<Final, std::tuple<Expanded...>, T, Args...> {
            using type = typename pack_expand_impl<Final, std::tuple<Expanded..., T>, Args...>::type;
        };

    }

    template<template <class...> class VariadicType, typename... TypePacks>
    using pack_expand = typename details::pack_expand_impl<VariadicType, std::tuple<>, TypePacks...>::type;

    template<typename... TypePacks>
    using pack_expand_tuple = pack_expand<std::tuple, TypePacks...>;

}

#endif /* FTYPEPACK_HPP */
