#ifndef FKERNEL_CONCEPTS_HPP
#define FKERNEL_CONCEPTS_HPP

#include <type_traits>

#include "inria/logic.hpp"

#include "Containers/FTreeCoordinate.hpp"

namespace scalfmm {
    namespace meta {
        namespace details {

            #ifdef declare_has_method
            #error declare_has_method macro is already defined... This needs fixing.
            #endif

            #define declare_has_method(name)                            \
                template<typename K, typename Ret, typename ... Args>   \
                struct has_##name {                                     \
                    template<typename...> using void_t = void;          \
                    template<typename k, typename = void_t<> >          \
                    struct check : std::false_type {};                  \
                    template<typename k>                                \
                    struct check<k, void_t<decltype(std::declval<k>(). name(std::declval<Args>()...))> > \
                        : std::true_type {};                            \
                    constexpr static const bool value = check<K>::value; \
                };

            declare_has_method(P2M);
            declare_has_method(P2L);
            declare_has_method(M2M);
            declare_has_method(M2P);
            declare_has_method(M2L);
            declare_has_method(L2L);
            declare_has_method(L2P);
            declare_has_method(P2P);

            #undef declare_has_method
        }

        template<typename Tree, typename Kernel>
        struct has_M2P {
            enum : bool {
                value =
                    details::has_M2P<Kernel, void,
                                     typename Tree::node_t::data_t::multipole_t*,
                                     typename Tree::node_t::symbolic_data_t*,
                                     typename Tree::node_t::particle_container_t*,
                                     typename Tree::node_t::symbolic_data_t*
                                     >::value
                    };
        };

        template<typename Tree, typename Kernel>
        struct has_P2L {
            enum : bool {
                value =
                    details::has_P2L<Kernel, void,
                                     typename Tree::node_t::data_t::local_expansion_t*,
                                     typename Tree::node_t::symbolic_data_t*,
                                     typename Tree::node_t::particle_container_t*,
                                     typename Tree::node_t::symbolic_data_t*
                                     >::value
                    };
        };

        template<typename Tree, typename Kernel>
        struct has_P2M {
            enum : bool {
                value =
                    details::has_P2M<Kernel, void,
                                     typename Tree::node_t::data_t::multipole_t*,
                                     typename Tree::node_t::symbolic_data_t*,
                                     typename Tree::node_t::particle_container_t*
                                     >::value
                    };
        };

        template<typename Tree, typename Kernel>
        struct has_M2M {
            enum : bool {
                value =
                    details::has_M2M<Kernel, void,
                                     typename Tree::node_t::data_t::multipole_t*,
                                     typename Tree::node_t::symbolic_data_t*,
                                     typename Tree::node_t::data_t::multipole_t**,
                                     typename Tree::node_t::symbolic_data_t**
                                     >::value
                    };
        };

                template<typename Tree, typename Kernel>
        struct has_M2L {
            enum : bool {
                value =
                    details::has_M2L<Kernel, void,
                                     typename Tree::node_t::data_t::local_expansion_t*,
                                     typename Tree::node_t::symbolic_data_t*,
                                     typename Tree::node_t::data_t::multipole_t**,
                                     typename Tree::node_t::symbolic_data_t**,
                                     int*, int
                                     >::value
                    };
        };

        template<typename Tree, typename Kernel>
        struct has_L2L {
            enum : bool {
                value =
                    details::has_L2L<Kernel, void,
                                     typename Tree::node_t::data_t::local_expansion_t*,
                                     typename Tree::node_t::symbolic_data_t*,
                                     typename Tree::node_t::data_t::local_expansion_t**,
                                     typename Tree::node_t::symbolic_data_t**
                                     >::value
                    };
        };

        template<typename Tree, typename Kernel>
        struct has_L2P {
            enum : bool {
                value =
                    details::has_L2P<Kernel, void,
                                     typename Tree::node_t::data_t::local_expansion_t*,
                                     typename Tree::node_t::symbolic_data_t*,
                                     typename Tree::node_t::particle_container_t*
                                     >::value
                    };
        };

        template<typename Tree, typename Kernel>
        struct has_P2P {
            enum : bool {
                value =
                    details::has_P2P<Kernel, void,
                                     FTreeCoordinate,
                                     typename Tree::node_t::particle_container_t*,
                                     typename Tree::node_t::particle_container_t*,
                                     typename Tree::node_t::particle_container_t**,
                                     int*, int
                                     >::value
                    };
        };

        template<typename Tree, typename Kernel>
        struct has_partial_P2P {
            enum : bool {
                value =
                    details::has_P2P<Kernel, void,
                                     FTreeCoordinate,
                                     typename Tree::node_t::particle_container_t*,
                                     typename Tree::node_t::particle_container_t*,
                                     typename Tree::node_t::particle_container_t**,
                                     int*, int, bool
                                     >::value
                    };
        };

        template<class Tree, class Kernel>
        using adaptive_compatible = inria::conjunction<
            has_P2M<Tree, Kernel>,
            has_M2M<Tree, Kernel>,
            has_M2L<Tree, Kernel>,
            has_L2L<Tree, Kernel>,
            has_L2P<Tree, Kernel>,
            has_P2P<Tree, Kernel>,
            has_M2P<Tree, Kernel>,
            has_P2L<Tree, Kernel>
            >;

        template<class Tree, class Kernel>
        constexpr bool check_adaptive_compatible() noexcept {
            #ifdef CHECK_CONCEPT_FAILURE
            #error "Macro CHECK_CONCEPT_FAILURE is already defined..."
            #endif
            #define CHECK_CONCEPT_FAILURE(NAME)  \
                static_assert(has_##NAME<Tree,Kernel>::value, #NAME);
            CHECK_CONCEPT_FAILURE(P2M);
            CHECK_CONCEPT_FAILURE(M2M);
            CHECK_CONCEPT_FAILURE(M2L);
            CHECK_CONCEPT_FAILURE(L2L);
            CHECK_CONCEPT_FAILURE(L2P);
            CHECK_CONCEPT_FAILURE(P2P);
            CHECK_CONCEPT_FAILURE(M2P);
            CHECK_CONCEPT_FAILURE(P2L);
            #undef CHECK_CONCEPT_FAILURE
            return adaptive_compatible<Tree,Kernel>::value;
        }
    }
}


#endif /* FKERNEL_CONCEPTS_HPP */
