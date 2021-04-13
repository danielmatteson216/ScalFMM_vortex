#ifndef SDT_ALIGN_POLYFILL_HPP_
#define SDT_ALIGN_POLYFILL_HPP_


namespace scalfmm {
    namespace details {

        template<class>
        struct sfinae_false : std::false_type {};

        template<class T>
        static auto test_align(int)
            -> sfinae_false<
            decltype( align(std::declval<std::size_t>(),
                            std::declval<std::size_t>(),
                            std::declval<T*&>(),
                            std::declval<std::size_t>()
                          ))
            >;

        template<class>
        static auto test_align(...) -> std::true_type;


        template<class T>
        struct has_not_align : decltype(test_align<T>(0)){};

        template<class T>
        using has_not_align_t = typename std::enable_if<has_not_align<T>::value, void>::type;

    }
}

namespace std {

    /** Shamelessly copied from https://gcc.gnu.org/bugzilla/show_bug.cgi?id=57350#c11 */
    template< typename T = scalfmm::details::has_not_align_t<void> >
    inline T* align( std::size_t alignment, std::size_t size,
                        T *&ptr, std::size_t &space ) {
	std::uintptr_t pn = reinterpret_cast< std::uintptr_t >( ptr );
	std::uintptr_t aligned = ( pn + alignment - 1 ) & - alignment;
	std::size_t padding = aligned - pn;
	if ( space < size + padding ) return nullptr;
	space -= padding;
	return ptr = reinterpret_cast< void * >( aligned );
    }
}

#endif
