#ifndef _SCALFMM_KERNEL_UTILITIES_HPP_
#define _SCALFMM_KERNEL_UTILITIES_HPP_

#include "inria/meta.hpp"

namespace scalfmm {

namespace details {

/**
 * \brief Check for kernel.setup(tree) existence
 */
template<class Kernel, class... Args>
struct has_setup {
    template<class K>
    constexpr static auto check(K* k, Args* ... args)
        -> decltype(k->setup(*args...), void(), true)
    { return true; }
    constexpr static bool check(...) {return false;}
    static constexpr bool value = check((Kernel*)0, (Args*)0 ...);
};

/**
 * \brief Call kernel.setup(tree)
 *
 * \param kernel Kernel to setup
 * \param tree   Tree to setup with
 *
 * \tparam Kernel Kernel type
 * \tparam Tree   Tree type
 */
template<class Kernel, class Tree>
void setup_kernel(has_setup<Kernel, Tree>, Kernel& kernel, Tree& tree) {
    kernel.setup(tree);
}

/**
 * \brief Call kernel.setup()
 *
 * \param kernel Kernel to setup
 *
 * \tparam Kernel Kernel type
 * \tparam Tree   Tree type, the tree is not used
 */
template<class Kernel, class Tree>
void setup_kernel(has_setup<Kernel>, Kernel& kernel, const Tree&) {
    kernel.setup();
}

/**
 * \brief No-op
 *
 * \tparam Args Unused arguments
 */
template<class... Args>
void setup_kernel(std::true_type, const Args&...) {}


} // close namespace [scalfmm]::details

/**
 * \brief Call kernel setup method with arguments if possible
 *
 * Current behaviour calls (first defined):
 *   - kernel.setup(tree)
 *   - kernel.setup()
 *
 * \param kernel Kernel to setup
 * \param tree   Tree to setup kernel (if needed)
 */
template<class Kernel, class Tree>
void setup_kernel(Kernel& kernel, Tree& tree) {
    using tag = inria::first_true_t<
        details::has_setup<Kernel, Tree>,
        details::has_setup<Kernel>,
        std::true_type
        >;
    details::setup_kernel(tag{}, kernel, tree);
}




namespace details {

/**
 * \brief Check for kernel.cleanup(tree) existence
 */
template<class Kernel, class... Args>
struct has_cleanup {
    template<class K>
    constexpr static auto check(K* k, Args* ... args)
        -> decltype(k->cleanup(*args...), void(), true)
    { return true; }
    constexpr static bool check(...) {return false;}
    static constexpr bool value = check((Kernel*)0, (Args*)0 ...);
};

/**
 * \brief Call kernel.cleanup(tree)
 *
 * \param kernel Kernel to cleanup
 * \param tree   Tree to cleanup with
 *
 * \tparam Kernel Kernel type
 * \tparam Tree   Tree type
 */
template<class Kernel, class Tree>
void cleanup_kernel(has_cleanup<Kernel, Tree>, Kernel& kernel, Tree& tree) {
    kernel.cleanup(tree);
}

/**
 * \brief Call kernel.cleanup()
 *
 * \param kernel Kernel to cleanup
 *
 * \tparam Kernel Kernel type
 * \tparam Tree   Tree type, the tree is not used
 */
template<class Kernel, class Tree>
void cleanup_kernel(has_cleanup<Kernel>, Kernel& kernel, const Tree&) {
    kernel.cleanup();
}

/**
 * \brief No-op
 *
 * \tparam Args Unused arguments
 */
template<class... Args>
void cleanup_kernel(std::true_type, const Args&...) {}


} // close namespace [scalfmm]::details

/**
 * \brief Call kernel cleanup method with arguments if possible
 *
 * Current behaviour calls (first defined):
 *   - kernel.cleanup(tree)
 *   - kernel.cleanup()
 *
 * \param kernel Kernel to cleanup
 * \param tree   Tree to cleanup kernel (if needed)
 */
template<class Kernel, class Tree>
void cleanup_kernel(Kernel& kernel, Tree& tree) {
    using tag = inria::first_true_t<
        details::has_cleanup<Kernel, Tree>,
        details::has_cleanup<Kernel>,
        std::true_type
        >;
    details::cleanup_kernel(tag{}, kernel, tree);
}












} // close namespace scalfmm

#endif /* _SCALFMM_KERNEL_UTILITIES_HPP_ */
