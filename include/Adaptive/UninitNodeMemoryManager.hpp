
#ifndef UNINITNODEMEMORYMANAGER_HPP
#define UNINITNODEMEMORYMANAGER_HPP

#include <functional>
#include <list>
#include <vector>

#include "UninitialisedMemoryProvider.hpp"

/**
 * \brief Manages uninitialised memory for nodes
 *
 * This class' intent is to provide block memory management for nodes. This
 * reduces the amount of memory allocations.
 *
 * The nodes' memory is handled per level. Memory for a level is not garanteed
 * to be contiguous, it is however contiguous by block.
 *
 * The memory requests are forwarded to UninitialisedMemoryProvider objects. A
 * vector of list of providers is kept to ensure constant access to a level. The
 * lists allow constant access to their last provider and allow adding new ones
 * without reallocating the entire list.
 *
 * Memory is alway requested from the last provider in a list. When the last
 * provider does not have enough memory available, a new one is appended to the
 * list.
 *
 * \tparam _Node Node type
 */
template<class _Node>
class UninitNodeMemoryManager {
    /**
     * \brief Policy for memory provider capacity
     *
     * \param level Level of the nodes in the tree
     * \param count Node count
     *
     * \return The provider size
     */
    std::function<std::size_t(std::size_t level, std::size_t count)>
    provider_size_policy = [](const std::size_t level, const std::size_t count) {
        enum {max_level = 4};
        return std::max(std::size_t(1) << 3*(level > max_level ? max_level : level), count);
    };

    /**
     * \brief memory provider lists
     *
     * The vector allows constant access to any level, the lists allow constant
     * access to the last element.
     */
    std::vector<std::list<UninitialisedMemoryProvider<_Node>>> providers;

    /**
     * \brief Fill the providers vector until given level is reached
     *
     * \param level Tree level
     */
    void add_until_level(const std::size_t level) noexcept {
        const std::size_t provider_count = this->providers.size();
        if(level < provider_count) {
            return;
        }
        this->providers.resize(level+1);
        for(std::size_t l = provider_count; l < level+1; ++l) {
            this->providers[l].emplace_back(provider_size_policy(level, 0));
        }
    }

public:

    /**
     * \brief Provide uninitialised memory for given node count and tree level
     *
     * \param level Node level in the tree
     * \param count Node count
     *
     * \return A pointer to the first node. Memory is uninitialised.
     */
    _Node* provide(const std::size_t level, const std::size_t count) noexcept {
        // Check that a provider for given level exists
        if(this->providers.size() < level+1) {
            this->add_until_level(level+1);
        }

        // Add a provider when needed
        if(! this->providers[level].back().can_provide(count)) {
            this->providers[level].emplace_back(provider_size_policy(level, count));
        }

        return this->providers[level].back().provide(count);
    }


};


#endif /* UNINITNODEMEMORYMANAGER_HPP */
