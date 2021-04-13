#ifndef _SCALFMM_IN_ORDER_NODE_ITERATOR_HPP_
#define _SCALFMM_IN_ORDER_NODE_ITERATOR_HPP_

#include "FNodeIterator.hpp"

namespace scalfmm {
    namespace tests {
        struct test_InOrderNodeIterator;
    }
}

template<class Node>
class FInOrderNodeIterator : public FNodeIterator<FInOrderNodeIterator<Node>, Node> {

    using base_t = FNodeIterator<FInOrderNodeIterator<Node>, Node>;
    using Direction = scalfmm::FNodeIteratorDirection;

    friend scalfmm::tests::test_InOrderNodeIterator;
    friend base_t;

public:

    using base_t::base_t;
    using base_t::operator++;
    using base_t::operator--;

    FInOrderNodeIterator& operator++() {
        return move_iterator<Direction::forwards>();
    }

    FInOrderNodeIterator& operator--() {
        return move_iterator<Direction::backwards>();
    }


private:

    /** Indicates the index of the middle child */
    template<Direction dir>
    typename base_t::index_t mid_idx() const {
        if(dir == Direction::forwards) { // Compile time optimisation will remove unused branch
            return (base_t::child_count/2)-1;
        } else {
            return (base_t::child_count/2);
        }
    }

    /** Moves the cursor to the first node
     *
     * Resets the iterator to a valid state.
     */
    template<Direction dir>
    FInOrderNodeIterator& move_to_boundary() {

        base_t::template goto_root<dir>();

        while(! this->_cur->is_leaf()) {
            base_t::template move_down<dir>();
        }
        return *this;
    }

    template<Direction dir>
    FInOrderNodeIterator& leaf_next() {
        // Leave current leaf
        if(!this->_cur_state.empty()) {
            base_t::move_up();
        }
        // If we visited half the parent's children, stay on the parent
        if (!this->_cur_state.empty() && mid_idx<dir>() == this->_cur_state.back()) {
            return *this;
        }
        // If we visited all the children of the node, we move up
        while(!this->_cur_state.empty() && base_t::template at_last_child<dir>()) {
            base_t::move_up();

            if (!this->_cur_state.empty() && mid_idx<dir>() == this->_cur_state.back()) {
                return *this;
            }
        }
        // Moving too far up means we got to the end of the tree
        if(this->_cur_state.empty()) {
            base_t::template set_end<dir>();
            return *this;
        }
        // Once we moved up enough, we move down to the next leaf
        while(! this->_cur->is_leaf()) {
            base_t::template move_down<dir>();
        }

        return *this;
    }

    template<Direction dir>
    FInOrderNodeIterator& node_next() {
        // If we are on an internal node, we move down to the children
        if(dir == Direction::forwards) { // Compile time optimisation will remove unused branch
            this->_cur_state.back() = static_cast<char>(mid_idx<dir>()+1);
        } else {
            this->_cur_state.back() = static_cast<char>(mid_idx<dir>()-1);
        }
        this->_cur = this->_cur->getChild((std::size_t)this->_cur_state.back());
        this->_cur_state.push_back(base_t::template first_idx<dir>());

        while(! this->_cur->is_leaf()) {
            base_t::template move_down<dir>();
        }
        return *this;
    }


    template<Direction dir>
    FInOrderNodeIterator& move_iterator() {
        // If before the beginning, move to the first node
        if(base_t::template is_other_end<dir>()) {
            return move_to_boundary<dir>();
        }
        // Don't move if already past the end
        if(base_t::template is_end<dir>()) {
            return *this;
        }
        if(this->_cur->is_leaf()) {
            return leaf_next<dir>();
        } else {
            return node_next<dir>();
        }
    }
};

#endif
