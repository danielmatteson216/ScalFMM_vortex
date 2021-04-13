#ifndef _SCALFMM____ORDER_NODE_ITERATOR_HPP_
#define _SCALFMM____ORDER_NODE_ITERATOR_HPP_

#include "FNodeIterator.hpp"

#define NAMESPACE_START namespace scalfmm { namespace details { namespace adaptive_iterator
#define NAMESPACE_END } }

namespace scalfmm {
    namespace tests {
        struct test_PostOrderNodeIterator;
        struct test_PreOrderNodeIterator;
    }
}


namespace scalfmm {
namespace details {
namespace adaptive_iterator {

    enum class NodeIteratorOrder {preorder, postorder};

    template<class Node, NodeIteratorOrder it_type>
    class __OrderNodeIterator : public FNodeIterator<__OrderNodeIterator<Node, it_type>, Node> {

        friend scalfmm::tests::test_PostOrderNodeIterator;
        friend scalfmm::tests::test_PreOrderNodeIterator;

        using base_t = FNodeIterator<__OrderNodeIterator<Node, it_type>, Node>;
        using Direction = FNodeIteratorDirection;

    public:

        friend base_t;
        using base_t::base_t;
        using base_t::operator++;
        using base_t::operator--;

        __OrderNodeIterator& operator++() {
            if(it_type == NodeIteratorOrder::preorder) {
                return move_iterator_prefix<Direction::forwards>();
            } else {
                return move_iterator_postfix<Direction::forwards>();
            }
        }

        __OrderNodeIterator& operator--() {
            if(it_type == NodeIteratorOrder::preorder) {
                return move_iterator_postfix<Direction::backwards>();
            } else {
                return move_iterator_prefix<Direction::backwards>();
            }
        }


    protected:


        /** Moves the cursor to the first node
         *
         * Resets the iterator to a valid state.
         */
        template<Direction dir>
        __OrderNodeIterator& move_to_boundary() {

            base_t::template goto_root<dir>();

            if( (it_type == NodeIteratorOrder::preorder
                 && dir == Direction::backwards)
                || (it_type == NodeIteratorOrder::postorder
                    && dir == Direction::forwards) )
            {
                while(! this->_cur->is_leaf()) {
                    base_t::template move_down<dir>();
                }
            }
            return *this;
        }


        ////////////// POSTORDER TRAVERSAL /////////////////
        template<Direction dir>
        __OrderNodeIterator& move_iterator_postfix() {
            if(base_t::template is_other_end<dir>()) {
                return move_to_boundary<dir>();
            }
            if(base_t::template is_end<dir>()) {
                return *this;
            }
            return next_postfix<dir>();
        }

        template<Direction dir>
        __OrderNodeIterator& next_postfix() {
            // Leave current leaf
            if(!this->_cur_state.empty()) {
                base_t::move_up();
            }
            // If we visited all the children of the node, we move up
            if(!this->_cur_state.empty() && base_t::template at_last_child<dir>()) {
                if(dir == Direction::backwards) {
                    --(this->_cur_state.back());
                } else {
                    ++(this->_cur_state.back());
                }
                return *this;
            }
            // Moving too far up means we got to the end of the tree
            if(this->_cur_state.empty()) {
                base_t::template set_end<dir>();
                return *this;
            }
            // Once we moved up enough, we move down to the next node
            while(! this->_cur->is_leaf()) {
                base_t::template move_down<dir>();
            }
            return *this;
        }

        ////////////// PREORDER TRAVERSAL //////////////////
        template<Direction dir>
        __OrderNodeIterator& move_iterator_prefix() {
            if(base_t::template is_other_end<dir>()) {
                return move_to_boundary<dir>();
            }
            if(base_t::template is_end<dir>()) {
                return *this;
            }
            if(this->_cur->is_leaf()) {
                return leaf_next_prefix<dir>();
            } else {
                return node_next_prefix<dir>();
            }
        }

        template<Direction dir>
        __OrderNodeIterator& leaf_next_prefix() {
            // Leave current leaf
            if(!this->_cur_state.empty()) {
                base_t::move_up();
            }
            // If we visited all the children of the node, we move up
            while(!this->_cur_state.empty() && base_t::template at_last_child<dir>()) {
                base_t::move_up();
            }
            // Moving too far up means we got to the end of the tree
            if(this->_cur_state.empty()) {
                base_t::template set_end<dir>();
                return *this;
            }
            // Once we moved up enough, we move down to the next node
            base_t::template move_down<dir>();
            return *this;
        }

        template<Direction dir>
        __OrderNodeIterator& node_next_prefix() {
            base_t::template move_down<dir>();
            return *this;
        }


        ////////////////////////////////////////////////////////////
    };

} } }


#undef NAMESPACE_START
#undef NAMESPACE_END


template<class Node>
using FPreOrderNodeIterator = scalfmm::details::adaptive_iterator::
    __OrderNodeIterator<Node,scalfmm::details::adaptive_iterator::
                        NodeIteratorOrder::preorder>;
template<class Node>
using FPostOrderNodeIterator = scalfmm::details::adaptive_iterator::
    __OrderNodeIterator<Node,scalfmm::details::adaptive_iterator::
                        NodeIteratorOrder::postorder>;


#endif
