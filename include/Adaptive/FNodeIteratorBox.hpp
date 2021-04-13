#ifndef _SCALFMM_NODE_ITERATOR_BOX_
#define _SCALFMM_NODE_ITERATOR_BOX_

#include <iterator>

#include "FNodeIterator.hpp"

template <class Node, template <class> class NodeIterator>
class FNodeIteratorBox {
    using node_t = Node;
public:

    using value_type      = node_t;
    using reference       = node_t&;
    using const_reference = const node_t&;
    using iterator        = NodeIterator<node_t>;
    using const_iterator  = NodeIterator<const node_t>;
    using difference_type = typename iterator::difference_type;
    using size_type = std::size_t;


    node_t* _root = nullptr;

    FNodeIteratorBox(node_t* root) : _root(root) { }

    iterator begin() {
        return iterator(_root);
    }

    iterator end() {
        return iterator(_root, scalfmm::FNodeIteratorPosition::end);
    }

    const_iterator begin() const {
        return const_iterator(_root, scalfmm::FNodeIteratorPosition::begin);
    }

    const_iterator end() const {
        return const_iterator(_root, scalfmm::FNodeIteratorPosition::end);
    }

    const_iterator cbegin() const {
        return this->begin();
    }

    const_iterator cend() const {
        return this->cend();
    }

};

#endif
