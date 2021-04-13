#ifndef _SCALFMM_NODE_ITERATOR_HPP_
#define _SCALFMM_NODE_ITERATOR_HPP_

#include <vector>

#include "Utils/FConstFuncs.hpp"

namespace scalfmm {
namespace tests {
struct test_NodeIterator;
}

enum class FNodeIteratorPosition {begin, end, reverse_end};
enum class FNodeIteratorDirection {forwards, backwards};

}

template<class Derived, class Node>
class FNodeIterator {
public:

    friend struct scalfmm::tests::test_NodeIterator;

    using node_t = Node;

    // Iterator types
    using value_type = node_t;
    using difference_type = std::ptrdiff_t;
    using pointer = node_t*;
    using reference = node_t&;
    using iterator_category = std::bidirectional_iterator_tag;

    constexpr static std::size_t child_count = node_t::child_count;

protected:
    using state_container = std::vector<char>;
    using index_t = typename state_container::value_type;


    /// Root node for the graph walk
    /** \warning This may not be the first node of the walk. */
    node_t* _root = nullptr;
    /// Current node for the graph walk
    node_t* _cur  = nullptr;
    /** Current state in the walk
     * Contains the indexes of the next children to visit per level
     * in the currently visited tree branch.
     */
    state_container _cur_state;
    /// Special past_the_end iterator status
    bool _is_end = false;
    /// Special past_the_reverse_end iterator status
    bool _is_rend = false;

public:
    /// Main constructor
    /**
     * \param root The graph root node before walking accross it.
     *
     * \warning This constructor must be completed by the derived class.
     * It does not set #_cur nor #_cur_state to the right value. */
    FNodeIterator(node_t* root, scalfmm::FNodeIteratorPosition p = scalfmm::FNodeIteratorPosition::begin) :
        _root(root),
        _is_end(p == scalfmm::FNodeIteratorPosition::end),
        _is_rend(p == scalfmm::FNodeIteratorPosition::reverse_end)
    {
        if(!_is_end && ! _is_rend) {
            static_cast<Derived*>(this)->template move_to_boundary<scalfmm::FNodeIteratorDirection::forwards>();
        }
    }

    /// Default constructor
    /** Such a constructed iterator is not in a valid state. */
    FNodeIterator() = default;
    /// Copy constructor
    FNodeIterator(const FNodeIterator& other) = default;
    /// Copy assignment
    FNodeIterator& operator=(const FNodeIterator& other) = default;
    /// Move constructor
    FNodeIterator(FNodeIterator&& other) = default;
    /// Move assignment
    FNodeIterator& operator=(FNodeIterator&& other) = default;

    /// Destructor
    ~FNodeIterator() = default;

    /// Swap operation
    void swap(FNodeIterator&& other) noexcept {
        #define _MY_SWAP_(a) {auto tmp = a; a = other.a; other.a = tmp;}
        _MY_SWAP_(_root);
        _MY_SWAP_(_cur);
        _MY_SWAP_(_is_end);
        _MY_SWAP_(_is_rend);
        _cur_state.swap(other._cur_state);
        #undef _MY_SWAP_
    }


    /// Dereference operator
    node_t& operator*() {
        return *_cur;
    }
    /// Dereference const operator
    const node_t& operator*() const {
        return *_cur;
    }

    /// Pointer dereference operator
    node_t* operator->() {
        return _cur;
    }
    /// Pointer const dereference operator
    const node_t* operator->() const {
        return _cur;
    }

    /// Prefix increment
    //virtual FNodeIterator& operator++() = 0;

    /// Postfix increment
    Derived operator++(int) {
        Derived it(static_cast<Derived&>(*this));
        ++static_cast<Derived&>(*this);
        return it;
    }

    /// Prefix decrement
    //virtual FNodeIterator& operator--() = 0;

    /// Postfix decrement
    Derived operator--(int) {
        Derived it(static_cast<Derived&>(*this));
        --static_cast<Derived&>(*this);
        return it;
    }

    /// Equality operator
    friend bool operator==(const FNodeIterator& lhs, const FNodeIterator& rhs) {
        #define _TEST_EQ(name) (lhs.name == rhs.name)
        return
            _TEST_EQ(_root) &&
            _TEST_EQ(_cur) &&
            _TEST_EQ(_is_end) &&
            _TEST_EQ(_is_rend);
        #undef _TEST_EQ
    }

    /// Inequality operator
    friend bool operator!=(const FNodeIterator& lhs, const FNodeIterator& rhs) {
        return ! (lhs == rhs);
    }

    friend Derived operator+=(Derived lhs, difference_type n) {
        if(n < 0) {
            for(; n != 0; ++n, --lhs);
        }
        if(n > 0) {
            for(; n != 0; --n, ++lhs);
        }
        return lhs;
    }

    friend Derived operator-=(const Derived& lhs, difference_type n) {
        return operator+=(lhs, -n);
    }

    friend Derived operator+(const Derived& lhs, difference_type n) {
        return lhs += n;
    }

    friend Derived operator-(const Derived& lhs, difference_type n) {
        return operator+=(lhs, -n);
    }

protected:

    // Methods used in all subclasses

    template<scalfmm::FNodeIteratorDirection dir>
    void goto_root() {
        this->_is_rend = false;
        this->_is_end = false;
        this->_cur = this->_root;
        this->_cur_state.clear();
        this->_cur_state.push_back(first_idx<dir>());
    }

    /** Returns the end flag corresponding to the direction */
    template<scalfmm::FNodeIteratorDirection dir>
    bool is_end() {
        return dir == scalfmm::FNodeIteratorDirection::forwards ? this->_is_end : this->_is_rend;
    }

    /** Returns the end flag corresponding to the oposite direction */
    template<scalfmm::FNodeIteratorDirection dir>
    bool is_other_end() {
        return dir == scalfmm::FNodeIteratorDirection::forwards ? this->_is_rend : this->_is_end;
    }

    /** Sets the end flag according to the direction */
    template<scalfmm::FNodeIteratorDirection dir>
    void set_end() {
        if(dir == scalfmm::FNodeIteratorDirection::forwards) {
            this->_is_end = true;
        } else {
            this->_is_rend = true;
        }
    }

    /** Returns the first index to insert in #_cur_state */
    template<scalfmm::FNodeIteratorDirection dir>
    index_t first_idx() const {
        return dir == scalfmm::FNodeIteratorDirection::forwards ? -1 : (index_t)Node::child_count;
    }


    /** Indicates whether the last child of a node has been visited */
    template<scalfmm::FNodeIteratorDirection dir>
    bool at_last_child() const {
        if(dir == scalfmm::FNodeIteratorDirection::forwards) { // Compile time optimisation will remove unused branch
            return Node::child_count-1 == this->_cur_state.back();
        } else {
            return 0 == this->_cur_state.back();
        }
    }

    /** Moves the cursor to its parent
     *
     * \warning Undefined behaviour when the cursor points to a leaf.
     */
    void move_up() {
        this->_cur_state.pop_back();
        this->_cur = this->_cur->getParent();
    }

    /** Moves the cursor to the next child to visit
     *
     * \warning Undefined behaviour when the cursor points to a leaf.
     */
    template<scalfmm::FNodeIteratorDirection dir>
    void move_down() {
        if(dir == scalfmm::FNodeIteratorDirection::forwards) { // Compile time optimisation will remove unused branch
            ++(this->_cur_state.back());
        } else {
            --(this->_cur_state.back());
        }
        this->_cur = this->_cur->getChild((std::size_t) this->_cur_state.back());
        this->_cur_state.push_back(first_idx<dir>());
    }



};

template<class Derived, class Node>
void swap(FNodeIterator<Derived, Node>& it1, FNodeIterator<Derived, Node>& it2) {
    it1.swap(it2);
}

#endif
