#ifndef SCALFMM_TREE_HPP_
#define SCALFMM_TREE_HPP_

#include <unordered_set>
#include <functional>

#include "FBox.hpp"
#include "FNode.hpp"
#include "FNodeIteratorBox.hpp"
#include "FInOrderNodeIterator.hpp"
#include "FPrePostOrderNodeIterator.hpp"
#include "UninitNodeMemoryManager.hpp"

#include "inria/meta.hpp"

/**
 * \brief Extract the particle type from the container
 *
 * \tparam Range Container type
 */
template<class Range, class = void>
struct particle_type {
    /// Particle type
    using type = typename Range::value_type;
};

/**
 * \brief Extract the particle type from the container
 *
 * \tparam Range Container type
 */
template<class Range>
struct particle_type<Range, inria::void_t<typename Range::particle_t>> {
    /// Particle type
    using type = typename Range::particle_t;
};


/**
 * \brief Adaptive FMM tree
 *
 * \author Quentin Khan <quentin.khan@inria.fr>
 *
 * Implements the adaptive tree to be used with an adaptive FMM algorithm.
 *
 * \tparam _ParticleContainer The class use to store particles. Expected to
 * expose a value_type type definition that exposes the particle type.
 * \tparam _NodeData The cell data class that is used by the FMM kernel.
 */
template<
    class _ParticleContainer,
    class _NodeData >
class FTree {
public:
    /// Internal node structure
    using node_t = FNode<FTree, _ParticleContainer, _NodeData> ;
    /// Particle container type
    using particle_container_t = _ParticleContainer;
    /// Particle type
    using particle_t = typename particle_type<_ParticleContainer>::type;
    /// Floating point numbers type
    using FReal = typename particle_t::FReal;
    /// Particle position type
    using position_t = typename particle_t::position_t;
    /// Space dimension count
    constexpr static std::size_t Dim = position_t::Dim;
    /// Box type use to slice space
    using box_t = FBox<position_t>;
    /// List type to store leaves
    using leaf_list_t = std::unordered_set<node_t*>;
    /// In order tree traversal mock container type
    using in_order_walk_t   = FNodeIteratorBox<node_t, FInOrderNodeIterator>;
    /// Pre-order tree traversal mock container type
    using pre_order_walk_t  = FNodeIteratorBox<node_t, FPreOrderNodeIterator>;
    /// Post-order tree traversal mock container type
    using post_order_walk_t = FNodeIteratorBox<node_t, FPostOrderNodeIterator>;

    using LeafClass = node_t;

    friend node_t;

    /// Level storage type
    using level_t = typename node_t::level_t;
    /// Morton index storage type
    using morton_index_t = typename node_t::morton_index_t;

    /// Print particles when printing tree or not
    bool print_particles = false;
private:
    /**
     * \brief Tree maximum heigh
     * \note Current morton index implementation prevents heigths greater than
     * 20.
     */
    level_t _max_height = 20;
    /// Maximum particle per leaf density
    std::size_t _leaf_max_particle_count = 50;
    /// Tree space bounding box
    box_t _box;
    /// Tree root node
    node_t* _root = nullptr;
    /// Tree leaf list
    leaf_list_t _leaves;
    /// Tree height
    level_t _height = 0;

    UninitNodeMemoryManager<node_t> node_memory_manager;

public:
    /**
     * \brief Swaps two trees
     */
    void swap(FTree& tree) {
        using std::swap;
        swap(node_memory_manager, tree.node_memory_manager);
        swap(_max_height, tree._max_height);
        swap(_box, tree._box);
        swap(_root, tree._root);
        swap(_leaves, tree._leaves);
        swap(_leaf_max_particle_count, tree._leaf_max_particle_count);

        if(tree.root() && root()) {
            tree.root()->setTree(&(root()->getTree()));
        }
        if(root()) {
            root()->setTree(this);
        }
    }

    /**
     * \brief Builds a tree
     *
     * \param box Tree bounding box
     */
    FTree(box_t box_) : _box(box_) {
        _root = this->node_memory_manager.provide(0,1);
        new(this->_root) node_t(this);
    }

    /**
     * \brief Deleted copy constructor
     */
    FTree(FTree&) = delete;

    /**
     * \brief Move constructor
     */
    FTree(FTree&& tree) {
        swap(tree);
    }

    /**
     * \brief Deleted copy-assignment operator
     */
    FTree& operator=(FTree&) = delete;

    /**
     * \brief Move assignment
     */
    FTree& operator=(FTree&& tree) {
        swap(tree);
    }

    /**
     * \brief Destructor
     */
    ~FTree() {
        _root->~node_t();
    }

    /**
     * \brief Tree height
     * \return The current tree height
     *
     * \b Complexity \f$O(1)\f$
     */
    level_t height() const noexcept {
        return _height;
    }

    level_t getHeight() const noexcept {
        return _height;
    }

private:
    /// Use internally to set height
    void set_height(level_t new_height) {
        _height = new_height;
    }
public:

    /**
     * \brief Maximum height accessor
     * \return Tree maximum height
     */
    level_t max_height() const {
        return _max_height;
    }

    /**
     * \brief Maximum height setter
     * \param height New maximum height
     */
    void max_height(const level_t& new_height) {
        _max_height = new_height;
    }

    /**
     * \brief Maximum leaf particle density accessor
     * \return Tree maximum leaf particle density
     */
    std::size_t leaf_max_particle_count() const {
        return _leaf_max_particle_count;
    }
    /**
     * \brief Maximum leaf particle density accessor
     * \param count New maximum density
     */
    void leaf_max_particle_count(std::size_t count) {
        _leaf_max_particle_count = count;
    }

    /**
     * \brief Tree root accessor
     * \return Root pointer
     */
    node_t* root() {
        return _root;
    }

    /**
     * \brief Tree root const accessor
     * \return Root pointer
     */
    const node_t* root() const {
        return _root;
    }

    /**
     * \brief Bounding box accessor
     * \return Tree bounding box
     */
    const box_t& box() const {
        return _box;
    }

    position_t getBoxCenter() const {
        return _box.center();
    }

    FReal getBoxWidth() const {
        return _box.width(0);
    }

    /**
     * \brief Leaf list accessor
     * \return Tree leaf list
     */
    leaf_list_t& leaves() {
        return _leaves;
    }

    /**
     * \brief Leaf list const accessor
     * \return Tree leaf list
     */
    const leaf_list_t& leaves() const {
        return _leaves;
    }

    /**
     * \brief Builds a ghost leaf at specified location
     *
     * Builds the minimal subtree needed so that the leaf at given level with
     * given morton index exists. Nodes created to reach the given leaf are
     * marked as owned by given owner.
     *
     * \param level Leaf level
     * \param m_idx Leaf morton index
     * \param owner_id Id of process owning the leaf
     */
    void build_ghost_leaf(level_t l, morton_index_t m, unsigned owner_id) {
        this->root()->build_ghost_leaf(l, m, owner_id);
    }

    /**
     * \brief Builds a ghost leaf at specified location
     *
     * Builds the minimal subtree needed so that the leaf at given level with
     * given morton index exists. Nodes created to reach the given leaf are
     * marked as owned by given owner.
     *
     * \param node_info Node information
     * \param owner_id  Id of process owning the leaf
     *
     * \tparam NodeInfo Node information type.
     *
     * The node_info object must have two free functions defined:
     *
     *   - `level(node_info)`
     *   - `morton_index(node_info)`
     */
    template<class NodeInfo>
    void build_ghost_leaf(const NodeInfo& node_info, unsigned owner_id) {
        level_t l = static_cast<level_t>(level(node_info));
        morton_index_t m = static_cast<morton_index_t>(morton_index(node_info));
        this->root()->build_ghost_leaf(l, m, owner_id);
    }

    /**
     * \brief Builds ghost leaves from a range of locations
     *
     * Builds the minimal subtree needed so that the leaf at given level with
     * given morton index exists. Nodes created to reach the given leaf are
     * marked as owned by given owner.
     *
     * \param r Node information range
     * \param owner_id  Id of process owning the leaf
     *
     * \tparam NodeInfoRange Range of node information type.
     *
     * The NodeInfo objects must have two free functions defined:
     *
     *   - `level(node_info)`
     *   - `morton_index(node_info)`
     */
    template<class NodeInfoRange>
    void build_ghost_subtree(const NodeInfoRange& r, unsigned owner_id) {
        for(const auto& node_info : r) {
            this->build_ghost_leaf(node_info, owner_id);
        }
    }


    /**
     * \brief In order walk accessor
     */
    in_order_walk_t in_order_walk() {
        return in_order_walk_t(root());
    }

    /**
     * \brief Pre-order walk accessor
     */
    pre_order_walk_t pre_order_walk() {
        return pre_order_walk_t(root());
    }

    /**
     * \brief Post-order walk accessor
     */
    post_order_walk_t post_order_walk() {
        return post_order_walk_t(root());
    }

    /**
     * \brief Proxy call for FNode#insert applied to #root
     */
    void insert(const particle_t& p) {
        root()->insert(p);
    }

    template<typename... Args>
    void insert(const position_t& position, const Args&... args) {
        root()->insert(position, args...);
    }

    /**
     * \brief Proxy call for FNode#extract appliced to #root
     */
    void extract(const particle_t& p) {
        root()->extract(p);
    }

    /**
     * \brief Updates the underlying graph after moving some particles
     *
     * After an FMM run, particles may have moved out of their leaf's bounding
     * box. This extracts the particles from the leaf and inserts them again at
     * the right place.
     *
     */
    void reshape() {
        set_height(0);
        root()->reshape();
    }

    friend std::ostream& operator<<(std::ostream& os, const FTree& tree) {
        os << *(tree.root());
        return os;
    }


    // Sadly, some ugly functions are needed for compatibility with the rest of
    // the lib...


    template<typename F>
    void forEachLeaf(F&& func) {
        for(auto& leaf : this->leaves()) {
            func(leaf);
        }
    }

    template<typename F>
    void forEachLeaf(F&& func) const {
        for(auto& leaf : this->leaves()) {
            func(leaf);
        }
    }


};
#endif
