#ifndef SCALFMM_NODE_HPP_
#define SCALFMM_NODE_HPP_

#include <algorithm>
#include <cmath>
#include <cassert>
#include <functional>
#include <memory>
#include <ostream>
#include <iomanip>
#include <stdexcept>
#include <string>
#include <vector>
#include <unordered_set>

#include "Utils/FOstreamTuple.hpp"
#include "inria/ostream_joiner.hpp"
#include "inria/detection_idiom.hpp"

#include "FBox.hpp"
#include "Utils/FPoint.hpp"
#include "Utils/FConstFuncs.hpp"
#include "Utils/FOstreamTuple.hpp"

#include "Components/FSymbolicData.hpp"
#include "Components/FBasicCell.hpp"

namespace scalfmm {

namespace tests {
struct test_Node;
struct test_NodeIterator;
struct test_InOrderNodeIterator;
struct test_PreOrderNodeIterator;
struct test_PostOrderNodeIterator;
} // close namespace [scalfmm]::tests


namespace meta {
/**
 * \brief Sink in type for SFINAE purpose
 *
 * Derives from std::true_type.
 *
 * \tparam Ts Type parameter pack
 */
template<typename... Ts>
struct exist_t : std::true_type {};

/**
 * \brief Compile check that `T::push(Args...)` exists
 *
 * This does not check the return type
 */
template<typename T, typename... Args>
class has_push {
    /**
     * \brief Overload found if `T::push(Args...)` exists
     *
     * If `T::push(Args...)` exists, this overload is found when
     * computing #value at compile time.
     *
     * \tparam U used as T
     * \tparam unnamed Default parameter dependent on U, fails
     * compilation if `U::push(Args...)` does not exist.
     */
    template<typename U,
             typename std::enable_if<
                 exist_t<decltype(std::declval<U>().push(std::declval<Args>()...))>::value
                 >::type* = nullptr>
    static std::true_type get(U);

    /// SFINAE fallback
    static std::false_type get(...);
public:
    /// True if `T::push(Args...)` exists, false otherwise
    constexpr static const bool value = decltype(get(std::declval<T>()))::value;
};


/**
 * \brief Compile check that `T::push_back(Args...)` exists
 *
 * This does not check the return type
 */
template<typename T, typename... Args>
class has_push_back {
    /**
     * \brief Overload found if `T::push_back(Args...)` exists
     *
     * If `T::push_back(Args...)` exists, this overload is found when
     * computing #value at compile time.
     *
     * \tparam U used as T
     * \tparam unnamed Default parameter dependent on U, fails
     * compilation if `U::push_back(Args...)` does not exist.
     */
    template<typename U,
             typename std::enable_if<
                 exist_t<decltype(std::declval<U>().push_back(std::declval<Args>()...))>::value
                 >::type* = nullptr>
    static std::true_type get(U);

    /// SFINAE fallback
    static std::false_type get(...);
public:
    /// True if `T::push_back(Args...)` exists, false otherwise
    constexpr static const bool value = decltype(get(std::declval<T>()))::value;
};

} // close [scalfmm]::meta


namespace sfinae {

template<bool value, template<class...> class Checker, typename... Args>
using use_if = typename std::enable_if<value == Checker<Args...>::value >::type*;


/**
 * \brief Exists (as void*) if T derives from Base
 *
 * Type alias for use in SFINAE. Is a void* if Derived is derived from
 * Base, does not exist otherwise.
 *
 * \tparam Base Wanted base type
 * \tparam T Type to check against Base for derivation
 */
template<typename Base, typename T>
using derived_from = use_if<true, std::is_base_of, Base,T >;
}

struct fmt {
    enum : long {TERM, OBJ};

    /// Get the stream flag index for the output format
    static const int& node_os_format_id() {
        static const int i = std::ios_base::xalloc();
        return i;
    }

    /// Get the stream flag index for the particle output in terminal format
    static const int& node_os_particle_id() {
        static const int i = std::ios_base::xalloc();
        return i;
    }

    static std::ostream& term(std::ostream& os) {
        os.iword(node_os_format_id()) = TERM;
        return os;
    };

    static std::ostream& obj(std::ostream& os) {
        os.iword(node_os_format_id()) = OBJ;
        return os;
    };

    static std::ostream& parts(std::ostream& os) {
        os.iword(node_os_particle_id()) = true;
        return os;
    };

    static std::ostream& no_parts(std::ostream& os) {
        os.iword(node_os_particle_id()) = false;
        return os;
    };


};

}

struct NodeEmptyData {};



/**
 * \brief Tree node implementation
 */
template<class _Tree, class _ParticleContainer, class NodeData>
class FNode {
public:
    /// Tree type this class belongs to
    using tree_t = _Tree;
    /// Type used to represent real numbers
    using FReal = typename tree_t::FReal;
    /// Space dimension count
    constexpr static const std::size_t Dim = _Tree::Dim;
    /// Child count if the node is not a leaf
    constexpr static const std::size_t child_count = 1 << Dim;

    /// Node position type
    using position_t = typename tree_t::position_t;
    /// Node bounding box type
    using box_t =  typename tree_t::box_t;
    /// Interaction lists type
    using interaction_list_t = std::unordered_set<FNode*>;
    /// Children array type
    using child_node_array_t = std::array<FNode*, child_count>;
    /// Particle container type
    using particle_container_t = _ParticleContainer;
    /**
     * \brief Particle type
     * The particle must satisfy the following conditions:
     *   - Default constructible
     *   - Constructible from particle_container_t::value_type
     *   - position_t position() const method must exist
     *
     */
    using particle_t = typename tree_t::particle_t;
    /// Node data structure
    using data_t = NodeData;

    /**
     * \brief Node structural data
     * Keeps data about the node that may be read by kernels or algorithms.
     */
    using symbolic_data_t
    = typename extract_symbolic_data_t_or_fallback_to_default<FSymbolicData,data_t>::type;

    static_assert(models_symbolic_data<symbolic_data_t>::value,
                  "The symbolic_data_t type does not model the required symbolic data interface.");

    /// Level type alias
    using level_t = typename std::remove_reference<
        decltype(std::declval<symbolic_data_t>().getLevel())>::type;

    /// Morton index type alias
    using morton_index_t = typename std::remove_reference<
        decltype(std::declval<symbolic_data_t>().getMortonIndex())>::type;

    constexpr static int not_ghost_tag = ~0;

private:

    friend struct scalfmm::tests::test_Node;
    friend struct scalfmm::tests::test_NodeIterator;
    friend struct scalfmm::tests::test_InOrderNodeIterator;
    friend struct scalfmm::tests::test_PreOrderNodeIterator;
    friend struct scalfmm::tests::test_PostOrderNodeIterator;
    friend tree_t;

    /// Children array, filled with nullptr is Node is a leaf
    child_node_array_t _children{{}};
    /// Node parent, nullptr if the node is a root
    FNode* _parent = nullptr;
    /// Node spatial bounding box
    box_t _box;
    /// Particle container
    particle_container_t _p_container;
    /// Node data
    data_t _data = data_t();

    /// Tree the Node belongs to
    tree_t* _tree;

    /// Indicates whether node is a leaf
    bool _is_leaf = true;
    /// Owner process when node is a ghost
    int ghost_owner_ = ~0;

    /// Node data that may be of use to algorithms and kernels
    symbolic_data_t _symbolic_data;

public:
    /// Near-field leaves interaction list
    interaction_list_t U;
    /// Mid-field node interaction list
    interaction_list_t V;
    /// Near-field node interaction list
    interaction_list_t W;
    /// Mid-field leaf interaction list
    interaction_list_t X;

    /**
     * \brief Constructor called from parent
     *
     * This constructor is meant to be called from a parent node.
     *
     * \param parent The parent node
     * \param child_index The index of this node in the parent children array
     */
    FNode(FNode& parent, const std::size_t& child_index) :
        _parent(&parent),
        _box   (parent.getBox().center(), parent.getBox().corner(child_index) ),
        _tree  (parent._tree ),
        _symbolic_data{static_cast<int>(parent.getLevel()+1), (parent.getIndex() << Dim) + child_index}
    {
        if (child_index >= child_count) {
            throw std::invalid_argument(std::string("Wrong child index in node contructor: got ")
                                        + std::to_string(child_index)
                                        + std::string(", expected at most ")
                                        + std::to_string(Dim)
                                        + ".");
        }

        this->common_init();
    }

    /**
     * \brief Root constructor called from tree
     *
     * This constructor must be called when building the tree root.
     */
    FNode(tree_t* tree) :
        _box(tree->box()),
        _tree(tree),
        _symbolic_data{}
    {
        tree->leaves().insert(this);
        // `this` belongs to its own U list, not done in the other constructors
        // because managed by split or fuse
        this->U.insert(this);
        this->common_init();
    }

    /**
     * \brief Deleted default constructor
     */
    FNode() = delete;

    /**
     * \brief Deleted copy constructor
     */
    FNode(const FNode& other) = delete;

    /**
     * \brief Deleted copy operator
     */
    FNode& operator=(const FNode& other) = delete;

    /**
     * \brief Deleted move constructor
     */
    FNode(FNode&& other) = delete;

    /**
     * \brief Deleted move operator
     */
    FNode& operator=(FNode&& other) = delete;

    /**
     * \brief Deleted destructor
     */
    ~FNode() {
        this->delete_children();
    }


    /// Data accessor
    data_t* getData() noexcept {
        return &_data;
    }
    /// Data const accessor
    const data_t* getData() const noexcept {
        return &_data;
    }

    /// Symbolic data accessor
    symbolic_data_t& getSymbolicData() noexcept {
        return this->_symbolic_data;
    }
    /// Symbolic data const accessor
    const symbolic_data_t& getSymbolicData() const noexcept {
        return this->_symbolic_data;
    }


    /// Children container accessor
    child_node_array_t& getChildren() noexcept {
        return _children;
    }
    /// Children const container accessor
    const child_node_array_t& getChildren() const noexcept {
        return _children;
    }

    /**
     * \brief Child container accessor
     *
     * \param index Child index
     */
    FNode* getChild(const std::size_t& index) noexcept {
        return getChildren()[index];
    }
    /**
     * \brief Child container const accessor
     *
     * \param index Child index
     */
    const FNode* getChild(const std::size_t& index) const noexcept {
        return getChildren()[index];
    }

    /// Parent accessor
    FNode* getParent() noexcept {
        return _parent;
    }
    /// Parent const accessor
    const FNode* getParent() const noexcept {
        return _parent;
    }

    /// Depth accessor
    [[gnu::deprecated("Use getLevel for interface consistency.")]]
    level_t getDepth() const noexcept {
        return _symbolic_data.getLevel();
    }
    /// Level accessor
    level_t getLevel() const noexcept {
        return _symbolic_data.getLevel();
    }

    /// Morton index accessor
    morton_index_t getIndex() const noexcept {
        return _symbolic_data.m_idx;
    }

    /// Check whether node is a ghost
    bool is_ghost() const {
        return this->ghost_owner_ != not_ghost_tag;
    }

    /// Owner process rank
    int ghost_owner() const {
        return this->ghost_owner_;
    }

    void set_ghost_owner(int owner_id) {
        this->ghost_owner_ = owner_id;
    }

    /// Tree accessor
    tree_t& getTree() noexcept {
        return *_tree;
    }
    /// Tree const accessor
    const tree_t& getTree() const noexcept {
        return *_tree;
    }

    /// Box const accessor
    const box_t& getBox() const noexcept {
        return _box;
    }

    /// Particle container accessor
    particle_container_t* getParticleContainer() noexcept {
        return &_p_container;
    }
    /// Particle container accessor
    const particle_container_t* getParticleContainer() const noexcept {
        return &_p_container;
    }

    /// Particle container accessor
    particle_container_t* getTargets() noexcept {
        return &_p_container;
    }
    /// Particle container accessor
    const particle_container_t* getTargets() const noexcept {
        return &_p_container;
    }

    /// Particle count for the container
    std::size_t getParticleCount() const noexcept {
        return getParticleContainer()->size();
    }

    /**
     * \brief Find out whether this node and the 'other' node are adjacent
     *
     * The nodes are assumed to belong to the same tree.
     *
     * To check whether nodes are adjacent, on each axis, the distance between
     * the nodes' center is compared to the sum of their half diameter. For at
     * least one of the axes, the two must be equal. For the others, the
     * distance must be less than or equal to the sum. This ensures that a node
     * is not adjacent to one of its descendants.
     *
     * \param other The node to test adjacency with.
     *
     * \return true if this FNode and the 'other' FNode are adjacent.
     */
    bool is_adjacent(const FNode& other) const noexcept {
        // Sum of the half side lengh of the two nodes boxes.
        // Boxes are cubes, we only need one side.
        FReal centers_distance = getBox().center()[0] - getBox().c1()[0]
            + other.getBox().center()[0] - other.getBox().c1()[0];
        // Used to check that the other box isn't overlapping with this box
        bool one_axis_is_at_exact_distance = false;

        position_t my_center = getBox().center();
        position_t other_center = other.getBox().center();

        for(std::size_t i = 0; i < Dim; ++i) {
	  FReal distance = std::abs(my_center[i] - other_center[i]);
            if( Ffeq(distance, centers_distance) ) {
                one_axis_is_at_exact_distance = true;
            } else if(distance > centers_distance) {
                return false;
            }
        }

        return one_axis_is_at_exact_distance;
    }

    /**
     * \brief Tests whether this node and the 'other' node are adjacent
     *
     * \return true if the nodes are adjacent.
     */
    bool is_adjacent(const FNode* other) const noexcept {
        if(nullptr == other) {
            return false;
        } else if (this == other){
            return true;
        }

        return is_adjacent(*other);
    }

    /**
     * \brief Tests whether this node is a leaf
     *
     * \return true if this node is a leaf.
     */
    bool is_leaf() const noexcept {
        return _is_leaf;
    }

private:
    /**
     * \brief Push particle in container if it has compatible 'push'
     *
     * \tparam T Container type
     * \tparam scalfmm::sfinae::use_if SFINAE default argument. Compilation of
     * this overload fails if the condition isn't true
     *
     * \param container Container to push into
     * \param p Particle to push
     */
    template<typename T, class... Args,
             scalfmm::sfinae::use_if<true, scalfmm::meta::has_push, T, Args...> = nullptr>
    void particle_push(T& container, const Args&... args) {
        container.push(args...);
    }

    /**
     * \brief Push particle in container if it has compatible `push_back` and
     * no `push`
     *
     * \tparam T Container type
     * \tparam scalfmm::sfinae::use_if SFINAE default argument
     * \tparam scalfmm::sfinae::use_if SFINAE default argument
     *
     * \param container Container to push into
     * \param p Particle to push
     */
    template<typename T,
             scalfmm::sfinae::use_if<true, scalfmm::meta::has_push_back, T, particle_t> = nullptr,
             scalfmm::sfinae::use_if<false, scalfmm::meta::has_push, T, particle_t> = nullptr>
    void particle_push(T& container, const particle_t& p) {
        container.push_back(p);
    }

public:
    /**
     * \brief Inserts a particle in the tree rooted at node
     *
     * Pushes a particle in the node particle container if it is a leaf. If it
     * isn't, the particle is forwarded to the relevant child.
     *
     * \note The `push_back` method of #particle_container_t is used to insert in a leaf
     */
    void insert(const particle_t& p) {
        if(! is_leaf()) {
            std::size_t child_index = box_t::space_filling_curve_t::index(p.position(), getBox().center());
            getChild(child_index)->insert(p);
        } else {
            this->set_ghost_owner(not_ghost_tag);
            this->getParticleContainer()->reserve(this->getTree().leaf_max_particle_count());
            particle_push(*(this->getParticleContainer()), p);
            if(getParticleContainer()->size() > getTree().leaf_max_particle_count()) {
                split();
            }
        }
    }

    /**
     * \brief Inserts a particle in the tree rooted at node
     *
     * Emplaces a particle in the node particle container if it is a leaf. If it
     * isn't, the particle is forwarded to the relevant child.
     *
     * \note The `emplace_back` method of #particle_container_t is used to insert in a leaf
     */
    template<typename... Args>
    void insert(const position_t& position, Args&&... args) {
        if(! is_leaf()) {
            std::size_t child_index = box_t::space_filling_curve_t::index(position, getBox().center());
            getChild(child_index)->insert(position, args...);
        } else {
            this->set_ghost_owner(not_ghost_tag);
            this->getParticleContainer()->reserve(this->getTree().leaf_max_particle_count());
            particle_push(*(this->getParticleContainer()), position, args...);
            if(getParticleContainer()->size() > getTree().leaf_max_particle_count()) {
                split();
            }
        }
    }

    /**
     * \brief Extracts a particle from a leaf
     *
     * If the nodes is not a leaf, does nothing.
     *
     * \return If in a leaf, the corresponding particle, otherwise a default
     * constructed particle_t
     */
    particle_t extract(const std::size_t& idx) {
        if(getParticleContainer()) {
            particle_t p(*(getParticleContainer()->begin() + idx));
            getParticleContainer()->erase(getParticleContainer()->begin() + idx);
            return p;
        }
        return particle_t();
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
    void build_ghost_leaf(level_t level, morton_index_t m_idx, int owner_id) {
        constexpr unsigned mask = ~((~0u) << Dim);
        auto* n = this;

        while(n->getLevel() != level) {
            if(n->is_leaf()) {
                n->split();
                for(FNode* child : n->_children) {
                    child->set_ghost_owner(owner_id);
                }
            }

            auto child_idx = (m_idx >> (level - n->getLevel() - 1)) & mask;
            n = n->getChild(child_idx);
        }

        n->set_ghost_owner(owner_id);
    }


    /**
     * \brief Splits or fuses nodes after tree modifications
     *
     * Iterates over all the sub-tree particles and splits the leafs that
     * hold too many particles or fuses those that hold too few.
     *
     * This may be necessary when the tree parameters have been changed
     * (leaf max particle count modified) or when particles have moved and
     * been reinserted in the tree to be in the right boxes.
     */
    void reshape() {
        if(is_leaf()) {
            if(getParticleCount() > getTree().leaf_max_particle_count()) {
                split();
            } else {
                update_tree_height(); // split will update the tree height otherwise
            }
        } else {
            std::size_t p_count = 0;
            for(auto child : getChildren()) {
                child->reshape();
                p_count += child->getParticleCount();
            }

            if(p_count <= getTree().leaf_max_particle_count()) {
                fuse();
            }
        }
    }

    void split_to_height(const std::size_t new_height) {
        assert(new_height < 30); // Trees cannot go this high, means an underflow happened
        if(new_height <= 1) {
            return;
        }

        if(this->is_leaf()) {
            this->split();
            for(auto& child : this->getChildren()) {
                child->split_to_height(new_height - 1);
            }
        }
    }

    /**
     * \brief Applies a function to the node and it descendants in order
     */
    void for_each_in_order(std::function<void(FNode*)> lambda) {
        std::size_t idx = 0;
        if(! is_leaf()) {
            for(; idx < getChildren().size()/2; ++idx) {
                getChildren()[idx]->for_each_in_order(lambda);
            }
        }
        lambda(this);
        if(! is_leaf()) {
            for(; idx < getChildren().size(); ++idx) {
                getChildren()[idx]->for_each_in_order(lambda);
            }
        }
    }

    /**
     * \brief Applies a function to the node and it descendants post order (children first)
     */
    void for_each_post_order(std::function<void(FNode*)> lambda) {
        if(! is_leaf()) {
            for(std::size_t idx = 0; idx < getChildren().size(); ++idx) {
                getChildren()[idx]->for_each_post_order(lambda);
            }
        }
        lambda(this);
    }

    /**
     * \brief Applies a function to the node and it descendants pre order (parent first)
     */
    void for_each_pre_order(std::function<void(FNode*)> lambda) {
        lambda(this);
        if(! is_leaf()) {
            for(std::size_t idx = 0; idx < getChildren().size(); ++idx) {
                getChildren()[idx]->for_each_pre_order(lambda);
            }
        }
    }

    /**
     * \brief Equality test operator
     */
    bool operator==(const FNode& other) const {
        return other.getParent() == getParent()
            && other.getLevel() == getLevel()
            && other.getIndex() == getIndex();
    }

private:


    /**
     * \brief Initialization of FNode cell data
     *
     * Encloses the init functions that are used to initializes cell data
     * according to their base class.
     */
    struct cell_data_initializer {
        /**
         * \brief Initializer when CellData is derived from FBasicCell
         *
         * \tparam CellData Type of the node data
         * \tparam Unnamed SFINAE type to check that CellData is derived from
         * FBasicCell
         */
        template<typename CellData, scalfmm::sfinae::derived_from<FBasicCell, CellData> = nullptr >
        static void init(CellData* data, FNode* const node) {
            data->setMortonIndex(node->getIndex());
            FTreeCoordinate coord(node->getIndex());
            data->setCoordinate(coord);
            data->setLevel(node->getLevel());
        }

        /**
         * \brief Catch-all initializer, no-op
         */
        static void init(...) {}

    };

    /**
     * \brief Common initialization done by constructors
     */
    void common_init() {
        cell_data_initializer::init(this->getData(), this);
    }

    /**
     * \brief Tree setter
     */
    void setTree(tree_t* t) {
        _tree = t;
        if(! is_leaf()) {
            for(FNode*& child : getChildren()) {
                child->setTree(t);
            }
        }
    }

    /**
     * \brief Updates the tree height when called from a leaf
     *
     * \warning Must be called from a leaf
     */
    void update_tree_height() {
        assert(is_leaf() == true);
        if(getLevel()+1 > getTree().height())
            getTree().set_height(getLevel() + 1);
    }

    /**
     * \brief Creates or allocates the particle container
     */
    void create_particle_container() {
        _p_container.clear();
        _p_container.reserve(this->getTree().leaf_max_particle_count());
    }

    /**
     * \brief Deletes the particle container to save space when it is not needed
     */
    void delete_particle_container() {
        _p_container.clear();
    }

    /**
     * \brief Allocates this node's children
     */
    void create_children() {
        std::size_t idx = 0;
        // Remove this node from tree leaf list
        getTree().leaves().erase(this);
        // Create the children, add them to tree leaf list
        FNode* tmp = this->getTree().node_memory_manager.provide(this->getLevel()+1, child_count);
        for(FNode*& child : getChildren()) {
            child = new(tmp+idx) FNode(*this, idx);
            getTree().leaves().insert(child);
            ++idx;
        }
        // Update tree height from child
        getChild(0)->update_tree_height();
        // Remove leaf status from this node
        _is_leaf = false;
    }

    /*
     * \brief Deletes this node's children
     */
    void delete_children() {
        // Remove children from tree leaf list, free them
        for(FNode*& child : getChildren()) {
            if(child) {
                getTree().leaves().erase(child);
                child->~FNode();
            }
        }
        std::fill_n(this->getChildren().data(), child_count, nullptr);
        // Insert this node in tree leaf list
        getTree().leaves().insert(this);
        // Set leaf status for this node
        _is_leaf = true;
    }

    /**
     * \brief Adds children to a leaf node
     *
     * Adds children to a leaf node and redistributes its particles among
     * the newly created leaves.
     */
    void split() {
        assert(this->is_leaf());

        if(getLevel()+1 > getTree().max_height()) {
            // TODO: log that there were too many particles
            return;
        }

        create_children();

        // Update interaction lists
        this->U.erase(this);
        for(FNode* child : getChildren()) {
            // Children are adjacent to each other
            child->U.insert(getChildren().begin(), getChildren().end());
            // Find where to put U list items in child
            for(FNode* u_item : this->U) {
                if(child->is_adjacent(u_item)) {
                    // Adjacent to child
                    child->U.insert(u_item);
                    u_item->U.insert(child);
                } else if(u_item->getLevel() < child->getLevel()) {
                    // Adjacent to parent (this) but not to child
                    child->X.insert(u_item);
                    u_item->W.insert(child);
                } else { // u_item->getLevel() >= child->getLevel()
                    // Find ancestor of u_item that is adjacent to child
                    while(u_item->getLevel() > child->getLevel()) {
                        if(child->is_adjacent(u_item->getParent())) {
                            // Parent is adjacent -> W list
                            child->W.insert(u_item);
                            u_item->X.insert(child);
                            break;
                        } else {
                            u_item = u_item->getParent();
                        }
                    }
                    if(u_item->getLevel() == child->getLevel()) {
                        // No adjacent ancestor -> add a neighbour
                        child->V.insert(u_item);
                        u_item->V.insert(child);
                    }
                }
            }

            // Find where to put W list items in child
            for(FNode* w_item : this->W) {
                // Find first ancestor of w_item that is adjacent to child
                // not needed, done in U list treatment, only check parent
                if(child->getLevel() < w_item->getLevel()) {
                    if(child->is_adjacent(w_item->getParent())) {
                        child->W.insert(w_item);
                        w_item->X.insert(child);
                    }
                } else if(child->getLevel() == w_item->getLevel()) {
                    // No adjacent ancestor -> add a neighbour
                    child->V.insert(w_item);
                    w_item->V.insert(child);
                }
            }
        }
        // Remove this from other lists
        for(FNode* w_item : this->W) {
            w_item->X.erase(this);
        }
        for(FNode* u_item : this->U) {
            u_item->U.erase(this);
        }
        // Clear leaf-only lists
        this->U.clear();
        this->W.clear();

        move_particles_to_children();
    }


    /**
     * \brief Splits a particle tuple into a postion and attributes
     *
     * \param t Tuple to split
     *
     * \tparam Tuple Tuple type, std::get<I>(t) must extract the I_th element
     */
    template<std::size_t... Is, std::size_t... Js, class Tuple>
    static auto get_position_and_attributes_impl(
        const Tuple& t,
        inria::index_sequence<Is...>, inria::index_sequence<Js...>)
    {
        return std::make_pair(position_t(std::get<Is>(t)...),
                              std::make_tuple(std::get<Js+Dim>(t)...));
    }

    /**
     * \brief Splits a particle tuple into a postion and attributes
     *
     * \param t Tuple to split
     *
     * \tparam Tuple Tuple type, std::get<I>(t) must extract the I_th element
     */
    template<class Tuple>
    static auto get_position_and_attributes(const Tuple& t)
    {
        return get_position_and_attributes_impl(
            t,
            inria::make_index_sequence<Dim>{} ,
            inria::make_index_sequence<std::tuple_size<Tuple>::value - Dim>{});
    }

    template<typename Node, typename Multipole>
    static Multipole *getMultipoleDataFromNode(Node *n) {
        return n == nullptr ? nullptr
            : &(n->getData()->getMultipoleData());
    }

    template<typename Node, typename SymbolicData>
    static SymbolicData *getSymbolicDataFromNode(Node *n) {
        return n == nullptr ? nullptr
            : &(n->getSymbolicData());
    }

    template<typename Node, typename LocalExpansion>
    static LocalExpansion *getLocalExpansionDataFromNode(Node *n) {
        return n == nullptr ? nullptr
            : &(n->getData()->getLocalExpansionData());
    }

    /**
     * \brief Reinsert particle in tree after spliting a container
     *
     * The #insert member function cannot be used to reinsert particles because
     * the value_type of the container may not be a particle.
     *
     * This overload is used when the value_type is a tuple.
     *
     * \param p Tuple representation of the particle to insert
     */
    template<class... Types>
    void reinsert(std::tuple<Types...> p) {
        reinsert_impl(get_position_and_attributes(p),
                      inria::make_index_sequence<sizeof...(Types)-Dim>{});
    }

    /**
     * \brief Implementation of particle reinsertion
     *
     * \param pos_attrs The position and attributes tuple of the particle.
     *
     * \tparam Pair
     * \parblock
     * A std::pair like object that holds two attributes
     *
     *   - `first` is the position of the particle
     *   - `second` is a tuple containing the particle attributes
     * \endparblock
     *
     */
    template<class Pair, std::size_t... Is>
    void reinsert_impl(const Pair& pos_attrs, inria::index_sequence<Is...>) {
        this->insert(pos_attrs.first, std::get<Is>(pos_attrs.second)...);
    }


    /**
     * \brief Reinsert particle in tree after spliting a container
     *
     * The #insert member function cannot be used to reinsert particles because
     * the value_type of the container may not be a particle.
     *
     * This overload is used when the value_type is a particle.
     *
     * \param p Particle to insert
     */
    void reinsert(const particle_t& p) {
        this->insert(p);
    }


    /**
     * \brief Reinsert particles once children have been added
     *
     * The particle container is deleted afterwards.
     */
    void move_particles_to_children() {
        for(auto&& p : *getParticleContainer()) {
            this->reinsert(p);
        }
        delete_particle_container();
    }


    /**
     * \brief Fuses the children nodes until this node becomes a leaf
     */
    void fuse() {
        if(is_leaf()) {
            return; // In a leaf, there's nothing to do
        }

        for(FNode* child : getChildren()) {
            child->fuse(); // Fuse children into leaves
        }

        create_particle_container();
        _is_leaf = true;

        // Remove children from U lists
        for(FNode* child_1 : getChildren()) {
            for(FNode* child_2 : getChildren()) {
                child_1->U.erase(child_2);
            }
        }
        // Use the children interaction lists to update this one
        for(FNode* child : getChildren()) {
            // Child U list items get into this U list
            for(FNode* u_item : child->U) {
                // Remove child from u_item U list
                u_item->U.erase(child);
                // Add this and u_item in each other's U lists
                this->U.insert(u_item);
                u_item->U.insert(this);
            }
            child->U.clear();

            // Child X list items get into the this U list
            for(FNode* x_item : child->X) {
                // Remove child from x_item W list
                x_item->W.erase(child);
                // Add this and x_item in each other's U lists
                this->U.insert(x_item);
                x_item->U.insert(this);
            }

            // Child W items get into this W list
            for(FNode* w_item : child->W) {
                // Remove child from w_item X list
                w_item->X.erase(child);
                // Add this and w_item in each other's W and X lists
                // when w_item is not adjacent to this
                if(! is_adjacent(w_item)) {
                    this->W.insert(w_item);
                    w_item->X.insert(this);
                }
            }

            // Child V list items get into this W list
            for(FNode* v_item : child->V) {
                // Remove child from the v_item V list
                v_item->V.erase(child);
                // Add this and v_item in each other's W and X lists
                // when v_item is not adjacent to this
                if(! is_adjacent(v_item)) {
                    this->W.insert(v_item);
                    v_item->X.insert(this);
                }
            }

            for(auto&& p : *(child->getParticleContainer())) {
                insert(p);
            }
        }
        // This belongs to its own U list
        this->U.insert(this);
        // Clear the children after everything is updated
        delete_children();

    }

    void terminal_print(std::ostream& os, bool print_particles) const {

        // Indentation
        if(! this->getParent()) {
            os << "└─ ";
        } else {
            std::string indent;
            const FNode* parent = this->getParent();
            while(parent) {
                if((parent->getIndex() & (child_count-1)) == child_count-1
                   || ! parent->getParent()) {
                    indent = "   " + indent;
                } else {
                    indent = "│  " + indent;
                }
                parent = parent->getParent();
            }
            if((this->getIndex() & (child_count-1)) == child_count-1) {
                os << indent << "└─ ";
            } else{
                os << indent << "├─ ";
            }
        }

        // Node data and recurrence
        auto index_width =
            std::setw(static_cast<int>(
                          std::log10(1 << (Dim * (this->getLevel()+1)))));

        if(! this->is_leaf()) { // Internal node
            os << "node " << index_width << this->getIndex() << ":\n";
            for(auto&& child : this->getChildren()) {
                child->terminal_print(os, print_particles);
            }
        } else { // Leaf
            os << "leaf "<< index_width << this->getIndex() <<": ";
            const auto& p_container = *(this->getParticleContainer());
            auto count_width =
                std::setw(static_cast<int>(
                              std::log10(this->getTree().leaf_max_particle_count())));
            os << count_width << p_container.size() << " particle";
            os << (p_container.size() > 1 ? "s" : "");

            if(print_particles && p_container.size() > 0) {
                os << " [";
                auto it = inria::make_ostream_joiner(os, ", ");
                std::transform(p_container.begin(), p_container.end(),
                               it, scalfmm::tuple_out);
                os << "]";
            }
            os << '\n';
        }
    }

    void obj_print(std::ostream& os) const {
        if(Dim != 3 && Dim != 2) {
            return;
        }
        // When not on a leaf, recurse to leaf
        if(! this->is_leaf()) {
            for(auto& child : getChildren()) {
                child->obj_print(os);
            }
        } else {
            // Print each of the box vertices
            const box_t& box = this->getBox();
            // note, child_count is the box's corner count
            for(std::size_t i = 0; i < child_count; ++i) {
                const position_t p = box.corner(i);
                os << "v";
                for(auto c : p) {
                    os << ' ' << c;
                }
                os << '\n';
            }

            // Print the box faces
            if(Dim == 3) {
                os << "f -8 -7 -5 -6\n";
                os << "f -8 -7 -3 -4\n";
                os << "f -6 -5 -1 -2\n";
                os << "f -4 -3 -1 -2\n";
                os << "f -7 -5 -1 -3\n";
                os << "f -8 -6 -2 -4\n";
                os << '\n';
            } else if(Dim == 2) {
                os << "f -4 -3 -1 -2\n";
            }
        }
    }

public:
    friend std::ostream& operator<<(std::ostream& os, const FNode& node) {
        if(os.iword(scalfmm::fmt::node_os_format_id()) == scalfmm::fmt::OBJ) {
            node.obj_print(os);
        } else {
            bool p_parts = os.iword(scalfmm::fmt::node_os_particle_id());
            node.terminal_print(os, p_parts);
        }
        return os;
    }

    friend morton_index_t morton_index(const FNode& n) {
        return n.getIndex();
    }

    friend level_t level(const FNode& n) {
        return n.getLevel();
    }
};


#endif
