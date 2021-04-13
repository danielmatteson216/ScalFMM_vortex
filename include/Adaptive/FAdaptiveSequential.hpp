#ifndef SCALFMM_SEQUENTIAL_ALGO_HPP_
#define SCALFMM_SEQUENTIAL_ALGO_HPP_

#include <algorithm>
#include <cmath> // Used to round box differences
#include <functional>
#include <map>
#include <vector>


#include "Core/FCoreCommon.hpp"
#include "Containers/FTreeCoordinate.hpp"
#include "Utils/FAlgorithmTimers.hpp"

#include "Kernels/FKernelConcepts.hpp"

#include "kernel_utilities.hpp"

template<
    class _Tree, class _Kernel,
    class = inria::require<
        scalfmm::meta::adaptive_compatible<_Tree,_Kernel>
        >
    >
class FAdaptiveSequential : public FAlgorithmInterface, public FAlgorithmTimers {
public:
    using tree_t = _Tree;
    using kernel_t = _Kernel;
    using FReal = typename tree_t::FReal;

    using node_t = typename tree_t::node_t;

    using multipole_t = typename node_t::data_t::multipole_t;
    using local_expansion_t = typename node_t::data_t::local_expansion_t;
    using symbolic_data_t = typename node_t::symbolic_data_t;

private:
    tree_t&   _tree;
    kernel_t& _kernel;

public:

    FAdaptiveSequential(tree_t* tree, kernel_t* kernel) :
        _tree(*tree),
        _kernel(*kernel) {
    }

    std::string name() const override {
        return "Sequential adaptive algorithm";
    }

    std::string description() const override {
        return "";
    }

    using FAlgorithmInterface::execute;

    /** \brief Run specific steps of the algorithm
     *
     * \param operations Specifies the algorithm operations to run, see
     * FFmmOperations.
     */
    void execute(const unsigned int operations) override {
        this->run(operations);
    }

    void run(int operations = FFmmNearAndFarFields) {

        scalfmm::setup_kernel(this->_kernel, this->_tree);

        if(operations & FFmmP2M) {
            // 1. source to up, P2M
            Timers["P2M"].tic();
            this->source_to_up();;
            Timers["P2M"].tac();
            std::cout << "    P2M: " << Timers["P2M"].elapsed() << std::endl;
        }

        if(operations & FFmmM2M) {
            // 2. up to up, M2M
            Timers["M2M"].tic();
            this->up_to_up();
            Timers["M2M"].tac();
            std::cout << "    M2M: " << Timers["M2M"].elapsed() << std::endl;
        }

        if(operations & FFmmM2L) {
            // 3a V-list, M2L
            Timers["M2L"].tic();
            this->v_list_step();
            Timers["M2L"].tac();
            std::cout << "    M2L: " << Timers["M2L"].elapsed() << std::endl;
        }

        if(operations & FFmmP2L) {
            // 3b X-list, P2L
            Timers["P2L"].tic();
            this->x_list_step();
            Timers["P2L"].tac();
            std::cout << "    P2L: " << Timers["P2L"].elapsed() << std::endl;
        }

        if(operations & FFmmL2L) {
            // 4. down to down, L2L
            Timers["L2L"].tic();
            this->down_to_down();
            Timers["L2L"].tac();
            std::cout << "    L2L: " << Timers["L2L"].elapsed() << std::endl;
        }

        if(operations & FFmmM2P) {
        // 5a W-list, M2P
            Timers["M2P"].tic();
            this->w_list_step();
            Timers["M2P"].tac();
            std::cout << "    M2P: " << Timers["M2P"].elapsed() << std::endl;
        }

        if(operations & FFmmL2P) {
            // 5b down to target, L2P
            Timers["L2P"].tic();
            this->down_to_target();
            Timers["L2P"].tac();
            std::cout << "    L2P: " << Timers["L2P"].elapsed() << std::endl;
        }

        if(operations & FFmmP2P) {
            // A. U-list, P2P
            Timers["P2P"].tic();
            this->u_list_step();
            Timers["P2P"].tac();
            std::cout << "    P2P: " << Timers["P2P"].elapsed() << std::endl;
        }
    }

    // P2M
    void source_to_up() {
        for(node_t* leaf : _tree.leaves()) {
            if(leaf->getParticleContainer()->size() != 0) {
                _kernel.P2M(&(leaf->getData()->getMultipoleData()),
                            &(leaf->getSymbolicData()),
                            leaf->getParticleContainer());
            }
        }
    }

    // M2M
    void up_to_up() {
        for(node_t& node : _tree.post_order_walk()) {
            if(! node.is_leaf()) {
                // Fill array of child data
                std::array<multipole_t*, node_t::child_count> child_multipoles {};
                std::transform(std::begin(node.getChildren()), std::end(node.getChildren()),
                               child_multipoles.begin(),
                               [](node_t* n) {
                                   return n == nullptr ? nullptr
                                       : &(n->getData()->getMultipoleData());
                               });

                std::array<symbolic_data_t*, node_t::child_count> child_symbolics {};
                std::transform(std::begin(node.getChildren()), std::end(node.getChildren()),
                               child_symbolics.begin(),
                               [](node_t* n) {
                                   return n == nullptr ? nullptr
                                       : &(n->getSymbolicData());
                               });

                // Call kernel module
                _kernel.M2M(&(node.getData()->getMultipoleData()),
                            &(node.getSymbolicData()),
                            child_multipoles.data(),
                            child_symbolics.data()
                    );
            }
        }
    }

    // M2L
    void v_list_step() {
        std::vector<multipole_t*> v_item_multipoles;
        std::vector<symbolic_data_t*> v_item_symbolics;
        std::vector<int> v_item_indices;

        for(node_t& node : _tree.in_order_walk()) {
            if(node.is_leaf() && node.getParticleContainer()->size() == 0) {
                continue;
            }

            v_item_multipoles.clear();
            v_item_symbolics.clear();
            v_item_indices.clear();

            // Needed to compute offset between boxes
            for(node_t* v_item : node.V) {
                if(v_item->is_leaf()
                   && v_item->getParticleContainer()->size() == 0) {
                    continue;
                }
                v_item_indices.push_back(compute_box_offset_index(&node, v_item, 3));
                v_item_multipoles.push_back(&(v_item->getData()->getMultipoleData()));
                v_item_symbolics.push_back(&(v_item->getSymbolicData()));
            }

            // Call kernel M2L operator
            _kernel.M2L(&(node.getData()->getLocalExpansionData()),
                        &(node.getSymbolicData()),
                        v_item_multipoles.data(),
                        v_item_symbolics.data(),
                        v_item_indices.data(),
                        static_cast<int>(v_item_multipoles.size())
                );
        }
    }

    // P2L
    void x_list_step() {
        /* NOTE: the X list and W list are complementary: if A is in X(B) then B
         * is in W(A).
         *
         * We loop over the leaves first to detect the empty ones early on.
         */
        for(node_t* leaf : _tree.leaves()) {
            if(leaf->getParticleContainer()->size() > 0) {
                for(node_t* w_item : leaf->W) {
                    if(w_item->is_leaf() && w_item->getParticleContainer()->size() == 0) {
                        continue;
                    }
                    _kernel.P2L(&(w_item->getData()->getLocalExpansionData()),
                                &(w_item->getSymbolicData()),
                                leaf->getParticleContainer(),
                                &(leaf->getSymbolicData())
                        );
                }
            }
        }
    }

    // L2L
    void down_to_down() {
        for(node_t& node : _tree.pre_order_walk()) {
            if(! node.is_leaf()) {

                std::array<local_expansion_t*, node_t::child_count> child_local_expansions {};
                std::transform(std::begin(node.getChildren()), std::end(node.getChildren()),
                               child_local_expansions.begin(),
                               [](node_t* n) {
                                   return n == nullptr ? nullptr
                                       : &(n->getData()->getLocalExpansionData());
                               });

                std::array<symbolic_data_t*, node_t::child_count> child_symbolics {};
                std::transform(std::begin(node.getChildren()), std::end(node.getChildren()),
                               child_symbolics.begin(),
                               [](node_t* n) {
                                   return n == nullptr ? nullptr
                                       : &(n->getSymbolicData());
                               });

                _kernel.L2L(&(node.getData()->getLocalExpansionData()),
                            &(node.getSymbolicData()),
                            child_local_expansions.data(),
                            child_symbolics.data()
                    );
            }
        }
    }

    // M2P
    void w_list_step() {
        for(node_t* leaf : _tree.leaves()) {
            if(leaf->getParticleContainer()->size() > 0) {
                for(node_t* w_item : leaf->W) {
                    if(w_item->is_leaf() && w_item->getParticleContainer()->size() == 0) {
                        continue;
                    }
                    _kernel.M2P(&(w_item->getData()->getMultipoleData()),
                                &(w_item->getSymbolicData()),
                                leaf->getParticleContainer(),
                                &(leaf->getSymbolicData())
                        );
                }
            }
        }
    }

    // L2P
    void down_to_target() {
        for(node_t* leaf : _tree.leaves()) {
            if(leaf->getParticleContainer()->size() != 0) {
                _kernel.L2P(&(leaf->getData()->getLocalExpansionData()),
                            &(leaf->getSymbolicData()),
                            leaf->getParticleContainer());
            }
        }
    }


    /** \brief Direct computation step (P2P)
     *
     * For each tree leaf, a direct computation is done with the adjacent
     * leaves.
     *
     * The symmetric computation kernels expect to receive the list of adjacent
     * leaves, stored in a leaf's U list, sorted by offset index.
     *
     * The offset index is of a node compared to another follows the
     * numerotation of the adjacent nodes, on each axis. In 3D, this woud give,
     * with L the reference leaf:
     *
     *          x          x            x
     *        0 1 2  |   9 10 11  |  18 19 20
     *     y  3 4 5  |  12  L 14  |  21 22 23
     *        6 7 8  |  15 16 17  |  24 25 26
     *         z=0        z=1          z=2
     *
     * When two node that are to be compared are not on the same level, the
     * lowest one's ancestor that is at the highest one's level is used to
     * compute the index.
     *
     * To statisfy this condition, the leaves have to be sorted before being
     * given to the kernel P2P method.
     *
     */
    void u_list_step() {
        using container_t = typename node_t::particle_container_t;

        // Vectors to pass to kernel P2P
        std::vector<container_t*> u_item_source_particle_containers;
        std::vector<int> u_item_indices;

        for(node_t* leaf : _tree.leaves()) {

            container_t* const leaf_source_particle_container =
                leaf->getParticleContainer();
            container_t* const leaf_target_particle_container =
                leaf->getParticleContainer();

            // Skip empty leaves
            if( leaf_source_particle_container->size() == 0
                && leaf_target_particle_container->size() == 0) {
                continue;
            }

            u_item_source_particle_containers.clear();
            u_item_indices.clear();

            for(node_t* u_item : leaf->U) {
                // The kernels do not consider leaf to be adjacent to
                // itself. Skip empty u_items
                if(u_item == leaf
                   || u_item->getParticleContainer()->size() == 0)
                {
                    continue;
                }
                u_item_source_particle_containers.push_back(u_item->getParticleContainer());
                u_item_indices.push_back(compute_box_offset_index(leaf, u_item, 1));
            }

            // Call P2P on vectors data
            _kernel.P2P(
                FTreeCoordinate(MortonIndex(leaf->getIndex())),
                leaf_target_particle_container,
                leaf_source_particle_container,
                u_item_source_particle_containers.data(),
                u_item_indices.data(),
                static_cast<int>(u_item_source_particle_containers.size())
                );
        }
    }


    /** Compute resulting index in pre-computation array.
     *
     * Kernels precompute some operators according to node offsets. This method
     * returns the `other_node` index in the precomputation array.
     *
     * The index has the following form, with d the space dimension:
     *
     * x+n * 7^d + y+n * 7^(d-1) + ... + z+n
     *
     * Because an `other_node` is at most `n` boxes away in every direction from a
     * `node`, in order to get a positive index, we add `n` to the offset in every
     * direction.
     *
     * \param node The target node
     * \param other_node The interacting node
     * \param n The longest possible distance between node and other_node in
     * terms of boxes on an axis
     *
     * \warning If the nodes are not at the same level, the lowest one's
     * ancestor at the highest one's level is used.
     */
    int compute_box_offset_index(node_t* node, node_t* other_node, const std::size_t n) {
        while(node->getLevel() > other_node->getLevel()) {
            node = node->getParent();
        }
        while(node->getLevel() < other_node->getLevel()) {
            other_node = other_node->getParent();
        }

        typename node_t::FReal box_width = node->getBox().width(0);

        auto center_offset = other_node->getBox().center() - node->getBox().center();
        std::size_t other_node_index = 0;
        for(std::size_t i = 0; i < node_t::Dim; ++i) {
            other_node_index = other_node_index * (2 * n + 1)
                + static_cast<unsigned long>(std::lround(center_offset[i] / box_width))
                + n;
        }

        return static_cast<int>(other_node_index);
    }


};

#endif
