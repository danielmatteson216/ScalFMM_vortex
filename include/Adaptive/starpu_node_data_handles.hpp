/**
 * \brief StarPU handle management adaptive algorithms
 * \file
 *
 * \author Quentin Khan
 */

#ifndef _SCALFMM_STARPU_NODE_DATA_HANDLES_HPP_
#define _SCALFMM_STARPU_NODE_DATA_HANDLES_HPP_

// @FUSE_STARPU

#include <starpu.h>

#include <utility>

/**
 * \brief Data handle manager for StarPU dependency management
 *
 * Holds data handles for the data parts of a node.
 *
 * \note The class is marked final because its destructor, which manages the
 * handles unregistration is not marked virtual.
 */
struct node_data_handles final {

    /// StarPU data handle for node symbolic data
    starpu_data_handle_t symbolic;
    /// StarPU data handle for node multipole expansion data
    starpu_data_handle_t multipole;
    /// StarPU data handle for node local expansion data
    starpu_data_handle_t local_exp;
    /// StarPU data handle for node particle data, used for leaf nodes
    starpu_data_handle_t particles;

    /// Flag to distinguish nodes that use the additionnal #particles handle
    bool registered_particles = false;

    /**
     * \brief Deleted copy constructor
     *
     * Copy is forbidden to avoid double unregistration on object destruction.
     */
    node_data_handles(const node_data_handles&) = delete;
    /**
     * \brief Deleted copy assignment operator
     *
     * Copy is forbidden to avoid double unregistration on object destruction.
     */
    node_data_handles& operator=(const node_data_handles&) = delete;

    /**
     * \brief Move constructor
     */
    node_data_handles(node_data_handles&& other) {
        *this = std::move(other);
    }
    /**
     * \brief Move assignment operator
     *
     * \param other Data handles to swap with
     *
     * \implementation Swaps the contents of this and other.
     */
    node_data_handles& operator=(node_data_handles&& other) {
        using std::swap;
        swap(this->symbolic, other.symbolic);
        swap(this->multipole, other.multipole);
        swap(this->local_exp, other.local_exp);
        swap(this->particles, other.particles);
        return *this;
    }

    /**
     * \brief Build and register the data handles for a node
     *
     * \param n Node to register data for
     *
     * \tparam node_t Node type
     */
    template<class node_t>
    node_data_handles(node_t& n) : registered_particles(false) {
        starpu_variable_data_register(
            &(this->symbolic),
            STARPU_MAIN_RAM,
            reinterpret_cast<uintptr_t>(&(n.getSymbolicData())),
            sizeof(n.getSymbolicData()));
        starpu_variable_data_register(
            &(this->multipole),
            STARPU_MAIN_RAM,
            reinterpret_cast<uintptr_t>(&(n.getData()->getMultipoleData())),
            sizeof(typename node_t::data_t::multipole_t));
        starpu_variable_data_register(
            &(this->local_exp),
            STARPU_MAIN_RAM,
            reinterpret_cast<uintptr_t>(&(n.getData()->getLocalExpansionData())),
            sizeof(typename node_t::data_t::local_expansion_t));
        if(n.is_leaf()) {
            starpu_variable_data_register(
                &(this->particles),
                STARPU_MAIN_RAM,
                reinterpret_cast<uintptr_t>(n.getParticleContainer()),
                sizeof(*(n.getParticleContainer())));
            this->registered_particles = true;
        }
    }

    /**
     * \brief Unregister all handles
     */
    ~node_data_handles() {
        starpu_data_unregister(this->symbolic);
        starpu_data_unregister(this->multipole);
        starpu_data_unregister(this->local_exp);
        if(this->registered_particles) {
            starpu_data_unregister(this->particles);
            this->registered_particles = false;
        }
    }
};



#endif /* _SCALFMM_STARPU_NODE_DATA_HANDLES_HPP_ */
