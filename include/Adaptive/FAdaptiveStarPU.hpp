#ifndef SCALFMM_STARPU_ALGO_HPP_
#define SCALFMM_STARPU_ALGO_HPP_

//@FUSE_STARPU

#include <algorithm>
#include <cmath> // Used to round box differences
#include <functional>
#include <map>
#include <memory>
#include <vector>
#include <unordered_map>

#include <starpu.h>

#include "Core/FCoreCommon.hpp"
#include "Containers/FTreeCoordinate.hpp"
#include "Utils/FAlgorithmTimers.hpp"

#include "Kernels/FKernelConcepts.hpp"

#include "kernel_utilities.hpp"
#include "starpu_node_data_handles.hpp"

template<class _Tree, class _Kernel,
         class = inria::require<
             scalfmm::meta::adaptive_compatible<_Tree,_Kernel>
             >
         >
class FAdaptiveStarPU : public FAlgorithmInterface, public FAlgorithmTimers {
public:
    using tree_t = _Tree;
    using kernel_t = _Kernel;
    using FReal = typename tree_t::FReal;

    using node_t = typename tree_t::node_t;

private:
    using symbolic_data_t = typename tree_t::node_t::symbolic_data_t;
    using multipole_t = typename tree_t::node_t::data_t::multipole_t;
    using local_expansion_t = typename tree_t::node_t::data_t::local_expansion_t;
    using container_t = typename tree_t::node_t::particle_container_t;

    tree_t&   _tree;

    kernel_t& _ref_kernel;
    /// Vector of kernels, one per thread
    std::vector<std::unique_ptr<kernel_t>> _kernels;

    starpu_codelet P2M_cl;
    starpu_codelet M2M_cl;
    starpu_codelet M2L_cl;
    starpu_codelet L2L_cl;
    starpu_codelet L2P_cl;
    starpu_codelet P2P_cl;
    starpu_codelet M2P_cl;
    starpu_codelet P2L_cl;

    std::unordered_map<node_t*, node_data_handles> _data_handles;

public:

    FAdaptiveStarPU(tree_t* tree, kernel_t* kernel) :
        _tree(*tree),
        _ref_kernel(*kernel),
        _kernels()
    {}

    std::string name() const override {
        return "StarPU adaptive algorithm";
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

        scalfmm::setup_kernel(this->_ref_kernel, this->_tree);
        this->setup_starpu();

        auto run = [this](const char* name, void (FAdaptiveStarPU::*op)()) {
            Timers[name].tic();
            (this->*op)();;
            Timers[name].tac();
            std::cout << "  " << name << ": " << Timers[name].elapsed() << std::endl;
        };

        if(operations & FFmmP2P) {
            run("P2P", &FAdaptiveStarPU::u_list_step);
        }

        if(operations & FFmmP2M) {
            run("P2M", &FAdaptiveStarPU::source_to_up);
        }

        if(operations & FFmmP2L) {
            run("P2L", &FAdaptiveStarPU::x_list_step);
        }

        if(operations & FFmmM2M) {
            run("M2M", &FAdaptiveStarPU::up_to_up);
        }

        if(operations & FFmmM2L) {
            run("M2L", &FAdaptiveStarPU::v_list_step);
        }

        if(operations & FFmmM2P) {
            run("M2P", &FAdaptiveStarPU::w_list_step);
        }

        if(operations & FFmmL2L) {
            run("L2L", &FAdaptiveStarPU::down_to_down);
        }

        if(operations & FFmmL2P) {
            run("L2P", &FAdaptiveStarPU::down_to_target);
        }

        this->cleanup_starpu();
        scalfmm::cleanup_kernel(this->_ref_kernel, this->_tree);
    }


private:

    template<class K>
    auto fuse_kernel_results(K* ref, K* other) -> decltype(ref->fuse_results(*other)) {
        return ref->fuse_results(*other);
    }

    void fuse_kernel_results(...) {}

    starpu_codelet setup_worker_kernel_cl;

    static void setup_worker_kernel(void**, void* cl_arg) {
        FAdaptiveStarPU* algo = (FAdaptiveStarPU*) cl_arg;
        auto& ref_kernel = algo->_ref_kernel;
        auto& kernels = algo->_kernels;

        int id = starpu_worker_get_id();
        std::cout << "initialising worker " << id << '\n';
        kernels[id].reset(new kernel_t(ref_kernel));
        scalfmm::setup_kernel(*kernels[id].get(), algo->_tree);
    }

    void setup_starpu() {
        if(0 != starpu_init(NULL)) {
            std::cerr << "StarPU could not be initialized.\n";
            std::exit(EXIT_FAILURE);
        }

        this->_kernels.resize(starpu_worker_get_count());
        std::cerr << "StarPU worker count: " << starpu_worker_get_count() << "\n";

        starpu_codelet_init(&this->setup_worker_kernel_cl);
        setup_worker_kernel_cl.cpu_funcs[0] = setup_worker_kernel;
        for(unsigned int i = 0; i < starpu_worker_get_count(); ++i) {
            starpu_task* task = starpu_task_create();
            task->cl = &this->setup_worker_kernel_cl;
            task->cl_arg = this;
            task->execute_on_a_specific_worker = true;
            task->workerid = i;
            task->synchronous = true;
            this->submit_starpu_task(task);
        }

        this->init_P2M_codelet();
        this->init_M2M_codelet();
        this->init_M2L_codelet();
        this->init_L2L_codelet();
        this->init_L2P_codelet();
        this->init_P2P_codelet();
        this->init_M2P_codelet();
        this->init_P2L_codelet();

        for(auto& node : _tree.pre_order_walk()) {
            this->_data_handles.emplace(&node, node);
        }
    }

    starpu_codelet cleanup_worker_kernel_cl;

    static void cleanup_worker_kernel(void**, void* cl_arg) {
        FAdaptiveStarPU* algo = (FAdaptiveStarPU*) cl_arg;
        auto& ref_kernel = algo->_ref_kernel;
        auto& kernels = algo->_kernels;

        int id = starpu_worker_get_id();
        std::cout << "cleaning up worker " << id << '\n';
        algo->fuse_kernel_results(&ref_kernel, kernels[id].get());
        scalfmm::cleanup_kernel(*kernels[id].get(), algo->_tree);
    }


    void cleanup_starpu() {
        starpu_task_wait_for_all();
        this->_data_handles.clear();

        starpu_codelet_init(&this->cleanup_worker_kernel_cl);
        cleanup_worker_kernel_cl.cpu_funcs[0] = cleanup_worker_kernel;
        for(unsigned int i = 0; i < starpu_worker_get_count(); ++i) {
            starpu_task* task = starpu_task_create();
            task->cl = &this->cleanup_worker_kernel_cl;
            task->cl_arg = this;
            task->execute_on_a_specific_worker = true;
            task->workerid = i;
            task->synchronous = true;
            this->submit_starpu_task(task);
        }

        starpu_shutdown();
    }

    void submit_starpu_task(starpu_task* task) {
        const auto res = starpu_task_submit(task);
        if(0 != res) {
            std::cerr << task->cl->name << " was not submitted...\n";
            if(res == -ENODEV) {
                std::cerr << "No device found to execute task.\n";
            }
        }
    }

    static kernel_t* get_starpu_worker_kernel(FAdaptiveStarPU* algo) {
        int worker_id = starpu_worker_get_id();
        return algo->_kernels[worker_id].get();
    }


    /**
     * \brief Allocate dynamic handles and mode buffers if needed
     *
     * \warning `task->nbuffers` must have been set
     *
     * \return A tuple holding the handles and modes buffers. In no allocation
     * took place, those are equal to `task->handles` and `task->modes`.
     */
    std::tuple<starpu_data_handle_t*, starpu_data_access_mode*>
    allocate_dynamic_handles_if_needed(starpu_task*& task) {
        starpu_data_handle_t* handles = task->handles;
        starpu_data_access_mode* modes = task->modes;

        if(task->nbuffers > STARPU_NMAXBUFS) {
            // Allocate dynamic StaPU dynamic handles
            // TODO: check that those are freed by StarPU
            task->dyn_handles = (starpu_data_handle_t*)
                std::malloc(task->nbuffers * sizeof(starpu_data_handle_t));
            task->dyn_modes = (starpu_data_access_mode*)
                std::malloc(task->nbuffers * sizeof(starpu_data_access_mode));
            handles = task->dyn_handles;
            modes = task->dyn_modes;
        }
        return std::make_tuple(handles, modes);
    }


    // P2M
    void source_to_up() {
        for(node_t* leaf : this->_tree.leaves()) {
            if(leaf->getParticleContainer()->size() != 0) {
                node_data_handles& leaf_handles = this->_data_handles.at(leaf);

                starpu_task* task = starpu_task_create();
                task->cl = &(this->P2M_cl);
                task->handles[0] = leaf_handles.multipole;
                task->handles[1] = leaf_handles.symbolic;
                task->handles[2] = leaf_handles.particles;
                task->cl_arg = this;
                task->cl_arg_size = sizeof(this);

                this->submit_starpu_task(task);
            }
        }
    }

    void init_P2M_codelet() {
        starpu_codelet_init(&this->P2M_cl);
        this->P2M_cl.name = "P2M";
        this->P2M_cl.cpu_funcs[0] = P2M_cpu;
        this->P2M_cl.cpu_funcs_name[0] = {const_cast<char*>("P2M task")};
        this->P2M_cl.nbuffers = 3; // {0: leaf multipole, 1: leaf symb, 2: leaf particles}
        this->P2M_cl.modes[0] = STARPU_RW;
        this->P2M_cl.modes[1] = STARPU_R;
        this->P2M_cl.modes[2] = STARPU_R;
    }

    static void P2M_cpu(void** buffers, void* cl_arg) {
        multipole_t* leaf_multipole =
            (multipole_t*) STARPU_VARIABLE_GET_PTR(buffers[0]);

        symbolic_data_t* leaf_symb_data =
            (symbolic_data_t*) STARPU_VARIABLE_GET_PTR(buffers[1]);

        container_t* particle_container =
            (container_t*) STARPU_VARIABLE_GET_PTR(buffers[2]);

        auto algo = (FAdaptiveStarPU<tree_t, kernel_t>*) cl_arg;
        kernel_t* ker = get_starpu_worker_kernel(algo);
        ker->P2M(leaf_multipole, leaf_symb_data, particle_container);
    }


    // M2M
    void up_to_up() {
        for(node_t& node : _tree.post_order_walk()) {
            if(node.is_leaf()) {
                continue;
            }

            starpu_task* task = starpu_task_create();
            task->cl = &(this->M2M_cl);
            task->cl_arg = this;
            task->cl_arg_size = sizeof(this);

            auto children = node.getChildren();
            // Setup handle array to be used and buffer count
            task->nbuffers = 2;
            // Get total handle count: mutipole + symbolic per child
            for(auto child : children) {
                task->nbuffers += (child != nullptr) * 2;
            }
            starpu_data_handle_t* handles = nullptr;
            starpu_data_access_mode* modes = nullptr;
            std::tie(handles, modes) = this->allocate_dynamic_handles_if_needed(task);

            // Set parent multipole handle
            handles[0] = this->_data_handles.at(&node).symbolic;
            modes[0] = STARPU_R;
            handles[1] = this->_data_handles.at(&node).multipole;
            modes[1] = STARPU_RW;

            // Set children multipole handles
            for(std::size_t i = 0, j = 2; i < node_t::child_count; ++i) {
                if(nullptr != children[i]) {
                    handles[j] = this->_data_handles.at(children[i]).symbolic;
                    modes[j] = STARPU_R;
                    handles[j+1] = this->_data_handles.at(children[i]).multipole;
                    modes[j+1] = STARPU_R;
                    j += 2;
                }
            }

            this->submit_starpu_task(task);
        }
    }

    void init_M2M_codelet() {
        starpu_codelet_init(&this->M2M_cl);
        this->M2M_cl.name = "M2M";
        this->M2M_cl.cpu_funcs[0] = M2M_cpu;
        this->M2M_cl.cpu_funcs_name[0] = {const_cast<char*>("M2M task")};
        this->M2M_cl.nbuffers = STARPU_VARIABLE_NBUFFERS;
    }

    static void M2M_cpu(void** buffers, void* cl_arg) {
        const starpu_task* const current_task =  starpu_task_get_current();
        std::size_t buffer_count = current_task->nbuffers;

        const symbolic_data_t* node_symbolic =
            (symbolic_data_t*) STARPU_VARIABLE_GET_PTR(buffers[0]);
        multipole_t* node_multipole =
            (multipole_t*) STARPU_VARIABLE_GET_PTR(buffers[1]);

        multipole_t* child_multipoles[node_t::child_count] = {};
        symbolic_data_t* child_symbolics[node_t::child_count] = {};

        // Children buffer indices start at 2
        for(std::size_t i = 0, j = 2, k = 3;
            k < buffer_count;
            ++i, j += 2, k += 2)
        {
            child_symbolics[i] =
                (symbolic_data_t*) STARPU_VARIABLE_GET_PTR(buffers[j]);
            child_multipoles[i] =
                (multipole_t*) STARPU_VARIABLE_GET_PTR(buffers[k]);
        }

        auto algo = (FAdaptiveStarPU<tree_t, kernel_t>*) cl_arg;
        kernel_t* ker = get_starpu_worker_kernel(algo);
        ker->M2M(node_multipole, node_symbolic,
                 child_multipoles, child_symbolics);
    }


    // M2L
    void v_list_step() {
        for(node_t& node : _tree.in_order_walk()) {
            // Empty leaves have an empty multipole, we skip them
            if(node.is_leaf() && node.getParticleContainer()->size() == 0) {
                continue;
            }

            starpu_task* task = starpu_task_create();
            task->cl = &(this->M2L_cl);
            task->cl_arg = this;
            task->cl_arg_size = sizeof(this);

            // Returns true if the V item must be used
            auto pick_item = [](const node_t* v) {
                return ! v->is_leaf() || v->getParticleContainer()->size() != 0;
            };
            auto v_item_count = std::count_if(node.V.begin(), node.V.end(), pick_item);

            // V-items symbolic data and multipole, plus the local expansion and
            // its symbolic data.
            task->nbuffers = 2 + (2 * static_cast<int>(v_item_count));

            starpu_data_handle_t* handles = nullptr;
            starpu_data_access_mode* modes = nullptr;
            std::tie(handles, modes) = this->allocate_dynamic_handles_if_needed(task);

            handles[0] = this->_data_handles.at(&node).symbolic;
            handles[1] = this->_data_handles.at(&node).local_exp;
            modes[0] = STARPU_R;
            modes[1] = STARPU_RW;

            std::size_t j = 2;
            for(node_t* v_item : node.V) {
                if(! pick_item(v_item)) {
                    continue;
                }
                handles[j+0] = this->_data_handles.at(v_item).symbolic;
                handles[j+1] = this->_data_handles.at(v_item).multipole;
                modes[j+0] = STARPU_R;
                modes[j+1] = STARPU_RW;
                j += 2;
            }

            this->submit_starpu_task(task);
        }
    }

    void init_M2L_codelet() {
        starpu_codelet_init(&this->M2L_cl);
        this->M2L_cl.name = "M2L";
        this->M2L_cl.cpu_funcs[0] = M2L_cpu;
        this->M2L_cl.cpu_funcs_name[0] = {const_cast<char*>("M2L task")};
        this->M2L_cl.nbuffers = STARPU_VARIABLE_NBUFFERS;
    }


    // buffer order: symbolic local_exp then symbolic multipole ....
    static void M2L_cpu(void** buffers, void* cl_arg) {
        const starpu_task* const current_task =  starpu_task_get_current();
        std::size_t buffer_count = current_task->nbuffers;

        const symbolic_data_t* node_symbolic =
            (symbolic_data_t*) STARPU_VARIABLE_GET_PTR(buffers[0]);
        local_expansion_t* node_local_exp =
            (local_expansion_t*) STARPU_VARIABLE_GET_PTR(buffers[1]);


        std::vector<const multipole_t*> multipoles;
        std::vector<const symbolic_data_t*> symbolics;
        std::vector<int> offset_indices;
        // two buffers per V item, first is for the local exp
        const std::size_t vec_init_size = (buffer_count-2) / 2;
        multipoles.reserve(vec_init_size);
        symbolics.reserve(vec_init_size);
        offset_indices.reserve(vec_init_size);

        for(std::size_t i = 2; i < buffer_count; i += 2) {
            const symbolic_data_t* v_item_symb = (symbolic_data_t*) STARPU_VARIABLE_GET_PTR(buffers[i]);
            symbolics.emplace_back(v_item_symb);
            multipoles.emplace_back((multipole_t*) STARPU_VARIABLE_GET_PTR(buffers[i+1]));
            offset_indices.emplace_back(compute_box_offset_index(*node_symbolic, *v_item_symb, 3));
        }

        auto algo = (FAdaptiveStarPU<tree_t, kernel_t>*) cl_arg;
        kernel_t* ker = get_starpu_worker_kernel(algo);
        ker->M2L(node_local_exp, node_symbolic,
                 multipoles.data(), symbolics.data(),
                 offset_indices.data(), static_cast<int>(multipoles.size()));
    }


    // P2L
    void x_list_step() {
        /* NOTE: the X list and W list are complementary: if A is in X(B) then B
         * is in W(A).
         *
         * We loop over the leaves first to detect the empty ones early on.
         */
        for(node_t* leaf : _tree.leaves()) {
            if(leaf->getParticleContainer()->size() == 0) {
                continue;
            }
            for(node_t* w_item : leaf->W) {
                if(w_item->is_leaf() && w_item->getParticleContainer()->size() == 0) {
                    continue;
                }

                starpu_task* task = starpu_task_create();
                task->cl = &(this->P2L_cl);
                task->handles[0] = this->_data_handles.at(w_item).symbolic;
                task->handles[1] = this->_data_handles.at(w_item).local_exp;
                task->handles[2] = this->_data_handles.at(leaf).symbolic;
                task->handles[3] = this->_data_handles.at(leaf).particles;
                task->cl_arg = this;
                task->cl_arg_size = sizeof(this);

                this->submit_starpu_task(task);
            }

        }
    }

    void init_P2L_codelet() {
        starpu_codelet_init(&this->P2L_cl);
        this->P2L_cl.name = "P2L";
        this->P2L_cl.cpu_funcs[0] = P2L_cpu;
        this->P2L_cl.cpu_funcs_name[0] = {const_cast<char*>("P2L task")};
        this->P2L_cl.nbuffers = 4; // {node symb, node local exp, leaf symb, leaf particles}
        this->P2L_cl.modes[0] = STARPU_R;
        this->P2L_cl.modes[1] = STARPU_RW;
        this->P2L_cl.modes[2] = STARPU_R;
        this->P2L_cl.modes[3] = STARPU_R;
    }

    static void P2L_cpu(void** buffers, void* cl_arg) {
        symbolic_data_t* node_symb_data =
            (symbolic_data_t*) STARPU_VARIABLE_GET_PTR(buffers[0]);

        local_expansion_t* node_local_exp =
            (local_expansion_t*) STARPU_VARIABLE_GET_PTR(buffers[1]);

        symbolic_data_t* leaf_symb_data =
            (symbolic_data_t*) STARPU_VARIABLE_GET_PTR(buffers[2]);

        container_t* leaf_particles =
            (container_t*) STARPU_VARIABLE_GET_PTR(buffers[3]);

        auto algo = (FAdaptiveStarPU<tree_t, kernel_t>*) cl_arg;
        kernel_t* ker = get_starpu_worker_kernel(algo);
        ker->P2L(node_local_exp, node_symb_data,
                 leaf_particles, leaf_symb_data);
    }


    // L2L
    void down_to_down() {
        for(node_t& node : _tree.pre_order_walk()) {
            if(node.is_leaf()) {
                continue;
            }

            starpu_task* task = starpu_task_create();
            task->cl = &(this->L2L_cl);
            task->cl_arg = this;
            task->cl_arg_size = sizeof(this);

            auto children = node.getChildren();
            // Setup handle array to be used and buffer count
            task->nbuffers = 2;
            // Get total handle count: mutipole + symbolic per child
            for(auto child : children) {
                task->nbuffers += (child != nullptr) * 2;
            }
            starpu_data_handle_t* handles = nullptr;
            starpu_data_access_mode* modes = nullptr;
            std::tie(handles, modes) = this->allocate_dynamic_handles_if_needed(task);

            // Set parent multipole handle
            handles[0] = this->_data_handles.at(&node).symbolic;
            modes[0] = STARPU_R;
            handles[1] = this->_data_handles.at(&node).local_exp;
            modes[1] = STARPU_R;

            // Set children multipole handles
            for(std::size_t i = 0, j = 2; i < node_t::child_count; ++i) {
                if(nullptr != children[i]) {
                    handles[j] = this->_data_handles.at(children[i]).symbolic;
                    modes[j] = STARPU_R;
                    handles[j+1] = this->_data_handles.at(children[i]).local_exp;
                    modes[j+1] = STARPU_RW;
                    j += 2;
                }
            }

            // Submit task
            this->submit_starpu_task(task);
        }
    }

    void init_L2L_codelet() {
        starpu_codelet_init(&this->L2L_cl);
        this->L2L_cl.name = "L2L";
        this->L2L_cl.cpu_funcs[0] = L2L_cpu;
        this->L2L_cl.cpu_funcs_name[0] = {const_cast<char*>("L2L task")};
        this->L2L_cl.nbuffers = STARPU_VARIABLE_NBUFFERS;
    }

    static void L2L_cpu(void** buffers, void* cl_arg) {
        const starpu_task* const current_task =  starpu_task_get_current();
        std::size_t buffer_count = current_task->nbuffers;

        const symbolic_data_t* parent_symbolic =
            (symbolic_data_t*) STARPU_VARIABLE_GET_PTR(buffers[0]);
        local_expansion_t* parent_local_exp =
            (local_expansion_t*) STARPU_VARIABLE_GET_PTR(buffers[1]);

        const symbolic_data_t* child_symbolics[node_t::child_count] = {};
        local_expansion_t* child_local_exps[node_t::child_count] = {};

        // Children buffer indices start at 2
        for(std::size_t i = 0, j = 2, k = 3;
            k < buffer_count;
            ++i, j += 2, k += 2)
        {
            child_symbolics[i] =
                (symbolic_data_t*) STARPU_VARIABLE_GET_PTR(buffers[j]);
            child_local_exps[i] =
                (local_expansion_t*) STARPU_VARIABLE_GET_PTR(buffers[k]);
        }

        auto algo = (FAdaptiveStarPU<tree_t, kernel_t>*) cl_arg;
        kernel_t* ker = get_starpu_worker_kernel(algo);
        ker->L2L(parent_local_exp, parent_symbolic,
                 child_local_exps, child_symbolics);
    }


    // M2P
    void w_list_step() {
        for(node_t* leaf : _tree.leaves()) {
            if(leaf->getParticleContainer()->size() == 0) {
                continue;
            }
            for(node_t* w_item : leaf->W) {
                if(w_item->is_leaf() && w_item->getParticleContainer()->size() == 0) {
                    continue;
                }

                starpu_task* task = starpu_task_create();
                task->cl = &(this->M2P_cl);
                task->handles[0] = this->_data_handles.at(w_item).symbolic;
                task->handles[1] = this->_data_handles.at(w_item).multipole;
                task->handles[2] = this->_data_handles.at(leaf).symbolic;
                task->handles[3] = this->_data_handles.at(leaf).particles;
                task->cl_arg = this;
                task->cl_arg_size = sizeof(this);

                this->submit_starpu_task(task);
            }
        }
    }

    void init_M2P_codelet() {
        starpu_codelet_init(&this->M2P_cl);
        this->M2P_cl.name = "M2P";
        this->M2P_cl.cpu_funcs[0] = M2P_cpu;
        this->M2P_cl.cpu_funcs_name[0] = {const_cast<char*>("M2P task")};
        this->M2P_cl.nbuffers = 4; // {node symb, node multipole, leaf symb, leaf particles}
        this->M2P_cl.modes[0] = STARPU_R;
        this->M2P_cl.modes[1] = STARPU_R;
        this->M2P_cl.modes[2] = STARPU_R;
        this->M2P_cl.modes[3] = STARPU_RW;
    }

    static void M2P_cpu(void** buffers, void* cl_arg) {
        symbolic_data_t* node_symb_data =
            (symbolic_data_t*)STARPU_VARIABLE_GET_PTR(buffers[0]);

        multipole_t* node_multipole =
            (multipole_t*) STARPU_VARIABLE_GET_PTR(buffers[1]);

        symbolic_data_t* leaf_symb_data =
            (symbolic_data_t*)STARPU_VARIABLE_GET_PTR(buffers[2]);

        container_t* leaf_particles =
            (container_t*) STARPU_VARIABLE_GET_PTR(buffers[3]);

        auto algo = (FAdaptiveStarPU<tree_t, kernel_t>*) cl_arg;
        kernel_t* ker = get_starpu_worker_kernel(algo);
        ker->M2P(node_multipole, node_symb_data,
                 leaf_particles, leaf_symb_data);
    }


    // L2P
    void down_to_target() {
        for(node_t* leaf : _tree.leaves()) {
            if(leaf->getParticleContainer()->size() == 0) {
                continue;
            }

            node_data_handles& leaf_handles = this->_data_handles.at(leaf);

            starpu_task* task = starpu_task_create();
            task->cl = &(this->L2P_cl);
            task->handles[0] = leaf_handles.local_exp;
            task->handles[1] = leaf_handles.symbolic;
            task->handles[2] = leaf_handles.particles;
            task->cl_arg = this;
            task->cl_arg_size = sizeof(this);

            this->submit_starpu_task(task);
        }
    }

    void init_L2P_codelet() {
        starpu_codelet_init(&this->L2P_cl);
        this->L2P_cl.name = "L2P";
        this->L2P_cl.cpu_funcs[0] = L2P_cpu;
        this->L2P_cl.cpu_funcs_name[0] = {const_cast<char*>("L2P task")};
        this->L2P_cl.nbuffers = 3; // {0: leaf local exp, 1: leaf symb, 2: leaf particles}
        this->L2P_cl.modes[0] = STARPU_R;
        this->L2P_cl.modes[1] = STARPU_R;
        this->L2P_cl.modes[2] = STARPU_RW;
    }

    static void L2P_cpu(void** buffers, void* cl_arg) {
        local_expansion_t* leaf_local_exp =
            (local_expansion_t*) STARPU_VARIABLE_GET_PTR(buffers[0]);

        symbolic_data_t* leaf_symb_data =
            (symbolic_data_t*) STARPU_VARIABLE_GET_PTR(buffers[1]);

        container_t* particle_container =
            (container_t*) STARPU_VARIABLE_GET_PTR(buffers[2]);

        auto algo = (FAdaptiveStarPU<tree_t, kernel_t>*) cl_arg;
        kernel_t* ker = get_starpu_worker_kernel(algo);
        ker->L2P(leaf_local_exp, leaf_symb_data, particle_container);
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
        for(node_t* leaf : _tree.leaves()) {
            auto* const leaf_source_particle_container =
                leaf->getParticleContainer();
            auto* const leaf_target_particle_container =
                leaf->getParticleContainer();
            // Skip empty leaves
            if( leaf_source_particle_container->size() == 0
                && leaf_target_particle_container->size() == 0) {
                continue;
            }

            starpu_task* task = starpu_task_create();
            task->cl = &(this->P2P_cl);
            task->cl_arg = this;
            task->cl_arg_size = sizeof(this);
            task->nbuffers = 3;

            auto pick_item = [leaf](const node_t* u){
                return u != leaf && u->getParticleContainer()->size() != 0;
            };

            for(node_t* u_item : leaf->U) {
                // The kernels do not consider leaf to be adjacent to itself
                if(! pick_item(u_item)) {
                    continue;
                }
                task->nbuffers += 2;
            }

            task->dyn_handles = (starpu_data_handle_t*)
                std::malloc(task->nbuffers * sizeof(starpu_data_handle_t));
            task->dyn_modes = (starpu_data_access_mode*)
                std::malloc(task->nbuffers * sizeof(starpu_data_access_mode));

            task->dyn_handles[0] = this->_data_handles.at(leaf).symbolic;
            task->dyn_handles[1] = this->_data_handles.at(leaf).particles;
            task->dyn_handles[2] = this->_data_handles.at(leaf).particles;

            task->dyn_modes[0] = STARPU_R;
            task->dyn_modes[1] = STARPU_RW;
            task->dyn_modes[2] = STARPU_R;

            std::size_t j = 3;
            for(node_t* u_item : leaf->U) {
                if(! pick_item(u_item)) {
                    continue;
                }
                task->dyn_handles[j+0] = this->_data_handles.at(u_item).symbolic;
                task->dyn_handles[j+1] = this->_data_handles.at(u_item).particles;
                task->dyn_modes[j+0] = STARPU_R;
                task->dyn_modes[j+1] = STARPU_RW;
                j += 2;
            }

            this->submit_starpu_task(task);
        }
    }

    void init_P2P_codelet() {
        starpu_codelet_init(&this->P2P_cl);
        this->P2P_cl.name = "P2P";
        this->P2P_cl.cpu_funcs[0] = P2P_cpu;
        this->P2P_cl.cpu_funcs_name[0] = {const_cast<char*>("P2P task")};
        this->P2P_cl.nbuffers = STARPU_VARIABLE_NBUFFERS;
    }

    static void P2P_cpu(void** buffers, void* cl_arg) {
        const starpu_task* const current_task =  starpu_task_get_current();
        std::size_t buffer_count = current_task->nbuffers;

        const symbolic_data_t* node_symbolic =
            (symbolic_data_t*) STARPU_VARIABLE_GET_PTR(buffers[0]);
        container_t* node_target_particles =
            (container_t*) STARPU_VARIABLE_GET_PTR(buffers[1]);
        container_t* node_source_particles =
            (container_t*) STARPU_VARIABLE_GET_PTR(buffers[2]);

        std::vector<container_t*> u_source_particles;
        std::vector<int> offset_indices;
        // two buffers per V item, first is for the local exp
        u_source_particles.reserve((buffer_count-3) / 2);
        offset_indices.reserve((buffer_count-3) / 2);

        for(std::size_t i = 3; i < buffer_count; i += 2) {
            const symbolic_data_t* u_item_symb = (symbolic_data_t*) STARPU_VARIABLE_GET_PTR(buffers[i]);
            offset_indices.emplace_back(compute_box_offset_index(*node_symbolic, *u_item_symb, 1));
            u_source_particles.emplace_back((container_t*) STARPU_VARIABLE_GET_PTR(buffers[i+1]));
        }

        auto algo = (FAdaptiveStarPU<tree_t, kernel_t>*) cl_arg;
        kernel_t* ker = get_starpu_worker_kernel(algo);
        ker->P2P(
            FTreeCoordinate(node_symbolic->m_idx),
            node_target_particles, node_source_particles,
            u_source_particles.data(),
            offset_indices.data(),
            static_cast<int>(u_source_particles.size()));
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
    static int compute_box_offset_index(
        const typename node_t::symbolic_data_t& node,
        const typename node_t::symbolic_data_t& other_node,
        const std::size_t n)
    {
        MortonIndex node_m_idx = node.m_idx;
        MortonIndex other_node_m_idx = other_node.m_idx;
        if(node.depth > other_node.depth) {
            node_m_idx >>= ((node.depth - other_node.depth) * node_t::Dim);
        } else if(node.depth < other_node.depth) {
            other_node_m_idx >>= ((other_node.depth - node.depth) * node_t::Dim);
        }
	FTreeCoordinate node_offset = FTreeCoordinate(other_node_m_idx) - FTreeCoordinate(node_m_idx);
	//        FTreeCoordinate node_offset(FTreeCoordinate(other_node_m_idx) - FTreeCoordinate(node_m_idx));
        std::size_t other_node_index = 0;
        for(std::size_t i = 0; i < node_t::Dim; ++i) {
            other_node_index = other_node_index * (2 * n + 1)
                + node_offset[i]
                + n;
        }
        return static_cast<int>(other_node_index);
    }


};

#endif
