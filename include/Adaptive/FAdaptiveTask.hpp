#ifndef SCALFMM_ADAPTIVE_TASK_ALGO_HPP_
#define SCALFMM_ADAPTIVE_TASK_ALGO_HPP_

#include <algorithm>
#include <numeric>
#include <cmath> // Used to round box differences
#include <functional>
#include <map>
#include <vector>
#include <list>
#include <array>
#include <type_traits>
#include <memory>
#include <sstream>

#include <omp.h>

#include <unistd.h>
#include "Utils/FGlobal.hpp"
#include "Core/FCoreCommon.hpp"
#include "Containers/FTreeCoordinate.hpp"
#include "Utils/FAlgorithmTimers.hpp"

#include "Kernels/FKernelConcepts.hpp"

#include "kernel_utilities.hpp"
#include "inria/logic.hpp"
template<class _Tree, class _Kernel,
         class = inria::require<
             scalfmm::meta::adaptive_compatible<_Tree,_Kernel>
             >
         >
class FAdaptiveTask : public FAlgorithmInterface, public FAlgorithmTimers {
public:
    using tree_t = _Tree;
    using kernel_t = _Kernel;
    using FReal = typename tree_t::FReal;

    using node_t = typename tree_t::node_t;

    using multipole_t = typename node_t::data_t::multipole_t;
    using local_expansion_t = typename node_t::data_t::local_expansion_t;
    using symbolic_data_t = typename node_t::symbolic_data_t;

protected:
    /// Tree
    tree_t&   _tree;
    /// Vector of kernels, one per thread
    std::vector<std::unique_ptr<kernel_t>> _kernels;

    /**
     * \brief Data pool for mock dependency
     *
     * Generates addresses that can be used as OpenMP task dependencies. To do
     * so, a list of char arrays (buckets) is created. Each char is a potential
     * dependency. When an array is exhausted, an new one is appended to the
     * list.
     */
    struct mock_dependency {
        /// Size of a bucket
        enum : std::size_t {bucket_size = 512};
        /// List of buckets. A bucket is implented as an array
        std::list<std::array<char, bucket_size>> pool;
        /// Index of the next dependency to hand out in the current bucket
        std::size_t mock_dependency_index = bucket_size;

        /// Yield a new dependency
        char* next() {
            // When the bucket is empty, add a new one
            if(bucket_size == mock_dependency_index) {
                pool.emplace_back();
                mock_dependency_index = 0;
            }
            return (char*) (&(this->pool.back()[mock_dependency_index++]));
        }
    } mock_dep;

    /**
     * \brief Dependency type names structure
     *
     * Each member of the enum identifies a dependency type name.
     */
    enum class dep_t {
        M,   ///< Multipole
        L,   ///< Local expansion
        P_s, ///< Source particle
        P_t, ///< Target particle
    };


    /**
     * \brief Get pointer to node ressource handling given dependency type
     *
     * OpenMP task dependencies do not need to be correlated with real data,
     * they only need to point to the same area in memory. This method maps a
     * dependency type to an address inside the given node object, cast to
     * `const char`.
     *
     * \tparam I Dependency type
     *
     * \param node Node to which the dependency applies
     *
     * \return A pointer to memory inside the node.
     */
    template<dep_t I>
    const char* get_dependency(node_t* node) const {
        if(I == dep_t::M) {
            return (const char*) node->getData();
        } else if(I == dep_t::L) {
            using type = typename std::decay<decltype((node->getSymbolicData()))>::type;
            return (const char*) const_cast<type*>(&(node->getSymbolicData()));
        } else if(I == dep_t::P_s) {
            return (const char*) node;
        } else if(I == dep_t::P_t) {
            return (const char*) node->getParticleContainer();
        }
        return nullptr;
    }


public:

    /**
     * \brief Constructor
     *
     * \param tree The tree to run over
     * \param kernel The kernel to use
     */
  FAdaptiveTask(tree_t* tree, kernel_t* kernel) :
            _tree(*tree)
    {
        // Place kernel objects near their threads in memory to avoid NUMA
        // latency
        int threads = 1;
        #pragma omp parallel
        #pragma omp master 
          threads = omp_get_num_threads();

        _kernels.resize(threads);

        #pragma omp parallel num_threads(threads)
        {
          #pragma omp critical
          {
            _kernels[omp_get_thread_num()] = std::make_unique<kernel_t>(*kernel);
          }
        }
    }

    std::string name() const override {
        return "Task adaptive algorithm";
    }

    std::string description() const override {
        int threads = 1;
        #pragma omp parallel shared(threads)
        {
            #pragma omp single nowait
            {
                threads = omp_get_num_threads();
            }
        }
        return std::string("threads: ") + std::to_string(threads);
    }

    using FAlgorithmInterface::execute;

    /**
     * \brief Run specific steps of the algorithm
     *
     * \param operations Specifies the algorithm operations to run, see
     * FFmmOperations.
     */
    void execute(const unsigned int operations) override {
      //    this->run(operations);
    }

    void run(int operations) {

        #pragma omp parallel
        {
            scalfmm::setup_kernel(*(this->_kernels[omp_get_thread_num()].get()),
                                  this->_tree);
            #pragma omp barrier

            #pragma omp single nowait
            {
                if(operations & FFmmP2M) {
                    // 1. source to up, P2M
                    Timers["P2M"].tic();
                    this->source_to_up();;
                    Timers["P2M"].tac();
                    std::cout << "    P2M: " << Timers["P2M"].elapsed() << std::endl;
                }

                if(operations & FFmmP2L) {
                    // 3b X-list, P2L
                    Timers["P2L"].tic();
                    this->x_list_step();
                    Timers["P2L"].tac();
                    std::cout << "    P2L: " << Timers["P2L"].elapsed() << std::endl;
                }

                if(operations & FFmmP2P) {
                    // A. U-list, P2P
                    Timers["P2P"].tic();
                    this->u_list_step();
                    Timers["P2P"].tac();
                    std::cout << "    P2P: " << Timers["P2P"].elapsed() << std::endl;
                }

                if(operations & FFmmM2M) {
                    // 2. up to up, M2M
                    Timers["M2M"].tic();
                    this->up_to_up();
                    Timers["M2M"].tac();
                    std::cout << "    M2M: " << Timers["M2M"].elapsed() << std::endl;
                }

                if(operations & FFmmM2P) {
                    // 5a W-list, M2P
                    Timers["M2P"].tic();
                    this->w_list_step();
                    Timers["M2P"].tac();
                    std::cout << "    M2P: " << Timers["M2P"].elapsed() << std::endl;
                }

                if(operations & FFmmM2L) {
                    // 3a V-list, M2L
                    Timers["M2L"].tic();
                    this->v_list_step();
                    Timers["M2L"].tac();
                    std::cout << "    M2L: " << Timers["M2L"].elapsed() << std::endl;
                }

                if(operations & FFmmL2L) {
                    // 4. down to down, L2L
                    Timers["L2L"].tic();
                    this->down_to_down();
                    Timers["L2L"].tac();
                    std::cout << "    L2L: " << Timers["L2L"].elapsed() << std::endl;
                }

                if(operations & FFmmL2P) {
                    // 5b down to target, L2P
                    Timers["L2P"].tic();
                    this->down_to_target();
                    Timers["L2P"].tac();
                    std::cout << "    L2P: " << Timers["L2P"].elapsed() << std::endl;
                }



                std::cout << "    tasks creation: "
                          << std::accumulate(Timers.begin(), Timers.end(), 0.,
                                             [](double & res, const typename decltype(Timers)::value_type & t)
                                             {return res + t.second.cumulated();})
                          << '\n';
            }
        }
    }

    /** \brief Walk through leaves and queue P2M tasks */
    void source_to_up() {
        for(node_t* leaf : _tree.leaves()) {
            this->create_P2M_task(leaf);
        }
        // #pragma omp taskwait
    }

    /**
     * \brief Create and queue a P2M task
     *
     * \param leaf The P2M leaf
     */
    void create_P2M_task(node_t* leaf) {
            const char* ps_dep = get_dependency<dep_t::P_s>(leaf);(void)ps_dep;
            const char* m_dep = get_dependency<dep_t::M>(leaf);(void)m_dep;

            #pragma omp task firstprivate(leaf)                        \
                depend(in: ps_dep[:1])                                 \
                depend(inout: m_dep[:1])
            {
                const int thread_num = omp_get_thread_num();
                _kernels[thread_num]->P2M(
                    &(leaf->getData()->getMultipoleData()),
                    &(leaf->getSymbolicData()),
                    leaf->getParticleContainer());
            }
    }


    /** \brief Walk through tree and queue M2M tasks */
    void up_to_up() {
        for(node_t& n : _tree.post_order_walk()) {
            if(! n.is_leaf()) {
                create_M2M_task(&n);
            }
        }
        // #pragma omp taskwait
    }

    /**
     * \brief Create and queue an M2M task
     *
     * \param node An internal node
     *
     * \warning node is assumed not to be a leaf.
     */
    void create_M2M_task(node_t* node) {
        assert(! node->is_leaf());
        // Setup task dependencies
        // children data
        const char* children_dep[node_t::child_count] = {};
        for(node_t* child : node->getChildren()) {
            children_dep[child->getIndex() & (node_t::child_count-1)]
                = get_dependency<dep_t::M>(child);
        }
        // node data
        const char* parent_dep = get_dependency<dep_t::M>(node); (void) parent_dep;

        #pragma omp task                        \
            depend(in:                          \
                   children_dep[0][:1],         \
                   children_dep[1][:1],         \
                   children_dep[2][:1],         \
                   children_dep[3][:1],         \
                   children_dep[4][:1],         \
                   children_dep[5][:1],         \
                   children_dep[6][:1],         \
                   children_dep[7][:1])         \
            depend(out: parent_dep[:1])
        {
            const int thread_num = omp_get_thread_num();

                std::array<multipole_t*, node_t::child_count> child_multipoles {};
                std::transform(std::begin(node->getChildren()), std::end(node->getChildren()),
                               child_multipoles.begin(),
                               node_t::template getMultipoleDataFromNode<node_t, multipole_t>);

                std::array<symbolic_data_t*, node_t::child_count> child_symbolics {};
                std::transform(std::begin(node->getChildren()), std::end(node->getChildren()),
                               child_symbolics.begin(),
                               node_t::template getSymbolicData<node_t, symbolic_data_t>);

            // Call kernel module
                this->_kernels[thread_num]->M2M(
                    &(node->getData()->getMultipoleData()),
                    &(node->getSymbolicData()),
                    child_multipoles.data(),
                    child_symbolics.data()
                    );
        }
    }



    /** \brief Walk through tree and queue M2L tasks */
    void v_list_step() {
        for(node_t& n : _tree.in_order_walk()) {
            if(n.is_leaf() && n.getParticleContainer()->size() == 0) {
                continue;
            }
            create_M2L_task(&n);
        }
        // #pragma omp taskwait
    }


    /**
     * \brief Create and queue an M2L task
     *
     * \param node A tree node
     */
    void create_M2L_task(node_t* node) {

        // Generate task dependencies
        // There cannot be more than 7^Dim cells involved in a M2L, in 3D, this is 343
        const char* task_deps[343];
        const char* data_dep = get_dependency<dep_t::L>(node);(void) data_dep;
        std::size_t idx_dep = 0;
        // Add existing dependencies
        for(node_t* v_item : node->V) {
            if(v_item->is_leaf()
               && v_item->getParticleContainer()->size() == 0) {
                continue;
            }
            task_deps[idx_dep] = get_dependency<dep_t::M>(v_item);
            ++idx_dep;
        }
        // Add mock dependencies, these are generated on the fly and used
        // only once, that way they can never stop a task from starting
        while(idx_dep < 343) {
            task_deps[idx_dep] = this->mock_dep.next();
            ++idx_dep;
        }

        #pragma omp task                                                \
            depend(in:                                                  \
                   task_deps[0][:1],                                    \
                   task_deps[1][:1],                                    \
                   task_deps[2][:1],                                    \
                   task_deps[3][:1],                                    \
                   task_deps[4][:1],                                    \
                   task_deps[5][:1],                                    \
                   task_deps[6][:1],                                    \
                   task_deps[7][:1],                                    \
                   task_deps[8][:1],                                    \
                   task_deps[9][:1],                                    \
                   task_deps[10][:1],                                   \
                   task_deps[11][:1],                                   \
                   task_deps[12][:1],                                   \
                   task_deps[13][:1],                                   \
                   task_deps[14][:1],                                   \
                   task_deps[15][:1],                                   \
                   task_deps[16][:1],                                   \
                   task_deps[17][:1],                                   \
                   task_deps[18][:1],                                   \
                   task_deps[19][:1],                                   \
                   task_deps[20][:1],                                   \
                   task_deps[21][:1],                                   \
                   task_deps[22][:1],                                   \
                   task_deps[23][:1],                                   \
                   task_deps[24][:1],                                   \
                   task_deps[25][:1],                                   \
                   task_deps[26][:1],                                   \
                   task_deps[27][:1],                                   \
                   task_deps[28][:1],                                   \
                   task_deps[29][:1],                                   \
                   task_deps[30][:1],                                   \
                   task_deps[31][:1],                                   \
                   task_deps[32][:1],                                   \
                   task_deps[33][:1],                                   \
                   task_deps[34][:1],                                   \
                   task_deps[35][:1],                                   \
                   task_deps[36][:1],                                   \
                   task_deps[37][:1],                                   \
                   task_deps[38][:1],                                   \
                   task_deps[39][:1],                                   \
                   task_deps[40][:1],                                   \
                   task_deps[41][:1],                                   \
                   task_deps[42][:1],                                   \
                   task_deps[43][:1],                                   \
                   task_deps[44][:1],                                   \
                   task_deps[45][:1],                                   \
                   task_deps[46][:1],                                   \
                   task_deps[47][:1],                                   \
                   task_deps[48][:1],                                   \
                   task_deps[49][:1],                                   \
                   task_deps[50][:1],                                   \
                   task_deps[51][:1],                                   \
                   task_deps[52][:1],                                   \
                   task_deps[53][:1],                                   \
                   task_deps[54][:1],                                   \
                   task_deps[55][:1],                                   \
                   task_deps[56][:1],                                   \
                   task_deps[57][:1],                                   \
                   task_deps[58][:1],                                   \
                   task_deps[59][:1],                                   \
                   task_deps[60][:1],                                   \
                   task_deps[61][:1],                                   \
                   task_deps[62][:1],                                   \
                   task_deps[63][:1],                                   \
                   task_deps[64][:1],                                   \
                   task_deps[65][:1],                                   \
                   task_deps[66][:1],                                   \
                   task_deps[67][:1],                                   \
                   task_deps[68][:1],                                   \
                   task_deps[69][:1],                                   \
                   task_deps[70][:1],                                   \
                   task_deps[71][:1],                                   \
                   task_deps[72][:1],                                   \
                   task_deps[73][:1],                                   \
                   task_deps[74][:1],                                   \
                   task_deps[75][:1],                                   \
                   task_deps[76][:1],                                   \
                   task_deps[77][:1],                                   \
                   task_deps[78][:1],                                   \
                   task_deps[79][:1],                                   \
                   task_deps[80][:1],                                   \
                   task_deps[81][:1],                                   \
                   task_deps[82][:1],                                   \
                   task_deps[83][:1],                                   \
                   task_deps[84][:1],                                   \
                   task_deps[85][:1],                                   \
                   task_deps[86][:1],                                   \
                   task_deps[87][:1],                                   \
                   task_deps[88][:1],                                   \
                   task_deps[89][:1],                                   \
                   task_deps[90][:1],                                   \
                   task_deps[91][:1],                                   \
                   task_deps[92][:1],                                   \
                   task_deps[93][:1],                                   \
                   task_deps[94][:1],                                   \
                   task_deps[95][:1],                                   \
                   task_deps[96][:1],                                   \
                   task_deps[97][:1],                                   \
                   task_deps[98][:1],                                   \
                   task_deps[99][:1],                                   \
                   task_deps[100][:1],                                  \
                   task_deps[101][:1],                                  \
                   task_deps[102][:1],                                  \
                   task_deps[103][:1],                                  \
                   task_deps[104][:1],                                  \
                   task_deps[105][:1],                                  \
                   task_deps[106][:1],                                  \
                   task_deps[107][:1],                                  \
                   task_deps[108][:1],                                  \
                   task_deps[109][:1],                                  \
                   task_deps[110][:1],                                  \
                   task_deps[111][:1],                                  \
                   task_deps[112][:1],                                  \
                   task_deps[113][:1],                                  \
                   task_deps[114][:1],                                  \
                   task_deps[115][:1],                                  \
                   task_deps[116][:1],                                  \
                   task_deps[117][:1],                                  \
                   task_deps[118][:1],                                  \
                   task_deps[119][:1],                                  \
                   task_deps[120][:1],                                  \
                   task_deps[121][:1],                                  \
                   task_deps[122][:1],                                  \
                   task_deps[123][:1],                                  \
                   task_deps[124][:1],                                  \
                   task_deps[125][:1],                                  \
                   task_deps[126][:1],                                  \
                   task_deps[127][:1],                                  \
                   task_deps[128][:1],                                  \
                   task_deps[129][:1],                                  \
                   task_deps[130][:1],                                  \
                   task_deps[131][:1],                                  \
                   task_deps[132][:1],                                  \
                   task_deps[133][:1],                                  \
                   task_deps[134][:1],                                  \
                   task_deps[135][:1],                                  \
                   task_deps[136][:1],                                  \
                   task_deps[137][:1],                                  \
                   task_deps[138][:1],                                  \
                   task_deps[139][:1],                                  \
                   task_deps[140][:1],                                  \
                   task_deps[141][:1],                                  \
                   task_deps[142][:1],                                  \
                   task_deps[143][:1],                                  \
                   task_deps[144][:1],                                  \
                   task_deps[145][:1],                                  \
                   task_deps[146][:1],                                  \
                   task_deps[147][:1],                                  \
                   task_deps[148][:1],                                  \
                   task_deps[149][:1],                                  \
                   task_deps[150][:1],                                  \
                   task_deps[151][:1],                                  \
                   task_deps[152][:1],                                  \
                   task_deps[153][:1],                                  \
                   task_deps[154][:1],                                  \
                   task_deps[155][:1],                                  \
                   task_deps[156][:1],                                  \
                   task_deps[157][:1],                                  \
                   task_deps[158][:1],                                  \
                   task_deps[159][:1],                                  \
                   task_deps[160][:1],                                  \
                   task_deps[161][:1],                                  \
                   task_deps[162][:1],                                  \
                   task_deps[163][:1],                                  \
                   task_deps[164][:1],                                  \
                   task_deps[165][:1],                                  \
                   task_deps[166][:1],                                  \
                   task_deps[167][:1],                                  \
                   task_deps[168][:1],                                  \
                   task_deps[169][:1],                                  \
                   task_deps[170][:1],                                  \
                   task_deps[171][:1],                                  \
                   task_deps[172][:1],                                  \
                   task_deps[173][:1],                                  \
                   task_deps[174][:1],                                  \
                   task_deps[175][:1],                                  \
                   task_deps[176][:1],                                  \
                   task_deps[177][:1],                                  \
                   task_deps[178][:1],                                  \
                   task_deps[179][:1],                                  \
                   task_deps[180][:1],                                  \
                   task_deps[181][:1],                                  \
                   task_deps[182][:1],                                  \
                   task_deps[183][:1],                                  \
                   task_deps[184][:1],                                  \
                   task_deps[185][:1],                                  \
                   task_deps[186][:1],                                  \
                   task_deps[187][:1],                                  \
                   task_deps[188][:1],                                  \
                   task_deps[189][:1],                                  \
                   task_deps[190][:1],                                  \
                   task_deps[191][:1],                                  \
                   task_deps[192][:1],                                  \
                   task_deps[193][:1],                                  \
                   task_deps[194][:1],                                  \
                   task_deps[195][:1],                                  \
                   task_deps[196][:1],                                  \
                   task_deps[197][:1],                                  \
                   task_deps[198][:1],                                  \
                   task_deps[199][:1],                                  \
                   task_deps[200][:1],                                  \
                   task_deps[201][:1],                                  \
                   task_deps[202][:1],                                  \
                   task_deps[203][:1],                                  \
                   task_deps[204][:1],                                  \
                   task_deps[205][:1],                                  \
                   task_deps[206][:1],                                  \
                   task_deps[207][:1],                                  \
                   task_deps[208][:1],                                  \
                   task_deps[209][:1],                                  \
                   task_deps[210][:1],                                  \
                   task_deps[211][:1],                                  \
                   task_deps[212][:1],                                  \
                   task_deps[213][:1],                                  \
                   task_deps[214][:1],                                  \
                   task_deps[215][:1],                                  \
                   task_deps[216][:1],                                  \
                   task_deps[217][:1],                                  \
                   task_deps[218][:1],                                  \
                   task_deps[219][:1],                                  \
                   task_deps[220][:1],                                  \
                   task_deps[221][:1],                                  \
                   task_deps[222][:1],                                  \
                   task_deps[223][:1],                                  \
                   task_deps[224][:1],                                  \
                   task_deps[225][:1],                                  \
                   task_deps[226][:1],                                  \
                   task_deps[227][:1],                                  \
                   task_deps[228][:1],                                  \
                   task_deps[229][:1],                                  \
                   task_deps[230][:1],                                  \
                   task_deps[231][:1],                                  \
                   task_deps[232][:1],                                  \
                   task_deps[233][:1],                                  \
                   task_deps[234][:1],                                  \
                   task_deps[235][:1],                                  \
                   task_deps[236][:1],                                  \
                   task_deps[237][:1],                                  \
                   task_deps[238][:1],                                  \
                   task_deps[239][:1],                                  \
                   task_deps[240][:1],                                  \
                   task_deps[241][:1],                                  \
                   task_deps[242][:1],                                  \
                   task_deps[243][:1],                                  \
                   task_deps[244][:1],                                  \
                   task_deps[245][:1],                                  \
                   task_deps[246][:1],                                  \
                   task_deps[247][:1],                                  \
                   task_deps[248][:1],                                  \
                   task_deps[249][:1],                                  \
                   task_deps[250][:1],                                  \
                   task_deps[251][:1],                                  \
                   task_deps[252][:1],                                  \
                   task_deps[253][:1],                                  \
                   task_deps[254][:1],                                  \
                   task_deps[255][:1],                                  \
                   task_deps[256][:1],                                  \
                   task_deps[257][:1],                                  \
                   task_deps[258][:1],                                  \
                   task_deps[259][:1],                                  \
                   task_deps[260][:1],                                  \
                   task_deps[261][:1],                                  \
                   task_deps[262][:1],                                  \
                   task_deps[263][:1],                                  \
                   task_deps[264][:1],                                  \
                   task_deps[265][:1],                                  \
                   task_deps[266][:1],                                  \
                   task_deps[267][:1],                                  \
                   task_deps[268][:1],                                  \
                   task_deps[269][:1],                                  \
                   task_deps[270][:1],                                  \
                   task_deps[271][:1],                                  \
                   task_deps[272][:1],                                  \
                   task_deps[273][:1],                                  \
                   task_deps[274][:1],                                  \
                   task_deps[275][:1],                                  \
                   task_deps[276][:1],                                  \
                   task_deps[277][:1],                                  \
                   task_deps[278][:1],                                  \
                   task_deps[279][:1],                                  \
                   task_deps[280][:1],                                  \
                   task_deps[281][:1],                                  \
                   task_deps[282][:1],                                  \
                   task_deps[283][:1],                                  \
                   task_deps[284][:1],                                  \
                   task_deps[285][:1],                                  \
                   task_deps[286][:1],                                  \
                   task_deps[287][:1],                                  \
                   task_deps[288][:1],                                  \
                   task_deps[289][:1],                                  \
                   task_deps[290][:1],                                  \
                   task_deps[291][:1],                                  \
                   task_deps[292][:1],                                  \
                   task_deps[293][:1],                                  \
                   task_deps[294][:1],                                  \
                   task_deps[295][:1],                                  \
                   task_deps[296][:1],                                  \
                   task_deps[297][:1],                                  \
                   task_deps[298][:1],                                  \
                   task_deps[299][:1],                                  \
                   task_deps[300][:1],                                  \
                   task_deps[301][:1],                                  \
                   task_deps[302][:1],                                  \
                   task_deps[303][:1],                                  \
                   task_deps[304][:1],                                  \
                   task_deps[305][:1],                                  \
                   task_deps[306][:1],                                  \
                   task_deps[307][:1],                                  \
                   task_deps[308][:1],                                  \
                   task_deps[309][:1],                                  \
                   task_deps[310][:1],                                  \
                   task_deps[311][:1],                                  \
                   task_deps[312][:1],                                  \
                   task_deps[313][:1],                                  \
                   task_deps[314][:1],                                  \
                   task_deps[315][:1],                                  \
                   task_deps[316][:1],                                  \
                   task_deps[317][:1],                                  \
                   task_deps[318][:1],                                  \
                   task_deps[319][:1],                                  \
                   task_deps[320][:1],                                  \
                   task_deps[321][:1],                                  \
                   task_deps[322][:1],                                  \
                   task_deps[323][:1],                                  \
                   task_deps[324][:1],                                  \
                   task_deps[325][:1],                                  \
                   task_deps[326][:1],                                  \
                   task_deps[327][:1],                                  \
                   task_deps[328][:1],                                  \
                   task_deps[329][:1],                                  \
                   task_deps[330][:1],                                  \
                   task_deps[331][:1],                                  \
                   task_deps[332][:1],                                  \
                   task_deps[333][:1],                                  \
                   task_deps[334][:1],                                  \
                   task_deps[335][:1],                                  \
                   task_deps[336][:1],                                  \
                   task_deps[337][:1],                                  \
                   task_deps[338][:1],                                  \
                   task_deps[339][:1],                                  \
                   task_deps[340][:1],                                  \
                   task_deps[341][:1],                                  \
                   task_deps[342][:1]                                   \
                ) depend(inout: data_dep[:1]) firstprivate(node)
        {

            const int thread_num = omp_get_thread_num();

            std::vector<multipole_t*> v_item_multipoles;
            std::vector<symbolic_data_t*> v_item_symbolics;
            std::vector<int> v_item_indices;

            // Needed to compute offset between boxes
            for(node_t* v_item : node->V) {
                if(v_item->is_leaf()
                   && v_item->getParticleContainer()->size() == 0) {
                    continue;
                }
                v_item_indices.push_back(compute_box_offset_index(node, v_item, 3));
                v_item_multipoles.push_back(&(v_item->getData()->getMultipoleData()));
                v_item_symbolics.push_back(&(v_item->getSymbolicData()));
            }

            // Call kernel M2L operator
            this->_kernels[thread_num]->M2L(
                &(node->getData()->getLocalExpansionData()),
                &(node->getSymbolicData()),
                v_item_multipoles.data(),
                v_item_symbolics.data(),
                v_item_indices.data(),
                static_cast<int>(v_item_multipoles.size())
                );
        }
    }

    /** \brief Walk through leaves and queue P2L tasks */
    void x_list_step() {
        /* NOTE: the X list and W list are complementary: if A is in X(B) then B
         * is in W(A).
         *
         * We loop over the leaves first to detect the empty ones early on.
         */
        for(node_t* leaf : _tree.leaves()) {
            if(leaf->getParticleContainer()->size() > 0) {
                this->create_P2L_task(leaf);
            }
        }
        // #pragma omp taskwait
    }


    /**
     * \brief Create and queue a P2L task
     *
     * \param leaf A tree leaf
     */
    void create_P2L_task(node_t* leaf) {
        for(node_t* w_item : leaf->W) {
            if(w_item->is_leaf() && w_item->getParticleContainer()->size() == 0) {
                continue;
            }

            const char* w_dep = get_dependency<dep_t::L>(w_item);(void)w_dep;
            const char* ps_dep = get_dependency<dep_t::P_s>(leaf);(void)ps_dep;

            #pragma omp task firstprivate(leaf, w_item) depend(in: ps_dep[:1]) depend(inout: w_dep[:1])
            {
                const int thread_num = omp_get_thread_num();
                this->_kernels[thread_num]->P2L(
                    &(w_item->getData()->getLocalExpansionData()),
                    &(w_item->getSymbolicData()),
                    leaf->getParticleContainer(),
                    &(leaf->getSymbolicData())
                    );
            }
        }
    }

    /** Walk through tree and queue L2L tasks */
    void down_to_down() {
        for(node_t& n : _tree.pre_order_walk()) {
            if(! n.is_leaf()) {
                this->create_L2L_task(&n);
            }
        }
        // #pragma omp taskwait
    }


    /**
     * \brief Create and queue an L2L task
     *
     * \param node An internal tree node
     *
     * \node node is assumed not to be a leaf
     */
    void create_L2L_task(node_t* node) {
        assert(! node->is_leaf());

        const char* parent_dep = get_dependency<dep_t::L>(node);(void)parent_dep;
        const char* children_dep[node_t::child_count];
        for(std::size_t i = 0; i < node_t::child_count; ++i) {
            children_dep[i] = get_dependency<dep_t::L>(node->getChild(i));
        }

        #pragma omp task                        \
            depend(in: parent_dep[:1])          \
            depend(inout:                       \
                   children_dep[0][:1],         \
                   children_dep[1][:1],         \
                   children_dep[2][:1],         \
                   children_dep[3][:1],         \
                   children_dep[4][:1],         \
                   children_dep[5][:1],         \
                   children_dep[6][:1],         \
                   children_dep[7][:1]          \
                )
        {
            const int thread_num = omp_get_thread_num();

                std::array<local_expansion_t*, node_t::child_count> child_local_expansions {};
                std::transform(std::begin(node->getChildren()), std::end(node->getChildren()),
                               child_local_expansions.begin(),
                               node_t::template getLocalExpansionDataFromNode<node_t, local_expansion_t>);
                std::array<symbolic_data_t*, node_t::child_count> child_symbolics {};
                std::transform(std::begin(node->getChildren()), std::end(node->getChildren()),
                               child_symbolics.begin(),
                               node_t::template getSymbolicDataFromNode<node_t, symbolic_data_t>);

                this->_kernels[thread_num]->L2L(
                    &(node->getData()->getLocalExpansionData()),
                    &(node->getSymbolicData()),
                    child_local_expansions.data(),
                    child_symbolics.data()
                    );
        }
    }



    /** \brief Walk through the leaves an queue M2P tasks */
    void w_list_step() {
        for(node_t* leaf : _tree.leaves()) {
            if(leaf->getParticleContainer()->size() > 0) {
                create_M2P_task(leaf);
            }
        }
        // #pragma omp taskwait
    }

    /**
     * \brief Create and queue an M2P task
     *
     * \param leaf A leaf
     */
    void create_M2P_task(node_t* leaf) {
        for(node_t* w_item : leaf->W) {
            if(w_item->is_leaf() && w_item->getParticleContainer()->size() == 0) {
                continue;
            }

            const char* w_dep = get_dependency<dep_t::M>(w_item); (void)w_dep;
            const char* pt_dep = get_dependency<dep_t::P_t>(leaf); (void)pt_dep;

            #pragma omp task depend(inout: pt_dep[:1]) depend(in: w_dep[:1])
            {
                const int thread_num = omp_get_thread_num();
                this->_kernels[thread_num]->M2P(
                    &(w_item->getData()->getMultipoleData()),
                    &(w_item->getSymbolicData()),
                    leaf->getParticleContainer(),
                    &(leaf->getSymbolicData())
                    );
            }
        }
    }


    /** \brief Walk through the leaves and queue L2P tasks */
    void down_to_target() {
        for(node_t* leaf : _tree.leaves()) {
            if(leaf->getParticleContainer()->size() != 0) {
                this->create_L2P_task(leaf);
            }
        }
        // #pragma omp taskwait
    }

    /**
     * \brief Create and queue an L2P task
     *
     * \param leaf A leaf
     */
    void create_L2P_task(node_t* leaf) {
        const char* data_dep = get_dependency<dep_t::L>(leaf); (void)data_dep;
        const char* pt_dep = get_dependency<dep_t::P_t>(leaf); (void)pt_dep;
        #pragma omp task depend(inout: pt_dep[:1]) depend(in: data_dep[:1])
        {
            const int thread_num = omp_get_thread_num();
            this->_kernels[thread_num]->L2P(
                &(leaf->getData()->getLocalExpansionData()),
                &(leaf->getSymbolicData()),
                leaf->getParticleContainer());
        }
    }


    /** \brief Walk through the leaves and queue P2P tasks */
    void u_list_step() {
        for(node_t* leaf : _tree.leaves()) {
            this->create_P2P_task(leaf);
        }
        // #pragma omp taskwait
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
    void create_P2P_task(node_t* leaf) {
        using container_t = typename node_t::particle_container_t;

        container_t* const leaf_source_particle_container =
            leaf->getParticleContainer();
        container_t* const leaf_target_particle_container =
            leaf->getParticleContainer();

        // Skip empty leaves
        if( leaf_source_particle_container->size() == 0
            && leaf_target_particle_container->size() == 0) {
            return;
        }

        auto it = leaf->U.begin();
        bool do_inner = true;

        while(it != leaf->U.end()) {
            constexpr const std::size_t max_task_size = 27;
            std::size_t i = 0;
            auto first = it;

            const char* task_deps[max_task_size];
            const char* leaf_s_dep = get_dependency<dep_t::P_s>(leaf);(void)leaf_s_dep;
            const char* leaf_t_dep = get_dependency<dep_t::P_t>(leaf);(void)leaf_t_dep;

            while((it != leaf->U.end()) && (i < max_task_size)) {
                if((*it) == leaf || (*it)->getParticleContainer()->size() == 0) {
                    ++it;
                    continue;
                }

                task_deps[i] = get_dependency<dep_t::P_t>(*it);
                ++i;
                ++it;
            }

            while(i < max_task_size) {
                task_deps[i] = this->mock_dep.next();
                ++i;
            }

            #pragma omp task                                            \
                firstprivate(leaf, leaf_source_particle_container, leaf_target_particle_container, first, it, do_inner) \
                depend(inout:                                           \
                       leaf_s_dep[:1],                                  \
                       leaf_t_dep[:1],                                  \
                       task_deps[0][:1],                                \
                       task_deps[1][:1],                                \
                       task_deps[2][:1],                                \
                       task_deps[3][:1],                                \
                       task_deps[4][:1],                                \
                       task_deps[5][:1],                                \
                       task_deps[6][:1],                                \
                       task_deps[7][:1],                                \
                       task_deps[8][:1],                                \
                       task_deps[9][:1],                                \
                       task_deps[10][:1],                               \
                       task_deps[11][:1],                               \
                       task_deps[12][:1],                               \
                       task_deps[13][:1],                               \
                       task_deps[14][:1],                               \
                       task_deps[15][:1],                               \
                       task_deps[16][:1],                               \
                       task_deps[17][:1],                               \
                       task_deps[18][:1],                               \
                       task_deps[19][:1],                               \
                       task_deps[20][:1],                               \
                       task_deps[21][:1],                               \
                       task_deps[22][:1],                               \
                       task_deps[23][:1],                               \
                       task_deps[24][:1],                               \
                       task_deps[25][:1],                               \
                       task_deps[26][:1]                                \
                    )
            {
                const int thread_num = omp_get_thread_num();
                // Vectors to be filled after sort and passed to kernel P2P
                std::vector<container_t*> u_item_source_particle_containers;
                std::vector<int> u_item_indices;

                auto last = it;

                while(first != last) {
                    // Skip empty u_items
                    if((*first) == leaf || (*first)->getParticleContainer()->size() == 0) {
                        ++first;
                        continue;
                    }
                    u_item_source_particle_containers.push_back((*first)->getParticleContainer());
                    u_item_indices.push_back(compute_box_offset_index(leaf, *first, 1));
                    ++first;
                }
                // Call P2P on vectors data
                this->_kernels[thread_num]
                    ->P2P(
                        FTreeCoordinate(MortonIndex(leaf->getIndex())),
                        leaf_target_particle_container,
                        leaf_source_particle_container,
                        u_item_source_particle_containers.data(),
                        u_item_indices.data(),
                        static_cast<int>(u_item_source_particle_containers.size()),
                        do_inner
                        );
            }
            do_inner = false;

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
