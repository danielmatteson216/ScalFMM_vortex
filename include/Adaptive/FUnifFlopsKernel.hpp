#include <cstdlib>
#include <array>

#include "Utils/FConstFuncs.hpp"
#include "Containers/FTreeCoordinate.hpp"

/** \brief Kernel to evaluate the cost of an FMM run using the uniform kernel
 *
 * This kernel does no actual computation, it knows the theoretical cost of each
 * operator and sums the costs as the algorithm calls them.
 *
 * \tparam ORDER Uniform interpolation kernel order
 * \tparam NVALS Physical value count associated to each particle
 * \tparam MatrixFlopsKernelClass
 * \parblock
 * Flops counting class for the matrix kernel used.
 *
 * This class is expected to provide the same functions as those used in the
 * FUnifKernel MatrixKernelClass but that return the Flops count instead of
 * performing the calculation.
 *
 * The class is expected to be copy constructible.
 * \endparblock
 * \tparam CellClass Data storage class for a tree node
 * \tparam ContainerClass Class used to store particles
 *
 */
template<std::size_t ORDER,
         std::size_t NVALS,
         class MatrixFlopsKernelClass,
         class CellClass,
         class ContainerClass>
struct FUnifFlopsKernel {

    /// Matrix Flops kernel
    MatrixFlopsKernelClass _matrixFlopsKernel;

    /// Class representing a particle
    using particle_t = typename ContainerClass::particle_t;
    /// Number of children for a node
    constexpr static std::size_t child_count = Fpow(2, particle_t::position_t::Dim);

    /// FMM operator names
    struct operators {
        enum {P2M, M2M, M2L, P2L, L2L, L2P, P2P, M2P, ALL, size};
        std::array<std::string, size> name {{"P2M", "M2M", "M2L", "P2L", "L2L", "L2P", "P2P", "M2P", "ALL"}};
    };

    /// Flops counter array for each operator
    std::array<std::size_t, operators::size> _flops;
    /// Call counter array for each operator
    std::array<std::size_t, operators::size> _calls;

    /// Default constructor
    FUnifFlopsKernel() : FUnifFlopsKernel(MatrixFlopsKernelClass())
    {}

    /// Constructor with matrix kernel
    template<class MK>
    FUnifFlopsKernel(MK&& mk) :
        _matrixFlopsKernel(mk)
    {
        for(std::size_t i = 0; i < operators::size; ++i) {
            _flops[i] = 0;
            _calls[i] = 0;
        }
    }

    /** \brief Setup the kernel before running, no-op here
     *
     * \tparam Tree FMM tree
     *
     * \param unnamed Unused, tree to setup for
     */
    template<typename Tree>
    void setup(Tree&) {
        // TODO Add precomputation cost ?
    };

    /** \brief Particle to multipole operator
     *
     * Creates a leaf multipole expantion from its particles.
     *
     * \param leaf_data Leaf data
     * \param source_particle_container Leaf particle container
     */
    template<class SymbolicData>
    void P2M(typename CellClass::multipole_t* const,
             const SymbolicData* const,
             const ContainerClass* const source_particle_container)
    {
        std::size_t flops = source_particle_container->size() * (3 * Fpow(ORDER,3) + 3 * 3 * ORDER * (ORDER-1));
        _flops[operators::P2M] += flops;
        _flops[operators::ALL] += flops;
        _calls[operators::P2M] += 1;
    }


    /** \brief Multipole to multipole
     *
     * Propagates children multipoles to their father.
     *
     * \param node_data Parent multipole node data
     * \param child_data Array of pointer to children node data
     * \param unnamed Unused, level of the parent
     */
    template<class SymbolicData>
    void M2M(typename CellClass::multipole_t* const FRestrict,
             const SymbolicData* const,
             const typename CellClass::multipole_t*const FRestrict *const FRestrict,
             const SymbolicData* const [])
    {
        for(std::size_t idx = 0 ; idx < child_count ; ++idx){
            std::size_t flops = 3 * Fpow(ORDER,3) *(2*ORDER-1);
            _flops[operators::M2M] += flops;
            _flops[operators::ALL] += flops;
            _calls[operators::M2M] += 1;
        }
    }


    /** \brief Multipole to local development
     *
     * Computes a multipole influence on a far field local development.
     *
     * \param node_data Local development node data
     * \param v_item_data Multipole node data pointer array
     * \param unnamed Multipole nodes positions relative to node_data in terms
     * of boxes at nodes level
     * \param v_item_data_size Size of v_item_data
     * \param unnamed Unused, level of the nodes
     *
     * \note All nodes are at the same level
     */
    template<class SymbolicData>
    void M2L(typename CellClass::local_expansion_t * const FRestrict ,
             const SymbolicData* const ,
             const typename CellClass::multipole_t * const FRestrict [],
             const SymbolicData* const FRestrict [],
             const int [],
             const int v_item_data_size)
    {
        for(int idx = 0 ; idx < v_item_data_size ; ++idx){
            std::size_t flops = Fpow(ORDER, 3);
            _flops[operators::M2L] += flops;
            _flops[operators::ALL] += flops;
            _calls[operators::M2L] += 1;
        }
    }


    /** \brief Local development to local development
     *
     * Transfers a parent local development to it children.
     *
     * \param node_data Parent node data
     * \param child_data Child node data pointer array
     * \param unnamed Unused, level of the parent node
     */
    template<class SymbolicData>
    void L2L(const typename CellClass::local_expansion_t * const FRestrict,
             const SymbolicData* const,
             typename CellClass::local_expansion_t * FRestrict * const FRestrict,
             const SymbolicData* const [])
    {
        for(std::size_t idx = 0 ; idx < child_count ; ++idx){
            std::size_t flops = 3 * Fpow(ORDER,3) *(2*ORDER-1);
            _flops[operators::L2L] += flops;
            _flops[operators::ALL] += flops;
            _calls[operators::L2L] += 1;
        }
    }


    /** \brief Particle to local development
     *
     * Direct computation of a particle influence on a local development.
     *
     * \param node_data Local developmet node data
     * \param source_particle_container Particle container
     */
    template<class SymbolicData>
    void P2L(typename CellClass::local_expansion_t* const,
             const SymbolicData * const,
             const ContainerClass* const source_particle_container,
             const SymbolicData * const)
    {
        std::size_t flops = source_particle_container->size() * NVALS * _matrixFlopsKernel.evaluate();
        _flops[operators::P2L] += flops;
        _flops[operators::ALL] += flops;
        _calls[operators::P2L] += 1;
    }


    /** \brief Local development to particle
     *
     * Transfer of a local development to the underlying particles.
     *
     * \param leaf_data Leaf node data
     * \param target_particle_container Leaf particle container
     */
    template<class SymbolicData>
    void L2P(const typename CellClass::local_expansion_t* const,
             const SymbolicData* const,
             ContainerClass* const target_particle_container)
    {
        std::size_t flops = target_particle_container->size() * (4 * Fpow(ORDER,3) + 9 * ORDER * (ORDER-1));
        _flops[operators::L2P] += flops;
        _flops[operators::ALL] += flops;
        _calls[operators::L2P] += 1;
    }


    /** \brief Multipole to particle
     *
     * Direct computation of a multipole interaction with particles.
     *
     * \param node_data Multipole node data
     * \param target_particle_container Target particle container
     */
    template<class SymbolicData>
    void M2P(const typename CellClass::multipole_t* const,
             const SymbolicData* const,
             ContainerClass* const target_particle_container,
             const SymbolicData * const)
    {
        std::size_t flops = target_particle_container->size() * NVALS * _matrixFlopsKernel.evaluateBlockAndDerivative();
        _flops[operators::M2P] += flops;
        _flops[operators::ALL] += flops;
        _calls[operators::M2P] += 1;
    }


    /** \brief Particle to particle
     *
     * Direct computation of interactions between particles of a leaf and those
     * of its neighbours.
     *
     * \param unnamed Unused, coordinates of the target node
     * \param target_particle_container Target node target particle container
     * \param source_particle_container Target node source particle container
     * \param u_item_source_particle_container Source particle container pointer array
     * \param unnamed Unused, source particle container array positions relative to current leaf
     * \param u_item_count u_item_source_particle_container size
     */
    void P2P(const FTreeCoordinate& ,
             ContainerClass* const target_particle_container,
             const ContainerClass* const source_particle_container,
             ContainerClass* const /*u_item_source_particle_container*/[],
             const int /*positions*/[],
             const int /*u_item_count*/) {

        std::size_t flops = source_particle_container->size() * target_particle_container->size()
            * _matrixFlopsKernel.evaluateBlockAndDerivative();
        _flops[operators::P2P] += flops;
        _flops[operators::ALL] += flops;
        _calls[operators::P2P] += 1;

    }
};
