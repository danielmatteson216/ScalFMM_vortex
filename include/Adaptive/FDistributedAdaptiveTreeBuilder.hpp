/**
 * \brief Functions to fill a distributed adaptive tree
 * \file
 *
 * \author Quentin Khan
 */

#ifndef _FDISTRIBUTEDADAPTIVETREEBUILDER_HPP_
#define _FDISTRIBUTEDADAPTIVETREEBUILDER_HPP_

#include "FTree.hpp"
#include "inria/algorithm/distributed/mpi.hpp"
#include "inria/algorithm/distributed/sort.hpp"
#include "inria/linear_tree/balance_tree.hpp"
#include "inria/span.hpp"


/**
 * \brief Fill an adaptive tree from a distributed range
 *
 * \warning The `particles` distributed range must already match the
 * `local_linear_tree`.
 *
 * The resulting tree is a Locally Essential Tree (LET). It consists in the
 * nodes for the local particles plus the remote ghosts that are needed to
 * compute the local expansions of the local nodes.
 *
 * \param conf MPI configuration
 * \param tree Adaptive tree
 * \param local_linear_tree Linear tree
 * \param particles Distributed particle container
 * \param tree_inserter
 * \parblock Function object that inserts a particle in a tree. Must have the
 * following signature:
 *
 * ~~~{.cpp}
 * tree_inserter(Tree&, AssignableRange::value_type);
 * ~~~
 * \endparblock
 *
 * \tparam Tree Tree type
 * \tparam AssignableRange Range with an `assign` member function
 * \tparam LinearTree
 * \tparam TreeInserter
 *
 */
template<class Tree, class AssignableRange, class LinearTree, class TreeInserter>
void fill_distributed_adaptive_tree(
    inria::mpi_config conf,
    Tree& tree,
    LinearTree& local_linear_tree,
    AssignableRange& particles,
    TreeInserter&& tree_inserter
    )
{
    using inria::linear_tree::node::level;
    using inria::linear_tree::node::morton_index;

    // MPI communicator
    auto comm = conf.comm;

    // 1) Gather all processes linear trees to create a common base for the tree
    //
    // The global linear tree is the base of all processes' tree. To avoid many
    // point-to-point communications, it is gathered on the first process then
    // broadcast to the others.

    // Sizes of local linear tree
    std::vector<int> tree_sizes(comm.size());
    tree_sizes[comm.rank()] = static_cast<int>(local_linear_tree.size());
    comm.allgather(tree_sizes.data(), 1, MPI_INT);
    // Offsets for final reception buffer, the last element is the global linear
    // tree leaf count
    std::vector<int> tree_offsets(comm.size());
    std::partial_sum(begin(tree_sizes), end(tree_sizes), begin(tree_offsets));

    // note: this type is a vector
    LinearTree global_linear_tree(tree_offsets.back());

    // Gather global linear tree on root process, then broadcast it

    inria::mpi::datatype_commit_guard guard (
        inria::mpi::get_datatype<typename LinearTree::value_type>()
        );

    if(comm.rank() == 0) {
        std::copy(begin(local_linear_tree), end(local_linear_tree),
                  begin(global_linear_tree));

        std::vector<inria::mpi::request> reqs;
        reqs.reserve(comm.size() - 1);
        for(int i = 1; i < comm.size(); ++i) {
            auto req = comm.irecv(global_linear_tree.data() + tree_offsets[i-1],
                                  tree_sizes[i], guard.datatype, i,
                                  conf.base_tag + 471);
            reqs.push_back(req);
        }
        inria::mpi::request::waitall(reqs);
    } else {
        comm.send(local_linear_tree.data(),
                  static_cast<int>(local_linear_tree.size()),
                  guard.datatype, 0, conf.base_tag + 471);
    }

    comm.bcast(global_linear_tree.data(),
               static_cast<int>(global_linear_tree.size()),
               guard.datatype, 0);

    // 2) Insert ghost nodes in the tree
    //
    // Once the global linear tree is known, the tree is built to match it. All
    // the nodes created during this step are considered ghosts belonging to the
    // process that sent the corresponding leaf.
    //
    // We use the tree_offsets to track which process sent which leaf.

    int owner_rank = 0, count = 0;
    for(auto it = begin(global_linear_tree); it != end(global_linear_tree); ++it, ++count) {
        while(count >= tree_offsets[owner_rank]) {
            ++owner_rank;
        }
        tree.build_ghost_leaf(*it, owner_rank);
    }

    // 3) Insert particles in the tree
    //
    // The local particles are inserted into the tree.
    //
    // TODO: refine the linear tree (after creating the global linear tree and
    // before building ghost leaves) to avoid particle container reallocation
    // when splitting nodes.

    for(auto&& p : particles) {
        tree_inserter(tree, p);
    }

    #ifndef NDEBUG
    // Check that local ghost don't exist
    for(auto& node : tree.pre_order_walk()) {
        if(node.ghost_owner() == conf.comm.rank()) {
            std::cerr << inria::linear_tree::node::make_info(node) << ' ';
            std::cerr << node.is_leaf() << '\n';
            std::cerr << node.ghost_owner() << '\n';
        }
        assert(node.ghost_owner() != conf.comm.rank() && "Ghost should not have local process id");
    }
    #endif // NDEBUG

    // 3) For each process, find out which nodes have to be sent
    //
    // The local nodes which hold ghost nodes in their interaction lists must be
    // sent to the corresponding ghost owners. Because the interaction lists all
    // hold some symmetry [1], the remotes nodes that we will send and receive
    // are guaranteed to be usefull.
    //
    // [1] A in U_B <=> B in U_A
    //     A in V_B <=> B in V_A
    //     A in W_B <=> B in X_A

    using node_info = inria::linear_tree::node::info<Tree::Dim>;

    std::map<int, std::vector<node_info>> to_send;

    for(auto& node : tree.in_order_walk()) {
        if(! node.is_ghost()) {
            auto add_if_ghost = [&](typename Tree::node_t* ext) {
                if(ext->is_ghost()) {
                    node_info ni = inria::linear_tree::node::make_info(node);
                    int key = ext->ghost_owner();
                    if(to_send[key].empty() || to_send[key].back() != ni) {
                        to_send[key].push_back(ni);
                    }
                }
            };
            std::for_each(begin(node.U), end(node.U), add_if_ghost);
            std::for_each(begin(node.V), end(node.V), add_if_ghost);
            std::for_each(begin(node.W), end(node.W), add_if_ghost);
            std::for_each(begin(node.X), end(node.X), add_if_ghost);
        }
    }

    /////// DEBUG // remove if useless
    #if 0
    {
        using std::begin;
        using std::end;

        std::cerr << "Non ghost nodes: ";
        auto walk = tree.in_order_walk();
        std::cerr << std::accumulate(begin(walk), end(walk) , 0,
                                     [](int r, const typename Tree::node_t& n) {
                                         return r + ! n.is_ghost();
                                     }) << '\n';

        std::cerr << "Send lists: \n";
        for(auto& ni : to_send) {
            std::cerr << "  " << ni.first << ": " << ni.second.size() << '\n';
        }
    }
    #endif

    // 4) Send the nodes found in step 3.
    //
    // We know that if we send at least one node to a process, this process will
    // send us a node too.
    //
    // To setup the reception buffer for each process, a first communication
    // round is done where only the node count is sent. Then the minimal
    // information about the nodes is sent (level, morton index).

    std::map<int,int> recv_sizes;
    std::vector<inria::mpi::request> reqs;
    std::vector<int> send_sizes(to_send.size());
    int accu = 0;

    enum {SIZE_TAG, COMM_TAG};

    // Send and receive sizes

    for(auto&& proc_val : to_send) {
        int rank = proc_val.first;
        std::vector<node_info>& data = proc_val.second;
        reqs.emplace_back(
            comm.irecv(&recv_sizes[rank], 1, MPI_INT, rank, conf.base_tag + SIZE_TAG));
        send_sizes[accu] = static_cast<int>(data.size());

        /////// DEBUG
        std::cerr << "Sending to " << rank << " size = " << send_sizes[accu] << '\n';

        reqs.emplace_back(
            comm.isend(&send_sizes[accu], 1, MPI_INT, rank, conf.base_tag + SIZE_TAG));
        ++accu;
    }

    inria::mpi::request::waitall(reqs);
    reqs.clear();

    /////// DEBUG // remove if useless
    #if 0
    {
        std::cerr << "Receive sizes: \n";
        for(auto& ni : recv_sizes) {
            std::cerr << "  " << ni.first << ": " << ni.second << '\n';
        }
    }
    #endif

    // Setup reception buffer
    //
    // All the processes communication sizes are accumulated. The buffer is
    // split into sub-spans for each communication.

    std::unique_ptr<node_info[]> recv_buffer = std::unique_ptr<node_info[]>(
        new node_info[std::accumulate(begin(recv_sizes), end(recv_sizes), 0,
                                      [](const int& r, const std::pair<int,int>& e){
                                          return r + e.second;
                                      })]);

    std::map<int, inria::span<node_info>> recv_spans;

    guard = inria::mpi::datatype_commit_guard(inria::mpi::get_datatype<node_info>());

    // Start the send and reception asynchronous tasks.

    accu = 0;
    for(auto& proc_val : recv_sizes) {
        int rank = proc_val.first;
        int size = proc_val.second;
        recv_spans[rank] = inria::span<node_info>(recv_buffer.get() + accu, size);
        reqs.emplace_back(
            comm.irecv(recv_spans[rank].data(),
                       static_cast<int>(recv_spans[rank].size()),
                       guard.datatype, rank, conf.base_tag + COMM_TAG));
        accu += size;
    }

    for(auto& proc_val : to_send) {
        int rank = proc_val.first;
        node_info* buffer = proc_val.second.data();
        int size = static_cast<int>(proc_val.second.size());
        reqs.emplace_back(
            comm.isend(buffer, size, guard.datatype, rank, conf.base_tag+COMM_TAG));
    }

    inria::mpi::request::waitall(reqs);
    reqs.clear();

    // 5) Build the remote nodes as ghosts

    for(auto&& proc_val : recv_spans) {
        tree.build_ghost_subtree(proc_val.second, proc_val.first);
    }

}







#endif /* _FDISTRIBUTEDADAPTIVETREEBUILDER_HPP_ */
