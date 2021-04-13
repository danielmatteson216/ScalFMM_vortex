/**
 * This file contain function to manage the FGroupLinearTree and build a
 * GroupTree with LET
 * The LET is the Local Essential Tree
 * The LET is the symbolic information of leaf for P2P M2L and M2M operation
 *
 * @author benjamin.dufoyer@inria.fr
 */
// ==== CMAKE =====
// @FUSE_MPI
// ================
//

#ifndef _FDISTRIBUTED_GROUPTREE_BUILDER_HPP_
#define _FDISTRIBUTED_GROUPTREE_BUILDER_HPP_

#include <cmath>
#include <algorithm>
#include <stdint.h>
#include <limits.h>

#include "inria/algorithm/distributed/mpi.hpp"


// Define a MPI type for std::size_t
#if SIZE_MAX == UCHAR_MAX
   #define my_MPI_SIZE_T MPI_UNSIGNED_CHAR
#elif SIZE_MAX == USHRT_MAX
   #define my_MPI_SIZE_T MPI_UNSIGNED_SHORT
#elif SIZE_MAX == UINT_MAX
   #define my_MPI_SIZE_T MPI_UNSIGNED
#elif SIZE_MAX == ULONG_MAX
   #define my_MPI_SIZE_T MPI_UNSIGNED_LONG
#elif SIZE_MAX == ULLONG_MAX
   #define my_MPI_SIZE_T MPI_UNSIGNED_LONG_LONG
#else
   #error "FDistributedGroupTreeBuilder.hpp: MPI_TYPE what is happening here?"
#endif



namespace dstr_grp_tree_builder{

/**
 * fill_new_linear_tree this function fill the new linear tree with the value
 * of the current linear tree
 * The current linear tree is keep intact.
 * @author benjamin.dufoyer@inria.fr
 * @param  source               old linear tree
 * @param  destination          new linear tree
 * @param  nb_leaf_recev_left   number of leaf recev from left
 * @param  nb_leaf_recev_right  number of leaf recev from right
 * @param  nb_leaf_send_left    number of leaf send to left
 * @param  nb_leaf_send_right   number of leaf send to right
 */
template<class node_t>
void fill_new_linear_tree(
    std::vector<node_t>* source,
    std::vector<node_t>* destination,
    unsigned nb_leaf_recev_left,
    unsigned nb_leaf_recev_right,
    unsigned nb_leaf_send_left,
    unsigned nb_leaf_send_right
){
    unsigned min_copy = nb_leaf_send_left;
    unsigned max_copy = (unsigned)source->size() - nb_leaf_send_right;
    unsigned min_destination = nb_leaf_recev_left;
    unsigned max_destination = (unsigned)destination->size() - nb_leaf_recev_right;

    // find the smallest interval
    unsigned destination_interval = max_destination-min_destination;
    unsigned source_interval = max_copy-min_copy;
    if(source_interval < destination_interval){
        memcpy(&destination->data()[min_destination],&source->data()[min_copy],
	       sizeof(node_t)*source_interval);
    } else {
        memcpy(&destination->data()[min_destination],&source->data()[min_copy],
	       sizeof(node_t)*destination_interval);
    }
}

/**
 * This function redistribute leaf between proc according to groupSize
 * She modify the linear_tree put in parameter
 * She define a new vector avec she swap the pointer
 *
 * [RESTRICTION] : the number of leaf is very important, we cannot be in case
 *                  where a proc havn't leaf
 * the leaf on proc according to the groupsize
 * on a distributed linear tree
 * @author benjamin.dufoyer@inria.fr
 * @param  conf         MPI conf
 * @param  linear_tree  current distributed linear tree
 * @param  group_size   size of the group of leaf
 */
 template<class node_t>
 void parrallel_build_block(const inria::mpi_config& conf,
                            std::vector<node_t>* linear_tree,
                            const int& group_size)
{
    // define usefull variables
     int  nb_local_leaf  = (int)linear_tree->size();
     const int  nb_proc        = conf.comm.size();
     std::vector<int> array_global_nb_leaf(nb_proc,0);
     //int* array_global_nb_leaf  = (int *)malloc(sizeof(int) * nb_proc); //nb leaf
     const int  my_rank        = conf.comm.rank();
     // Check if i have leaf on my proc
     FAssertLF(nb_local_leaf > 0);
     // Distribute the local number of leaf to every process
     conf.comm.allgather(&nb_local_leaf,
                     1,
                     MPI_INT,
                     array_global_nb_leaf.data(),
                     1,
                     MPI_INT);

     int nb_global_leaf = 0;
     for(int i = 0 ; i < nb_proc ; i++)
         nb_global_leaf += array_global_nb_leaf[i];

     int nb_global_group = nb_global_leaf / group_size;
     int nb_local_group  = nb_global_group / nb_proc;
     int nb_leaf_needed  = nb_local_group * group_size;
     // Check if we habe enought leafs for every proc
     if( (nb_leaf_needed*(nb_proc-1)) > nb_global_leaf ){
         std::cout << " nb_leaf_needed: " << nb_leaf_needed << std::endl;
         std::cout << " nb_global_leaf: " << nb_global_leaf << std::endl;
         std::cout << " res:            " << (nb_leaf_needed*(nb_proc-1)) << std::endl;
     }
     FAssertLF( (nb_leaf_needed*(nb_proc-1)) < nb_global_leaf );  // OC: Pourquoi cela ? Ne suffit il pas de faire un exit dans le if ??

     struct message_info{
         int process_rank;
         int nb_leaf;
     };
     // We stock the future interaction in 2 vectors
     std::vector<message_info> interaction_send;
     std::vector<message_info> interaction_recev;

     // The number of leaf send and revev from left
     // it's used to fill the new linear_tree
     int nb_leaf_recev_left  = 0;
     int nb_leaf_recev_right = 0;
     int nb_leaf_send_right  = 0;
     int nb_leaf_send_left   = 0;
     //
     // COMPUTE FOR LEFT PROCESS
     // Check to know if the current proc need to send leaf
     // The compute is from left to right because it's the right process
     // who don't have a fix number of particle
     //
     // OC: Ne peut-on mettre une topologie 1d dans le communicateur pour simplifier le code
     //
     if(!my_rank == 0){ //The first process don't have computation on his left  OC: Execpt in periodic
         for(int i = 1 ; i < my_rank ; ++i ){
             array_global_nb_leaf[i] += array_global_nb_leaf[i-1];
         }
         // Check if on left process need leaf or have too many leaf
         // The idea is, if the total of leaf is a multiple of nb_leaf_needed, no communication is needed
         int nb_leaf_to_send = array_global_nb_leaf[my_rank-1] - (nb_leaf_needed * my_rank);
         if(nb_leaf_to_send != 0 ){
             message_info mess;
             mess.process_rank = my_rank-1;
             //Check if the left process have too much leaf
             if(nb_leaf_to_send > 0){
                 mess.nb_leaf = nb_leaf_to_send;
                 interaction_recev.push_back(mess);
                 //Update the array global with future value for the current proc
                 array_global_nb_leaf[my_rank] += nb_leaf_to_send;
                 nb_leaf_recev_left += nb_leaf_to_send;
             // Else the left process don't have leaf to create block_size
             } else {
                 nb_leaf_to_send = -nb_leaf_to_send;
                 mess.nb_leaf = nb_leaf_to_send;
                 interaction_send.push_back(mess);
                 //Update the array global with future value for the current proc
                 array_global_nb_leaf[my_rank] -= nb_leaf_to_send;
                 nb_leaf_send_left += nb_leaf_to_send;
             }
         }
     }

     // COMPUTE FOR RIGHT PROCESS
     // every proc compute his number of leaf after the first communication
     int nb_leaf_to_send = array_global_nb_leaf[my_rank] - nb_leaf_needed;
     // The last proc don't have proc on his right
     if( (my_rank+1) != nb_proc) {
         if(nb_leaf_to_send != 0){
             message_info mess;
             mess.process_rank = my_rank+1;
             //Check if the current process DON'T have too much leaf
             if(nb_leaf_to_send < 0){
                 nb_leaf_to_send = -nb_leaf_to_send;
                 mess.nb_leaf = nb_leaf_to_send;
                 interaction_recev.push_back(mess);
                 // Else the left process don't have leaf to form block_size
                 //Update the array global with future value for the current proc
                 array_global_nb_leaf[my_rank] += nb_leaf_to_send;
                 nb_leaf_recev_right += nb_leaf_to_send;
             } else {
                 mess.nb_leaf = nb_leaf_to_send;
                 interaction_send.push_back(mess);
                 //Update the array global with future value for the current proc
                 array_global_nb_leaf[my_rank] -= nb_leaf_to_send;
                 nb_leaf_send_right += nb_leaf_to_send;
             }
         }
     }

     // Now we have 2 vectors with all interaction with other process
     // in the first we will post every recev message
     // in a second time we post every send message

     // We declare the new linear_tree who have the new number of cell
     std::vector<node_t>* new_linear_tree = new std::vector<node_t>(array_global_nb_leaf[my_rank]);
     //Posting receiv message
     //For every interaction
     // Array of request to know the status
     inria::mpi::request tab_mpi_status[interaction_recev.size()];
     for(unsigned i = 0 ; i < interaction_recev.size(); i++ ){
         // Size of buffer
         int size_recev = (int) (sizeof(node_t)*interaction_recev[i].nb_leaf);
         // Compute the pointer to write result
         unsigned start = 0;
         if(my_rank < interaction_recev[i].process_rank){
             start = (unsigned)new_linear_tree->size() - interaction_recev[i].nb_leaf;
         }
         // Sending request
         tab_mpi_status[i] =
             conf.comm.irecv(&new_linear_tree->data()[start],
                         size_recev,
                         MPI_CHAR,
                         interaction_recev[i].process_rank,1);
     }

     ////Posting sending message
     for(unsigned i = 0 ; i < (unsigned)interaction_send.size(); ++i ){
         int sizeToSend = (int)sizeof(node_t)*interaction_send[i].nb_leaf;
         // Compute the pointer to send cell
         unsigned start = 0;
         if(my_rank < interaction_send[i].process_rank){
             start = (unsigned)linear_tree->size() - interaction_send[i].nb_leaf;
         }

         //sending leaf
         conf.comm.isend(&linear_tree->data()[start],
                        sizeToSend,
                        MPI_CHAR,
                        interaction_send[i].process_rank,1);
     }

     // Filling vector with the old linear_tree
     // The function need to know how many cell the current process send on
     // the right and on the left because the MPI request write on the same
     //  pointer
     fill_new_linear_tree(linear_tree,
                          new_linear_tree,
                          nb_leaf_recev_left,
                          nb_leaf_recev_right,
                          nb_leaf_send_left,
                          nb_leaf_send_right);

     // waiting for the send of all MPI request
     // usefull as buffer are local in the procedure
     inria::mpi::request::waitall(interaction_recev.size(),tab_mpi_status);
     //

     //free(array_global_nb_leaf);
     // swaping linear_tree pointer
     std::swap(*linear_tree,*new_linear_tree);

 }

/**
 * The idea of this function is to know the repartition of particle
 * morton index on every proc
 * Every proc send he first and his last particle morton index with a MPI
 * allgather
 * @author benjamin.dufoyer@inria.fr
 * @param   conf MPI
 * @param   particle : vector where particles a stock, they will be sorted by
 *                     morton index BEFORE calling function
 * @param   particle_repartition : a empty vector with the size of the
 *                                 number  of process
 */

 template<class type1_t,
          class type2_t>
 void share_particle_division(
     const inria::mpi_config& conf,
     std::pair<type1_t,type2_t>& my_pair,
     std::vector<std::pair<type1_t,type2_t>>& particle_index_distribution
 ){
     conf.comm.allgather(
         &my_pair,
         sizeof(my_pair),
         MPI_CHAR,
         particle_index_distribution.data(),
         sizeof(my_pair),
         MPI_CHAR);
 }

/**
 * this function is a surcharge of the previons function, she call the previous
 * function, it's just to use her more easly and callable with juste a pair
 * of index and no with juste a particle container
 * @author benjamin.dufoyer@inria.fr
 * @param  conf                     conf_MPI
 * @param  particle                 particle container
 * @param  particle_index_distribution     vector to stock the particle distribution
 */
template<class particle_t,
         class type1_t,
         class type2_t>
void share_particle_division(
    const inria::mpi_config& conf,
    std::vector<particle_t>& particle,
    std::vector<std::pair<type1_t,type2_t>>& particle_index_distribution)
{
    FAssertLF(particle_index_distribution.size() == (unsigned)conf.comm.size());
    FAssertLF(particle.size() > 0);

    std::pair<type1_t,type2_t> my_idx;
    my_idx.first = particle.front().morton_index;
    my_idx.second = particle.back().morton_index;
    // Distribute the local min max morton_index to every process
    share_particle_division(conf,my_idx,particle_index_distribution);
}

/**
 * This function sort and remove duplicate data from a vector of MortonIndex
 * @author benjamin.dufoyer@inria.fr
 * @param  data_to_modify       The MortonIndex's vector
 * @param  nb_data              number of data in the vector
 * @return                      Vector of MortonIndex
 */
std::vector<MortonIndex> sort_and_delete_duplicate_data(
    std::vector<MortonIndex> data_to_modify,
    unsigned nb_data
){
    if(nb_data != 0) {

        // Sort every morton index
        FQuickSort<MortonIndex>::QsSequential(data_to_modify.data(),nb_data);

        // Compute the number of different morton index
        // to allocate vector
        unsigned    nb_leaf = 1;
        MortonIndex last_m_idx = data_to_modify[0];
        for(unsigned i = 1 ; i < nb_data ; ++i){
            if(last_m_idx != data_to_modify[i]){
                last_m_idx = data_to_modify[i];
                ++nb_leaf;
            }
        }
        // Alloc the returned vector
        std::vector<MortonIndex> leaf_needed(nb_leaf,0);
        // Fill the returned vector
        MortonIndex current_idx = 1;
        last_m_idx = data_to_modify[0];
        leaf_needed[0] = data_to_modify[0];
        for(unsigned i = 1 ; i < nb_data ; ++i){
            if(last_m_idx != data_to_modify[i]){
                last_m_idx = data_to_modify[i];
                leaf_needed[current_idx] = data_to_modify[i];
                ++current_idx;
            }
            if((unsigned)current_idx == nb_leaf)
                break;
        }
        return leaf_needed;
    } else {
        std::vector<MortonIndex> leaf_needed(0,0);
        return leaf_needed;
    }
}



/**
 * This function compute the morton index of every leaf needed for the P2P
 * First we compute every morton index needed for every leaf
 * We sort the result
 * And to finish we remove multiple iteration
 * @author benjamin.dufoyer@inria.fr
 * @param   GroupOctreeClass Actual group tree
 * @param global_min_m_idx min morton index of the simulation
 * @param global_max_m_idx max morton index of the simulation
 * @return  std::vector with morton_index of cell who is needed to compute P2P
 */
 template<class GroupOctreeClass>
std::vector<MortonIndex> get_leaf_P2P_interaction(
    GroupOctreeClass& tree,
    const MortonIndex& global_min_m_idx,
    const MortonIndex& global_max_m_idx,
    const MortonIndex& local_min_m_idx,
    const MortonIndex& local_max_m_idx
){
    // 26 is the number of neigbors of one cell
   std::vector<MortonIndex> external_interaction(tree.getTotalNbLeaf()*26,0); //OC: Tableau tres grand
    // Reset interactions
    // idx to know where we are in the vector
    unsigned idx_vector= 0;
    // First leaf level
    {
        // We iterate on all particle group // OC: Local on the each proc ?
        for(int idxGroup = 0 ; idxGroup < tree.getNbParticleGroup() ; ++idxGroup){
            // get the particle group
            // it's ugly but, if i use template, it's not convenient
            auto* containers = tree.getParticleGroup(idxGroup);
            {
                // Iterate on every leaf
                for(int leafIdx = 0;
                    leafIdx < containers->getNumberOfLeavesInBlock();
                    ++leafIdx){
                    // Getting the morton index of the leaf
                    const MortonIndex mindex = containers->getLeafMortonIndex(leafIdx);

                    // Define array to receive computed information
                    MortonIndex interactionsIndexes[26];
                    int interactionsPosition[26];
                    //Getting coordinate of the current leaf
                    FTreeCoordinate coord(mindex);
                    // Getting neigbors
                    int counter = coord.getNeighborsIndexes(tree.getHeight(),interactionsIndexes,interactionsPosition);

                    // Iterate on every neighbors
                    for(int idxInter = 0 ; idxInter < counter ; ++idxInter){
                        // Check if the current proc already have the leaf
                        if(interactionsIndexes[idxInter] >=  local_min_m_idx &&  interactionsIndexes[idxInter] <= local_max_m_idx){
                            // do nothing
                        } else {
                            // Check if the leaf exist
                            if(interactionsIndexes[idxInter] >= global_min_m_idx && interactionsIndexes[idxInter] <= global_max_m_idx ){
                                external_interaction[idx_vector] = interactionsIndexes[idxInter];
                                ++idx_vector;
                            }
                        }
                    }
                }
            }
        }
    }
    return (sort_and_delete_duplicate_data(external_interaction,idx_vector));
}


/**
 * This function compute the leaf needed for the M2L operation
 * We take every leaf of the tree, get her parent, get the neigbors of
 * the parents and take every child of the parent's neighbors.
 * After we sort the result and remove duplicate data
 * @author benjamin.dufoyer@inria.fr
 * @param global_min_m_idx  global min morton index
 * @param global_max_m_idx  global max morton index
 * @param local_min_m_idx   local min morton index
 * @param local_max_m_idx   local max morton index
 * @param tree              local group tree
 * @param dim               optionnal parameter to compute in other dimension
 * @return vector with leaf needed for the M2M operation
 */
template<class GroupOctreeClass>
std::vector<MortonIndex> get_leaf_M2L_interaction_at_level(
        const MortonIndex& global_min_m_idx,
        const MortonIndex& global_max_m_idx,
        const MortonIndex& local_min_m_idx,
        const MortonIndex& local_max_m_idx,
        int level,
        GroupOctreeClass& tree,
        int dim = 3)
{

    // idx to fill the vector
    unsigned idx_vector = 0;
    // All External leaf
    std::vector<MortonIndex> external_interaction(tree.getNbCellGroupAtLevel(level)*tree.getNbElementsPerBlock()*189,0);
    // iterate on the group
    for(int idxGroup = 0 ; idxGroup < tree.getNbCellGroupAtLevel(level) ; ++idxGroup){
        auto* containers = tree.getCellGroup(level,idxGroup);
        MortonIndex curr_m_idx;
        // +1 is to pass the test at the first try
        for(int leafIdx = 0;
                leafIdx < containers->getNumberOfCellsInBlock();
                ++leafIdx){
            // Getting the current morton index
            curr_m_idx  = containers->getCellMortonIndex(leafIdx);
            // Compute coordinate
            MortonIndex interactionsIndexes[189];
            int interactionsPosition[189];
            FTreeCoordinate coord(curr_m_idx);
            // Getting neigbors of the father
            int counter = coord.getInteractionNeighbors(level,interactionsIndexes,interactionsPosition);
            for(int idxNeighbor = 0 ; idxNeighbor < counter ; ++idxNeighbor){
                MortonIndex tmp = interactionsIndexes[idxNeighbor];
                if( tmp    >= global_min_m_idx
                    && tmp <= global_max_m_idx)
                {
                    if(tmp >= local_min_m_idx &&
                        tmp <= local_max_m_idx){
                            // do nothing
                    } else {
                        //Stock the leaf
                        external_interaction[idx_vector] = tmp;
                        ++idx_vector;
                    }
                }
            } // end for neigbors
        } // end for leaf
    } // end for group
    // if we have leaf in the vector
    return (sort_and_delete_duplicate_data(external_interaction,idx_vector));
}


/**
 * The goal of this function is to concat the two vector
 * and erase duplicate data
 * first we compute the number of leaf
 * second we fill the vector with erase duplicate data
 * @author benjamin.dufoyer@inria.fr
 * @param  leaf_P2P vector of P2P leaf
 * @param  leaf_M2L vector of M2L leaf
 * @return vector with morton index of leaf needed
 */
std::vector<MortonIndex> concat_M2L_P2P(
    std::vector<MortonIndex>& leaf_P2P,
    std::vector<MortonIndex>& leaf_M2L)
{
    // Checking if one of vector is empty
    if(leaf_P2P.size() == 0)
        return leaf_M2L;
    if(leaf_M2L.size() == 0)
        return leaf_P2P;
    // Compute the total number of leaf needed
    unsigned idx_P2P = 0;
    unsigned idx_M2L = 0;
    std::size_t nb_leaf = 0;
    while(idx_P2P != leaf_P2P.size() && idx_M2L != leaf_M2L.size()){
        MortonIndex curr_P2P = leaf_P2P[idx_P2P];
        MortonIndex curr_M2M = leaf_M2L[idx_M2L];
        if(curr_P2P != curr_M2M){
            if(curr_P2P < curr_M2M){
                ++idx_P2P;
            } else {
                ++idx_M2L;
            }
        } else {
            ++idx_P2P;
            ++idx_M2L;
        }
        ++nb_leaf;
    }
    if(idx_P2P == leaf_P2P.size()){
        nb_leaf += (leaf_M2L.size()) - idx_M2L;
    } else {
        nb_leaf += (leaf_P2P.size()) - idx_P2P;
    }
    // Allocate the vector
    std::vector<MortonIndex> leaf_needed(nb_leaf,-1);
    idx_P2P = 0;
    idx_M2L = 0;
    std::size_t idx_leaf = 0;
    // fill the vector
    while(idx_P2P != leaf_P2P.size() && idx_M2L != leaf_M2L.size()){
        MortonIndex curr_P2P = leaf_P2P[idx_P2P];
        MortonIndex curr_M2M = leaf_M2L[idx_M2L];
        if(curr_P2P != curr_M2M){
            if(curr_P2P < curr_M2M){
                leaf_needed[idx_leaf] = leaf_P2P[idx_P2P];
                ++idx_P2P;
            } else {
                leaf_needed[idx_leaf] = leaf_M2L[idx_M2L];
                ++idx_M2L;
            }
            //++nb_leaf;
        } else {
            leaf_needed.at(idx_leaf) = leaf_M2L[idx_M2L];
            ++idx_P2P;
            ++idx_M2L;
        }
        ++idx_leaf;
    }
  //   std::cout << "  idx_leaf  " << idx_leaf << " nb_leaf "  << nb_leaf <<std::endl;
    // Copy the rest of leaf with a memcpy
    if(idx_leaf < nb_leaf){
  //       std::cout << "  MEMCOPY " << std::endl;
        void* destination =  &leaf_needed.data()[idx_leaf];
        void* source;
        std::size_t num = 0;
        if(idx_P2P == leaf_P2P.size()){
 //           std::cout << "    M2L " <<idx_M2L << "  " << leaf_M2L[idx_M2L]<< " "
  //                    << leaf_M2L.size() -1 << "  " << leaf_M2L[leaf_M2L.size() -1 ]
 //                        << "  "<<  leaf_M2L.size() << " " << leaf_M2L.size() - idx_M2L <<std::endl;
            source = &leaf_M2L[idx_M2L];
            num = sizeof(MortonIndex)*(leaf_M2L.size() - idx_M2L);
        } else {
            source = &leaf_P2P[idx_P2P];
            num = sizeof(MortonIndex)* (leaf_P2P.size() - idx_P2P);
        }
        memcpy(destination,source,num);
    }
    return leaf_needed;

}



/**
 * This function compute the number of leaf needed on every proc,
 * she check if leaf exit, if a leaf doesn't exist this function
 * remove here at the end of the function
 * @author benjamin.dufoyer@inria.fr
 * @param  needed_leaf vector with the morton index of every leaf
 * needed for my proc
 * @param  particle_repartition vector with the repartation of particle on
 * every proc
 * @param  conf conf MPI
 * @return vector<std::size_t> with a size nb_proc*nb_proc
 * it's a interaction matrix
 */
std::vector<std::vector<std::size_t>> get_matrix_interaction(
    std::vector<MortonIndex>& needed_leaf,
    std::vector<std::pair<MortonIndex,MortonIndex>>& particle_distribution,
    const inria::mpi_config& conf)
{
    // Getting MPI Info
    const int  nb_proc        = conf.comm.size();
    const int  my_rank        = conf.comm.rank();
    // Alloc interaction matrix
    std::vector<std::vector<std::size_t>> matrix_interaction(2,std::vector<std::size_t>(nb_proc,0));
    std::vector<std::size_t> global_matrix_interaction(nb_proc,0);
    // Initialise idx on particle_distribution
    size_t idx_part = 0;
    // Iterate on every leaf to know where she is
    MortonIndex max_morton_index = 0;
    if(needed_leaf.size() > 0)
        max_morton_index = needed_leaf[needed_leaf.size()-1]+1;
    // iterate on every mortonIndex
    for(unsigned idx_leaf = 0; idx_leaf < needed_leaf.size(); ++idx_leaf){
        MortonIndex current_leaf = needed_leaf[idx_leaf];
        // if she is on the current proc
        if(current_leaf >= particle_distribution[idx_part].first
        && current_leaf <= particle_distribution[idx_part].second){
            if(idx_part == (unsigned)my_rank){
                needed_leaf[idx_leaf] = max_morton_index;
            } else {
                matrix_interaction[0][idx_part] += 1;
            }
        } else {
            // While the current leaf is not on the good interval
            while(idx_part < particle_distribution.size() && particle_distribution[idx_part].second < current_leaf){
                idx_part += 1;
            }
            if(idx_part == particle_distribution.size())
                break;
            if(particle_distribution[idx_part].first > current_leaf){
                // in this case the leaf is not in interval, so she doesn't exist
                needed_leaf[idx_leaf] = max_morton_index;
            } else {
                // In the case it's a normal case, we juste increment the
                // number of leaf send at the proc idx_part
                if(idx_part == (unsigned)my_rank){
                    needed_leaf[idx_leaf] = max_morton_index;
                } else {
                    matrix_interaction[0][idx_part] += 1;
                }
            }
        }
    }
    // i don't need to send to me
    matrix_interaction[0][my_rank] = 0;
    // now we have the number of leaf to send at every proc
    // we proceed a AlltoAll to share this information at every proc
    conf.comm.alltoall(matrix_interaction[0].data(),
                       1,
                       my_MPI_SIZE_T,
                       matrix_interaction[1].data(),
                       1,
                       my_MPI_SIZE_T);
    // removing bad leaf
    needed_leaf.erase(std::remove(needed_leaf.begin(),needed_leaf.end(),max_morton_index),needed_leaf.end());

    return {begin(matrix_interaction),end(matrix_interaction)};
}

/**
 * This function return the number of block at node level
 * This algo is different than the computation at leaf level, because
 * it's only the proc who have the smallest rank who have the attribution of
 * the block
 * @author benjamin.dufoyer@inria.fr
 * @param  tree                 local GroupTree
 * @param  node_needed          List of needed node
 * @param  nb_node              Number of node needed in the array
 * @param  level                Level of the node
 * @return                      Vector of index of block
 */
 template<class GroupOctreeClass>
 std::vector<MortonIndex> get_nb_block_from_node(GroupOctreeClass& tree,
              MortonIndex* node_needed,
              std::size_t nb_node,
              int level)
 {
     int idx_vect = 0 ;
     std::vector<int> block_to_send(tree.getNbCellGroupAtLevel(level),0);
     unsigned idx_node = 0;
     // iterate on every group
     for(unsigned idx_group = 0; idx_group < (unsigned)tree.getNbCellGroupAtLevel(level) ;++idx_group){
         // if the current block hasnt been already send
        auto* containers = tree.getCellGroup(level,idx_group);
        // check if we have check every node
        if(idx_node == nb_node){
            break;
        }
        // while the morton index of the current node is not high
        while(idx_node < nb_node && node_needed[idx_node] < containers->getStartingIndex()){
            ++idx_node;
        }
        // while the current morton index is in the block
        while(idx_node < nb_node && node_needed[idx_node] < containers->getEndingIndex()){
             // if the container have the current morton index
             // keep the block and go out of the while
             if(containers->isInside(node_needed[idx_node])){
                 block_to_send[idx_vect] = idx_group;
                 ++idx_vect;
                 ++idx_node;
                 break;
             }
             ++idx_node;
        }
        if(idx_node == nb_node){
            break;
        }
     }
     return {block_to_send.begin(),block_to_send.begin()+idx_vect};
 }

 /*
  * In this function know the block needed by every proc
  * So we need need to compute the parents of every block and send the
  * number of parent's block to every proc so we have
  *
  * 1st Step : Compute the parents block
  * 2nd Step : Send and recev the number of parents block per level
  * @author benjamin.dufoyer@inria.fr
  * @param  nb_block_to_receiv number of block to receiv by every proc per level
  * @param  list_of_block_to_send list of block to send for every proc per level
  * @param  tree         local group octree
 * @param  conf          MPI CONF
  */
 template<class GroupOctreeClass>
 void send_get_number_of_block_node_level(
     std::vector<MortonIndex>& vect_recv,
     std::vector<std::vector<std::size_t>> global_matrix_interaction,
     std::size_t& nb_msg_recv,
     GroupOctreeClass& tree,
     std::vector<std::pair<int,int>>& nb_block_to_receiv,
     std::vector<std::pair<int,std::vector<MortonIndex>>>& block_to_send,
     int level,
     const inria::mpi_config& conf
 )
 {
     int idx_status = 0;
     int idx_proc   = 0;
     inria::mpi::request tab_mpi_status[nb_msg_recv];

     // Post the number reception of the number of block
     for(unsigned i = 0; i < global_matrix_interaction[0].size() ; ++i)
     {
         // If we have interaction with this proc
         if(global_matrix_interaction[0][i] != 0){
             nb_block_to_receiv[idx_status].first = idx_proc;
             tab_mpi_status[idx_status] = conf.comm.irecv(
                 &nb_block_to_receiv[idx_status].second,
                     1,
                     MPI_INT,
                     idx_proc,1
                 );
                 idx_status += 1;
         }
         idx_proc += 1;
     }
     idx_proc = 0;
     std::size_t idx_vect = 0;
     idx_status = 0;
     // Posting every send message
     for(unsigned i = 0; i < global_matrix_interaction[1].size() ; ++i){
         // If we have interaction with this proc
         if(global_matrix_interaction[1][i] != 0){
             // Compute the number of leaf
             int  nb_block;
             block_to_send[idx_status].second = get_nb_block_from_node(
                                 tree,
                                 &vect_recv.data()[idx_vect],
                                 global_matrix_interaction[1][i],
                                 level);
             nb_block = (int)block_to_send[idx_status].second.size();
             block_to_send[idx_status].first = i;
             // send the number of leaf
             conf.comm.isend(
                 &nb_block,
                 1,
                 MPI_INT,
                 i,1
             );
             idx_vect += global_matrix_interaction[1][i];
             idx_status += 1;
         }
         idx_proc += 1;
     }
     // Waiting for all request
     if(nb_msg_recv != 0 ){
         inria::mpi::request::waitall(nb_msg_recv,tab_mpi_status);
     }
 }


/**
 * This function compute the leaf needed by every proc
 * @author benjamin.dufoyer@inria.fr
 * @param  needed_leaf      List of leaf needed by me
 * @param  global_matrix_interaction
 * @param  nb_msg_recv      number of recev message
 * @param  nb_leaf_recv     number of recev leaf
 * @param  conf             MPI conf
 * @return                  Vector with all leaf needed by other proc
 */
std::vector<MortonIndex> send_get_leaf_morton(
    std::vector<MortonIndex>&     needed_leaf,
    std::vector<std::vector<std::size_t>>&     global_matrix_interaction,
    std::size_t&                  nb_msg_recv,
    std::size_t&                  nb_leaf_recv,
    const inria::mpi_config& conf)
{
    // allocate tab of mpi status for synchronisazion
    std::vector<inria::mpi::request> tab_mpi_status(nb_msg_recv);

    // Initialiser the reception vector
    std::vector<std::size_t> vect_recv(nb_leaf_recv,0);

    // Initialize every variable
    int    idx_status = 0;
    unsigned idx_vect   = 0;
    int    idx_proc   = 0;
    // Posting every recv message
    for(unsigned i = 0; i < global_matrix_interaction[1].size() ; ++i ){
        if(global_matrix_interaction[1][i] != 0){
            std::size_t nb_leaf = global_matrix_interaction[1][i];
            tab_mpi_status[idx_status] = conf.comm.irecv(
                &vect_recv[idx_vect],
                int(nb_leaf*sizeof(MortonIndex)),
                MPI_CHAR,
                i,1
            );
            idx_vect += (unsigned)nb_leaf;
            idx_status+= 1;
        }
        idx_proc += 1;
    }

    // Posting every send message
    idx_proc = 0;
    idx_vect = 0;
    for(unsigned i = 0; i < global_matrix_interaction[0].size() ; ++i){
        if(global_matrix_interaction[0][i] != 0){
            std::size_t nb_leaf = global_matrix_interaction[0][i];
            conf.comm.isend(
                &needed_leaf[idx_vect],
                int(nb_leaf*sizeof(MortonIndex)),
                MPI_CHAR,
                i,1
            );
            idx_vect += (unsigned)nb_leaf;
        }
        idx_proc += 1;
    }
    if(nb_msg_recv != 0 ){
        inria::mpi::request::waitall(idx_status,tab_mpi_status.data());
    }
    conf.comm.barrier();
    return{begin(vect_recv),end(vect_recv)};

}


struct particle_symbolic_block{
    int idx_global_block;
    FSize nb_particles;
    std::vector<FSize> nb_particle_per_leaf;
    friend
    std::ostream& operator<<(std::ostream& os, const particle_symbolic_block& n) {
        return os << "--->  nb particle " << n.nb_particles << "<--";
    }
};

struct cell_symbolic_block{
    int idx_global_block;
    MortonIndex start_index;
    MortonIndex end_index;
    int nb_leaf_in_block;
    std::vector<MortonIndex> m_idx_in_block;

    friend
    std::ostream& operator<<(std::ostream& os, const cell_symbolic_block& n) {
        return os << "--> n_block : " << n.idx_global_block << " start : " << n.start_index << " end : " << n.end_index << " nb_leaf " << n.nb_leaf_in_block  << "<--";
    }
};

/**
 * This function echange symbolic information of block of a GroupTree
 *
 * @author benjamin.dufoyer@inria.fr
 * @param  nb_block_to_receiv It's a vector of pair of int,int the first int
 *                            is for the number of the proc who send me block
 *                            the second int is the number of block he send to
 *                            me
 * @param  block_to_send      it's a vector of pair, the first is the number of
 *                            the proc who i send the block, the second is the
 *                            list of the block needed by the proc
 * @param  nb_msg_recv        it's the number of message i will receiv
 * @param  tree              it's the GroupTree where block are stock
 * @param  conf                it's the MPI conf
 */

template<class GroupOctreeClass>
std::pair<std::vector<cell_symbolic_block>,std::vector<particle_symbolic_block>> exchange_block(
    std::vector<std::pair<int,int>>& nb_block_to_receiv,
    std::vector<std::pair<int,std::vector<MortonIndex>>>& block_to_send,
    GroupOctreeClass& tree,
    int level,
    const inria::mpi_config& conf
)
{
    struct sending_cell_structure{
        int idx_global_block;
        MortonIndex start_index;
        MortonIndex end_index;
        int nb_leaf_in_block;
    };

    int my_rank =  conf.comm.rank();
    bool leaf_level = ( level == tree.getHeight() -1 );
    int block_size = tree.getNbElementsPerBlock();
    // declaration of the array of MPI status for synchro
    unsigned nb_message_recv  = 0;
    unsigned nb_block_to_recv = 0;
    for(unsigned i = 0 ; i <  nb_block_to_receiv.size() ;++i ){
        if(nb_block_to_receiv[i].second != 0 && nb_block_to_receiv[i].first != my_rank){
            // computing of the number of message and the number of block to recv
            ++nb_message_recv;
            nb_block_to_recv += nb_block_to_receiv[i].second;
        }
    }
    if(leaf_level)
        nb_message_recv = nb_message_recv+ (nb_message_recv*2);
    // compute the total number of message
    std::vector<inria::mpi::request> tab_mpi_status(nb_message_recv);

    // Compute the total number of block i will send
    std::size_t total_size = 0;
    for(unsigned i = 0 ; i < block_to_send.size(); ++i){
        total_size += block_to_send[i].second.size();
    }

    // Buffer to send the cell structure
    std::vector<sending_cell_structure>             cell_symb_to_send(total_size);

    // buffer to send the morton index
    std::vector<size_t>                morton_index_send(total_size*block_size,0);
    // buffer to send particles block
    std::vector<FSize>  nb_particle_per_leaf(0,0);
    std::vector<unsigned> particle_symb_to_send(0,0);

    if(leaf_level){
        nb_particle_per_leaf.resize(total_size*block_size,0);
        particle_symb_to_send.resize(total_size);
    }


    std::size_t idx_vect_to_send = 0;
    std::size_t idx_m_idx        = 0;
    // Filling the buffer of block
    for(unsigned i = 0 ; i < block_to_send.size(); ++i){
        for(unsigned j = 0 ; j < block_to_send[i].second.size() ; ++j){
            if(block_to_send[i].first != my_rank){
                auto* container = tree.getCellGroup(level ,((int)block_to_send[i].second[j]));
                sending_cell_structure block_to_add{
                    container->getIdxGlobal(),
                    container->getStartingIndex(),
                    container->getEndingIndex(),
                    container->getNumberOfCellsInBlock()
                };
                // Get all morton index of the block
                for(int k = 0 ; k < container->getNumberOfCellsInBlock(); ++k ){
                    morton_index_send[idx_m_idx+k] = container->getCellMortonIndex(k);
                }
                for(int k = container->getNumberOfCellsInBlock(); k  < block_size; ++k ){
                    morton_index_send[idx_m_idx+k] = container->getCellMortonIndex(container->getNumberOfCellsInBlock()-1);
                }
                // add the block to the vector
                cell_symb_to_send[idx_vect_to_send] = block_to_add;
                if(leaf_level){
                    // get the particle container associated
                    auto* container_particle =  tree.getParticleGroup(((int)block_to_send[i].second[j]));
                    particle_symb_to_send[idx_vect_to_send] = container_particle->getIdxGlobal();
                    // iterate on every leaf
                    for(int k = 0 ; k < container_particle->getNumberOfLeavesInBlock(); ++k){
                        // stock the number of particles in the leaf
                        nb_particle_per_leaf[idx_m_idx+k] = container_particle->getNbParticlesInLeaf(k);
                    }
                }
                idx_m_idx += block_size;
                ++idx_vect_to_send;
            }
        }
    }
    // Now i have my vector(s) to send all of my blocks

    // the first vector will contain all of cell_block and the send all of
    // particle block
    std::vector<sending_cell_structure> symbolic_block_rcv(nb_block_to_recv);

    int size_of_vect = nb_block_to_recv*block_size;
    std::vector<FSize>       nb_part_leaf(0,0);
    std::vector<unsigned>    idx_global_particle_rcv(0,0);
    if(leaf_level){
        nb_part_leaf.resize(size_of_vect,0);
        idx_global_particle_rcv.resize(nb_block_to_recv,0);
    }


    int idx_status = 0;
    unsigned offset_block = 0;
    unsigned offset_m_idx = 0;
    // Posting recv
    for(unsigned i = 0; i < nb_block_to_receiv.size(); ++i)
    {
        if(nb_message_recv == 0)
            break;
        // Setting parameter
        int source   = nb_block_to_receiv[i].first;
        int nb_block = nb_block_to_receiv[i].second;
        if(nb_block != 0 && source != my_rank){
            // Posting reveiv message
            tab_mpi_status[idx_status] =
            conf.comm.irecv(
                &symbolic_block_rcv[offset_block],
                int(nb_block*sizeof(sending_cell_structure)),
                MPI_CHAR,
                source,1
            );
            idx_status += 1;

            // if it's the leaf level, i need to recv the particle block
            if(leaf_level){
                tab_mpi_status[idx_status] =
                conf.comm.irecv(
                    &nb_part_leaf[offset_m_idx],
                    int((nb_block*block_size*sizeof(FSize))),
                    MPI_CHAR,
                    source,3
                );
                idx_status += 1;

                tab_mpi_status[idx_status] =
                conf.comm.irecv(
                    &idx_global_particle_rcv[offset_block],
                    nb_block,
                    MPI_UNSIGNED,
                    source,4
                );
                idx_status += 1;

            }
            offset_block += nb_block;
            offset_m_idx += (nb_block*block_size);
        }
    }
    FAssertLF(idx_status == (int)nb_message_recv);

    // post sending message
    offset_block = 0;
    offset_m_idx = 0;
    for(unsigned i = 0 ; i < block_to_send.size(); ++i){
        // Posting send message
        int nb_block    = (int)block_to_send[i].second.size();
        int destination = block_to_send[i].first;

        if(nb_block != 0 && destination != my_rank){
            // Setting parameters
            conf.comm.isend(
                &cell_symb_to_send[offset_block],
                int(nb_block*sizeof(sending_cell_structure)),
                MPI_CHAR,
                destination,1
            );

            if(leaf_level){
                conf.comm.isend(
                    &nb_particle_per_leaf[offset_m_idx],
                    int(nb_block*block_size*sizeof(FSize)),
                    MPI_CHAR,
                    destination,3
                );
                conf.comm.isend(
                    &particle_symb_to_send[offset_block],
                    nb_block,
                    MPI_UNSIGNED,
                    destination,4
                );
            }
            offset_block = (offset_block+nb_block);
            offset_m_idx = (offset_m_idx + (nb_block*block_size));
        }

    }
    // Waiting for all request
    inria::mpi::request::waitall(idx_status,tab_mpi_status.data());
    // Sending morton idx
    if(leaf_level)
        nb_message_recv /= 2;
    idx_status = 0;
    std::vector<size_t> m_idx_to_recv(size_of_vect,0);
    inria::mpi::request tab_status[nb_message_recv];
    offset_block = 0;
    offset_m_idx = 0;
    for(unsigned i = 0; i < nb_block_to_receiv.size(); ++i)
    {
        if(nb_message_recv == 0)
            break;
        // Setting parameter
        int source   = nb_block_to_receiv[i].first;
        int nb_block = nb_block_to_receiv[i].second;
        if(nb_block != 0 && source != my_rank){
            // Posting reveiv message
            tab_status[idx_status] =
            conf.comm.irecv(
                &m_idx_to_recv.data()[offset_m_idx],
                int(nb_block*block_size),
                my_MPI_SIZE_T,
                source,2
            );
            idx_status += 1;
            offset_m_idx = (offset_m_idx + (nb_block*block_size));
        }
    }

    offset_block = 0;
    offset_m_idx = 0;
    for(unsigned i = 0 ; i < block_to_send.size(); ++i){
        // Posting send message
        int nb_block    = (int)block_to_send[i].second.size();
        int destination = block_to_send[i].first;

        if(nb_block != 0 && destination != my_rank){
            conf.comm.isend(
                &morton_index_send.data()[offset_m_idx],
                int(nb_block*block_size),
                my_MPI_SIZE_T,
                destination,2
            );
            offset_m_idx = (offset_m_idx + (nb_block*block_size));
        }
    }
    if(nb_message_recv > 0)
        inria::mpi::request::waitall(idx_status,tab_status);
    conf.comm.barrier();

    if(nb_message_recv > 0){
        std::pair<std::vector<cell_symbolic_block>,
                std::vector<particle_symbolic_block>> pair_return;
        pair_return.first.resize(symbolic_block_rcv.size());

        if(leaf_level)
            pair_return.second.resize(symbolic_block_rcv.size());
        else
            pair_return.second.resize(0);

        int nb_leaf_before_me = 0;

        for(unsigned i = 0 ; i < symbolic_block_rcv.size() ; ++i){
            // filling symbolique information
            cell_symbolic_block new_block{
                symbolic_block_rcv[i].idx_global_block,
                symbolic_block_rcv[i].start_index,
                symbolic_block_rcv[i].end_index,
                symbolic_block_rcv[i].nb_leaf_in_block
            };
            // filling morton index vector
            new_block.m_idx_in_block.clear();
            new_block.m_idx_in_block.insert(
                new_block.m_idx_in_block.begin(),
                m_idx_to_recv.begin()+nb_leaf_before_me,
        m_idx_to_recv.begin()+(nb_leaf_before_me+new_block.nb_leaf_in_block));
                for (size_t nb = 0; nb < new_block.m_idx_in_block.size(); nb++) {
            if(new_block.m_idx_in_block[nb] > symbolic_block_rcv[i].end_index || new_block.m_idx_in_block[nb] < symbolic_block_rcv[i].start_index){
                std::cout << "ERROR" << i << '\n';
                std::cout << new_block.m_idx_in_block[nb] << " " ;
                std::cout << symbolic_block_rcv[i].end_index << " ";
                std::cout << symbolic_block_rcv[i].start_index << '\n';
                std::cout << "nb_leaf_in_block "<< symbolic_block_rcv[i].nb_leaf_in_block << "\n";
                std::cout << m_idx_to_recv.size() << std::endl;
                for(int idx = nb_leaf_before_me ; idx < nb_leaf_before_me+block_size; ++idx ){
                    std::cout << " " << m_idx_to_recv[i] ;
                }
                std::cout << std::endl;
            }
        }
        // adding to the vector
        pair_return.first[i] = new_block;
            if(leaf_level){
                particle_symbolic_block new_p_block;
                new_p_block.idx_global_block =idx_global_particle_rcv[i];
                new_p_block.nb_particle_per_leaf.clear();
                new_p_block.nb_particle_per_leaf.insert(
                new_p_block.nb_particle_per_leaf.begin(),
                nb_part_leaf.begin()+nb_leaf_before_me,
                nb_part_leaf.begin()+(nb_leaf_before_me+new_block.nb_leaf_in_block));

                FSize nb_particles_in_block = 0;

                for(unsigned j = 0 ; j <  new_p_block.nb_particle_per_leaf.size() ; ++j){
                    nb_particles_in_block += new_p_block.nb_particle_per_leaf[j];
                }

                new_p_block.nb_particles = nb_particles_in_block;
                pair_return.second[i] = new_p_block;
            }
            nb_leaf_before_me += block_size;
        }
        return {pair_return.first,pair_return.second};
    }
    std::pair<std::vector<cell_symbolic_block>,
            std::vector<particle_symbolic_block>> pair_return;
    pair_return.first.resize(0);
    pair_return.second.resize(0);
    return {pair_return.first,pair_return.second};
}



/**
 * Compute the block needed by other proc, this function stock the idx of
 * block into the current_level variable
 * This function have 2 step, first we compute the number of parents of
 * the block needed a the level -1
 * after we allocate the vector and fill me
 * @author benjamin.dufoyer@inria.fr
 * @param  under_level      block for the under level
 * @param  current_level    block for the current level
 * @param  level            level of the tree
 * @param  tree             Local GroupTree
 */
template<class GroupOctreeClass>
void compute_block_node_level(
    std::vector<std::pair<int,std::vector<std::size_t>>>& under_level,
    std::vector<std::pair<int,std::vector<std::size_t>>>& current_level,
    int level,
    GroupOctreeClass& tree
){
    FAssertLF(under_level.size() == current_level.size() );
    // Iterate on every interaction of the under level
    for(unsigned i = 0 ; i < under_level.size() ; ++i){
        // Init variables for the search
        std::size_t start_m_idx = tree.getCellGroup(level+1,int(tree.getNbCellGroupAtLevel(level+1)-1))->getEndingIndex() ;
        std::size_t last_m_idx  = 0;
        std::size_t current_min_m_idx;
        std::size_t current_max_m_idx;
        int nb_block = 0;
        int last_block = -1;
        // Iterate on every block to find his parents
        for(unsigned j = 0 ; j < under_level[i].second.size() ; j++){
            // get the morton index at level
            current_min_m_idx = (tree.getCellGroup(level+1,
                                    (int)under_level[i].second[j])->getStartingIndex() >> 3);
            current_max_m_idx = (tree.getCellGroup(level+1,
                                    (int)under_level[i].second[j])->getEndingIndex() >> 3);
            // If The parents of this block is unknow, it's a new block to add
            if(current_min_m_idx < start_m_idx || current_max_m_idx > last_m_idx){
                // iterate on all group of the current level to find the
                // good block
                for(int idx_block = 0 ;  idx_block < tree.getNbCellGroupAtLevel(level) ; ++idx_block){
                    // get the current block
                    auto* container = tree.getCellGroup(level,idx_block);
                    // if the morton index is inside
                    if(container->isInside(current_min_m_idx) ||
                       container->isInside(current_max_m_idx) ){
                        // update the number of block
                        // if it's not the last block
                        if(idx_block != last_block){
                            nb_block += 1;
                            last_block = idx_block;
                            // update idx
                            if((std::size_t)container->getStartingIndex() < start_m_idx){
                                start_m_idx = container->getStartingIndex();
                            }
                            if((std::size_t)container->getEndingIndex() > last_m_idx){
                                last_m_idx = container->getEndingIndex();
                            }
                        }
                        // stop seeking for this idx
                        break;
                    }
                }
            }
        }
        // Modify current_level
        current_level[i].first = under_level[i].first;
        current_level[i].second.resize(nb_block);
        // reinit variables
        unsigned idx_vect = 0;
        start_m_idx = tree.getCellGroup(level+1,int(tree.getNbCellGroupAtLevel(level+1)-1))->getEndingIndex();
        last_m_idx  = 0;
        last_block = -1;
        // fill the current level with the same technic
        for(unsigned j = 0 ; j < under_level[i].second.size() ; j++){
            // get the morton index at level
            current_min_m_idx = (tree.getCellGroup(level+1,
                                (int)under_level[i].second[j])->getStartingIndex() >> 3);
            current_max_m_idx = (tree.getCellGroup(level+1,
                                (int)under_level[i].second[j])->getEndingIndex() >> 3);
            // The parents of this block is unknow
            if(current_min_m_idx < start_m_idx || current_max_m_idx > last_m_idx){
                // iterate on all group of the current level
                for(int idx_block = 0 ;  idx_block < tree.getNbCellGroupAtLevel(level) ; ++idx_block){
                    // get the current block
                    auto* container = tree.getCellGroup(level,idx_block);
                    if(container->isInside(current_min_m_idx)  ||
                       container->isInside(current_max_m_idx)){
                        if(idx_block != last_block){
                            // update the number of block
                            current_level[i].second[idx_vect] = idx_block;
                            idx_vect += 1;
                            last_block = idx_block;
                            // update idx
                            if((std::size_t)container->getStartingIndex() < start_m_idx){
                                start_m_idx = container->getStartingIndex();
                            }
                            if((std::size_t)container->getEndingIndex() > last_m_idx){
                                last_m_idx = container->getEndingIndex();
                            }
                        }
                        break;
                    }
                }
            }
        }
    }
}



/**
 * this function return a vector with the symbolic information of the needed
 * block to build a LET Group tree
 *
 * @author benjamin.dufoyer@inria.fr
 * @param  needed_leaf list of the leaf needed
 * @param  global_matrix_interaction matrix of interaction
 * @param  tree local group tree
 * @param  conf MPI conf
 */
template<class GroupOctreeClass>
std::pair<std::vector<cell_symbolic_block>,std::vector<particle_symbolic_block>> send_get_symbolic_block_at_level(
    std::vector<MortonIndex>&               needed_leaf,
    std::vector<std::vector<size_t>>&       matrix_interaction,
    GroupOctreeClass&                       tree,
    int                                     level,
    const inria::mpi_config&                conf
){
    std::size_t nb_msg_recv = 0;
    std::size_t nb_leaf_recv = 0;
    std::size_t nb_msg_send = 0;

    // Getting the number of sended message
    for(unsigned i = 0 ; i < matrix_interaction[0].size() ; ++i){
        if(matrix_interaction[0][i] > 0 )
            ++nb_msg_send;
    }

    // Getting the number of recv message and the number of leaf
    for(unsigned i = 0 ; i < matrix_interaction[1].size() ; ++i){
        if(matrix_interaction[1][i] > 0){
            ++nb_msg_recv;
            nb_leaf_recv += matrix_interaction[1][i];
        }
    }


    ////////////////////////////////////////////////////////////
    /// FIRST STEP
    /// Getting the list of leaf needed by every proc
    ////////////////////////////////////////////////////////////
    std::vector<MortonIndex> vect_recv =
        send_get_leaf_morton(
            needed_leaf,
            matrix_interaction,
            nb_msg_recv,
            nb_leaf_recv,
            conf);
    // free needed_leaf
    std::vector<MortonIndex>().swap(needed_leaf);

    ////////////////////////////////////////////////////////////
    // SECOND STEP
    // Compute the block to send to other proc
    // And send the number of block sended
    ////////////////////////////////////////////////////////////
    // Init variable to stock
    std::vector<std::pair<int,int>> nb_block_to_receiv(nb_msg_send);
    std::vector<std::pair<int,std::vector<MortonIndex>>>
                                list_of_block_to_send(nb_msg_recv);

    send_get_number_of_block_node_level(
        vect_recv,
        matrix_interaction,
        nb_msg_send,
        tree,
        nb_block_to_receiv,
        list_of_block_to_send,
        level,
        conf);

    std::vector<MortonIndex>().swap(vect_recv);




    ////////////////////////////////////////////////////////////
    /// THIRD STEP
    /// Getting the list of leaf needed by every proc
    ////////////////////////////////////////////////////////////
    return exchange_block(
              nb_block_to_receiv,
              list_of_block_to_send,
              tree,
              level,
              conf);
}


/**
 * This algorithm compute the global index of every block in the local tree
 * @author benjamin.dufoyer@inria.fr
 * @param  tree [description]
 * @param  conf [description]
 */
 template<class GroupOctreeClass>
 int set_cell_group_global_index_at(
           GroupOctreeClass&     tree,
           int                   level,
           int                   nb_block_under_level,
     const inria::mpi_config&    conf,
           bool                  particle = false
 ){
     int nb_proc = conf.comm.size();
     int my_rank = conf.comm.rank();
     int nb_block_before_me = 0;
     int my_nb_block;
     if(!particle) {
        my_nb_block = tree.getNbCellGroupAtLevel(level);
    } else {
        my_nb_block = tree.getNbParticleGroup();
     }
     // get the number of block at the under level
     if(my_rank == 0){
         conf.comm.recv(
             &nb_block_before_me,
             1,
             MPI_INT,
             nb_proc-1,level);
     } else if( my_rank == nb_proc-1){
         conf.comm.send(
             &nb_block_under_level,
             1,
             MPI_INT,
             0,level
         );
     }

     if(nb_proc != 0){
         // get the number of block before me
         if(my_rank != 0){
             conf.comm.recv(
                 &nb_block_before_me,
                 1,
                 MPI_INT,
                 my_rank-1,level
             );
         }
         // send the number of block before me with my number
         if(my_rank != (nb_proc-1) ){
             int nb_block_after_me = my_nb_block + nb_block_before_me;
             conf.comm.send(
                 &nb_block_after_me,
                 1,
                 MPI_INT,
                 my_rank+1,level
             );
         }
     }
     // Now i have the total number of block before me, i will compute
     // the idex of all of my block at this level
     for(int idx_group = 0 ; idx_group < my_nb_block ;++idx_group){
         if(!particle){
            auto* container = tree.getCellGroup(level,idx_group);
            container->setIdxGlobal(nb_block_before_me);
        } else {
            auto* container = tree.getParticleGroup(idx_group);
            container->setIdxGlobal(int(nb_block_before_me));
        }
        ++nb_block_before_me;
     }
     return nb_block_before_me;
}

/**
 * This function launch the computaition of the flobal index of every
 * group at every level
 * @author benjamin.dufoyer@inria.fr
 * @param  tree         local group tree
 * @param  conf         MPI conf
 * @param  level_min    [OPTIONNAL] minimum level
 */
template<class GroupOctreeClass>
int set_cell_group_global_index(
          GroupOctreeClass&     tree,
    const inria::mpi_config&    conf,
          int                   level_min = 1
){
    int nb_proc = conf.comm.size();
    if(nb_proc > 1){
        // Can be a task
        int nb_block_before_me = 0;
        // set the idx global on the particle block
        nb_block_before_me = set_cell_group_global_index_at(tree,0,nb_block_before_me,conf,true);

        for(int i = tree.getHeight()-1; i >= level_min ; --i){
            nb_block_before_me = set_cell_group_global_index_at(tree,i,nb_block_before_me,conf);
        }
        conf.comm.bcast(
            &nb_block_before_me,
            1,
            MPI_INT,
            nb_proc-1
        );
        return nb_block_before_me;
    } else {
        int idx_global = 0;
        for(int i = 0 ; i < tree.getNbParticleGroup(); ++i){
            tree.getParticleGroup(i)->setIdxGlobal(idx_global);
            ++idx_global;
        }
        for(int i = tree.getHeight()-1; i >= 1 ; --i){
            for(int j = 0; j < tree.getNbCellGroupAtLevel(i) ; ++j ){
                tree.getCellGroup(i,j)->setIdxGlobal(idx_global);
                ++idx_global;
            }
        }
        return idx_global;
    }
}


/**
 * This function add the blocks for the M2M operation
 * 1) Compute the min and max morton index of my distribution at the level
 * 2) Check if i have this morton index at the upper level
 * 3) Share this information at my neighboor
 * 4) Send the morton index needed and post recv of block
 * 5) Send the block to my neihboor if needed
 * 6) Add block to the tree
 *
 * [RESTRICTION] You need to add the LET to the tree BEFORE calling this
 * function
 *
 * @author benjamin.dufoyer@inria.fr
 * @param  tree     LET GroupTree
 * @param  conf     MPI conf
 * @param  level    Level to check
 */
template<class GroupOctreeClassWithLET>
void send_get_block_M2M_at_level(
            GroupOctreeClassWithLET&   tree,
    const   inria::mpi_config&  conf,
            int                 level
){
    // structure for sending message
    struct sending_cell_structure_M2M{
        int idx_global_block;
        MortonIndex start_index;
        MortonIndex end_index;
        int nb_leaf_in_block;
        int idx_global_particle_block = 0;
    };
    // boolean to know if we are at the leaf level
    bool leaf_level = (tree.getHeight()-1 == level);
    // get the block_size
    int block_size = tree.getNbElementsPerBlock();
    // Prepare buffer for sending to other proc
    sending_cell_structure_M2M block_needer_min{-1,0,0,0};
    sending_cell_structure_M2M block_needer_max{-1,0,0,0};
    std::vector<MortonIndex> m_idx_min(block_size,0);
    std::vector<MortonIndex> m_idx_max(block_size,0);
    std::vector<FSize>       nb_particle_min(0,0);
    std::vector<FSize>       nb_particle_max(0,0);
    //  IDEA can be a task
    //  Compute the minimum morton index of my distribution
    //  iterate on every group
    for(int idx_group = 0 ; idx_group < tree.getNbCellGroupAtLevel(level);++idx_group){
        auto* container = tree.getCellGroup(level,idx_group);
        // if the block is mine
        if(container->isMine()){
            // get the symbolic information of the block
            block_needer_min.idx_global_block = container->getIdxGlobal();
            block_needer_min.start_index = container->getStartingIndex();
            block_needer_min.end_index = container->getEndingIndex();
            block_needer_min.nb_leaf_in_block = container->getNumberOfCellsInBlock();
            // get every morton index
            for(int idx_cell = 0 ; idx_cell < container->getNumberOfCellsInBlock() ; ++idx_cell){
                m_idx_min[idx_cell] = container->getCellMortonIndex(idx_cell);
            }
            // if it's the leaf level we need to send the number of particule
            // too
            if(leaf_level){
                auto* container_particle = tree.getParticleGroup(idx_group);
                block_needer_min.idx_global_particle_block = container_particle->getIdxGlobal();
                nb_particle_min.resize(block_size,0);
                for(int idx_cell = 0 ; idx_cell < container_particle->getNumberOfLeavesInBlock() ; ++idx_cell){
                    nb_particle_min[idx_cell] = container_particle->getNbParticlesInLeaf(idx_cell);
                }
            }
            // break the loop
            break;
        }
    }
    // IDEA can be a task
    // compute the maximum morton index of my distribution
    // iterate on every groups
    for(int idx_group = tree.getNbCellGroupAtLevel(level)-1 ; idx_group >= 0; --idx_group){
        auto* container = tree.getCellGroup(level,idx_group);
        // if the block is Mine
        if(container->isMine()){
            // stock symbolic information
            block_needer_max.idx_global_block = container->getIdxGlobal();
            block_needer_max.start_index = container->getStartingIndex();
            block_needer_max.end_index = container->getEndingIndex();
            block_needer_max.nb_leaf_in_block = container->getNumberOfCellsInBlock();
            // get every morton index
            for(int idx_cell = 0 ; idx_cell < container->getNumberOfCellsInBlock() ; ++idx_cell){
                m_idx_max[idx_cell] = container->getCellMortonIndex(idx_cell);
            }
            // if it's the leaf level we need to send to number of particule
            // per leaf
            if(leaf_level){
                auto* container_particle = tree.getParticleGroup(idx_group);
                block_needer_max.idx_global_particle_block = container_particle->getIdxGlobal();
                nb_particle_max.resize(block_size,0);
                for(int idx_cell = 0 ; idx_cell < container_particle->getNumberOfLeavesInBlock() ; ++idx_cell){
                    nb_particle_max[idx_cell] = container_particle->getNbParticlesInLeaf(idx_cell);
                }
            }
            break;
        }
    }
    // compute the max and the min morton Index att the upper level

    // Now we have our max and our min at the current level
    // Now we want to check if we have the parents of our min and our max
    bool flag_min = false;
    bool flag_max = false;

    // MPI info
    int nb_proc = conf.comm.size();
    int my_rank = conf.comm.rank();

    // reception buffer
    // Symbolic block buffer
    sending_cell_structure_M2M buffer_right_neighbor{-1,0,0,0};
    sending_cell_structure_M2M buffer_left_neighbor{-1,0,0,0};
    // Morton index buffer
    std::vector<MortonIndex> buffer_m_idx_right(block_size,0);
    std::vector<MortonIndex> buffer_m_idx_left(block_size,0);
    // number of particle buffer
    std::vector<FSize>       buffer_nb_particle_right(0,0);
    std::vector<FSize>       buffer_nb_particle_left(0,0);

    // flag for neighboot
    bool flag_right_neighboor = false;
    bool flag_left_neighboor  = false;

    // array of request
    int nb_message = 0;
    inria::mpi::request tab_mpi_status[12];

    // if i'm 0, i don't need a block from left
    if(my_rank == 0){
        flag_min = true;
        flag_left_neighboor = true;
    }
    // if i'm the last proc, i don't need block from right
    if(my_rank == nb_proc-1){
        flag_max = true;
        flag_right_neighboor = true;
    }

    // Now we need to send to the neighboor if we need a block, and recv if
    // he need block

    // First send and recv from left

    if(my_rank != 0){
        /////////////////////////////////////////////
        //// SENDING
        /////////////////////////////////////////////
        // Send symbolic information of my min block
        tab_mpi_status[nb_message] =
            conf.comm.isend(
                &block_needer_min,
                sizeof(sending_cell_structure_M2M),
                MPI_CHAR,
                my_rank-1,1);
        ++nb_message;
        // send the morton index of the min block
        tab_mpi_status[nb_message] =
            conf.comm.isend(
                m_idx_min.data(),
                int(sizeof(MortonIndex)*block_size),
                MPI_CHAR,
                my_rank-1,2);
        ++nb_message;
        // if it's the leaf level
        if(leaf_level){
            // send the number of particle of the particle block attached
            tab_mpi_status[nb_message] =
                conf.comm.isend(
                    nb_particle_min.data(),
                    int(sizeof(FSize)*block_size),
                    MPI_CHAR,
                    my_rank-1,3);
                    ++nb_message;
        }
        // recev the symbolic block from my left neighbor
        tab_mpi_status[nb_message] =
            conf.comm.irecv(
                &buffer_left_neighbor,
                sizeof(sending_cell_structure_M2M),
                MPI_CHAR,
                my_rank-1,1);
        ++nb_message;
        // recv the morton index of the block send by my left neighbor
        tab_mpi_status[nb_message] =
            conf.comm.irecv(
                buffer_m_idx_left.data(),
                int(sizeof(MortonIndex)*block_size),
                MPI_CHAR,
                my_rank-1,2);
        ++nb_message;
        // if it's the leaf level
        if(leaf_level){
            // need to recev the number of particle of the particle attached
            buffer_nb_particle_left.resize(block_size,0);
            tab_mpi_status[nb_message] =
                conf.comm.irecv(
                    buffer_nb_particle_left.data(),
                    int(sizeof(FSize)*block_size),
                    MPI_CHAR,
                    my_rank-1,3);
                    ++nb_message;
        }
    }
    // Send and recv from right
    if(my_rank != nb_proc-1){
        // send my block max
        tab_mpi_status[nb_message] =
            conf.comm.isend(
                &block_needer_max,
                sizeof(sending_cell_structure_M2M),
                MPI_CHAR,
                my_rank+1,1);
        ++nb_message;
        // send the morton index of the right block
        tab_mpi_status[nb_message] =
            conf.comm.isend(
                m_idx_max.data(),
                int(sizeof(MortonIndex)*block_size),
                MPI_CHAR,
                my_rank+1,2);
        ++nb_message;
        // if it's the leaf level
        if(leaf_level){
            // send the number of particle of the particle block attached
            tab_mpi_status[nb_message] =
                conf.comm.isend(
                    nb_particle_max.data(),
                    int(sizeof(FSize)*block_size),
                    MPI_CHAR,
                    my_rank+1,3);
            ++nb_message;
        }

        // recv the block from the right
        tab_mpi_status[nb_message] =
            conf.comm.irecv(
                &buffer_right_neighbor,
                sizeof(sending_cell_structure_M2M),
                MPI_CHAR,
                my_rank+1,1);
        ++nb_message;
        // recv the morton index of the right
        tab_mpi_status[nb_message] =
            conf.comm.irecv(
                buffer_m_idx_right.data(),
                int(sizeof(MortonIndex)*block_size),
                MPI_CHAR,
                my_rank+1,2);
        ++nb_message;
        // if it's leaf level
        if(leaf_level){
            // recv number of particle
            buffer_nb_particle_right.resize(block_size,0);
            tab_mpi_status[nb_message] =
                conf.comm.irecv(
                    buffer_nb_particle_right.data(),
                    int(sizeof(FSize)*block_size),
                    MPI_CHAR,
                    my_rank+1,3);
            ++nb_message;
        }

    }
    // Wait all request
    inria::mpi::request::waitall(nb_message,tab_mpi_status);

    // Now we have the min and the max block at the level L
    // But now we need to send the block needer of this block
    // to insert task on him with starPU

    // buffer to recv blocks
    sending_cell_structure_M2M block_from_left;
    sending_cell_structure_M2M block_from_right;

    // buffer to recv morton index
    std::vector<MortonIndex> m_idx_from_left(block_size,-1);
    std::vector<MortonIndex> m_idx_from_right(block_size,-1);

    // buffer to send morton index
    std::vector<MortonIndex> m_idx_to_send_right(block_size,-1);
    std::vector<MortonIndex> m_idx_to_send_left(block_size,-1);
    // buffer to send symbolic information
    sending_cell_structure_M2M block_to_left{-1,0,0,0};
    sending_cell_structure_M2M block_to_right{-1,0,0,0};


    nb_message = 0;
    // IDEA can be a task
    // we post the recv for the left block
    if(!flag_min){
        // posting reception of the block
        tab_mpi_status[nb_message] =
            conf.comm.irecv(
                &block_from_left,
                sizeof(sending_cell_structure_M2M),
                MPI_CHAR,
                my_rank-1,2);
        ++nb_message;
        tab_mpi_status[nb_message] =
            conf.comm.irecv(
                m_idx_from_left.data(),
                int(sizeof(MortonIndex)*block_size),
                MPI_CHAR,
                my_rank-1,3);
        ++nb_message;
    }
    // IDEA can be a task
    // we post the recv for the right block
    if(!flag_max){
        // posting the reception buffer
        tab_mpi_status[nb_message] =
            conf.comm.irecv(
                &block_from_right,
                sizeof(sending_cell_structure_M2M),
                MPI_CHAR,
                my_rank+1,2);
        ++nb_message;
        tab_mpi_status[nb_message] =
            conf.comm.irecv(
                m_idx_from_right.data(),
                int(sizeof(MortonIndex)*block_size),
                MPI_CHAR,
                my_rank+1,3);
        ++nb_message;
    }

    // IDEA Can be a task
    // if i need to send to the right
    if(!flag_right_neighboor){
        bool flag = false;
        // seeking the first block who is mine at the upper level
        for(int i = (tree.getNbCellGroupAtLevel(level-1)-1) ; i >= 0  ; --i){
            auto* container = tree.getCellGroup(level-1,i);
            // if the block is mine
            if(container->isMine()){
                // stock symbolic information
                block_to_right.idx_global_block = container->getIdxGlobal();
                block_to_right.start_index = container->getStartingIndex();
                block_to_right.end_index = container->getEndingIndex();
                block_to_right.nb_leaf_in_block = container->getNumberOfCellsInBlock();
                // stock the morton index of the block
                for(int idx_cell = 0 ; idx_cell < container->getNumberOfCellsInBlock(); ++idx_cell ){
                    m_idx_to_send_right[idx_cell] = container->getCellMortonIndex(idx_cell);
                }
                // put the flag on true
                flag = true;
                // send the 2 buffer
                tab_mpi_status[nb_message] =
                conf.comm.isend(
                    &block_to_right,
                    sizeof(sending_cell_structure_M2M),
                    MPI_CHAR,
                    my_rank+1,2);
                ++nb_message;
                tab_mpi_status[nb_message] =
                conf.comm.isend(
                    &m_idx_to_send_right.data()[0],
                    int(sizeof(MortonIndex)*block_size),
                    MPI_CHAR,
                    my_rank+1,3);
                ++nb_message;
                break;
            }
        }
        // we don't have block at the upper level
        if(!flag){
            // send fake block
            tab_mpi_status[nb_message] =
            conf.comm.isend(
                &block_to_right,
                sizeof(sending_cell_structure_M2M),
                MPI_CHAR,
                my_rank+1,2);
            ++nb_message;
            tab_mpi_status[nb_message] =
            conf.comm.isend(
                m_idx_to_send_right.data(),
                int(sizeof(MortonIndex)*block_size),
                MPI_CHAR,
                my_rank+1,3);
            ++nb_message;
        }
    }
    // IDEA Can be a task

    if(!flag_left_neighboor){
        bool flag = false;
        // seek the first block who is mine
        for(int i = 0 ; i < tree.getNbCellGroupAtLevel(level-1) ; ++i){
            auto* container = tree.getCellGroup(level-1,i);
            // send the first left block who is mine
            if(container->isMine()){
                // stock symbolic information
                block_to_left.idx_global_block = container->getIdxGlobal();
                block_to_left.start_index = container->getStartingIndex();
                block_to_left.end_index = container->getEndingIndex();
                block_to_left.nb_leaf_in_block = container->getNumberOfCellsInBlock();
                // stock morton index
                for(int idx_cell = 0 ; idx_cell < container->getNumberOfCellsInBlock(); ++idx_cell ){
                    m_idx_to_send_left[idx_cell] = container->getCellMortonIndex(idx_cell);
                }
                // put the flag on true
                flag = true;
                // send block
                tab_mpi_status[nb_message] =
                conf.comm.isend(
                    &block_to_left,
                    sizeof(sending_cell_structure_M2M),
                    MPI_CHAR,
                    my_rank-1,2);
                ++nb_message;
                tab_mpi_status[nb_message] =
                conf.comm.isend(
                    m_idx_to_send_left.data(),
                    int(sizeof(MortonIndex)*block_size),
                    MPI_CHAR,
                    my_rank-1,3);
                ++nb_message;
                break;
            }
        }
        // we don't have block at the upper level
        if(!flag){
            // send fake block
            tab_mpi_status[nb_message] =
            conf.comm.isend(
                &block_to_left,
                sizeof(sending_cell_structure_M2M),
                MPI_CHAR,
                my_rank-1,2);
            ++nb_message;
            tab_mpi_status[nb_message] =
            conf.comm.isend(
                m_idx_to_send_left.data(),
                int(sizeof(MortonIndex)*block_size),
                MPI_CHAR,
                my_rank-1,3);
            ++nb_message;
        }
    }


    // Wait for the send/recv
    if(nb_message > 0)
        inria::mpi::request::waitall(nb_message,tab_mpi_status);

    // now i have the block needed for the M2M
    // now we need to add this block

    // We add the block, if the idx_global_block is -1, the block
    // is invalid so we don't need him, and we don't need to add him to the tree
    if(!flag_min && block_from_left.idx_global_block != -1){
        tree.insert_block(block_from_left,m_idx_from_left,level-1);
    }
    if(!flag_max && block_from_right.idx_global_block != -1){
        tree.insert_block(block_from_right,m_idx_from_right,level-1);
    }
    if(!flag_right_neighboor && buffer_right_neighbor.idx_global_block != -1){
        tree.insert_block(buffer_right_neighbor,buffer_m_idx_right,level,&buffer_nb_particle_right);
    }
    if(!flag_left_neighboor && buffer_left_neighbor.idx_global_block != -1){
        tree.insert_block(buffer_left_neighbor,buffer_m_idx_left,level,&buffer_nb_particle_left);
    }
}

/**
 * This function exchange blocks with neighbors proc
 * The left proc have the first block
 * The right proc have the last block
 *
 * The blocks send are the block who have the boolean "isMine" on 1 on the
 * sender
 *
 * @author benjamin.dufoyer@inria.fr
 * @param  tree         The group tree
 * @param  conf
 * @param  level_min    [OPTIONNAL] minimum level to apply this function
 */
template<class GroupOctreeClass>
void send_get_block_M2M(
            GroupOctreeClass&   tree,
    const   inria::mpi_config&  conf,
            int                 level_min = 1
){
    int nb_proc = conf.comm.size();
    // if we have less than 1 proc, we don't need to exchange block
    if(nb_proc > 1){
        // get the M2M block at every level
        for(int i = tree.getHeight()-1 ; i > level_min ; --i){
            send_get_block_M2M_at_level(tree,conf,i);
        }
    }
}


}


#endif /*_FDISTRIBUTED_GROUPTREE_BUILDER_HPP_*/
