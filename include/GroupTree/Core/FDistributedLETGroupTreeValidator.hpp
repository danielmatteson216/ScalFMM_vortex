// ==== CMAKE ====
// @FUSE_MPI
// ================
//


#ifndef _FDISTRIBUTED_LET_GROUPTREE_VALIDATOR_
#define _FDISTRIBUTED_LET_GROUPTREE_VALIDATOR_


#include "inria/algorithm/distributed/mpi.hpp"


namespace dstr_grp_tree_vldr{

/**
 * This function check the level of the LetGroupTree to check if we forget
 * a group
 * The principe is simple. We compute every interaction of every cell at the
 * level, we check if we have the morton index of the interaction in our tree
 * If we don't have this index, we send a request to the proc who have this
 * index to check if he exist, if he exist, it's a error
 * @author benjamin.dufoyer@inria.fr
 * @param  tree  localGroupTree
 * @param  level level to check
 * @param  conf  conf MPI
 * @return true if it's ok, false if we forget a group
 */
template<class GroupOctreeClass>
bool validate_group_tree_at_level(
    GroupOctreeClass& tree,
    int level,
    const inria::mpi_config& conf
){
    // MPI information
    const int nb_proc = conf.comm.size();
    const int my_rank = conf.comm.rank();

    // Compute my min and my max morton index at the level
    MortonIndex min_morton_index_at_level = 0;
    MortonIndex max_morton_index_at_level = 0;
    for(int i = 0 ; i < tree.getNbCellGroupAtLevel(level) ; ++i){
        auto* container =  tree.getCellGroup(level,i);
        if(container->isMine()){
            min_morton_index_at_level = container->getStartingIndex();
            break;
        }
    }
    for(int i = tree.getNbCellGroupAtLevel(level)-1; i >= 0 ; --i){
        auto* container =  tree.getCellGroup(level,i);
        if(container->isMine()){
            max_morton_index_at_level = container->getEndingIndex();
            break;
        }
    }

    // Sharing my interval and getting interval from all proc
    std::pair<MortonIndex,MortonIndex> my_interval(min_morton_index_at_level,max_morton_index_at_level);
    std::vector<std::pair<MortonIndex,MortonIndex>> all_interval(nb_proc);
    conf.comm.allgather(&my_interval,
                        sizeof(my_interval),
                        MPI_CHAR,
                        all_interval.data(),
                        sizeof(my_interval),
                        MPI_CHAR);

    // if i have 1 block or more
    // Get all MortonIndex for interaction
    std::vector<MortonIndex> morton_index_not_in_tree(0);
    if(my_interval.second != 0){
        // vector to stock all MortonIndex
        std::vector<MortonIndex> external_interaction(tree.getNbCellGroupAtLevel(level)*tree.getNbElementsPerBlock()*189,0);
        unsigned idx_vector = 0;
        // iterate on every group
        for(int idx_group = 0 ; idx_group < tree.getNbCellGroupAtLevel(level) ; ++idx_group){
            // get the current group
            auto* container = tree.getCellGroup(level,idx_group);
            if(container->isMine()){
                // iterate on every cell
                for(int cell_idx = 0;
                        cell_idx < container->getNumberOfCellsInBlock();
                        ++cell_idx){
                            // Getting the current morton index
                            MortonIndex curr_m_idx  = container->getCellMortonIndex(cell_idx);
                            MortonIndex interactionsIndexes[189];
                            int interactionsPosition[189];
                            FTreeCoordinate coord(curr_m_idx);
                            // Getting neigbors of the father
                            int counter = coord.getInteractionNeighbors(level,interactionsIndexes,interactionsPosition);
                            for(int idx_neighbor = 0 ; idx_neighbor < counter ; ++idx_neighbor){
                                MortonIndex tmp = interactionsIndexes[idx_neighbor];
                                if(tmp >= min_morton_index_at_level && tmp < max_morton_index_at_level){
                                    // do nothing, it's my interval
                                } else {
                                    //Stock the index
                                    external_interaction[idx_vector] = tmp;
                                    ++idx_vector;
                                }
                            } // end for neigbors
                        } // end for leaf
                    } // end for group
                }
        if(idx_vector > 0){
            FQuickSort<MortonIndex>::QsSequential(external_interaction.data(),idx_vector);
            // vector to have all mortonIndex with no duplicate data
            std::vector<MortonIndex> morton_needed(0);
            MortonIndex last_morton_index = -1;
            for(unsigned i = 0 ; i < idx_vector ; ++i){
                if(external_interaction[i] != last_morton_index){
                    morton_needed.push_back(external_interaction[i]);
                    last_morton_index = external_interaction[i];
                }
            }
            // free the old vector
            std::vector<MortonIndex>().swap(external_interaction);
            // vector to stock morton index who are not in the tree

            for(unsigned i = 0 ; i < morton_needed.size(); ++i ){
                bool flag = false;
                MortonIndex current_morton_index = morton_needed[i];
                for(int j = 0 ; j < tree.getNbCellGroupAtLevel(level); ++j){
                    auto* container = tree.getCellGroup(level,j);
                    if(!container->isMine()){
                        if(container->isInside(current_morton_index) || container->getEndingIndex() == current_morton_index){
                            flag =true;
                            break;
                        }
                    }
                }
                // if we are here, we don't have the interaction
                if(!flag)
                    morton_index_not_in_tree.push_back(current_morton_index);
            }
        }
    }

    // Now we have all morton index who is not in our tree
    std::vector<unsigned> nb_message_to_send(nb_proc,0);
    std::vector<unsigned> nb_message_to_recev(nb_proc,0);
    for(unsigned i = 0 ; i < morton_index_not_in_tree.size() ;++i ){
        for(unsigned j = 0 ; j < all_interval.size() ; ++j){
            MortonIndex min = all_interval[j].first;
            MortonIndex max = all_interval[j].second;
            if(morton_index_not_in_tree[i] >= min &&  morton_index_not_in_tree[i] <= max ){
                nb_message_to_send[j] += 1;
                break;
            }
        }
    }

    // Send the number of morton index we will send
    conf.comm.alltoall(nb_message_to_send.data(),
                        1,
                        MPI_UNSIGNED,
                        nb_message_to_recev.data(),
                        1,
                        MPI_UNSIGNED);

    // Compute the number of message and the number of morton index
    int nb_morton_index = 0;
    int nb_message =0;
    for(unsigned i = 0 ; i < nb_message_to_recev.size() ; ++i){
        nb_morton_index += nb_message_to_recev[i];
        if(nb_message_to_recev[i] > 0){
            ++nb_message;
        }
        if(nb_message_to_send[i] > 0){
            ++nb_message;
        }
    }

    // declare the reception buffer
    std::vector<MortonIndex> morton_recv(nb_morton_index,0);
    // tab of MPI request to wait the completion
    inria::mpi::request tab_mpi_status[nb_message];

    int idx_message =0;
    unsigned offset = 0;
    // post all reception
    for(unsigned i = 0 ; i < nb_message_to_recev.size(); ++i){
        if(nb_message_to_recev[i] > 0){
            unsigned nb_m_idx = nb_message_to_recev[i];
            tab_mpi_status[idx_message] = conf.comm.irecv(&morton_recv[offset],
                            int(nb_m_idx*sizeof(MortonIndex)),
                            MPI_CHAR,
                            i,1);
            ++idx_message;
            offset += nb_m_idx;
        }
    }

    offset = 0;
    // post all send message
    for(unsigned i = 0 ; i < nb_message_to_send.size(); ++i){
        if(nb_message_to_send[i] > 0){
            unsigned nb_m_idx = nb_message_to_send[i];
            tab_mpi_status[idx_message] = conf.comm.isend(&morton_index_not_in_tree[offset],
                            int(nb_m_idx*sizeof(MortonIndex)),
                            MPI_CHAR,
                            i,1);
            ++idx_message;
            offset += nb_m_idx;
        }
    }

    // Wait all request
    inria::mpi::request::waitall(idx_message,tab_mpi_status);

    offset = 0 ;
    bool flag = true;
    for(unsigned i = 0 ; i < nb_message_to_recev.size() ; ++i ){
        unsigned nb_morton_index_2 = nb_message_to_recev[i];
        for(unsigned j = 0 ; j < nb_morton_index_2 ; ++j ){
            MortonIndex current_idx = morton_recv[j+offset];
            for(int k = 0; k < tree.getNbCellGroupAtLevel(level);++k){
                auto* container = tree.getCellGroup(level,k);
                if(container->isMine()){
                    if(container->isInside(current_idx)){
                        std::cout << " [Error][level "<<level << "] " << current_idx << " on " << my_rank << " Not transfered to " << i << std::endl;
                        flag = false;
                    }
                }
            }
        }
        offset += nb_morton_index_2;
    }

    // return the flag
    return flag;

}

/**
 * This function check every level of the LetGroupTree to know if we forget
 * a group
 * @author benjamin.dufoyer@inria.fr
 * @param  tree local group tree + let
 * @param  conf MPI cong
 * @return true if the tree is ok
 */
template<class GroupOctreeClass>
bool validate_group_tree(
    GroupOctreeClass& tree,
    const inria::mpi_config& conf
){
    bool res = true;
    // check every level
    for(int i = tree.getHeight()-1 ; i > 0 ; --i){
        res = validate_group_tree_at_level(tree,i,conf);
        // if the current level is not good
        if(!res)
            break;
    }
    return res;
}


}

#endif
