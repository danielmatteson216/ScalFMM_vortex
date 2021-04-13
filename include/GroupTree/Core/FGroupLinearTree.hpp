#ifndef _FGROUP_LINEAR_TREE_HPP_
#define _FGROUP_LINEAR_TREE_HPP_

#include <vector>
#include "Utils/FGlobal.hpp"
#include "Utils/FLog.hpp"
//#include "inria/algorithm/distributed/mpi.hpp"
//using FReal = double;

template<class node_t>
class FGroupLinearTree {



protected:

    int block_size;   //<
    int nb_block;     //<

    // Copy of the MPI conf 
    const inria::mpi_config mpi_conf;  //<

    std::vector<node_t>*                           linear_tree;         //<
    std::vector<std::pair<MortonIndex,MortonIndex>> index_particle_distribution;  //<
    bool unknow_index_particle_distribution = true;  //<

public:

////////////////////////////////////////////////
// constructor
////////////////////////////////////////////////

    /**
     * FBlockedLinearTree  Constructor of blocked linear tree
     * @author benjamin.dufoyer@inria.fr
     * @param  in_block_size    Block size needed
     * @param  in_linear_tree   Linear tree
     * @param  in_box_center    Box Center of particle container
     * @param  in_box_width     Box Width of particle container
     * @warning We copy the MPI comm because O3 compilation fail the utest
     */
    FGroupLinearTree(const inria::mpi_config conf):
        mpi_conf(conf),
        index_particle_distribution(conf.comm.size())
    {
        linear_tree = new std::vector<node_t>[1];
    }

////////////////////////////////////////////////////
// Function of initialisation
////////////////////////////////////////////////////

    /**
     * This function create a blocked linear tree from the current distributed
     * linear tree
     * This function stock the linear tree with his adress
     * @author benjamin.dufoyer@inria.fr
     * @param  in_linear_tree linear tree
     * @param  in_block_size  block size
     */
    void create_local_group_linear_tree(
            std::vector<node_t>* in_linear_tree,
            int in_block_size
    ){
        this->create(in_linear_tree,in_block_size);
    }

    /**
     * this function create a blocked linear tree from the current distributed
     * linear tree and she redistribute block according to the block size
     * the function stock the linear tree with his adress
     * @author benjamin.dufoyer@inria.fr
     * @param  in_linear_tree linear tree
     * @param  in_block_size  blocksize needed
     */
    void create_global_group_linear_tree(
            std::vector<node_t>* in_linear_tree,
            int in_block_size
    ){
        this->create(in_linear_tree,in_block_size);
        this->redistribute_block();
    }

    void create(
        std::vector<node_t>* in_linear_tree,
        int in_block_size
    ){

        this->block_size    = in_block_size;
        this->linear_tree   = in_linear_tree;
        this->nb_block      = (int)in_linear_tree->size()/in_block_size;
        if(this->linear_tree->size()%this->block_size != 0)
            this->nb_block += 1;
    }

////////////////////////////////////////////////
// destructor
////////////////////////////////////////////////

    ~FGroupLinearTree(){
        linear_tree = nullptr;
    }

////////////////////////////////////////////////
// Function
////////////////////////////////////////////////

    /**
     * redistribute_block redistribute leaf of the linear_tree with the good
     * block size. For N proc, N-1 proc have the same number of leaf,
     * the rest is for the proc N
     * @author benjamin.dufoyer@inria.fr
     */
    void redistribute_block(){

        dstr_grp_tree_builder::parrallel_build_block(
                        this->mpi_conf,
                        this->linear_tree,
                        this->block_size);
        //Update nb_block
        if(this->linear_tree->size()%block_size == 0)
            this->nb_block = (int)this->linear_tree->size()/block_size;
        else
            this->nb_block = (int)this->linear_tree->size()/block_size+1;

    }

    size_t get_nb_leaf() const{
        return this->linear_tree->size();
    }

    int get_nb_block() const{
        return this->nb_block;
    }

    int get_block_size() const{
        return this->block_size;
    }

    /**
     * get_block_size_at return the block size of the number of the block
     * placed in parametter,
     * [INFO] first block is 0
     * [INFO] last block is this->nb_block-1
     * @author benjamin.dufoyer@inria.fr
     * @param  num_block number of the block
     * @return size of the block
     */
    int get_block_size_at(int num_block) const{
        FAssertLF(num_block < this->nb_block);
        int size;
        if(num_block == this->nb_block-1){
            size = this->linear_tree->size() - ((this->nb_block-1)*this->block_size);
        } else {
            size = this->block_size;
        }
        return size;
    }

    /**
     * get_leaf_at return the leaf at the position placed in parameter
     * @author benjamin.dufoyer@inria.fr
     * @param  position position of the leaf
     * @return the leaf
     */
    node_t get_leaf_at(int position){
        return this->linear_tree->at(position);
    }

    /**
     * get_leaf_at return the leaf at the position placed in parameter
     * @author benjamin.dufoyer@inria.fr
     * @param  position position of the leaf
     * @return the leaf
     */
    node_t at(int position){
        return this->get_leaf_at(position);
    }

    size_t get_leaf_level() const{
        return this->linear_tree->back().level;
    }

    size_t get_tree_height() const{
        return this->get_leaf_level();
    }

    size_t get_first_morton_index() const{
        return this->linear_tree->front().morton_index;
    }

    size_t get_last_morton_index() const{
        return this->linear_tree->back().morton_index;
    }


    const inria::mpi_config get_mpi_conf() const{
        return mpi_conf;
    }


    /**
     * This function print the information of this current class
     * @author benjamin.dufoyer@inria.fr
     */
    void print_info_tree(){
        std::cout << " nb_leaf : " << this->linear_tree->size() << std::endl;
        std::cout << " nb_block : " << nb_block << std::endl;
        std::cout << " block_size : " << block_size << std::endl;
        for(unsigned i = 0 ; i < this->linear_tree->size() ; i++){
            std::cout << linear_tree->data()[i] << std::endl;
        }
    }

    /**
     * This function return the linear tree associated with this calss
     * @author benjamin.dufoyer@inria.fr
     */
    std::vector<node_t>* get_tree(){
        return this->linear_tree;
    }

    /**
     * This function fill the index_particle_distribution with a call of the function
     * from FDistributedGroupTreeBuilder
     * @author benjamin.dufoyer@inria.fr
     * @param  particle_container [description]
     */
    template<class particle_t>
    void set_index_particle_distribution( std::vector<particle_t>& particle_container)
    {
        unknow_index_particle_distribution = false;
        if(this->mpi_conf.comm.size() > 1){
            dstr_grp_tree_builder::share_particle_division(
                  this->mpi_conf,
                  particle_container,
                  this->index_particle_distribution);
        } else {
            this->index_particle_distribution.resize(1);
            std::pair<MortonIndex,MortonIndex> my_distrib;
            my_distrib.first  = particle_container.front().morton_index;
            my_distrib.second = particle_container.back().morton_index;
            this->index_particle_distribution[0] = my_distrib;
        }
    }

    /**ad
     * this function do a update of the current particle distribution
     * and the pair is put in parameter
     * @author benjamin.dufoyer@inria.fr
     * @param  new_distrib [description]
     */
    void update_index_particle_distribution(std::pair<MortonIndex,MortonIndex> new_distrib){
        unknow_index_particle_distribution = false;
        if(this->mpi_conf.comm.size() > 1){
            dstr_grp_tree_builder::share_particle_division(
                this->mpi_conf,
                new_distrib,
                this->index_particle_distribution);
        } else {
            this->index_particle_distribution.resize(1);
            this->index_particle_distribution[0] = new_distrib;
        }
    }

    /**
     * this function return a pointer of the total particule repartition
     * @author benjamin.dufoyer@inria.fr
     */
    std::vector<std::pair<MortonIndex,MortonIndex>>
    get_index_particle_distribution(){
        // TO get the particle repartition, you will compute it before
        FAssert(!unknow_index_particle_distribution);
        return this->index_particle_distribution;
    }


  std::vector<MortonIndex> get_index_particle_distribution_implicit(){

    std::vector<MortonIndex> distribution( (this->index_particle_distribution.size()*2) /*+2*/,-1); // Pouruoi +2 OC ?
    if(this->mpi_conf.comm.size() == 0){
      for(unsigned i = 1; i < distribution.size() ; ++i ){
	distribution[i] = this->index_particle_distribution[0].second;
      }
    } 
    else {
  //    int idx_vect = 0 ;
   //   distribution[0] =  ;
      distribution[1] = this->index_particle_distribution[0].second ;

      for(unsigned i = 1 ; i < this->index_particle_distribution.size() ; ++i){
          distribution[2*i]   = this->index_particle_distribution[i-1].second ;
          distribution[2*i+1] = this->index_particle_distribution[i].second;
      }
//      int idx_vect = static_cast<int>(2*this->index_particle_distribution.size() );
//      ///////////// TO REMOVE ???
//      distribution[idx_vect] = this->index_particle_distribution[this->index_particle_distribution.size()-1].second;
//      ++idx_vect;
//      distribution[idx_vect] = this->index_particle_distribution[this->index_particle_distribution.size()-1].second;
    }
    return distribution;
  }

    /**
     * this function return the particle distribution for a rank of proc
     * put in parameter
     * @author benjamin.dufoyer@inria.fr
     * @param  proc_rank rank of the proc
     * @return a pair of morton index
     */
    std::pair<MortonIndex,MortonIndex> get_index_particle_distribution_at(unsigned proc_rank){
        // TO get the particle repartition, you will compute it before
        FAssert(!unknow_index_particle_distribution);
        FAssert(proc_rank < this->index_particle_distribution.size());
        return this->index_particle_distribution.data()[proc_rank];
    }

    /**
     * return the max morton index of the left proc
     * it's needed for the distributed construction of the groupTree
     * [RESTRICTION] you will compute the repartition before calling this function
     * @author benjamin.dufoyer@inria.fr
     * @return max left morton index, -1 if he doesn't exist
     */
    MortonIndex get_left_limit(){
        // Check if the repartition have been compute
        FAssert(!unknow_index_particle_distribution);
        int my_rank = this->mpi_conf.comm.rank();
        MortonIndex left_limit = -1;
        if(my_rank != 0){
            left_limit = static_cast<MortonIndex>(this->index_particle_distribution[my_rank-1].second);
        }
        return left_limit;
    }


    /**
     * This function is used to show the FGroupLinearTee more easly
     * @author benjamin.dufoyer@inria.fr
     */
    friend
    std::ostream& operator<<(std::ostream& os, const FGroupLinearTree& n) {
    return os << "--> Number of leaf : " << n.get_nb_leaf()
            << "\n first leaf : "      << n.get_first_morton_index()
            << "\n last  leaf : "      << n.get_last_morton_index()
            << "\n block_size "          << n.get_block_size()
            << "\n number of block : "   << n.get_nb_block();
    }


};

#endif //_FGROUP_LINEAR_TREE_HPP_
