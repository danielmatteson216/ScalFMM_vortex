// See LICENCE file at project root

// ==== CMAKE =====
// @FUSE_BLAS
// @FUSE_MPI
// ================

#include "Utils/FGlobal.hpp"
// include algo for linear tree
#include "inria/algorithm/distributed/mpi.hpp"
#include "inria/linear_tree/balance_tree.hpp"
// tree class
#include "GroupTree/Core/FGroupTree.hpp"
// symbolic data
#include "Components/FSymbolicData.hpp"
// cell class
#include "Kernels/Chebyshev/FChebCell.hpp"
// parameter
#include "Utils/FParameters.hpp"
#include "Utils/FParameterNames.hpp"
// GroupParticleContianer
#include "GroupTree/Core/FP2PGroupParticleContainer.hpp"
// file loader
#include "Files/FMpiFmaGenericLoader.hpp"
// FBox
#include "Adaptive/FBox.hpp"
// Group linear tree
#include "GroupTree/Core/FGroupLinearTree.hpp"
// Function for GroupLinearTree
#include "GroupTree/Core/FDistributedGroupTreeBuilder.hpp"
#include "GroupTree/Core/FDistributedLETGroupTreeValidator.hpp"
#include <memory>

#include "FUTester.hpp"

static const int ORDER = 6;
using FReal               = double;

struct particle_t {
    using position_t = FPoint<FReal>;
    position_t pos;
    FReal phi;

    std::size_t morton_index;
    const auto& position() const {
        return pos;
    }
    const FPoint<FReal>& getPosition(){
        return pos;
    }
    const auto& physicalValue() const{
        return phi;
    }
    const auto& getPositions() const {
        return pos;
    }
    int weight() const { return 1;}
    friend constexpr std::size_t morton_index(const particle_t& p) {
        return p.morton_index;
    }
};

class TestLetGroupTree : public FUTesterMpi<TestLetGroupTree>{


    template<class GroupCellClass
             ,class GroupCellUpClass
             ,class GroupCellDownClass
             ,class GroupCellSymbClass
             ,class GroupContainerClass
             ,class GroupOctreeClass>
    void RunTest(){
        const int TreeHeight = 5;
        const int level = TreeHeight-1;
        const int groupSize = 32;
        // Definition of the particle type


        // Load the Quentin MPI
        inria::mpi::communicator mpi_comm(app.global().getComm());

        // Selection of the file
        const std::string parFile( (sizeof(FReal) == sizeof(float))?
                                       "Test/DirectFloatbfma":
                                       "test20k.fma");
        std::string filename(SCALFMMDataPath+parFile);

        // Load the file
        FMpiFmaGenericLoader<FReal> loader(filename, app.global());

        // declare vector to stock particle
        std::vector<particle_t> myParticles(loader.getMyNumberOfParticles());


        const std::size_t max_level = sizeof(particle_t::morton_index) * 8 / 3;
        // define a box, used in the sort
        const FBox<FPoint<FReal>> box{loader.getBoxWidth(),loader.getCenterOfBox()};

        // iterate on all of my particles
        for(FSize idxPart = 0; idxPart <loader.getMyNumberOfParticles();++idxPart){
            particle_t tmp;
            // get the current particles
            loader.fillParticle(&tmp.pos,&tmp.phi);
            // set the morton index of the current particle at the max_level
            tmp.morton_index =  inria::linear_tree::get_morton_index(
                                tmp.pos, box, max_level);
            // set the weight of the particle
            tmp.phi = 0.1;
            //  add the particle to my vector of particle
            myParticles[idxPart] = tmp;
        }
        // Sort particules
        inria::sort(mpi_comm,myParticles,
                [](const auto& p1, const auto& p2) {
                    return p1.morton_index < p2.morton_index;
                });
        // Create linear tree
        auto linear_tree = inria::linear_tree::create_balanced_linear_tree_at_level(
            mpi_comm,
            level,
            box,                                                     myParticles);
        // Create empty instance of group linear tree
        FGroupLinearTree<decltype(linear_tree)::value_type>group_linear_tree{mpi_comm};


        // Fill the group linear tree
        group_linear_tree.create_local_group_linear_tree(
            &linear_tree,
            groupSize);

        // Redistribute particule according to the linear tree
        inria::linear_tree::redistribute_particles(mpi_comm,
                                               linear_tree,
                                               myParticles);


       // Modify the Morton Index to accord him to the level
       for(unsigned i = 0 ; i < myParticles.size(); ++i){
            myParticles.at(i).morton_index = inria::linear_tree::get_morton_index(
                myParticles.at(i).pos, box, level);
        }

        group_linear_tree.set_index_particle_distribution(myParticles);


        GroupOctreeClass  localGroupTree = GroupOctreeClass::template get_block_tree_instance<GroupCellSymbClass, GroupCellUpClass, GroupCellDownClass, GroupContainerClass>(TreeHeight,
                                             groupSize,
                                             loader.getCenterOfBox(),
                                             loader.getBoxWidth());
         std::cout << " Creating Tree" << std::endl;

        localGroupTree.create_tree(group_linear_tree,myParticles);
        std::cout << "Tree construit" << std::endl;
        int nb_block = dstr_grp_tree_builder::set_cell_group_global_index(localGroupTree,mpi_comm);
        // Check if we have no duplication of global index
        std::vector<bool> group_index_checker(nb_block,false);
        // Check on particle group
        for(int i = 0 ; i < localGroupTree.getNbParticleGroup() ; ++i){
            auto* container = localGroupTree.getParticleGroup(i);
            int idx_global = container->getIdxGlobal();
            uassert(group_index_checker[idx_global] == false);
            group_index_checker[idx_global] = true;
        }
        // Check on cell group
        for(int j = 0 ; j < localGroupTree.getHeight() ; ++j){
            for(int i = 0 ; i < localGroupTree.getNbCellGroupAtLevel(j); ++i){
                auto* container = localGroupTree.getCellGroup(j,i);
                int idx_global = container->getIdxGlobal();
                uassert(group_index_checker[idx_global] == false);
                group_index_checker[idx_global] = true;
            }
        }
        // Create LET
        localGroupTree.create_LET(group_linear_tree);
        // launch the let checker
        bool flag = dstr_grp_tree_vldr::validate_group_tree(localGroupTree,mpi_comm);
        // if the LET is correct, flag is true
        uassert(flag);
    }



    void TestLet(){
        using GroupCellClass      = FChebCell<FReal, ORDER>;
        using GroupCellUpClass    = typename GroupCellClass::multipole_t;
        using GroupCellDownClass  = typename GroupCellClass::local_expansion_t;
        using GroupCellSymbClass  = FSymbolicData;
        using GroupContainerClass = FP2PGroupParticleContainer<FReal>;
        using GroupOctreeClass    = FGroupTree<FReal,
                                        GroupCellSymbClass,
                                        GroupCellUpClass,
                                        GroupCellDownClass, GroupContainerClass, 1, 4, FReal>;

        RunTest<GroupCellClass,
                GroupCellUpClass,
                GroupCellDownClass,
                GroupCellSymbClass,
                GroupContainerClass,
                GroupOctreeClass>();
    }

    void SetTests(){
        AddTest(&TestLetGroupTree::TestLet,"Test the building of the LET ");
    }

public:
    TestLetGroupTree(int argc, char ** argv):
    FUTesterMpi(argc,argv){
    }

};

TestClassMpi(TestLetGroupTree);
