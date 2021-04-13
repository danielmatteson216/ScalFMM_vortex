// See LICENCE file at project root
#ifndef FTESTKERNELS_HPP
#define FTESTKERNELS_HPP


#include <iostream>

#include "FAbstractKernels.hpp"
#include "../Containers/FOctree.hpp"
#include "../Utils/FGlobal.hpp"



/**
 * @example testFmmAlgorithm
 * It illustrates how the FTestKernels to validate a FMM core algorithm.
 */


/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class AbstractKernels
* @brief This kernels is a virtual kernels to validate that the fmm algorithm is correctly done on particles.
*
* It should use FTestCell and FTestParticle.
* A the end of a computation, the particles should host then number of particles in the simulation (-1).
*/
template< class CellClass, class ContainerClass>
class FTestKernels  : public FAbstractKernels<CellClass,ContainerClass> {
public:
    /** Default destructor */
    virtual ~FTestKernels(){
    }

    /** Before upward */
    template<class Symb>
    void P2M(typename CellClass::multipole_t* const leaf_multipole,
             const Symb* const /* leaf_symbolic_data */,
             const ContainerClass* const particles)
    {
        // the pole represents all particles under
        leaf_multipole->set(leaf_multipole->get() + particles->getNbParticles());
    }

    template<class Symb>
    void M2M(typename CellClass::multipole_t * const node_multipole,
             const Symb* const /*node_symbolic_data*/ ,
             const typename CellClass::multipole_t * const * const child_multipoles,
             const Symb* const /* child_symbolic_data */ [])
    {
        // A parent represents the sum of the child
        for(int idx = 0 ; idx < 8 ; ++idx){
            if(child_multipoles[idx]){
                node_multipole->set(node_multipole->get() + child_multipoles[idx]->get());
            }
        }
    }

    /** Before Downward */
    template<class Symb>
    void M2L(typename CellClass::local_expansion_t* const node_local_expansion,
             const Symb* const /*node_symbolic_data*/,
             const typename CellClass::multipole_t* const v_item_multipoles[],
             const Symb* const /*v_item_symbolic_data*/[],
             const int /*position*/[],
             const int v_item_data_size)
    {
        // The pole is impacted by what represent other poles
        for(int idx = 0 ; idx < v_item_data_size ; ++idx){
            node_local_expansion->set(node_local_expansion->get()
                                      + v_item_multipoles[idx]->get());
        }
    }

    /** During Downward */
    template<class Symb>
    void L2L(const typename CellClass::local_expansion_t* const node_local_exp,
             const Symb* const /*node_symbolic_data*/,
             typename CellClass::local_expansion_t** const child_local_exps,
             const Symb* const /*child_symbolic_data*/[])
    {
        // Each child is impacted by the father
        for(int idx = 0 ; idx < 8 ; ++idx){
            if(child_local_exps[idx]){
                child_local_exps[idx]->set(node_local_exp->get()
                                           + child_local_exps[idx]->get());
            }
        }

    }

    /** After Downward */
    template<class Symb>
    void L2P(const typename CellClass::local_expansion_t* const leaf_local_exp,
             const Symb* const /* target_symbolic_data */,
             ContainerClass* const particles)
    {
        // The particles is impacted by the parent cell
        long long int*const particlesAttributes = particles->getDataDown();
        for(FSize idxPart = 0 ; idxPart < particles->getNbParticles() ; ++idxPart){
            particlesAttributes[idxPart] += leaf_local_exp->get();
        }
    }


    /** After Downward */
    void P2P(const FTreeCoordinate& ,
                 ContainerClass* const /*FRestrict*/ targets,
             const ContainerClass* const  /*FRestrict*/ sources,
                 ContainerClass* const directNeighborsParticles[], const int /*positions*/[], const int inSize) override {
        // Each particles targeted is impacted by the particles sources
        long long int inc = sources->getNbParticles();
        if(targets == sources){
            inc -= 1;
        }
        for(int idx = 0 ; idx < inSize ; ++idx){
                inc += directNeighborsParticles[idx]->getNbParticles();
        }

        long long int*const particlesAttributes = targets->getDataDown();
        for(FSize idxPart = 0 ; idxPart < targets->getNbParticles() ; ++idxPart){
            particlesAttributes[idxPart] += inc;
        }
    }

    void P2POuter(const FTreeCoordinate& /*inLeafPosition*/,
             ContainerClass* const FRestrict targets,
             ContainerClass* const directNeighborsParticles[], const int neighborPositions[],
             const int inSize) override {
        long long int inc = 0;

        for(int idx = 0 ; idx < inSize ; ++idx){
                inc += directNeighborsParticles[idx]->getNbParticles();
        }

        long long int*const particlesAttributes = targets->getDataDown();
        for(FSize idxPart = 0 ; idxPart < targets->getNbParticles() ; ++idxPart){
            particlesAttributes[idxPart] += inc;
        }
    }

    /** After Downward */
    void P2PRemote(const FTreeCoordinate& ,
                   ContainerClass* const FRestrict targets, const ContainerClass* const FRestrict /*sources*/,
                 const ContainerClass* const directNeighborsParticles[], const int /*positions*/[], const int inSize) override  {
        // Each particles targeted is impacted by the particles sources
        long long int inc = 0;
        for(int idx = 0 ; idx < inSize ; ++idx){
                inc += directNeighborsParticles[idx]->getNbParticles();
        }

        long long int*const particlesAttributes = targets->getDataDown();
        for(FSize idxPart = 0 ; idxPart < targets->getNbParticles() ; ++idxPart){
            particlesAttributes[idxPart] += inc;
        }
    }
};


/** This function test the octree to be sure that the fmm algorithm
  * has worked completly.
  */
template< class OctreeClass, class CellClass, class ContainerClass, class LeafClass>
void ValidateFMMAlgo(OctreeClass* const tree){
    std::cout << "Check Result\n";
    const int TreeHeight = tree->getHeight();
    long long int NbPart = 0;
    { // Check that each particle has been summed with all other
        tree->forEachCellLeaf([&](CellClass* cell, LeafClass* leaf){
            if(cell->getMultipoleData().get() != leaf->getSrc()->getNbParticles() ){
                    std::cout << "Problem P2M : " << cell->getMultipoleData().get() <<
                                 " (should be " << leaf->getSrc()->getNbParticles() << ")\n";
            }
            NbPart += leaf->getSrc()->getNbParticles();
        });
    }
    { // Ceck if there is number of NbPart summed at level 1
        typename OctreeClass::Iterator octreeIterator(tree);
        octreeIterator.moveDown();
        long long int res = 0;
        do{
            res += octreeIterator.getCurrentCell()->getMultipoleData().get();
        } while(octreeIterator.moveRight());
        if(res != NbPart){
            std::cout << "Problem M2M at level 1 : " << res << "\n";
        }
    }
    { // Ceck if there is number of NbPart summed at level 1
        typename OctreeClass::Iterator octreeIterator(tree);
        octreeIterator.gotoBottomLeft();
        for(int idxLevel = TreeHeight - 1 ; idxLevel > 1 ; --idxLevel ){
            long long int res = 0;
            do{
                res += octreeIterator.getCurrentCell()->getMultipoleData().get();
            } while(octreeIterator.moveRight());
            if(res != NbPart){
                std::cout << "Problem M2M at level " << idxLevel << " : " << res << "\n";
            }
            octreeIterator.moveUp();
            octreeIterator.gotoLeft();
        }
    }
    { // Check that each particle has been summed with all other
        typename OctreeClass::Iterator octreeIterator(tree);
        octreeIterator.gotoBottomLeft();
        do{
            const bool isUsingTsm = (octreeIterator.getCurrentListTargets() != octreeIterator.getCurrentListSrc());
            const long long int* dataDown = octreeIterator.getCurrentListTargets()->getDataDown();
            for(FSize idxPart = 0 ; idxPart < octreeIterator.getCurrentListTargets()->getNbParticles() ; ++idxPart){
                if( (!isUsingTsm && dataDown[idxPart] != NbPart - 1) ||
                    (isUsingTsm && dataDown[idxPart] != NbPart) ){
                    std::cout << "Problem L2P + P2P : " << dataDown[idxPart] << ", " <<
                                 " NbPart : " << NbPart << ", " <<
                                 " ( Index " << octreeIterator.getCurrentGlobalIndex() << ")\n";
                }
            }
        } while(octreeIterator.moveRight());
    }

    std::cout << "Done\n";
}



#endif //FTESTKERNELS_HPP
