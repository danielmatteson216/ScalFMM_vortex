// See LICENCE file at project root

// ==== CMAKE =====
// @FUSE_MPI
// ================

#include "Utils/FMpi.hpp"
#include "Utils/FTic.hpp"

#include "Containers/FOctree.hpp"
#include "Containers/FVector.hpp"
#include "Utils/FParameters.hpp"
#include "Utils/FGlobal.hpp"

#include "Components/FSimpleLeaf.hpp"

#include "Utils/FPoint.hpp"

#include "Components/FTestCell.hpp"
#include "Components/FTestKernels.hpp"
#include "Components/FTestParticleContainer.hpp"

#include "Core/FFmmAlgorithmThreadProcTsm.hpp"
#include "Core/FFmmAlgorithmThreadTsm.hpp"

#include "Files/FRandomLoader.hpp"
#include "Files/FMpiTreeBuilder.hpp"

#include "Components/FBasicKernels.hpp"

#include "Utils/FLeafBalance.hpp"

#include "Utils/FParameterNames.hpp"
#include "Extensions/FExtendCellType.hpp"

#include <iostream>
#include <cstdio>
#include <cstdlib>


/** This program show an example of use of the fmm threaded + mpi algo
 * it also check that each particles is impacted each other particles
 */

/////////////////////////////////////////////////////////////////////////////
// Test function
/////////////////////////////////////////////////////////////////////////////

// Check if tree is built correctly
template<class OctreeClass>
void ValidateTree(OctreeClass& realTree,
                  OctreeClass& treeValide){

    typename OctreeClass::Iterator octreeIteratorValide(&treeValide);
    octreeIteratorValide.gotoBottomLeft();

    typename OctreeClass::Iterator octreeIterator(&realTree);
    octreeIterator.gotoBottomLeft();

    while(octreeIteratorValide.getCurrentGlobalIndex() != octreeIterator.getCurrentGlobalIndex()){
        if(octreeIteratorValide.moveRight() == false){
            std::cout << "Error the real tree smaller than the parallel one\n";
            std::cout << "Valide tree stop at " << octreeIteratorValide.getCurrentGlobalIndex() << "\n";
            std::cout << "Other at " << octreeIterator.getCurrentGlobalIndex() << "\n";
            return;
        }
    }

    std::cout << "The tree starts at " << octreeIteratorValide.getCurrentGlobalIndex() << "\n";

    while(true){
        if(octreeIteratorValide.getCurrentGlobalIndex() != octreeIterator.getCurrentGlobalIndex()){
            std::cout << "Error the trees do not have the same idx.\n";
            std::cout << "Correct one is " << octreeIteratorValide.getCurrentGlobalIndex() << "\n";
            std::cout << "Incorrect one is " << octreeIterator.getCurrentGlobalIndex() << "\n";
            return;
        }

        if(octreeIteratorValide.getCurrentListSrc()->getNbParticles() != octreeIterator.getCurrentListSrc()->getNbParticles()){
            std::cout << "Error the trees do not have the same number of particles at idx " << octreeIteratorValide.getCurrentGlobalIndex() << ".\n";
            std::cout << "Correct one is " << octreeIteratorValide.getCurrentListSrc()->getNbParticles() << "\n";
            std::cout << "Incorrect one is " << octreeIterator.getCurrentListSrc()->getNbParticles() << "\n";
            return;
        }

        if(octreeIterator.moveRight() == false){
            break;
        }

        if(octreeIteratorValide.moveRight() == false){
            std::cout << "Error the real tree smaller than the parallel one\n";
        }
    }

    std::cout << "The tree stops at " << octreeIteratorValide.getCurrentGlobalIndex() << "\n";
}



/** This function tests the octree to be sure that the fmm algorithm
 * has worked completly.
 */
template<class OctreeClass, class ContainerClass, class FmmClassProc>
void ValidateFMMAlgoProc(OctreeClass* const badTree,
                         OctreeClass* const valideTree,
                         FmmClassProc* const fmm){
    std::cout << "The working interval is from " << fmm->getWorkingInterval(badTree->getHeight()-1).leftIndex << "\n";
    std::cout << "The working interval is to " << fmm->getWorkingInterval(badTree->getHeight()-1).rightIndex << "\n";
    std::cout << "\tValidate L2L M2M...\n";
    const int OctreeHeight = badTree->getHeight();
    {
        typename OctreeClass::Iterator octreeIterator(badTree);
        octreeIterator.gotoBottomLeft();

        typename OctreeClass::Iterator octreeIteratorValide(valideTree);
        octreeIteratorValide.gotoBottomLeft();

        for(int level = OctreeHeight - 1 ; level > 0 && fmm->hasWorkAtLevel(level) ; --level){

            while(octreeIteratorValide.getCurrentGlobalIndex() != octreeIterator.getCurrentGlobalIndex()) {
                octreeIteratorValide.moveRight();
            }

            while(octreeIteratorValide.getCurrentGlobalIndex() != fmm->getWorkingInterval(level).leftIndex){
                octreeIteratorValide.moveRight();
                octreeIterator.moveRight();
            }

            FSize countCheck = 0;
            do{
                if(octreeIterator.getCurrentGlobalIndex() != octreeIteratorValide.getCurrentGlobalIndex()){
                    std::cout << "Problem Error index are not equal! " << octreeIterator.getCurrentGlobalIndex() << " " << octreeIteratorValide.getCurrentGlobalIndex() << std::endl;
                }
                else{
                    if(octreeIterator.getCurrentCell()->getDataUp() != octreeIteratorValide.getCurrentCell()->getDataUp()){
                        std::cout << "Problem M2M error at level " << level << " up bad " << octreeIterator.getCurrentCell()->getDataUp()
                                  << " good " << octreeIteratorValide.getCurrentCell()->getDataUp() << " index " << octreeIterator.getCurrentGlobalIndex() << std::endl;
                    }
                    if(octreeIterator.getCurrentCell()->getDataDown() != octreeIteratorValide.getCurrentCell()->getDataDown()){
                        std::cout << "Problem L2L error at level " << level << " down bad " << octreeIterator.getCurrentCell()->getDataDown()
                                  << " good " << octreeIteratorValide.getCurrentCell()->getDataDown() << " index " << octreeIterator.getCurrentGlobalIndex() << std::endl;
                    }
                }
                ++countCheck;
            } while(octreeIteratorValide.moveRight() && octreeIterator.moveRight());

            // Check that each particle has been summed with all other

            octreeIterator.moveUp();
            octreeIterator.gotoLeft();

            octreeIteratorValide.moveUp();
            octreeIteratorValide.gotoLeft();
        }
    }
    std::cout << "\tValidate L2P P2P...\n";
    {
        FSize NbPart = 0;
        FSize NbLeafs = 0;
        { // Check that each particle has been summed with all other
            typename OctreeClass::Iterator octreeIterator(valideTree);
            octreeIterator.gotoBottomLeft();
            do{
                NbPart += octreeIterator.getCurrentListSrc()->getNbParticles();
                ++NbLeafs;
            } while(octreeIterator.moveRight());
        }
        {
            // Check that each particle has been summed with all other
            typename OctreeClass::Iterator octreeIterator(badTree);
            octreeIterator.gotoBottomLeft();

            do {
                const bool isUsingTsm = (octreeIterator.getCurrentListTargets() != octreeIterator.getCurrentListSrc());

                ContainerClass* container = (octreeIterator.getCurrentListTargets());
                const long long int*const dataDown = container->getDataDown();

                for(FSize idxPart = 0 ; idxPart < container->getNbParticles() ; ++idxPart){
                    // If a particles has been impacted by less than NbPart - 1 (the current particle)
                    // there is a problem
                    if( (!isUsingTsm && dataDown[idxPart] != NbPart - 1) ||
                            (isUsingTsm && dataDown[idxPart] != NbPart) ){
                        std::cout << "Problem L2P + P2P, value on particle is : " << dataDown[idxPart] <<
                                     " at pos " << idxPart << " index is " << octreeIterator.getCurrentGlobalIndex() << "\n";
                    }
                }
            } while( octreeIterator.moveRight());
        }
    }
    std::cout << "\tValidate P2M...\n";
    {
        {
            // Check that each particle has been summed with all other
            typename OctreeClass::Iterator octreeIterator(badTree);
            octreeIterator.gotoBottomLeft();

            do {
                if(octreeIterator.getCurrentListSrc()->getNbParticles() != octreeIterator.getCurrentCell()->getDataUp()){
		  std::cout << "P2M Problem nb part " 
			    << octreeIterator.getCurrentListSrc()->getNbParticles()
			    << " data up "<< octreeIterator.getCurrentCell()->getDataUp()
			    <<std::endl; 
                }
            } while( octreeIterator.moveRight() );
        }
    }
    std::cout << "\tValidate Particles...\n";
    {
        // Check that each particle has been summed with all other
        typename OctreeClass::Iterator octreeIterator(badTree);
        octreeIterator.gotoBottomLeft();

        typename OctreeClass::Iterator valideOctreeIterator(valideTree);
        valideOctreeIterator.gotoBottomLeft();
        while(valideOctreeIterator.getCurrentGlobalIndex() != octreeIterator.getCurrentGlobalIndex()){
            valideOctreeIterator.moveRight();
        }

        do {
            if(valideOctreeIterator.getCurrentGlobalIndex() != octreeIterator.getCurrentGlobalIndex()){
                printf("Problem Do not have the same index valide %lld invalide %lld \n",
                       valideOctreeIterator.getCurrentGlobalIndex(), octreeIterator.getCurrentGlobalIndex());
                break;
            }

            if(octreeIterator.getCurrentListTargets()->getNbParticles() != valideOctreeIterator.getCurrentListTargets()->getNbParticles()){
                printf("Problem Do not have the same number of particle at leaf id %lld, valide %lld invalide %lld \n",
                       octreeIterator.getCurrentGlobalIndex(), valideOctreeIterator.getCurrentListTargets()->getNbParticles(),
                       octreeIterator.getCurrentListTargets()->getNbParticles());
            }
            else {
                ContainerClass* container = (octreeIterator.getCurrentListTargets());
                const long long int*const dataDown = container->getDataDown();

                ContainerClass* containerValide = (valideOctreeIterator.getCurrentListTargets());
                const long long int*const dataDownValide = containerValide->getDataDown();

                for(FSize idxPart = 0 ; idxPart < container->getNbParticles() ; ++idxPart){
                    // If a particles has been impacted by less than NbPart - 1 (the current particle)
                    // there is a problem
                    if( dataDown[idxPart] != dataDownValide[idxPart]){
                        std::cout << "Problem on leaf " << octreeIterator.getCurrentGlobalIndex() <<
                                     " part " << idxPart << " valide data down " << dataDownValide[idxPart] <<
                                     " invalide " << dataDown[idxPart] << "\n";
                        std::cout << "Data down for leaf is: valide " << valideOctreeIterator.getCurrentCell()->getDataDown()
                                  << " invalide " << octreeIterator.getCurrentCell()->getDataDown()
                                  << " size is: valide " <<  valideOctreeIterator.getCurrentListTargets()->getNbParticles()
                                  << " invalide " << octreeIterator.getCurrentListTargets()->getNbParticles() << std::endl;
                    }
                }
            }

        }while( octreeIterator.moveRight() && valideOctreeIterator.moveRight());
    }
    std::cout << "\tDone!\n";
}


/** To print an octree
 * used to debug and understand how the values were passed
 */
template<class OctreeClass>
void print(OctreeClass* const valideTree){
    typename OctreeClass::Iterator octreeIterator(valideTree);
    for(int idxLevel = valideTree->getHeight() - 1 ; idxLevel > 1 ; --idxLevel ){
        do{
            std::cout << "[" << octreeIterator.getCurrentGlobalIndex() << "] up:" << octreeIterator.getCurrentCell()->getDataUp() << " down:" << octreeIterator.getCurrentCell()->getDataDown() << "\t";
        } while(octreeIterator.moveRight());
        std::cout << "\n";
        octreeIterator.gotoLeft();
        octreeIterator.moveDown();
    }
}



/////////////////////////////////////////////////////////////////////
// Main
/////////////////////////////////////////////////////////////////////

class FTestCellTsm: public FTestCell , public FExtendCellType{
};

// Simply create particles and try the kernels
int main(int argc, char ** argv){
    /////////////////////////////////////////////////////////////////////
    // Define the classes to use
    /////////////////////////////////////////////////////////////////////

    typedef double FReal;

    typedef FTestCellTsm                  CellClass;
    typedef FTestParticleContainer<FReal>     ContainerClass;

    typedef FSimpleLeaf<FReal, ContainerClass >                     LeafClass;
    typedef FOctree<FReal, CellClass, ContainerClass , LeafClass >  OctreeClass;
    typedef FTestKernels< CellClass, ContainerClass >         KernelClass;

    typedef FFmmAlgorithmThreadTsm<OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass >     FmmClass;
    typedef FFmmAlgorithmThreadProcTsm<OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass >     FmmClassProc;

    FHelpDescribeAndExit(argc, argv,
                         "Test FMM distributed algorithm by counting the nb of interactions each particle receive.",
                         FParameterDefinitions::OctreeHeight, FParameterDefinitions::OctreeSubHeight,
                         FParameterDefinitions::InputFile);
    ///////////////////////What we do/////////////////////////////
    std::cout << ">> This executable has to be used to test the FMM algorithm.\n";
    //////////////////////////////////////////////////////////////

    FMpi app( argc, argv);

    const int NbLevels = FParameters::getValue(argc,argv,FParameterDefinitions::OctreeHeight.options, 5);
    const int SizeSubLevels = FParameters::getValue(argc,argv,FParameterDefinitions::OctreeSubHeight.options, 3);
    FTic counter;
    const FSize nbParticles = 100000;
    FRandomLoader<FReal> loader(nbParticles, 1.0, FPoint<FReal>(0,0,0), nbParticles);
    if(!loader.isOpen()) throw std::runtime_error("Particle file couldn't be opened!");

    std::vector<FPoint<FReal>> allParticles(nbParticles);
    for(int idxPart = 0 ; idxPart < nbParticles ; ++idxPart){
        loader.fillParticle(&allParticles[idxPart]);
    }

    const double chunckSize = static_cast<double>(nbParticles+app.global().processCount()-1)/double(app.global().processCount());
    const FSize offsetParticles = static_cast<FSize>(chunckSize*static_cast<double>(app.global().processId()));
    const FSize nbPartsForMe = std::min(nbParticles, FSize(chunckSize*double(app.global().processId()+1)))
                                - offsetParticles;
    const FReal boxWidth = loader.getBoxWidth();
    const FPoint<FReal> centerOfBox = loader.getCenterOfBox();

    std::cout << "Simulation properties :\n";
    std::cout << "Nb Particles For me " << nbPartsForMe << "\n";
    std::cout << "Box Width : " << boxWidth << "\n";
    std::cout << "Box Center : " << centerOfBox << "\n";
    std::cout << "offsetParticles " << offsetParticles << ":\n";
    std::cout << "nbPartsForMe " << nbPartsForMe << ":\n";

    // The real tree to work on
    OctreeClass realTree(NbLevels, SizeSubLevels,boxWidth,centerOfBox);

    FSize totalNbSources = 0;
    FSize totalNbTargets = 0;

    if( app.global().processCount() != 1){
        //////////////////////////////////////////////////////////////////////////////////
        // Build tree from mpi loader
        //////////////////////////////////////////////////////////////////////////////////
        std::cout << "Build Tree ..." << std::endl;
        counter.tic();

        struct TestParticle{
            FPoint<FReal> position;
            const FPoint<FReal>& getPosition(){
                return position;
            }
            bool isSource;
        };

        TestParticle* particles = new TestParticle[nbPartsForMe];
        memset(particles, 0, sizeof(TestParticle) * nbPartsForMe);
        FSize nbSources = 0;
        FSize nbTargets = 0;
        for(FSize idxPart = 0 ; idxPart < nbPartsForMe ; ++idxPart){
            FPoint<FReal> position = allParticles[idxPart + offsetParticles];
            particles[idxPart].position = position;
            particles[idxPart].isSource = ((idxPart + offsetParticles)&1);
            if(particles[idxPart].isSource){
                nbSources += 1;
            }
            else{
                nbTargets += 1;
            }
        }


        std::cout << "I init " << nbPartsForMe << " particles\n";
        std::cout << "I init " << nbSources << " sources\n";
        std::cout << "I init " << nbTargets << " sources\n";

        totalNbSources = app.global().allReduceSum(nbSources);
        totalNbTargets = app.global().allReduceSum(nbTargets);

        std::cout << "totalNbSources " << totalNbSources << " sources\n";
        std::cout << "totalNbTargets " << totalNbTargets << " sources\n";

        FVector<TestParticle> finalParticles;
        FLeafBalance balancer;
        FMpiTreeBuilder< FReal,TestParticle >::DistributeArrayToContainer(app.global(),particles,
                                                                    nbPartsForMe,
                                                                    realTree.getBoxCenter(),
                                                                    realTree.getBoxWidth(),realTree.getHeight(),
                                                                    &finalParticles, &balancer);
        std::cout << "I have now " << finalParticles.getSize() << " particles\n";

        for(int idx = 0 ; idx < finalParticles.getSize(); ++idx){
            realTree.insert(finalParticles[idx].position,
                            finalParticles[idx].isSource? FParticleType::source: FParticleType::target);
        }

        delete[] particles;

        counter.tac();
        std::cout << "Done  " << "(" << counter.elapsed() << "s)." << std::endl;

        //////////////////////////////////////////////////////////////////////////////////
    }
    else{
        FAssertLF(offsetParticles == 0);
        for(FSize idxPart = 0 ; idxPart < nbPartsForMe ; ++idxPart){
            const bool isSource = idxPart&1;
            realTree.insert(allParticles[idxPart],
                            isSource? FParticleType::source: FParticleType::target);
            if(isSource){
                totalNbSources += 1;
            }
            else{
                totalNbTargets += 1;
            }
        }
    }
    FAssertLF(totalNbSources + totalNbTargets == loader.getNumberOfParticles(),
              totalNbSources, "+", totalNbTargets, "==", loader.getNumberOfParticles());

    //////////////////////////////////////////////////////////////////////////////////
    // Create real tree
    //////////////////////////////////////////////////////////////////////////////////

    OctreeClass treeValide(NbLevels, SizeSubLevels,boxWidth,centerOfBox);
    {
        for(FSize idxPart = 0 ; idxPart < nbParticles ; ++idxPart){
            bool isSource = idxPart&1;
            treeValide.insert(allParticles[idxPart],
                              isSource? FParticleType::source: FParticleType::target);
        }
    }

    //////////////////////////////////////////////////////////////////////////////////
    // Check particles in tree
    //////////////////////////////////////////////////////////////////////////////////
    std::cout << "Validate tree ..." << std::endl;
    counter.tic();

    ValidateTree(realTree, treeValide);

    counter.tac();
    std::cout << "Done  " << "(" << counter.elapsed() << "s)." << std::endl;

    //////////////////////////////////////////////////////////////////////////////////

    std::cout << "Working parallel particles ..." << std::endl;
    counter.tic();

    KernelClass kernels;

    FmmClassProc algo(app.global(),&realTree,&kernels);
    algo.execute();

    counter.tac();
    std::cout << "Done  " << "(@Algorithm Particles = " << counter.elapsed() << "s)." << std::endl;

    //////////////////////////////////////////////////////////////////////////////////

    std::cout << "Working sequential particles ..." << std::endl;
    counter.tic();

    FmmClass algoValide(&treeValide,&kernels);
    algoValide.execute();

    counter.tac();
    std::cout << "Done  " << "(@Algorithm Particles = " << counter.elapsed() << "s)." << std::endl;

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    std::cout << "Checking data ..." << std::endl;
    counter.tic();

    ValidateFMMAlgoProc<OctreeClass,ContainerClass, FmmClassProc>(&realTree,&treeValide,&algo);

    counter.tac();
    std::cout << "Done  " << "(" << counter.elapsed() << "s)." << std::endl;

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    return 0;
}
