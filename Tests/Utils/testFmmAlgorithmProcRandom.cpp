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

//#include "Core/FFmmAlgorithmProcMpi.hpp"
#include "Core/FFmmAlgorithmThreadProc.hpp"
#include "Core/FFmmAlgorithmThread.hpp"

#include "Files/FRandomLoader.hpp"
#include "Files/FMpiTreeBuilder.hpp"

#include "Components/FBasicKernels.hpp"

#include "Utils/FLeafBalance.hpp"

#include "Utils/FParameterNames.hpp"

#include <iostream>
#include <cstdio>
#include <cstdlib>




/////////////////////////////////////////////////////////////////////
// Define the classes to use
/////////////////////////////////////////////////////////////////////

typedef double FReal;

typedef FTestCell                  CellClass;
typedef FTestParticleContainer<FReal>     ContainerClass;

typedef FSimpleLeaf<FReal, ContainerClass >                     LeafClass;
typedef FOctree<FReal, CellClass, ContainerClass , LeafClass >  OctreeClass;
typedef FTestKernels< CellClass, ContainerClass >         KernelClass;

typedef FFmmAlgorithmThread<OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass >     FmmClass;
typedef FFmmAlgorithmThreadProc<OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass >     FmmClassProc;

/////////////////////////////////////////////////////////////////////
// Main
/////////////////////////////////////////////////////////////////////

// Simply create particles and try the kernels
int main(int argc, char ** argv){
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

    const FSize NbParticles = FParameters::getValue(argc,argv,FParameterDefinitions::NbParticles.options, FSize(10000));

    FRandomLoader<FReal> loader(NbParticles, 1.0, FPoint<FReal>(0,0,0), app.global().processId());
    if(!loader.isOpen()) throw std::runtime_error("Particle file couldn't be opened!");

    const FReal boxWidth = loader.getBoxWidth();
    const FPoint<FReal> centerOfBox = loader.getCenterOfBox();

    std::cout << "Simulation properties :\n";
    std::cout << "Nb Particles For me " << NbParticles << "\n";
    std::cout << "Box Width : " << boxWidth << "\n";
    std::cout << "Box Center : " << centerOfBox << "\n";

    // The real tree to work on
    OctreeClass realTree(NbLevels, SizeSubLevels,boxWidth,centerOfBox);
    {
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
        };

        TestParticle* particles = new TestParticle[NbParticles];
        memset(particles, 0, sizeof(TestParticle) * NbParticles);
        for(FSize idxPart = 0 ; idxPart < NbParticles ; ++idxPart){
            FPoint<FReal> position;
            loader.fillParticle(&position);
            particles[idxPart].position = position;
        }
        std::cout << "Creating  " << "(" << counter.tacAndElapsed() << "s)." << std::endl;
        counter.tic();

        FVector<TestParticle> finalParticles;
        FLeafBalance balancer;
        FMpiTreeBuilder< FReal,TestParticle >::DistributeArrayToContainer(app.global(),particles,
                                                                          NbParticles,
                                                                          realTree.getBoxCenter(),
                                                                          realTree.getBoxWidth(),realTree.getHeight(),
                                                                          &finalParticles, &balancer);
        std::cout << "Sorting  " << "(" << counter.tacAndElapsed() << "s)." << std::endl;
        std::cout << "I have now " << finalParticles.getSize() << " particles\n";
        counter.tic();

        for(int idx = 0 ; idx < finalParticles.getSize(); ++idx){
            realTree.insert(finalParticles[idx].position);
        }

        delete[] particles;

        std::cout << "Inserting  " << "(" << counter.tacAndElapsed() << "s)." << std::endl;
    }
    //////////////////////////////////////////////////////////////////////////////////


    //////////////////////////////////////////////////////////////////////////////////
    // Check particles in tree
    //////////////////////////////////////////////////////////////////////////////////

    std::cout << "Working parallel particles ..." << std::endl;
    counter.tic();

    KernelClass kernels;

    FmmClassProc algo(app.global(),&realTree,&kernels);
    algo.execute();

    counter.tac();
    std::cout << "Done  " << "(@Algorithm Particles = " << counter.elapsed() << "s)." << std::endl;

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////


    {
        std::cout << "Testing data ..." << std::endl;
        const FSize totalNbParticles = app.global().allReduceSum(NbParticles);
        std::cout << "totalNbParticles " << totalNbParticles << std::endl;

        counter.tic();

        bool hasWorked = true;

        // Check that each particle has been summed with all other
        typename OctreeClass::Iterator octreeIterator(&realTree);
        octreeIterator.gotoBottomLeft();

        do {
            ContainerClass* container = (octreeIterator.getCurrentListTargets());
            const long long int*const dataDown = container->getDataDown();

            for(FSize idxPart = 0 ; idxPart < container->getNbParticles() ; ++idxPart){
                // If a particles has been impacted by less than NbPart - 1 (the current particle)
                // there is a problem
                if( dataDown[idxPart] != totalNbParticles-1){
                    hasWorked = false;
                }
            }

        }while( octreeIterator.moveRight() );

        counter.tac();
        std::cout << "Done  " << "(" << counter.elapsed() << "s)." << std::endl;
        std::cout << "Res  " << (hasWorked?"Worked":"Failed") << std::endl;
    }

    return 0;
}
