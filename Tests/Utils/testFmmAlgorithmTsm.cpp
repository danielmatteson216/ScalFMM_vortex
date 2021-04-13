// See LICENCE file at project root

#include <iostream>

#include <cstdio>
#include <cstdlib>

#include "Utils/FTic.hpp"

#include "Containers/FOctree.hpp"
#include "Containers/FVector.hpp"
#include "Utils/FParameters.hpp"

#include "Components/FTypedLeaf.hpp"

#include "Utils/FPoint.hpp"

#include "Components/FTestCell.hpp"
#include "Components/FTestKernels.hpp"

#include "Extensions/FExtendCellType.hpp"

#include "Core/FFmmAlgorithmTsm.hpp"
#include "Core/FFmmAlgorithmThreadTsm.hpp"

#include "Components/FBasicKernels.hpp"

#include "Files/FRandomLoader.hpp"

#include "Components/FTestParticleContainer.hpp"

#include "Utils/FParameterNames.hpp"

/** This program show an example of use of
  * the fmm basic algo
  * it also check that each particles is impacted each other particles
  */
class FTestCellTsm: public FTestCell , public FExtendCellType{
};

// Simply create particles and try the kernels
int main(int argc, char ** argv){
    FHelpDescribeAndExit(argc, argv,
                         "Test FMM TSM (target source model) algorithm by counting the nb of interactions each particle receive.",
                         FParameterDefinitions::OctreeHeight, FParameterDefinitions::OctreeSubHeight,
                         FParameterDefinitions::NbParticles);

    typedef double FReal;
    typedef FTestCellTsm                 CellClassTyped;
    typedef FTestParticleContainer<FReal>       ContainerClassTyped;

    typedef FTypedLeaf< FReal, ContainerClassTyped >                      LeafClassTyped;
    typedef FOctree<FReal, CellClassTyped, ContainerClassTyped , LeafClassTyped >  OctreeClassTyped;
    typedef FTestKernels< CellClassTyped, ContainerClassTyped >          KernelClassTyped;

    typedef FFmmAlgorithmThreadTsm<OctreeClassTyped, CellClassTyped, ContainerClassTyped, KernelClassTyped, LeafClassTyped > FmmClassTyped;
    ///////////////////////What we do/////////////////////////////
    std::cout << ">> This executable has to be used to test the FMM algorithm.\n";
    //////////////////////////////////////////////////////////////

    const int NbLevels = FParameters::getValue(argc,argv,FParameterDefinitions::OctreeHeight.options, 5);
    const int SizeSubLevels = FParameters::getValue(argc,argv,FParameterDefinitions::OctreeSubHeight.options, 3);
    const FSize NbPart = FParameters::getValue(argc,argv,FParameterDefinitions::NbParticles.options, FSize(2000000));
    FTic counter;

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    FRandomLoaderTsm<FReal> loader(NbPart, 1, FPoint<FReal>(0.5,0.5,0.5), 1);
    OctreeClassTyped tree(NbLevels, SizeSubLevels,loader.getBoxWidth(),loader.getCenterOfBox());

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    std::cout << "Creating " << NbPart << " particles ..." << std::endl;
    counter.tic();

    {
        FPoint<FReal> particlePosition;
        FParticleType particleType;
        for(FSize idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
            loader.fillParticle(&particlePosition, &particleType);
            tree.insert(particlePosition, particleType);
        }
    }

    counter.tac();
    std::cout << "Done  " << "(@Creating and Inserting Particles = " << counter.elapsed() << "s)." << std::endl;


    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    std::cout << "Working on particles ..." << std::endl;
    counter.tic();

    KernelClassTyped kernels;

    FmmClassTyped algo(&tree,&kernels);
    algo.execute();

    counter.tac();
    std::cout << "Done  " << "(@Algorithm = " << counter.elapsed() << "s)." << std::endl;

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    ValidateFMMAlgo<OctreeClassTyped, CellClassTyped, ContainerClassTyped, LeafClassTyped>(&tree);

    return 0;
}



