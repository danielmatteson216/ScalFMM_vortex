// ===================================================================================// See LICENCE file at project root

#include <iostream>
#include <iomanip>

#include <cstdio>  //printf
#include <cstdlib>
#include <cstring>  //memset

#include <cmath>
#include <algorithm>
#include <string>

#include  "ScalFmmConfig.h"
#include "Utils/FTic.hpp"
#include "Utils/FMath.hpp"
#include "Utils/FParameters.hpp"
#include "Utils/FParameterNames.hpp"
#include "Files/FIOVtk.hpp"

#include "Containers/FOctree.hpp"
#include "Files/FRandomLoader.hpp"

#include "Kernels/P2P/FP2PR.hpp"
#include "Kernels/P2P/FP2PParticleContainerIndexed.hpp"

// Simply create particles and try the kernels
int main(int argc, char ** argv){
    FHelpDescribeAndExit(argc, argv,
                         ">> This executable test the efficiency of the computation of the P2P",
                         FParameterDefinitions::NbParticles);

    const FSize nbParticles = FParameters::getValue(argc, argv, FParameterDefinitions::NbParticles.options, 1000);
    std::cout << "Test with " << nbParticles << " particles." << std::endl;

    //////////////////////////////////////////////////////////

    typedef double FReal;

    FRandomLoader<FReal> loader(nbParticles*2);

    FTic timer;
    FP2PParticleContainer<FReal> leaf1;
    for(FSize idxPart = 0 ; idxPart < nbParticles ; ++idxPart){
        FPoint<FReal> pos;
        loader.fillParticle(&pos);
        leaf1.push(pos, 1.0);
    }

    FP2PParticleContainer<FReal> leaf2;
    for(FSize idxPart = 0 ; idxPart < nbParticles ; ++idxPart){
        FPoint<FReal> pos;
        loader.fillParticle(&pos);
        leaf2.push(pos, 1.0);
    }
    FP2PParticleContainer<FReal> * const pleaf2 = &leaf2;
    std::cout << "Timer taken to create and insert the particles = " << timer.tacAndElapsed() << "s" << std::endl;

    //////////////////////////////////////////////////////////

    std::cout << "Double pricision:" <<  std::endl;

    timer.tic();
    FP2PRT<double>::FullMutual<FP2PParticleContainer<FReal>>( &leaf1, &pleaf2, 1);
    timer.tac();
    std::cout << "Timer taken by FullMutual = " << timer.elapsed() << "s" << std::endl;

    timer.tic();
    FP2PRT<double>::FullRemote<FP2PParticleContainer<FReal>>( &leaf1, &pleaf2, 1);
    timer.tac();
    std::cout << "Timer taken by FullRemote = " << timer.elapsed() << "s" << std::endl;

    //////////////////////////////////////////////////////////

    return 0;
}
