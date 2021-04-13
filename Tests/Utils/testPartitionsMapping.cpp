// See LICENCE file at project root

// ==== CMAKE =====
// @FUSE_MPI
// ================

#include <iostream>

#include <cstdio>
#include <cstdlib>


#include "Kernels/Rotation/FRotationCell.hpp"
#include "Kernels/Rotation/FRotationKernel.hpp"

#include "Components/FSimpleLeaf.hpp"
#include "Kernels/P2P/FP2PParticleContainerIndexed.hpp"

#include "Utils/FParameters.hpp"
#include "Utils/FMemUtils.hpp"
#include "Utils/FPartitionsMapping.hpp"

#include "Containers/FOctree.hpp"
#include "Containers/FVector.hpp"

#include "Files/FRandomLoader.hpp"
#include "Files/FMpiTreeBuilder.hpp"

#include "Core/FFmmAlgorithm.hpp"
#include "Core/FFmmAlgorithmThread.hpp"
#include "Core/FFmmAlgorithmThreadProc.hpp"

#include "Utils/FLeafBalance.hpp"

#include "Utils/FParameterNames.hpp"

/**
 * This program runs the FMM Algorithm Distributed with the Rotation kernel
 */

// Simply create particles and try the kernels
int main(int argc, char* argv[])
{
    FHelpDescribeAndExit(argc, argv,
                         "Test with MPI the chebyshev FMM and compare it to the direct computation for debugging purpose.",
                         FParameterDefinitions::NbParticles, FParameterDefinitions::OctreeHeight,
                         FParameterDefinitions::OctreeSubHeight, FParameterDefinitions::NbThreads);

    typedef double FReal;

    FMpi app(argc,argv);

    const FSize nbParticles       = FParameters::getValue(argc,argv, FParameterDefinitions::NbParticles.options, 10000000ULL);
    const unsigned int TreeHeight    = FParameters::getValue(argc, argv, FParameterDefinitions::OctreeHeight.options, 5);
    FTic time;

    std::cout << ">> This executable has to be used to test Proc Rotation Algorithm. \n";

    // init particles position and physical value
    struct TestParticle{
        FPoint<FReal> position;
        const FPoint<FReal>& getPosition(){
            return position;
        }
    };

    // open particle file
    std::cout << "Creating : " << nbParticles << "\n" << std::endl;
    FRandomLoader<FReal> loader(nbParticles, 1.0, FPoint<FReal>(0,0,0), app.global().processId());

    time.tic();
    std::unique_ptr<TestParticle[]> particles(new TestParticle[loader.getNumberOfParticles()]);
    for(FSize idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
        loader.fillParticle(&particles[idxPart].position);
    }

    FPartitionsMapping<FReal> map(app.global());

    FVector<FPartitionsMapping<double>::TestParticle<0> > finalParticles = map.distributeParticles<0>(loader.getNumberOfParticles(),
                                                                   loader.getCenterOfBox(),
                                                                   loader.getBoxWidth(), TreeHeight,
                                                                   [&](const int idx, FPoint<FReal>* position, std::array<FReal, 0>* /*val*/){
        position->setPosition(particles[idx].position.getX(),
                              particles[idx].position.getY(),
                              particles[idx].position.getZ());
    });

    // Test every particles exist
    {
        std::unique_ptr<int[]> count(new int[app.global().processCount() * nbParticles]);
        memset(count.get(), 0, sizeof(int) * app.global().processCount() * nbParticles);
        for(FSize part = 0 ; part < finalParticles.getSize() ; ++part){
            const FSize idx = finalParticles[part].initialProcOwner*nbParticles + finalParticles[part].localIndex;
            FAssertLF(count[idx] == 0)
            count[idx] = 1;
        }

        FMpi::Assert( MPI_Allreduce(MPI_IN_PLACE, count.get(), 
				    app.global().processCount()*static_cast<int>(nbParticles), 
				    MPI_INT, MPI_SUM,
				    app.global().getComm()), __LINE__);

        for(FSize part = 0 ; part < app.global().processCount()*nbParticles ; ++part){
            FAssertLF(count[part] == 1);
        }
    }

    // Test to send data
    {
        std::unique_ptr<std::array<FReal, 2>[]> toDistr(new std::array<FReal, 2>[nbParticles]);
        for(FSize idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
            toDistr[idxPart][0] = app.global().processId();
            toDistr[idxPart][1] = FReal(idxPart);
        }

        std::unique_ptr<std::array<FReal, 2>[]> res = map.distributeData<2>(toDistr.get());

        for(FSize part = 0 ; part < finalParticles.getSize() ; ++part){
            if(int(res[part][0]) != finalParticles[part].initialProcOwner)
                    std::cout << "[" << app.global().processId() << "]Res proc is " << res[part][0] << " should be " << finalParticles[part].initialProcOwner << std::endl;
            if(int(res[part][1]) != finalParticles[part].localIndex)
                    std::cout << "[" << app.global().processId() << "]Res localidx is " << res[part][1] << " should be " << finalParticles[part].localIndex << std::endl;
        }

        std::unique_ptr<std::array<FReal, 2>[]> resBack = map.getResultingData<2>(res.get());

        for(FSize part = 0 ; part < loader.getNumberOfParticles() ; ++part){
            if(int(resBack[part][0]) != app.global().processId())
                    std::cout << "[" << app.global().processId() << "]ResBack proc is " << resBack[part][0] << " should be " << app.global().processId() << std::endl;
            if(int(resBack[part][1]) != part)
                    std::cout << "[" << app.global().processId() << "]ResBack localidx is " << resBack[part][1] << " should be " << part << std::endl;
        }
    }

    return 0;
}
