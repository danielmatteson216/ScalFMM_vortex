// See LICENCE file at project root

// ==== CMAKE =====
// @FUSE_MPI
// ================

#include "ScalFmmConfig.h"
#include <cstdlib>
#include <string.h>
#include <stdexcept>
#include <algorithm>
#include <vector>

#include "FUTester.hpp"

#include "Utils/FMpi.hpp"
#include "Containers/FVector.hpp"

#include "Files/FRandomLoader.hpp"
#include "Utils/FLeafBalance.hpp"
#include "Containers/FTreeCoordinate.hpp"
#include "Containers/FCoordinateComputer.hpp"

#include "Utils/FQuickSortMpi.hpp"
#include "Utils/FBitonicSort.hpp"
#include "Files/FMpiStaticTreeBuilder.hpp"
#include "Core/FCoreCommon.hpp"

#include "Utils/FPoint.hpp"
#include "Utils/FMath.hpp"


class TestMpiTreeBuilder :  public FUTesterMpi< class TestMpiTreeBuilder> {

    template <class FReal>
    struct TestParticle{
        FSize indexInFile;
        FPoint<FReal> position;
        FReal physicalValue;

        const FPoint<FReal>& getPosition()const{
            return position;
        }

        bool operator==(const TestParticle& other){
            return indexInFile == other.indexInFile
                    && position.getX() == other.position.getX()
                    && position.getY() == other.position.getY()
                    && position.getZ() == other.position.getZ()
                    && physicalValue == other.physicalValue;
        }
    };

    template <class FReal, int localNbParticlesLocal>
    void RunTest(){
        const int totalNbParticles = localNbParticlesLocal*app.global().processCount();
        std::unique_ptr<TestParticle<FReal>[]> globalParticles(new TestParticle<FReal>[totalNbParticles]);

        FRandomLoader<FReal> loader(totalNbParticles, 1., FPoint<FReal>(0,0,0), 0);

        for(int idxPart = 0 ; idxPart < totalNbParticles ; ++idxPart){
            loader.fillParticle(&globalParticles[idxPart].position);
            globalParticles[idxPart].physicalValue = FReal(idxPart);
            globalParticles[idxPart].indexInFile = idxPart;
        }

        for(int treeHeigt = 3 ; treeHeigt < 8 ; ++treeHeigt){
            FVector<TestParticle<FReal>> finalParticles;
            FMpiStaticTreeBuilder< FReal,TestParticle<FReal> >::DistributeArrayToContainer(app.global(),
                                                                    &globalParticles[localNbParticlesLocal*app.global().processId()],
                                                                    localNbParticlesLocal,
                                                                    loader.getCenterOfBox(),
                                                                    loader.getBoxWidth(),treeHeigt,
                                                                    &finalParticles);

            {
                assert(finalParticles.getSize() <= std::numeric_limits<int>::max());
                int myNbParticles = int(finalParticles.getSize());
                int nbParticlesAfter = 0;
                FMpi::MpiAssert(MPI_Allreduce(&myNbParticles,
                                &nbParticlesAfter,
                                1,
                                MPI_INT,
                                MPI_SUM,
                                app.global().getComm()), __LINE__);
                uassert(nbParticlesAfter == totalNbParticles);
            }

            std::unique_ptr<int[]> exists(new int[totalNbParticles]());

            for(int idxPart = 0 ;idxPart < finalParticles.getSize() ; ++idxPart){
                uassert(finalParticles[idxPart] == globalParticles[finalParticles[idxPart].indexInFile]);
                uassert(exists[finalParticles[idxPart].indexInFile] == 0);
                exists[finalParticles[idxPart].indexInFile] += 1;
            }

            std::unique_ptr<int[]> existsAll(new int[totalNbParticles]());

            FMpi::MpiAssert(MPI_Allreduce(exists.get(),
                            existsAll.get(),
                            totalNbParticles,
                            MPI_INT,
                            MPI_SUM,
                            app.global().getComm()), __LINE__);


            for(int idxPart = 0 ; idxPart < totalNbParticles ; ++idxPart){
                uassert(existsAll[idxPart] == 1);
            }
        }
    }

    /** If memstas is running print the memory used */
    void PostTest() {
        if( FMemStats::controler.isUsed() ){
            std::cout << app.global().processId() << "-> Memory used at the end " << FMemStats::controler.getCurrentAllocated()
                      << " Bytes (" << FMemStats::controler.getCurrentAllocatedMB() << "MB)\n";
            std::cout << app.global().processId() << "-> Max memory used " << FMemStats::controler.getMaxAllocated()
                      << " Bytes (" << FMemStats::controler.getMaxAllocatedMB() << "MB)\n";
            std::cout << app.global().processId() << "-> Total memory used " << FMemStats::controler.getTotalAllocated()
                      << " Bytes (" << FMemStats::controler.getTotalAllocatedMB() << "MB)\n";
        }
    }

    void SetTests(){
        AddTest(&TestMpiTreeBuilder::RunTest<double, 1000>,"Generate particles and distribute them");
        AddTest(&TestMpiTreeBuilder::RunTest<float, 1000>,"Generate particles and distribute them");

        AddTest(&TestMpiTreeBuilder::RunTest<double, 1200>,"Generate particles and distribute them");
        AddTest(&TestMpiTreeBuilder::RunTest<float, 1200>,"Generate particles and distribute them");
    }

public:
    TestMpiTreeBuilder(int argc,char ** argv) : FUTesterMpi(argc,argv){
    }


};

TestClassMpi(TestMpiTreeBuilder);
