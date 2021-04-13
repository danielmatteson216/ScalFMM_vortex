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

#include "Files/FFmaGenericLoader.hpp"
#include "Files/FMpiFmaGenericLoader.hpp"
#include "Utils/FLeafBalance.hpp"
#include "Containers/FTreeCoordinate.hpp"
#include "Containers/FCoordinateComputer.hpp"

#include "Utils/FQuickSortMpi.hpp"
#include "Utils/FBitonicSort.hpp"
#include "Files/FMpiTreeBuilder.hpp"
#include "Core/FCoreCommon.hpp"

#include "Utils/FPoint.hpp"
#include "Utils/FMath.hpp"


class TestMpiTreeBuilder :  public FUTesterMpi< class TestMpiTreeBuilder> {

    template <class FReal>
    struct TestParticle{
        MortonIndex mindex;
        FSize indexInFile;
        FPoint<FReal> position;
        FReal physicalValue;

        const FPoint<FReal>& getPosition()const{
            return position;
        }

        bool operator<(const TestParticle& rhs)const{
            return rhs.mindex < this->mindex || (rhs.mindex == this->mindex && this->indexInFile < rhs.indexInFile);
        }

    };

    void RunTest(){
        typedef double FReal;
        //
        // Load particles
        //
        if(sizeof(FReal) == sizeof(float) ) {
            std::cerr << "No input data available for Float "<< std::endl;
            exit(EXIT_FAILURE);
        }
        const std::string parFile( (sizeof(FReal) == sizeof(float))?
                                       "Test/DirectFloatbfma":
                                       "UTest/DirectDouble.bfma");
        //Let the choice there to test
        std::string filename(SCALFMMDataPath+parFile);
        //std::string filename("../Data/unitCubeXYZQ100.bfma");

        int TreeHeight = 3;

        FMpiFmaGenericLoader<FReal> loader(filename,app.global());
        if(!loader.isOpen()) throw std::runtime_error("Particle file couldn't be opened!") ;

        //Get the needed informations
        const FReal boxWidth = loader.getBoxWidth();
        const FReal boxWidthAtLeafLevel = boxWidth/FReal(1 << (TreeHeight - 1));

        const FPoint<FReal> centerOfBox = loader.getCenterOfBox();
        const FPoint<FReal> boxCorner   = centerOfBox - boxWidth/2;
        //Now, we sort again the particles with MPI QuickSort
        const FSize idxStart = loader.getStart();
        FAssertLF(idxStart + loader.getMyNumberOfParticles() <= loader.getNumberOfParticles());

        FMpiTreeBuilder<FReal,TestParticle<FReal>>::IndexedParticle * arrayToBeSorted = new FMpiTreeBuilder<FReal,TestParticle<FReal>>::IndexedParticle[loader.getMyNumberOfParticles()];
        //Copy the TestParticles into an array of indexedParticle
        for(FSize i=0 ; i<loader.getMyNumberOfParticles() ; ++i){
            //Fill automatically position AND physicalValue attributes
            loader.fillParticle(&(arrayToBeSorted[i].particle.position),&(arrayToBeSorted[i].particle.physicalValue));

            //We store the index in the file
            arrayToBeSorted[i].particle.indexInFile = i + idxStart;

            //Build temporary TreeCoordinate
            FTreeCoordinate host;
            host.setX( FCoordinateComputer::GetTreeCoordinate<FReal>( arrayToBeSorted[i].particle.getPosition().getX() - boxCorner.getX(), boxWidth, boxWidthAtLeafLevel, TreeHeight ));
            host.setY( FCoordinateComputer::GetTreeCoordinate<FReal>( arrayToBeSorted[i].particle.getPosition().getY() - boxCorner.getY(), boxWidth, boxWidthAtLeafLevel, TreeHeight ));
            host.setZ( FCoordinateComputer::GetTreeCoordinate<FReal>( arrayToBeSorted[i].particle.getPosition().getZ() - boxCorner.getZ(), boxWidth, boxWidthAtLeafLevel, TreeHeight ));

            //Set Morton index from Tree Coordinate
            arrayToBeSorted[i].particle.mindex = host.getMortonIndex();
            arrayToBeSorted[i].index = arrayToBeSorted[i].particle.mindex;
        }
        FMpiTreeBuilder<FReal,TestParticle<FReal>>::IndexedParticle* outputArray = nullptr;
        FSize outputSize;
        FQuickSortMpi<FMpiTreeBuilder<FReal,TestParticle<FReal>>::IndexedParticle,MortonIndex,FSize>::QsMpi(arrayToBeSorted,loader.getMyNumberOfParticles(),&outputArray,&outputSize,app.global());

        if(app.global().processId() == 0){
            FSize allparts;
            MPI_Reduce(&outputSize, &allparts, 1, FMpi::GetType(outputSize), MPI_SUM, 0, app.global().getComm());
            FAssertLF(allparts == loader.getNumberOfParticles());
        }
        else{
            MPI_Reduce(&outputSize, nullptr, 1, FMpi::GetType(outputSize), MPI_SUM, 0, app.global().getComm());
        }

        int* allPartsCount = new int[loader.getNumberOfParticles()];
        memset(allPartsCount, 0, sizeof(int)*loader.getNumberOfParticles());

        for(int idxPart = 0 ; idxPart < outputSize ; ++idxPart){
            FAssertLF(allPartsCount[outputArray[idxPart].particle.indexInFile] == 0);
            allPartsCount[outputArray[idxPart].particle.indexInFile] += 1;
        }

        if(app.global().processId() == 0){
            int* allPartsCountReduced = new int[loader.getNumberOfParticles()];
            MPI_Reduce(allPartsCount, allPartsCountReduced, int(loader.getNumberOfParticles()), FMpi::GetType(*allPartsCount), MPI_SUM, 0, app.global().getComm());

            for(int idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
                FAssertLF(allPartsCountReduced[idxPart] == 1);
            }

            delete[] allPartsCountReduced;
        }
        else{
            MPI_Reduce(allPartsCount, nullptr, int(loader.getNumberOfParticles()), FMpi::GetType(*allPartsCount), MPI_SUM, 0, app.global().getComm());
        }

        MortonIndex leftright[2] = {-1, -1};

        if(outputSize){
            leftright[0] = outputArray[0].particle.mindex;
            leftright[1] = outputArray[0].particle.mindex;

            for(int idxPart = 1 ; idxPart < outputSize ; ++idxPart){
                leftright[0] = std::min(leftright[0], outputArray[idxPart].particle.mindex);
                leftright[1] = std::max(leftright[1], outputArray[idxPart].particle.mindex);
            }
        }

        if(app.global().processId() == 0){
            MortonIndex* allintervals = new MortonIndex[app.global().processCount()*2];
            MPI_Gather(leftright, 2, FMpi::GetType(*leftright),
                       allintervals, 2, FMpi::GetType(*allintervals), 0, app.global().getComm());

            MortonIndex currentLimit = -1;

            for(int idxProc = 0 ; idxProc < app.global().processCount() ; ++idxProc){
                FAssertLF(allintervals[idxProc*2] != -1 || allintervals[idxProc*2+1] == -1);
                if(allintervals[idxProc*2] != -1){
                    FAssertLF(allintervals[idxProc*2] <= allintervals[idxProc*2+1]);

                    if(idxProc && allintervals[idxProc*2-1] != -1){
                        FAssertLF(allintervals[idxProc*2-1] < allintervals[idxProc*2]);
                    }
                    if(idxProc != app.global().processCount()-1 && allintervals[idxProc*2+1] != -1){
                        FAssertLF(allintervals[idxProc*2] < allintervals[idxProc*2+1]);
                    }

                    FAssertLF(currentLimit < allintervals[idxProc*2]);
                    currentLimit = allintervals[idxProc*2+1];
                }
            }

            delete[] allintervals;
        }
        else{
            MPI_Gather(leftright, 2, FMpi::GetType(*leftright),
                       nullptr, 2, FMpi::GetType(*leftright), 0, app.global().getComm());
        }


        delete[] allPartsCount;
        delete[] outputArray;
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
        AddTest(&TestMpiTreeBuilder::RunTest,"Load a File, sort it, merge it, and Equalize it (4 steps)");
    }

public:
    TestMpiTreeBuilder(int argc,char ** argv) : FUTesterMpi(argc,argv){
    }


};

TestClassMpi(TestMpiTreeBuilder);
