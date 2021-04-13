// See LICENCE file at project root
#ifndef FMPISTATICTREEBUILDER_H
#define FMPISTATICTREEBUILDER_H

#include "../Utils/FMpi.hpp"
#include "../Utils/FQuickSortMpi.hpp"
#include "../Utils/FBitonicSort.hpp"
#include "../Utils/FTic.hpp"
#include "../Utils/FEnv.hpp"

#include "../Utils/FMemUtils.hpp"

#include "../Containers/FVector.hpp"

#include "../Utils/FLeafBalance.hpp"
#include "../Utils/FEqualize.hpp"

#include "../Containers/FCoordinateComputer.hpp"

/**
 * This class manage the loading of particles for the mpi version.
 * It work in several steps.
 * First it load the data from a file or an array and sort them amon the MPI processes.
 * Then, it carrefully manage if a leaf is shared by multiple processes.
 * Finally it balances the data using an external interval builder.
 *
 */
template<class FReal, class ParticleClass>
class FMpiStaticTreeBuilder{
public:
    /**
     * A particle may not have a MortonIndex Method (set/get morton index)
     * But in this algorithm they are sorted based on their morton indexes.
     * So an IndexedParticle is storing a real particle + its index.
     */
    struct IndexedParticle{
    public:
        MortonIndex index;
        ParticleClass particle;

        operator MortonIndex() const {
            return this->index;
        }
    };

    //////////////////////////////////////////////////////////////////////////
    // The builder function
    //////////////////////////////////////////////////////////////////////////

    template <class ContainerClass>
    static void DistributeArrayToContainer(const FMpi::FComm& communicator, const ParticleClass originalParticlesArray[], const FSize originalNbParticles,
                                           const FPoint<FReal>& boxCenter, const FReal boxWidth, const int treeHeight,
                                           ContainerClass* particleSaver){
        // Allocate the particles array
        std::unique_ptr<IndexedParticle[]> originalParticlesCore(new IndexedParticle[originalNbParticles]);
        FMemUtils::memset(originalParticlesCore.get(), 0, sizeof(IndexedParticle) * originalNbParticles);

        FPoint<FReal> boxCorner(boxCenter - (boxWidth/2));
        FTreeCoordinate host;
        const FReal boxWidthAtLeafLevel = boxWidth / FReal(1 << (treeHeight - 1) );

        FLOG(FTic counterTime);

        MortonIndex minMaxIndexes[2];
        minMaxIndexes[0] = std::numeric_limits<decltype(MortonIndex())>::max();
        minMaxIndexes[1] = std::numeric_limits<decltype(MortonIndex())>::min();

        // Fill the array and compute the morton index
        for(FSize idxPart = 0 ; idxPart < originalNbParticles ; ++idxPart){
            originalParticlesCore[idxPart].particle = originalParticlesArray[idxPart];
            host.setX( FCoordinateComputer::GetTreeCoordinate<FReal>( originalParticlesCore[idxPart].particle.getPosition().getX() - boxCorner.getX(), boxWidth, boxWidthAtLeafLevel,
                                           treeHeight ));
            host.setY( FCoordinateComputer::GetTreeCoordinate<FReal>( originalParticlesCore[idxPart].particle.getPosition().getY() - boxCorner.getY(), boxWidth, boxWidthAtLeafLevel,
                                           treeHeight ));
            host.setZ( FCoordinateComputer::GetTreeCoordinate<FReal>( originalParticlesCore[idxPart].particle.getPosition().getZ() - boxCorner.getZ(), boxWidth, boxWidthAtLeafLevel,
                                           treeHeight ));

            originalParticlesCore[idxPart].index = host.getMortonIndex();

            minMaxIndexes[0] = std::min(host.getMortonIndex(),minMaxIndexes[0]);
            minMaxIndexes[1] = std::max(host.getMortonIndex(),minMaxIndexes[1]);
        }

        FQuickSort<IndexedParticle>::QsOmp(originalParticlesCore.get(), originalNbParticles);

        MortonIndex globalMinMaxIndexes[2];
        FMpi::MpiAssert(MPI_Allreduce(&minMaxIndexes[0],
                            &globalMinMaxIndexes[0],
                            1,
                            MPI_LONG_LONG,
                            MPI_MIN,
                            communicator.getComm()), __LINE__);

        FMpi::MpiAssert(MPI_Allreduce(&minMaxIndexes[1],
                            &globalMinMaxIndexes[1],
                            1,
                            MPI_LONG_LONG,
                            MPI_MAX,
                            communicator.getComm()), __LINE__);


        const int nb_processes = communicator.processCount();
        const int my_rank = communicator.processId();

        const MortonIndex intervalSize = 1 + globalMinMaxIndexes[1] - globalMinMaxIndexes[0];
        const MortonIndex stepNbCells = (intervalSize+nb_processes-1)/nb_processes;

        std::vector<int> nb_items_to_send(nb_processes, 0);
        for(FSize idxPart = 0 ; idxPart < originalNbParticles ; ++idxPart){
            nb_items_to_send[(originalParticlesCore[idxPart].index - globalMinMaxIndexes[0])/stepNbCells] += 1;
        }

        std::vector<int> offset_items_to_send(nb_processes+1, 0);
        for(int idxProc = 0; idxProc < nb_processes ; ++idxProc){
            assert(std::numeric_limits<int>::max()-offset_items_to_send[idxProc] >= nb_items_to_send[idxProc]);
            offset_items_to_send[idxProc+1] = offset_items_to_send[idxProc] + nb_items_to_send[idxProc];
        }

        std::vector<int> nb_items_to_sendrecv_all(nb_processes*nb_processes);
        FMpi::MpiAssert(MPI_Allgather(const_cast<int*>(nb_items_to_send.data()), nb_processes, MPI_INT,
                                  nb_items_to_sendrecv_all.data(), nb_processes, MPI_INT,
                                  communicator.getComm()), __LINE__);

        int total_to_recv = 0;
        std::vector<int> nb_items_to_recv(nb_processes, 0);
        std::vector<int> offset_items_to_recv(nb_processes+1, 0);
        for(int idx_proc = 0 ; idx_proc < nb_processes ; ++idx_proc){
            const int nbrecv = nb_items_to_sendrecv_all[idx_proc*nb_processes + my_rank];
            assert(static_cast<long long int>(total_to_recv) + static_cast<long long int>(nbrecv) <= std::numeric_limits<int>::max());
            total_to_recv += nbrecv;
            nb_items_to_recv[idx_proc] = nbrecv;
            assert(static_cast<long long int>(nb_items_to_recv[idx_proc]) + static_cast<long long int>(offset_items_to_recv[idx_proc]) <= std::numeric_limits<int>::max());
            offset_items_to_recv[idx_proc+1] = nb_items_to_recv[idx_proc]
                                                    + offset_items_to_recv[idx_proc];
        }


        std::unique_ptr<IndexedParticle[]> out_to_recv(new IndexedParticle[total_to_recv]);

        {
            // Convert to byte!
            const int sizeOfParticle = sizeof(IndexedParticle);

            for(int& val : nb_items_to_send){
                assert(std::numeric_limits<int>::max()-sizeOfParticle >= val);
                val *= sizeOfParticle;
            }
            for(int& val : offset_items_to_send){
                assert(std::numeric_limits<int>::max()-sizeOfParticle >= val);
                val *= sizeOfParticle;
            }
            for(int& val : nb_items_to_recv){
                assert(std::numeric_limits<int>::max()-sizeOfParticle >= val);
                val *= sizeOfParticle;
            }
            for(int& val : offset_items_to_recv){
                assert(std::numeric_limits<int>::max()-sizeOfParticle >= val);
                val *= sizeOfParticle;
            }

            FMpi::MpiAssert(MPI_Alltoallv(originalParticlesCore.get(), const_cast<int*>(nb_items_to_send.data()),
                                  const_cast<int*>(offset_items_to_send.data()), MPI_BYTE, out_to_recv.get(),
                                  const_cast<int*>(nb_items_to_recv.data()), const_cast<int*>(offset_items_to_recv.data()), MPI_BYTE,
                                  communicator.getComm()), __LINE__);
        }

        for(FSize idxPart = 0 ; idxPart < total_to_recv ; ++idxPart){
            particleSaver->push(out_to_recv[idxPart].particle);
        }

#ifdef SCALFMM_USE_LOG
        /** To produce stats after the Equalize phase  */
        {
            const FSize finalNbParticles = particleSaver->getSize();

            if(communicator.processId() != 0){
                FMpi::MpiAssert(MPI_Gather(const_cast<FSize*>(&finalNbParticles),1,FMpi::GetType(finalNbParticles),nullptr,
                                           1,FMpi::GetType(finalNbParticles),0,communicator.getComm()), __LINE__);
            }
            else{
                const int nbProcs = communicator.processCount();
                std::unique_ptr<FSize[]> nbPartsPerProc(new FSize[nbProcs]);

                FMpi::MpiAssert(MPI_Gather(const_cast<FSize*>(&finalNbParticles),1,FMpi::GetType(finalNbParticles),nbPartsPerProc.get(),
                                           1,FMpi::GetType(finalNbParticles),0,communicator.getComm()), __LINE__);

                FReal averageNbParticles = 0;
                FSize minNbParticles = finalNbParticles;
                FSize maxNbParticles = finalNbParticles;

                for(int idxProc = 0 ; idxProc < nbProcs ; ++idxProc){
                    maxNbParticles = FMath::Max(maxNbParticles, nbPartsPerProc[idxProc]);
                    minNbParticles = FMath::Min(minNbParticles, nbPartsPerProc[idxProc]);
                    averageNbParticles += FReal(nbPartsPerProc[idxProc]);
                }
                averageNbParticles /= float(nbProcs);

                printf("Particles Distribution: End of Equalize Phase : \n \t Min number of parts : %lld \n \t Max number of parts : %lld \n \t Average number of parts : %e \n",
                       minNbParticles,maxNbParticles,averageNbParticles);
            }
        }
#endif
    }


};

#endif // FMPITREEBUILDER_H
