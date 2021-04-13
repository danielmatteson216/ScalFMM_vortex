// See LICENCE file at project root

// ==== CMAKE =====
// @FUSE_MPI
// ================


#ifndef FMPIFMAGENERICLOADER_HPP
#define FMPIFMAGENERICLOADER_HPP

#include <cstdlib>
#include <vector>

#include "Utils/FMpi.hpp"
#include "Files/FFmaGenericLoader.hpp"

template <class FReal>
class FMpiFmaGenericLoader : public FFmaGenericLoader<FReal> {
protected:
    using FFmaGenericLoader<FReal>::nbParticles;
    using FFmaGenericLoader<FReal>::file;
    using FFmaGenericLoader<FReal>::typeData;

  FSize myNbOfParticles;     //Number of particles that the calling process will manage
  MPI_Offset idxParticles;   //
  FSize start;               // number of my first parts in file
  size_t headerSize;
public:
    FMpiFmaGenericLoader(const std::string inFilename,const FMpi::FComm& comm, const bool /*useMpiIO*/ = false)
        : FFmaGenericLoader<FReal>(inFilename),myNbOfParticles(0),idxParticles(0),headerSize(0)
  {
    FSize startPart = comm.getLeft(nbParticles);
    FSize endPart   = comm.getRight(nbParticles);
    this->start = startPart;
    this->myNbOfParticles = endPart-startPart;
    std::cout << " startPart "<< startPart << " endPart " << endPart<<std::endl;
    std::cout << "Proc " << comm.processId() << " will hold " << myNbOfParticles << std::endl;

    if(this->binaryFile) {
        //This is header size in bytes
        //   MEANING :      sizeof(FReal)+nbAttr, nb of parts, boxWidth+boxCenter
        headerSize = sizeof(int)*2 + sizeof(FSize) + sizeof(FReal)*4;
        //To this header size, we had the parts that belongs to proc on my left
        file->seekg(headerSize + startPart*typeData[1]*sizeof(FReal));
    } else {
        // First finish to read the current line
        file->ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        for(int i = 0; i < startPart; ++i) {
            file->ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        }
    }
  }

  ~FMpiFmaGenericLoader(){
  }

  FSize getMyNumberOfParticles() const{
    return myNbOfParticles;
  }

  FSize getStart() const{
    return start;
  }

  /**
   * Given an index, get the one particle from this index
   */
  void fill1Particle(FReal*datas,FSize indexInFile){
    file->seekg(headerSize+(int(indexInFile)*FFmaGenericLoader<FReal>::getNbRecordPerline()*sizeof(FReal)));
    file->read((char*) datas,FFmaGenericLoader<FReal>::getNbRecordPerline()*sizeof(FReal));
  }

};
/**
 *
 * \brief Writes a set of distributed particles to an FMA formated file.
 *
 * The file may be in ASCII or binary mode. The example below shows how to use the class.
 *
 * \code
 * // Instanciate the writer with a binary fma file (extension .bfma).
 * \endcode
 * ----------------------------------------
 * FMA is a simple format to store particles in a file. It is organized as follow.
 *
 * file
 */
template <class FReal>
class FMpiFmaGenericWriter : public FFmaGenericWriter<FReal> {

protected:
  const FMpi* _parallelManager ;
  bool _writeDone ;
  int _headerSize;
  int _nbDataTowritePerRecord ;  //< number of data to write for one particle
  FSize _numberOfParticles ;     //< number of particle (global) to write in the file
  using FFmaGenericWriter<FReal>::file;
  MPI_File _mpiFile;              //< MPI pointer on data file (write mode)

public:
    /**
     * This constructor opens a file to be written to.
     *
     * - The opening mode is guessed from the file extension : `.fma` will open
     * in ASCII mode, `.bfma` will open in binary mode.
     *
     * @param filename the name of the file to open.
     */
  FMpiFmaGenericWriter(const std::string inFilename,  const FMpi& para)  : FFmaGenericWriter<FReal>(inFilename),
                       _parallelManager(&para),_writeDone(false),_headerSize(0),_nbDataTowritePerRecord(8),_numberOfParticles(0)
  {
    if ( ! this->isBinary()){
        std::cout << "FMpiFmaGenericWriter only works with binary file (.bfma)." << std::endl;
        std::exit(EXIT_FAILURE);
      }
    int fileIsOpen = MPI_File_open( _parallelManager->global().getComm(), const_cast<char*>(inFilename.c_str()),
                                    MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &_mpiFile );
    // Is it open?
    if(fileIsOpen != MPI_SUCCESS){
        std::cout << "Cannot create parallel file, FMpiFmaGenericWriter constructeur abort." << std::endl;
        std::exit(EXIT_FAILURE);
        return;
      }
  }
    /**
     * Writes the header of FMA file.
     *
     * Should be used if we write the particles with writeArrayOfReal method
     *
     * @param centerOfBox      The center of the Box (FPoint<FReal> class)
     * @param boxWidth         The width of the box
     * @param nbParticles      Number of particles in the box (or to save)
     * @param dataType         Size of the data type of the values in particle
     * @param nbDataPerRecord  Number of record/value per particle
     */
    void writeHeader(const FPoint<FReal> &centerOfBox,const FReal &boxWidth, const FSize &nbParticles,
                     const unsigned int  dataType, const unsigned int  nbDataPerRecord) {
//      * \code
//      *   DatatypeSize  Number_of_record_per_line
//      *   NB_particles  half_Box_width  Center_X  Center_Y  Center_Z
//      *   Particle_values
//      * \endcode
          _headerSize = 0 ;
      _nbDataTowritePerRecord = nbDataPerRecord ;
      _numberOfParticles     = nbParticles ;
      if(_parallelManager->global().processId()==0){
          int sizeType=0 ;
            int ierr = 0 ;
          MPI_Datatype mpiFSize_t = _parallelManager->GetType(nbParticles) ;
          MPI_Datatype mpiFReal_t = _parallelManager->GetType(boxWidth) ;
  //
          unsigned int typeFReal[2]  = {sizeof(FReal) , nbDataPerRecord};

          ierr =MPI_File_write_at(_mpiFile, 0, &typeFReal, 2, MPI_INT,  MPI_STATUS_IGNORE);
          MPI_Type_size(MPI_INT, &sizeType) ;
          _headerSize += sizeType*2 ;
          ierr =MPI_File_write_at(_mpiFile, _headerSize, &nbParticles, 1, mpiFSize_t,  MPI_STATUS_IGNORE);
          MPI_Type_size(mpiFSize_t, &sizeType) ;
          _headerSize += sizeType*1 ;

         FReal boxSim[4] = {boxWidth ,centerOfBox.getX() , centerOfBox.getX() , centerOfBox.getX() } ;

          ierr =MPI_File_write_at(_mpiFile, _headerSize, &boxSim[0], 4, mpiFReal_t,  MPI_STATUS_IGNORE);
         MPI_Type_size(mpiFReal_t, &sizeType) ;
         _headerSize += sizeType*4 ;
         // Build the header offset
         std::cout << " headerSize "<<   _headerSize << std::endl;
          }
         MPI_Bcast(&_headerSize,1,MPI_INT,0,_parallelManager->global().getComm());
          std::cout << "  _headerSize  " <<  _headerSize  <<std::endl;

    }
  ~FMpiFmaGenericWriter(){
      MPI_File_close(&_mpiFile );
  }

    /**
     *  Write all for all particles the position, physical values, potential and forces
     *
     * @param myOctree the octree
     * @param nbParticlesnumber of particles
     * @param mortonLeafDistribution the morton distribution of the leaves (this is a vecor of size 2* the number of MPI processes
     *
     */
    template <class OCTREECLASS>
    void writeDistributionOfParticlesFromGroupedOctree( OCTREECLASS &myOctree, const FSize& nbParticles, const std::vector<MortonIndex> &mortonLeafDistribution){
      //
      // Write the header
      int sizeType = 0,ierr = 0 ;
      FReal tt =0.0 ;
      MPI_Datatype mpiFSize_t = _parallelManager->GetType(nbParticles) ;
      MPI_Datatype mpiFReal_t = _parallelManager->GetType(tt) ;
      MPI_Type_size(mpiFReal_t, &sizeType) ;
      int myRank = _parallelManager->global().processId()  ;
      _headerSize = 0 ;
      //
      unsigned int typeFReal[2]  = {sizeof(FReal) , static_cast<unsigned int>(_nbDataTowritePerRecord)};
      if(myRank==0){
          ierr =MPI_File_write_at(_mpiFile, 0, &typeFReal, 2, MPI_INT,  MPI_STATUS_IGNORE);
        }
      MPI_Type_size(MPI_INT, &sizeType) ;
      _headerSize += sizeType*2 ;
      if(myRank==0){
          ierr =MPI_File_write_at(_mpiFile, _headerSize, &nbParticles, 1, mpiFSize_t,  MPI_STATUS_IGNORE);
        }
      MPI_Type_size(mpiFSize_t, &sizeType) ;
      _headerSize += sizeType*1 ;
      auto centerOfBox =myOctree.getBoxCenter()  ;
      FReal boxSim[4] = {myOctree.getBoxWidth()*0.5 , centerOfBox.getX() , centerOfBox.getX() , centerOfBox.getX() } ;

      if(myRank==0){
          ierr =MPI_File_write_at(_mpiFile, _headerSize, &boxSim[0], 4, mpiFReal_t,  MPI_STATUS_IGNORE);
        }
      if(ierr >0){
          std::cerr << "Error during the construction of the header in FMpiFmaGenericWriter::writeDistributionOfParticlesFromOctree"<<std::endl;
        }
      MPI_Type_size(mpiFReal_t, &sizeType) ;
      _headerSize += sizeType*4 ;
    //
    // Construct the local number of particles on my process
    FSize nbLocalParticles =0 ,maxPartLeaf =0;
    MortonIndex starIndex = mortonLeafDistribution[2*myRank], endIndex =  mortonLeafDistribution[2*myRank+1];
    myOctree.template forEachCellLeaf<typename OCTREECLASS::LeafClass_T >(
        [&](typename OCTREECLASS::GroupSymbolCellClass_T* gsymb ,
                     typename OCTREECLASS::GroupCellUpClass_T*   /* gmul */,
                     typename OCTREECLASS::GroupCellDownClass_T* /* gloc */,
                     typename OCTREECLASS::LeafClass_T * leafTarget
                     )
    {
        if (! (gsymb->getMortonIndex() < starIndex || gsymb->getMortonIndex() > endIndex)) {
            auto n =  leafTarget->getNbParticles();
            nbLocalParticles += n;
            maxPartLeaf = std::max(maxPartLeaf,n);
          }
      }
    );
    std::vector<FReal> particles(maxPartLeaf*_nbDataTowritePerRecord);
    // Build the offset for eaxh processes
    FSize before=0;  // Number of particles before me (rank < myrank)
    MPI_Scan(&nbLocalParticles,&before,1,mpiFSize_t,MPI_SUM,_parallelManager->global().getComm());
    before -= nbLocalParticles ;
    MPI_Offset offset = _headerSize + sizeType*_nbDataTowritePerRecord*before;
    //
    // Write particles in file
    myOctree.template forEachCellLeaf<typename OCTREECLASS::LeafClass_T >(
        [&](typename OCTREECLASS::GroupSymbolCellClass_T* gsymb ,
                     typename OCTREECLASS::GroupCellUpClass_T*   /* gmul */,
                     typename OCTREECLASS::GroupCellDownClass_T* /* gloc */,
                     typename OCTREECLASS::LeafClass_T * leafTarget
                     )
    {
        if (! (gsymb->getMortonIndex() < starIndex || gsymb->getMortonIndex() > endIndex)) {
            const FSize nbPartsInLeaf = leafTarget->getNbParticles();
            const FReal*const posX = leafTarget->getPositions()[0];
            const FReal*const posY = leafTarget->getPositions()[1];
            const FReal*const posZ = leafTarget->getPositions()[2];
            const FReal*const physicalValues = leafTarget->getPhysicalValues();
            const FReal*const forceX = leafTarget->getForcesX();
            const FReal*const forceY = leafTarget->getForcesY();
            const FReal*const forceZ = leafTarget->getForcesZ();
            const FReal*const potential = leafTarget->getPotentials();
            for (int i=0, k=0 ; i < nbPartsInLeaf ;++i,k+=_nbDataTowritePerRecord ) {
                particles[k] = posX[i];  particles[k+1] = posY[i];  particles[k+2] = posZ[i];
                particles[k+3] = physicalValues[i]; particles[k+4] = potential[i];
                particles[k+5] = forceX[i];  particles[k+6] = forceY[i];  particles[k+7] = forceZ[i];
              }
            MPI_File_write_at(_mpiFile, offset, particles.data(), static_cast<int>(_nbDataTowritePerRecord*nbPartsInLeaf), mpiFReal_t,  MPI_STATUS_IGNORE);
            offset+=sizeType*_nbDataTowritePerRecord*nbPartsInLeaf;
          }
      }
    );

    MPI_File_close(&_mpiFile );

    }

    /**
     *  Write all for all particles the position, physical values, potential and forces
     *
     * @param myOctree the octree
     * @param nbParticlesnumber of particles
     * @param nbLocalParticles number of local particles (on the MPI processus
     * @param mortonLeafDistribution the morton distribution of the leaves (this is a vector of size 2* the number of MPI processes
     *
     */
    template <class OCTREECLASS>
    void writeDistributionOfParticlesFromOctree( OCTREECLASS &myOctree,
                                                 const FSize& nbParticles,
                                                 const FSize& nbLocalParticles,
                                                 const std::vector<MortonIndex> &mortonLeafDistribution){
      //
      // Write the header
      int sizeType = 0,ierr = 0 ;
      FReal tt{}  ;

      MPI_Datatype mpiFSize_t = _parallelManager->GetType(nbParticles) ;
      MPI_Datatype mpiFReal_t = _parallelManager->GetType(tt) ;
      MPI_Type_size(mpiFReal_t, &sizeType) ;
      int myRank = _parallelManager->global().processId()  ;
      _headerSize = 0 ;
      //
      unsigned int typeFReal[2]  = {sizeof(FReal) , static_cast<unsigned int>(_nbDataTowritePerRecord)};
      if(myRank==0){
          ierr =MPI_File_write_at(_mpiFile, 0, &typeFReal, 2, MPI_INT,  MPI_STATUS_IGNORE);
        }
      MPI_Type_size(MPI_INT, &sizeType) ;
      _headerSize += sizeType*2 ;
      if(myRank==0){
          ierr = MPI_File_write_at(_mpiFile, _headerSize, const_cast<FSize*>(&nbParticles), 1, mpiFSize_t,  MPI_STATUS_IGNORE);
        }
      MPI_Type_size(mpiFSize_t, &sizeType) ;
      _headerSize += sizeType*1 ;
      auto centerOfBox =myOctree.getBoxCenter()  ;
      FReal boxSim[4] = {myOctree.getBoxWidth()*0.5 , centerOfBox.getX() , centerOfBox.getX() , centerOfBox.getX() } ;

      if(myRank==0){
          ierr =MPI_File_write_at(_mpiFile, _headerSize, &boxSim[0], 4, mpiFReal_t,  MPI_STATUS_IGNORE);
        }
      if(ierr >0){
          std::cerr << "Error during the construction of the header in FMpiFmaGenericWriter::writeDistributionOfParticlesFromOctree"<<std::endl;
        }
      MPI_Type_size(mpiFReal_t, &sizeType) ;
      _headerSize += sizeType*4 ;
    //
    // Construct the local number of particles on my process
    FSize maxPartLeaf =0, nn =0;
    MortonIndex starIndex = mortonLeafDistribution[2*myRank], endIndex =  mortonLeafDistribution[2*myRank+1];
  //  myOctree.template forEachCellLeaf<typename OCTREECLASS::LeafClass_T >(
          myOctree.forEachCellLeaf(
        [&](typename OCTREECLASS::CellClassType* cell,
                     typename OCTREECLASS::LeafClass_T * leaf)
    {
        if (! (cell->getMortonIndex() < starIndex || cell->getMortonIndex() > endIndex)) {
            auto n =  leaf->getSrc()->getNbParticles();
            maxPartLeaf = std::max(maxPartLeaf,n);
            nn += n ;
          }
      }
    );
          std::cout << "  nn " << nn << "  should be " << nbLocalParticles << std::endl;
    std::vector<FReal> particles(maxPartLeaf*_nbDataTowritePerRecord);
    // Build the offset for eaxh processes
    FSize before=0;  // Number of particles before me (rank < myrank)
    MPI_Scan(const_cast<FSize*>(&nbLocalParticles),&before,1,mpiFSize_t,MPI_SUM,_parallelManager->global().getComm());
    before -= nbLocalParticles ;
    MPI_Offset offset = _headerSize + sizeType*_nbDataTowritePerRecord*before;
    //
    // Write particles in file
    myOctree.forEachCellLeaf(
        [&](typename OCTREECLASS::CellClassType* cell,
          typename OCTREECLASS::LeafClass_T * leaf )
    {
        if (! (cell->getMortonIndex() < starIndex || cell->getMortonIndex() > endIndex)) {
            const FSize nbPartsInLeaf = leaf->getTargets()->getNbParticles();
            const FReal*const posX = leaf->getTargets()->getPositions()[0];
            const FReal*const posY = leaf->getTargets()->getPositions()[1];
            const FReal*const posZ = leaf->getTargets()->getPositions()[2];
            const FReal*const physicalValues = leaf->getTargets()->getPhysicalValues();
            const FReal*const forceX = leaf->getTargets()->getForcesX();
            const FReal*const forceY = leaf->getTargets()->getForcesY();
            const FReal*const forceZ = leaf->getTargets()->getForcesZ();
            const FReal*const potential = leaf->getTargets()->getPotentials();
            for (int i=0, k=0 ; i < nbPartsInLeaf ;++i,k+=_nbDataTowritePerRecord ) {
                particles[k] = posX[i];  particles[k+1] = posY[i];  particles[k+2] = posZ[i];
                particles[k+3] = physicalValues[i]; particles[k+4] = potential[i];
                particles[k+5] = forceX[i];  particles[k+6] = forceY[i];  particles[k+7] = forceZ[i];
              }
            MPI_File_write_at(_mpiFile, offset, particles.data(),
                              static_cast<int>(_nbDataTowritePerRecord*nbPartsInLeaf),
                              mpiFReal_t,  MPI_STATUS_IGNORE);
            offset+=sizeType*_nbDataTowritePerRecord*nbPartsInLeaf;
          }
      }
    );
#ifdef TODO
#endif
    MPI_File_close(&_mpiFile );


    }

//    /**
//     *  Write an array of data in a file Fill
//     *
//     * @param dataToWrite array of particles of type FReal
//     * @param nbData number of data per particle
//     * @param N number of particles
//     *
//     *   The size of the array is N*nbData
//     *
//     *   example
//     * \code
//     * FmaRParticle * const particles = new FmaRParticle[nbParticles];
//     * memset(particles, 0, sizeof(FmaRParticle) * nbParticles) ;
//     * ...
//     * FFmaGenericWriter<FReal> writer(filenameOut) ;
//     * Fwriter.writeHeader(Centre,BoxWith, nbParticles,*particles) ;
//     * Fwriter.writeArrayOfReal(particles, nbParticles);
//     * \endcode
//     */
//    void writeArrayOfReal(const FReal *dataToWrite, const FSize nbData, const FSize N){
//      /*
//      if(! _writeDone){
//          FSize previousNumberofParticles;
//          MPI_Scan(&N,&previousNumberofParticles,1,_parallelManager->GetType(N),MPI_SUM,_parallelManager->global().getComm());
//          FSize offset= previousNumberofParticles-N;
//          //To this header size, we had the parts that belongs to proc on my left
//          this->skipHeaderAndPart(offset) ;
//          FFmaGenericWriter<FReal>::writeArrayOfReal(dataToWrite,4,N) ;
//          //

//          std::cout <<" node " << _parallelManager->global().processId() << "Npart " << N << "  before Me" << previousNumberofParticles-N<< std::endl;
//          _writeDone  = true;
//        }
//      else {
//          std::cerr << " The writeArrayOfReal should be call only once !!!! "<< std::endl;
//          std::exit(EXIT_FAILURE);
//        }
//        */
//   }
//private:
//     void skipHeaderAndPart(const FSize &numberOfParticleToSkip){
//       if(this->binaryFile) {
//           //This is header size in bytes
//           //   MEANING :      sizeof(FReal)+nbAttr, nb of parts, boxWidth+boxCenter
//           _headerSize = sizeof(int)*2 + sizeof(FSize) + sizeof(FReal)*4;
//          file->seekg(_headerSize+numberOfParticleToSkip* 4/*FFmaGenericWriter<FReal>::getNbRecordPerline()*/*sizeof(FReal), std::ios::beg);

//       } else {
//           // First finish to read the current line
//           file->ignore(std::numeric_limits<std::streamsize>::max(), '\n');
//           for(int i = 0; i < numberOfParticleToSkip; ++i) {
//               file->ignore(std::numeric_limits<std::streamsize>::max(), '\n');
//           }
//       }
//     }



} ;

#endif //FMPIFMAGENERICLOADER_HPP
