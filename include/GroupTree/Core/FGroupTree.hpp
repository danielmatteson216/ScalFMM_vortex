// ==== CMAKE =====
//
// ================

//
#ifndef FGROUPTREE_HPP
#define FGROUPTREE_HPP
#include <vector>
#include <functional>

#include "Utils/FAssert.hpp"
#include "Utils/FPoint.hpp"
#include "Utils/FQuickSort.hpp"
#include "Containers/FTreeCoordinate.hpp"
#include "Containers/FCoordinateComputer.hpp"
#include "FGroupOfCells.hpp"
#include "FGroupOfParticles.hpp"
#include "FGroupAttachedLeaf.hpp"
#include "Kernels/P2P/FP2PParticleContainer.hpp"
#ifdef SCALFMM_USE_MPI
#include "FDistributedGroupTreeBuilder.hpp"
#endif

template <class FReal, class SymbolCellClass, class PoleCellClass, class LocalCellClass,
          class GroupAttachedLeafClass, unsigned NbSymbAttributes, unsigned NbAttributesPerParticle, class AttributeClass = FReal>
class FGroupTree {
public:
  typedef GroupAttachedLeafClass     BasicAttachedClass;   // Leaf
  typedef FGroupOfParticles<FReal, NbSymbAttributes, NbAttributesPerParticle,AttributeClass>     ParticleGroupClass;
  typedef FGroupOfCells<SymbolCellClass, PoleCellClass, LocalCellClass> CellGroupClass;
  typedef SymbolCellClass  GroupSymbolCellClass_T ;
  typedef LocalCellClass   GroupCellDownClass_T ;
  typedef PoleCellClass    GroupCellUpClass_T ;
  typedef GroupAttachedLeafClass     LeafClass_T;   // Leaf

protected:
  //< height of the tree (1 => only the root)
  const int _treeHeight;
  //< max number of cells in a block
  const int _nbElementsPerBlock;
  //< all the blocks of the tree
  std::vector<CellGroupClass*>* _cellBlocksPerLevel;
  //< all the blocks of leaves
  std::vector<ParticleGroupClass*> _particleBlocks;

  //< the space system center
  const FPoint<FReal> boxCenter;
  //< the space system corner (used to compute morton index)
  const FPoint<FReal> boxCorner;
  //< the space system width
  const FReal boxWidth;
  //< the width of a box at width level
  const FReal boxWidthAtLeafLevel;

public:
  typedef typename std::vector<CellGroupClass*>::iterator           CellGroupIterator;
  typedef typename std::vector<CellGroupClass*>::const_iterator     CellGroupConstIterator;
  typedef typename std::vector<ParticleGroupClass*>::iterator       ParticleGroupIterator;
  typedef typename std::vector<ParticleGroupClass*>::const_iterator ParticleGroupConstIterator;

  /** This constructor create a blocked octree from a usual octree
     * The cell are allocated as in the usual octree (no copy constructor are called!)
     * Once allocated each cell receive its morton index and tree coordinate.
     * No blocks are allocated at level 0.
     */
  template<class OctreeClass>
  FGroupTree()
  {}
  template<class OctreeClass>
  FGroupTree(const int in_treeHeight, const int in_nbElementsPerBlock, OctreeClass*const inOctreeSrc)
    : _treeHeight(in_treeHeight), _nbElementsPerBlock(in_nbElementsPerBlock), _cellBlocksPerLevel(nullptr),
      boxCenter(inOctreeSrc->getBoxCenter()), boxCorner(inOctreeSrc->getBoxCenter(),-(inOctreeSrc->getBoxWidth()/2)),
      boxWidth(inOctreeSrc->getBoxWidth()), boxWidthAtLeafLevel(inOctreeSrc->getBoxWidth()/FReal(1<<(in_treeHeight-1))){

    _cellBlocksPerLevel = new std::vector<CellGroupClass*>[_treeHeight];

    // Iterate on the tree and build
    typename OctreeClass::Iterator octreeIterator(inOctreeSrc);
    octreeIterator.gotoBottomLeft();

    { // First leaf level, we create leaves and cells groups
      const int idxLevel = _treeHeight-1;
      typename OctreeClass::Iterator avoidGotoLeft = octreeIterator;
      // For each cell at this level
      do {
          typename OctreeClass::Iterator blockIteratorInOctree = octreeIterator;
          // Move the iterator per _nbElementsPerBlock (or until it cannot move right)
          int sizeOfBlock = 1;
          FSize nbParticlesInGroup = octreeIterator.getCurrentLeaf()->getSrc()->getNbParticles();
          while(sizeOfBlock < _nbElementsPerBlock && octreeIterator.moveRight()){
              sizeOfBlock += 1;
              nbParticlesInGroup += octreeIterator.getCurrentLeaf()->getSrc()->getNbParticles();
            }

          // Create a block with the apropriate parameters
          CellGroupClass*const newBlock = new CellGroupClass(blockIteratorInOctree.getCurrentGlobalIndex(),
                                                             octreeIterator.getCurrentGlobalIndex()+1,
                                                             sizeOfBlock);
          FGroupOfParticles<FReal, NbSymbAttributes, NbAttributesPerParticle, AttributeClass>*const newParticleBlock = new FGroupOfParticles<FReal, NbSymbAttributes, NbAttributesPerParticle, AttributeClass>(blockIteratorInOctree.getCurrentGlobalIndex(),
                                                                                                                                                                                                               octreeIterator.getCurrentGlobalIndex()+1,
                                                                                                                                                                                                               sizeOfBlock, nbParticlesInGroup);

          // Initialize each cell of the block
          int cellIdInBlock = 0;
          size_t nbParticlesOffsetBeforeLeaf = 0;
          while(cellIdInBlock != sizeOfBlock){
              const MortonIndex newNodeIndex = blockIteratorInOctree.getCurrentCell()->getMortonIndex();
              const FTreeCoordinate newNodeCoordinate = blockIteratorInOctree.getCurrentCell()->getCoordinate();
              // Add cell
              newBlock->newCell(newNodeIndex, cellIdInBlock);

              SymbolCellClass& symbolic = newBlock->getSymbolic(cellIdInBlock);
              symbolic.setMortonIndex(newNodeIndex);
              symbolic.setCoordinate(newNodeCoordinate);
              symbolic.setLevel(idxLevel);

              // Add leaf
              nbParticlesOffsetBeforeLeaf = newParticleBlock->newLeaf(newNodeIndex, cellIdInBlock,
                                                                      blockIteratorInOctree.getCurrentLeaf()->getSrc()->getNbParticles(),
                                                                      nbParticlesOffsetBeforeLeaf);

              BasicAttachedClass attachedLeaf = newParticleBlock->template getLeaf<BasicAttachedClass>(cellIdInBlock);
              attachedLeaf.copyFromContainer(blockIteratorInOctree.getCurrentLeaf()->getSrc(), 0);

              cellIdInBlock += 1;
              blockIteratorInOctree.moveRight();
            }

          // Keep the block
          newBlock->declare_mine();
          _cellBlocksPerLevel[idxLevel].push_back(newBlock);
          _particleBlocks.push_back(newParticleBlock);

          // If we can move right then add another block
        } while(octreeIterator.moveRight());

      avoidGotoLeft.moveUp();
      octreeIterator = avoidGotoLeft;
    }

    // For each level from heigth - 2 to 1
    for(int idxLevel = _treeHeight-2; idxLevel > 0 ; --idxLevel){
        typename OctreeClass::Iterator avoidGotoLeft = octreeIterator;
        // For each cell at this level
        do {
            typename OctreeClass::Iterator blockIteratorInOctree = octreeIterator;
            // Move the iterator per _nbElementsPerBlock (or until it cannot move right)
            int sizeOfBlock = 1;
            while(sizeOfBlock < _nbElementsPerBlock && octreeIterator.moveRight()){
                sizeOfBlock += 1;
              }

            // Create a block with the apropriate parameters
            CellGroupClass*const newBlock = new CellGroupClass(blockIteratorInOctree.getCurrentGlobalIndex(),
                                                               octreeIterator.getCurrentGlobalIndex()+1,
                                                               sizeOfBlock);
            // Initialize each cell of the block
            int cellIdInBlock = 0;
            while(cellIdInBlock != sizeOfBlock){
                const MortonIndex newNodeIndex = blockIteratorInOctree.getCurrentCell()->getMortonIndex();
                const FTreeCoordinate newNodeCoordinate = blockIteratorInOctree.getCurrentCell()->getCoordinate();
                newBlock->newCell(newNodeIndex, cellIdInBlock);

                SymbolCellClass& symbolic = newBlock->getSymbolic(cellIdInBlock);
                symbolic.setMortonIndex(newNodeIndex);
                symbolic.setCoordinate(newNodeCoordinate);
                symbolic.setLevel(idxLevel);

                cellIdInBlock += 1;
                blockIteratorInOctree.moveRight();
              }

            // Keep the block
            newBlock->declare_mine();
            _cellBlocksPerLevel[idxLevel].push_back(newBlock);

            // If we can move right then add another block
          } while(octreeIterator.moveRight());

        avoidGotoLeft.moveUp();
        octreeIterator = avoidGotoLeft;
      }
  }

  /**
     * This constructor create a group tree from a particle container index.
     * The morton index are computed and the particles are sorted in a first stage.
     * Then the leaf level is done.
     * Finally the other leve are proceed one after the other.
     * It should be easy to make it parallel using for and tasks.
     * If no limite give inLeftLimite = -1
     */
  template<class ParticleContainer>
  FGroupTree(const int in_treeHeight, const FReal inBoxWidth, const FPoint<FReal>& inBoxCenter,
             const int in_nbElementsPerBlock, ParticleContainer* inParticlesContainer,
             const bool particlesAreSorted = false, MortonIndex inLeftLimite = -1):
    _treeHeight(in_treeHeight),_nbElementsPerBlock(in_nbElementsPerBlock),_cellBlocksPerLevel(nullptr),
    boxCenter(inBoxCenter), boxCorner(inBoxCenter,-(inBoxWidth/2)), boxWidth(inBoxWidth),
    boxWidthAtLeafLevel(inBoxWidth/FReal(1<<(in_treeHeight-1)))
  {
    _cellBlocksPerLevel = new std::vector<CellGroupClass*>[_treeHeight];

    MortonIndex* currentBlockIndexes = new MortonIndex[_nbElementsPerBlock];
    // First we work at leaf level
    {
      // Build morton index for particles
      struct ParticleSortingStruct{
        FSize originalIndex;
        MortonIndex mindex;
      };
      // Convert position to morton index
      const FSize nbParticles = inParticlesContainer->getNbParticles();
      ParticleSortingStruct* particlesToSort = new ParticleSortingStruct[nbParticles];
      {
        const FReal* xpos = inParticlesContainer->getPositions()[0];
        const FReal* ypos = inParticlesContainer->getPositions()[1];
        const FReal* zpos = inParticlesContainer->getPositions()[2];

        for(FSize idxPart = 0 ; idxPart < nbParticles ; ++idxPart){
            const FTreeCoordinate host = FCoordinateComputer::GetCoordinateFromPositionAndCorner<FReal>(this->boxCorner, this->boxWidth,
                                                                                                        _treeHeight,
                                                                                                        FPoint<FReal>(xpos[idxPart], ypos[idxPart], zpos[idxPart]) );
            const MortonIndex particleIndex = host.getMortonIndex();
            particlesToSort[idxPart].mindex = particleIndex;
            particlesToSort[idxPart].originalIndex = idxPart;
          }
      }

      // Sort if needed
      if(particlesAreSorted == false){
          FQuickSort<ParticleSortingStruct, FSize>::QsOmp(particlesToSort, nbParticles, [](const ParticleSortingStruct& v1, const ParticleSortingStruct& v2){
              return v1.mindex <= v2.mindex;
            });
        }

      FAssertLF(nbParticles == 0 || inLeftLimite < particlesToSort[0].mindex);
      // Convert to block
      const int idxLevel = (_treeHeight - 1);
      FSize* nbParticlesPerLeaf = new FSize[_nbElementsPerBlock];
      FSize firstParticle = 0;
      // We need to proceed each group in sub level
      while(firstParticle != nbParticles){
          int sizeOfBlock = 0;
          FSize lastParticle = firstParticle;
          // Count until end of sub group is reached or we have enough cells
          while(sizeOfBlock < _nbElementsPerBlock && lastParticle < nbParticles){
              if(sizeOfBlock == 0 || currentBlockIndexes[sizeOfBlock-1] != particlesToSort[lastParticle].mindex){
                  currentBlockIndexes[sizeOfBlock] = particlesToSort[lastParticle].mindex;
                  nbParticlesPerLeaf[sizeOfBlock]  = 1;
                  sizeOfBlock += 1;
                }
              else{
                  nbParticlesPerLeaf[sizeOfBlock-1] += 1;
                }
              lastParticle += 1;
            }
          while(lastParticle < nbParticles && currentBlockIndexes[sizeOfBlock-1] == particlesToSort[lastParticle].mindex){
              nbParticlesPerLeaf[sizeOfBlock-1] += 1;
              lastParticle += 1;
            }

          // Create a group
          CellGroupClass*const newBlock = new CellGroupClass(currentBlockIndexes[0],
              currentBlockIndexes[sizeOfBlock-1]+1,
              sizeOfBlock);
          FGroupOfParticles<FReal, NbSymbAttributes, NbAttributesPerParticle, AttributeClass>*const newParticleBlock = new FGroupOfParticles<FReal, NbSymbAttributes, NbAttributesPerParticle, AttributeClass>(currentBlockIndexes[0],
              currentBlockIndexes[sizeOfBlock-1]+1,
              sizeOfBlock, lastParticle-firstParticle);

          // Init cells
          size_t nbParticlesOffsetBeforeLeaf = 0;
          FSize offsetParticles = firstParticle;
          for(int cellIdInBlock = 0; cellIdInBlock != sizeOfBlock ; ++cellIdInBlock){
              newBlock->newCell(currentBlockIndexes[cellIdInBlock], cellIdInBlock);

              SymbolCellClass& symbolic = newBlock->getSymbolic(cellIdInBlock);
              symbolic.setMortonIndex(currentBlockIndexes[cellIdInBlock]);
              FTreeCoordinate coord;
              coord.setPositionFromMorton(currentBlockIndexes[cellIdInBlock]);
              symbolic.setCoordinate(coord);
              symbolic.setLevel(idxLevel);

              // Add leaf
              nbParticlesOffsetBeforeLeaf = newParticleBlock->newLeaf(currentBlockIndexes[cellIdInBlock], cellIdInBlock,
                                                                      nbParticlesPerLeaf[cellIdInBlock], nbParticlesOffsetBeforeLeaf);

              BasicAttachedClass attachedLeaf = newParticleBlock->template getLeaf<BasicAttachedClass>(cellIdInBlock);
              // Copy each particle from the original position
              for(FSize idxPart = 0 ; idxPart < nbParticlesPerLeaf[cellIdInBlock] ; ++idxPart){
                  attachedLeaf.setParticle(idxPart, particlesToSort[idxPart + offsetParticles].originalIndex, inParticlesContainer);
                }
              offsetParticles += nbParticlesPerLeaf[cellIdInBlock];
            }

          // Keep the block
          newBlock->declare_mine();
          _cellBlocksPerLevel[idxLevel].push_back(newBlock);
          _particleBlocks.push_back(newParticleBlock);

          sizeOfBlock = 0;
          firstParticle = lastParticle;
        }
      delete[] nbParticlesPerLeaf;
      delete[] particlesToSort;
    }

    // For each level from heigth - 2 to 1
    for(int idxLevel = _treeHeight-2; idxLevel > 0 ; --idxLevel){
        inLeftLimite = (inLeftLimite == -1 ? inLeftLimite : (inLeftLimite>>3));

        CellGroupConstIterator iterChildCells = _cellBlocksPerLevel[idxLevel+1].begin();
        const CellGroupConstIterator iterChildEndCells = _cellBlocksPerLevel[idxLevel+1].end();

        // Skip blocks that do not respect limit
        while(iterChildCells != iterChildEndCells
              && ((*iterChildCells)->getEndingIndex()>>3) <= inLeftLimite){
            ++iterChildCells;
          }
        // If lower level is empty or all blocks skiped stop here
        if(iterChildCells == iterChildEndCells){
            break;
          }

        MortonIndex currentCellIndex = (*iterChildCells)->getStartingIndex();
        if((currentCellIndex>>3) <= inLeftLimite) currentCellIndex = ((inLeftLimite+1)<<3);
        int sizeOfBlock = 0;

        // We need to proceed each group in sub level
        while(iterChildCells != iterChildEndCells){
            // Count until end of sub group is reached or we have enough cells
            while(sizeOfBlock < _nbElementsPerBlock && iterChildCells != iterChildEndCells ){
                if((sizeOfBlock == 0 || currentBlockIndexes[sizeOfBlock-1] != (currentCellIndex>>3))
                   && (*iterChildCells)->exists(currentCellIndex)){
                    currentBlockIndexes[sizeOfBlock] = (currentCellIndex>>3);
                    sizeOfBlock += 1;
                    currentCellIndex = (((currentCellIndex>>3)+1)<<3);
                  }
                else{
                    currentCellIndex += 1;
                  }
                // If we are at the end of the sub group, move to next
                while(iterChildCells != iterChildEndCells && (*iterChildCells)->getEndingIndex() <= currentCellIndex){
                    ++iterChildCells;
                    // Update morton index
                    if(iterChildCells != iterChildEndCells && currentCellIndex < (*iterChildCells)->getStartingIndex()){
                        currentCellIndex = (*iterChildCells)->getStartingIndex();
                      }
                  }
              }

            // If group is full
            if(sizeOfBlock == _nbElementsPerBlock || (sizeOfBlock && iterChildCells == iterChildEndCells)){
                // Create a group
                CellGroupClass*const newBlock = new CellGroupClass(currentBlockIndexes[0],
                    currentBlockIndexes[sizeOfBlock-1]+1,
                    sizeOfBlock);
                // Init cells
                for(int cellIdInBlock = 0; cellIdInBlock != sizeOfBlock ; ++cellIdInBlock){
                    newBlock->newCell(currentBlockIndexes[cellIdInBlock], cellIdInBlock);

                    SymbolCellClass& symbolic = newBlock->getSymbolic(cellIdInBlock);
                    symbolic.setMortonIndex(currentBlockIndexes[cellIdInBlock]);
                    FTreeCoordinate coord;
                    coord.setPositionFromMorton(currentBlockIndexes[cellIdInBlock]);
                    symbolic.setCoordinate(coord);
                    symbolic.setLevel(idxLevel);
                  }

                // Keep the block
                newBlock->declare_mine();
                _cellBlocksPerLevel[idxLevel].push_back(newBlock);

                sizeOfBlock = 0;
              }
          }
      }
    delete[] currentBlockIndexes;
  }

  /**
     * This constructor create a group tree from a particle container index.
     * The morton index are computed and the particles are sorted in a first stage.
     * Then the leaf level is done.
     * Finally the other leve are proceed one after the other.
     * It should be easy to make it parallel using for and tasks.
     * If no limite give inLeftLimite = -1
     * The cover ration is the minimum pourcentage of cell that should
     * exist in a group (0 means no limite, 1 means the block must be dense)
     * oneParent should be turned on if it is better to have one block parent
     * per sublock (in case of have the cost of FMM that increase with the level
     * this could be an asset).
     */
  template<class ParticleContainer>
  FGroupTree(const int in_treeHeight, const FReal inBoxWidth, const FPoint<FReal>& inBoxCenter,
             const int in_nbElementsPerBlock, ParticleContainer* inParticlesContainer,
             const bool particlesAreSorted, const bool oneParent,
             const FReal inCoverRatio = 0.0, MortonIndex inLeftLimite = -1):
    _treeHeight(in_treeHeight),_nbElementsPerBlock(in_nbElementsPerBlock),_cellBlocksPerLevel(nullptr),
    boxCenter(inBoxCenter), boxCorner(inBoxCenter,-(inBoxWidth/2)), boxWidth(inBoxWidth),
    boxWidthAtLeafLevel(inBoxWidth/FReal(1<<(in_treeHeight-1)))
  {

    FAssertLF(inCoverRatio == 0.0 || oneParent == true, "If a ratio is choosen oneParent should be turned on");
    const bool userCoverRatio = (inCoverRatio != 0.0);

    _cellBlocksPerLevel = new std::vector<CellGroupClass*>[_treeHeight];

    MortonIndex* currentBlockIndexes = new MortonIndex[_nbElementsPerBlock];
    // First we work at leaf level
    {
      // Build morton index for particles
      struct ParticleSortingStruct{
        FSize originalIndex;
        MortonIndex mindex;
      };
      // Convert position to morton index
      const FSize nbParticles = inParticlesContainer->getNbParticles();
      ParticleSortingStruct* particlesToSort = new ParticleSortingStruct[nbParticles];
      {
        const FReal* xpos = inParticlesContainer->getPositions()[0];
        const FReal* ypos = inParticlesContainer->getPositions()[1];
        const FReal* zpos = inParticlesContainer->getPositions()[2];

        for(FSize idxPart = 0 ; idxPart < nbParticles ; ++idxPart){
            const FTreeCoordinate host = FCoordinateComputer::GetCoordinateFromPositionAndCorner<FReal>(this->boxCorner, this->boxWidth,
                                                                                                        _treeHeight,
                                                                                                        FPoint<FReal>(xpos[idxPart], ypos[idxPart], zpos[idxPart]) );
            const MortonIndex particleIndex = host.getMortonIndex();
            particlesToSort[idxPart].mindex = particleIndex;
            particlesToSort[idxPart].originalIndex = idxPart;
          }
      }

      // Sort if needed
      if(particlesAreSorted == false){
          FQuickSort<ParticleSortingStruct, FSize>::QsOmp(particlesToSort, nbParticles, [](const ParticleSortingStruct& v1, const ParticleSortingStruct& v2){
              return v1.mindex <= v2.mindex;
            });
        }

      FAssertLF(nbParticles == 0 || inLeftLimite < particlesToSort[0].mindex);

      // Convert to block
      const int idxLevel = (_treeHeight - 1);
      int* nbParticlesPerLeaf = new int[_nbElementsPerBlock];
      int firstParticle = 0;
      // We need to proceed each group in sub level
      while(firstParticle != nbParticles){
          int sizeOfBlock = 0;
          int lastParticle = firstParticle;
          // Count until end of sub group is reached or we have enough cells
          while(sizeOfBlock < _nbElementsPerBlock && lastParticle < nbParticles
                && (userCoverRatio == false
                    || sizeOfBlock == 0
                    || currentBlockIndexes[sizeOfBlock-1] == particlesToSort[lastParticle].mindex
                    || (FReal(sizeOfBlock+1)/FReal(particlesToSort[lastParticle].mindex-particlesToSort[firstParticle].mindex)) >= inCoverRatio)){
              if(sizeOfBlock == 0 || currentBlockIndexes[sizeOfBlock-1] != particlesToSort[lastParticle].mindex){
                  currentBlockIndexes[sizeOfBlock] = particlesToSort[lastParticle].mindex;
                  nbParticlesPerLeaf[sizeOfBlock]  = 1;
                  sizeOfBlock += 1;
                }
              else{
                  nbParticlesPerLeaf[sizeOfBlock-1] += 1;
                }
              lastParticle += 1;
            }
          while(lastParticle < nbParticles && currentBlockIndexes[sizeOfBlock-1] == particlesToSort[lastParticle].mindex){
              nbParticlesPerLeaf[sizeOfBlock-1] += 1;
              lastParticle += 1;
            }

          // Create a group
          CellGroupClass*const newBlock = new CellGroupClass(currentBlockIndexes[0],
              currentBlockIndexes[sizeOfBlock-1]+1,
              sizeOfBlock);
          FGroupOfParticles<FReal, NbSymbAttributes, NbAttributesPerParticle, AttributeClass>*const newParticleBlock = new FGroupOfParticles<FReal, NbSymbAttributes, NbAttributesPerParticle, AttributeClass>(currentBlockIndexes[0],
              currentBlockIndexes[sizeOfBlock-1]+1,
              sizeOfBlock, lastParticle-firstParticle);

          // Init cells
          size_t nbParticlesOffsetBeforeLeaf = 0;
          int offsetParticles = firstParticle;
          for(int cellIdInBlock = 0; cellIdInBlock != sizeOfBlock ; ++cellIdInBlock){
              newBlock->newCell(currentBlockIndexes[cellIdInBlock], cellIdInBlock);

              SymbolCellClass& symbolic = newBlock->getSymbolic(cellIdInBlock);
              symbolic.setMortonIndex(currentBlockIndexes[cellIdInBlock]);
              FTreeCoordinate coord;
              coord.setPositionFromMorton(currentBlockIndexes[cellIdInBlock]);
              symbolic.setCoordinate(coord);
              symbolic.setLevel(idxLevel);

              // Add leaf
              nbParticlesOffsetBeforeLeaf = newParticleBlock->newLeaf(currentBlockIndexes[cellIdInBlock], cellIdInBlock,
                                                                      nbParticlesPerLeaf[cellIdInBlock], nbParticlesOffsetBeforeLeaf);

              BasicAttachedClass attachedLeaf = newParticleBlock->template getLeaf<BasicAttachedClass>(cellIdInBlock);
              // Copy each particle from the original position
              for(FSize idxPart = 0 ; idxPart < nbParticlesPerLeaf[cellIdInBlock] ; ++idxPart){
                  attachedLeaf.setParticle(idxPart, particlesToSort[idxPart + offsetParticles].originalIndex, inParticlesContainer);
                }
              offsetParticles += nbParticlesPerLeaf[cellIdInBlock];
            }

          // Keep the block
          newBlock->declare_mine();
          _cellBlocksPerLevel[idxLevel].push_back(newBlock);
          _particleBlocks.push_back(newParticleBlock);

          sizeOfBlock = 0;
          firstParticle = lastParticle;
        }
      delete[] nbParticlesPerLeaf;
      delete[] particlesToSort;
    }


    // For each level from heigth - 2 to 1
    for(int idxLevel = _treeHeight-2; idxLevel > 0 ; --idxLevel){
        inLeftLimite = (inLeftLimite == -1 ? inLeftLimite : (inLeftLimite>>3));

        CellGroupConstIterator iterChildCells = _cellBlocksPerLevel[idxLevel+1].begin();
        const CellGroupConstIterator iterChildEndCells = _cellBlocksPerLevel[idxLevel+1].end();

        // Skip blocks that do not respect limit
        while(iterChildCells != iterChildEndCells
              && ((*iterChildCells)->getEndingIndex()>>3) <= inLeftLimite){
            ++iterChildCells;
          }
        // If lower level is empty or all blocks skiped stop here
        if(iterChildCells == iterChildEndCells){
            break;
          }

        MortonIndex currentCellIndex = (*iterChildCells)->getStartingIndex();
        if((currentCellIndex>>3) <= inLeftLimite) currentCellIndex = ((inLeftLimite+1)<<3);
        int sizeOfBlock = 0;

        if(oneParent == false){
            // We need to proceed each group in sub level
            while(iterChildCells != iterChildEndCells){
                // Count until end of sub group is reached or we have enough cells
                while(sizeOfBlock < _nbElementsPerBlock && iterChildCells != iterChildEndCells ){
                    if((sizeOfBlock == 0 || currentBlockIndexes[sizeOfBlock-1] != (currentCellIndex>>3))
                       && (*iterChildCells)->exists(currentCellIndex)){
                        currentBlockIndexes[sizeOfBlock] = (currentCellIndex>>3);
                        sizeOfBlock += 1;
                        currentCellIndex = (((currentCellIndex>>3)+1)<<3);
                      }
                    else{
                        currentCellIndex += 1;
                      }
                    // If we are at the end of the sub group, move to next
                    while(iterChildCells != iterChildEndCells && (*iterChildCells)->getEndingIndex() <= currentCellIndex){
                        ++iterChildCells;
                        // Update morton index
                        if(iterChildCells != iterChildEndCells && currentCellIndex < (*iterChildCells)->getStartingIndex()){
                            currentCellIndex = (*iterChildCells)->getStartingIndex();
                          }
                      }
                  }

                // If group is full
                if(sizeOfBlock == _nbElementsPerBlock || (sizeOfBlock && iterChildCells == iterChildEndCells)){
                    // Create a group
                    CellGroupClass*const newBlock = new CellGroupClass(currentBlockIndexes[0],
                        currentBlockIndexes[sizeOfBlock-1]+1,
                        sizeOfBlock);
                    // Init cells
                    for(int cellIdInBlock = 0; cellIdInBlock != sizeOfBlock ; ++cellIdInBlock){
                        newBlock->newCell(currentBlockIndexes[cellIdInBlock], cellIdInBlock);

                        SymbolCellClass& symbolic = newBlock->getSymbolic(cellIdInBlock);
                        symbolic.setMortonIndex(currentBlockIndexes[cellIdInBlock]);
                        FTreeCoordinate coord;
                        coord.setPositionFromMorton(currentBlockIndexes[cellIdInBlock]);
                        symbolic.setCoordinate(coord);
                        symbolic.setLevel(idxLevel);
                      }

                    // Keep the block
                    newBlock->declare_mine();
                    _cellBlocksPerLevel[idxLevel].push_back(newBlock);

                    sizeOfBlock = 0;
                  }
              }
          }
        else{
            // We need to proceed each group in sub level
            while(iterChildCells != iterChildEndCells){
                // We want one parent group per child group so we will stop the parent group
                // when we arrive to the same parent as lastChildIndex (which is lastChildIndex>>3)
                const MortonIndex lastChildIndex = ((*iterChildCells)->getEndingIndex()-1);
                // Count until end of sub group is reached or we passe the requested parent
                while( iterChildCells != iterChildEndCells
                       && (currentCellIndex>>3) <= (lastChildIndex>>3) ){
                    // Proceed until the requested parent
                    while(currentCellIndex != (*iterChildCells)->getEndingIndex()
                          && (currentCellIndex>>3) <= (lastChildIndex>>3) ){
                        if((*iterChildCells)->exists(currentCellIndex)){
                            currentBlockIndexes[sizeOfBlock] = (currentCellIndex>>3);
                            sizeOfBlock += 1;
                            currentCellIndex = (((currentCellIndex>>3)+1)<<3);
                          }
                        else{
                            currentCellIndex += 1;
                          }
                      }
                    // If we are at the end of the sub group, move to next (otherwise we have consume a part of it)
                    while(iterChildCells != iterChildEndCells && (*iterChildCells)->getEndingIndex() <= currentCellIndex){
                        ++iterChildCells;
                        // Update morton index
                        if(iterChildCells != iterChildEndCells && currentCellIndex < (*iterChildCells)->getStartingIndex()){
                            currentCellIndex = (*iterChildCells)->getStartingIndex();
                          }
                      }
                  }

                // If group is full
                if(sizeOfBlock){
                    // Create a group
                    CellGroupClass*const newBlock = new CellGroupClass(currentBlockIndexes[0],
                        currentBlockIndexes[sizeOfBlock-1]+1,
                        sizeOfBlock);
                    // Init cells
                    for(int cellIdInBlock = 0; cellIdInBlock != sizeOfBlock ; ++cellIdInBlock){
                        newBlock->newCell(currentBlockIndexes[cellIdInBlock], cellIdInBlock);

                        SymbolCellClass& symbolic = newBlock->getSymbolic(cellIdInBlock);
                        symbolic.setMortonIndex(currentBlockIndexes[cellIdInBlock]);
                        FTreeCoordinate coord;
                        coord.setPositionFromMorton(currentBlockIndexes[cellIdInBlock]);
                        symbolic.setCoordinate(coord);
                        symbolic.setLevel(idxLevel);
                      }

                    // Keep the block
                    newBlock->declare_mine();
                    _cellBlocksPerLevel[idxLevel].push_back(newBlock);

                    sizeOfBlock = 0;
                  }
              }
          }
      }
    delete[] currentBlockIndexes;
  }
    /**
     * Sequential Constructor of GroupTree
     * used to construct a duplicated Ggroup tree on all processes
     * @param[in] in_treeHeight size of the tree
     * @param[in] in_boxWidth   bow witdh
     * @param[in] in_boxCenter  box center
     * @param[in] in__nbElementsPerBlock block size
     * @param[in] inParticlesContainer  an array of particles
     * @param[out] blockSizeAtEachLevel  box width at leaf level
     * @param[in] particlesAreSorted  True if the particle are sorted
     */
  template<class ParticleContainer>
  FGroupTree(const int in_treeHeight, const FReal inBoxWidth, const FPoint<FReal>& inBoxCenter,
             const int in_nbElementsPerBlock, ParticleContainer* inParticlesContainer,
             std::vector<std::vector<int>> & blockSizeAtEachLevel,
             const bool particlesAreSorted = false):
    _treeHeight(in_treeHeight),_nbElementsPerBlock(in_nbElementsPerBlock),_cellBlocksPerLevel(nullptr),
    boxCenter(inBoxCenter), boxCorner(inBoxCenter,-(inBoxWidth/2)), boxWidth(inBoxWidth),
    boxWidthAtLeafLevel(inBoxWidth/FReal(1<<(in_treeHeight-1)))
  {
    _cellBlocksPerLevel = new std::vector<CellGroupClass*>[_treeHeight];

    MortonIndex* currentBlockIndexes = new MortonIndex[_nbElementsPerBlock];
    // First we work at leaf level
    {
      // Build morton index for particles
      struct ParticleSortingStruct{
        FSize originalIndex;
        MortonIndex mindex;
      };
      // Convert position to morton index
      const FSize nbParticles = inParticlesContainer->getNbParticles();
      ParticleSortingStruct* particlesToSort = new ParticleSortingStruct[nbParticles];
      {
        const FReal* xpos = inParticlesContainer->getPositions()[0];
        const FReal* ypos = inParticlesContainer->getPositions()[1];
        const FReal* zpos = inParticlesContainer->getPositions()[2];

        for(FSize idxPart = 0 ; idxPart < nbParticles ; ++idxPart){
            const FTreeCoordinate host = FCoordinateComputer::GetCoordinateFromPositionAndCorner<FReal>(this->boxCorner, this->boxWidth,
                                                                                                        _treeHeight,
                                                                                                        FPoint<FReal>(xpos[idxPart], ypos[idxPart], zpos[idxPart]) );
        //    const MortonIndex particleIndex         = host.getMortonIndex();
            particlesToSort[idxPart].mindex         = host.getMortonIndex();
            particlesToSort[idxPart].originalIndex = idxPart;
          }
      }


      // Sort if needed
      if(particlesAreSorted == false){
          FQuickSort<ParticleSortingStruct, FSize>::QsOmp(particlesToSort, nbParticles, [](const ParticleSortingStruct& v1, const ParticleSortingStruct& v2){
              return v1.mindex <= v2.mindex;
            });
        }

      // Convert to block
      const int idxLevel = (_treeHeight - 1);
      int idxBlock = 0;
      FSize* nbParticlesPerLeaf = new FSize[_nbElementsPerBlock];
      FSize firstParticle = 0;
      // We need to proceed each group in sub level
      while(firstParticle != nbParticles){
          int sizeOfBlock = 0;
          FSize lastParticle = firstParticle;
          // Count until end of sub group is reached or we have enough cells
          while(sizeOfBlock < blockSizeAtEachLevel[_treeHeight-1][idxBlock] && lastParticle < nbParticles){
              if(sizeOfBlock == 0 || currentBlockIndexes[sizeOfBlock-1] != particlesToSort[lastParticle].mindex){
                  currentBlockIndexes[sizeOfBlock] = particlesToSort[lastParticle].mindex;
                  nbParticlesPerLeaf[sizeOfBlock]  = 1;
                  sizeOfBlock += 1;
                }
              else{
                  nbParticlesPerLeaf[sizeOfBlock-1] += 1;
                }
              lastParticle += 1;
            }
          while(lastParticle < nbParticles && currentBlockIndexes[sizeOfBlock-1] == particlesToSort[lastParticle].mindex){
              nbParticlesPerLeaf[sizeOfBlock-1] += 1;
              lastParticle += 1;
            }

          // Create a group
          CellGroupClass*const newBlock = new CellGroupClass(currentBlockIndexes[0],
              currentBlockIndexes[sizeOfBlock-1]+1,
              sizeOfBlock);
          FGroupOfParticles<FReal, NbSymbAttributes, NbAttributesPerParticle, AttributeClass>*const newParticleBlock = new FGroupOfParticles<FReal, NbSymbAttributes, NbAttributesPerParticle, AttributeClass>(currentBlockIndexes[0],
              currentBlockIndexes[sizeOfBlock-1]+1,
              sizeOfBlock, lastParticle-firstParticle);

          /////////////////////////  TO REMOVE ?? //////////////
//          #include <iostream>
//	  using namespace std;
//	  if(currentBlockIndexes[sizeOfBlock-1]+1 == 511)
//	    cout << "Suricate" << endl;
	  /////////////////////////////////////////////////////

	  // Init cells
	  size_t nbParticlesOffsetBeforeLeaf = 0;
	  FSize offsetParticles = firstParticle;
	  for(int cellIdInBlock = 0; cellIdInBlock != sizeOfBlock ; ++cellIdInBlock){
	      newBlock->newCell(currentBlockIndexes[cellIdInBlock], cellIdInBlock);

	      SymbolCellClass& symbolic = newBlock->getSymbolic(cellIdInBlock);
	      symbolic.setMortonIndex(currentBlockIndexes[cellIdInBlock]);
	      FTreeCoordinate coord;
	      coord.setPositionFromMorton(currentBlockIndexes[cellIdInBlock]);
	      symbolic.setCoordinate(coord);
	      symbolic.setLevel(idxLevel);

	      // Add leaf
	      nbParticlesOffsetBeforeLeaf = newParticleBlock->newLeaf(currentBlockIndexes[cellIdInBlock], cellIdInBlock,
								      nbParticlesPerLeaf[cellIdInBlock], nbParticlesOffsetBeforeLeaf);

	      BasicAttachedClass attachedLeaf = newParticleBlock->template getLeaf<BasicAttachedClass>(cellIdInBlock);
	      // Copy each particle from the original position
	      for(FSize idxPart = 0 ; idxPart < nbParticlesPerLeaf[cellIdInBlock] ; ++idxPart){
		  attachedLeaf.setParticle(idxPart, particlesToSort[idxPart + offsetParticles].originalIndex, inParticlesContainer);
		}
	      offsetParticles += nbParticlesPerLeaf[cellIdInBlock];
	    }

	  // Keep the block
	  newBlock->declare_mine();
	  _cellBlocksPerLevel[idxLevel].push_back(newBlock);
	  _particleBlocks.push_back(newParticleBlock);

	  sizeOfBlock = 0;
	  firstParticle = lastParticle;
	  ++idxBlock;
	}
      delete[] nbParticlesPerLeaf;
      delete[] particlesToSort;
    }
//    MPI_Barrier(MPI_COMM_WORLD);

    // For each level from heigth - 2 to 1
    for(int idxLevel = _treeHeight-2; idxLevel > 0 ; --idxLevel){
        CellGroupConstIterator iterChildCells = _cellBlocksPerLevel[idxLevel+1].begin();
        const CellGroupConstIterator iterChildEndCells = _cellBlocksPerLevel[idxLevel+1].end();

        // If lower level is empty or all blocks skiped stop here
        if(iterChildCells == iterChildEndCells){
            break;
          }

        MortonIndex currentCellIndex = (*iterChildCells)->getStartingIndex();
        int sizeOfBlock = 0;
        int idxBlock    = 0;

        // We need to proceed each group in sub level
        while(iterChildCells != iterChildEndCells){

            // Count until end of sub group is reached or we have enough cells
            while(sizeOfBlock < blockSizeAtEachLevel[idxLevel][idxBlock] && iterChildCells != iterChildEndCells ){

                if((sizeOfBlock == 0 || currentBlockIndexes[sizeOfBlock-1] != (currentCellIndex>>3))
                   && (*iterChildCells)->exists(currentCellIndex)){
                    currentBlockIndexes[sizeOfBlock] = (currentCellIndex>>3);
                    sizeOfBlock += 1;
                    currentCellIndex = (((currentCellIndex>>3)+1)<<3);
                  }
                else{
                    currentCellIndex += 1;
                  }
                // If we are at the end of the sub group, move to next
                while(iterChildCells != iterChildEndCells && (*iterChildCells)->getEndingIndex() <= currentCellIndex){
                    ++iterChildCells;
                    // Update morton index
                    if(iterChildCells != iterChildEndCells && currentCellIndex < (*iterChildCells)->getStartingIndex()){
                        currentCellIndex = (*iterChildCells)->getStartingIndex();
                      }
                  }
              }

            // If group is full
            if(sizeOfBlock == blockSizeAtEachLevel[idxLevel][idxBlock] || (sizeOfBlock && iterChildCells == iterChildEndCells)){ //NOTE la seconde partie va sûrement sauter, car la taille est pré-calculée
                // Create a group
                CellGroupClass*const newBlock = new CellGroupClass(currentBlockIndexes[0],
                    currentBlockIndexes[sizeOfBlock-1]+1,
                    sizeOfBlock);
                // Init cells
                for(int cellIdInBlock = 0; cellIdInBlock != sizeOfBlock ; ++cellIdInBlock){
                    newBlock->newCell(currentBlockIndexes[cellIdInBlock], cellIdInBlock);

                    SymbolCellClass& symbolic = newBlock->getSymbolic(cellIdInBlock);
                    symbolic.setMortonIndex(currentBlockIndexes[cellIdInBlock]);
                    FTreeCoordinate coord;
                    coord.setPositionFromMorton(currentBlockIndexes[cellIdInBlock]);
                    symbolic.setCoordinate(coord);
                    symbolic.setLevel(idxLevel);
                  }

                // Keep the block
                newBlock->declare_mine();
                _cellBlocksPerLevel[idxLevel].push_back(newBlock);

                sizeOfBlock = 0;
                ++idxBlock;
              }
          }
      }
    delete[] currentBlockIndexes;
  }


  /**
     * Minimal Constructor of GroupTree
     * @author benjamin.dufoyer@inria.fr
     * @param in__treeHeight size of the tree
     * @param in__nbElementsPerBlock block size
     * @param in_boxCenter  box center
     * @param in_boxCorner  box cornet
     * @param in_boxWidth   bow witdh
     * @param in_boxWidthAtLeafLevel  box width at leaf level
     */
  FGroupTree(
      int in__treeHeight,
      int in__nbElementsPerBlock,
      FPoint<FReal> in_boxCenter,
      FPoint<FReal> in_boxCorner,
      FReal in_boxWidth,
      FReal in_boxWidthAtLeafLevel
      ):
    _treeHeight(in__treeHeight),
    _nbElementsPerBlock(in__nbElementsPerBlock),
    boxCenter(in_boxCenter),
    boxCorner(in_boxCorner),
    boxWidth(in_boxWidth),
    boxWidthAtLeafLevel(in_boxWidthAtLeafLevel)
  {
    this->_cellBlocksPerLevel = new std::vector<CellGroupClass*>[_treeHeight];

  }


  /**
     * get_block_tree_instance return a new instance of FGroupTree from
     * a blocked linear tree
     * @author benjamin.dufoyer@inria.fr
     * @param  blocked_linear_tree blocked linear tree
     * @return new FGroupTree
     */
  template<
      class GroupCellSymbClass,
      class GroupCellUpClass,
      class GroupCellDownClass,
      class GroupContainerClass
      >
  static FGroupTree get_block_tree_instance(
      int in_tree_height,
      int in_block_size,
      FPoint<FReal> in_box_center,
      FReal in_box_width
      ){
    // Compute every information to initialise the group tree
    FPoint<FReal> box_corner = FPoint<FReal>(in_box_center, -in_box_width/2);
    FReal box_width_at_leaf_level = in_box_width/FReal( 1<< (in_tree_height-1));

    // Return a new instance of a empty group tree
    return FGroupTree<FReal,GroupCellSymbClass,GroupCellUpClass, GroupCellDownClass, GroupContainerClass, NbSymbAttributes, NbAttributesPerParticle, FReal>(
          in_tree_height
          ,in_block_size
          ,in_box_center
          ,box_corner
          ,in_box_width
          ,box_width_at_leaf_level);
  }


  /////////////////////////////////////////////////////////
  // Function to init group tree
  /////////////////////////////////////////////////////////

  /**
     * create_tree this function fill the tree from blocked_linear_tree
     * She build the group tree from the bottom
     * @author benjamin.dufoyer@inria.fr
     * @param  in_lin_tree  Blocked linear tree
     * @param  particles    vector where particle are stock,
     *                      they will be sort BEFORE calling this function
     */
  template<class Group_Linear_tree,
           class Particle_Container
           >
  void create_tree(Group_Linear_tree& in_lin_tree,
                   const Particle_Container& particles){
    MortonIndex in_left_limit = in_lin_tree.get_left_limit();
    // Creation of the leaf level and groups of particle
    std::vector<MortonIndex> current_block_indexes = create_leaf_level(in_lin_tree,particles);
    // Creation of every level of the tree
    create_block_nodes_level(
          current_block_indexes,
          in_left_limit);
  }

  /** This function dealloc the tree by deleting each block */
  ~FGroupTree(){
    for(int idxLevel = 0 ; idxLevel < _treeHeight ; ++idxLevel){
        std::vector<CellGroupClass*>& levelBlocks = _cellBlocksPerLevel[idxLevel];
        for (CellGroupClass* block: levelBlocks){
            delete block;
          }
      }
    delete[] _cellBlocksPerLevel;
    for (ParticleGroupClass* block: _particleBlocks){
        delete block;
      }
  }

  ////////////////////////////////////////////////////////
  // Lambda function to apply to all member
  /////////////////////////////////////////////////////////

  /**
   * @brief forEachLeaf iterate on the leaf and apply the function
   * @param function
   */
  template<class ParticlesAttachedClass>
  void forEachLeaf(std::function<void(ParticlesAttachedClass*)> function){
    for (ParticleGroupClass* block: _particleBlocks){
        block->forEachLeaf(function);
      }
  }

  /**
   * @brief forEachLeaf iterate on the leaf and apply the generic function
   * @param function
   */
  template<class F>
  void forEachLeaf(F&& function){
    for (ParticleGroupClass* block: _particleBlocks){
        block->forEachLeaf(function);
      }
  }

  /**
   * @brief forEachMyLeaf iterate on the leaf and apply the function
   * @param function
   */
  template<class ParticlesAttachedClass>
  void forEachMyLeaf(std::function<void(ParticlesAttachedClass*)> function){
    for (ParticleGroupClass* block: _particleBlocks){
        if(block->isMine())
          block->forEachLeaf(function);
      }
  }

  /**
   * @brief forEachCell iterate on the cell and apply the function
   * @param function
   */
  void forEachCell(std::function<void(SymbolCellClass*,PoleCellClass*,LocalCellClass*)> function){
    for(int idxLevel = 0 ; idxLevel < _treeHeight ; ++idxLevel){
        std::vector<CellGroupClass*>& levelBlocks = _cellBlocksPerLevel[idxLevel];
        for (CellGroupClass* block: levelBlocks){
            block->forEachCell(function);
          }
      }
  }

  void forEachMyCell(std::function<void(SymbolCellClass*,PoleCellClass*,LocalCellClass*)> function){
    for(int idxLevel = 0 ; idxLevel < _treeHeight ; ++idxLevel){
        std::vector<CellGroupClass*>& levelBlocks = _cellBlocksPerLevel[idxLevel];
        for (CellGroupClass* block: levelBlocks){
            if(block->isMine())
              block->forEachCell(function);
          }
      }
  }

  /**
   * @brief forEachLeaf iterate on the cell and apply the function
   * @param function
   */
  void forEachCellWithLevel(std::function<void(SymbolCellClass*,PoleCellClass*,LocalCellClass*,const int)> function){
    for(int idxLevel = 0 ; idxLevel < _treeHeight ; ++idxLevel){
        std::vector<CellGroupClass*>& levelBlocks = _cellBlocksPerLevel[idxLevel];
        for (CellGroupClass* block: levelBlocks){
            block->forEachCell(function, idxLevel);
          }
      }
  }

  void forEachMyCellWithLevel(std::function<void(SymbolCellClass*,PoleCellClass*,LocalCellClass*,const int)> function){
    for(int idxLevel = 0 ; idxLevel < _treeHeight ; ++idxLevel){
        std::vector<CellGroupClass*>& levelBlocks = _cellBlocksPerLevel[idxLevel];
        for (CellGroupClass* block: levelBlocks){
            if(block->isMine())
              block->forEachCell(function, idxLevel);
          }
      }
  }

  /**
   * @brief forEachLeaf iterate on the cells that i own and apply the function on the leaves
   *
   *  We iterate on the owned group and for for each group we iterate on its cells.
   * We obtain the leaf through the morton index of the cell and we apply the function on it.
   *
   * @param function function to apply on each leaf
   *
   */
  template<class ParticlesAttachedClass>
  void forEachCellLeaf(std::function<void(SymbolCellClass*,PoleCellClass*,LocalCellClass*,ParticlesAttachedClass*)> function){
    CellGroupIterator       iterCells    = _cellBlocksPerLevel[_treeHeight-1].begin();
    const CellGroupIterator iterEndCells = _cellBlocksPerLevel[_treeHeight-1].end();

    ParticleGroupIterator       iterLeaves    = _particleBlocks.begin();
    const ParticleGroupIterator iterEndLeaves = _particleBlocks.end();
   // int count = 0 ;
     //  Iterate on Cell group and Leaf group
   while(iterCells != iterEndCells && iterLeaves != iterEndLeaves){
       //
       //  Iterate all Cells inside the group
        (*iterCells)->forEachCell(
              [&](SymbolCellClass* symb,
              PoleCellClass* mult,
              LocalCellClass* loc)
        {
         //  The group of cells is mine --> Leaf exist
//        std::cout << " forEachCellLeaf " << (*iterCells)->isMine() << " morton " << symb->getMortonIndex()
//                  << "   leafIdx " << (*iterLeaves)->getLeafIndex(symb->getMortonIndex())
//                  << "   count  " <<count <<std::endl;
        if ((*iterCells)->isMine() ) {

              //  get leafindex (position in the group) from the Morton index
              const int leafIdx = (*iterLeaves)->getLeafIndex(symb->getMortonIndex());
              //  Leaf exists and we apply function on it
              FAssertLF(leafIdx != -1);
              ParticlesAttachedClass aLeaf = (*iterLeaves)->template getLeaf <ParticlesAttachedClass>(leafIdx);
              FAssertLF(aLeaf.isAttachedToSomething());
              function(symb, mult, loc, &aLeaf);
            }
   //       ++count;
        });
        ++iterCells;
        ++iterLeaves;
      }

    FAssertLF(iterCells == iterEndCells && iterLeaves == iterEndLeaves);
  }


//  template<class ParticlesAttachedClass>
//  void forEachCellMyLeaf(std::function<void(SymbolCellClass*,PoleCellClass*,LocalCellClass*,ParticlesAttachedClass*)> function){
//    CellGroupIterator iterCells = _cellBlocksPerLevel[_treeHeight-1].begin();

//    const CellGroupIterator iterEndCells = _cellBlocksPerLevel[_treeHeight-1].end();

//    ParticleGroupIterator iterLeaves = _particleBlocks.begin();
//    const ParticleGroupIterator iterEndLeaves = _particleBlocks.end();

//    while(iterCells != iterEndCells && iterLeaves != iterEndLeaves){
//        if((*iterCells)->isMine()){
//            (*iterCells)->forEachCell(
//                  [&](SymbolCellClass* symb,
//                  PoleCellClass*       mult,
//                  LocalCellClass*      loc)
//            {
//              const int leafIdx = (*iterLeaves)->getLeafIndex(symb->getMortonIndex());
//              FAssertLF(leafIdx != -1);
//              ParticlesAttachedClass aLeaf = (*iterLeaves)->template getLeaf <ParticlesAttachedClass>(leafIdx);
//              FAssertLF(aLeaf.isAttachedToSomething());
//              function(symb, mult, loc, &aLeaf);
//            });
//          }
//        ++iterCells;
//        ++iterLeaves;
//      }

//    FAssertLF(iterCells == iterEndCells && iterLeaves == iterEndLeaves);
//  }


  /** @brief, for statistic purpose, display each block with number of
   * cell, size of header, starting index, and ending index
   */
  void printInfoBlocks(){
    std::cout << "Group Tree information:\n";
    std::cout << "\t Group Size = " << _nbElementsPerBlock << "\n";
    std::cout << "\t Tree height = " << _treeHeight << "\n";
    for(int idxLevel = 1 ; idxLevel < _treeHeight ; ++idxLevel){
        std::vector<CellGroupClass*>& levelBlocks = _cellBlocksPerLevel[idxLevel];
        std::cout << "Level " << idxLevel << ", there are " << levelBlocks.size() << " groups.\n";
        int idxGroup = 0;
        for (const CellGroupClass* block: levelBlocks){
            std::cout << "\t Group " << (idxGroup++);
          //  std::cout << "\t local " << std::boolalpha << block->isMine();
            std::cout << "\t Size = " << block->getNumberOfCellsInBlock();
            std::cout << "\t Starting Index = " << block->getStartingIndex();
            std::cout << "\t Ending Index = " << block->getEndingIndex();
         //   std::cout << "\t Global index  = " << block->getIdxGlobal();
            std::cout << "\t Ratio of usage = " <<
                         float(block->getNumberOfCellsInBlock())/float(block->getEndingIndex()-block->getStartingIndex()) << "\n";
          }
      }

    std::cout << "There are " << _particleBlocks.size() << " leaf-groups.\n";
    int idxGroup = 0;
    FSize totalNbParticles = 0;
    for (const ParticleGroupClass* block: _particleBlocks){
        std::cout << "\t Group " << (idxGroup++);

        std::cout << "\t Size = " << block->getNumberOfLeavesInBlock();
        std::cout << "\t Starting Index = " << block->getStartingIndex();
        std::cout << "\t Ending Index = " << block->getEndingIndex();
        std::cout << "\t Nb Particles = " << block->getNbParticlesInGroup();
        std::cout << "\t Global index  = " << block->getIdxGlobal();
        std::cout << "\t Ratio of usage = " << float(block->getNumberOfLeavesInBlock())/float(block->getEndingIndex()-block->getStartingIndex()) << "\n";
        totalNbParticles += block->getNbParticlesInGroup();
      }
    std::cout << "There are " << totalNbParticles << " particles.\n";
  }

  /////////////////////////////////////////////////////////
  // Algorithm function
  /////////////////////////////////////////////////////////

  int getHeight() const {
    return _treeHeight;
  }

  CellGroupIterator cellsBegin(const int inLevel){
    FAssertLF(inLevel < _treeHeight);
    return _cellBlocksPerLevel[inLevel].begin();
  }

  CellGroupConstIterator cellsBegin(const int inLevel) const {
    FAssertLF(inLevel < _treeHeight);
    return _cellBlocksPerLevel[inLevel].begin();
  }

  CellGroupIterator cellsEnd(const int inLevel){
    FAssertLF(inLevel < _treeHeight);
    return _cellBlocksPerLevel[inLevel].end();
  }

  CellGroupConstIterator cellsEnd(const int inLevel) const {
    FAssertLF(inLevel < _treeHeight);
    return _cellBlocksPerLevel[inLevel].end();
  }

  int getNbCellGroupAtLevel(const int inLevel) const {
    FAssertLF(inLevel < _treeHeight);
    return int(_cellBlocksPerLevel[inLevel].size());
  }

  CellGroupClass* getCellGroup(const int inLevel, const int inIdx){
    FAssertLF(inLevel < _treeHeight);
    FAssertLF(inIdx < int(_cellBlocksPerLevel[inLevel].size()));
    return _cellBlocksPerLevel[inLevel][inIdx];
  }

  const int getNbElementsPerBlock() const{
    return this->_nbElementsPerBlock;
  }

  const CellGroupClass* getCellGroup(const int inLevel, const int inIdx) const {
    FAssertLF(inLevel < _treeHeight);
    FAssertLF(inIdx < int(_cellBlocksPerLevel[inLevel].size()));
    return _cellBlocksPerLevel[inLevel][inIdx];
  }

  ParticleGroupIterator leavesBegin(){
    return _particleBlocks.begin();
  }

  ParticleGroupConstIterator leavesBegin() const {
    return _particleBlocks.begin();
  }

  ParticleGroupIterator leavesEnd(){
    return _particleBlocks.end();
  }

  ParticleGroupConstIterator leavesEnd() const {
    return _particleBlocks.end();
  }

  int getNbParticleGroup() const {
    return int(_particleBlocks.size());
  }

  ParticleGroupClass* getParticleGroup(const int inIdx){
    FAssertLF(inIdx < int(_particleBlocks.size()));
    return _particleBlocks[inIdx];
  }

  const ParticleGroupClass* getParticleGroup(const int inIdx) const {
    FAssertLF(inIdx < int(_particleBlocks.size()));
    return _particleBlocks[inIdx];
  }

  const FPoint<FReal> getBoxCenter() const{
    return this->boxCenter;
  }

  const FReal getBoxWidth() const{
    return this->boxWidth;
  }

  std::size_t getTotalNbLeaf() {
    std::size_t nbLeaf = 0;
    for(int i = 0 ; i < this->getNbParticleGroup();++i){
        nbLeaf += this->_particleBlocks[i]->getNumberOfLeavesInBlock();
      }
    return nbLeaf;
  }
  /**
     * RESTRICTION : The array will be initialise BEFORE
     * RESTRICTION : The morton index of particle will be at _treeHeight
     * get_number_of_particle compute the total number of
     * particle on every leaf, he just fill the array nb_particles_per_leaf
     * @author benjamin.dufoyer@inria.fr
     * @param  container    container of particle
     */
  template<class particle_t>
  void get_number_of_particle(const std::vector<particle_t>& container,
                              std::vector<std::size_t>& nb_particles_per_leaf){
    FAssert(container.size() != 0);
    int current_idx = 0;
    std::size_t old_m_index   = container.front().morton_index;
    std::size_t current_m_idx = old_m_index;
    for(std::size_t i = 0 ; i < container.size(); ++i){
        current_m_idx = container[i].morton_index;
        if(current_m_idx == old_m_index){
            nb_particles_per_leaf[current_idx] += 1;
          } else {
            current_idx += 1;
            nb_particles_per_leaf[current_idx] += 1;
            old_m_index = current_m_idx;
          }
      }
  }

  /**
     * create_leaf_level create the leaf level of the
     * Group tree from a blocked linear tree
     * @author benjamin.dufoyer@inria.fr
     * @param  in_lin_tree  Blocked linear tree
     * @param  particles    container of particle, will be a std::vector
     */

  template<class Blocked_Linear_tree,
           class Particle_Container>
  std::vector<MortonIndex> create_leaf_level(Blocked_Linear_tree& in_lin_tree,
                                             Particle_Container& particles)
  {
    // set parametter for the function
    const int idxLevel = this->_treeHeight-1;
    const int nb_block = in_lin_tree.get_nb_block();
    const int block_size = this->_nbElementsPerBlock;
    std::size_t in_nb_leaf = in_lin_tree.get_nb_leaf();
    auto tree = in_lin_tree.get_tree();
    // alloc the vector for the current block index
    // get the number of particle per leaf
    std::vector<MortonIndex> current_block_indexes(this->_nbElementsPerBlock,0);
    std::vector<std::size_t> nb_particle_per_leaf(in_nb_leaf,0);
    this->get_number_of_particle(particles,nb_particle_per_leaf);
    // put the particle in the FP2PParticleContainer
    FP2PParticleContainer<FReal> particle_container;
    for(unsigned i = 0 ; i < particles.size() ; ++i){
        particle_container.push(particles[i].position(), particles[i].physicalValue());
      }

    std::size_t leaf_number = 0;
    std::size_t leaf_number_min = 0;

    // Create every block
    std::size_t idx_particules = 0;
    for(int n_block = 0 ; n_block < nb_block ; ++n_block){
        // Compute the morton index for the first and the
        // last cell of the block
        unsigned size_of_block = 0;
        while(size_of_block < (unsigned)block_size
              && leaf_number < in_nb_leaf)
          {
            current_block_indexes[size_of_block] = tree->data()[leaf_number].morton_index;
            leaf_number += 1;
            size_of_block += 1;
          }

        CellGroupClass*const new_block = new CellGroupClass(current_block_indexes[0],
            current_block_indexes[size_of_block-1]+1, //+1 is need by the class
            size_of_block);
        size_t current_nb_particle = 0;
        for(size_t i = 0 ; i < size_of_block ; ++i){
            current_nb_particle += nb_particle_per_leaf[leaf_number_min+i];
          }
        FGroupOfParticles<
            FReal,
            NbSymbAttributes,
            NbAttributesPerParticle,
            AttributeClass>*const new_particle_block
            = new FGroupOfParticles<
            FReal,
            NbSymbAttributes,
            NbAttributesPerParticle,
            AttributeClass>
            (current_block_indexes[0],
            current_block_indexes[size_of_block-1]+1,
            size_of_block,
            current_nb_particle);

        // Initialise each cell of the block
        size_t nb_particles_offset_before_leaf = 0;
        for(unsigned cell_id_in_block = 0;  cell_id_in_block < size_of_block; ++cell_id_in_block)
          {
            // Adding cell into leaf block
            new_block->newCell(
                  current_block_indexes[cell_id_in_block],
                  cell_id_in_block);

            // Fill symbolic information of the block
            SymbolCellClass& symbolic =
                new_block->getSymbolic(cell_id_in_block);
            symbolic.setMortonIndex(current_block_indexes[cell_id_in_block]);
            FTreeCoordinate coord;
            coord.setPositionFromMorton(current_block_indexes[cell_id_in_block]);
            symbolic.setCoordinate(coord);
            symbolic.setLevel(idxLevel);

            // Adding cell into particle blockCells

            nb_particles_offset_before_leaf =
                new_particle_block->newLeaf(
                  current_block_indexes[cell_id_in_block],
                  cell_id_in_block,
                  FSize(nb_particle_per_leaf[leaf_number_min+cell_id_in_block]),
                nb_particles_offset_before_leaf
                );


            BasicAttachedClass attached_leaf =
                new_particle_block->template getLeaf<BasicAttachedClass>(cell_id_in_block);

            // Adding particle
            for(size_t idxPart = 0 ; idxPart <   nb_particle_per_leaf[leaf_number_min+cell_id_in_block] ; ++idxPart ){
                attached_leaf.setParticle(
                      idxPart,
                      idx_particules,
                      //nb_particles_offset_before_leaf+idxPart,
                      &particle_container);
                ++idx_particules;
              }
            // Setting the offset to don't use particle twice
            //offset_particles += nb_particle_per_leaf[idx_nb_particle_in_block];
            //idx_nb_particle_in_block += 1;
            // cell_id_in_block += 1;
          }
        leaf_number_min = leaf_number;
        new_block->declare_mine();
        //Stock the block cell and the block particles
        _cellBlocksPerLevel[idxLevel].push_back(new_block);
        _particleBlocks.push_back(new_particle_block);
        size_of_block = 0;
      }
    return {current_block_indexes.begin(),current_block_indexes.end()};

  }

  /**
    * create_level create every level
    * It's juste a factorisation from the Beregenger constructor
    * @author benjamin.dufoyer@inria.fr
    * @param currentBlockIndexes block repartition at leaf level
    * to construct
    * @param inLeftLimite left limit of block of the current proc
    * this parameter is not used with the blocked_linear_tree, he is here
    * to have compatibility with old constructor
    */
  void create_block_nodes_level(std::vector<MortonIndex>& currentBlockIndexes,
                                MortonIndex inLeftLimite = -1
      ){
    // Cronstruct every level
    for(int idxLevel = _treeHeight-2; idxLevel > 0 ; --idxLevel){
        inLeftLimite = (inLeftLimite == -1 ? inLeftLimite : (inLeftLimite>>3));

        CellGroupConstIterator iterChildCells          = _cellBlocksPerLevel[idxLevel+1].begin();
        const CellGroupConstIterator iterChildEndCells = _cellBlocksPerLevel[idxLevel+1].end();

        // Skip blocks that do not respect limit
        while(iterChildCells != iterChildEndCells
              && ((*iterChildCells)->getEndingIndex()>>3) <= inLeftLimite){
            ++iterChildCells;
          }
        // If lower level is empty or all blocks skiped stop here
        if(iterChildCells == iterChildEndCells){
            break;
          }

        MortonIndex currentCellIndex = (*iterChildCells)->getStartingIndex();
        if((currentCellIndex>>3) <= inLeftLimite) currentCellIndex = ((inLeftLimite+1)<<3);
        int sizeOfBlock = 0;

        // We need to proceed each group in sub level
        while(iterChildCells != iterChildEndCells){
            // Count until end of sub group is reached or we have enough cells
            while(sizeOfBlock < _nbElementsPerBlock && iterChildCells != iterChildEndCells ){
                if((sizeOfBlock == 0 || currentBlockIndexes[sizeOfBlock-1] != (currentCellIndex>>3))
                   && (*iterChildCells)->exists(currentCellIndex)){
                    currentBlockIndexes[sizeOfBlock] = (currentCellIndex>>3);
                    sizeOfBlock += 1;
                    currentCellIndex = (((currentCellIndex>>3)+1)<<3);
                  }
                else{
                    currentCellIndex += 1;
                  }
                // If we are at the end of the sub group, move to next
                while(iterChildCells != iterChildEndCells && (*iterChildCells)->getEndingIndex() <= currentCellIndex){
                    ++iterChildCells;
                    // Update morton index
                    if(iterChildCells != iterChildEndCells && currentCellIndex < (*iterChildCells)->getStartingIndex()){
                        currentCellIndex = (*iterChildCells)->getStartingIndex();
                      }
                  }
              }

            // If group is full
            if(sizeOfBlock == _nbElementsPerBlock || (sizeOfBlock && iterChildCells == iterChildEndCells)){
                // Create a group
                CellGroupClass*const newBlock = new CellGroupClass(currentBlockIndexes[0],
                    currentBlockIndexes[sizeOfBlock-1]+1,
                    sizeOfBlock);
                // Init cells
                for(int cellIdInBlock = 0; cellIdInBlock != sizeOfBlock ; ++cellIdInBlock){
                    newBlock->newCell(currentBlockIndexes[cellIdInBlock], cellIdInBlock);

                    SymbolCellClass& symbolic = newBlock->getSymbolic(cellIdInBlock);
                    symbolic.setMortonIndex(currentBlockIndexes[cellIdInBlock]);
                    FTreeCoordinate coord;
                    coord.setPositionFromMorton(currentBlockIndexes[cellIdInBlock]);
                    symbolic.setCoordinate(coord);
                    symbolic.setLevel(idxLevel);
                  }
                newBlock->declare_mine();
                // Keep the block
                _cellBlocksPerLevel[idxLevel].push_back(newBlock);
                sizeOfBlock = 0;
              }
          }
      }
  }



  /**
    * This function add all LET block put in parameter
    * She put block in order according to idx_global
    * She detect if we are at leaf level and create particle group
    *
    * @author benjamin.dufoyer@inria.fr
    * @param  block_to_insert pair of symbolic information of cellGroup and Part
    *                         group
    * @param  level           The level where are adding LET group
    */
  template<class particle_symbolic_block_t,
           class cell_symbolic_block_t>
  void add_LET_block(
      std::pair<std::vector<cell_symbolic_block_t>,
      std::vector<particle_symbolic_block_t>>& block_to_insert,
      int                                                level
      ){
    // Check if we are at leaf level
    bool leaf_level = ( level == ( _treeHeight - 1 ) );
    // Bind the vector of the pair
    std::vector<cell_symbolic_block_t> cell_to_insert = block_to_insert.first;
    std::vector<particle_symbolic_block_t> particle_to_insert = block_to_insert.second;
    // If we are at leaf level
    if(leaf_level){
        // Check if we have the same number of symoblic information of cellBlock
        // and of particleBlock
        FAssert(cell_to_insert.size() == particle_to_insert.size());
      } else {
        // Else check if the particle block is empty
        FAssert(particle_to_insert.size() == 0);
      }
    // if we have no block to insert, we don't need to continue this function
    if(cell_to_insert.size() == 0)
      return;
    // Get my local minimum index global
    int min_idx_global = this->getCellGroup(level,0)->getIdxGlobal();
    // Allocate vector of new block
    std::vector<CellGroupClass*> vect_block(cell_to_insert.size());

    // Fill the vector of new block
    unsigned block_at_begin = 0;
    // iterate on every cell
    for(unsigned i = 0; i < cell_to_insert.size(); ++i){
        // create new cell
        vect_block[i] = new CellGroupClass(
              cell_to_insert[i].start_index ,
              cell_to_insert[i].end_index,
              (int)cell_to_insert[i].nb_leaf_in_block );
        // set the global index of the cell
        vect_block[i]->setIdxGlobal(cell_to_insert[i].idx_global_block);
        // if the global index is less than the local idex, we need to
        // insert
        // the block at the beginning of the tree
        if(cell_to_insert[i].idx_global_block < min_idx_global){
            ++block_at_begin;
          }
        // init each cell of the new block
        for(unsigned j = 0; j < cell_to_insert[i].m_idx_in_block.size(); ++j){
            vect_block[i]->newCell(cell_to_insert[i].m_idx_in_block[j],j);
          }
      }
    // Add block at beginning of the level
    _cellBlocksPerLevel[level].insert(
          _cellBlocksPerLevel[level].begin(),
          vect_block.begin(),
          vect_block.begin()+block_at_begin);

    // Add block a the end of the level
    _cellBlocksPerLevel[level].insert(
          _cellBlocksPerLevel[level].end(),
          vect_block.begin()+block_at_begin,
          vect_block.end());
    // if we are at the leaf level
    if(leaf_level ){
        // init of the vector of particle
        std::vector<ParticleGroupClass*> vect_particle(particle_to_insert.size());

        block_at_begin = 0;
        // iterate on every symbolic particle group
        for(unsigned i = 0 ; i < particle_to_insert.size(); ++i ){
            // create a new particle group
            vect_particle[i] = new ParticleGroupClass(
                  cell_to_insert[i].start_index ,
                  cell_to_insert[i].end_index,
                  (int)cell_to_insert[i].nb_leaf_in_block,
                  particle_to_insert[i].nb_particles);
            // set the global index of the new particle group
            vect_particle[i]-> setIdxGlobal(particle_to_insert[i].idx_global_block);
            // if the current idx global block have a idx global smaller than
            // the global index in local
            if(cell_to_insert[i].idx_global_block < min_idx_global){
                ++block_at_begin;
              }
            size_t offset = 0;
            // init all leaf of the current particle group
            for(int j = 0; j < cell_to_insert[i].nb_leaf_in_block; ++j){
                offset = vect_particle[i]->newLeaf(
                      cell_to_insert[i].m_idx_in_block[j],
                      j,
                      particle_to_insert[i].nb_particle_per_leaf[j],
                      offset);
              }
          }
        // Add block at beginning of the level
        _particleBlocks.insert(
              _particleBlocks.begin(),
              vect_particle.begin(),
              vect_particle.begin()+block_at_begin);
        // Add block a the end of the level
        _particleBlocks.insert(
              _particleBlocks.end(),
              vect_particle.begin()+block_at_begin,
              vect_particle.end());
      }
  }


#ifdef SCALFMM_USE_MPI
  /**
     * This function compute and add the local essential tree (LET) at
     * the level.
     * We compute interaction for the P2P(if needed) and M2L. We communicate
     * other proc to get the GroupOfCell needed for building the LET
     * @author benjamin.dufoyer@inria.fr
     * @param  group_linear_tree        The group linear tree
     * @param  level                    The level to build the LET
     * @param  dim                      The dimension of Coordinate
     */
  template<class GroupLinearTree>
  void create_LET_at_level(
      GroupLinearTree&    group_linear_tree,
      int&                level,
      MortonIndex&        gmin,
      MortonIndex&        gmax,
      MortonIndex&        lmin,
      MortonIndex&        lmax,
      int                 dim = 3
      ){
    // stock in the variable if we are at the leaf level
    bool leaf_level = (this->getHeight()-1 == level);
    // update the morton index
    if(!leaf_level){
        gmin = gmin >> 3;
        gmax = gmax >> 3;
      }
    const MortonIndex global_min_m_idx = gmin;
    const MortonIndex global_max_m_idx = gmax;
    // Compute min and max local morton index at the level needed
    if(this->getNbCellGroupAtLevel(level) > 0){
        lmin = this->getCellGroup(level,0)->getStartingIndex();
        lmax = this->getCellGroup(level,this->getNbCellGroupAtLevel(level)-1)->getEndingIndex()-1;
      } else {
        lmin = -1;
        lmax = -1;
      }
    const MortonIndex local_min_m_idx = lmin;
    const MortonIndex local_max_m_idx = lmax;

    // declare variable, needed because we fill it in a if case
    std::vector<MortonIndex> leaf_P2P;
    if(leaf_level){
        // IDEA : can be a task
        // This function compute the leaf needed by the P2P operation
        // This function return a vector with all leaf needed
        // The P2P interaction is only needed at leaf level
        leaf_P2P = dstr_grp_tree_builder::get_leaf_P2P_interaction(
              *this,
              global_min_m_idx,
              global_max_m_idx,
              local_min_m_idx,
              local_max_m_idx);
      }

    // IDEA can be a task
    // This function compute the leaf needed by the M2L operation
    // This function return a vector with all leaf needed
    // get leaf M2L
    std::vector<MortonIndex> leaf_M2L =
        dstr_grp_tree_builder::get_leaf_M2L_interaction_at_level(
          global_min_m_idx,
          global_max_m_idx,
          local_min_m_idx,
          local_max_m_idx,
          level,
          *this,
          dim);
    std::vector<MortonIndex> needed_leaf;
    if(leaf_level){
        // this function return the concatenation of the leaf for the P2P and
        // the leaf for the M2L
        needed_leaf = dstr_grp_tree_builder::concat_M2L_P2P(leaf_P2P,leaf_M2L);
      } else {
        // if it's not the leaf level, we juste need the M2L
        needed_leaf = leaf_M2L;
        group_linear_tree.update_index_particle_distribution(
              std::pair<MortonIndex,MortonIndex>(local_min_m_idx
                                                 ,local_max_m_idx)
              );
      }
    // free memory
    // this call swap the current vector to a empty vector
    std::vector<MortonIndex>().swap(leaf_P2P);
    std::vector<MortonIndex>().swap(leaf_M2L);

    std::vector<std::pair<MortonIndex,MortonIndex>> index_particle_distribution =
        group_linear_tree.get_index_particle_distribution();

    // Get the interaction matrix
    // matrix[2][nproc]
    // first line for Morton index to Send
    // second line for Morton index to Recv
    std::vector<std::vector<size_t>> global_matrix_interaction = dstr_grp_tree_builder::get_matrix_interaction(
          needed_leaf,
          index_particle_distribution,
          group_linear_tree.get_mpi_conf());

    // Send and get leaf
    // Auto is used to get the block more easly
    // it's a std::pair<std::vector<cell_symbolic_block>,std::vector<particle_symbolic_block>>
    // block_t is a struct define on FDistributedGroupTreeBuilder.hpp
    auto let_block =
        dstr_grp_tree_builder::send_get_symbolic_block_at_level(
          needed_leaf,
          global_matrix_interaction,
          *this,
          level,
          group_linear_tree.get_mpi_conf());

    // free needed leaf
    std::vector<MortonIndex>().swap(needed_leaf);
    // free interaction matrix
    std::vector<std::vector<size_t>>().swap(global_matrix_interaction);


    //add the LET block to the tree
    this->add_LET_block(
          let_block,
          level);
  }

  /**
     * this function create the local essential tree at every level requested
     * by this function
     * @author benjamin.dufoyer@inria.fr
     * @param  group_linear_tree        The group linear tree
     * @param  level_min                The minimum level to build the LET
     * @param  dim                      The dimension of coordinate
     */
  template<class GroupLinearTree>
  void create_LET(
      GroupLinearTree&    group_linear_tree,
      int                 level_min = 2,
      int                 dim = 3
      ){
    // get the particle distribution
    std::vector<std::pair<MortonIndex,MortonIndex>> index_particle_distribution =
        group_linear_tree.get_index_particle_distribution();

    // compute the min and the max global morton index at the level needed
    // Compute min and max global morton index at the level needed
    // This variable is used to put value in const
    MortonIndex gmin = index_particle_distribution.front().first;
    MortonIndex gmax = index_particle_distribution.back().second;
    MortonIndex lmin = this->getParticleGroup(0)->getStartingIndex();
    MortonIndex lmax = this->getParticleGroup((this->getNbParticleGroup()-1) )->getEndingIndex();
    // if we have more than 1 proc
    if( group_linear_tree.get_mpi_conf().comm.size() != 1 ){
        // compute the LET at every level
        for(int i = this->_treeHeight-1 ; i >= level_min ; --i){
            //        std::cout << "Start creating LET at " << i << std::endl;
            this->create_LET_at_level(group_linear_tree,i,gmin,gmax,lmin,lmax,dim);
          }
      }
    dstr_grp_tree_builder::send_get_block_M2M(
          *this,
          group_linear_tree.get_mpi_conf()
          );
  }
#endif
  /**
     * IDEA une factorisation peut être faite avec la fonction d'ajout du LET
     * This function allow to insert 1 block at a level needed
     * @author benjamin.dufoyer@inria.fr
     * @param  block_to_insert          symbolique information of the block
     * @param  list_m_idx               List of Morton Index
     * @param  level                    Level to insert
     * @param  nb_particle_per_leaf     [OPTIONNAL] number of particle per leaf
     */
  template<class info_symb_cell_t>
  void insert_block(
      info_symb_cell_t&           block_to_insert,
      std::vector<MortonIndex>&   list_m_idx,
      int                         level,
      std::vector<FSize>*         nb_particle_per_leaf = nullptr
      ){
    // Check if we already have this block
    for(int i = 0 ; i < this->getNbCellGroupAtLevel(level); ++i){
        auto* container = this->getCellGroup(level,i);
        // break the loop if the globalIdx is too big
        if(container->getIdxGlobal() > block_to_insert.idx_global_block ){
            break;
          }
        if(container->getIdxGlobal() == block_to_insert.idx_global_block ){
            return;
          }
      }
    // Check if we are at the leaf level
    bool leaf_level = ( level == ( _treeHeight - 1 ) );
    // create the new block
    CellGroupClass* new_block = new CellGroupClass(
          block_to_insert.start_index,
          block_to_insert.end_index,
          block_to_insert.nb_leaf_in_block);
    // set the global idx to the new block
    new_block->setIdxGlobal(block_to_insert.idx_global_block);
    // init all cell of the new block
    for(int i = 0 ; i < block_to_insert.nb_leaf_in_block ; ++i){
        new_block->newCell(list_m_idx[i],i);
      }
    // if we are at leaf level
    if(leaf_level){
        MortonIndex min_global_idx = 0;
        MortonIndex max_global_idx = 0;
        int idx_min = 0;
        int idx_max = 0;
        // seek the min morton index of my blocks
        for(int i = 0 ; i < this->getNbParticleGroup() ; ++i){
            if(this->getCellGroup(level,i)->isMine()){
                min_global_idx = this->getParticleGroup(i)->getIdxGlobal()-1;
                idx_min = i;
                break;
              }
          }
        // seek the max morton index of my blocks
        for(int i = this->getNbParticleGroup()-1 ; i >= 0 ; --i){
            if(this->getCellGroup(level,i)->isMine()){
                max_global_idx = this->getParticleGroup(i)->getIdxGlobal()+1;
                idx_max = i;
                break;
              }
          }
        // compute the number of particle of this block
        FSize nb_particle = 0;
        for(unsigned i = 0; i < nb_particle_per_leaf->size(); ++i){
            nb_particle += nb_particle_per_leaf->data()[i];
          }
        // create the particle group
        ParticleGroupClass* new_block_p = new ParticleGroupClass(
              block_to_insert.start_index ,
              block_to_insert.end_index,
              (int)block_to_insert.nb_leaf_in_block,
              nb_particle);
        // set the global index of the particle group
        new_block_p->setIdxGlobal((int)min_global_idx);
        std::size_t offset = 0;
        // create all leaf of the particle group
        for(int  i = 0 ; i < block_to_insert.nb_leaf_in_block ; ++i){
            offset = new_block_p->newLeaf(
                  list_m_idx[i],
                  i,
                  nb_particle_per_leaf->data()[i],
                  offset
                  );
          }
        // insert the particle group at the good place
        if(this->getParticleGroup(idx_min)->getStartingIndex() > block_to_insert.start_index){
            new_block_p->setIdxGlobal((int)min_global_idx);
            _particleBlocks.insert(
                  _particleBlocks.begin()+idx_min,
                  new_block_p
                  );
          } else {
            new_block_p->setIdxGlobal((int)max_global_idx);
            _particleBlocks.insert(
                  _particleBlocks.begin()+idx_max+1,
                  new_block_p
                  );
          }

      }
    // if we need to put the new block at first
    // if we already have a block at this level
    if(this->getNbCellGroupAtLevel(level) > 0) {
        if(this->getCellGroup(level,0)->getIdxGlobal() > block_to_insert.idx_global_block){
            _cellBlocksPerLevel[level].insert(
                  _cellBlocksPerLevel[level].begin(),
                  new_block);
            return;
          }
        // if we don't have block at this level
      } else {
        _cellBlocksPerLevel[level].insert(
              _cellBlocksPerLevel[level].begin(),
              new_block);
        return;
      }
    // else find the place of the block
    // iterate on every block
    for(int idx_block = 0 ; idx_block < this->getNbCellGroupAtLevel(level) ; ++idx_block){
        auto* container = this->getCellGroup(level,idx_block);
        // if the block i want to insert is already here
        if(container->getIdxGlobal() == block_to_insert.idx_global_block ){
            return;
          }
        if(container->getIdxGlobal() >  block_to_insert.idx_global_block ){
            _cellBlocksPerLevel[level].insert(
                  _cellBlocksPerLevel[level].begin()+idx_block,
                  new_block);
            return;
          }
      }
    _cellBlocksPerLevel[level].insert(
          _cellBlocksPerLevel[level].end(),
          new_block);
  }

};




#endif // FGROUPTREE_HPP
