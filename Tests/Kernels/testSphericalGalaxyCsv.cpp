// See LICENCE file at project root

#include <iostream>

#include <cstdio>
#include <cstdlib>

#include "Utils/FTic.hpp"
#include "Utils/FParameters.hpp"

#include "Containers/FOctree.hpp"
#include "Containers/FVector.hpp"

#include "Core/FFmmAlgorithm.hpp"
#include "Core/FFmmAlgorithmThread.hpp"

#include "Kernels/Spherical/FSphericalKernel.hpp"
#include "Kernels/Spherical/FSphericalCell.hpp"

#include "Files/FTreeCsvSaver.hpp"
#include "Files/FFmaGenericLoader.hpp"
#include "Arranger/FOctreeArranger.hpp"

#include "Components/FSimpleLeaf.hpp"
#include "Components/FParticleType.hpp"

#include "Kernels/P2P/FP2PParticleContainer.hpp"

#include "Utils/FParameterNames.hpp"
#include "Arranger/FAbstractMover.hpp"


template <class FReal>
class VelocityContainer : public FP2PParticleContainer<FReal> {
    typedef FP2PParticleContainer<FReal> Parent;
    FVector<FPoint<FReal>> velocities;

public:
    template<typename... Args>
    void push(const FPoint<FReal>& inParticlePosition, const FPoint<FReal>& velocity, Args... args){
        Parent::push(inParticlePosition, args... );
        velocities.push(velocity);
    }

    const FVector<FPoint<FReal>>& getVelocities() const{
        return velocities;
    }

    FVector<FPoint<FReal>>& getVelocities() {
        return velocities;
    }

    void fillToCsv(const FSize partIdx, FReal values[4]) const {
        values[0] = Parent::getPositions()[0][partIdx];
        values[1] = Parent::getPositions()[1][partIdx];
        values[2] = Parent::getPositions()[2][partIdx];
        values[3] = Parent::getPotentials()[partIdx];
    }
};


template <class FReal>
class GalaxyLoader : public FFmaGenericLoader<FReal> {
public:
    GalaxyLoader(const char* const filename) : FFmaGenericLoader<FReal>(filename) {
    }

    void fillParticle(FPoint<FReal>* position, FReal* physivalValue, FPoint<FReal>* velocity){
        FReal x,y,z,data, vx, vy, vz;
        *(FFmaGenericLoader<FReal>::file) >> x >> y >> z >> data >> vx >> vy >> vz;
        position->setPosition(x,y,z);
        *physivalValue = (data);
        velocity->setPosition(vx,vy,vz);
    }
};

template<class FReal, class OctreeClass>
class GalaxyMover : public FAbstractMover<FReal, OctreeClass, VelocityContainer<FReal>>{
private:
    VelocityContainer<FReal> toStoreRemovedParts;

public:
    GalaxyMover() {
    }

    virtual ~GalaxyMover(){
    }

    /** To get the position of the particle at idx idxPart in leaf lf */
    void getParticlePosition(VelocityContainer<FReal>* lf, const FSize idxPart, FPoint<FReal>* particlePos){
        (*particlePos) = FPoint<FReal>(lf->getPositions()[0][idxPart],lf->getPositions()[1][idxPart],lf->getPositions()[2][idxPart]);
    }

    /** Remove a particle but keep it to reinsert it later*/
    void removeFromLeafAndKeep(VelocityContainer<FReal>* lf, const FPoint<FReal>& particlePos, const FSize idxPart,FParticleType /*type*/){
        std::array<typename VelocityContainer<FReal>::AttributesClass, VelocityContainer<FReal>::NbAttributes> particleValues;
        for(int idxAttr = 0 ; idxAttr < VelocityContainer<FReal>::NbAttributes ; ++idxAttr){
            particleValues[idxAttr] = lf->getAttribute(idxAttr)[idxPart];
        }

        toStoreRemovedParts.push(particlePos,lf->getVelocities()[idxPart],particleValues);

        lf->getVelocities().removeOne(idxPart);
        lf->removeParticles(&idxPart,1);
    }

    /** Reinsert the previously saved particles */
    void insertAllParticles(OctreeClass* tree){
        std::array<typename VelocityContainer<FReal>::AttributesClass, VelocityContainer<FReal>::NbAttributes> particleValues;

        for(FSize idxToInsert = 0; idxToInsert<toStoreRemovedParts.getNbParticles() ; ++idxToInsert){
            for(int idxAttr = 0 ; idxAttr < VelocityContainer<FReal>::NbAttributes ; ++idxAttr){
                particleValues[idxAttr] = toStoreRemovedParts.getAttribute(idxAttr)[idxToInsert];
            }
            const FPoint<FReal> particlePos(toStoreRemovedParts.getPositions()[0][idxToInsert],
                                     toStoreRemovedParts.getPositions()[1][idxToInsert],
                                     toStoreRemovedParts.getPositions()[2][idxToInsert]);

            tree->insert(particlePos, toStoreRemovedParts.getVelocities()[idxToInsert], particleValues);
        }

        toStoreRemovedParts.clear();
        toStoreRemovedParts.getVelocities().clear();
    }
};


// Simply create particles and try the kernels
int main(int argc, char ** argv){
    FHelpDescribeAndExit(argc, argv,
                         "Run a Spherical Harmonic (Old Implementation) FMM kernel with several time step.",
                         FParameterDefinitions::InputFile, FParameterDefinitions::OctreeHeight,
                         FParameterDefinitions::OctreeSubHeight, FParameterDefinitions::SHDevelopment,
                         FParameterDefinitions::DeltaT, FParameterDefinitions::OutputFile);

    typedef double FReal;
    typedef FSphericalCell<FReal>          CellClass;
    typedef VelocityContainer<FReal>  ContainerClass;

    typedef FSimpleLeaf<FReal, ContainerClass >                     LeafClass;
    typedef FOctree<FReal, CellClass, ContainerClass , LeafClass >  OctreeClass;
    typedef FSphericalKernel<FReal, CellClass, ContainerClass >   KernelClass;

    typedef FFmmAlgorithmThread<OctreeClass,  CellClass, ContainerClass, KernelClass, LeafClass > FmmClass;

    typedef GalaxyMover<FReal, OctreeClass> MoverClass;
    typedef FOctreeArranger<FReal,OctreeClass, ContainerClass, MoverClass> ArrangerClass;
    ///////////////////////What we do/////////////////////////////
    std::cout << ">> This executable has to be used to test Spherical algorithm.\n";
    //////////////////////////////////////////////////////////////

    const int NbLevels = FParameters::getValue(argc,argv,FParameterDefinitions::OctreeHeight.options, 6);
    const int SizeSubLevels = FParameters::getValue(argc,argv,FParameterDefinitions::OctreeSubHeight.options, 3);
    const FReal DT          = FParameters::getValue(argc,argv,FParameterDefinitions::DeltaT.options, FReal(0.1));
    const int DevP          = FParameters::getValue(argc,argv,FParameterDefinitions::SHDevelopment.options, 5);

    FSphericalCell<FReal>::Init(DevP);

    GalaxyLoader<FReal> loader(FParameters::getStr(argc,argv,FParameterDefinitions::InputFile.options, "../Data/galaxy.fma.tmp"));

    // -----------------------------------------------------

    OctreeClass tree(NbLevels, SizeSubLevels, loader.getBoxWidth(), loader.getCenterOfBox());

    // -----------------------------------------------------

    std::cout << "Creating & Inserting " << loader.getNumberOfParticles() << " particles ..." << std::endl;
    std::cout << "\tHeight : " << NbLevels << " \t sub-height : " << SizeSubLevels << std::endl;

    {
        FPoint<FReal> position, velocity;
        FReal physicalValue;

        for(FSize idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
            loader.fillParticle(&position, &physicalValue, &velocity);
            tree.insert(position, velocity, physicalValue);
        }
    }

    // -----------------------------------------------------

    KernelClass kernels( DevP, NbLevels, loader.getBoxWidth(), loader.getCenterOfBox());
    FmmClass algo( &tree, &kernels);
    ArrangerClass arranger(&tree);
    FTreeCsvSaver<FReal, OctreeClass, ContainerClass> saver(FParameters::getStr(argc,argv,FParameterDefinitions::OutputFile.options, "/tmp/test%d.csv"));

    for(int idx = 0; idx < 100 ; ++idx){
        algo.execute();
        { // update velocity and position
            OctreeClass::Iterator octreeIterator(&tree);
            octreeIterator.gotoBottomLeft();
            do{
                kernels.computeVelocity(octreeIterator.getCurrentListTargets(), DT);
                kernels.updatePosition(octreeIterator.getCurrentListTargets(), DT);
            } while(octreeIterator.moveRight());
        }
        // update tree and vtk
        arranger.rearrange();
        saver.exportTree(&tree);
    }

    // -----------------------------------------------------

    return 0;
}
