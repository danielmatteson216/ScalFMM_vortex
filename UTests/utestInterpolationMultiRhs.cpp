// See LICENCE file at project root

// ==== CMAKE =====
// @FUSE_BLAS
// @FUSE_FFT
// @SCALFMM_PRIVATE
// ==============
#include <array>

#include "ScalFmmConfig.h"
#include "Utils/FGlobal.hpp"

#include "Containers/FOctree.hpp"

#include "Files/FFmaGenericLoader.hpp"

#include "Core/FFmmAlgorithmThread.hpp"
#include "Core/FFmmAlgorithm.hpp"

#include "FUTester.hpp"

#include "Components/FSimpleLeaf.hpp"


#include "Kernels/Chebyshev/FChebCell.hpp"
#include "Kernels/Interpolation/FInterpMatrixKernel.hpp"
#include "Kernels/Chebyshev/FChebKernel.hpp"
#include "Kernels/Chebyshev/FChebSymKernel.hpp"

#include "Kernels/Uniform/FUnifCell.hpp"
#include "Kernels/Interpolation/FInterpMatrixKernel.hpp"
#include "Kernels/Uniform/FUnifKernel.hpp"


#include "Kernels/P2P/FP2PParticleContainerIndexed.hpp"
/*
  In this test we compare the spherical FMM results and the direct results.
 */


/** the test class
 *
 */
class TestInterpolationKernel : public FUTester<TestInterpolationKernel> {

    ///////////////////////////////////////////////////////////
    // The tests!
    ///////////////////////////////////////////////////////////

    template <class FReal, class CellClass, class ContainerClass, class KernelClass, class MatrixKernelClass,
              class LeafClass, class OctreeClass, class FmmClass, const int NVals>
    void RunTest()	{
        // Warning in make test the exec dir it Build/UTests
        // Load particles
        //
        // Load particles
        //
        if(sizeof(FReal) == sizeof(float) ) {
            std::cerr << "No input data available for Float "<< std::endl;
            exit(EXIT_FAILURE);
        }
        const std::string parFile( (sizeof(FReal) == sizeof(float))?
                                       "Test/DirectFloat.bfma":
                                       "UTest/DirectDouble.bfma");
        //
        std::string filename(SCALFMMDataPath+parFile);
        //
        FFmaGenericLoader<FReal> loader(filename);
        if(!loader.isOpen()){
            Print("Cannot open particles file.");
            uassert(false);
            return;
        }
        Print("Number of particles:");
        Print(loader.getNumberOfParticles());

        const int NbLevels        = 4;
        const int SizeSubLevels = 2;

        // std::cout << "\nInterpolation FMM (ORDER="<< ORDER << ") ... " << std::endl;

        // Create Matrix Kernel
        const MatrixKernelClass MatrixKernel; // FUKernelTester is only designed to work with 1/R, i.e. matrix kernel ctor takes no argument.
        //
        FSize nbParticles = loader.getNumberOfParticles() ;
        FmaRWParticle<FReal, 8,8>* const particles = new FmaRWParticle<FReal, 8,8>[nbParticles];

        loader.fillParticle(particles,nbParticles);

        // Create octree
        OctreeClass tree(NbLevels, SizeSubLevels, loader.getBoxWidth(), loader.getCenterOfBox());
        // Insert particle in the tree
        for(FSize idxPart = 0 ; idxPart < nbParticles ; ++idxPart){
            // Convert FReal[NVALS] to std::array<FReal,NVALS>
            std::array<FReal, (1+4*1)*NVals> physicalState;
            for(int idxVals = 0 ; idxVals < NVals ; ++idxVals){
                physicalState[0*NVals+idxVals]= particles[idxPart].getPhysicalValue();
                physicalState[1*NVals+idxVals]=0.0;
                physicalState[2*NVals+idxVals]=0.0;
                physicalState[3*NVals+idxVals]=0.0;
                physicalState[4*NVals+idxVals]=0.0;
            }
            // put in tree
            tree.insert(particles[idxPart].getPosition(), idxPart, physicalState);
        }


        // Run FMM
        Print("Fmm...");
        KernelClass kernels(NbLevels, loader.getBoxWidth(), loader.getCenterOfBox(),&MatrixKernel);
        FmmClass algo(&tree,&kernels);
        algo.execute();
        //
        FReal energy= 0.0 , energyD = 0.0 ;
        for(FSize idx = 0 ; idx < loader.getNumberOfParticles()  ; ++idx){
            energyD +=  particles[idx].getPotential()*particles[idx].getPhysicalValue() ;
        }
        //
        // Compare
        Print("Compute Diff...");
        FMath::FAccurater<FReal> potentialDiff[NVals];
        FMath::FAccurater<FReal> fx, fy, fz;
        {
            tree.forEachLeaf([&](LeafClass* leaf){
                //
                for(int idxVals = 0 ; idxVals < NVals ; ++idxVals){
                    const FReal* const physicalValues = leaf->getTargets()->getPhysicalValues(idxVals);
                    const FReal*const potentials = leaf->getTargets()->getPotentials(idxVals);
                    const FReal*const forcesX = leaf->getTargets()->getForcesX(idxVals);
                    const FReal*const forcesY = leaf->getTargets()->getForcesY(idxVals);
                    const FReal*const forcesZ = leaf->getTargets()->getForcesZ(idxVals);
                    const FSize nbParticlesInLeaf = leaf->getTargets()->getNbParticles();
                    const FVector<FSize>& indexes = leaf->getTargets()->getIndexes();

                    for(FSize idxPart = 0 ; idxPart < nbParticlesInLeaf ; ++idxPart){

                        const FSize indexPartOrig = indexes[idxPart];
                        //					std::cout << " index "<< indexPartOrig << "   "  << particles[indexPartOrig].getPotential() << "   " << potentials[idxPart] << std::endl;
                        potentialDiff[idxVals].add(particles[indexPartOrig].getPotential(),potentials[idxPart]);
                        //
                        fx.add(particles[indexPartOrig].getForces()[0],forcesX[idxPart]);
                        fy.add(particles[indexPartOrig].getForces()[1],forcesY[idxPart]);
                        fz.add(particles[indexPartOrig].getForces()[2],forcesZ[idxPart]);
                        //

                        energy   += potentials[idxPart]*physicalValues[idxPart];
                    }

                }
            });

        }
        delete[] particles;
        energy /=NVals;
        // Print for information
        double errorPotRL2=0.0, errorPotRMS=0.0;
        Print("Potential diff is = ");
        for(int idxVals = 0 ; idxVals < NVals ; ++idxVals){
            printf("   Charge: %d\n",		idxVals);
            printf("         Pot L2Norm     %e\n",potentialDiff[idxVals].getL2Norm());
            printf("         Pot RL2Norm   %e\n",potentialDiff[idxVals].getRelativeL2Norm());
            printf("         Pot RMSError   %e\n",potentialDiff[idxVals].getRMSError());
            errorPotRL2 = std::max(errorPotRL2, potentialDiff[idxVals].getRelativeL2Norm());
            errorPotRMS = std::max(errorPotRMS, potentialDiff[idxVals].getRMSError());
        }
        Print("Fx diff is = ");
        printf("         Fx L2Norm     %e\n",fx.getL2Norm());
        printf("         Fx RL2Norm   %e\n",fx.getRelativeL2Norm());
        printf("         Fx RMSError   %e\n",fx.getRMSError());
        Print("Fy diff is = ");
        printf("        Fy L2Norm     %e\n",fy.getL2Norm());
        printf("        Fy RL2Norm   %e\n",fy.getRelativeL2Norm());
        printf("        Fy RMSError   %e\n",fy.getRMSError());
        Print("Fz diff is = ");
        printf("        Fz L2Norm     %e\n",fz.getL2Norm());
        printf("        Fz RL2Norm   %e\n",fz.getRelativeL2Norm());
        printf("        Fz RMSError   %e\n",fz.getRMSError());
        FReal L2error = (fx.getRelativeL2Norm()*fx.getRelativeL2Norm() + fy.getRelativeL2Norm()*fy.getRelativeL2Norm()  + fz.getRelativeL2Norm() *fz.getRelativeL2Norm()  );
        printf(" Total L2 Force Error= %e\n",FMath::Sqrt(L2error)) ;
        printf("  Energy Error  =   %.12e\n",FMath::Abs(energy-energyD));
        printf("  Energy FMM    =   %.12e\n",FMath::Abs(energy));
        printf("  Energy DIRECT =   %.12e\n",FMath::Abs(energyD));

        // Assert
        const FReal MaximumDiffPotential = FReal(9e-3);
        const FReal MaximumDiffForces     = FReal(9e-2);

        Print("Test1 - Error Relative L2 norm Potential ");
        uassert(errorPotRL2 < MaximumDiffPotential);    //1
        Print("Test2 - Error RMS L2 norm Potential ");
        uassert(errorPotRMS< MaximumDiffPotential);  //2
        Print("Test3 - Error Relative L2 norm FX ");
        uassert(fx.getRelativeL2Norm()  < MaximumDiffForces);                       //3
        Print("Test4 - Error RMS L2 norm FX ");
        uassert(fx.getRMSError() < MaximumDiffForces);                      //4
        Print("Test5 - Error Relative L2 norm FY ");
        uassert(fy.getRelativeL2Norm()  < MaximumDiffForces);                       //5
        Print("Test6 - Error RMS L2 norm FY ");
        uassert(fy.getRMSError() < MaximumDiffForces);                      //6
        Print("Test7 - Error Relative L2 norm FZ ");
        uassert(fz.getRelativeL2Norm()  < MaximumDiffForces);                      //8
        Print("Test8 - Error RMS L2 norm FZ ");
        uassert(fz.getRMSError() < MaximumDiffForces);                                           //8
        Print("Test9 - Error Relative L2 norm F ");
        uassert(L2error              < MaximumDiffForces);                                            //9   Total Force
        Print("Test10 - Relative error Energy ");
        uassert(FMath::Abs(energy-energyD) /energyD< MaximumDiffPotential);                     //10  Total Energy


        // Compute multipole local rhs diff
        FMath::FAccurater<FReal> localDiff;
        FMath::FAccurater<FReal> multiPoleDiff;
        tree.forEachCell([&](CellClass* cell){
            for( int idxRhs = 1 ; idxRhs < NVals ; ++idxRhs){
                localDiff.add(cell->getLocalExpansionData().get(0), cell->getLocalExpansionData().get(idxRhs), cell->getLocalExpansionData().getVectorSize());
                multiPoleDiff.add(cell->getMultipoleData().get(0), cell->getMultipoleData().get(idxRhs), cell->getMultipoleData().getVectorSize());
            }
        });
        Print("Local diff is = ");
        Print(localDiff.getL2Norm());
        Print(localDiff.getInfNorm());
        Print("Multipole diff is = ");
        Print(multiPoleDiff.getL2Norm());
        Print(multiPoleDiff.getInfNorm());

        uassert(localDiff.getL2Norm()  < 1e-10);
        uassert(localDiff.getInfNorm() < 1e-10);
        uassert(multiPoleDiff.getL2Norm()  < 1e-10);
        uassert(multiPoleDiff.getInfNorm() < 1e-10);
    }

    /** If memstas is running print the memory used */
    void PostTest() {
        if( FMemStats::controler.isUsed() ){
            std::cout << "Memory used at the end " << FMemStats::controler.getCurrentAllocated()
                      << " Bytes (" << FMemStats::controler.getCurrentAllocatedMB() << "MB)\n";
            std::cout << "Max memory used " << FMemStats::controler.getMaxAllocated()
                      << " Bytes (" << FMemStats::controler.getMaxAllocatedMB() << "MB)\n";
            std::cout << "Total memory used " << FMemStats::controler.getTotalAllocated()
                      << " Bytes (" << FMemStats::controler.getTotalAllocatedMB() << "MB)\n";
        }
    }


    ///////////////////////////////////////////////////////////
    // Set the tests!
    ///////////////////////////////////////////////////////////


    /** TestUnifKernel */
    void TestUnifKernel(){
        typedef double FReal;
        const int NVals = 3;
        const unsigned int ORDER = 6 ;
        // run test
        typedef FInterpMatrixKernelR<FReal> MatrixKernelClass;


        typedef FP2PParticleContainerIndexed<FReal,1,1,NVals> ContainerClass;
        typedef FSimpleLeaf<FReal, ContainerClass >  LeafClass;
        typedef FUnifCell<FReal,ORDER,1,1,NVals> CellClass;
        typedef FOctree<FReal, CellClass,ContainerClass,LeafClass> OctreeClass;
        typedef FUnifKernel<FReal,CellClass,ContainerClass,MatrixKernelClass,ORDER,NVals> KernelClass;
        typedef FFmmAlgorithm<OctreeClass,CellClass,ContainerClass,KernelClass,LeafClass> FmmClass;

        RunTest<FReal,CellClass,ContainerClass,KernelClass,MatrixKernelClass,LeafClass,OctreeClass,FmmClass, NVals>();
    }

    /** TestChebSymKernel */
    void TestChebSymKernel(){
        typedef double FReal;
        const int NVals = 3;
        const unsigned int ORDER = 6;
        typedef FP2PParticleContainerIndexed<FReal,1,1,NVals> ContainerClass;
        typedef FSimpleLeaf<FReal, ContainerClass> LeafClass;
        typedef FInterpMatrixKernelR<FReal> MatrixKernelClass;
        typedef FChebCell<FReal,ORDER, 1, 1, NVals> CellClass;
        typedef FOctree<FReal, CellClass,ContainerClass,LeafClass> OctreeClass;
        typedef FChebSymKernel<FReal,CellClass,ContainerClass,MatrixKernelClass,ORDER, NVals> KernelClass;
        typedef FFmmAlgorithm<OctreeClass,CellClass,ContainerClass,KernelClass,LeafClass> FmmClass;
        // run test
        RunTest<FReal,CellClass,ContainerClass,KernelClass,MatrixKernelClass,LeafClass,OctreeClass,FmmClass, NVals>();
    }

    /** TestChebKernel */
    void TestChebKernel(){
        typedef double FReal;
        const int NVals = 3;
        const unsigned int ORDER = 6;
        typedef FP2PParticleContainerIndexed<FReal,1,1,NVals> ContainerClass;
        typedef FSimpleLeaf<FReal, ContainerClass> LeafClass;
        typedef FInterpMatrixKernelR<FReal> MatrixKernelClass;
        typedef FChebCell<FReal,ORDER, 1, 1, NVals> CellClass;
        typedef FOctree<FReal, CellClass,ContainerClass,LeafClass> OctreeClass;
        typedef FChebKernel<FReal,CellClass,ContainerClass,MatrixKernelClass,ORDER, NVals> KernelClass;
        typedef FFmmAlgorithm<OctreeClass,CellClass,ContainerClass,KernelClass,LeafClass> FmmClass;
        // run test
        RunTest<FReal,CellClass,ContainerClass,KernelClass,MatrixKernelClass,LeafClass,OctreeClass,FmmClass, NVals>();
    }

    ///////////////////////////////////////////////////////////
    // Set the tests!
    ///////////////////////////////////////////////////////////

    /** set test */
    void SetTests(){

        //AddTest(&TestInterpolationKernel::TestUnifKernel,"Test Lagrange/Uniform grid FMM");
        AddTest(&TestInterpolationKernel::TestChebSymKernel,"Test Symmetric Chebyshev Kernel with 16 small SVDs and symmetries");
        AddTest(&TestInterpolationKernel::TestChebKernel,"Test Chebyshev Kernel with 1 large SVD");

    }
};


// You must do this
TestClass(TestInterpolationKernel)
