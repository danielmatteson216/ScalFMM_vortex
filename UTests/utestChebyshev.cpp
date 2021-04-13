// See LICENCE file at project root

// ==== CMAKE =====
// @FUSE_BLAS
// ================

#include "FUKernelTester.hpp"

#include "Components/FSimpleLeaf.hpp"

#include "Kernels/Chebyshev/FChebCell.hpp"
#include "Kernels/Interpolation/FInterpMatrixKernel.hpp"
#include "Kernels/Chebyshev/FChebDenseKernel.hpp"
#include "Kernels/Chebyshev/FChebKernel.hpp"
#include "Kernels/Chebyshev/FChebSymKernel.hpp"

#include "Kernels/P2P/FP2PParticleContainerIndexed.hpp"
/*
  In this test we compare the Chebyschev fmm results and the direct results.
*/

/** the test class
 *
 */
class TestChebyshevDirect : public FUKernelTester<TestChebyshevDirect> {

  ///////////////////////////////////////////////////////////
  // Set the tests!
  ///////////////////////////////////////////////////////////


  /** TestChebDenseKernel */
  void TestChebDenseKernel(){
    typedef double FReal;
    const unsigned int ORDER = 6;
    typedef FP2PParticleContainerIndexed<FReal> ContainerClass;
    typedef FSimpleLeaf<FReal, ContainerClass> LeafClass;
    typedef FInterpMatrixKernelR<FReal> MatrixKernelClass;
    typedef FChebCell<FReal,ORDER> CellClass;
    typedef FOctree<FReal, CellClass,ContainerClass,LeafClass> OctreeClass;
    typedef FChebDenseKernel<FReal,CellClass,ContainerClass,MatrixKernelClass,ORDER> KernelClass;
    typedef FFmmAlgorithm<OctreeClass,CellClass,ContainerClass,KernelClass,LeafClass> FmmClass;
    // run test
    RunTest<FReal,CellClass,ContainerClass,KernelClass,MatrixKernelClass,LeafClass,OctreeClass,FmmClass>(
                           [&](int NbLevels, FReal boxWidth, FPoint<FReal> centerOfBox, const MatrixKernelClass *const MatrixKernel){
                             return std::unique_ptr<KernelClass>(new KernelClass(NbLevels, boxWidth, centerOfBox, MatrixKernel));
                           });
  }

  /** TestChebKernel */
  void TestChebKernel(){
    typedef double FReal;
    const unsigned int ORDER = 6;
    typedef FP2PParticleContainerIndexed<FReal> ContainerClass;
    typedef FSimpleLeaf<FReal, ContainerClass> LeafClass;
    typedef FInterpMatrixKernelR<FReal> MatrixKernelClass;
    typedef FChebCell<FReal,ORDER> CellClass;
    typedef FOctree<FReal, CellClass,ContainerClass,LeafClass> OctreeClass;
    typedef FChebKernel<FReal,CellClass,ContainerClass,MatrixKernelClass,ORDER> KernelClass;
    typedef FFmmAlgorithm<OctreeClass,CellClass,ContainerClass,KernelClass,LeafClass> FmmClass;
    // run test
    RunTest<FReal,CellClass,ContainerClass,KernelClass,MatrixKernelClass,LeafClass,OctreeClass,FmmClass>(
													 [&](int NbLevels, FReal boxWidth, FPoint<FReal> centerOfBox, const MatrixKernelClass *const MatrixKernel){
													   return std::unique_ptr<KernelClass>(new KernelClass(NbLevels, boxWidth, centerOfBox, MatrixKernel));
													 });
  }

  /** TestChebSymKernel */
  void TestChebSymKernel(){
    typedef double FReal;
    const unsigned int ORDER = 6;
    typedef FP2PParticleContainerIndexed<FReal> ContainerClass;
    typedef FSimpleLeaf<FReal, ContainerClass> LeafClass;
    typedef FInterpMatrixKernelR<FReal> MatrixKernelClass;
    typedef FChebCell<FReal,ORDER> CellClass;
    typedef FOctree<FReal, CellClass,ContainerClass,LeafClass> OctreeClass;
    typedef FChebSymKernel<FReal,CellClass,ContainerClass,MatrixKernelClass,ORDER> KernelClass;
    typedef FFmmAlgorithm<OctreeClass,CellClass,ContainerClass,KernelClass,LeafClass> FmmClass;
    // run test
    RunTest<FReal,CellClass,ContainerClass,KernelClass,MatrixKernelClass,LeafClass,OctreeClass,FmmClass>(
													 [&](int NbLevels, FReal boxWidth, FPoint<FReal> centerOfBox, const MatrixKernelClass *const MatrixKernel){
													   return std::unique_ptr<KernelClass>(new KernelClass(NbLevels, boxWidth, centerOfBox, MatrixKernel));
													 });
  }



  ///////////////////////////////////////////////////////////
  // Set the tests!
  ///////////////////////////////////////////////////////////

  /** set test */
  void SetTests(){
    AddTest(&TestChebyshevDirect::TestChebDenseKernel,"Test Chebyshev Kernel without compression.");
    AddTest(&TestChebyshevDirect::TestChebKernel,"Test Chebyshev Kernel with 1 large compression.");
    AddTest(&TestChebyshevDirect::TestChebSymKernel,"Test Chebyshev Kernel with 16 small SVDs and symmetries.");
  }
};


// You must do this
TestClass(TestChebyshevDirect)




