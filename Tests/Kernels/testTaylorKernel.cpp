// See LICENCE file at project root

#include <string>

#include "Utils/FPoint.hpp"
#include "Utils/FLog.hpp"
#include "Utils/FMath.hpp"

#include "Kernels/Taylor/FTaylorCell.hpp"
#include "Kernels/Taylor/FTaylorKernel.hpp"

#include "Kernels/P2P/FP2PParticleContainer.hpp"

#include "Components/FSimpleLeaf.hpp"

#include "Containers/FVector.hpp"
#include "Containers/FOctree.hpp"

#include "Core/FFmmAlgorithm.hpp"
#include "Core/FFmmAlgorithmThread.hpp"
#include "Core/FFmmAlgorithmTask.hpp"
#include "Utils/FParameterNames.hpp"

int main(int argc,char** argv){
    FHelpDescribeAndExit(argc, argv,
                         "Compile a Taylor Kernel (but do nothing).");
    const int P = 10;
    const int order = 1;
    typedef double FReal;
    FPoint<FReal> centerBox = FPoint<FReal>(0,0,0);

    typedef FTaylorCell<FReal, P,order> CellClass;
    typedef FP2PParticleContainer<FReal> ContainerClass;
    typedef FTaylorKernel<FReal,CellClass,ContainerClass,P,order> KernelClass;
    //typedef FSimpleLeaf<FReal, ContainerClass > LeafClass;
    //typedef FOctree<FReal, CellClass, ContainerClass , LeafClass > OctreeClass;

    KernelClass kernel(9,1.0,centerBox);

    return 0;
}
