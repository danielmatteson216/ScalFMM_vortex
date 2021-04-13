// ===================================================================================
// Ce LOGICIEL "ScalFmm" est couvert par le copyright Inria 20xx-2012.
// Inria détient tous les droits de propriété sur le LOGICIEL, et souhaite que
// la communauté scientifique l'utilise afin de le tester et de l'évaluer.
// Inria donne gracieusement le droit d'utiliser ce LOGICIEL. Toute utilisation
// dans un but lucratif ou à des fins commerciales est interdite sauf autorisation
// expresse et préalable d'Inria.
// Toute utilisation hors des limites précisées ci-dessus et réalisée sans l'accord
// expresse préalable d'Inria constituerait donc le délit de contrefaçon.
// Le LOGICIEL étant un produit en cours de développement, Inria ne saurait assurer
// aucune responsabilité et notamment en aucune manière et en aucun cas, être tenu
// de répondre d'éventuels dommages directs ou indirects subits par l'utilisateur.
// Tout utilisateur du LOGICIEL s'engage à communiquer à Inria ses remarques
// relatives à l'usage du LOGICIEL
// ===================================================================================

// ==== CMAKE =====
// @FUSE_FFT
// @FUSE_BLAS
// ================
// Keep in private GIT
// @SCALFMM_PRIVATE



#include <iostream>

#include <stdlib.h>
#include <time.h>

#include "Utils/FTic.hpp"
#include "Utils/FMath.hpp"
#include "Utils/FBlas.hpp" 

// 
#include "Utils/stdComplex.hpp"
#include "Utils/FDft.hpp"

#include "Utils/FParameterNames.hpp"


int main(int argc, char ** argv){
    FHelpDescribeAndExit(argc, argv, "Test the FFT (only the code is interesting).");

    typedef double FReal;
    const FReal FRandMax = FReal(RAND_MAX);

    FTic time;

    const unsigned int ORDER = 3;
    const unsigned int N=(2*ORDER-1);


    //////////////////////////////////////////////////////////////////////////////
    // Application of Convolution Theorem
    // Let there be a Circulant Toeplitz matrix K of size N

    FReal K[N*N];
    FReal Y[N];
    FReal X[N];

    for(unsigned int i=0; i<N; ++i) X[i]=0.0;
    for(unsigned int i=0; i<N; ++i) Y[i]=FReal(rand())/FRandMax;

    std::cout<< "Y: "<<std::endl;
    for(unsigned int i=0; i<N; ++i)
        std::cout<< Y[i] << ", ";
    std::cout<<std::endl;

    // Assemble Matrix K
    std::cout<< "Assemble Full Matrix K: ";
    time.tic();
    for(unsigned int i=0; i<N; ++i){
        for(unsigned int j=0; j<N; ++j){
            if(i>j) K[i*N+j]=i-j-1;
            else   K[i*N+j]=N+i-j-1;
        }
    }
    std::cout << "took " << time.tacAndElapsed() << "sec." << std::endl;

    std::cout<< "Circulant Toeplitz matrix K: "<<std::endl;
    for(unsigned int i=0; i<N; ++i){
        for(unsigned int j=0; j<N; ++j){
            std::cout<< K[i*N+j] << ", ";
        } std::cout<<std::endl;
    } std::cout<<std::endl;

    // Direct application of Circulant Matrix
    std::cout<< "Direct application of K: ";
    time.tic();
    for(unsigned int i=0; i<N; ++i){
        for(unsigned int j=0; j<N; ++j){
            X[i]+=K[i*N+j]*Y[j];
        }
    }
    std::cout << "took " << time.tacAndElapsed() << "sec." << std::endl;

    std::cout<< "X=KY: "<<std::endl;
    for(unsigned int i=0; i<N; ++i)
        std::cout<< X[i] << ", ";
    std::cout<<std::endl;

    // now compute via DFT and use convolution theorem
    stdComplex<FReal> FK[N];
    stdComplex<FReal> FY[N];
    stdComplex<FReal> FX[N];
    FReal iFX[N];

    // Init DFTor
    std::cout<< "Set DFT: ";
    time.tic();
    const int dim = 1;
    //FDft<FReal> Dft(N);// direct version (Beware! Ordering of output differs from REAL valued-FFT)
    FFftw<FReal,stdComplex<FReal>,dim> Dft(N);// fast version
    std::cout << "took " << time.tacAndElapsed() << "sec." << std::endl;

    // Initialize manually
    //  for(unsigned int s=0; s<N; ++s){
    //    FX[s] = FY[s] = FK[s] =stdComplex<FReal>(0.0,0.0); // init
    //    iFX[s]=0.0;
    //  }
    // ... or using Blas routines
    FBlas::c_setzero(N,reinterpret_cast<FReal*>(FX));
    FBlas::c_setzero(N,reinterpret_cast<FReal*>(FY));
    FBlas::c_setzero(N,reinterpret_cast<FReal*>(FK));
    FBlas::setzero(N,iFX);


    // Transform Y in Fourier space
    std::cout<< "Transform Y->FY: ";
    time.tic();
    Dft.applyDFT(Y,FY);
    std::cout << "took " << time.tacAndElapsed() << "sec." << std::endl;

    std::cout<< "FY: "<<std::endl;
    for(unsigned int i=0; i<N; ++i){
        std::cout<< FY[i] << ", ";
    }std::cout<<std::endl;

    // Transform first column of K
    FReal tK[N];
    for(unsigned int i=0; i<N; ++i) tK[i]=K[i*N]; // first column
    //  for(unsigned int i=0; i<N; ++i) tK[i]=K[i]; // first row
    std::cout<< "Transform tK->FK: ";
    time.tic();
    Dft.applyDFT(tK,FK);
    std::cout << "took " << time.tacAndElapsed() << "sec." << std::endl;

    //  std::cout<< "Transformed Matrix FK: "<<std::endl;
    //  for(unsigned int i=0; i<N; ++i){
    //    std::cout<< FK[i] << ", ";
    //  }std::cout<<std::endl;

    // Compute X=KY in Fourier space
    std::cout<< "Apply K in Fourier space (entrywise prod.): ";
    time.tic();
    for(unsigned int s=0; s<N; ++s){// TODO loop only over non zero entries
        FX[s]=FK[s];
        FX[s]*=FY[s];
    }
    std::cout << "took " << time.tacAndElapsed() << "sec." << std::endl;

    std::cout<< "Transformed Local Exp FX=FK:FY: "<<std::endl;
    for(unsigned int i=0; i<N; ++i){
        std::cout<< FX[i] << ", ";
    }std::cout<<std::endl;

    // Transform FX back in Physical space
    std::cout<< "Transform FX->iF(FX): ";
    time.tic();
    Dft.applyIDFTNorm(FX,iFX);
    std::cout << "took " << time.tacAndElapsed() << "sec." << std::endl;

    std::cout<< "iF(FK:FY): "<<std::endl;
    for(unsigned int i=0; i<N; ++i){
        std::cout<< iFX[i] << ", ";
    }std::cout<<std::endl;

    std::cout << "Relative error   = " << FMath::FAccurater<FReal>( X, iFX, N) << std::endl;


    return 0;
}
