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
// @FUSE_BLAS
// ================

#include <iostream>
#include <stdlib.h>
#include <stdexcept>

#include "Utils/FGlobal.hpp"
#include "Utils/FBlas.hpp"
#include "Utils/FTic.hpp"
#include "Utils/FParameterNames.hpp"

/**
 * Test functionality of C - interfaced LAPACK functions : GEQP3
 */

int main(int argc, char ** argv)
{
    FHelpDescribeAndExit(argc, argv, "Test the lapack compilation and linking (only the code is interesting).");

    /*
   * List of tested functions:
   * QR decomposition with pivoting: FBlas::geqp3()
   */
    
    FTic time;

    typedef double FReal;
    const FReal FRandMax = FReal(RAND_MAX);

    const unsigned int m = 10, n = 10;
    FReal* A = new FReal [m * n]; // matrix: column major ordering

    // A= LL^T ////////////////////////////////////
    // define symmetric definite positive matrix A
    for (unsigned int i=0; i<m; ++i) 
        for (unsigned int j=0; j<n; ++j)
            A[i+j*m]=FReal(rand())/FRandMax;

    // copy A in C
    FReal* C = new FReal [m * n]; // matrix: column major ordering
    for (unsigned int ii=0; ii<m; ++ii)
        for (unsigned int jj=0; jj<n; ++jj)
            C[ii + jj*m]=A[ii + jj*m];

    std::cout<<"\nA=["<<std::endl;
    for (unsigned int i=0; i<m; ++i) {
        for (unsigned int j=0; j<n; ++j)
            std::cout << A[i+j*m] << " ";
        std::cout<< std::endl;
    }
    std::cout<<"]"<<std::endl;

    // Workspace
    // init SVD
    const unsigned int minMN = std::min(m,n);
    const unsigned int maxMN = std::max(m,n);
    //const FSize LWORK = 2*4*minMN; // for square matrices
    const unsigned int LWORK = 2*std::max(3*minMN+maxMN, 5*minMN);
    FReal *const WORK = new FReal[LWORK];
    FReal* TAU = new FReal[n];

    // Pivot
    unsigned* jpiv = new unsigned[n]; 

    // perform Cholesky decomposition
    std::cout<<"\nQR decomposition ";
    FTic timeQR;
    //int INF = FBlas::geqrf(m, n, A, TAU, LWORK, WORK);
    int INFO = FBlas::geqp3(m, n, A, jpiv, TAU, LWORK, WORK);
     //FBlas::geqp3(m, n, A, jpiv);

    if(INFO==0) {std::cout<<"succeeded!"<<std::endl;}
    else {std::cout<<"failed!"<<std::endl;}



    std::cout << "took " << timeQR.tacAndElapsed() << "sec." << std::endl;

    std::cout<<"\nA_out=["<<std::endl;
    for (unsigned int i=0; i<m; ++i) {
        for (unsigned int j=0; j<n; ++j)
            std::cout << A[i+j*m] << " ";
        std::cout<<std::endl;
    }
    std::cout<<"]"<<std::endl;

    std::cout<<"\njpiv=["<<std::endl;
    for (unsigned int i=0; i<n; ++i) {
            std::cout << jpiv[i] << " ";
    }
    std::cout<<"]"<<std::endl;


    // get Q
    FTic timeGETQ;
    timeGETQ.tic();
    INFO = FBlas::orgqr(m, n, A, TAU, int(LWORK), WORK);
    if(INFO!=0) {
      std::stringstream stream;
      stream << INFO;
      delete [] A;
      delete [] C;
      delete [] jpiv;
      delete [] TAU;
      delete [] WORK;
      throw std::runtime_error("get Q failed! INFO=" + stream.str());
    }
    double tGETQ = timeGETQ.tacAndElapsed();
    std::cout << "... took @tGETQ = "<< tGETQ <<"\n";


    std::cout<<"\nQ=["<<std::endl;
    for (unsigned int i=0; i<m; ++i) {
        for (unsigned int j=0; j<n; ++j)
            std::cout << A[i+j*m] << " ";
        std::cout<< std::endl;
    }
    std::cout<<"]"<<std::endl;

    delete [] TAU;
    delete [] WORK;

    // TODO Verify QR facto! How? Form R then apply Q? How?

    delete [] A;
    delete [] C;
    delete [] jpiv;

    std::cout << "took " << time.tacAndElapsed() << "sec." << std::endl;

    return 0;
}
