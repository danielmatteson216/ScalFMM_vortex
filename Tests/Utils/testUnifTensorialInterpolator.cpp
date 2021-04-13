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
// ================
// Keep in private GIT
// @SCALFMM_PRIVATE

#include <iostream>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "Utils/FTic.hpp"
#include "Utils/FMath.hpp"
#include "Utils/FBlas.hpp"

#include "Containers/FVector.hpp"

#include "Utils/FAssert.hpp"
#include "Utils/FPoint.hpp"

#include "Utils/FParameterNames.hpp"

#include "Kernels/Uniform/FUnifInterpolator.hpp"
#include "Kernels/Interpolation/FInterpMatrixKernel_TensorialInteractions.hpp"
#include "Kernels/Uniform/FUnifTensor.hpp"

// Check DFT
#include "Utils/stdComplex.hpp"
#include "Utils/FDft.hpp"


#include "Kernels/P2P/FP2PParticleContainer.hpp"
#include "Components/FSimpleLeaf.hpp"




/**
 * In this file we compute the interactions (direct and Unif FM-approximate) for
 * a tensorial interaction kernel (R_ij) as well as the forces (comparison
 * with direct computation using R_ijk kernel).
 */

int main(int argc, char ** argv){
    FHelpDescribeAndExit(argc, argv, "Test the uniform tensorial interpolator (only the code is interesting)");

    typedef double FReal;

    typedef FInterpMatrixKernel_R_IJ<FReal> MatrixKernelClass;
    const double a = 0.0; // core width (Beware! if diff from 0. then Kernel should be NON HOMOGENEOUS !!!)

    const unsigned int ncmp = MatrixKernelClass::NCMP;
    const unsigned int nrhs = MatrixKernelClass::NRHS;
    const unsigned int nlhs = MatrixKernelClass::NLHS;
    const unsigned int npot = MatrixKernelClass::NPOT;

    typedef FP2PParticleContainer<FReal,nrhs,nlhs> ContainerClass;
    typedef FSimpleLeaf<FReal, ContainerClass> LeafClass;

    ///////////////////////What we do/////////////////////////////
    std::cout << "\nTask: Compute interactions between source particles in leaf Y and target\n";
    std::cout << " particles in leaf X. Compare the fast summation K ~ Lx K Ly' with the\n";
    std::cout << " direct computation.\n" << std::endl;
    //////////////////////////////////////////////////////////////

    const FReal FRandMax = FReal(RAND_MAX);
    FTic time;


    // Leaf size
    FReal width = FReal(3.723);

    ////////////////////////////////////////////////////////////////////
    LeafClass X;
    FPoint<FReal> cx(0., 0., 0.);
    const unsigned long M = 5000;
    std::cout << "Fill the leaf X of width " << width
              << " centered at cx=" << cx << " with M=" << M << " target particles" << std::endl;
    {

        for(unsigned long i=0; i<M; ++i){
            FReal x = (FReal(rand())/FRandMax - FReal(.5)) * width + cx.getX();
            FReal y = (FReal(rand())/FRandMax - FReal(.5)) * width + cx.getY();
            FReal z = (FReal(rand())/FRandMax - FReal(.5)) * width + cx.getZ();
            // PB: need to know the actual value of NRHS (=3 here)
            X.push(FPoint<FReal>(x, y, z), FReal(rand())/FRandMax, FReal(rand())/FRandMax, FReal(rand())/FRandMax);
        }
    }


    ////////////////////////////////////////////////////////////////////
    LeafClass Y;
    FPoint<FReal> cy(FReal(2.)*width, 0., 0.);
    //FPoint<FReal> cy(FReal(2.)*width, FReal(2.)*width, 0.);

    const unsigned long N = 5000;
    std::cout << "Fill the leaf Y of width " << width
              << " centered at cy=" << cy	<< " with N=" << N << " source particles" << std::endl;
    {
        for(unsigned long i=0; i<N; ++i){
            FReal x = (FReal(rand())/FRandMax - FReal(.5)) * width + cy.getX();
            FReal y = (FReal(rand())/FRandMax - FReal(.5)) * width + cy.getY();
            FReal z = (FReal(rand())/FRandMax - FReal(.5)) * width + cy.getZ();
            // PB: need to know the actual value of NRHS (=3 here)
            Y.push(FPoint<FReal>(x, y, z), FReal(rand())/FRandMax, FReal(rand())/FRandMax, FReal(rand())/FRandMax);
        }
    }



    ////////////////////////////////////////////////////////////////////
    // approximative computation
    const unsigned int ORDER = 6;
    const unsigned int nnodes = TensorTraits<ORDER>::nnodes;
    typedef FUnifInterpolator<FReal,ORDER,MatrixKernelClass> InterpolatorClass;
    InterpolatorClass S;
    MatrixKernelClass MatrixKernel;

    std::cout << "\nCompute interactions approximatively, interpolation order = " << ORDER << " ..." << std::endl;

    std::cout << "\nP2M ... " << std::flush;
    time.tic();
    // Anterpolate: W_n = \sum_j^N S(y_j,\bar y_n) * w_j
    FReal W[nrhs*nnodes]; // multipole expansion
    // tensorial case interpolate same Y for each component
    S.applyP2M(cy, width, W, Y.getSrc()); // the multipole expansions are set to 0 in S.applyP2M
    std::cout << "took " << time.tacAndElapsed() << "s" << std::endl;

    std::cout << "M2L ... " << std::flush;
    time.tic();
    // Multipole to local: F_m = \sum_n^L K(\bar x_m, \bar y_n) * W_n
    FPoint<FReal> rootsX[nnodes], rootsY[nnodes];
    FUnifTensor<FReal,ORDER>::setRoots(cx, width, rootsX);
    FUnifTensor<FReal,ORDER>::setRoots(cy, width, rootsY);

    FReal F[nlhs*nnodes]; // local expansion
    for (unsigned int i=0; i<nnodes*nlhs; ++i) F[i] = FReal(0.);

    for (unsigned int i=0; i<nnodes; ++i) {
        for (unsigned int j=0; j<nnodes; ++j){

            for (unsigned int idxLhs=0; idxLhs<nlhs; ++idxLhs){
                unsigned int idxRhs = idxLhs % npot;
                unsigned int d = MatrixKernel.getPosition(idxLhs);

                F[i+idxLhs*nnodes] += MatrixKernelClass(a,d).evaluate(rootsX[i], rootsY[j]) * W[j+idxRhs*nnodes];

            }
        }
    }
    std::cout << "took " << time.tacAndElapsed() << "s" << std::endl;

    //  std::cout<< "F via direct applyM2L: "<<std::endl;
    //  for (unsigned int d=0; d<nlhs; ++d){
    //    for (unsigned int i=0; i<nnodes; ++i)
    //      std::cout<< F[i+d*nnodes] << ", ";
    //    std::cout<<std::endl;
    //  }
    //  std::cout<<std::endl;

    ////////////////////////////////////////////////////////////////////////////
    // Store M2L in K and apply K
    FReal K[ncmp*nnodes*nnodes]; // local expansion
    for (unsigned int i=0; i<nnodes; ++i) {
        for (unsigned int j=0; j<nnodes; ++j){

            for (unsigned int d=0; d<ncmp; ++d){
                K[d*nnodes*nnodes+i*nnodes+j] = MatrixKernelClass(a,d).evaluate(rootsX[i], rootsY[j]);
            }

        }
    }
    std::cout<< "Apply M2L in usual sense: ";
    time.tic();
    for (unsigned int i=0; i<nnodes*nlhs; ++i) F[i] = FReal(0.);

    for (unsigned int i=0; i<nnodes; ++i)
        for (unsigned int j=0; j<nnodes; ++j){

            for (unsigned int idxLhs=0; idxLhs<nlhs; ++idxLhs){
                unsigned int idxRhs = idxLhs % npot;
                unsigned int d = MatrixKernel.getPosition(idxLhs);

                F[i+idxLhs*nnodes] += K[d*nnodes*nnodes+i*nnodes+j] * W[j+idxRhs*nnodes];

            }
        }

    time.tac();
    std::cout << "took " << time.elapsed() << "sec." << std::endl;

    //  std::cout<< "F via store and applyM2L: "<<std::endl;
    //  for (unsigned int d=0; d<nlhs; ++d){
    //    for (unsigned int i=0; i<nnodes; ++i)
    //      std::cout<< F[i+d*nnodes] << ", ";
    //    std::cout<<std::endl;
    //  }
    //  std::cout<<std::endl;

    /////////////////////////////////////////////////////////////////////////////////////
    // PB: Verify storage improvement works (indexing etc...)
    // 1) store circulant matrix
    const unsigned int rc = (2*ORDER-1)*(2*ORDER-1)*(2*ORDER-1);
    FReal C[ncmp*rc];

    typedef FUnifTensor<FReal,ORDER> TensorType;
    unsigned int node_diff[nnodes*nnodes];
    TensorType::setNodeIdsDiff(node_diff);
    unsigned int node_ids_pairs[rc][2];
    TensorType::setNodeIdsPairs(node_ids_pairs);

    unsigned int ido=0;

    for(unsigned int l=0; l<2*ORDER-1; ++l)
        for(unsigned int m=0; m<2*ORDER-1; ++m)
            for(unsigned int n=0; n<2*ORDER-1; ++n){

                for (unsigned int d=0; d<ncmp; ++d){

                    C[d*rc + ido]
                            = MatrixKernelClass(a,d).evaluate(rootsX[node_ids_pairs[ido][0]],
                            rootsY[node_ids_pairs[ido][1]]);
                }

                ido++;
            }

    //  // Display C (gathers every values of K that need to be stored,
    //  // corresponds to the first line of the padded matrix (reverse order?))
    //  std::cout<<"C="<<std::endl;
    //    for (unsigned int d=0; d<ncmp; ++d){
    //      for (unsigned int n=0; n<rc; ++n)
    //        std::cout<< C[n + d*rc] << ", ";
    //      std::cout<<std::endl;
    //    }
    //  std::cout<<std::endl;

    //////////////////////////////////////////////////////////////////////////////////////////////////////
    // K is a block Toeplitz matrix
    // i.e. a blockwise Toeplitz matrix where the block also have the Toeplitz structure.
    // e.g. for ORDER=3: K=[K_{1,1} K_{1,2} K_{1,3},  where K_{i,j}=[k11 k12 k13,
    //                      K_{2,1} K_{1,1} K_{1,2},                 k21 k11 k12,
    //                      K_{3,1} K_{2,1} K_{1,1}];                k31 k21 k11];
    // K is of size order^3 x order^3
    // (i.e. order^2 x order^2 Toeplitz blocks of size order x order),
    // K is very close to be Toeplitz itself and even circulant.
    // In order to actually embed K into a circulant matrix C one just
    // needs to insert (ORDER-1) extra lines/columns (to each block).

    //  std::cout<< "K=" <<std::endl;
    //  for (unsigned int i=0; i<nnodes; ++i){
    //    for (unsigned int j=0; j<nnodes; ++j){
    //      std::cout<< K[i*nnodes+j]<<", ";
    //    }
    //    std::cout<<std::endl;
    //  }
    //  std::cout<<std::endl;

    //  // Check matrix node_diff
    //  std::cout<< "node_diff=" <<std::endl;
    //  for (unsigned int i=0; i<nnodes; ++i){
    //    for (unsigned int j=0; j<nnodes; ++j){
    //      std::cout<< node_diff[i*nnodes+j] <<", ";
    //    }
    //    std::cout<<std::endl;
    //  }
    //  std::cout<<std::endl;

    //  // Expected ido for the (2*ORDER-1)^3x(2*ORDER-1)^3 circulant matrix
    //  for (unsigned int i=0; i<rc; ++i){
    //    for (unsigned int j=0; j<rc; ++j){
    //      if(i>j) std::cout<< i-j-1 << ", ";
    //      else std::cout<< rc+i-j-1 << ", ";
    //    } std::cout<<std::endl;
    //  } std::cout<<std::endl;

    //  // kernel evaluated at previous ido returns a circulant matrix
    //  for (unsigned int i=0; i<rc/2; ++i){
    //    for (unsigned int j=0; j<rc/2; ++j){
    //      if(i>j) std::cout<< C[i-j-1] << ", ";
    //      else std::cout<< C[rc+i-j-1] << ", ";
    //    } std::cout<<std::endl;
    //  } std::cout<<std::endl;

    // In 1D the Zero Padding consists in
    // inserting ORDER-1 zeros in the multipole exp
    // in order to apply the (ORDER+ORDER-1)x(ORDER+ORDER-1)
    // circulant matrix to it.
    // Let us extend it to the 3D case:
    FReal MultExp[nrhs*nnodes]; FReal PaddedMultExp[nrhs*rc];
    for (unsigned int i=0; i<nrhs*nnodes; ++i) MultExp[i]=W[i];
    FReal LocalExp[nlhs*nnodes]; FReal PaddedLocalExp[nlhs*rc];
    FBlas::setzero(nlhs*nnodes,LocalExp);
    FBlas::setzero(nlhs*rc,PaddedLocalExp);

    //  std::cout<< "Expected LocalExp: "<<std::endl;
    //  for (unsigned int d=0; d<nlhs; ++d){
    //    for (unsigned int i=0; i<nnodes; ++i)
    //      std::cout<< F[i+d*nnodes] << ", ";
    //    std::cout<<std::endl;
    //  }

    /////////////////////////////////////////////////////////////////////////////////////
    // Application of circulant Toeplitz system in PHYSICAL SPACE
    std::cout<< "Apply circulant M2L in Physical space: ";
    time.tic();
    for (unsigned int i=0; i<nnodes; ++i){

        // Pad Multipole Expansion with respect to current row
        FBlas::setzero(nrhs*rc,PaddedMultExp);
        for (unsigned int j=0; j<nnodes; ++j)
            for (unsigned int d=0; d<nrhs; ++d)
                PaddedMultExp[node_diff[i*nnodes+j] + d*rc]=MultExp[j + d*nnodes];

        //    std::cout<< "Padded MultExp for row i=" << i << ": "<<std::endl;
        //    for (unsigned int p=0; p<rc; ++p)
        //      std::cout<< PaddedMultExp[p] << ", ";
        //    std::cout<<std::endl;

        // Application of M2L in PHYSICAL SPACE
        for (unsigned int pj=0; pj<rc; ++pj)
            for (unsigned int idxLhs=0; idxLhs<nlhs; ++idxLhs){
                unsigned int idxRhs = idxLhs % npot;
                unsigned int d = MatrixKernel.getPosition(idxLhs);

                LocalExp[i + idxLhs*nnodes]+=C[pj + d*rc]*PaddedMultExp[pj + idxRhs*rc];

            }

    }// end i
    time.tac();
    std::cout << "took " << time.elapsed() << "sec." << std::endl;

    //  std::cout<< "LocalExp via product in PHYSICAL SPACE: "<<std::endl;
    //  for (unsigned int d=0; d<nlhs; ++d){
    //    for (unsigned int i=0; i<nnodes; ++i)
    //      std::cout<< LocalExp[i+d*nnodes] << ", ";
    //    std::cout<<std::endl;
    //  }
    //  std::cout<<std::endl;

    /////////////////////////////////////////////////////////////////////////////////////
    // Efficient application of the Toeplitz system in FOURIER SPACE

    // Init DFT
    const int dimfft = 1;
    //FDft Dft(rc); // direct version
    FFftw<FReal,stdComplex<FReal>,dimfft> Dft(rc); // fast version

    // Get first COLUMN of K and Store in T
    FReal T[ncmp*rc];
    // use permutations
    unsigned int perm[rc];
    for(unsigned int p=0; p<rc; ++p){
        if(p>0) perm[p]=p-1;
        else perm[p]=rc+p-1;
        //    std::cout << "perm["<< p << "]="<< perm[p] << std::endl;
    }

    for (unsigned int i=0; i<rc; ++i){
        // keep this lines commented to see what permutation accounts for:
        //  for (unsigned int j=0; j<rc; ++j){
        //      if(i>0) T[i]=C[i-0-1];
        //      else T[i]=C[rc+i-0-1];
        for (unsigned int d=0; d<ncmp; ++d)
            T[i + d*rc]=C[perm[i] + d*rc];
        //  }
    }

    //  std::cout<< "First column of C[rc x rc]: "<<std::endl;
    //  for (unsigned int p=0; p<rc; ++p)
    //    std::cout<< T[p] << ", ";
    //  std::cout<<std::endl;

    // Apply DFT to T
    stdComplex<FReal> FT[ncmp*rc];
    //  for (unsigned int n=0; n<rc; ++n) FT[n]=stdComplex<FReal>(0.0,0.0);
    FBlas::c_setzero(ncmp*rc,reinterpret_cast<FReal*>(FT));

    // if first COLUMN (T) of C is used
    for (unsigned int d=0; d<ncmp; ++d)
        Dft.applyDFT(T+d*rc,FT+d*rc);
    //  // if first ROW of C is used
    //  Dft.applyDFT(C,FT);

    stdComplex<FReal> FPMultExp[nrhs*rc];
    stdComplex<FReal> FPLocalExp[nlhs*rc];
    FReal PLocalExp[nlhs*rc];

    //for (unsigned int n=0; n<rc; ++n) FPLocalExp[n]=stdComplex<FReal>(0.0,0.0);
    FBlas::c_setzero(nlhs*rc,reinterpret_cast<FReal*>(FPLocalExp));

    FBlas::setzero(nlhs*rc,PLocalExp);

    // Pad physical MultExp
    FBlas::setzero(nrhs*rc,PaddedMultExp); //part of padding
    for (unsigned int idRhs=0; idRhs<nrhs; ++idRhs)
        for (unsigned int j=0; j<nnodes; ++j){
            // if first COLUMN (T) of C is used
            PaddedMultExp[node_diff[j*nnodes]+idRhs*rc]=MultExp[j+idRhs*nnodes];
            //    // if first ROW of C is used
            //    PaddedMultExp[node_diff[j]]=MultExp[j];
        }

    //    std::cout<< "Padded MultExp: "<<std::endl;
    //    for (unsigned int p=0; p<rc; ++p)
    //      std::cout<< PaddedMultExp[p] << ", ";
    //    std::cout<<std::endl;


    // Set transformed MultExp to 0
    //  for (unsigned int n=0; n<rc; ++n) FPMultExp[n]=stdComplex<FReal>(0.0,0.0);
    FBlas::c_setzero(nrhs*rc,reinterpret_cast<FReal*>(FPMultExp));

    // Transform PaddedMultExp
    for (unsigned int idxRhs=0; idxRhs<nrhs; ++idxRhs) // apply nrhs 1 dimensionnal FFT
        Dft.applyDFT(PaddedMultExp+idxRhs*rc,FPMultExp+idxRhs*rc);

    std::cout<< "Apply M2L in Fourier space: ";
    time.tic();

    // Application of M2L in FOURIER SPACE
    stdComplex<FReal> tmpFX;
    for (unsigned int idxLhs=0; idxLhs<nlhs; ++idxLhs){
        unsigned int idxRhs = idxLhs % npot;
        unsigned int d = MatrixKernel.getPosition(idxLhs);

        for (unsigned int pj=0; pj<rc; ++pj){
            tmpFX=FT[pj + d*rc];
            tmpFX*=FPMultExp[pj+idxRhs*rc];
            FPLocalExp[pj+idxLhs*rc]+=tmpFX; // add new contribution +RijYj
        }

    }
    time.tac();
    std::cout << "took " << time.elapsed() << "sec." << std::endl;

    //    std::cout<< "Transfo Padded LocalExp: "<<std::endl;
    //    for (unsigned int p=0; p<rc; ++p)
    //      std::cout<< FPLocalExp[p] << ", ";
    //    std::cout<<std::endl;

    for (unsigned int idxLhs=0; idxLhs<nlhs; ++idxLhs) // apply nrhs 1 dimensionnal FFT
        Dft.applyIDFTNorm(FPLocalExp+idxLhs*rc,PLocalExp+idxLhs*rc);

    //  std::cout<< "Padded LocalExp: "<<std::endl;
    //  for (unsigned int p=0; p<rc; ++p)
    //    std::cout<< PLocalExp[p] << ", ";
    //  std::cout<<std::endl;

    // Unpad
    for (unsigned int idLhs=0; idLhs<nlhs; ++idLhs)
        for (unsigned int j=0; j<nnodes; ++j){
            // if first COLUMN (T) of C is used
            LocalExp[j+idLhs*nnodes]=PLocalExp[node_diff[nnodes-j-1]+idLhs*rc];
            //    // if first ROW of C is used
            //    LocalExp[j]=PLocalExp[node_diff[j*nnodes]];
        }

    //  std::cout<< "Mask to be applied to Padded LocalExp: "<<std::endl;
    //  for (unsigned int j=0; j<nnodes; ++j)
    //    std::cout<< node_diff[nnodes-j-1] << ", ";
    //  std::cout<<std::endl;

    //  std::cout<< "LocalExp via product in FOURIER SPACE: "<<std::endl;
    //  for (unsigned int d=0; d<nlhs; ++d){
    //    for (unsigned int p=0; p<nnodes; ++p)
    //      std::cout<< LocalExp[p + d*nnodes] << ", ";
    //    std::cout<<std::endl;
    //  }
    //  std::cout<<std::endl;


    /////////////////////////////////////////////////////////////////////////////////////

    std::cout << "L2P (potential) ... " << std::flush;
    time.tic();
    // Interpolate p_i = \sum_m^L S(x_i,\bar x_m) * F_m
    S.applyL2P(cx, width, F, X.getTargets());
    std::cout << "took " << time.tacAndElapsed() << "s" << std::endl;

    std::cout << "L2P (forces) ... " << std::flush;
    time.tic();
    // Interpolate f_i = \sum_m^L P(x_i,\bar x_m) * F_m
    S.applyL2PGradient(cx, width, F, X.getTargets());
    std::cout << "took " << time.tacAndElapsed() << "s" << std::endl;

    ////////////////////////////////////////////////////////////////////
    // direct computation
    std::cout << "Compute interactions directly ..." << std::endl;
    time.tic();

    FReal** approx_f = new FReal* [npot];
    FReal**        f = new FReal* [npot];
    for (unsigned int i=0; i<npot; ++i){
        approx_f[i] = new FReal [M * 3];
        f[i] = new FReal [M * 3];
        FBlas::setzero(M*3, f[i]);
    }

    FReal** approx_p = new FReal* [npot];
    FReal**        p = new FReal* [npot];
    for (unsigned int i=0; i<npot; ++i){
        approx_p[i] = new FReal [M];
        p[i] = new FReal [M];
        FBlas::setzero(M, p[i]);
    }

    { // start direct computation
        unsigned int counter = 0;

        for(FSize idxPartX = 0 ; idxPartX < X.getSrc()->getNbParticles() ; ++idxPartX){
            const FPoint<FReal> x = FPoint<FReal>(X.getSrc()->getPositions()[0][idxPartX],
                    X.getSrc()->getPositions()[1][idxPartX],
                    X.getSrc()->getPositions()[2][idxPartX]);
            const FReal wx[nrhs] = {X.getSrc()->getPhysicalValues(0)[idxPartX],
                                    X.getSrc()->getPhysicalValues(1)[idxPartX],
                                    X.getSrc()->getPhysicalValues(2)[idxPartX]};

            for(FSize idxPartY = 0 ; idxPartY < Y.getSrc()->getNbParticles() ; ++idxPartY){
                const FPoint<FReal> y = FPoint<FReal>(Y.getSrc()->getPositions()[0][idxPartY],
                        Y.getSrc()->getPositions()[1][idxPartY],
                        Y.getSrc()->getPositions()[2][idxPartY]);
                const FReal wy[nrhs] = {Y.getSrc()->getPhysicalValues(0)[idxPartY],
                                        Y.getSrc()->getPhysicalValues(1)[idxPartY],
                                        Y.getSrc()->getPhysicalValues(2)[idxPartY]};

                //        // 1/R
                //        const FReal one_over_r = MatrixKernel.evaluate(x, y);
                //
                //        // potential
                //        p[counter] += one_over_r * wy;
                //        // force
                //        FPoint<FReal> force(y - x);
                //        force *= one_over_r*one_over_r*one_over_r;
                //        f[counter*3 + 0] += force.getX() * wx * wy;
                //        f[counter*3 + 1] += force.getY() * wx * wy;
                //        f[counter*3 + 2] += force.getZ() * wx * wy;

                FReal rij[ncmp];
                FReal rijk[ncmp][3];
                MatrixKernel.evaluateBlockAndDerivative(x,y,rij,rijk);

                // R,ij and (R,ij),k
                for (unsigned int i=0; i<npot; ++i) // sum all compo
                    for (unsigned int j=0; j<nrhs; ++j){
                        unsigned int d = MatrixKernel.getPosition(i*nrhs+j);
                        // potential
                        p[i][counter] += rij[d] * wy[j];
                        // force
                        FReal force_k;
                        for (unsigned int k=0; k<3; ++k){
                            // Convention in matrix kernel: R_ij(x-y), while R_ijk(y-x)
                            force_k = FReal(1.) * rijk[d][k];
                            f[i][counter*3 + k] += force_k * wx[j] * wy[j];
                        }
                    }

            }
            counter++;
        }
    } // end direct computation


    time.tac();
    std::cout << "Done in " << time.elapsed() << "sec." << std::endl;


    ////////////////////////////////////////////////////////////////////
    unsigned int counter = 0;
    for(FSize idxPartX = 0 ; idxPartX < X.getSrc()->getNbParticles() ; ++idxPartX){
        for (unsigned int i=0; i<npot; ++i){
            approx_p[i][counter] = X.getSrc()->getPotentials(i)[idxPartX];
            const FPoint<FReal> force = FPoint<FReal>(X.getSrc()->getForcesX(i)[idxPartX],
                                        X.getSrc()->getForcesY(i)[idxPartX],
                                        X.getSrc()->getForcesZ(i)[idxPartX]);
            approx_f[i][counter*3 + 0] = force.getX();
            approx_f[i][counter*3 + 1] = force.getY();
            approx_f[i][counter*3 + 2] = force.getZ();
        }
        counter++;
    }

    //  std::cout << "Check Potential, forceX, forceY, forceZ " << std::endl;
    //  for (unsigned int i=0; i<npot; ++i){
    //    std::cout<< "idxLhs="<< i << std::endl;
    //    for(FSize idxPart = 0 ; idxPart < 20 ; ++idxPart){
    //      std::cout << approx_p[i][idxPart]     << ", "<< p[i][idxPart] << "|| ";
    //      std::cout << approx_f[i][idxPart]     << ", "<< f[i][idxPart] << "|| ";
    //      std::cout << approx_f[i][idxPart+M]   << ", "<< f[i][idxPart+M] << "|| ";
    //      std::cout << approx_f[i][idxPart+2*M] << ", "<< f[i][idxPart+2*M] << "|| ";
    //      std::cout << std::endl;
    //    }
    //    std::cout << std::endl;
    //  }
    //  std::cout << std::endl;

    std::cout << "\nRelative Inf/L2 errors: " << std::endl;
    std::cout << "  Potential:" << std::endl;
    for(unsigned i = 0; i<npot;++i) {
        std::cout << "    " << i << ": "
                  << FMath::FAccurater<FReal>(p[i],approx_p[i],M).getRelativeInfNorm()<<", "
                  << FMath::FAccurater<FReal>(p[i],approx_p[i],M).getRelativeL2Norm()
                  << std::endl;
    }
    std::cout << std::endl;

    std::cout << "  Force:" << std::endl;
    for(unsigned i = 0; i<npot;++i) {
        std::cout << "    " << i << ": "
                  << FMath::FAccurater<FReal>(f[i],approx_f[i],3*M).getRelativeInfNorm()<<", "
                  << FMath::FAccurater<FReal>(f[i],approx_f[i],3*M).getRelativeL2Norm()
                  << std::endl;
    }
    std::cout << std::endl;

    // free memory
    delete [] approx_p;
    delete [] p;
    delete [] approx_f;
    delete [] f;


    return 0;
}
