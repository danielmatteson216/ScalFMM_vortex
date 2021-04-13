/* ----------------------------------------------------------------------
   Compile
   g++ -v -I /usr/include/openblas/ -lblas -llapack testing.cpp -o testing
   ---------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   ProTip (to self) ~ it helps Reading The Fine Manual
   ---------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   define lapack_int   int
   define lapack_complex_double   double _Complex
   define lapack_complex_double_real(z)   (creal(z))
   define lapack_complex_double_imag(z)   (cimag(z))

   void LAPACK_zgeqrf(
      lapack_int            *m,
      lapack_int            *n,
      lapack_complex_double *A,
      lapack_int            *lda,
      lapack_complex_double *tau,
      lapack_complex_double *work,
      lapack_int            *lwork,
      lapack_int            *info
   );

   A/R : On exit, the elements on and above the diagonal of the array contain
         the min(M,N)-by-N upper trapezoidal matrix R (R is upper triangular
         if m >= n);
         the elements below the diagonal, with the array TAU, represent the
         orthogonal matrix Q as a product of min(m,n) elementary reflectors.

   Q :   The matrix Q is represented as a product of elementary reflectors
            Q = H(1) H(2) . . . H(k), where k = min(m,n).
         Each H(i) has the form
            H(i) = I - tau * v * v'
         where tau is a real scalar, and v is a real vector with v(1:i-1) = 0
         and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i), and tau in
         TAU(i).

   ---------------------------------------------------------------------- */

#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <complex>
//#include <lapack.h>
#include <lapacke.h>


//*****************************************************************************
//****************************  WRITE MATRIX  *********************************
//*****************************************************************************
void write_matrix( const std::complex<double> * const A,
		               const int &m, const int &n,
		               char const * const str ) {
  std::cout << "--- ("
	          << str
	          << ") ------------------------------------------------------------"
	          << std::endl;

    for( unsigned i=0; i<m; i++ ) {
        for( unsigned j=0; j<n; j ++ ) {
          std::cout << std::setw(18)
  		              << std::setprecision(4)
  		              << std::fixed
  		              << A[i+j*m] << " ";
        }
        std::cout << std::endl;
    }

std::cout << "-----------------------------------------------------------------------------"
	        << std::endl << std::endl;
}
//*****************************************************************************





//*****************************************************************************
//****************************      MAIN      *********************************
//*****************************************************************************
int main() {

  int     m=5;    // [in]  The number of rows of the matrix A.  M >= 0
  int     n=5;    // [in[  The number of columns of the matrix A.  N >= 0.
  int     lda;    // [in]  The leading dimension of the array A.  LDA >= max(1,M).
  int     lwork;  // [in]  The dimension of the array WORK.
  int     info;   // [out]

  std::complex<double> *A;      // [in,out] COMPLEX*16 array, dimension (LDA,N)
  std::complex<double> *tau;    // [out]    COMPLEX*16 array, dimension (min(M,N))
  std::complex<double> *work;   // [out]    COMPLEX*16 array, dimension (MAX(1,LWORK)

  std::complex<double> *R;      // [other]  COMPLEX*16 array, dimension (LDA,N)

  // The leading dimension of the array A.  LDA >= max(1,M).
  lda = (m>1) ? m : 1;


  std::cout << "------------------------------------------------------------"
	    << std::endl
	    << "Data type size information: " << std::endl
	    << "sizeof(std::complex<double>) = " << sizeof(std::complex<double>)
	    << std::endl
	    << "sizeof(double) = " << sizeof(double)
	    << std::endl
	    << "------------------------------------------------------------"
	    << std::endl << std::endl;

  A    =  new std::complex<double>[m*n];
  R    =  new std::complex<double>[m*n];
  tau  =  new std::complex<double>[m];
  work =  new std::complex<double>[1];


// ------------------------------------------------------ BUILD A -------------
  // A is a Hilbert Matrix
  for( int i=0; i<m; i++ ) {
    for( int j=0; j<n; j++ ) {
      A[i+j*m] = std::complex<double>(1.0/(1.0+i+j), 0.0);
    }
  }
// ------------------------------------------------------ BUILD A -------------


  write_matrix( A, m, n, "A" );        //------------------->  PRINT A MATRIX

  // Workspace Query.
  lwork = -1;

info =   LAPACKE_zgeqrf_work( 0,m,
		             n,
		             reinterpret_cast <__complex__ double*>(A),
		             lda,
		             reinterpret_cast <__complex__ double*>(tau),
		             reinterpret_cast <__complex__ double*>(work),
		             lwork);

  std::cout << "info:    " << info    << std::endl
	    << "lwork:   " << lwork   << std::endl
	    << "work[0]: " << work[0] << std::endl
	    << "work[0]: " << real(work[0])
	    << " this is the optimal LWORK"
	    << std::endl
	    << "work[0]: " << imag(work[0])
	    << std::endl << std::endl;

  // Resize lwork to optimal
  lwork = real(work[0]);
  delete[] work;
  work =  new std::complex<double>[lwork];

  // A should be untouched
  write_matrix( A, m, n, "A" );        //------------------->  PRINT A MATRIX

  // Compute Call
info =   LAPACKE_zgeqrf_work( 0,m,
		             n,
		             reinterpret_cast <__complex__ double*>(A),
		             lda,
		             reinterpret_cast <__complex__ double*>(tau),
		             reinterpret_cast <__complex__ double*>(work),
		             lwork);


  if( info ) {
    std::cout << "info:    " << info    << std::endl;
    return( info );
  }


  write_matrix( A, m, n, "A, after zgeqrf" );


// ------------------------------------------------------ BUILD R -------------
  // R is the upper triangular part of A
  for( unsigned i=0; i<m; i++ ) {
    for( unsigned j=i; j<m; j++ ) {
      R[i+j*m] = A[i+j*m];
    }
  }
// ------------------------------------------------------ BUILD R -------------


  write_matrix( R, m, n, "R" );        //------------------->  PRINT R MATRIX
  write_matrix( tau, m, 1, "tau" );    //------------------->  PRINT TAU

  // We need some outer-product action to actually build the Q-matrix.
  // (zgemm) (zgerc)

  delete[] A, tau, work;
  return(0);
}
//*****************************************************************************
