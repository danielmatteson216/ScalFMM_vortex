// See LICENCE file at project root
#ifndef FBLAS_HPP
#define FBLAS_HPP

#include "ScalFmmConfig.h"
#include "FGlobal.hpp"
#include "FFortranMangling.hpp"

#ifndef SCALFMM_USE_BLAS
#error The BLAS header is included while SCALFMM_USE_BLAS is turned OFF
#endif

// This file interfaces the blas functions
// to enable a generic use.
// If no blas has been enabled in the cmake,
// the function will be empty


// for real
namespace scalfmm {
  const double D_ZERO =  0.0;
  const double D_ONE  =  1.0;
  const double D_MONE = -1.0;
  const float  S_ZERO =  0.0;
  const float  S_ONE  =  1.0;
  const float  S_MONE = -1.0;
  // for complex
  const double Z_ZERO[2] =  {0.0,0.0};
  const double Z_ONE[2]  =  {1.0,0.0};
  const double Z_MONE[2] =  {-1.0,0.0};
  const float  C_ZERO[2] =  {0.0,0.0};
  const float  C_ONE[2]  =  {1.0,0.0};
  const float  C_MONE[2] =  {-1.0,0.0};

  //const double D_PREC = 1e-16;

  const unsigned N_ONE = 1;
  const int N_MONE = -1;
  const char JOB_STR[] = "NTOSVULCR";
}

extern "C"
{
  // double //////////////////////////////////////////////////////////
  // blas 1
  double Fddot(const unsigned*, const double*, const unsigned*, const double*, const unsigned*);
  void   Fdscal(const unsigned*, const double*, const double*, const unsigned*);
  void   Fdcopy(const unsigned*, const double*, const unsigned*, double*, const unsigned*);
  void   Fdaxpy(const unsigned*, const double*, const double*, const unsigned*, double*, const unsigned*);
  // blas 2
  void   Fdgemv(const char*, const unsigned*, const unsigned*, const double*,
		  const double*, const unsigned*, const double*, const unsigned*,
	      const double*, double*, const unsigned*);
  // blas 3
  void Fdgemm(const char*, const char*, const unsigned*, const unsigned*,
	      const unsigned*, const double*, double*, const unsigned*,
	      double*, const unsigned*, const double*, double*,	const unsigned*);
  // lapack
  void Fdgesvd(const char*, const char*, const unsigned*, const unsigned*,
	       double*, const unsigned*, double*, double*, const unsigned*,
	       double*, const unsigned*, double*, const unsigned*, int*);
  void Fdgeqrf(const unsigned*, const unsigned*, double*, const unsigned*,
	       double*, double*, const unsigned*, int*);

  void Fdgeqp3(const unsigned*, const unsigned*, double*, const unsigned*, /*TYPE OF JPIV*/ unsigned*,
	       double*, double*, const unsigned*, int*);
  void Fdorgqr(const unsigned*, const unsigned*, const unsigned*,
	       double*, const unsigned*, double*, double*, const unsigned*, int*);
  void Fdormqr(const char*, const char*,
               const unsigned*, const unsigned*, const unsigned*,
	       const double*, const unsigned*, 
               double*, double*, const unsigned*, 
               double*, const unsigned*, int*);
  void Fdpotrf(const char*, const unsigned*, double*, const unsigned*, int*);

  // single //////////////////////////////////////////////////////////
  // blas 1
  float Fsdot(const unsigned*, const float*, const unsigned*,	const float*, const unsigned*);
  void Fsscal(const unsigned*, const float*, const float*, const unsigned*);
  void Fscopy(const unsigned*, const float*, const unsigned*,	float*, const unsigned*);
  void Fsaxpy(const unsigned*, const float*, const float*, const unsigned*, float*, const unsigned*);
  // blas 2
  void Fsgemv(const char*, const unsigned*, const unsigned*, const float*,
	      const float*, const unsigned*, const float*, const unsigned*,
	      const float*, float*, const unsigned*);
  // blas 3
  void Fsgemm(const char*, const char*, const unsigned*, const unsigned*,
	      const unsigned*, const float*, float*, const unsigned*,
	      float*, const unsigned*, const float*, float*, const unsigned*);
  // lapack
  void Fsgesvd(const char*, const char*, const unsigned*, const unsigned*,
	       float*, const unsigned*, float*, float*, const unsigned*,
	       float*, const unsigned*, float*, const unsigned*, int*);
  void Fsgeqrf(const unsigned*, const unsigned*, float*, const unsigned*,
	       float*, float*, const unsigned*, int*);
  void Fsgeqp3(const unsigned*, const unsigned*, float*, const unsigned*, /*TYPE OF JPIV*/ unsigned*,
	       float*, float*, const unsigned*, int*);
  void Fsorgqr(const unsigned*, const unsigned*, const unsigned*,
	       float*, const unsigned*, float*, float*, const unsigned*, int*);
  void Fsormqr(const char*, const char*,
               const unsigned*, const unsigned*, const unsigned*,
	       const float*, const unsigned*, 
               float*, float*, const unsigned*, 
               float*, const unsigned*, int*);
  void Fspotrf(const char*, const unsigned*, float*, const unsigned*, int*);

  // double complex //////////////////////////////////////////////////
  // blas 1
  double Fzdot(const unsigned*, const double*, const unsigned*, const double*, const unsigned*);  //added
  void Fzscal(const unsigned*, const double*, const double*, const unsigned*);
  void Fzcopy(const unsigned*, const double*, const unsigned*, double*, const unsigned*);
  void Fzaxpy(const unsigned*, const double*, const double*, const unsigned*, double*, const unsigned*);
  // blas 2
  void Fzgemv(const char*, const unsigned*, const unsigned*, const double*,
	      const double*, const unsigned*, const double*, const unsigned*,
	      const double*, double*, const unsigned*);
  // blas 3
  void Fzgemm(const char*, const char*, const unsigned*, const unsigned*,
	      const unsigned*, const double*, double*, const unsigned*,
	      double*, const unsigned*, const double*, double*, const unsigned*);
  void Fzgesvd(const char*, const char*, const unsigned*, const unsigned*,
	       double*, const unsigned*, double*, double*, const unsigned*,
	       double*, const unsigned*, double*,   unsigned*,  double*,   int*);  //change of input argument types (#13) WAS: int*  IS: unsigned*

  void Fzgeqrf(const unsigned*, const unsigned*, double*, const unsigned*,
	       double*, double*, const unsigned*, int*);
  void Fzgeqp3(const unsigned*, const unsigned*, double*, const unsigned*,/*TYPE OF JPIV*/ unsigned*,
	       double*, double*, const unsigned*, int*);

  void Fzpotrf(const char*, const unsigned*, double*, const unsigned*, int*);
  void Fzungqr(const unsigned*, const unsigned*, const unsigned*,
	       double*, const unsigned*, double*, double*, const unsigned*, int*);		// added zungqr
  // single complex //////////////////////////////////////////////////
  // blas 1
  void Fcscal(const unsigned*, const float*, const float*, const unsigned*);
  void Fccopy(const unsigned*, const float*, const unsigned*,	float*, const unsigned*);
  void Fcaxpy(const unsigned*, const float*, const float*, const unsigned*, float*, const unsigned*);
  // blas 2
  void Fcgemv(const char*, const unsigned*, const unsigned*, const float*,
	      const float*, const unsigned*, const float*, const unsigned*,
	      const float*, float*, const unsigned*);
  // blas 3
  void Fcgemm(const char*, const char*, const unsigned*, const unsigned*,
	      const unsigned*, const float*, float*, const unsigned*,
	      float*, const unsigned*, const float*, float*, const unsigned*);
  void Fcgeqrf(const unsigned*, const unsigned*, float*, const unsigned*,
	       float*, float*, const unsigned*, int*);
  void Fcgeqp3(const unsigned*, const unsigned*, float*, const unsigned*, /*TYPE OF JPIV*/ unsigned*,
	       float*, float*, const unsigned*, int*);    
  void Fcpotrf(const char*, const unsigned*, float*, const unsigned*, int*);


}


namespace FBlas {

  // copy
  inline void copy(const unsigned n, double* orig, double* dest)
  {	Fdcopy(&n, orig, &scalfmm::N_ONE, dest, &scalfmm::N_ONE);	}
  inline void copy(const unsigned n, const double* orig, double* dest)
  {	Fdcopy(&n, orig, &scalfmm::N_ONE, dest, &scalfmm::N_ONE);	}
  inline void copy(const unsigned n, float* orig, float* dest)
  {	Fscopy(&n, orig, &scalfmm::N_ONE, dest, &scalfmm::N_ONE);	}
  inline void copy(const unsigned n, const float* orig, float* dest)
  {   Fscopy(&n, orig, &scalfmm::N_ONE, dest, &scalfmm::N_ONE); }
  inline void c_copy(const unsigned n, double* orig, double* dest)
  {	Fzcopy(&n, orig, &scalfmm::N_ONE, dest, &scalfmm::N_ONE);	}
  inline void c_copy(const unsigned n, const double* orig, double* dest)
  {	Fzcopy(&n, orig, &scalfmm::N_ONE, dest, &scalfmm::N_ONE);	}
  inline void c_copy(const unsigned n, float* orig, float* dest)
  {	Fccopy(&n, orig, &scalfmm::N_ONE, dest, &scalfmm::N_ONE);	}
  inline void c_copy(const unsigned n, const float* orig, float* dest)
  {   Fccopy(&n, orig, &scalfmm::N_ONE, dest, &scalfmm::N_ONE); }

  // copy (variable increment)
  inline void copy(const unsigned n, double* orig, const unsigned inco, double* dest, const unsigned incd)
  {	Fdcopy(&n, orig, &inco, dest, &incd);	}
  inline void copy(const unsigned n, float* orig, const unsigned inco, float* dest, const unsigned incd)
  {	Fscopy(&n, orig, &inco, dest, &incd);	}
  inline void c_copy(const unsigned n, double* orig, const unsigned inco, double* dest, const unsigned incd)
  {	Fzcopy(&n, orig, &inco, dest, &incd);	}
  inline void c_copy(const unsigned n, float* orig, const unsigned inco, float* dest, const unsigned incd)
  {	Fccopy(&n, orig, &inco, dest, &incd);	}

  // scale
  inline void scal(const unsigned n, const double d, double* const x)
  {	Fdscal(&n, &d, x, &scalfmm::N_ONE); }
  inline void scal(const unsigned n, const float d, float* const x)
  {	Fsscal(&n, &d, x, &scalfmm::N_ONE); } 
  inline void c_scal(const unsigned n, const double* d, double* const x) //changed the input arg type for c_scal
  {	Fzscal(&n, d, x, &scalfmm::N_ONE); }
  inline void c_scal(const unsigned n, const float d, float* const x)
  {	Fcscal(&n, &d, x, &scalfmm::N_ONE); }

  // scale (variable increment)
  inline void scal(const unsigned n, const double d, double* const x, const unsigned incd)
  {	Fdscal(&n, &d, x, &incd); }
  inline void scal(const unsigned n, const float d, float* const x, const unsigned incd)
  {	Fsscal(&n, &d, x, &incd); }
  inline void c_scal(const unsigned n, const double d, double* const x, const unsigned incd)
  {	Fzscal(&n, &d, x, &incd); }
  inline void c_scal(const unsigned n, const float d, float* const x, const unsigned incd)
  {	Fcscal(&n, &d, x, &incd); }

  // set zero
  inline void setzero(const unsigned n, double* const x)
  {	for (unsigned i=0; i<n; ++i) x[i] = 0.0; }
  inline void setzero(const unsigned n, float* const x)
  {	for (unsigned i=0; i<n; ++i) x[i] = 0.0f; }
  inline void c_setzero(const unsigned n, double* const x)
  {	for (unsigned i=0; i<n; ++i) x[i*2] = x[i*2+1] = 0.0; }
  inline void c_setzero(const unsigned n, float* const x)
  {	for (unsigned i=0; i<n; ++i) x[i*2] = x[i*2+1] = 0.0f; }

  // y += x
  inline void add(const unsigned n, double* const x, double* const y)
  {	Fdaxpy(&n, &scalfmm::D_ONE, x, &scalfmm::N_ONE, y, &scalfmm::N_ONE);	}
  inline void add(const unsigned n, float* const x, float* const y)
  {	Fsaxpy(&n, &scalfmm::S_ONE, x, &scalfmm::N_ONE, y, &scalfmm::N_ONE);	}
  inline void c_add(const unsigned n, float* const x, float* const y)
  {	Fcaxpy(&n, scalfmm::C_ONE, x, &scalfmm::N_ONE, y, &scalfmm::N_ONE);	}
  inline void c_add(const unsigned n, double* const x,double* const y)
  {	Fzaxpy(&n, scalfmm::Z_ONE, x, &scalfmm::N_ONE, y, &scalfmm::N_ONE);	}

  // y += d x
  inline void axpy(const unsigned n, const double d, const double* const x, double* const y)
  {	Fdaxpy(&n, &d, x, &scalfmm::N_ONE, y, &scalfmm::N_ONE);	}
  inline void axpy(const unsigned n, const float d, const float* const x, float* const y)
  {	Fsaxpy(&n, &d, x, &scalfmm::N_ONE, y, &scalfmm::N_ONE);	}
  inline void c_axpy(const unsigned n, const float* d, const float* const x, float* const y)
  {	Fcaxpy(&n, d, x, &scalfmm::N_ONE, y, &scalfmm::N_ONE);	}
  inline void c_axpy(const unsigned n, const double* d, const double* const x, double* const y)
  {	Fzaxpy(&n, d, x, &scalfmm::N_ONE, y, &scalfmm::N_ONE);	}



  //	// y = d Ax
  //	inline void gemv(const unsigned m, const unsigned n, double d, double* A, double *x, double *y)
  //	{	cblas_dgemv(CblasColMajor, CblasNoTrans, m, n, d, A, m, x, scalfmm::N_ONE, scalfmm::D_ZERO, y, scalfmm::N_ONE); }
  //	inline void gemv(const unsigned m, const unsigned n, float d, float* A, float *x, float *y)
  //	{	cblas_sgemv(CblasColMajor, CblasNoTrans, m, n, d, A, m, x, scalfmm::N_ONE, scalfmm::S_ZERO, y, scalfmm::N_ONE); }
  // y = d Ax
  inline void gemv(const unsigned m, const unsigned n, double d, double* A, double *x, double *y)
  {	Fdgemv(scalfmm::JOB_STR, &m, &n, &d, A, &m, x, &scalfmm::N_ONE, &scalfmm::D_ZERO, y, &scalfmm::N_ONE); }
  inline void gemv(const unsigned m, const unsigned n, float d, float* A, float *x, float *y)
  {	Fsgemv(scalfmm::JOB_STR, &m, &n, &d, A, &m, x, &scalfmm::N_ONE, &scalfmm::S_ZERO, y, &scalfmm::N_ONE); }
  inline void c_gemv(const unsigned m, const unsigned n, float* d, float* A, float *x, float *y)
  {	Fcgemv(scalfmm::JOB_STR, &m, &n, d, A, &m, x, &scalfmm::N_ONE, scalfmm::C_ZERO, y, &scalfmm::N_ONE); }
  inline void c_gemv(const unsigned m, const unsigned n, double* d, double* A, double *x, double *y)
  {	Fzgemv(scalfmm::JOB_STR, &m, &n, d, A, &m, x, &scalfmm::N_ONE, scalfmm::Z_ZERO, y, &scalfmm::N_ONE); }

  //	// y += d Ax
  //	inline void gemva(const unsigned m, const unsigned n, double d, double* A, double *x, double *y)
  //	{	cblas_dgemv(CblasColMajor, CblasNoTrans, m, n, d, A, m, x, scalfmm::N_ONE, scalfmm::D_ONE, y, scalfmm::N_ONE); }
  //	inline void gemva(const unsigned m, const unsigned n, float d, float* A, float *x, float *y)
  //	{	cblas_sgemv(CblasColMajor, CblasNoTrans, m, n, d, A, m, x, scalfmm::N_ONE, scalfmm::S_ONE, y, scalfmm::N_ONE); }
  // y += d Ax
  inline void gemva(const unsigned m, const unsigned n, double d, double* A, double *x, double *y)
  {	Fdgemv(scalfmm::JOB_STR, &m, &n, &d, A, &m, x, &scalfmm::N_ONE, &scalfmm::D_ONE, y, &scalfmm::N_ONE);	}
  inline void gemva(const unsigned m, const unsigned n, float d, float* A, float *x, float *y)
  {	Fsgemv(scalfmm::JOB_STR, &m, &n, &d, A, &m, x, &scalfmm::N_ONE, &scalfmm::S_ONE, y, &scalfmm::N_ONE);	}
  inline void c_gemva(const unsigned m, const unsigned n, const float* d, const float* A, const float *x, float *y)
  {	Fcgemv(scalfmm::JOB_STR, &m, &n, d, A, &m, x, &scalfmm::N_ONE, scalfmm::C_ONE, y, &scalfmm::N_ONE);	}
  inline void c_gemva(const unsigned m, const unsigned n, const double* d, const double* A, const double *x, double *y)
  {	Fzgemv(scalfmm::JOB_STR, &m, &n, d, A, &m, x, &scalfmm::N_ONE, scalfmm::Z_ONE, y, &scalfmm::N_ONE);	}

  //	// y = d A^T x
  //	inline void gemtv(const unsigned m, const unsigned n, double d, double* A, double *x, double *y)
  //	{ cblas_dgemv(CblasColMajor, CblasTrans, m, n, d, A, m, x, scalfmm::N_ONE, scalfmm::D_ZERO, y, scalfmm::N_ONE); }
  //	inline void gemtv(const unsigned m, const unsigned n, float d, float* A, float *x, float *y)
  //	{	cblas_sgemv(CblasColMajor, CblasTrans, m, n, d, A, m, x, scalfmm::N_ONE, scalfmm::S_ZERO, y, scalfmm::N_ONE); }
  // y = d A^T x
  inline void gemtv(const unsigned m, const unsigned n, double d, double* A, double *x, double *y)
  {	Fdgemv(scalfmm::JOB_STR+1, &m, &n, &d, A, &m, x, &scalfmm::N_ONE, &scalfmm::D_ZERO, y, &scalfmm::N_ONE); }
  inline void gemtv(const unsigned m, const unsigned n, float d, float* A, float *x, float *y)
  {	Fsgemv(scalfmm::JOB_STR+1, &m, &n, &d, A, &m, x, &scalfmm::N_ONE, &scalfmm::S_ZERO, y, &scalfmm::N_ONE); }
  inline void c_gemtv(const unsigned m, const unsigned n, float* d, float* A, float *x, float *y)
  {	Fcgemv(scalfmm::JOB_STR+1, &m, &n, d, A, &m, x, &scalfmm::N_ONE, scalfmm::C_ZERO, y, &scalfmm::N_ONE); }
  inline void c_gemtv(const unsigned m, const unsigned n, double* d, double* A, double *x, double *y)
  {	Fzgemv(scalfmm::JOB_STR+1, &m, &n, d, A, &m, x, &scalfmm::N_ONE, scalfmm::Z_ZERO, y, &scalfmm::N_ONE); }
  inline void c_gemhv(const unsigned m, const unsigned n, float* d, float* A, float *x, float *y)
  {	Fcgemv(scalfmm::JOB_STR+7, &m, &n, d, A, &m, x, &scalfmm::N_ONE, scalfmm::C_ZERO, y, &scalfmm::N_ONE); } // hermitian transposed
  inline void c_gemhv(const unsigned m, const unsigned n, double* d, double* A, double *x, double *y)
  {	Fzgemv(scalfmm::JOB_STR+7, &m, &n, d, A, &m, x, &scalfmm::N_ONE, scalfmm::Z_ZERO, y, &scalfmm::N_ONE); } // hermitian transposed

  //	// y += d A^T x
  //	inline void gemtva(const unsigned m, const unsigned n, double d, double* A, double *x, double *y)
  //	{	cblas_dgemv(CblasColMajor, CblasTrans, m, n, d, A, m, x, scalfmm::N_ONE, scalfmm::D_ONE, y, scalfmm::N_ONE); }
  //	inline void gemtva(const unsigned m, const unsigned n, float d, float* A, float *x, float *y)
  //	{	cblas_sgemv(CblasColMajor, CblasTrans, m, n, d, A, m, x, scalfmm::N_ONE, scalfmm::S_ONE, y, scalfmm::N_ONE); }
  // y += d A^T x
  inline void gemtva(const unsigned m, const unsigned n, double d, double* A, double *x, double *y)
  {	Fdgemv(scalfmm::JOB_STR+1, &m, &n, &d, A, &m, x, &scalfmm::N_ONE, &scalfmm::D_ONE, y, &scalfmm::N_ONE);	}
  inline void gemtva(const unsigned m, const unsigned n, float d, float* A, float *x, float *y)
  {	Fsgemv(scalfmm::JOB_STR+1, &m, &n, &d, A, &m, x, &scalfmm::N_ONE, &scalfmm::S_ONE, y, &scalfmm::N_ONE);	}
  inline void c_gemtva(const unsigned m, const unsigned n, float* d, float* A, float *x, float *y)
  {	Fcgemv(scalfmm::JOB_STR+1, &m, &n, d, A, &m, x, &scalfmm::N_ONE, scalfmm::C_ONE, y, &scalfmm::N_ONE);	}
  inline void c_gemtva(const unsigned m, const unsigned n, double* d, double* A, double *x, double *y)
  {	Fzgemv(scalfmm::JOB_STR+1, &m, &n, d, A, &m, x, &scalfmm::N_ONE, scalfmm::Z_ONE, y, &scalfmm::N_ONE); }
  inline void c_gemhva(const unsigned m, const unsigned n, float* d, float* A, float *x, float *y)
  {	Fcgemv(scalfmm::JOB_STR+7, &m, &n, d, A, &m, x, &scalfmm::N_ONE, scalfmm::C_ONE, y, &scalfmm::N_ONE);	} // hermitian transposed
  inline void c_gemhva(const unsigned m, const unsigned n, double* d, double* A, double *x, double *y)
  {	Fzgemv(scalfmm::JOB_STR+7, &m, &n, d, A, &m, x, &scalfmm::N_ONE, scalfmm::Z_ONE, y, &scalfmm::N_ONE);	} // hermitian transposed




  // C = d A B, A is m x p, B is p x n
  inline void gemm(unsigned m, unsigned p, unsigned n, double d,
		   double* A, unsigned ldA, double* B, unsigned ldB, double* C, unsigned ldC)
  {	Fdgemm(scalfmm::JOB_STR, scalfmm::JOB_STR, &m, &n, &p, &d, A, &ldA, B, &ldB, &scalfmm::D_ZERO, C, &ldC);	}
  inline void gemm(unsigned m, unsigned p, unsigned n, float d,
		   float* A, unsigned ldA, float* B, unsigned ldB, float* C, unsigned ldC)
  {	Fsgemm(scalfmm::JOB_STR, scalfmm::JOB_STR, &m, &n, &p, &d, A, &ldA, B, &ldB, &scalfmm::S_ZERO, C, &ldC);	}
  inline void c_gemm(const unsigned m, const unsigned p, const unsigned n, const float* d,
		     float* A, const unsigned ldA, float* B, const unsigned ldB, float* C, const unsigned ldC)
  {
    Fcgemm(scalfmm::JOB_STR, scalfmm::JOB_STR, &m, &n, &p, d, A, &ldA, B, &ldB, scalfmm::C_ZERO, C, &ldC);	}
  inline void c_gemm(const unsigned m, const unsigned p, const unsigned n, const double* d,
		     double* A, const unsigned ldA, double* B, const unsigned ldB, double* C, const unsigned ldC)
  {
    Fzgemm(scalfmm::JOB_STR, scalfmm::JOB_STR, &m, &n, &p, d, A, &ldA, B, &ldB, scalfmm::Z_ZERO, C, &ldC);	}

  // C += d A B, A is m x p, B is p x n
  inline void gemma(unsigned m, unsigned p, unsigned n, double d,
		    double* A, unsigned ldA, double* B, unsigned ldB,	double* C, unsigned ldC)
  {	Fdgemm(scalfmm::JOB_STR, scalfmm::JOB_STR, &m, &n, &p, &d, A, &ldA, B, &ldB, &scalfmm::D_ONE, C, &ldC); }
  inline void gemma(unsigned m, unsigned p, unsigned n, float d,
		    float* A, unsigned ldA, float* B, unsigned ldB,	float* C, unsigned ldC)
  {	Fsgemm(scalfmm::JOB_STR, scalfmm::JOB_STR, &m, &n, &p, &d, A, &ldA, B, &ldB, &scalfmm::S_ONE, C, &ldC); }
  inline void c_gemma(unsigned m, unsigned p, unsigned n, float* d,
		      float* A, unsigned ldA, float* B, unsigned ldB,	float* C, unsigned ldC)
  {	Fcgemm(scalfmm::JOB_STR, scalfmm::JOB_STR, &m, &n, &p, d, A, &ldA, B, &ldB, scalfmm::C_ONE, C, &ldC); }
  inline void c_gemma(unsigned m, unsigned p, unsigned n, double* d,
		      double* A, unsigned ldA, double* B, unsigned ldB,	double* C, unsigned ldC)
  {	Fzgemm(scalfmm::JOB_STR, scalfmm::JOB_STR, &m, &n, &p, d, A, &ldA, B, &ldB, scalfmm::Z_ONE, C, &ldC); }

  // C = d A^T B, A is m x p, B is m x n
  inline void gemtm(unsigned m, unsigned p, unsigned n, double d,
		    double* A, unsigned ldA, double *B, unsigned ldB,	double* C, unsigned ldC)
  {	Fdgemm(scalfmm::JOB_STR+1, scalfmm::JOB_STR, &p, &n, &m, &d, A, &ldA, B, &ldB, &scalfmm::D_ZERO, C, &ldC);	}
  inline void gemtm(unsigned m, unsigned p, unsigned n, float d,
		    float* A, unsigned ldA, float *B, unsigned ldB,	float* C, unsigned ldC)
  {	Fsgemm(scalfmm::JOB_STR+1, scalfmm::JOB_STR, &p, &n, &m, &d, A, &ldA, B, &ldB, &scalfmm::S_ZERO, C, &ldC);	}
  inline void c_gemtm(unsigned m, unsigned p, unsigned n, float* d,
		      float* A, unsigned ldA, float *B, unsigned ldB,	float* C, unsigned ldC)
  {	Fcgemm(scalfmm::JOB_STR+1, scalfmm::JOB_STR, &p, &n, &m, d, A, &ldA, B, &ldB, scalfmm::C_ZERO, C, &ldC);	}
  inline void c_gemtm(unsigned m, unsigned p, unsigned n, double* d,
		      double* A, unsigned ldA, double *B, unsigned ldB,	double* C, unsigned ldC)
  {	Fzgemm(scalfmm::JOB_STR+1, scalfmm::JOB_STR, &p, &n, &m, d, A, &ldA, B, &ldB, scalfmm::Z_ZERO, C, &ldC);	}
  inline void c_gemhm(unsigned m, unsigned p, unsigned n, float* d, // hermitialn transposed
		      float* A, unsigned ldA, float *B, unsigned ldB,	float* C, unsigned ldC)
  {	Fcgemm(scalfmm::JOB_STR+7, scalfmm::JOB_STR, &p, &n, &m, d, A, &ldA, B, &ldB, scalfmm::C_ZERO, C, &ldC);	}
  inline void c_gemhm(unsigned m, unsigned p, unsigned n, double* d, // hermitian transposed
		      double* A, unsigned ldA, double *B, unsigned ldB,	double* C, unsigned ldC)
  {	Fzgemm(scalfmm::JOB_STR+7, scalfmm::JOB_STR, &p, &n, &m, d, A, &ldA, B, &ldB, scalfmm::Z_ZERO, C, &ldC);	}

  // C += d A^T B, A is m x p, B is m x n
  inline void gemtma(unsigned m, unsigned p, unsigned n, double d,
		     double* A, unsigned ldA, double *B, unsigned ldB, double* C, unsigned ldC)
  {	Fdgemm(scalfmm::JOB_STR+1, scalfmm::JOB_STR, &p, &n, &m, &d, A, &ldA, B, &ldB, &scalfmm::D_ONE, C, &ldC); }
  inline void gemtma(unsigned m, unsigned p, unsigned n, float d,
		     float* A, unsigned ldA, float *B, unsigned ldB, float* C, unsigned ldC)
  {	Fsgemm(scalfmm::JOB_STR+1, scalfmm::JOB_STR, &p, &n, &m, &d, A, &ldA, B, &ldB, &scalfmm::S_ONE, C, &ldC); }
  inline void c_gemtma(unsigned m, unsigned p, unsigned n, float* d,
		       float* A, unsigned ldA, float *B, unsigned ldB, float* C, unsigned ldC)
  {	Fcgemm(scalfmm::JOB_STR+1, scalfmm::JOB_STR, &p, &n, &m, d, A, &ldA, B, &ldB, scalfmm::C_ONE, C, &ldC); }
  inline void c_gemtma(unsigned m, unsigned p, unsigned n, double* d,
		       double* A, unsigned ldA, double *B, unsigned ldB, double* C, unsigned ldC)
  {	Fzgemm(scalfmm::JOB_STR+1, scalfmm::JOB_STR, &p, &n, &m, d, A, &ldA, B, &ldB, scalfmm::Z_ONE, C, &ldC); }
  inline void c_gemhma(unsigned m, unsigned p, unsigned n, float* d, // hermitian transposed
		       float* A, unsigned ldA, float *B, unsigned ldB, float* C, unsigned ldC)
  {	Fcgemm(scalfmm::JOB_STR+7, scalfmm::JOB_STR, &p, &n, &m, d, A, &ldA, B, &ldB, scalfmm::C_ONE, C, &ldC); }
  inline void c_gemhma(unsigned m, unsigned p, unsigned n, double* d, // hermitian transposed
		       double* A, unsigned ldA, double *B, unsigned ldB, double* C, unsigned ldC)
  {	Fzgemm(scalfmm::JOB_STR+7, scalfmm::JOB_STR, &p, &n, &m, d, A, &ldA, B, &ldB, scalfmm::Z_ONE, C, &ldC); }
	
	
  // C = d A B^T, A is m x p, B is n x p
  inline void gemmt(unsigned m, unsigned p, unsigned n, double d,
		    double* A, unsigned ldA, double *B, unsigned ldB, double* C, unsigned ldC)
  {	Fdgemm(scalfmm::JOB_STR, scalfmm::JOB_STR+1, &m, &n, &p, &d, A, &ldA, B, &ldB, &scalfmm::D_ZERO, C, &ldC);	}
  inline void gemmt(unsigned m, unsigned p, unsigned n, float d,
		    float* A, unsigned ldA, float *B, unsigned ldB,	float* C, unsigned ldC)
  {	Fsgemm(scalfmm::JOB_STR, scalfmm::JOB_STR+1, &m, &n, &p, &d, A, &ldA, B, &ldB, &scalfmm::S_ZERO, C, &ldC);	}
  inline void c_gemmt(unsigned m, unsigned p, unsigned n, float d,
		      float* A, unsigned ldA, float *B, unsigned ldB, float* C, unsigned ldC)
  {	Fcgemm(scalfmm::JOB_STR, scalfmm::JOB_STR+1, &m, &n, &p, &d, A, &ldA, B, &ldB, scalfmm::C_ZERO, C, &ldC); }
  inline void c_gemmt(unsigned m, unsigned p, unsigned n, double d,
		      double* A, unsigned ldA, double *B, unsigned ldB,	double* C, unsigned ldC)
  {	Fzgemm(scalfmm::JOB_STR, scalfmm::JOB_STR+1, &m, &n, &p, &d, A, &ldA, B, &ldB, scalfmm::Z_ZERO, C, &ldC);	}
  inline void c_gemmh(unsigned m, unsigned p, unsigned n, float d, // hermitian transposed
		      float* A, unsigned ldA, float *B, unsigned ldB, float* C, unsigned ldC)
  {	Fcgemm(scalfmm::JOB_STR, scalfmm::JOB_STR+7, &m, &n, &p, &d, A, &ldA, B, &ldB, scalfmm::C_ZERO, C, &ldC); }
  inline void c_gemmh(unsigned m, unsigned p, unsigned n, double d, // hermitian transposed
		      double* A, unsigned ldA, double *B, unsigned ldB,	double* C, unsigned ldC)
  {	Fzgemm(scalfmm::JOB_STR, scalfmm::JOB_STR+7, &m, &n, &p, &d, A, &ldA, B, &ldB, scalfmm::Z_ZERO, C, &ldC);	}

  // C += d A B^T, A is m x p, B is n x p
  inline void gemmta(unsigned m, unsigned p, unsigned n, double d,
		     double* A, unsigned ldA, double *B, unsigned ldB, double* C, unsigned ldC)
  {	Fdgemm(scalfmm::JOB_STR, scalfmm::JOB_STR+1, &m, &n, &p, &d, A, &ldA, B, &ldB, &scalfmm::D_ONE, C, &ldC); }
  inline void gemmta(unsigned m, unsigned p, unsigned n, float d,
		     float* A, unsigned ldA, float *B, unsigned ldB, float* C, unsigned ldC)
  {	Fsgemm(scalfmm::JOB_STR, scalfmm::JOB_STR+1, &m, &n, &p, &d, A, &ldA, B, &ldB, &scalfmm::S_ONE, C, &ldC); }
  inline void c_gemmta(unsigned m, unsigned p, unsigned n, float* d,
		       float* A, unsigned ldA, float *B, unsigned ldB, float* C, unsigned ldC)
  {	Fcgemm(scalfmm::JOB_STR, scalfmm::JOB_STR+1, &m, &n, &p, d, A, &ldA, B, &ldB, scalfmm::C_ONE, C, &ldC); }
  inline void c_gemmta(unsigned m, unsigned p, unsigned n, double* d,
		       double* A, unsigned ldA, double *B, unsigned ldB, double* C, unsigned ldC)
  {	Fzgemm(scalfmm::JOB_STR, scalfmm::JOB_STR+1, &m, &n, &p, d, A, &ldA, B, &ldB, scalfmm::Z_ONE, C, &ldC); }
  inline void c_gemmha(unsigned m, unsigned p, unsigned n, float* d, // hermitian transposed
		       float* A, unsigned ldA, float *B, unsigned ldB, float* C, unsigned ldC)
  {	Fcgemm(scalfmm::JOB_STR, scalfmm::JOB_STR+7, &m, &n, &p, d, A, &ldA, B, &ldB, scalfmm::C_ONE, C, &ldC); }
  inline void c_gemmha(unsigned m, unsigned p, unsigned n, double* d, // hermitian transposed
		       double* A, unsigned ldA, double *B, unsigned ldB, double* C, unsigned ldC)
  {	Fzgemm(scalfmm::JOB_STR, scalfmm::JOB_STR+7, &m, &n, &p, d, A, &ldA, B, &ldB, scalfmm::Z_ONE, C, &ldC); }


  // singular value decomposition
  //
  inline int gesvd(unsigned m, unsigned n, double* A, double* S, double* VT, unsigned ldVT,
		   unsigned nwk, double* wk)
  {
    int INF;
    Fdgesvd(scalfmm::JOB_STR+2, scalfmm::JOB_STR+3, &m, &n, A, &m, S, A, &m, VT, &ldVT,	wk, &nwk, &INF);
    return INF;
  }
  
  
  
  
  
  // c_gesvd
 //------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  //
  //    A = U * SIGMA * conjugate-transpose(V)
  // scalfmm::JOB_STR+2 = 'O':  the first min(m,n) columns of U (the left singular vectors) are overwritten on the array A;
  inline int c_gesvd(unsigned m, unsigned n, double* A, double* S, double* VT, unsigned ldVT,
		     unsigned nwk, double* wk,double* rwk)
  {
    int INF;
    Fzgesvd(scalfmm::JOB_STR+2, scalfmm::JOB_STR+3, &m, &n, A, &m, S, A, &m, VT, &ldVT,	wk, &nwk, rwk, &INF);
    return INF;
  }

 // gesvd
 //------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
 
  inline int gesvd(unsigned m, unsigned n, float* A, float* S, float* VT, unsigned ldVT,
		   unsigned nwk, float* wk)
  {
    int INF;
    Fsgesvd(scalfmm::JOB_STR+2, scalfmm::JOB_STR+3, &m, &n, A, &m, S, A, &m, VT, &ldVT,	wk, &nwk, &INF);
    return INF;
  }
   
  //------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  
  
  

  // singular value decomposition (SO)
  inline int gesvdSO(unsigned m, unsigned n, double* A, double* S, double* U, unsigned ldU,
		     unsigned nwk, double* wk)
  {
    int INF;
    Fdgesvd(scalfmm::JOB_STR+3, scalfmm::JOB_STR+2, &m, &n, A, &m, S, U, &m, A, &ldU, wk, &nwk, &INF);
    return INF;
  }
  inline int gesvdSO(unsigned m, unsigned n, float* A, float* S, float* U, unsigned ldU,
		     unsigned nwk, float* wk)
  {
    int INF;
    Fsgesvd(scalfmm::JOB_STR+3, scalfmm::JOB_STR+2, &m, &n, A, &m, S, U, &m, A, &ldU, wk, &nwk, &INF);
    return INF;
  }

  // singular value decomposition (AA)
  inline int gesvdAA(unsigned m, unsigned n, double* A, double* S, double* U, unsigned ldU,
		     unsigned nwk, double* wk)
  {
    int INF;
    Fdgesvd("A", "A", &m, &n, A, &m, S, U, &m, A, &ldU, wk, &nwk, &INF);
    return INF;
  }
  inline int gesvdAA(unsigned m, unsigned n, float* A, float* S, float* U, unsigned ldU,
		     unsigned nwk, float* wk)
  {
    int INF;
    Fsgesvd("A", "A", &m, &n, A, &m, S, U, &m, A, &ldU, wk, &nwk, &INF);
    return INF;
  }

  // Scalar product v1'*v2
  inline double scpr(const unsigned n, const double* const v1, const double* const v2)
  {	return Fddot(&n, v1, &scalfmm::N_ONE, v2, &scalfmm::N_ONE); }
  inline double c_scpr(const unsigned n, const double* const v1, const double* const v2)
  {	return Fzdot(&n, v1, &scalfmm::N_ONE, v2, &scalfmm::N_ONE); }  
  inline float scpr(const unsigned n, const float* const v1, const float* const v2)
  {	return Fsdot(&n, v1, &scalfmm::N_ONE, v2, &scalfmm::N_ONE);	}



  // QR factorisation
  inline int geqrf(const unsigned m, const unsigned n, double* A, double* tau, unsigned nwk, double* wk)
  {
    int INF;
    Fdgeqrf(&m, &n, A, &m, tau, wk, &nwk, &INF);
    return INF;
  }
 
  inline int geqrf(const unsigned m, const unsigned n, float* A, float* tau, unsigned nwk, float* wk)
  {
    int INF;
    Fsgeqrf(&m, &n, A, &m, tau, wk, &nwk, &INF);
    return INF;
  }

  // QR factorisation with column pivoting
  inline int geqp3(const unsigned m, const unsigned n, double* A, unsigned* jpiv, double* tau, unsigned nwk, double* wk)
  {
    int INF;
    Fdgeqp3(&m, &n, A, &m, jpiv, tau, wk, &nwk, &INF);
    return INF;
  }
  inline int geqp3(const unsigned m, const unsigned n, float* A, unsigned* jpiv, float* tau, unsigned nwk, float* wk)
  {
    int INF;
    Fsgeqp3(&m, &n, A, &m, jpiv, tau, wk, &nwk, &INF);
    return INF;
  }	
  inline int c_geqrf(const unsigned m, const unsigned n, float* A, float* tau, unsigned nwk, float* wk)
  {
    int INF;
    Fcgeqrf(&m, &n, A, &m, tau, wk, &nwk, &INF);
    return INF;
  }
	
  inline int c_geqrf(const unsigned m, const unsigned n, double* A, double* tau, unsigned nwk, double* wk)
  {
    int INF;
    Fzgeqrf(&m, &n, A, &m, tau, wk, &nwk, &INF);
    return INF;
  }
  inline int c_geqp3(const unsigned m, const unsigned n, float* A, unsigned* jpiv, float* tau, unsigned nwk, float* wk)
  {
    int INF;
    Fcgeqp3(&m, &n, A, &m, jpiv, tau, wk, &nwk, &INF);
    return INF;
  }
    
  inline int c_geqp3(const unsigned m, const unsigned n, double* A, unsigned* jpiv, double* tau, unsigned nwk, double* wk)
  {
    int INF;
    Fzgeqp3(&m, &n, A, &m, jpiv, tau, wk, &nwk, &INF);
    return INF;
  }

  // return full of Q-Matrix (QR factorization) in A
  inline int orgqr_full(const unsigned m, const unsigned n, double* A, double* tau, unsigned nwk, double* wk)
  {
    int INF;
    Fdorgqr(&m, &m, &n, A, &m, tau, wk, &nwk, &INF);
    return INF;
  }
  inline int orgqr_full(const unsigned m, const unsigned n, float* A, float* tau, unsigned nwk, float* wk)
  {
    int INF;
    Fsorgqr(&m, &m, &n, A, &m, tau, wk, &nwk, &INF);
    return INF;
  }
  // return the leading n columns of Q-Matrix (QR factorization) in A
  inline int orgqr(const unsigned m, const unsigned n, double* A, double* tau, unsigned nwk, double* wk)
  {
    int INF;
    Fdorgqr(&m, &n, &n, A, &m, tau, wk, &nwk, &INF);
    return INF;
  }
  inline int c_orgqr(const unsigned m, const unsigned n, double* A, double* tau, unsigned nwk, double* wk)			//added
  {
    int INF;
    Fzungqr(&m, &n, &n, A, &m, tau, wk, &nwk, &INF);
    return INF;
  }
  inline int orgqr(const unsigned m, const unsigned n, float* A, float* tau, unsigned nwk, float* wk)
  {
    int INF;
    Fsorgqr(&m, &n, &n, A, &m, tau, wk, &nwk, &INF);
    return INF;
  }



  // apply Q-Matrix (from QR factorization) to C
  // LEFT: Q(^T)C
  inline int left_ormqr(const char* TRANS, const unsigned m, const unsigned n, const double* A, double* tau, double* C, unsigned nwk, double* wk)
  {
    int INF;
    Fdormqr("L", TRANS, &m, &n, &m, A, &m, tau, C, &m, wk, &nwk, &INF);
    return INF;
  }

  inline int left_ormqr(const char* TRANS, const unsigned m, const unsigned n, const float* A, float* tau, float* C, unsigned nwk, float* wk)
  {
    int INF;
    Fsormqr("L", TRANS, &m, &n, &m, A, &m, tau, C, &m, wk, &nwk, &INF);
    return INF;
  }
  // RIGHT: CQ(^T)
  inline int right_ormqr(const char* TRANS, const unsigned m, const unsigned n, const double* A, double* tau, double* C, unsigned nwk, double* wk)
  {
    int INF;
    Fdormqr("R", TRANS, &m, &n, &n, A, &n, tau, C, &m, wk, &nwk, &INF);
    return INF;
  }


  // Cholesky decomposition: A=LL^T (if A is symmetric definite positive)
  inline int potrf(const unsigned m, double* A, const unsigned n)
  { 
    int INF;  
    Fdpotrf("L", &m, A, &n, &INF);
    return INF;
  }
  inline int potrf(const unsigned m, float* A, const unsigned n)
  { 
    int INF;  
    Fspotrf("L", &m, A, &n, &INF);
    return INF;
  }

} // end namespace FCBlas

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////


//#else
//enum CBLAS_ORDER {CblasRowMajor=101, CblasColMajor=102};
//enum CBLAS_TRANSPOSE {CblasNoTrans, CblasTrans, CblasConjTrans};
//
//template <typename T>
//void cblas_gemv(const CBLAS_ORDER order ,
//								const CBLAS_TRANSPOSE TransA , const int M , const int N ,
//								const void *alpha , const void *A , const int lda ,
//								const void *X , const int incX , const void *beta ,
//								void *Y , const int incY){
//}
//template <typename T>
//void cblas_dotu_sub( const int N , const void *X , const int incX ,
//										 const void *Y , const int incY , void *dotu){
//}


#endif //FBLAS_HPP

