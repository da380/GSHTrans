#ifndef GSH_TRANS_BLAS_GUARD_H
#define GSH_TRANS_BLAS_GUARD_H

// The one BLAS entry point this library uses, and the row-major wrapper the
// matrix kernel calls.
//
// Present only where GSHTRANS_HAVE_BLAS is defined; the header is conditional
// in its entirety, so including it without a BLAS is not an error and gives
// nothing. That is the shape RadialSplineDerivative.h and RadialResample.h
// already use for the optional Interpolation dependency.
//
// -- The Fortran interface, declared here rather than through cblas.h.
// find_package(BLAS) finds a library, not a header, and where cblas.h lives
// varies between OpenBLAS, MKL and Accelerate -- so depending on it would
// make the build fragile in exchange for an argument order. The Fortran
// symbol is the one thing every BLAS agrees on.

#ifdef GSHTRANS_HAVE_BLAS

#include <concepts>
#include <cstddef>

extern "C" {
void dgemm_(const char* transa, const char* transb, const int* m,
            const int* n, const int* k, const double* alpha, const double* a,
            const int* lda, const double* b, const int* ldb,
            const double* beta, double* c, const int* ldc);

void sgemm_(const char* transa, const char* transb, const int* m,
            const int* n, const int* k, const float* alpha, const float* a,
            const int* lda, const float* b, const int* ldb, const float* beta,
            float* c, const int* ldc);
}

namespace GSHTrans {

namespace BlasDetails {

// The precisions a BLAS has. There are four routines -- s, d, c and z -- and
// no wider one: `long double` is not part of the interface and never has
// been, so a grid over it cannot use the matrix kernel however it is built.
// GaussLegendreGrid refuses that combination at construction, where Real is
// known, rather than failing to compile inside a call.
template <typename T>
concept BlasReal = std::same_as<T, float> || std::same_as<T, double>;

inline void Gemm(const char* transa, const char* transb, int m, int n, int k,
                 double alpha, const double* a, int lda, const double* b,
                 int ldb, double beta, double* c, int ldc) {
  dgemm_(transa, transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
}

inline void Gemm(const char* transa, const char* transb, int m, int n, int k,
                 float alpha, const float* a, int lda, const float* b, int ldb,
                 float beta, float* c, int ldc) {
  sgemm_(transa, transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
}

// C(rows x cols) = alpha * A(rows x inner) * B(inner x cols) + beta * C, with
// every operand row-major and its own leading dimension.
//
// BLAS is column-major, so this is the standard reversal: a row-major array
// read column-major is its own transpose, and (A B)^T = B^T A^T. Passing B
// first, A second, and exchanging the roles of rows and cols therefore
// computes the row-major product with no data movement at all. Written once
// here because getting it wrong is easy and getting it wrong silently is
// easier -- the shapes still conform when the operands are square.
template <typename Real>
void RowMajorGemm(int rows, int cols, int inner, Real alpha, const Real* a,
                  int lda, const Real* b, int ldb, Real beta, Real* c,
                  int ldc) {
  Gemm("N", "N", cols, rows, inner, alpha, b, ldb, a, lda, beta, c, ldc);
}

// C(rows x cols) = alpha * A^T * B(inner x cols) + beta * C, where the array
// `a` holds A itself -- (inner x rows), row-major, leading dimension lda --
// and it is A's transpose that multiplies.
//
// This is what lets one stored Wigner matrix serve both directions. The
// forward transform multiplies by D and the inverse by D^T, and neither
// copies nor stores a second matrix: the transpose flag is the whole of the
// difference between them.
template <typename Real>
void RowMajorGemmTransposed(int rows, int cols, int inner, Real alpha,
                            const Real* a, int lda, const Real* b, int ldb,
                            Real beta, Real* c, int ldc) {
  Gemm("N", "T", cols, rows, inner, alpha, b, ldb, a, lda, beta, c, ldc);
}

}  // namespace BlasDetails

}  // namespace GSHTrans

#endif  // GSHTRANS_HAVE_BLAS

#endif  // GSH_TRANS_BLAS_GUARD_H
