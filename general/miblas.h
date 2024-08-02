/*
code developed by
    Rong-Qiang He (rqhe@ruc.edu.cn, RUC, China)
date 2013 - 2017
*/

#ifndef _MIBLAS_H_
#define _MIBLAS_H_

#define MKL_Complex8  std::complex<float>
#define MKL_Complex16 std::complex<double>
#include <mkl.h>

// dot product lhs^+ dot rhs
inline double DOT(const Vec<double> &lhs, const Vec<double> &rhs)
{
#ifdef _CHECK_DIMENSION_MATCH_
	ASSERT_EQ(lhs.size(), rhs.size());
#endif
	double *X = lhs.p();
	double *Y = rhs.p();
	MKL_INT N = (MKL_INT)lhs.size();
	MKL_INT incX = 1;
	MKL_INT incY = 1;
	return cblas_ddot(N, X, incX, Y, incY);
}

// dot product lhs^+ dot rhs
inline MKL_Complex16 DOT(const Vec<MKL_Complex16> &lhs, const Vec<MKL_Complex16> &rhs)
{
#ifdef _CHECK_DIMENSION_MATCH_
	ASSERT_EQ(lhs.size(), rhs.size());
#endif
	MKL_Complex16 *X = lhs.p();
	MKL_Complex16 *Y = rhs.p();
	MKL_INT N = (MKL_INT)lhs.size();
	MKL_INT incX = 1;
	MKL_INT incY = 1;
	MKL_Complex16 dotc;
	cblas_zdotc_sub(N, (void *)X, incX, (void *)Y, incY, (void *)(&dotc));
	return dotc;
}

// matrix-vector multiplication y = a * x
inline void MUL(Vec<double> &y, const Mat<double> &a, const Vec<double> &x)
{
#ifdef _CHECK_DIMENSION_MATCH_
	ASSERT_EQ2(y.size(), a.nrows(), a.ncols(), x.size());
#endif
	if (!y.size()) return;
	if (!x.size()) { y = 0.; return; }
	MKL_INT M = (MKL_INT)a.nrows();
	MKL_INT N = (MKL_INT)a.ncols();
	double *A = a.p();
	double *X = x.p();
	double *Y = y.p();
	double alpha = 1.;
	double beta = 0.;
	MKL_INT incX = 1;
	MKL_INT incY = 1;
	cblas_dgemv(CblasRowMajor, CblasNoTrans, M, N, alpha, A, N, X, incX, beta, Y, incY);
}

// matrix-vector multiplication y = a * x
inline void MUL(Vec<MKL_Complex16> &y, const Mat<MKL_Complex16> &a, const Vec<MKL_Complex16> &x)
{
#ifdef _CHECK_DIMENSION_MATCH_
	ASSERT_EQ2(y.size(), a.nrows(), a.ncols(), x.size());
#endif
	if (!y.size()) return;
	if (!x.size()) { y = 0.; return; }
	MKL_INT M = (MKL_INT)a.nrows();
	MKL_INT N = (MKL_INT)a.ncols();
	MKL_Complex16 *A = a.p();
	MKL_Complex16 *X = x.p();
	MKL_Complex16 *Y = y.p();
	MKL_Complex16 alpha = 1.;
	MKL_Complex16 beta = 0.;
	MKL_INT incX = 1;
	MKL_INT incY = 1;
	cblas_zgemv(CblasRowMajor, CblasNoTrans, 
		M, N, (void *)(&alpha), (void *)A, N, (void *)X, incX, (void *)(&beta), (void *)Y, incY);
}

// matrix-matrix multiplication a = b * c
inline void MUL(Mat<double> &a, const Mat<double> &b, const Mat<double> &c)
{
#ifdef _CHECK_DIMENSION_MATCH_
	ASSERT_EQ3(a.nrows(), b.nrows(), b.ncols(), c.nrows(), c.ncols(), a.ncols());
#endif
	if (!a.nrows() || !a.ncols()) return;
	if (!b.ncols()) { a = 0.; return; }
	MKL_INT M = (MKL_INT)a.nrows();
	MKL_INT N = (MKL_INT)a.ncols();
	MKL_INT P = (MKL_INT)b.ncols();
	double alpha = 1.;
	double beta = 0.;
	double *A = b.p();
	double *B = c.p();
	double *C = a.p();
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 
		M, N, P, alpha, A, P, B, N, beta, C, N);
}

// matrix-matrix multiplication a = b * c
inline void MUL(Mat<MKL_Complex16> &a, const Mat<MKL_Complex16> &b, const Mat<MKL_Complex16> &c)
{
#ifdef _CHECK_DIMENSION_MATCH_
	ASSERT_EQ3(a.nrows(), b.nrows(), b.ncols(), c.nrows(), c.ncols(), a.ncols());
#endif
	if (!a.nrows() || !a.ncols()) return;
	if (!b.ncols()) { a = 0.; return; }
	MKL_INT M = (MKL_INT)a.nrows();
	MKL_INT N = (MKL_INT)a.ncols();
	MKL_INT P = (MKL_INT)b.ncols();
	MKL_Complex16 *A = b.p();
	MKL_Complex16 *B = c.p();
	MKL_Complex16 *C = a.p();
	MKL_Complex16 alpha = 1.;
	MKL_Complex16 beta = 0.;
	cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 
		M, N, P, (void *)(&alpha), (void *)A, P, (void *)B, N, (void *)(&beta), (void *)C, N);
}

inline VecCmplx operator*(const MatReal &a, const VecCmplx &x) { return cmplx(a * real(x), a * imag(x)); }
inline VecCmplx operator*(const MatCmplx &a, const VecReal &x) { return a * cmplx(x); }


//	matrix-matrix multiplication a = alpha * b * c + beta * a
inline void MUL(MKL_Complex16 beta, Mat<MKL_Complex16> &a, MKL_Complex16 alpha, const Mat<MKL_Complex16> &b, const Mat<MKL_Complex16> &c)
{
#ifdef _CHECK_DIMENSION_MATCH_
	ASSERT_EQ3(a.nrows(), b.nrows(), b.ncols(), c.nrows(), c.ncols(), a.ncols());
#endif
	if (!a.nrows() || !a.ncols()) return;
	if (!b.ncols()) { a = 0.; return; }
	MKL_INT M = (MKL_INT)a.nrows();
	MKL_INT N = (MKL_INT)a.ncols();
	MKL_INT P = (MKL_INT)b.ncols();
	MKL_Complex16 *A = b.p();
	MKL_Complex16 *B = c.p();
	MKL_Complex16 *C = a.p();
	cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
		M, N, P, (void *)(&alpha), (void *)A, P, (void *)B, N, (void *)(&beta), (void *)C, N);
}

#endif /* _MIBLAS_H_ */
