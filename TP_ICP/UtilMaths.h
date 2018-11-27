////////////////////////////////////////////////////////////////////////
//  UtilsMaths.h -- Interface of routines on matrices
// 
// Modified : Raoul de Charette --- 2/2010
//

#ifndef UTILMATHS_H
#define UTILMATHS_H


#include <stdio.h>
#include <math.h>

#define TRUE 1
#define FALSE 0

#define ASSERT(x)	((x==1)? 1:0)	

#define SQR(x)			((x) * (x))
#define SIGN2(x,y)	    ((y)<0? -fabs(x):fabs(x))
#define SWAP(type, x, y)	{ type temp = (x); (x) = (y); (y) = temp; }
#define PT_PT_DIST(Pt1, Pt2)    sqrt(SQR(Pt1[0] - Pt2[0]) + \
                                     SQR(Pt1[1] - Pt2[1]) + \
                                     SQR(Pt1[2] - Pt2[2]))

// Matrix operation routines


//  This routine permutes the elements of the real symmetric matrix 'Mat'
//  so that all the larger elements in the matrix are relocated towards
//  the lower-right portion of the matrix.
//  This is done to minimize round-off errors in the tred2 and tqli algorithms,
//	for which the scatter matrix is being computed as an input.
//	"Order" is a 3 element array that will keep track of the order
//	changes in the matrix. The ordered matrix is returned in 'Mat'.
void Permute_Matrix(double** Mat, int* Order, int n);


//	Here are the routines for computing eigen values and eigen
//	vectors.
//	The routines are from the book "Numerical Recipes in C".
//	Include:
//			tred2(a,n,d,e)
//			tqli(d,e,n,z)

//	tred2 Householder reduction of a real, symmetric matrix a[1..n][1..n].
//	On output, a is replaced by the orthogonal matrix q effecting the
//	transformation. d[1..n] returns the diagonal elements of the
//	tridiagonal matrix, and e[1..n] the off-diagonal elements, with
//	e[1]=0.
//  (Notes:	*Here n=3 (simplified version of the NumRecInC tred2 and tqli algorithms),
//			*in the book, the index for array starts from 1 whereas in C, index should start from zero.
//			In order to cope with this, all indices [i] have been replaced by [i-1].)
void tred2(double** a, double* d, double* e, int n);

//	QL algo with implicit shift, to determine the eigenvalues and
//	eigenvectors of a real,symmetric  tridiagonal matrix, or of a real,
//	symmetric matrix previously reduced by algo tred2.
//	On input , d[1..n] contains the diagonal elements of the tridiagonal
//	matrix. On output, it returns the eigenvalues. The vector e[1..n]
//	inputs the subdiagonal elements of the tridiagonal matrix, with e[1]
//	arbitrary. On output e is destroyed. If the eigenvectors of a
//	tridiagonal matrix are desired, the matrix z[1..n][1..n] is input
//	as the identity matrix. If the eigenvectors of a matrix that has
//	been reduced by tred2 are required, then z is input as the matrix
//	output by tred2. In either case, the kth column of z returns the
//	normalized eigenvector corresponding to d[k].
void tqli(double* d, double* e, double** z, int n);

// Resolution of the system Ax=b  by the Gauss method
//
void RegularGaussj(const int n, double* b, double** A);

// Given a positive-definite symmetric matrix a[1..n][1..n], this routine constructs its Cholesky
// decomposition, A= L L^T.
int choldc(double **a, int n, double *p);

// Solves the set of linear equations a.x=b, where a is a positive-definite symmetric matrix
// a[1..n][1..n] and p[1..n] are input as the output of the routine choldc.
void cholsl(double **a, int n, double *p, double *b, double *x);

// The following function computes the eigen elements of a 3x3
// real symmetric matrix. The results are ordered as follows:
// EigenValues  = [lambda_max, lambda_mid, lambda_min]
// and EigenVectors = [v_max, v_mid, v_min].
void SymMatr_EigenAnalysis(double** Matrix, double* EigenVals, double** EigenVects, int n);

//Given a matrix a[1..m][1..n], this routine computes its singular value decomposition, A =
//U·W·V^t. The matrix U replaces a on output. The diagonal matrix of singular values W is output
//as a vector w[1..n]. The matrix V (not the transpose V^t ) is output as v[1..n][1..n].
void svdcmp(float **a, int m, int n, float w[], float **v);
//Computes (a2 + b2)1/2 without destructive underflow or overflow.
float pythag(float a, float b);

#endif  //UTILMATHS_H
