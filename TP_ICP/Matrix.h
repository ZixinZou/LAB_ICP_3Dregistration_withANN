//
//  Matrix.h -- Matrix class
//
//
// Author : thomas Chaperon --- 04/2001
//
// Modified : xavier Brun --- 11/2004
// Modified : Raoul de Charette --- 2/2010
//////////////////////////////////////////////////////////////////////

#ifndef MATRIX_H
#define MATRIX_H

#include "UtilMaths.h"
#include "Vector.h"
class Vector;

//#include <vector>
using namespace std;


class Matrix {
public:
	// constructors, destructor
	//
	Matrix();
	Matrix(const int& aNbRows, const int& aNbCols, const bool& aIsSymmetric = false);
	//rectangular matrix, symmetric (if square) or not. Default is non-symmetric.
	Matrix(const int& aSize, const bool& aIsSymmetric = false); // square matrix (symmetric or not)
	Matrix(const vector<Vector>& aColVectors); // copy constructor from a list of column vectors
	Matrix(const Matrix& aM);// copy constructor
	virtual ~Matrix();


	void	SetToZero();
	void	Resize(const int& aNbRows, const int& aNbCols);

	int		GetNbRow() const { return mNbRow; }
	int		GetNbCol() const { return mNbCol; }
	int		GetSize() const { return mSize; } // size of a square matrix (0 if the matrix is not square)


	double& operator() (const int& aRow, const int& aCol) { return mElt[aRow*mNbCol + aCol]; }
	// Fetch a reference inside the Matrix
	const double& operator() (const int& aRow, const int& aCol) const { return mElt[aRow*mNbCol + aCol]; }
	// Fetch a const reference inside the Matrix


	const double&	GetElt(const int& aRow, const int& aCol) const { return mElt[aRow*mNbCol + aCol]; }
	void	SetElt(const int& aRow, const int& aCol, const double& aK) { mElt[aRow*mNbCol + aCol] = aK; }


	Vector	GetRow(const int& aRow) const;
	Vector	GetCol(const int& aCol) const;
	void	SetRow(const int& aRow, const Vector& aV);
	void	SetCol(const int& aCol, const Vector& aV);


	// CAUTION this only works if the matrix has already been constructed
	void	operator=(const Matrix& aM);

	void	operator+=(const Matrix& aM);
	void	operator-=(const Matrix& aM);
	void	operator*=(const double& aK);
	void	operator/=(const double& aK);


	friend bool	operator==(const Matrix& aM, const Matrix& aN);
	friend Matrix	operator-(const Matrix& aM);
	friend Matrix	operator+(const Matrix& aM, const Matrix& aN);
	friend Matrix	operator-(const Matrix& aM, const Matrix& aN);
	friend Matrix	operator*(const double& aK, const Matrix& aM);
	friend Matrix	operator*(const Matrix& aM, const double& aK);
	friend Matrix	operator/(const Matrix& aM, const double& aK);
	friend Matrix	operator*(const Matrix& aM, const Matrix& aN);
	// computes aM * aN
	// CAUTION the number of columns in aM must equal the number of rows in aN
	// NB This produces a (aM.GetNbRow(),aN.GetNbCol()) matrix
	friend Vector	operator*(const Matrix& aM, const Vector& aU);
	// computes aM * aU
	// CAUTION the size of aU must equal the number of columns in aM


	bool IsSquare() const { return mIsSquare; }
	bool IsSymmetric() const { return mIsSymmetric; }



	void ComputeEigenElts(Vector& EigenVals, Matrix& EigenVects) const;
	void ComputeSVD(Matrix& U, Matrix& W, Matrix& V) const;


protected:
	Matrix & Copy(const Matrix& aM);
	// Copy the given matrix 


	//
	// Members
	//
protected:
	vector<double>	mElt;	// coefficients of the matrix
	int				mNbRow;	// number of rows	
	int				mNbCol; // number of columns

	bool			mIsSquare;// true or false whether the matrix is a square matrix
	bool			mIsSymmetric;// true or false whether the (square) matrix is symmetric
	int				mSize;// size of the matrix in case it is square, 0 otherwise
};

Matrix	Transpose(const Matrix& aM);
// transpose of a (n,p) matrix (NB: produces a (p,n) matrix)

Matrix	abT(const Vector& aU, const Vector& aV);
// a b^T product, where a is of size n and b is of size p
// NB: this produces an (n,p) matrix

Matrix	aaT(const Vector& aU);
// u u^T product 
//(NB: produces a square symmetric matrix)

Matrix mmT(const Matrix& aM);
// M M^T product 
//(NB: produces a symmetric matrix)

Matrix mTm(const Matrix& aM);
// M^T M product 
//(NB: produces a symmetric matrix)

Matrix	Inverse_3(const Matrix& aM);
// Inverse of a 3 by 3 matrix by simple computation of cofactors and determinant

Matrix	Inverse_Sym(const Matrix& aM);
// Inverse of a real symmetric matrix by Cholesky decomposition
// CAUTION this algorithm only applies to symmetric positive-definite matrices.
// NB uses routines 'choldc' and 'cholsl' of Numerical Recipes in C, Press et al. 1992

Matrix	Inverse(const Matrix& aM);
// so far only for symmetric matrices or 3 by 3 matrices

double Determinant_3(const Matrix& aM);
// Determinant of a 3 by 3 matrix

double MinorDeterminant_3(const int RowTop, const int ColLeft, const int RowBottom, const int ColRight, Matrix& aM);
// Minor determinant of a 3 by 3 matrix (hence the determinant of a 2 by 2 matrix)

double Determinant_Sym(const Matrix& aM);
// Determinant of a real symmetric matrix by computing the product of its eigen-values

double Determinant(const Matrix& aM);
// so far only for symmetric matrices or 3 by 3 matrices


#endif// UNEWMATRIX_H
