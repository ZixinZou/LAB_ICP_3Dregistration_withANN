 //
 //  Matrix.cpp --	Matrix class
 //
 //
 // Author : thomas Chaperon --- 04/2001
 //
 // Modified : xavier Brun --- 11/2004
 // Morified : Taha Ridene ---   2005-2007
 // Modified : Raoul de Charette --- 2/2010
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "Matrix.h"

/////////////////////////////////////////////////////////////////////////////
// Constructors, destructors

Matrix::Matrix()
{
}
  
	
Matrix::Matrix(const int& aNbRows, const int& aNbCols, const bool& aIsSymmetric)
{
	mNbRow = aNbRows;
	mNbCol = aNbCols;
	mElt.resize(aNbRows*aNbCols);

	if(mNbCol== mNbRow)
	{
		mIsSquare	= true;
		mSize		= mNbCol;
	}
	else
	{
		mIsSquare	= false;
		mSize		= 0;
	}
	mIsSymmetric = aIsSymmetric;

	// Set every element in the matrix to 0
	SetToZero();
}
	
Matrix::Matrix(const int& aSize, const bool& aIsSymmetric)
{
	mNbRow	= aSize;
	mNbCol	= aSize;
	mElt.resize(aSize*aSize);
	mSize	= aSize;
	mIsSquare	= true;
	mIsSymmetric = aIsSymmetric;

	// Set every element in the matrix to 0
	SetToZero();
}


Matrix::Matrix(const vector<Vector>& aColVectors)
{
	int i;

	// Check whether the list is empty or not
	ASSERT(aColVectors.size());

	// Check whether the vectors all have the same size
	for(i=1;i<(int)aColVectors.size();i++)
		ASSERT(aColVectors[i].Size() == aColVectors[0].Size());


	mNbRow = aColVectors[0].Size();
	mNbCol = (int)aColVectors.size();
	mElt.resize(mNbRow*mNbCol);
	if(mNbRow==mNbCol)
	{
		mSize = mNbRow;
		mIsSquare = true;
	}
	else
	{
		mSize = 0;
		mIsSquare = false;
	}
	mIsSymmetric = false;

	// Set the columns of the matrix
	for(i=0;i<mNbCol;i++)
		SetCol(i, aColVectors[i]);
}

Matrix::Matrix(const Matrix& aM) 
{ 
	Copy(aM);
}


Matrix::~Matrix()
{
	mElt.clear();
}

/////////////////////////////////////////////////////////////////////////////
 
void Matrix::SetToZero()
{
	// Set every element in the matrix to 0
	int i,j;
	for(i=0;i<mNbRow;i++)
		for(j=0;j<mNbCol;j++)
			SetElt(i,j, 0.);
}
/////////////////////////////////////////////////////////////////////////////

 	
void Matrix::Resize(const int& aNbRows, const int& aNbCols)
{	
	mElt.clear();

	mNbRow = aNbRows;
	mNbCol = aNbCols;
	mElt.resize(aNbRows*aNbCols);

	// Set every element in the matrix to 0
	SetToZero();
}
/////////////////////////////////////////////////////////////////////////////


Matrix& Matrix::Copy(const Matrix& aM)
{
	mNbRow = aM.GetNbRow();
	mNbCol = aM.GetNbCol();
	if(mNbCol== mNbRow)
	{
		mIsSquare	= true;
		mSize		= mNbCol;
	}
	else
	{
		mIsSquare	= false;
		mSize		= 0;
	}
	mIsSymmetric = aM.IsSymmetric();

	mElt.clear();
	mElt.resize(aM.GetNbCol() * aM.GetNbRow());
	int i,j;
	for(i=0;i<GetNbRow();i++)
	{
		for(j=0;j<GetNbCol();j++)
			SetElt(i,j, aM.GetElt(i,j));
	}
	return *this;
}

/////////////////////////////////////////////////////////////////////////////
void Matrix::operator=(const Matrix& aM)
{
	Copy(aM);
}

/////////////////////////////////////////////////////////////////////////////
 
void Matrix::operator +=(const Matrix& aM)
{
	mIsSymmetric = (mIsSymmetric && aM.IsSymmetric());

	int i,j;
	for(i=0;i<mNbRow;i++)
		for(j=0;j<mNbCol;j++) 			
			SetElt(i,j, GetElt(i,j) + aM.GetElt(i,j));
}
/////////////////////////////////////////////////////////////////////////////
 
void Matrix::operator -=(const Matrix& aM)
{
	mIsSymmetric = (mIsSymmetric && aM.IsSymmetric());

	int i,j;
	for(i=0;i<mNbRow;i++)
		for(j=0;j<mNbCol;j++) 			
			SetElt(i,j, GetElt(i,j) - aM.GetElt(i,j));
}
/////////////////////////////////////////////////////////////////////////////
 
void Matrix::operator*=(const double& aK)
{
	int i,j;
	for(i=0;i<mNbRow;i++)
		for(j=0;j<mNbCol;j++) 			
			SetElt(i,j, aK * GetElt(i,j));
}
/////////////////////////////////////////////////////////////////////////////
 
void Matrix::operator/=(const double& aK)
{
	ASSERT(aK);
	int i,j;
	for(i=0;i<mNbRow;i++)
		for(j=0;j<mNbCol;j++) 			
			SetElt(i,j, GetElt(i,j) / aK);
}
/////////////////////////////////////////////////////////////////////////////

bool operator==(const Matrix& aM, const Matrix& aN)
{
	ASSERT((aM.GetNbRow() == aN.GetNbRow()) && (aM.GetNbCol() == aN.GetNbCol()));
	int i,j;
	for(i=0;i<aM.GetNbRow();i++)
		for(j=0;j<aM.GetNbCol();j++) 			
			if(aM.GetElt(i,j) != aN.GetElt(i,j)) return false;
	return true;
}
/////////////////////////////////////////////////////////////////////////////

Matrix	operator-(const Matrix& aM)
{
	Matrix P(aM.GetNbRow(), aM.GetNbCol(), aM.IsSymmetric());

	int i,j;
	for(i=0;i<P.GetNbRow();i++)
		for(j=0;j<P.GetNbCol();j++) 			
			P.SetElt(i,j, -aM.GetElt(i,j));
	return P;
}
/////////////////////////////////////////////////////////////////////////////
	
 
Matrix	operator+(const Matrix& aM , const Matrix& aN)
{
	ASSERT((aM.GetNbRow() == aN.GetNbRow()) && (aM.GetNbCol() == aN.GetNbCol()));

	int i,j;
	Matrix P(aM.GetNbRow(), aM.GetNbCol(), (aM.IsSymmetric() && aN.IsSymmetric()) );
	for(i=0;i<P.GetNbRow();i++)
		for(j=0;j<P.GetNbCol();j++) 			
			P.SetElt(i,j, aM.GetElt(i,j) + aN.GetElt(i,j));
	return P;
}
/////////////////////////////////////////////////////////////////////////////

 
Matrix	operator-(const Matrix& aM , const Matrix& aN)
{
	ASSERT((aM.GetNbRow() == aN.GetNbRow()) && (aM.GetNbCol() == aN.GetNbCol()));
	int i,j;
	Matrix P(aM.GetNbRow(), aM.GetNbCol(), (aM.IsSymmetric() && aN.IsSymmetric()));
	for(i=0;i<P.GetNbRow();i++)
		for(j=0;j<P.GetNbCol();j++) 			
			P.SetElt(i,j, aM.GetElt(i,j) - aN.GetElt(i,j));
	return P;
}
/////////////////////////////////////////////////////////////////////////////
// multiplication by a scalar
 
Matrix	operator*(const double& aK, const Matrix& aM)
{
	int i,j;
	Matrix P(aM.GetNbRow(), aM.GetNbCol(), aM.IsSymmetric());
	for(i=0;i<P.GetNbRow();i++)
		for(j=0;j<P.GetNbCol();j++) 			
			P.SetElt(i,j, aK * aM.GetElt(i,j));
	return P;
}
/////////////////////////////////////////////////////////////////////////////
// division by a scalar
 
Matrix	operator/(const double& aK, const Matrix& aM)
{
	int i,j;
	Matrix P(aM.GetNbRow(), aM.GetNbCol(), aM.IsSymmetric());
	for(i=0;i<P.GetNbRow();i++)
		for(j=0;j<P.GetNbCol();j++) 			
			P.SetElt(i,j, aM.GetElt(i,j) / aK);
	return P;
}
/////////////////////////////////////////////////////////////////////////////
// matrix product

 
Matrix	operator*(const Matrix& aM, const Matrix& aN)
{
	ASSERT( (aM.GetNbCol() == aN.GetNbRow()) );

	int i,j,k;
	Matrix P(aM.GetNbRow(), aN.GetNbCol(), (aM.IsSymmetric() && aN.IsSymmetric()));
	double EltTmp;
	for(i=0;i<P.GetNbRow();i++)
	{
		for(j=0;j<P.GetNbCol();j++) 			
		{
			EltTmp = 0.;
			for(k=0;k<aM.GetNbCol();k++)	EltTmp += aM.GetElt(i,k) * aN.GetElt(k,j);		
			P.SetElt(i,j, EltTmp);
		}
	}
	return P;
}

/////////////////////////////////////////////////////////////////////////////

 
Vector	operator*(const Matrix& aM, const Vector& aU)
{
	ASSERT( (aU.Size() == aM.GetNbCol()) );

	int i,j;
	double EltTmp;
	Vector V( aM.GetNbRow() );	
	for(i=0;i<V.Size();i++)
	{
		EltTmp =0.;
		for(j=0;j<aM.GetNbCol();j++)
#if(1)
			EltTmp += aM.GetElt(i,j) * aU.Elt[j];
#endif
#if(0)
			EltTmp += aM.GetElt(i,j) * aU.GetElt(j);
#endif
#if(1)
		V.Elt[i] = EltTmp;
#endif
#if(0)
		V.SetElt(i, EltTmp);
#endif
	}
	return V;
}
/////////////////////////////////////////////////////////////////////////////


 	
Vector	Matrix::GetRow(const int& aRow) const
{	
	Vector V(mNbCol);
	int i;
	for(i=0;i<mNbCol;i++)
		V.Elt[i] = GetElt(aRow,i);
	return V;
}
/////////////////////////////////////////////////////////////////////////////
	
 	
Vector	Matrix::GetCol(const int& aCol) const
{	
	Vector V(mNbRow);
	int i;
	for(i=0;i<mNbRow;i++)
		V.Elt[i] = GetElt(i, aCol);
	return V;
}
/////////////////////////////////////////////////////////////////////////////
 	
void Matrix::SetRow(const int& aRow, const Vector& aV)	
{
	ASSERT( (aV.Size() == mNbCol) );
	int i;
	for(i=0;i<mNbCol;i++)
		SetElt(aRow,i, aV.Elt[i]);
}
/////////////////////////////////////////////////////////////////////////////
		
 	
void	Matrix::SetCol(const int& aCol, const Vector& aV)
{
	ASSERT( (aV.Size() == mNbRow) );
	int i;
	for(i=0;i<mNbRow;i++)
		SetElt(i,aCol, aV.Elt[i]);
}
/////////////////////////////////////////////////////////////////////////////

 
Matrix Transpose(const Matrix& aM)
{
	if(aM.IsSymmetric())
	{
		//Matrix* P = new Matrix(aM);
		//return *P;
		return aM;
	}
	
	Matrix P(aM.GetNbCol(), aM.GetNbRow());		
	int i,j;
	for(i=0;i<P.GetNbRow();i++)
	{	for(j=0;j<P.GetNbCol();j++) 			
			P.SetElt(i,j, aM.GetElt(j,i));
	}	
	return P;
}
/////////////////////////////////////////////////////////////////////////////
	
Matrix	abT(const Vector& aU, const Vector& aV)
{
	Matrix M(aU.Size(), aV.Size());
	int i,j;
	for(i=0;i<M.GetNbRow();i++)
	{
		for(j=0;j<M.GetNbCol();j++)
			M.SetElt(i,j, aU.Elt[i] * aV.Elt[j]);
	}
	return M;
}
/////////////////////////////////////////////////////////////////////////////

Matrix	aaT(const Vector& aU)
{
	Matrix M(aU.Size(), true); // define M as square and symmetric
	int i,j;
	for(i=0;i<M.GetNbRow();i++)
	{
		for(j=0;j<M.GetNbCol();j++)
			M.SetElt(i,j, aU.Elt[i] * aU.Elt[j]);
	}
	return M;
}
///////////////////////////////////////////////////////////////////


// this routine computes M M^T where M is a non-necessarily square matrix.
// It produces a square symmetric matrix of size the number of rows in M
Matrix mmT(const Matrix& aM)
{
	Matrix P(aM.GetNbRow(), true); // define P as square symmetric
	int i,j;
	for(i=0;i<P.GetSize();i++)
		for(j=0;j<P.GetSize();j++)
			P.SetElt(i,j, aM.GetRow(i) * aM.GetRow(j));// vector product
	return P;			
}
///////////////////////////////////////////////////////////////////

// this routine computes M^T M where M is a non-necessarily square matrix
// it produces a square symmetric matrix of size the number of columns in M
Matrix mTm(const Matrix& aM)
{
	Matrix P(aM.GetNbCol(), true); // define P as square symmetric
	int i,j;
	for(i=0;i<P.GetSize();i++)
		for(j=0;j<P.GetSize();j++)
			P.SetElt(i, j, aM.GetCol(i) * aM.GetCol(j)); // vector product
	return P;			
}
///////////////////////////////////////////////////////////////////

// Matrix inversion of a 3 by 3 matrix via the matrix of cofactors

Matrix Inverse_3(const Matrix& aM)
{
	ASSERT(aM.IsSquare() && (aM.GetSize() == 3));

	double Det = Determinant_3(aM);
	// Check whether the matrix is invertible
	ASSERT(Det != 0); 

	// Compute M^-1 = Det^-1 * Transpose(Com(M))

	Matrix P(3);
	P.SetElt(0,0,  (aM(1,1)*aM(2,2)-aM(2,1)*aM(1,2))/ Det);
	P.SetElt(0,1, -(aM(0,1)*aM(2,2)-aM(2,1)*aM(0,2))/ Det);
	P.SetElt(0,2,  (aM(0,1)*aM(1,2)-aM(0,2)*aM(1,1))/ Det);

	P.SetElt(1,0, -(aM(1,0)*aM(2,2)-aM(2,0)*aM(1,2))/ Det);
	P.SetElt(1,1,  (aM(0,0)*aM(2,2)-aM(2,0)*aM(0,2))/ Det);
	P.SetElt(1,2, -(aM(0,0)*aM(1,2)-aM(1,0)*aM(0,2))/ Det);

	P.SetElt(2,0,  (aM(1,0)*aM(2,1)-aM(2,0)*aM(1,1))/ Det);
	P.SetElt(2,1, -(aM(0,0)*aM(2,1)-aM(2,0)*aM(0,1))/ Det);
	P.SetElt(2,2,  (aM(0,0)*aM(1,1)-aM(1,0)*aM(0,1))/ Det);
		// NB actually Det^-1 * Com(Transpose(M)) is computed above.
		// This is equivalent since: det(A) = A * (Com(A))^T = (Com(A))^T A
		// thus Com(A^T) = (Com(A))^T

	return P;
}
///////////////////////////////////////////////////////////////////
// Matrix inversion by Cholesky decomposition
// CAUTION this algorithm only applies to symmetric positive-definite matrices.

Matrix	Inverse_Sym(const Matrix& aM)
{
	// Check whether the matrix is symmetric
	ASSERT( aM.IsSymmetric() );

	int  i,j;
	int n = aM.GetSize();

	Matrix P(n, true); // define P as symmetric

	// Copie avec conversion pour utiliser les fonctions choldc et cholsl de Numerical Recipes in C (dans Globals.cpp)
	double** Mtmp;
	Mtmp = new double*[n];
	for(i=0;i<n;i++)
		Mtmp[i] = new double[n];

	for(i=0;i<n;i++)
		for(j=0;j<n;j++)
			Mtmp[i][j] = aM.GetElt(i,j);// en fait seuls les éléments au dessus de la diagonale suffiraient ici (cf choldc)

	double* Dtmp= new double[n]; // diagonale



	// décomposition de la matrice sous la forme de Cholesky A= L L^T (L: triangulaire inférieure)
	// NB: lors de cette étape,
	// la partie triangulaire inferieure de Mtmp est remplacée par la matrice L
	// (sauf la diagonale qui est stockée dans Dtmp);
	// la partie triangulaire supérieure est laissée intacte.

	if (choldc(Mtmp, n, Dtmp)==0) 
	{  
		// matrice non inversible... Que faire??
		// TODO! treat this case!
	}
	else 
	{	
		double* ColonneTmp = new double[n];
		double* b_i = new double[n];

		for(i=0;i<n;i++)
		{
			for(j=0;j<n;j++) // vecteur colonne i de la matrice identit?de taille n
			{
				if(j==i)
					b_i[j] = 1;
				else 
					b_i[j] = 0;
			}

			// résolution du système
			cholsl(Mtmp, n, Dtmp, b_i, ColonneTmp);

			// stockage de la colonne dans la matrice P
			for(j=0;j<n;j++)
			{
				P.SetElt(j,i, ColonneTmp[j]);
			}
		}

		delete []ColonneTmp;
		ColonneTmp	= NULL;
		
		delete []b_i;
		b_i	= NULL;
	}


	for(i=0;i<n;i++)
	{
		delete Mtmp[i];
		Mtmp[i] = NULL;
	}
	delete []Mtmp;
	Mtmp	= NULL;


	delete []Dtmp;
	Dtmp = NULL;
	

	return P;
}
///////////////////////////////////////////////////////////////////

Matrix Inverse(const Matrix& aM)
{	
	ASSERT(aM.IsSquare());

	// so far only for symmetric matrices or 3 by 3 matrices
	ASSERT((aM.GetSize() == 3)||(aM.IsSymmetric())); 

	if(aM.GetSize() == 3)
		return Inverse_3(aM);
	else
		return Inverse_Sym(aM);
}

///////////////////////////////////////////////////////////////////

// The following function computes the eigen elements of a
// real symmetric matrix. The results are ordered as follows:
// EigenValues  = [lambda_max, ... , lambda_min]
// and EigenVectors = [v_max, ... , v_min].
void Matrix::ComputeEigenElts(Vector& EigenVals, Matrix& EigenVects) const
{
	// This works only for real symmetric matrices
	ASSERT( IsSymmetric() );

	int n;
	n = GetSize();

	int i,j;
	double** lMatrix;
	lMatrix = new double*[n];
	for(i=0;i<n;i++)
		lMatrix[i] = new double[n];

	for(i=0;i<n;i++)
		for(j=0;j<n;j++)
			lMatrix[i][j] = GetElt(i,j);


	int*	Order	= new int[n];
	int*    Order2	= new int[n];
	double*	Diagonal_Matrix		= new double[n];
	double*	Diagonal_Matrix2	= new double[n];
	double*	Off_Diagonal_Matrix = new double[n];
	double	max, min;

	// First preprocess the matrix to put all the larger elements in the
	// bottom  right corner (in order to avoid round-off errors in tred2 and tqli).
	// Order keeps track of the order changes.
	//
	Permute_Matrix(lMatrix, Order, n);

	tred2(lMatrix, Diagonal_Matrix, Off_Diagonal_Matrix, n);
	tqli(Diagonal_Matrix, Off_Diagonal_Matrix, lMatrix, n);

	//	Find the right order.
	min = Diagonal_Matrix[0];
	max = Diagonal_Matrix[0];


	//  Trouver la valeur propre minimum
	for (i=0; i<n; i++)
	{
		Order2[i] = i;
		if (Diagonal_Matrix[i] < min)		
			min = Diagonal_Matrix[i];
		Diagonal_Matrix2[i] = Diagonal_Matrix[i];
	}

	// On classe suivant les valeurs propres décroissantes
	min = min - 1.0;
	for (i=0;i<n;i++) 
	{
		max = min;
		for (j=0;j<n;j++) 
		{
			if (Diagonal_Matrix2[j] > max) 
			{
				max			= Diagonal_Matrix2[j];
				Order2[i]	= j;
			}
		}
		Diagonal_Matrix2[Order2[i]] = min;
	}

	// EigenVals et EigenVects sont classés suivant l'ordre décroissant des valeurs propres
	EigenVects.Resize(n, n);
	EigenVals.Resize(n);
	for (i=0;i<n;i++)
	{
		EigenVals.Elt[i] = Diagonal_Matrix[Order2[i]];
		for (j=0;j<n;j++) 
		{
			EigenVects.SetElt(Order[i], j, lMatrix[i][Order2[j]] );
		}
	}

	delete[] Order;
	delete[] Order2;
	delete[] Diagonal_Matrix;
	delete[] Diagonal_Matrix2;
	delete[] Off_Diagonal_Matrix;
	for(i=0;i<n;i++)
		delete[] lMatrix[i];
	delete[] lMatrix;
}

/////////////////////////////////////////////////////////////////////////////
//	The following function compute the Singular Value Decomposition 
//	de la matrice m*n : U*W*V^t
//  U : matrice de taille m*n orthogonal
//  W : matrice diagonale positive
//  V : matrice n*n orthogonal
void Matrix::ComputeSVD(Matrix& U,Matrix& W, Matrix& V) const
{
	int m,n;
	m = GetNbRow();
	n = GetNbCol();

	int i,j;
	float** lu;
	lu = new float*[m];
	for(i=0;i<m;i++)
		lu[i] = new float[n];

	for(i=0;i<m;i++)
		for(j=0;j<n;j++)
			lu[i][j] = (float)GetElt(i,j);

	float* lw=new float[n];

	float** lv;
	lv = new float*[n];
	for(i=0;i<n;i++)
		lv[i] = new float[n];

	svdcmp(lu,m,n,lw,lv);

	U.SetToZero();
	W.SetToZero();
	V.SetToZero();

	for (i=0;i<m;i++)
		for (j=0;j<n;j++) U.SetElt(i,j,lu[i][j]);

	for (i=0;i<n;i++)
	{
		W.SetElt(i,i,lw[i]);
		for (j=0;j<n;j++)
			V.SetElt(i,j,lv[i][j]);
	}

	for(i=0;i<m;i++)
		delete[] lu[i];
	delete[] lu;

	for(i=0;i<n;i++)
		delete[] lv[i];
	delete[] lv;

	delete[] lw;
}
/////////////////////////////////////////////////////////////////////////////

double Determinant_3(const Matrix& aM)
{
	// Check whether the matrix is 3 by 3
	ASSERT(aM.IsSquare() && (aM.GetSize() == 3));
	double Det;

	Det = aM(0,0)*aM(1,1)*aM(2,2) +aM(0,1)*aM(1,2)*aM(2,0) +aM(0,2)*aM(1,0)*aM(2,1)
			-aM(2,0)*aM(1,1)*aM(0,2) -aM(2,1)*aM(1,2)*aM(0,0) -aM(2,2)*aM(1,0)*aM(0,1);
	return Det;
}
/////////////////////////////////////////////////////////////////////////////

double MinorDeterminant_3(const int RowTop, const int ColLeft, const int RowBottom, const int ColRight, Matrix& aM)
{
	// Check whether the matrix is 3 by 3
	ASSERT(aM.IsSquare() && (aM.GetSize() == 3));

	return aM(RowTop,ColLeft)*aM(RowBottom,ColRight) - aM(RowBottom,ColLeft)*aM(RowTop,ColRight);
}

/////////////////////////////////////////////////////////////////////////////

double Determinant_Sym(const Matrix& aM)
{
	// Check whether the matrix is symmetric
	ASSERT(aM.IsSymmetric());

	double Det;

	// Compute the determinant as the product of the eigen-values
	//
	Matrix EigenVects(aM.GetSize());
	Vector EigenVals(aM.GetSize());
	aM.ComputeEigenElts(EigenVals, EigenVects);

	Det = 1.;
	for(int i=0;i<EigenVals.Size(); i++)
		Det *= EigenVals.Elt[i];
	return Det;
}
/////////////////////////////////////////////////////////////////////////////

double Determinant(const Matrix& aM)
{	
	// Check whether the matrix is square
	ASSERT(aM.IsSquare());

	// so far only for symmetric matrices or 3 by 3 matrices
	ASSERT((aM.GetSize() == 3)||(aM.IsSymmetric())); 

	if(aM.GetSize() == 3)
		return Determinant_3(aM);
	else
		return Determinant_Sym(aM);
}
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////


