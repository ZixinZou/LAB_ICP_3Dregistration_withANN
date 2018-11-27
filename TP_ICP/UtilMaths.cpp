//////////////////////////////////////////////////////////////////////
//  UtilMaths.cpp 
//
// Modified : Raoul de Charette --- 2/2010
//
//////////////////////////////////////////////////////////////////////
#include "stdafx.h"
//////////////////////////////////////////////////////////////////////
#include "UtilMaths.h"
//////////////////////////////////////////////////////////////////////

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

#define MAX(a,b) ((b) > (a) ? (b) : (a))

#define MIN(a,b) ((a) < (b) ? (a) : (b))

#define SQR(x)			((x) * (x))

#define PT_PT_DIST(Pt1, Pt2)    sqrt(SQR(Pt1[0] - Pt2[0]) + \
                                     SQR(Pt1[1] - Pt2[1]) + \
                                     SQR(Pt1[2] - Pt2[2]))


//////////////////////////////////////////////////////////////////////
//  This routine permutes the elements of the real symmetric matrix 'Mat'
//  so that all the larger elements in the matrix are relocated towards
//  the lower-right portion of the matrix.
//  This is done to minimize round-off errors in the tred2 and tqli algorithms,
//	for which the scatter matrix is being computed as an input.
//	"Order" is a 3 element array that will keep track of the order
//	changes in the matrix. The ordered matrix is returned in 'Mat'.
void Permute_Matrix(double** Mat, int* Order, int n)
{
	int i,j,k,ord=0;
	double a=0.0;
	double*		diag	= new double[n];
	double*		col		= new double[n];
	double*		row		= new double[n];

	for(i=0;i<n;i++)
		{
		diag[i]=col[i]=row[i]=0.0;
		Order[i]=i;
		}

	for(k=0;k<n;k++) diag[k] = Mat[k][k];

 	for (j=2;j<=n;j++) {
		i=j-1;
		a=diag[j-1];
			while (i > 0 && diag[i-1] > a) {
			// permute elements diag[i] and diag[i-1],
			// as well as columns Mat[.][i] and Mat[.][i-1]
			// and rows Mat[i][.] and Mat[i-1][.]
			a = diag[i];
			ord = Order[i];
			for(k=0;k<n;k++) col[k] = Mat[k][i]; // current column of index i

			diag[i] = diag[i-1];
			diag[i-1] = a;
			Order[i] = Order[i-1];
			Order[i-1] = ord;
			for(k=0;k<n;k++) // permute columns
				{
				Mat[k][i] = Mat[k][i-1];
				Mat[k][i-1] = col[k];
				}
			for(k=0;k<n;k++) row[k] = Mat[i][k]; // current row of index i
			for(k=0;k<n;k++)
				{
				Mat[i][k] = Mat[i-1][k];
				Mat[i-1][k] = row[k];
				}
			i--;
			}
		}

	delete[] diag;
	delete[] row;
	delete[] col;

}
//////////////////////////////////////////////////////////////////////



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
void tred2(double** a,double* d, double* e,int n)
{
  int		l,k,i,j;
  double	scale,hh,h,g,f;

	for(i=n;i>=2;i--)
	{
	l=i-1;
	h=scale=0.0;
	if(l>1)
		{
		for(k=1;k<=l;k++)
			scale+=fabs(a[i-1][k-1]);
		if(scale==0.0)		// skip transformation
			e[i-1]=a[i-1][l-1];
		else
			{
			for(k=1;k<=l;k++)
				{
				a[i-1][k-1]/=scale;	// use scaled a's for transformation.
				h+=a[i-1][k-1]*a[i-1][k-1];	// form sigma in h.
				}
			f=a[i-1][l-1];
			g=f>0? -sqrt(h):sqrt(h);
			e[i-1]=scale*g;
			h-=f*g;	// now h is equation (11.2.4)
			a[i-1][l-1]=f-g;	// store u in the ith row of a.
			f=0.0;
			for(j=1;j<=l;j++)
				{
				a[j-1][i-1]=a[i-1][j-1]/h; // store u/H in ith column of a.
				g=0.0;	// form an element of A.u in g
				for(k=1;k<=j;k++)
					g+=a[j-1][k-1]*a[i-1][k-1];
				for(k=j+1;k<=l;k++)
					g+=a[k-1][j-1]*a[i-1][k-1];
				e[j-1]=g/h; // form element of p in temorarliy unused element of e.
				f+=e[j-1]*a[i-1][j-1];
				}
			hh=f/(h+h);	// form K, equation (11.2.11)
			for(j=1;j<=l;j++) // form q and store in e overwriting p.
				{
				f=a[i-1][j-1]; // Note that e[l]=e[i-1] survives
				e[j-1]=g=e[j-1]-hh*f;
				for(k=1;k<=j;k++) // reduce a, equation (11.2.13)
					a[j-1][k-1]-=(f*e[k-1]+g*a[i-1][k-1]);
				}
			}
		}
	else
		e[i-1]=a[i-1][l-1];
	d[i-1]=h;
	}

  //	For computing eigenvector.
  d[0]=0.0;
  e[0]=0.0;

  for(i=1;i<=n;i++) // begin accumulating of transfomation matrices
	{
	l=i-1;
	if(d[i-1]) // this block skipped when i=1
		{
		for(j=1;j<=l;j++)
			{
			g=0.0;
			for(k=1;k<=l;k++) // use u and u/H stored in a to form P.Q
				g+=a[i-1][k-1]*a[k-1][j-1];
			for(k=1;k<=l;k++)
				a[k-1][j-1]-=g*a[k-1][i-1];
			}
		}
	d[i-1]=a[i-1][i-1];
	a[i-1][i-1]=1.0; // reset row and column of a to identity matrix for next iteration
	for(j=1;j<=l;j++)
		a[j-1][i-1]=a[i-1][j-1]=0.0;
	}
}
//////////////////////////////////////////////////////////////////////

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
void tqli(double* d,double* e,double** z,int n)
{
  int		m,l,iter,i,k;
  double	s,r,p,g,f,dd,c,b;

  for(i=2;i<=n;i++)
	e[i-2]=e[i-1];	// convenient to renumber the elements of e
  e[n-1]=0.0;
  for(l=1;l<=n;l++)
	{
	iter=0;
	do
		{
		for(m=l;m<=n-1;m++)
			{
			//	Look for a single small subdiagonal element
			//	to split the matrix.
			dd=fabs(d[m-1])+fabs(d[m]);
			if(fabs(e[m-1])+dd == dd)
				break;
			}
		if(m!=l)
			{
			if(iter++ == 30)
				{
//				printf("\nToo many iterations in TQLI");
                                // TODO: show corresponding error message!!
				}
			g=(d[l]-d[l-1])/(2.0*e[l-1]); /* form shift */
			r=sqrt((g*g)+1.0);
			g=d[m-1]-d[l-1]+e[l-1]/(g+SIGN2(r,g)); /* this is dm-ks */
			s=c=1.0;
			p=0.0;
			for(i=m-1;i>=l;i--)
				{
				//	A plane rotation as in the original
				//	QL, followed by Givens rotations to
				//	restore tridiagonal form.
				f=s*e[i-1];
				b=c*e[i-1];
				if(fabs(f) >= fabs(g))
					{
					c=g/f;
					r=sqrt((c*c)+1.0);
					e[i]=f*r;
					c*=(s=1.0/r);
					}
				else
					{
					s=f/g;
					r=sqrt((s*s)+1.0);
					e[i]=g*r;
					s*=(c=1.0/r);
					}
				g=d[i]-p;
				r=(d[i-1]-g)*s+2.0*c*b;
				p=s*r;
				d[i]=g+p;
				g=c*r-b;
				for(k=1;k<=n;k++)
					{
					//	Form eigenvectors
					f=z[k-1][i];
					z[k-1][i]=s*z[k-1][i-1]+c*f;
					z[k-1][i-1]=c*z[k-1][i-1]-s*f;
					}
				}
			d[l-1]=d[l-1]-p;
			e[l-1]=g;
			e[m-1]=0.0;
			}
		}while(m != l);
	}
}
//////////////////////////////////////////////////////////////////////

// Resolution du systeme : Ax=b  par la méthode de Gauss
void RegularGaussj(const int n, double* b, double** A)
{
	int *indxc,*indxr,*ipiv;
	int i,icol,irow,j,k,l,ll;
	double big,dum,pivinv;

	indxc = new int[n];
	indxr = new int[n];
	ipiv  = new int[n];

	for (j=1;j<=n;j++) ipiv[j-1]=0;
	for (i=1;i<=n;i++)
        {
		big=0.0;
		for (j=1;j<=n;j++)
			if (ipiv[j-1] != 1)
				for (k=1;k<=n;k++)
                                {
					if (ipiv[k-1] == 0)
                                        {
						if (fabs(A[j-1][k-1]) >= big)
                                                {
							big=(double)fabs(A[j-1][k-1]);
							irow=j;
							icol=k;
						}
					}
                                        else if (ipiv[k-1] > 1)
                                        {
                                        // printf("\n GAUSSJ: Singular Matrix-1");
                                        //TODO: show error message!!
                                        }
				}
		++(ipiv[icol-1]);
		if (irow != icol) {
			for (l=1;l<=n;l++) SWAP(double,A[irow-1][l-1],A[icol-1][l-1]);
			SWAP(double,b[irow-1],b[icol-1]);
		}
		indxr[i-1]=irow;
		indxc[i-1]=icol;
		if (A[icol-1][icol-1] == 0.0)
                        //printf("\n GAUSSJ: Singular Matrix-2");
                        //TODO:show error message!!
		pivinv=(double)(1.0/A[icol-1][icol-1]);
		A[icol-1][icol-1]=1.0;
		for (l=1;l<=n;l++) A[icol-1][l-1] *= pivinv;
		b[icol-1] *= pivinv;
		for (ll=1;ll<=n;ll++)
			if (ll != icol) {
				dum=A[ll-1][icol-1];
				A[ll-1][icol-1]=0.0;
				for (l=1;l<=n;l++) A[ll-1][l-1] -= A[icol-1][l-1]*dum;
				b[ll-1] -= b[icol-1]*dum;
			}
	}
	for (l=n;l>=1;l--) {
		if (indxr[l-1] != indxc[l-1])
			for (k=1;k<=n;k++)
				SWAP(double,A[k-1][indxr[l-1]-1],A[k-1][indxc[l-1]-1]);
	}

	delete	ipiv;
	delete	indxr;
	delete	indxc;
}
//////////////////////////////////////////////////////////////////////

// Given a positive-definite symmetric matrix a[1..n][1..n], this routine constructs its Cholesky
// decomposition, A= L L^T. On input, only the upper triangle of a need be given; it is not
// modified. The Cholesky factor L is returned in the lower triangle of a, except for its diagonal
// elements, which are returned in p[1..n]
// (Numerical Recipes in C, 2nd edition, Press et al. Cambridge  University Press, 1992)

// Complexity: n^3/6

// Since in C++ they should vary from 0 to n-1, all the indices [k] have been changed by [k-1]  (TC, 05/02/2001)

int choldc(double **a,int n,double *p)
{
  int             i, j, k;
  double           sum;

  for (i = 1; i <= n; i++)
  {
    for (j = i; j <= n; j++)
	{
      for (sum = a[i-1][j-1], k = i - 1; k >= 1; k--)
        sum -= a[i-1][k-1] * a[j-1][k-1];
      if (i == j)
	  {
        if (sum <= 0.0)
		{
          return 0; // Not OK
        }
        p[i-1] = sqrt(sum);
      } else
        a[j-1][i-1] = sum / p[i-1];
    }
  }
  return 1; // OK
}
//////////////////////////////////////////////////////////////////////




void cholsl(double **a,int n,double *p,double *b,double *x)
// Solves the set of linear equations a.x=b, where a is a positive-definite symmetric matrix
// a[1..n][1..n] and p[1..n] are input as the output of the routine choldc. Only the lower
// triangle of a is accessed. b[1..n] is input as the right-hand side vector. The solution vector is
// returned in x[1..n]. a, n, and p are not modified and can be left in place for successive calls
// with different right-hand side b. b is not modified unless you identify b and x in the calling
// sequence, which is allowed.
// (Numerical Recipes in C, 2nd edition, Press et al. Cambridge  University Press, 1992)

// Since in C++ they should vary from 0 to n-1, all the indices [k] have been changed by [k-1]  (TC, 05/02/2001)
{
  int i, k;
  double sum;

  for (i = 1; i <= n; i++) // Solve L.y = b, storing y in x
  {
    for (sum = b[i-1], k = i - 1; k >= 1; k--)
      sum -= a[i-1][k-1] * x[k-1];
    x[i-1] = sum / p[i-1];
  }
  for (i = n; i >= 1; i--) // Solve L^T.x = y
  {
    for (sum = x[i-1], k = i + 1; k <= n; k++)
      sum -= a[k-1][i-1] * x[k-1];
    x[i-1] = sum / p[i-1];
  }
}


////////////////////////////////////////////////////////////////////////////

// The following function computes the eigen elements of a 3x3
// real symmetric matrix. The results are ordered as follows:
// EigenValues  = [lambda_max, lambda_mid, lambda_min]
// and EigenVectors = [v_max, v_mid, v_min].
void SymMatr_EigenAnalysis(double** Matrix, double* EigenVals, double** EigenVects,int n)
{
	int*	Order	= new int[n];
	int*    Order2	= new int[n];
	double*	Diagonal_Matrix		= new double[n];
	double*	Diagonal_Matrix2	= new double[n];
	double*	Off_Diagonal_Matrix = new double[n];
	double	max, min;
	int	i,j;

	// First preprocess the matrix to put all the larger elements in the
	// bottom  right corner (in order to avoid round-off errors in tred2 and tqli).
	// Order keeps track of the order changes.
	//
	Permute_Matrix(Matrix, Order,n);

	tred2(Matrix,Diagonal_Matrix,Off_Diagonal_Matrix,n);
	tqli(Diagonal_Matrix,Off_Diagonal_Matrix,Matrix,n);

	//	Find the right order.
	min = Diagonal_Matrix[0];
	max = Diagonal_Matrix[0];


	//  Trouver la valeur propre minimum
	for (i=0; i<n; i++)
	{
		Order2[i] = i;
		if (Diagonal_Matrix[i]<min)		min = Diagonal_Matrix[i];
		Diagonal_Matrix2[i] = Diagonal_Matrix[i];
	}

	// On classe suivant les valeurs propres décroissantes
	min = min - 1.0;
	for (i=0;i<n;i++) {
		max = min;
		for (j=0;j<n;j++) {
			if (Diagonal_Matrix2[j]>max) {
				max			= Diagonal_Matrix2[j];
				Order2[i]	= j;
			}
		}
		Diagonal_Matrix2[Order2[i]] = min;
	}

	// EigenVals et EigenVects sont classés suivant l'ordre décroissant des valeurs propres
	for (i=0;i<n;i++) {
		EigenVals[i] = Diagonal_Matrix[Order2[i]];
		for (j=0;j<n;j++) {
			EigenVects[Order[i]][j] = Matrix[i][Order2[j]];
		}
	}

	delete Order;
	delete Order2;
	delete Diagonal_Matrix;
	delete Diagonal_Matrix2;
	delete Off_Diagonal_Matrix;
}
//////////////////////////////////////////////////////////////////////

//Given a matrix a[1..m][1..n], this routine computes its singular value decomposition, A =
//U·W·V^t. The matrix U replaces a on output. The diagonal matrix of singular values W is output
//as a vector w[1..n]. The matrix V (not the transpose V T ) is output as v[1..n][1..n].
void svdcmp(float **a, int m, int n, float w[], float **v)
{
	bool flag;
	int i,its,j,jj,k,l,nm;
	float anorm,c,f,g,h,s,scale,x,y,z;

	float *rv1 = new float[n];
	g=scale=anorm=0.0;
	for (i=0;i<n;i++) {
		l=i+2;
		rv1[i]=scale*g;
		g=s=scale=0.0;
		if (i < m) {
			for (k=i;k<m;k++) scale += fabs(a[k][i]);
			if (scale != 0.0) {
				for (k=i;k<m;k++) {
					a[k][i] /= scale;
					s += a[k][i]*a[k][i];
				}
				f=a[i][i];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a[i][i]=f-g;
				for (j=l-1;j<n;j++) {
					for (s=0.0,k=i;k<m;k++) s += a[k][i]*a[k][j];
					f=s/h;
					for (k=i;k<m;k++) a[k][j] += f*a[k][i];
				}
				for (k=i;k<m;k++) a[k][i] *= scale;
			}
		}
		w[i]=scale *g;
		g=s=scale=0.0;
		if (i+1 <= m && i+1 != n) {
			for (k=l-1;k<n;k++) scale += fabs(a[i][k]);
			if (scale != 0.0) {
				for (k=l-1;k<n;k++) {
					a[i][k] /= scale;
					s += a[i][k]*a[i][k];
				}
				f=a[i][l-1];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a[i][l-1]=f-g;
				for (k=l-1;k<n;k++) rv1[k]=a[i][k]/h;
				for (j=l-1;j<m;j++) {
					for (s=0.0,k=l-1;k<n;k++) s += a[j][k]*a[i][k];
					for (k=l-1;k<n;k++) a[j][k] += s*rv1[k];
				}
				for (k=l-1;k<n;k++) a[i][k] *= scale;
			}
		}
		anorm=MAX(anorm,(fabs(w[i])+fabs(rv1[i])));
	}
	for (i=n-1;i>=0;i--) {
		if (i < n-1) {
			if (g != 0.0) {
				for (j=l;j<n;j++)
					v[j][i]=(a[i][j]/a[i][l])/g;
				for (j=l;j<n;j++) {
					for (s=0.0,k=l;k<n;k++) s += a[i][k]*v[k][j];
					for (k=l;k<n;k++) v[k][j] += s*v[k][i];
				}
			}
			for (j=l;j<n;j++) v[i][j]=v[j][i]=0.0;
		}
		v[i][i]=1.0;
		g=rv1[i];
		l=i;
	}
	for (i=MIN(m,n)-1;i>=0;i--) {
		l=i+1;
		g=w[i];
		for (j=l;j<n;j++) a[i][j]=0.0;
		if (g != 0.0) {
			g=(float)(1.0/g);
			for (j=l;j<n;j++) {
				for (s=0.0,k=l;k<m;k++) s += a[k][i]*a[k][j];
				f=(s/a[i][i])*g;
				for (k=i;k<m;k++) a[k][j] += f*a[k][i];
			}
			for (j=i;j<m;j++) a[j][i] *= g;
		} else for (j=i;j<m;j++) a[j][i]=0.0;
		++a[i][i];
	}
	for (k=n-1;k>=0;k--) {
		for (its=0;its<30;its++) {
			flag=true;
			for (l=k;l>=0;l--) {
				nm=l-1;
				if (fabs(rv1[l])+anorm == anorm) {
					flag=false;
					break;
				}
				if (fabs(w[nm])+anorm == anorm) break;
			}
			if (flag) {
				c=0.0;
				s=1.0;
				for (i=l;i<k+1;i++) {
					f=s*rv1[i];
					rv1[i]=c*rv1[i];
					if (fabs(f)+anorm == anorm) break;
					g=w[i];
					h=pythag(f,g);
					w[i]=h;
					h=(float)(1.0/h);
					c=g*h;
					s = -f*h;
					for (j=0;j<m;j++) {
						y=a[j][nm];
						z=a[j][i];
						a[j][nm]=y*c+z*s;
						a[j][i]=z*c-y*s;
					}
				}
			}
			z=w[k];
			if (l == k) {
				if (z < 0.0) {
					w[k] = -z;
					for (j=0;j<n;j++) v[j][k] = -v[j][k];
				}
				break;
			}
			if (its == 29) printf("no convergence in 30 svdcmp iterations");
			x=w[l];
			nm=k-1;
			y=w[nm];
			g=rv1[nm];
			h=rv1[k];
			f=(float)(((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y));
			g=pythag(f,1.0);
			f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
			c=s=1.0;
			for (j=l;j<=nm;j++) {
				i=j+1;
				g=rv1[i];
				y=w[i];
				h=s*g;
				g=c*g;
				z=pythag(f,h);
				rv1[j]=z;
				c=f/z;
				s=h/z;
				f=x*c+g*s;
				g=g*c-x*s;
				h=y*s;
				y *= c;
				for (jj=0;jj<n;jj++) {
					x=v[jj][j];
					z=v[jj][i];
					v[jj][j]=x*c+z*s;
					v[jj][i]=z*c-x*s;
				}
				z=pythag(f,h);
				w[j]=z;
				if (z) {
					z=(float)(1.0/z);
					c=f*z;
					s=h*z;
				}
				f=c*g+s*y;
				x=c*y-s*g;
				for (jj=0;jj<m;jj++) {
					y=a[jj][j];
					z=a[jj][i];
					a[jj][j]=y*c+z*s;
					a[jj][i]=z*c-y*s;
				}
			}
			rv1[l]=0.0;
			rv1[k]=f;
			w[k]=x;
		}
	}

	delete []rv1;
}



//Computes (a2 + b2)1/2 without destructive underflow or overflow.
float pythag(float a, float b)
{
	float absa,absb;
	absa=fabs(a);
	absb=fabs(b);
	if (absa > absb) return absa*sqrt(1.0f+SQR(absb/absa));
	else return (absb == 0.0f ? 0.0f : absb*sqrt(1.0f+SQR(absa/absb)));
}
