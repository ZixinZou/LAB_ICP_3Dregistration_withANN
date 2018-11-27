 //
 //  Vector.h -- Vector class
 //
 //
 // Author : thomas Chaperon --- 04/2001
 //
 // Modified : xavier Brun --- 11/2004
 // Modified : Raoul de Charette --- 2/2010
//////////////////////////////////////////////////////////////////////
#ifndef VECTOR_H
#define VECTOR_H



#include<vector>

class Vector
{

public:
	// constructors, destructors
	//
	Vector();
	// the following  method allocates and initializes an n vector
	// whose coefficients are all set to zero.
	Vector(const int& n);
	Vector(const Vector& aV); // copy constructor
	~Vector();


	void	Init(const int& n);
	void	SetToZero();

	
	int     Size() const             { return (int)Elt.size(); }
	void    Resize(const int& aN)   { Elt.resize(aN);}
	void	Normalize();
	
	double&         operator[] (const int aIndex) { return Elt[aIndex]; }
	const double&   operator[] (const int aIndex) const { return Elt[aIndex]; }

	void	operator=(const Vector& aV);
	void	operator=(const float *aV);
	void	operator+=(const Vector& aV);
	void	operator+=(const float *aV);
	void	operator-=(const Vector& aV);
	void	operator*=(const double& aK);
                // multiplication by a scalar
	void	operator/=(const double& aK);
                // division by scalar

	 


	friend bool     operator==(const Vector& aM, const Vector& aN);
                //tests whether the two vectors are equal
                // nb: to be equal they must have the same size!

	friend Vector	operator-(const Vector& aV);
                //opposite

	friend Vector	operator+(const Vector& aM, const Vector& aN);
	friend Vector	operator-(const Vector& aM, const Vector& aV);
	friend Vector	operator*(const double& aK, const Vector& aV);
                // multiplication by a scalar
	friend Vector	operator*(const Vector& aV, const double& aK);
                // multiplication by a scalar
	friend Vector	operator/(const Vector& aV, const double& aK);
                // division by a scalar

	friend double	operator*(const Vector& aM, const Vector& aN);
                // Euclidean scalar product of R^n

	friend Vector operator^(const Vector& aM, const Vector& aN);
				//produit vectoriel

	//
	// Members
	//
public:
	std::vector<double> Elt;
};



#endif
