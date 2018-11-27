 //
 //  Vector.cpp -- Vector class
 //
 //
 // Author : thomas Chaperon --- 04/2001
 //
 // Modified : xavier Brun --- 11/2004
 // Modified : Raoul de Charette --- 2/2010
//////////////////////////////////////////////////////////////////////
#include "stdafx.h"
#include <math.h>

#include "Vector.h"
///////////////////////////////////////////////////////////////////////////////

// Constructors
Vector::Vector()
{

}

Vector::Vector(const int& n)
{
	Init(n);
}

	
Vector::Vector(const Vector& aV)
{
	Elt = aV.Elt;
}
	
// Destructor
Vector::~Vector()
{	
	Elt.clear();
}
///////////////////////////////////////////////////////////////////////////////

void Vector::Init(const int& n)
{
	int i;
	Elt.resize(n);

	for(i=0;i<n;i++)
	{
		Elt[i] = 0.0;
	}
}

///////////////////////////////////////////////////////////////////////////////

void Vector::SetToZero()
{
	int i;

	for(i=0;i<Size();i++)
	{
		Elt[i] = 0.0;
	}
}

///////////////////////////////////////////////////////////////////////////////


void Vector::operator=(const Vector& aV)
{
	int i;
	Resize(aV.Size());
	for(i=0;i<Size();i++)
	{
		Elt[i] = aV.Elt[i];
	}
}

///////////////////////////////////////////////////////////////////////////////


void Vector::operator=(const float *aV)
{
	int i;
	for(i=0;i<Size();i++)
	{
		Elt[i] = aV[i];
	}
}


///////////////////////////////////////////////////////////////////////////////


void Vector::operator+=(const Vector& aV)
{
	int i;
	for(i=0;i<Size();i++)
	{
		Elt[i] += aV.Elt[i];
	}
}

///////////////////////////////////////////////////////////////////////////////


void Vector::operator+=(const float *aV)
{
	int i;
	for(i=0;i<Size();i++)
	{
		Elt[i] += aV[i];
	}
}


///////////////////////////////////////////////////////////////////////////////

void Vector::operator-=(const Vector& aV)
{
	int i;
	for(i=0;i<Size();i++)
	{
		Elt[i] -= aV.Elt[i];
	}
}

///////////////////////////////////////////////////////////////////////////////


void Vector::operator*=(const double& aK)
{
	int i;
	for(i=0;i<Size();i++)
	{
		Elt[i] *= aK;
	}
}

///////////////////////////////////////////////////////////////////////////////


void Vector::operator/=(const double& aK)
{
//	ASSERT(aK);
	int i;
	for(i=0;i<Size();i++)
	{
		Elt[i] /= aK;
	}
}
///////////////////////////////////////////////////////////////////////////////

void Vector::Normalize()
{
	int i;
	float norm = 0.0;
	for(i=0;i<Size();i++)
	{
		norm += (float)(Elt[i]*Elt[i]);
	}
	for(i=0;i<Size();i++)
	{
		Elt[i] /= sqrt(norm);
	}

}


///////////////////////////////////////////////////////////////////////////////

bool	operator==(const Vector& aM, const Vector& aN)
{
//	ASSERT((aM.Size() == aN.Size())); //otherwise the question makes no sense 
	int i;
	for(i=0;i<aM.Size();i++)
	{
		if(aM.Elt[i] != aN.Elt[i]) return false ;
	}
	return true;
}

///////////////////////////////////////////////////////////////////////////////


Vector	operator-(const Vector& aM)
{
	Vector P(aM.Size());
	int i;
	for(i=0;i<aM.Size();i++)
	{
		P.Elt[i] = aM.Elt[i];
	}
	return P;
}

///////////////////////////////////////////////////////////////////////////////


Vector	operator+(const Vector& aM , const Vector& aN)
{
	Vector P(aM.Size());
	int i;
	for(i=0;i<aM.Size();i++)
	{
		P.Elt[i] = aM.Elt[i] + aN.Elt[i];
	}
	return P;
}

///////////////////////////////////////////////////////////////////////////////

Vector	operator-(const Vector& aM , const Vector& aN)
{
	Vector P(aM.Size());
	int i;
	for(i=0;i<aM.Size();i++)
	{
		P.Elt[i] = aM.Elt[i] - aN.Elt[i];
	}
	return P;
}

///////////////////////////////////////////////////////////////////////////////


Vector	operator*(const Vector& aM, const double& aK)
{
	Vector P(aM.Size());
	int i;
	for(i=0;i<aM.Size();i++)
	{
		P.Elt[i] = aK * aM.Elt[i];
	}
	return P;
}

///////////////////////////////////////////////////////////////////////////////


Vector	operator/(const Vector& aM, const double& aK)
{
//	ASSERT(aK);
	Vector P(aM.Size());
	int i;
	for(i=0;i<aM.Size();i++)
	{
		P.Elt[i] =  aM.Elt[i] / aK;
	}
	return P;
}

///////////////////////////////////////////////////////////////////////////////


double	operator*(const Vector& aM, const Vector& aN)
{
	double p = 0.0;

	int i;
	for(i=0;i<aM.Size();i++)
	{
		p += aM.Elt[i] * aN.Elt[i];
	}
	return p;
}


///////////////////////////////////////////////////////////////////////////////

Vector operator^(const Vector& aM, const Vector& aN)
{
	Vector P(aM.Size());
	if ( (aM.Size() == 3) && (aN.Size() == 3) )
	{

		P.Elt[0] =  aM.Elt[1]*aN.Elt[2] - aM.Elt[2]*aN.Elt[1];
		P.Elt[1] =  aM.Elt[2]*aN.Elt[0] - aM.Elt[0]*aN.Elt[2];
		P.Elt[2] =  aM.Elt[0]*aN.Elt[1] - aM.Elt[1]*aN.Elt[0];
		
	}
	
	return P;	
}

///////////////////////////////////////////////////////////////////////////////


