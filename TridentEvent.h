//This is a class for an event for the Trident reaction
//It will be used to generate an MC generator and is based on the four-vector class
//Jonathan Curtis
//02/12/2016

#ifndef EVENT_H
#define EVENT_H

#include <cmath>
#include <iostream>

/////
//3-vector class 
/////
class Vector{
protected:
	double x;	//the x component of the vector
	double y; 	//the y component
	double z;	//the z component 

public:
	Vector(double, double, double);	//This creates a three-vector given components as x,y,z
	void Print();	//Prints the vector to cout

	//Operators
	double operator *(const Vector&);	//The normal dot-product	
	Vector operator +(const Vector&);	//This adds two three vectors component-wise
	friend Vector operator*(const double,const Vector&);	//This is for scalar multiplication from the left 
	Vector operator -(const Vector& rh);	//This subtracts two three vectors component wise 
};


//////
//4-vector
//////
class FourVector{
protected:
	double t;	//the time component (x0)
	Vector v;	//The spatial components (x1,x2,x3) 

public:
	FourVector(double,Vector);	//This creates a four-vector given time component and spatial vector
	void Print();	//Prints the vector to cout

	double operator *(const FourVector&);	//The Lorentz invaraint contraction of two fourvectors
	FourVector operator +(const FourVector&);	//This adds two fourvectors component-wise
	friend FourVector operator*(const double,const FourVector&);	//This is for scalar multiplication from the left 
	FourVector operator -(const FourVector& rh);	//This subtracts two fourvectors component wise 
};
/*

//////
//Trident events
//////
class TridentEvent{
public:
	static FourVector p1;	//Incoming neutrino energy. It will be considered a parameter for all events 
	static FourVector P;	//Incoming nucleon energy. Also a parameter for all events 

	Event(ThreeVector,ThreeVector,ThreeVector);	//An event given by outgoing 3-momenta
	//There are four but one is fixed by conservation laws 
	//The given ones are outgoing muon, positron, and nucleon. The fixed one is the outgoing neutrino
	//Four vectors are also fixed by mass shell requiring 
	//E = sqrt(p^2 + m^2)
	//These will be generated from the given three vectors 

	void Print();	//This prints the event to cout 

protected:
	FourVector pn;	//Outgoing neutrino momentum
	FourVector pm;	//Outgoing muon momentum
	FourVector pe;	//Outgoin positron momentum
	FourVector PN;	//Outgoing nucleon momentum
};
*/

double getEnergy(Vector, double);	//This computes the energy of a momentum p and mass m

#endif