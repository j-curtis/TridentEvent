//This is a class for an event for the Trident reaction
//It will be used to generate an MC generator and is based on the four-vector class
//Jonathan Curtis
//02/12/2016

#ifndef EVENT_H
#define EVENT_H

#include <cmath>
#include <iostream>
#include <fstream>

//////
//GLOBAL CONSTANTS
//////
namespace constant
{
extern const double muon_mass;	//mass of muon
extern const double posi_mass;	//mass of positron/electron
extern const double proton_mass;//mass of proton
extern const double neutron_mass;	//mass of neutron
extern const double nuclear_num;//Mass of nucleus (number of protons)
extern const double nuclear_mass;//Mass of nucleus with given nuclear number 
extern const int 	muon_PDG;	//PDG code of muon
extern const int 	posi_PDG;	//PDG code of positron
extern const int 	nue_PDG;	//PDG code of electron neutrino
extern const int 	num_PDG;	//PDG code of muon neutrino
extern const int 	nuke_PDG;	//PDG code of selected nucleus 
}
/////
//3-VECTOR CLASS
/////
class Vector{
protected:
	double x;	//the x component of the vector
	double y; 	//the y component
	double z;	//the z component 

public:
	Vector(double, double, double);	//This creates a three-vector given components as x,y,z
	Vector();	//This generates a null (zero vector)
	void print();	//Prints the vector to cout
	friend class TridentEvent;	//This is so we can properly access members to print the events

	//Operators
	double operator *(const Vector&);	//The normal dot-product	
	Vector operator +(const Vector&);	//This adds two three vectors component-wise
	friend Vector operator*(const double,const Vector&);	//This is for scalar multiplication from the left 
	Vector operator -(const Vector& rh);	//This subtracts two three vectors component wise 
};


//////
//4-VECTOR CLASS
//////
class FourVector{
protected:
	double t;	//the time component (x0)
	Vector v;	//The spatial components (x1,x2,x3) 

public:
	FourVector(double,Vector);	//This creates a four-vector given time component and spatial vector
	FourVector();	//This is a constructor for a 0 four-vector 
	void print();	//Prints the vector to cout
	friend class TridentEvent;	//We need this to be a friend class so we can properly access components for printing the events 

	double operator *(const FourVector&);	//The Lorentz invaraint contraction of two fourvectors
	FourVector operator +(const FourVector&);	//This adds two fourvectors component-wise
	friend FourVector operator*(const double,const FourVector&);	//This is for scalar multiplication from the left 
	FourVector operator -(const FourVector& rh);	//This subtracts two fourvectors component wise 
};

double getEnergy(Vector, double);	//This computes the energy of a momentum p and mass m

//////
//TRIDENT EVENT CLASS 
//////
class TridentEvent{
public:
	static FourVector Pn0;	//Incoming neutrino energy. It will be considered a parameter for all events 
	static FourVector PN0;	//Incoming nucleon energy. Also a parameter for all events 

	TridentEvent(Vector,Vector,Vector);	//An event given by outgoing 3-momenta
	//There are four but one is fixed by conservation laws 
	//The given ones are outgoing muon, positron, and nucleon. The fixed one is the outgoing neutrino
	//Four vectors are also fixed by mass shell requiring 
	//E = sqrt(p^2 + m^2)
	//These will be generated from the given three vectors 

	std::ostream& printTo(std::ostream&);	//This prints the event to the given ostream

protected:
	FourVector Pn;	//Outgoing neutrino momentum
	FourVector Pm;	//Outgoing muon momentum
	FourVector Pe;	//Outgoin positron momentum
	FourVector PN;	//Outgoing nucleon momentum
};

//////
//MISC
//////
int calcPDG(int,int);	//This computes the PDG code for a nucleon of given atomic number and mass
#endif












