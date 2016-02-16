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
extern const double nu_incomingE;	//Incoming neutrino beam energy 
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
	static double QuadProd(FourVector,FourVector,FourVector,FourVector); //Returns the "quad product" of four four-vectors defined in the appendix A
};

//////
//TRIDENT EVENT CLASS 
//////
/*We label four-vectors for momentum with a capital P
We have 
P0 = initial nucleon momentum = P in Lovseth
P1 = initial neutrino momentum = p1 in Lovseth
Pf = final nucleon momentum = P' in Lovseth
P2 = final neutrino momentum = p2 in Lovseth
P3 = positron momentum = p3 in Lovseth
P4 = muon momentum = p4 in Lovseth
*/
class TridentEvent{
public:
	static FourVector P0;	//Incoming neutrino energy. It will be considered a parameter for all events 
	static FourVector P1;	//Incoming nucleon energy. Also a parameter for all events 
	static double TargetMass;	//The mass of the target nucleon
	static double TargetMass2;	//The mass of the target nucleon squared

	TridentEvent(Vector,Vector,Vector);	//An event given by outgoing 3-momenta
	//Arguments are 
	//P4, P3, P'
	//There are four but one is fixed by conservation laws 
	//The given ones are outgoing muon, positron, and nucleon. The fixed one is the outgoing neutrino
	//Four vectors are also fixed by mass shell requiring 
	//E = sqrt(p^2 + m^2)
	//These will be generated from the given three vectors 

	std::ostream& printTo(std::ostream&);	//This prints the event to the given ostream

protected:
	//These are the raw variables that are measureable 
	FourVector P2;	//Outgoing neutrino momentum
	FourVector P4;	//Outgoing muon momentum
	FourVector P3;	//Outgoin positron momentum
	FourVector Pf;	//Outgoing nucleon momentum

	//////
	//DERIVED VARIABLES
	//////
	//These are all computed on construction
	FourVector q;	//The q = (PN - PN0) of the reaction 
	double q2;	//The q^2 = (PN - PN0)^2 of the reacition
	double EnergyNormalizations;	//These are the 1/(E2E3E4E') that appear in the denominator of the differential cross section 
	double Nucleonx;	//The variable "x = q^2/4M^2" needed in the nuclear cross section 

	//We now try to build the leptonic matrix element
	//We first compute the D3 and D4 functions 
	double SixProd();		//Returns the six-way product for this event defined in the appendix A of Lovseth

	double D3;	//This is defined to be q^2 - 2q*p3
	double D4;	//Same, but with p4

	double calcLeptonMatPP();	//This compute the P*L*P matrix contraction for the lepton matrix element
	double LeptonMatPP;		//The result of the calculation of the P*L*P lepton matrix contraction. We compute at construction
};

//////
//MISC
//////

double getEnergy(Vector, double);	//This computes the energy of a momentum p and mass m

int calcPDG(int,int);	//This computes the PDG code for a nucleon of given atomic number and mass

#endif












