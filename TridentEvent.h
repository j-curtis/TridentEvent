//This is a class for an event for the Trident reaction
//It will be used to generate an MC generator and is based on the four-vector class
//Jonathan Curtis
//02/12/2016

#ifndef EVENT_H
#define EVENT_H

#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>	
#include <ctime>
#include <cstring>

//////
//GLOBAL CONSTANTS
//////
namespace constant
{
extern const double muon_mass;	//mass of muon
extern const double posi_mass;	//mass of positron/electron
extern const double proton_mass;//mass of proton
extern const double neutron_mass;	//mass of neutron
extern const int nuclear_A;//Number of nucleons
extern const int nuclear_Z;	//Number of protons
extern const double nuclear_mass;//Mass of nucleus (GeV) with given nuclear number 
extern const int 	muon_PDG;	//PDG code of muon
extern const int 	posi_PDG;	//PDG code of positron
extern const int 	nue_PDG;	//PDG code of electron neutrino
extern const int 	num_PDG;	//PDG code of muon neutrino
extern const int 	nuke_PDG;	//PDG code of selected nucleus 
extern const double nu_incomingE;	//Incoming neutrino beam energy 
extern const double Fermimeter;	//This is 1 Fermi (F, fm, 10^(-15) m) in units where hbar = c = 1
extern const double RMS_nuclear_radius;	//This is the RMS nuclear radius, used in the form-factor computation
extern const double mu_T_proton;		//The factor of mu_T for the proton used in the single nucleon interactions 
extern const double mu_T_neutron;		//Same but for the neutron instead of the proton
extern const double z_proton;			//This is one if it is a proton and zero for neutron
extern const double z_neutron;			//This is one for protons and zero for neutrons
extern const double p_fermi;			//This is the constant p_f used in the fermi exclusion factor
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

	friend double magnitude(Vector);	//This returns the magnitude of the passed vector
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
	double getSpatialMagnitude();	//This returns the magnitude of the spatial part of a four-vector

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
	static FourVector P;	//Incoming neutrino energy. It will be considered a parameter for all events 
	static FourVector P1;	//Incoming nucleon energy. Also a parameter for all events 
	
	TridentEvent(Vector,Vector,Vector);	//An event given by outgoing 3-momenta
	//Arguments are 
	//P4, P3, P'
	//There are four but one is fixed by conservation laws 
	//The given ones are outgoing muon, positron, and nucleon. The fixed one is the outgoing neutrino
	//Four vectors are also fixed by mass shell requiring 
	//E = sqrt(p^2 + m^2)
	//These will be generated from the given three vectors 

	void printTo(std::ostream&);	//This prints the event to the given ostream

	friend double calcDiffXC(TridentEvent);	//This externally computes the differential cross section of a given event 

protected:
	//These are the raw variables that are measureable 
	FourVector P2;	//Outgoing neutrino momentum
	FourVector P4;	//Outgoing muon momentum
	FourVector P3;	//Outgoin positron momentum
	FourVector Pf;	//Outgoing nucleon momentum

	//We make some energies protected members so they are visible to the friend function calcDiffXC
	double E2;
	double E3;
	double E4;
	double Ef;
};

//////
//MC SAMPLER CLASS
//////
/*
This class generates a sample of events via MC sampling
It will first compute the total cross section given an area of phase space
Then it will generate a sample according to 
dP = dsigma/sigma 
Finally, it will print these events to a given txt file
*/
class MCSampler{
public:
	MCSampler();	//This is just a generic constructor  
	
	void setP3Range(double,double);	//This sets the volume in P3 space that we integrate over.  
	//first parameter is max magnitude 
	//Second parameter is max opening angle 
	//We integrate symmetrically around the z axis in a cone 

	void setP4Range(double,double);	//Same but for P4 space 
	void setPfRange(double,double);	//Same but for Pf space 

protected:
	double maxP3_r;	//max |P3| 
	double maxP3_theta;	//max P3 opening angle 
	double maxP4_r;	//Same but for P4
	double maxP4_theta;
	double maxPf_r;
	double maxPf_theta;
	
};


//////
//MISC
//////

double getEnergy(Vector, double);	//This computes the energy of a momentum p and mass m

int calcPDG(int,int);	//This computes the PDG code for a nucleon of given atomic number and mass

double calcFormFactorConstant(int);	//This computes the form-factor constant for a given nuclear number

double calcNuclearMass(int,int);	//This computes the mass of a nucleas with Z protons and A nucleons

double fermiExclusionFactor(double);	//This computes the fermi exclusion principle factor used in the total cross section calculation

Vector sphericalVector(double,double,double);	//This generates a vector given a radius and polar and azimuthal angle

double getRandRange(double,double);	//This returns a randomly selected double from the given range 
#endif












