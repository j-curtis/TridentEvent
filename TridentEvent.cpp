//The definition file for the TridentEvent class
//Jonathan Curtis
//02/12/2016

#include "TridentEvent.h"

//We use mass/energy/momentum units of GeV
const double muon_mass = .1056583715;
const double posi_mass = .0005109989; 
const double proton_mass = .938272046;

/*
3-vector class methods 
*/

Vector::Vector(double arg1, double arg2, double arg3): x(arg1),y(arg2),z(arg3) {}

double Vector::operator *(const Vector& rh){ return ( (x*rh.x) + (y*rh.y) + (z*rh.z) ); }

Vector Vector::operator +(const Vector& rh){ return Vector( x+rh.x , y+rh.y , z+rh.z ); }

Vector operator *(const double c, const Vector& rh){ return Vector( c*rh.x , c*rh.y , c*rh.z );	}

Vector Vector::operator -(const Vector& rh){ return ( (*this) + ( (-1.0)*rh ) ); }

void Vector::Print(){
	std::cout<<"("<<x<<", "<<y<<", "<<z<<")";
}

//This computes the on-shell energy of a particle with momentum p and mass m
double getEnergy(Vector p, double m){
	double E2 = (p*p) + m*m;
	return std::sqrt(E2);
}

/*
4-vector methods 
*/

FourVector::FourVector(double arg0, Vector arg1): t(arg0), v(arg1) {}

double FourVector::operator *(const FourVector& rh){
	double sumv = (v)*(rh.v); 

	double sumt = (t)*(rh.t);

	return sumt - sumv;
}

FourVector FourVector::operator +(const FourVector& rh){ return FourVector( t+rh.t , v+rh.v ); }

FourVector operator *(const double c, const FourVector& rh){ return FourVector( c*rh.t , c*rh.v ); }

FourVector FourVector::operator -(const FourVector& rh){ return ( (*this) + ( (-1.0)*rh ) ); }

void FourVector::Print(){
	//Crude print method for a four vector 
	std::cout<<"("<<t<<", ";
	v.Print();
	std::cout<<")";
}

/*
TridentEvent methods 
*/
/*
TridentEvent::TridentEvent(ThreeVector arg1, ThreeVector arg2, ThreeVector arg3){
	//The arguments supplied are the outgoing 
	//muon momentum
	//electron momentum
	//nucleon momentum

	//First we compute the energies of the given three-vectors and masses 
	//We use a massless neutrino 
	//We assume an outgoing nucleon of 12 proton masses (carbon atom)

	//muon energy
	double Em = getEnergy(arg1,muon_mass);

	//posi energy
	double Ee = getEnergy(arg2,posi_mass);

	//nucleon energy 
	double EN = getEnergy(arg3,12.0*proton_mass);

	//Now we generate the fourvectors 
	pm = FourVector(arg1,Em);
	pe = FourVector(arg2,Ee);
	pN = FourVector(arg3,EN);

}*/