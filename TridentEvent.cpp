//The definition file for the TridentEvent class
//Jonathan Curtis
//02/12/2016

#include "TridentEvent.h"

//GLOBAL CONSTANT DEFINITIONS
const double 	constant::muon_mass = .1056583715;	//[GeV]
const double 	constant::posi_mass = .0005109989; 	//[GeV]
const double 	constant::proton_mass = .938272046;	//[GeV]
const double 	constant::neutron_mass = .9395654133;	//[GeV]
const double 	constant::nuclear_num = 12.0;		//[1]
const double 	constant::nuclear_mass = constant::nuclear_num*constant::proton_mass;	//[GeV]
const int 		constant::muon_PDG = 13;	//PDG 
const int 		constant::posi_PDG = -11;	//PDG 
const int 		constant::nue_PDG = 12;	//PDG
const int 		constant::num_PDG = 14;	//PDG 
const int 		constant::nuke_PDG = calcPDG(6,12);	//PDG generated from the atomic number/mass


/*
3-VECTOR CLASS DEFINITIONS
*/

Vector::Vector(double arg1, double arg2, double arg3): x(arg1),y(arg2),z(arg3) {}

Vector::Vector(): x(0.0), y(0.0), z(0.0) {}

double Vector::operator *(const Vector& rh){ return ( (x*rh.x) + (y*rh.y) + (z*rh.z) ); }

Vector Vector::operator +(const Vector& rh){ return Vector( x+rh.x , y+rh.y , z+rh.z ); }

Vector operator *(const double c, const Vector& rh){ return Vector( c*rh.x , c*rh.y , c*rh.z );	}

Vector Vector::operator -(const Vector& rh){ return ( (*this) + ( (-1.0)*rh ) ); }

void Vector::print(){
	std::cout<<"("<<x<<", "<<y<<", "<<z<<")";
}

//This computes the on-shell energy of a particle with momentum p and mass m
double getEnergy(Vector p, double m){
	double E2 = (p*p) + m*m;
	return std::sqrt(E2);
}

/*
4-VECTOR CLASS DEFINITIONS
*/

FourVector::FourVector(double arg0, Vector arg1): t(arg0), v(arg1) {}

FourVector::FourVector(): t(0.0), v(){}

double FourVector::operator *(const FourVector& rh){
	double sumv = (v)*(rh.v); 

	double sumt = (t)*(rh.t);

	return sumt - sumv;
}

FourVector FourVector::operator +(const FourVector& rh){ return FourVector( t+rh.t , v+rh.v ); }

FourVector operator *(const double c, const FourVector& rh){ return FourVector( c*rh.t , c*rh.v ); }

FourVector FourVector::operator -(const FourVector& rh){ return ( (*this) + ( (-1.0)*rh ) ); }

void FourVector::print(){
	//Crude print method for a four vector 
	std::cout<<"("<<t<<", ";
	v.print();
	std::cout<<")";
}

/*
TRIDENT EVENT CLASS DEFINITIONS 
*/

//We need to initialize all elements that don't have a default constructor. We then give them the appropriate valeus latter 
TridentEvent::TridentEvent(Vector arg1, Vector arg2, Vector arg3): Pm(),Pe(),PN(),Pn() {
	//The arguments supplied are the outgoing 
	//muon momentum
	//electron momentum
	//nucleon momentum

	//First we compute the energies of the given three-vectors and masses 
	//We use a massless neutrino 
	//We assume an outgoing nucleon of 12 proton masses (carbon atom)

	//muon energy
	double Em = getEnergy(arg1,constant::muon_mass);

	//posi energy
	double Ee = getEnergy(arg2,constant::posi_mass);

	//nucleon energy 
	double EN = getEnergy(arg3,constant::nuclear_mass);

	//Now we generate the fourvectors 
	Pm = FourVector(Em,arg1);
	Pe = FourVector(Ee,arg2);
	PN = FourVector(EN,arg3);

	//Now we generate the four-vector for the outgoing neutrino using conservation laws 
	Pn = ( (Pn0 + PN0) - (Pm + Pe + PN) );
}

std::ostream& TridentEvent::printTo(std::ostream& o){
	//This prints in the format accepted by the MINERvA file submissiomn format
	//We list an event as a block of particles 
	//Each line (for each particle) is formatted as
	//PDG X Y Z PX PY PZ E
	//We print as 
	//muon
	//positron
	//nucleon
	//neutrino 
	//We start each vertex at the origin for simplicity
	o<<constant::muon_PDG<<" 0.0 0.0 0.0 "<<Pm.v.x<<" "<<Pm.v.y<<" "<<Pm.v.z<<" "<<Pm.t<<std::endl;
	o<<constant::posi_PDG<<" 0.0 0.0 0.0 "<<Pe.v.x<<" "<<Pe.v.y<<" "<<Pe.v.z<<" "<<Pe.t<<std::endl;
	o<<constant::nuke_PDG<<" 0.0 0.0 0.0 "<<PN.v.x<<" "<<PN.v.y<<" "<<PN.v.z<<" "<<PN.t<<std::endl;
	o<<constant::nue_PDG <<" 0.0 0.0 0.0 "<<Pn.v.x<<" "<<Pn.v.y<<" "<<Pn.v.z<<" "<<Pn.t<<std::endl;
	return o;
}

//Initialize incoming neutrino vector 
FourVector TridentEvent::Pn0(10.0, Vector(0.0, 0.0, 10.0) );
FourVector TridentEvent::PN0(constant::nuclear_mass, Vector());


//////
//MISC
//////
int calcPDG(int Z,int A){
	//The PDG code for a nucleon (with isomer level 0 and no strange quarks)
	int digit_1 = 0;	//The isomer level
	int digit_234 = A;	//The digits 2-4 are the atomic mass (nucleon number)
	int digit_567 = Z;	//the digits 5-6 are the atomic number (number of protons)
	int digit_8910 = 100;	//the first three digits are 10 and then L=0 (number of strange baryons)

	//We insert appropriate factors of ten to get the code the correct number of digits
	return 1e0*digit_1 + 1e1*digit_234 + 1e4*digit_567 + 1e7*digit_8910;
}



