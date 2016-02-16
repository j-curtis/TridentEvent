//The definition file for the TridentEvent class
//Jonathan Curtis
//02/12/2016

#include "TridentEvent.h"

//GLOBAL CONSTANT DEFINITIONS
const double 	constant::muon_mass = .1056583715;	//[GeV]
const double 	constant::posi_mass = .0005109989; 	//[GeV]
const double 	constant::proton_mass = .938272046;	//[GeV]
const double 	constant::neutron_mass = .9395654133;	//[GeV]
const int 		constant::nuclear_A = 12;		//[1]
const int 		constant::nuclear_Z = 6;		//[1]
const double 	constant::nuclear_mass = calcNuclearMass(constant::nuclear_Z,constant::nuclear_A);	//[GeV]
const int 		constant::muon_PDG = 13;	//PDG 
const int 		constant::posi_PDG = -11;	//PDG 
const int 		constant::nue_PDG = 12;	//PDG
const int 		constant::num_PDG = 14;	//PDG 
const int 		constant::nuke_PDG = calcPDG(constant::nuclear_Z,constant::nuclear_A);	//PDG generated from the atomic number/mass
const double 	constant::nu_incomingE = 8.0;	//[GeV]
const double 	constant::Fermimeter = 1.0/.1973269788;	//[1/GeV] 1 Fermimeter = 10^-15 m in units of 1/GeV where hbar = c = 1
const double 	constant::RMS_nuclear_radius = calcFormFactorConstant(constant::nuclear_A);//[1/GeV] The rms radius of the nucleus in 1/GeV used in the form factor calculation
const double 	constant::mu_T_proton = 2.792;	//[1] The factor of mu_T for protons 
const double	constant::mu_T_neutron = -1.913;//[1] The factor of mu_T for neutrons
const double	constant::z_proton = 1.0;		//[1] This is one for protons and zero for neutrons
const double 	constant::z_neutron = 0.0;		//[1] This is one for protons and zero for neutrons
const double 	constant::p_fermi = .235;		//[GeV] The constant p_f used in Lovseth equation 20 for the Fermi exclusion factor 
/*
VECTOR CLASS DEFINITIONS
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

/*
4-VECTOR CLASS DEFINITIONS
*/

FourVector::FourVector(double arg0, Vector arg1): t(arg0), v(arg1) {}

FourVector::FourVector(): t(0.0), v(){}

double FourVector::operator *(const FourVector& rh){
	double sumv = (v)*(rh.v); 

	double sumt = (t)*(rh.t);

	return -sumt + sumv;
}

FourVector FourVector::operator +(const FourVector& rh){ return FourVector( t+rh.t , v+rh.v ); }

FourVector operator *(const double c, const FourVector& rh){ return FourVector( c*rh.t , c*rh.v ); }

FourVector FourVector::operator -(const FourVector& rh){ return ( (*this) + ( (-1.0)*rh ) ); }

double FourVector::QuadProd(FourVector v1, FourVector v2, FourVector v3, FourVector v4){
	//In the appendix A1 they introduce the notation 
	//(abcd) = (a.b)(c.d) - (a.c)(b.d) + (a.d)(b.c)
	//This computes this quadruple product 
	double term1 = (v1*v2)*(v3*v4);	//the first term 
	double term2 = (v1*v3)*(v2*v4);	//the second term
	double term3 = (v1*v4)*(v2*v3);	//the third term

	return term1 - term2 + term3;
}


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
TridentEvent::TridentEvent(Vector arg1, Vector arg2, Vector arg3): P4(),P3(),Pf(),P2(),q() {
	//The arguments supplied are the outgoing 
	//muon momentum
	//electron momentum
	//nucleon momentum

	//First we compute the energies of the given three-vectors and masses 
	//We use a massless neutrino 
	//We assume an outgoing nucleon of 12 proton masses (carbon atom)

	//muon energy
	double E4 = getEnergy(arg1,constant::muon_mass);

	//posi energy
	double E3 = getEnergy(arg2,constant::posi_mass);

	//nucleon energy 
	double Ef = getEnergy(arg3,constant::nuclear_mass);

	//Now we generate the fourvectors 
	P4 = FourVector(E4,arg1);
	P3 = FourVector(E3,arg2);
	Pf = FourVector(Ef,arg3);

	//Now we generate the four-vector for the outgoing neutrino using conservation laws 
	P2 = ( (P1 + P0) - (P4 + P3 + Pf) );


	//Now we precompute some of the kinematic variables 
	q = (P0 - Pf);	//The q of the reaction
	q2 = q*q;		//The q^2 is easy to compute 
	EnergyNormalizations = 1.0/(E4*E3*Ef*(P2.t));	//The 1/(E2 E3 E4 E') that appears in the denominator of the cross section calculations
	Nucleonx = q2/(4.0 * TridentEvent::TargetMass2);	//The variable "x = q^2 / 4M^2" needed in the single nucleon scattering case 
	NucleonG = std::pow(1.0+4.88*Nucleonx,-2);	//The variable G = 1/(1+4.88x)^2 used in the single nucleon cross section (Lovseth 11)
	//Now we compute the leptonic matrix elements 
	D3 = q2 - 2.0*(q*P3);
	D4 = q2 - 2.0*(q*P4);

	LeptonMatPP = calcLeptonMatPP();	//We construct the P*L*P lepton matrix element and store it for further use 
	LeptonMatTr = calcLeptonMatTr();	//We construct the L_ab delta_ab lepton matrix element contraction and store it for later use 
	FormFactor = std::exp(-std::pow(constant::RMS_nuclear_radius,2)*q2/6.0);	//The nuclear form-factor for this given q^2

	//Coherent matrix element contraction is 4 F^2(q^2) P*L*P
	coherentMatrixEl = 4.0*std::pow(FormFactor,2.0)*LeptonMatPP;

	//proton matrix element 
	protonMatrixEl = 4.0*NucleonG*NucleonG*( LeptonMatPP*( std::pow(constant::mu_T_proton*Nucleonx,2) + constant::z_proton ) ) 
					+ NucleonG*NucleonG	*LeptonMatTr*q2*constant::mu_T_proton*constant::mu_T_proton;

	neutronMatrixEl = 4.0*NucleonG*NucleonG*( LeptonMatPP*( std::pow(constant::mu_T_neutron*Nucleonx,2) + constant::z_neutron ) ) 
					+ NucleonG*NucleonG	*LeptonMatTr*q2*constant::mu_T_neutron*constant::mu_T_neutron;	

	FermiExclusionFactor = calcFermiExclusionFactor();

	//In Lovseth equation 19 they give the formula
	//dsigma = disgmacoherent + xclusionfactor*((A-Z)dsigmaneutron + Z dsigmaproton)
	diffXC = (coherentMatrixEl + FermiExclusionFactor*( (constant::nuclear_A-constant::nuclear_Z)*neutronMatrixEl + constant::nuclear_Z*protonMatrixEl))
		*(1.0/(EnergyNormalizations*q2*q2));
}	

double TridentEvent::SixProd(){
	//Returns the six-way product defined in the appendix
	//This is 
	//2p2*p4(p4p1p3P) + 2p1*p3(p4p2p3P)+m4^2(p2p1p3P)+m3^2(p4p2p1P)
	double term1 = 2.0*(P2*P4)*FourVector::QuadProd(P4,P1,P3,P0);
	double term2 = 2.0*(P1*P3)*FourVector::QuadProd(P4,P2,P3,P0);
	double term3 = constant::muon_mass * constant::muon_mass * FourVector::QuadProd(P2,P1,P3,P0);
	double term4 = constant::posi_mass * constant::posi_mass * FourVector::QuadProd(P4,P2,P1,P0);

	return term1+term2+term3+term4;
}

double TridentEvent::calcLeptonMatPP(){
	//Lovseth Equation (A1)
	//The equation is extremely long
	//We break it up into multiple pieces and construct the total equation from these 
	double sum = 0.0;

	double line1 = ( (P2*P4)/(D3*D3) )*(
		4.0*(P0*P3)*(P1*q)*(q*P0) 
		-TargetMass2*q2*(P1*P3)
		+2.0*TargetMass2*(q*P1)*(q*P3)
		-2.0*(P0*P3)*q2*(P1*P0)
		);

	//line 2 is just line 1 with P2,P4 replacing P1,P3
	double line2 = ( (P1*P3)/(D4*D4) )*(
		4.0*(P0*P4)*(P2*q)*(q*P0) 
		-TargetMass2*q2*(P2*P4)
		+2.0*TargetMass2*(q*P2)*(q*P4)
		-2.0*(P0*P4)*q2*(P2*P0)
		);

	double line3 = (1.0/(D3*D4))*(
			2.0*TargetMass2*(
			  (P4*q)*( (P2*P1)*(P3*q)  - (P2*P3)*(P1*q) ) 
			- (P2*q)*( (P4*P1)*(P3*q)  - (P4*P3)*(P1*q) )  
			+ 	  q2*( (P2*P3)*(P2*P1) - (P4*P3)*(P2*P1))
			)

			- 2.0*q*P0*SixProd() 
			- 2.0*q2*( (P0*P4)*(P2*P1)*(P3*P0) - (P0*P4)*(P2*P3)*(P1*P0) - (P0*P3)*(P1*P4)*(P2*P0) + (P0*P2)*(P4*P3)*(P1*P0) )
		);

	double line4 = 4.0 
	*( (P1*P3)*(P4*P0/D4 - P3*P0/D3) + FourVector::QuadProd(P1,P3,P0,q) )
	*( (P2*P4)*(P3*P0/D3 - P4*P0/D4) + FourVector::QuadProd(P2,P4,P0,q) );

	sum = line1 + line2 - line3 - line4;

	return std::pow(2.0,8)*sum;
}

double TridentEvent::calcLeptonMatTr(){
	//Lovseth Equation (A2)
	//The equation is extremely long
	//We break it up into multiple pieces and construct the total equation from these 
	double sum = 0.0;

	double term1 = 	((P1*P3)/(D4*D4))*( 2.0*(q*P2)*(constant::muon_mass*constant::muon_mass - q*P4)
					+ (P2*P4)*(q2-2.0*constant::muon_mass*constant::muon_mass));

	double term2 = 	((P2*P4)/(D3*D3))*( 2.0*(q*P1)*(constant::posi_mass*constant::posi_mass - q*P3)
					+ (P1*P3)*(q2-2.0*constant::posi_mass*constant::posi_mass));

	double term3 = 	((1.0)/(D3*D4))*(
					q2*FourVector::QuadProd(P3,P1,P2,P4) 
					+4.0*(P2*P4)*(P1*P3)*(P3*P4)
					-2.0*(P2*P4)*FourVector::QuadProd(P1,P3,P4,q)
					-2.0*(P1*P3)*FourVector::QuadProd(P2,P4,P3,q)
					);

	sum = term1 + term2 - term3;

	return std::pow(2.0,9)*sum;
}

double TridentEvent::calcFermiExclusionFactor(){
	//This is the function X(|q|) used in Lovseth equation 20
	if( std::sqrt(q.v*q.v) >= 2.0*constant::p_fermi){
		return 1.0;
	}
	else{
		return 1.5*(std::sqrt(q.v*q.v)/(2.0*constant::p_fermi)) - .5*std::pow( std::sqrt(q.v*q.v)/(2.0*constant::p_fermi),3);
	}
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
	o<<constant::muon_PDG<<" 0.0 0.0 0.0 "<<P4.v.x<<" "<<P4.v.y<<" "<<P4.v.z<<" "<<P4.t<<std::endl;
	o<<constant::posi_PDG<<" 0.0 0.0 0.0 "<<P3.v.x<<" "<<P3.v.y<<" "<<P3.v.z<<" "<<P3.t<<std::endl;
	o<<constant::nuke_PDG<<" 0.0 0.0 0.0 "<<Pf.v.x<<" "<<Pf.v.y<<" "<<Pf.v.z<<" "<<Pf.t<<std::endl;
	o<<constant::nue_PDG <<" 0.0 0.0 0.0 "<<P2.v.x<<" "<<P2.v.y<<" "<<P2.v.z<<" "<<P2.t<<std::endl;
	return o;
}

double TridentEvent::getDiffXC(){ return diffXC; }

//Initializing static member variables
//incoming neutrino momentum vector 
FourVector TridentEvent::P1(constant::nu_incomingE, Vector(0.0, 0.0, constant::nu_incomingE) );
//incoming nucleon momentum vector
FourVector TridentEvent::P0(constant::nuclear_mass, Vector());
//The target nucleon mass 
double TridentEvent::TargetMass = constant::nuclear_mass;	 
double TridentEvent::TargetMass2 = constant::nuclear_mass*constant::nuclear_mass;


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

//This computes the on-shell energy of a particle with momentum p and mass m
double getEnergy(Vector p, double m){
	double E2 = (p*p) + m*m;
	return std::sqrt(E2);
}

//This computes the form-factor constant in units of 1/GeV 
double calcFormFactorConstant(int A){
	return (.58 + .82*std::pow(A,1.0/3.0) )*constant::Fermimeter;
}

double calcNuclearMass(int Z, int A){
	return constant::proton_mass*Z + (A-Z)*constant::neutron_mass;
}

