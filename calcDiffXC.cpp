//This file is solely for the implementation of the cross section calculation
//Jonathan Curtis
//02/18/2016

#include "TridentEvent.h"

double calcDiffXC(TridentEvent event){
	//This computes the differential cross section for a given event
	//It has access to the protected members of the event, which are the outgoing momenta
	//We fist define some useful kinematic variables
	FourVector q, p1, p2, p3, p4, p;
	p1 = (TridentEvent::P1);	//The incoming neutrino momentum
	p2 = event.P2;	//Outgoing neutino momentum
	p3 = event.P3;	//Outgoing positron momentum
	p4 = event.P4; 	//Outgoing muon momentum
	p = TridentEvent::P;	//Incoming nucleon momentum

	//Energies used in normalization factors in the differential cross section
	double E2 = event.E2; 
	double E3 = event.E3;	
	double E4 = event.E4;
	double Ef = event.Ef;

	q = p - event.Pf;	//The momentum transfer	
	double q2 = q*q;	//The magnitude of the momentum transfer
	double vec_q_mag = q.getSpatialMagnitude();	//The magnitude of the spatial part of q

	double M = constant::nuclear_mass;	//The mass of the target nucleon
	double m3 = constant::posi_mass;	//These masses are used occasionally in the matrix elements
	double m4 = constant::muon_mass;

	double M2 = M*M;	//The mass squared of the target nucleon
	double m32 = m3*m3;	//The positron mass squared
	double m42 = m4*m4;	//Same for muon

	double x = q2/(4.0*M2);	//This variable is used in the nuclear matrix element for single nucleon reactions 
	double G = 1.0/std::pow(1.0+4.88*x,2);	//This is derived from x and is also used in the single nucleon cross sections 

	double mup = constant::mu_T_proton;		//The factor of mu_T for protons which is used in the single nucleon interactions 
	double mun = constant::mu_T_neutron;	//Same but for mu_T of the neutron

	double zp = constant::z_proton;	//The charge of the nucleon (used in the single nucleon interactions)
	double zn = constant::z_neutron;//Same but for neutron

	double GEp = zp*G;	//The proton electric form factor used in incoherent single nucleon scattering
	double GMp = mup*G;	//This is magnetic proton form factor used in the proton incoherent interactions
	double GEn = zn*G;	//The electric form factor for the neutron, used in incoherent reactions
	double GMn = mun*G;	//This is used in the neutron incoherent interactions

	double a = constant::RMS_nuclear_radius;	//This is used in the form factor calculation 
	double Fq2 = std::exp(- std::pow(a,2)*q2/6.0);	//This is F(q^2) which is the form factor for the coherent nuclear interaction
													//We use the form e^(-a^2q^2/6) which is valid in the low q2 regime

	double Z = double(constant::nuclear_Z);	//The number of protons (needed for weighting the various differential cross sections)
	double N = double(constant::nuclear_A-constant::nuclear_Z);	//The number of neutrons " "

	double Xq = fermiExclusionFactor(vec_q_mag);	//The exclusion principle weighting factor

	double D3 = q2 - 2.0*q*p3;	//Used in the leptonic matrix element calculations
	double D4 = q2 - 2.0*q*p4;	//Same but with p4 

	//This is the six-way product of four-vector that appears in the matrix element calculations 
	double sixprod 	= 	2.0*(p2*p4)*FourVector::QuadProd(p4,p1,p3,p) 
					+ 	2.0*(p1*p3)*FourVector::QuadProd(p4,p2,p3,p)
					+	m42*FourVector::QuadProd(p2,p1,p3,p)
					+	m32*FourVector::QuadProd(p4,p2,p1,p);

	//Now for the matrix element contractions 
	//We have two contractions that matter 
	//P*L*P and Tr(L)
	//We start with Tr(L)
	//Lovseth equation (A2)
	double TrL = std::pow(2.0,9)*(
				((p1*p3)/(D4*D4))*( 2.0*(q*p2)*(m42 - q*p4) + (p2*p4)*(q2-2.0*m42))
			+	((p2*p4)/(D3*D3))*( 2.0*(q*p1)*(m32 - q*p3) + (p1*p3)*(q2-2.0*m32))
			-	(1.0/(D3*D4))*(	q2*FourVector::QuadProd(p3,p1,p2,p4) 
								+ 4.0*(p2*p4)*(p1*p3)*(p3*p4)
								- 2.0*(p3*p4)*FourVector::QuadProd(p1,p3,p4,q)
								- 2.0*(p1*p3)*FourVector::QuadProd(p2,p4,p3,q)
								)
			);


	//Now for the P.L.P contraction
	//Lovseth equation (A1)
	double PLP = std::pow(2.0,8)*( 
		(4.0*(p*p3)*(p1*q)*(q*p) - M2*q2*(p1*p3) + 2.0*M2*(q*p1)*(q*p3) - 2.0*(p*p3)*(q2)*(p1*p) )*(p2*p4)/(D3*D3)
	+	(4.0*(p*p4)*(p2*q)*(q*p) - M2*q2*(p2*p4) + 2.0*M2*(q*p2)*(q*p4) - 2.0*(p*p4)*(q2)*(p2*p) )*(p1*p3)/(D4*D4)
	-	4.0*( (p1*p3)*(p4*p/D4-p3*p/D3) + FourVector::QuadProd(p1,p3,p,q)/D3 )*( (p2*p4)*(p3*p/D3-p4*p/D4) + FourVector::QuadProd(p2,p4,p,q)/D4 )
	-	(2.0*M2*((p2*q)*((p2*p1)*(p3*q)-(p2*p3)*(p1*q))-(p2*q)*((p4*p1)*(p3*q)-(p4*p3)*(p1*q))+q2*((p2*p3)*(p4*p1)-(p4*p3)*(p2*p1))) 
		-2.0*(q*p)*sixprod
		-2.0*q2*((p*p4)*(p2*p1)*(p3*p)-(p*p4)*(p2*p3)*(p1*p)-(p*p3)*(p1*p4)*(p2*p)+(p*p2)*(p4*p3)*(p1*p))
		)/(D3*D4)
		);

	//Now we contract the hadronic and leptonic matrix element contracitons 
	//There is a coherent term, a proton term, and a neutron term
	//We compute the coherent term first 
	//There are also factors of energy and propagators that go into the term
	//We have 
	//dXC_coh = 1/q2^2 *1/(E2 E3 E4 E') * L_ab H_ab(coherent)
	//With the effective part of the coherent matrix element being
	//H_ab(coherent) = 4Fq2^2 P_a P_b

	double dXC_coh = 4.0*(1.0/(q2*q2))*(1.0/(E2*E3*E4*Ef))*Fq2*Fq2*PLP;

	//Now we compute the proton/neutron contribution
	//The single nucleon contribution is 
	//1/Energies 1/q2^2 (4 PLP (GE^2 + x GM^2)/(1+x) + q2*GM^2 TrL)
	/*
	double dXC_pro = (1.0/(q2*q2))*(1.0/(E2*E3*E4*Ef))*( 4.0*PLP*(GEp*GEp + x*GMp*GMp)/(1.0+x) + q2*GMp*GMp*TrL );
	double dXC_neu = (1.0/(q2*q2))*(1.0/(E2*E3*E4*Ef))*( 4.0*PLP*(GEn*GEn + x*GMn*GMn)/(1.0+x) + q2*GMn*GMn*TrL );	
	*/

	//The equation for total differential cross section is 
	//dXC = dXC_coherent + X(q)*(Z*dXC_pro + N*dXC_neu)
	//We currently test with just the coherent part 
	double dXC = dXC_coh;
	return dXC;
}

double fermiExclusionFactor(double q){
	double p_f = constant::p_fermi;	//The constant used in the exclusion principle factor

	double y = .5*q/p_f;	//This is the variable used in the calculations 

	if(y>1.0){
		return 1.0;
	}
	else{
		return 1.5*y - .5*y*y*y;
	}
}






