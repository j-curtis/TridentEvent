//This is to test out all the classes in the TridentEvent package 
//Jonathan Curtis
//02/13/2016

#include "TridentEvent.h"

//All energies/momenta/masses are in units of GeV (c = 1)

int main(){
	Vector p4(0.0,.197,3.7);
	Vector p3(86.4e-3,0.0,.704);
	Vector pf(0.0,0.0,5.0e-3);

	TridentEvent Event(p4,p3,pf);
	Event.printTo(std::cout);
	std::cout<<std::endl;
	double dXC = calcDiffXC(Event);
	std::cout<<dXC<<std::endl;
	return 0;
}