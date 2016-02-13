//This is to test out all the classes in the TridentEvent package 
//Jonathan Curtis
//02/13/2016

#include "TridentEvent.h"

//All energies/momenta/masses are in units of GeV (c = 1)

int main(){
	Vector p1(12.0,3.0,-4.0);
	Vector p2(1.0,1.0,1.0);
	Vector p3(0.0,0.0,10.0);

	TridentEvent event(p1,p2,p3);

	std::cout<<event.printTo(std::cout);

	return 0;
}