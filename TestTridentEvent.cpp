//This is to test out all the classes in the TridentEvent package 
//Jonathan Curtis
//02/13/2016

#include "TridentEvent.h"
#include <ctime>

//All energies/momenta/masses are in units of GeV (c = 1)

int main(){
	Vector p1(12.0,3.0,-4.0);
	Vector p2(1.0,1.0,1.0);
	Vector p3(0.0,0.0,10.0);

	std::clock_t begin = std::clock();

	//We do a time trial to see how long it takes to generate 10^6 event samples 
	int num_sample = 1e6;

	for(int i = 0; i < num_sample; i++){
		TridentEvent Event(p1,p2,p3);
	}

	std::clock_t end = std::clock();

	double elapsed_secs = double(end - begin)/CLOCKS_PER_SEC;

	std::cout<<"Elapsed time: "<<elapsed_secs<<"(s)"<<std::endl;

	return 0;
}