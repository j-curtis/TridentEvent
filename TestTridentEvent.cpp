//This is to test out all the classes in the TridentEvent package 
//Jonathan Curtis
//02/13/2016

#include "TridentEvent.h"

int main(){
	Vector p1(12.0,3.0,-4.0);

	FourVector P1(13.0,p1);

	P1.Print();

	std::cout<<std::endl;
	std::cout<<P1*P1<<std::endl;

	(P1+P1).Print();
	std::cout<<std::endl;

	(P1-(3.0*P1)).Print();
	std::cout<<std::endl;

	return 0;
}