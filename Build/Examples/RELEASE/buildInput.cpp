#include <iostream>
#include <cmath>
using namespace std;

int main() {
	
	//int points = 10; //
	//int points = 16; //diff = 6
	//int points = 22; //diff = 6 
	//int points = 31; //diff = 9 
	//int points = 44; //diff = 13 
	//int points = 62; //diff = 18
	//int points = 88; //diff = 26
	//int points = 124; //diff = 36 
	//int points = 176; //diff = 52
	//int points = 182;  // -------------
	//int points = 183;  // --------------
	//int points = 248; 
	//int points = 351; // 123904
	//int points = 497; // 248004
	//int points = 548; // 301401
	int points = 641; // 721*721 =  /how high can we go?
	
	
	int numofpoints = (points + 1)*(points + 1);
	std::cout << "numofparticles  = " << numofpoints<< "     ,points: " << points << std::endl;
	
	double x[numofpoints];
	double y[numofpoints];
	double z[numofpoints];
	
	int counter = 0;
		for (double a = 0; a <= points; ++a)
		{		
			for (double b = 0; b <= points; ++b)
			{				
				x[counter] = (a/points);
				y[counter] = (0);
				z[counter] = (b/points);
				
				//std::cout << "x = " << x[counter]  << ", y = " << y[counter]  << ", z = " << z[counter] << std::endl;
				std::cout << x[counter]  << " " << y[counter]  << " " << z[counter] << " 0.01" <<std::endl;
				
				counter = counter + 1;
				//std::cout << "counter  = " << counter  << std::endl;
			}
		}
		
	
	
    return 0;
}