#include <iostream>
#include <cmath>
using namespace std;

int main() {

	double scalingfactor = 0.1;
	//int numofpoints = 121;  //10 --> 11 X 11
	//int numofpoints = 441;  //20 --> 21 X 21	
	//int numofpoints = 1681;	  //40 --> 41 X 41	
	//int numofpoints = 6561;	  //80 --> 81 X 81	
	int numofpoints = 25921;	  //160 --> 161 X 161		
	//float loopmax = sqrt(numofpoints) - 1;
	//std::cout << loopmax << std::endl;
	
	double x[numofpoints];
	double y[numofpoints];
	double z[numofpoints];
	
	int counter = 0;
		for (double a = 0; a <= 160; ++a)
		{		
			for (double b = 0; b <= 160; ++b)
			{				
				x[counter] = (scalingfactor*(a/16));
				y[counter] = (0);
				z[counter] = (scalingfactor*(b/16));
				
				//std::cout << "x = " << x[counter]  << ", y = " << y[counter]  << ", z = " << z[counter] << std::endl;
				std::cout << x[counter]  << " " << y[counter]  << " " << z[counter] << " 0.01" <<std::endl;
				
				counter = counter + 1;
				std::cout << "counter  = " << counter  << std::endl;
			}
		}
		
	
	
    return 0;
}