// See LICENCE file at project root
#include <iostream>
#include "Utils/FTic.hpp"
#include "Utils/FParameterNames.hpp"

#include <cstdlib>
#include <unistd.h>

/**
* Here we show an example of using FTic
*/

int main(int argc, char ** argv){
    FHelpDescribeAndExit(argc, argv,
                         "Only the code is interesting in order to understand the use of timers.");
    ///////////////////////What we do/////////////////////////////
    std::cout << ">> This executable is useless to execute.\n";
    std::cout << ">> It is only interesting to wath the code to understand\n";
    std::cout << ">> how to use FTic time counter.\n";
    //////////////////////////////////////////////////////////////
    {
	FTic counter;	
	counter.tic();
	usleep(1500000);
	//Sleep(1500); //on windows
	counter.tac();
	std::cout << counter.elapsed() << " (s)\n";
    }
    {
        FTic counter;
        usleep(1500000);
        //Sleep(1500); //on windows
        std::cout << counter.tacAndElapsed() << " (s)\n";
    }
    {
        FTic counter;
        usleep(1500000);
        //Sleep(1500); //on windows
        counter.tac();
        counter.tic();
        usleep(1500000);
        //Sleep(1500); //on windows
        std::cout << counter.tacAndElapsed() << " (s)\n";
        std::cout << counter.cumulated() << " (s)\n";
    }
    return 0;
}

