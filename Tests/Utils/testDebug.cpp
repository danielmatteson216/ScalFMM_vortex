// See LICENCE file at project root
#include "Utils/FLog.hpp"

#include "Utils/FParameterNames.hpp"

/**
* In this file we show how to use the debug module.
* Warning, in FGlobal.hpp (included in FLog.hpp) SCALFMM_USE_DEBUG might be undefined.
*/

int main(int argc, char ** argv){
    FHelpDescribeAndExit(argc, argv, "Just show how the debug and log stuff should be used.");
    ///////////////////////What we do/////////////////////////////
    std::cout << ">> This executable is useless to execute.\n";
    std::cout << ">> It is only interesting to wath the code to understand\n";
    std::cout << ">> to understand the FLog system.\n";
    //////////////////////////////////////////////////////////////
	// Print data simply
	FLOG( FLog::Controller << "Hello Wordl\n");

	// Print a variable (formated print)
	FLOG( int i = 50;)
	FLOG( FLog::Controller.writeVariableFromLine( "i", i, __LINE__, __FILE__););

	// Write a developer information
	FLOG( FLog::Controller.writeFromLine("Strange things are there!", __LINE__, __FILE__); )

	// Flush
	FLOG( FLog::Controller << FLog::Flush );

	// Change stream type
	FLOG( FLog::Controller.writeToFile("testDebug.out.temp"); )
	FLOG( FLog::Controller << "Hello Wordl 2 the return\n");

	return 0;
}



