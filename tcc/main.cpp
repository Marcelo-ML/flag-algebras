#include "LocalCutParamToSDP.h"

#include <iostream>
#include <string>

void interactive_main(){
	LocalCut::ParamToSDP constraint_handler;
	
	std::string filename;

	std::cout << "\nType the filename for the local cut constraints: ";
	std::cin >> filename;



	constraint_handler.addLCCondition(filename);

	std::cout << "\n\nType the filename for the semidefinite constraints: ";
	std::cin >> filename;
	std::cout << "\n\n";

	constraint_handler.addSDCondition(filename);


	double opt_val = constraint_handler.getOptimalValue();

	std::cout << "\n\nThe optimal value is: " << opt_val << std::endl;
}

int main(int argc, char **argv) {
	if (argc != 3 && argc != 4) {
		std::cout << "Usage: " 
		          << argv[0] << " "
		          << "<local constraints file> "
		          << "<semidefinite constraints file> "
		          << "[intermediate dat-s file (default: tmp.dat-s)]"
		          << std::endl;
		return 0;
	}

	LocalCut::ParamToSDP constraint_handler;
	if (argc == 4)
		constraint_handler.setTmpFilename(std::string(argv[3]));

	std::cout << "Adding local cut constraints..." << std::endl;
	constraint_handler.addLCCondition(std::string(argv[1]));
	std::cout << "Adding semidefinite constraints..." << std::endl;
	constraint_handler.addSDCondition(std::string(argv[2]));
	double opt_val = constraint_handler.getOptimalValue();

	std::cout << "The optimal value is: " << opt_val << std::endl;

	return 0;
}