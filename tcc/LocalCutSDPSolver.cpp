#ifndef LOCALCUTSDPSOLVER_CPP_INCLUDED
#define LOCALCUTSDPSOLVER_CPP_INCLUDED

#include "LocalCutSDPSolver.h"

double solve_sdpa_file(const std::string& filename) {
    // Assumes 'csdp' executable is in your system's PATH
    std::string command = "csdp " + filename;

    FILE* pipe = popen(command.c_str(), "r");
    if (!pipe) {
        throw std::runtime_error("popen() failed!");
    }

    char buffer[128];
    std::string result = "";
    bool solution_found = false;
    std::string primal_obj, dual_obj;
    
    while (fgets(buffer, sizeof(buffer), pipe) != nullptr) {
        result += buffer;
        std::string line = buffer;

        if (line.find("Primal objective value:") != std::string::npos) {
            // Extract the value by finding the substring after the label
            primal_obj = line.substr(line.find(":") + 1);
            solution_found = true;
        } else if (line.find("Dual objective value:") != std::string::npos) {
             dual_obj = line.substr(line.find(":") + 1);
        }
    }

    /* uncomment to print all info
    std::cout << "\nFull CSDP Output:\n" << result; // Uncomment to see everything CSDP prints
    */

    // Close the pipe and check the process exit status
    int status = pclose(pipe);
    if (status != 0) {
        std::cerr << "CSDP exited with non-zero status code: " << status << std::endl;
    }

    if (!solution_found) {
    	return -1.0;
    }

    /* uncomment to print essential info

    if (solution_found) {
        std::cout << "\nExternal CSDP Run Results:" << std::endl;
        std::cout << "Primal Objective: " << primal_obj;
        std::cout << "Dual Objective: " << dual_obj;
    } else {
        std::cout << "\nNo optimal solution found in external run output." << std::endl;
    }
    */
    

    return -std::stod(primal_obj);
}

#endif