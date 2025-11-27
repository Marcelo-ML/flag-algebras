#ifndef LOCALCUTSDPGENERATOR_H_INCLUDED
#define LOCALCUTSDPGENERATOR_H_INCLUDED

#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <iomanip>
#include <map>

namespace LocalCut {

    // Represents a matrix coefficient for the A_k variables in a constraint
    struct MatrixCoeff {
        unsigned long matrixIndex; // Which matrix variable A_1...A_s? (0-based)
        unsigned long row;         // Row in A_k
        unsigned long col;         // Col in A_k
        double value;              // Coefficient d_{ij}
    };

    class SDPToSDPA {
    private:
        unsigned long numAlphas;
        std::vector<unsigned long> matrixDimensions; // Dimensions of A_1, A_2, ...
        
        // Store constraints
        struct Constraint {
            double b; // RHS constant
            std::vector<double> c_coeffs; // Coefficients for alphas
            std::vector<MatrixCoeff> d_coeffs; // Coefficients for A matrices
        };
        std::vector<Constraint> constraints;

    public:
        SDPToSDPA(unsigned long);

        void addMatrixVariable(unsigned long);

        void addConstraint(double, const std::vector<double>&, const std::vector<MatrixCoeff>&);

        void writeToFile(const std::string&);

    private:
        void writeEntry(std::ofstream& file, unsigned long, unsigned long, unsigned long, unsigned long, double);
    };
}

#endif