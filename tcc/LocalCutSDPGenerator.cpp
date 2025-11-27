#ifndef LOCALCUTSDPGENERATOR_CPP_INCLUDED
#define LOCALCUTSDPGENERATOR_CPP_INCLUDED

#include "LocalCutSDPGenerator.h"

/**
 * A helper class to formulate a specific Semidefinite Program (SDP)
 * and export it to the SDPA sparse format (.dat-s) for solvers like CSDP.
 *
 * Problem Formulation:
 * Minimize: M - 0.08 * sum(alpha_i)
 * Subject to:
 * M >= b_k + sum(c_ki * alpha_i) + sum(trace(D_k * A))  (for each constraint k)
 * alpha >= 0
 * A_j >= 0 (PSD)
 *
 * Mapping to CSDP Primal Form (Max tr(CX) s.t. A(X)=b):
 * We Maximize: -M + 0.08 * sum(alpha_i)
 * Constraints: M - sum(c_ki * alpha_i) - sum(trace(D_k * A)) - s_k = b_k
 * Where s_k >= 0 are slack variables.
 */



LocalCut::SDPToSDPA::SDPToSDPA(unsigned long r_alphas) : numAlphas(r_alphas) {}


// Add a matrix variable A_k with dimension size x size
void LocalCut::SDPToSDPA::addMatrixVariable(unsigned long size) {
    matrixDimensions.push_back(size);
}

/**
 * Add a constraint: 
 * M >= b + sum(c_i * alpha_i) + sum(d_coeffs * A)
 * * Internally converts to:
 * M - sum(c_i * alpha_i) - sum(d * A) - slack = b
 */
void LocalCut::SDPToSDPA::addConstraint(double b, const std::vector<double>& c_coeffs, const std::vector<MatrixCoeff>& d_coeffs) {
    if (c_coeffs.size() != numAlphas) {
        std::cerr << "Error: Size of c_coeffs must match number of alphas.\n";
        return;
    }
    constraints.push_back({b, c_coeffs, d_coeffs});
}

// Write the problem to a .dat-s file
void LocalCut::SDPToSDPA::writeToFile(const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << "\n";
        return;
    }
    file << std::fixed << std::scientific << std::setprecision(15);

    unsigned long numConstraints = constraints.size();
    
    // --- 1. Header Information ---
    file << numConstraints << "\n"; // m (number of constraints)
    
    // Calculate number of blocks
    // Block 1: M (represented as M+ and M-, size 2 diag)
    // Block 2: Alpha (size numAlphas diag)
    // Blocks 3 to 2+s: Matrix variables A_k
    // Block Last: Slacks (size numConstraints diag)
    
    // handle the edge case where there are no alpha conditions
    bool hasAlphas = (numAlphas > 0);

    unsigned long numBlocks = 1 + hasAlphas + matrixDimensions.size() + 1; 
    file << numBlocks << "\n";

    // --- 2. Block Structure ---
    // Negative numbers indicate diagonal blocks, positive indicate dense blocks
    
    file << "-1 "; // Block 1: (Diagonal, size 1 since M >= 0)

    if (hasAlphas)
        file << "-" << numAlphas << " "; // Block 2: Alphas (Diagonal)
    
    // Matrix blocks
    for (unsigned long dim : matrixDimensions) {
        file << dim << " "; // Dense blocks
    }
    
    file << "-" << numConstraints << "\n"; // Slacks (Diagonal)

    // --- 3. RHS Vector (b) ---
    // CSDP 'c' vector is the RHS of the constraints
    for (const auto& cons : constraints) {
        file << cons.b << " ";
    }
    file << "\n";

    // --- 4. Matrix Definitions (Sparse Format) ---
    // Format: <Mat Index> <Block Index> <Row> <Col> <Value>
    // Mat Index 0 is Objective, 1..m are constraints
    

    // === MATRIX 0: OBJECTIVE FUNCTION ===
    // We want to Minimize: M - 0.08 * sum(alpha)
    // CSDP Maximizes tr(C*X). 
    // Maximize: -(M - 0.08*sum(alpha)) = -M + 0.08*sum(alpha)
    
    // Block 1 (M): Coeff is -1
    writeEntry(file, 0, 1, 1, 1, -1.0);

    // Block 2 (Alpha): Coeffs are +0.08
    for (unsigned long i = 1; i <= numAlphas; ++i) {
        writeEntry(file, 0, 2, i, i, 0.08);
    }
    // Other blocks (A_k, Slacks) have 0 coeff in objective

    // === MATRICES 1..m: CONSTRAINTS ===
    // Constraint k: M - sum(c * alpha) - sum(d * A) - slack_k = b
    
    for (unsigned long k = 0; k < numConstraints; ++k) {
        unsigned long matIdx = k + 1; // 1-based index for constraints

        // Block 1 (M): Coeff is 1
        writeEntry(file, matIdx, 1, 1, 1, 1.0);

        // Block 2 (Alpha): Coeffs -c_i
        const auto& c = constraints[k].c_coeffs;
        for (unsigned long i = 0; i < numAlphas; ++i) {
            if (c[i] != 0.0) {
                writeEntry(file, matIdx, 2, i + 1, i + 1, -c[i]);
            }
        }

        // Blocks 3.. (Matrices A_k): Coeffs -d
        // Note: Block index for A_j is 3 + j
        const auto& d = constraints[k].d_coeffs;
        for (const auto& coeff : d) {
            unsigned long blockIdx = 2 + hasAlphas + coeff.matrixIndex;
            // SDPA sparse format requires row <= col usually, but CSDP handles symmetric expansion.
            // It is safer to ensure we write the entry. 
            // If it's a dense block, standard is row, col.
            writeEntry(file, matIdx, blockIdx, coeff.row, coeff.col, -coeff.value);
        }

        // Block Last (Slacks): Coeff -1 for the k-th slack
        unsigned long slackBlockIdx = numBlocks;
        writeEntry(file, matIdx, slackBlockIdx, k + 1, k + 1, -1.0);
    }

    file.close();
    std::cout << "Successfully wrote SDPA file: " << filename << "\n";
}


void LocalCut::SDPToSDPA::writeEntry(std::ofstream& file, unsigned long mat, unsigned long blk, unsigned long r, unsigned long c, double val) {
    file << mat << " " << blk << " " << r << " " << c << " " << val << "\n";
}

using namespace LocalCut;

void example_usage() {
    // === EXAMPLE USAGE ===
    
    // 1. Setup: 2 alpha variables
    unsigned long r = 2; 
    SDPToSDPA sdp(r);

    // 2. Setup: 2 Matrix variables
    // A_0 is 2x2, A_1 is 3x3
    sdp.addMatrixVariable(2); 
    sdp.addMatrixVariable(3);

    // 3. Define Constraints
    
    // Constraint 1: M >= 5 + (2*a1 + 3*a2) + (1.5 * (A_0)1,1)
    // Converted form for solver: M - (2*a1 + 3*a2) - (1.5 * (A_0)1,1) >= 5
    
    double b1 = 5.0;
    std::vector<double> c1 = {2.0, 3.0}; // coeffs for alpha_1, alpha_2
    std::vector<MatrixCoeff> d1;
    // (A_0)_1,1 has coeff 1.5. Indices are 1-based for the user convenience here?
    // Let's use 1-based row/col for clarity in CSDP, or map internally.
    // The class writeEntry uses the passed integers directly. 
    // Usually CSDP dense blocks use 1-based indexing.
    d1.push_back({0, 1, 1, 1.5}); // Matrix 0 (A_0), row 1, col 1, val 1.5

    sdp.addConstraint(b1, c1, d1);

    // Constraint 2: M >= 2 + (0.5*a1) + (4 * (A_1)2,3)
    double b2 = 2.0;
    std::vector<double> c2 = {0.5, 0.0}; // alpha_2 has 0 coeff
    std::vector<MatrixCoeff> d2;
    d2.push_back({1, 2, 3, 4.0}); // Matrix 1 (A_1), row 2, col 3, val 4.0

    sdp.addConstraint(b2, c2, d2);

    // 4. Generate File
    sdp.writeToFile("problem.dat-s");
}

#endif