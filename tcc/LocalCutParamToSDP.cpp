#ifndef LOCALCUTPARAMTOSDP_CPP_INCLUDED
#define LOCALCUTPARAMTOSDP_CPP_INCLUDED

#include "LocalCutParamToSDP.h"

using namespace FlagAlgebra;
using namespace BigNum;
using std::vector;

namespace LocalCut {

ParamToSDP::ParamToSDP () { tmpFilename = "tmp.dat-s"; }

void ParamToSDP::setTmpFilename (const std::string& filename) { tmpFilename = std::string(filename); }

void ParamToSDP::addLCCondition (LocalCutCondition cond) {
	lcConditions.push_back(cond);
}

/* reads conditions from a file, each condition in the following format:
 * - Graph G as in Graph::readFromFile
 * - probability values assigned to the graph, as a bitmap
 * - their values p[0], p[1], ..., p[2^|G|-1]
 * For example:
 (3
 0 1
 0 2)
 0.0 0.0 1.0 1.0 0.5 0.5 0.5 1.0
 */

void ParamToSDP::addLCCondition (const std::string& filename) {

	FILE* file = std::fopen(filename.c_str(), "r");
	
	if (file == nullptr) {
		std::cerr << "Error opening the file" << std::endl;
		return;
	}


	Graph G;
	while (G.readFromFile(file) == 1) {
		vector<BigNum::frac> p(1l<<(G.size()),BigNum::frac(bn0));
		double x;
		for (unsigned long i=0; i<p.size(); ++i) {
			if (!std::fscanf(file, "%lf", &x)) {
				break;
			}
			p[i] = doubleToFrac(x);
		}
		addLCCondition({Graph(G),p});
	}

    std::fclose(file);

}

void ParamToSDP::addSDCondition (SDCondition cond) {
	sdConditions.push_back(cond);
}

void ParamToSDP::addSDCondition (const std::string& filename) {
	FILE* file = std::fopen(filename.c_str(), "r");
	
	if (file == nullptr) {
		std::cerr << "Error opening the file" << std::endl;
		return;
	}

	Graph G;
	unsigned long m;
	while (G.readFromFile(file) == 1) {
		if (std::fscanf(file, "%lu", &m)) {
			addSDCondition({G,m});
		}
	}

    std::fclose(file);
}

/* this is the most important function
 * it returns the SDPToSDPA with the conditions transformed into a solvable program
 */
SDPToSDPA ParamToSDP::getSDP() {
	unsigned long expSize = 0; // the minimum size of the expansion needed

	unsigned long r = lcConditions.size();
	unsigned long s = sdConditions.size();

	for (unsigned long i=0; i<r; ++i) {
		expSize = std::max(expSize, lcConditions[i].sigma.size() + 2);
	}

	for (unsigned long i=0; i<s; ++i) {
		expSize = std::max(expSize, 2 * sdConditions[i].m - sdConditions[i].pi.size());
	}

	std::cout << "Expansion in size " << expSize << std::endl;


	// STEP 0 - construct the sdp
	SDPToSDPA sdp(r);

	// STEP 1 - expand the flags from semidefinite conditions
	vector< container< Flag<Graph> > > enumFlags(s);
	for (unsigned long i=0; i<s; ++i) {
		morphism f(sdConditions[i].pi.size());
		for (unsigned long j=0; j<f.size(); ++j)
			f[j]=j;
		Flag<Graph> pi_flag(sdConditions[i].pi,f);

		enumFlags[i] = enumerateFlags<Graph>(sdConditions[i].m,pi_flag);
	}

	// STEP 1.1 - clear out the flags with triangles
	vector< container< Flag<Graph> > > enumFlagsNoTriangle(s);
	for (unsigned long i=0; i<s; ++i) {
		for (unsigned long j=0; j<enumFlags[i].size(); ++j) {
			if (!hasTriangle(enumFlags[i][j].object())) {
				enumFlagsNoTriangle[i].push_back(enumFlags[i][j]);
			}
		}
	}
	enumFlags = enumFlagsNoTriangle;

	// STEP 1.2 - the matrix sizes to the program
	for (unsigned long i=0; i<s; i++) {
		sdp.addMatrixVariable(enumFlags[i].size());
	}

	// STEP 2 - enumerate the graphs of size m
	container<Graph> objs_m = enumerateObjects<Graph>(expSize);
	
	vector<vbool> mat(2,vbool(2,false));
	mat[0][1] = mat[1][0] = true;
	Graph edge(mat);

	// STEP 3 - traverse the flags and add the constraints
	for (unsigned long i=0; i<objs_m.size(); ++i) {
		const Graph& F = objs_m[i];
		if (hasTriangle(F)) continue;

		// STEP 3.1 - determine b
		double b = fracToDouble(getCoefficient(F,edge));

		// STEP 3.2 - determine c
		vector<double> c(r,0.0);
		for (unsigned i=0; i<r; ++i) {
			BigNum::frac c_frac = 
				calculateLocalCutCoefficient(
					F,
					lcConditions[i].sigma,
					lcConditions[i].p
					);
			c[i] = fracToDouble(c_frac);
		}

		vector<MatrixCoeff> d(0);
		// STEP 3.3 - determine d
		for (unsigned matIdx=0; matIdx<s; ++matIdx) {
			unsigned long ell_k = enumFlags[matIdx].size();
			for (unsigned row=0; row<ell_k; ++row) {
				for (unsigned col=row; col<ell_k; ++col) {
					BigNum::frac d_frac = getCoefficientInProduct(
						F,
						enumFlags[matIdx][row],
						enumFlags[matIdx][col]
						);
					if (!d_frac.zero())
						d.push_back({matIdx, row+1, col+1, fracToDouble(d_frac)});
				}
			}
		}

		// STEP 3.4 - add the restriction
		sdp.addConstraint(b,c,d);
	}

	// STEP 4
	return sdp;
}

double ParamToSDP::getOptimalValue () {
	getSDP().writeToFile(tmpFilename);
	double opt_val = solve_sdpa_file(tmpFilename);
	return opt_val;
}

} // namespace LocalCut
#endif