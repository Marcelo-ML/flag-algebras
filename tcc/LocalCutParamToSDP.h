#ifndef LOCALCUTPARAMTOSDP_H_INCLUDED
#define LOCALCUTPARAMTOSDP_H_INCLUDED

#include "FlagGraph.h"
#include "Flag.h"
#include "bignum.h"

#include "LocalCutCoeff.h"
#include "LocalCutSDPGenerator.h"
#include "LocalCutSDPSolver.h"

#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <algorithm>
#include <cstdio>

namespace LocalCut {

	struct LocalCutCondition {
		const FlagAlgebra::Graph sigma;            // the type sigma
		const std::vector<BigNum::frac> p; // the probability attribution to each X subset of V(sigma)
                                       // represent
	};

	struct SDCondition {
		const FlagAlgebra::Graph pi;  // the type pi
		unsigned long m;  // m >= |pi| for the flag expansion before the semidefinite program
	};

	/* this class handles addition of local cut conditions and SDP creation */

	class ParamToSDP {
	public:
		ParamToSDP();

		void addLCCondition(LocalCutCondition);
		void addLCCondition(const std::string&);

		void addSDCondition(SDCondition);
		void addSDCondition(const std::string&);

		SDPToSDPA getSDP();
		double getOptimalValue();

	private:
		std::vector<LocalCutCondition> lcConditions;
		std::vector<SDCondition> sdConditions;
	};
}

#endif