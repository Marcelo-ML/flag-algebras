#ifndef LOCALCUTCOEFF_H_INCLUDED
#define LOCALCUTCOEFF_H_INCLUDED


#include "FlagGraph.h"
#include "Flag.h"
#include "bignum.h"

#include <vector>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <iostream>


namespace LocalCut {
    /* get [F] F_prime */
    BigNum::frac getCoefficient (
        const FlagAlgebra::Graph& F,
        const FlagAlgebra::Graph& F_prime
        );

    /* get [F] [[F_prime]] */
    BigNum::frac getCoefficient (
        const FlagAlgebra::Graph& F,
        const FlagAlgebra::Flag<FlagAlgebra::Graph>& F_prime
        );

    /* get the coefficient of F in [[F1 * F2]] */
    BigNum::frac getCoefficientInProduct (
        const FlagAlgebra::Graph& F,
        const FlagAlgebra::Flag<FlagAlgebra::Graph>& F1,
        const FlagAlgebra::Flag<FlagAlgebra::Graph>& F2
        );
	
    BigNum::frac calculateLocalCutCoefficient (
        const FlagAlgebra::Graph& F,                // o flag F
        const FlagAlgebra::Graph& sigma,            // o tipo sigma 
        const std::vector<BigNum::frac>& p  // p(X) para 2^k subconjuntos X
        );

}

bool hasTriangle ( const FlagAlgebra::Graph& G );

#endif