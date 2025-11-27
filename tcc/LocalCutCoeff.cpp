#ifndef LOCALCUTCOEFF_CPP_INCLUDED
#define LOCALCUTCOEFF_CPP_INCLUDED


#include "LocalCutCoeff.h"

using namespace FlagAlgebra;
using namespace LocalCut;
using namespace BigNum;
using std::vector;

#include <iostream>
#include <algorithm>

BigNum::frac LocalCut::getCoefficient (
    const Graph& F,
    const Graph& F_prime
    )
{
    unsigned long expSize = F.size();
    unsigned long baseSize = F_prime.size();

    const container<Graph>& obj_expSize = enumerateObjects<Graph>(expSize);
    const container<Graph>& obj_baseSize = enumerateObjects<Graph>(baseSize);

    unsigned long F_idx = 
        std::find(
            obj_expSize.begin(),
            obj_expSize.end(),
            F
        ) - obj_expSize.begin();
    unsigned long F_prime_idx =
        std::find(
            obj_baseSize.begin(),
            obj_baseSize.end(),
            F_prime
        ) - obj_baseSize.begin();

    const container< container<coeff_type> >& exp_matrix_obj = expandObjectEnumerate<Graph>(baseSize, expSize);

    for (unsigned long i=0; i<exp_matrix_obj[F_prime_idx].size(); ++i) {
        FlagAlgebra::coeff_type coeff_pair = exp_matrix_obj[F_prime_idx][i];
        if (coeff_pair.first == F_idx) {
            return coeff_pair.second;
        }
    }

    return BigNum::frac(bn0);
}

BigNum::frac LocalCut::getCoefficient (
    const Graph& F,
    const Flag<Graph>& F_prime
    )
{
    return F_prime.normFactor() * getCoefficient(F,F_prime.object());
}

BigNum::frac LocalCut::getCoefficientInProduct (
    const Graph& F,
    const Flag<Graph>& F1,
    const Flag<Graph>& F2
    )
{
    const unsigned long expSize = F.size();
    const unsigned long base1Size = F1.size();
    const unsigned long base2Size = F2.size();

    /* assuming F1 and F2 have the same type and the same typeMorphism */

    const Graph& unlabeledType = F1.unlabeledType();
    
    /* we have to copy the morphism */
    morphism f (F1.typeMorphism().size());
    for (unsigned long i=0; i<F1.typeMorphism().size(); ++i) {
        f[i] = (F1.typeMorphism())[i];
    }

    const container< container< container<coeff_type> > >& exp_prod_matrix =
        expandFlagProductEnumerate<Graph>(
            unlabeledType,
            base1Size,
            base2Size,
            expSize,
            f
        );

    const container< Flag<Graph> >& obj_base1Size = enumerateFlags<Graph>(base1Size,unlabeledType,f);
    const container< Flag<Graph> >& obj_base2Size = enumerateFlags<Graph>(base2Size,unlabeledType,f);
    const container< Flag<Graph> >& obj_expSize = enumerateFlags<Graph>(expSize,unlabeledType,f);

    unsigned long F1_idx = 
        std::find(
            obj_base1Size.begin(),
            obj_base1Size.end(),
            F1
        ) - obj_base1Size.begin();
    unsigned long F2_idx = 
        std::find(
            obj_base2Size.begin(),
            obj_base2Size.end(),
            F2
        ) - obj_base2Size.begin();

    BigNum::frac total_coefficient = BigNum::frac(bn0);

    for (unsigned long i=0; i<exp_prod_matrix[F1_idx][F2_idx].size(); ++i) {
        FlagAlgebra::coeff_type coeff_pair = exp_prod_matrix[F1_idx][F2_idx][i];
        const Flag<Graph>& F_prime = obj_expSize[coeff_pair.first];
        if (hasTriangle(F_prime.object())) continue;
        total_coefficient += coeff_pair.second * getCoefficient(F,F_prime);
    }

    return total_coefficient;
}


BigNum::frac LocalCut::calculateLocalCutCoefficient (
    const Graph& F,                // o flag F
    const Graph& sigma,            // o tipo sigma 
    const vector<BigNum::frac>& p  // p(X) para 2^k subconjuntos X (como mask)
    )
{
    const unsigned long k = sigma.size();
    const unsigned long k_prime = k+2;

    BigNum::frac total_coefficient(bn0);

    for (unsigned long idx_X = 0; idx_X < (1UL << k); ++idx_X) {
        for (unsigned long idx_Y = idx_X; idx_Y < (1UL << k); ++idx_Y){

            BigNum::frac p_X = p[idx_X], p_Y = p[idx_Y];
            BigNum::frac one(bn1);
            BigNum::frac p_coeff = p_X*p_Y + (one-p_X)*(one-p_Y);
            if (p_coeff.zero()) continue;

            /* build F_XY */

            vector<vbool> F_XY_adj_mat(k_prime,vbool(k_prime,false));

            /* add the sigma edges */
            for (unsigned long i=0; i<k; ++i) {
                for (unsigned long j=0; j<k; ++j) {
                    if (sigma.is_adj(i,j)) {
                        F_XY_adj_mat[i][j] = true;
                    }
                }
            }

            /* add the X and Y edges */
            for (unsigned long i_sig = 0; i_sig<k; ++i_sig) {
                if((idx_X >> i_sig) & 1) F_XY_adj_mat[k][i_sig] = F_XY_adj_mat[i_sig][k] = true;
                if((idx_Y >> i_sig) & 1) F_XY_adj_mat[k+1][i_sig] = F_XY_adj_mat[i_sig][k+1] = true;
            }

            /* add uv edge */
            F_XY_adj_mat[k][k+1] = F_XY_adj_mat[k+1][k] = true;

            /* check if the graph has triangles */
            if (hasTriangle(Graph(F_XY_adj_mat))) continue;

            /* build the morphism */
            morphism theta_F_XY(k);
            for (unsigned long i=0; i<k; ++i) theta_F_XY[i]=i;

            /* construct the flag */
            Flag<Graph> F_XY_flag(Graph(F_XY_adj_mat), theta_F_XY);

            total_coefficient += p_coeff * getCoefficient(F,F_XY_flag);
        }
    }

    return total_coefficient;
}

bool hasTriangle( const Graph& G ) {
    unsigned long n = G.size();

    for (unsigned long a=0;   a<n; ++a)
    for (unsigned long b=a+1; b<n; ++b)
    for (unsigned long c=b+1; c<n; ++c)
        if (G.is_adj(a,b) && G.is_adj(b,c) && G.is_adj(c,a))
            return true;

    return false;
}

#endif