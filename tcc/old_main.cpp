#include "LocalCutCoeff.h"
#include "LocalCutSDPGenerator.h"
#include "LocalCutSDPSolver.h"

#include <iostream>
#include <vector>

using namespace LocalCut;
using namespace FlagAlgebra;
using namespace BigNum;
using std::vector;
using std::cout;
using std::endl;
using std::cerr;

void test_basic_cut(){
	vector<vbool> mat(4,vbool(4,false));
	mat[0][1] = mat[1][0] = mat[0][2] = mat[2][0] = true;
	Graph F(mat);

	mat = vector<vbool>(1,vbool(1,false));
	Graph sigma(mat);

	vector<BigNum::frac> p(2); p[0]=BigNum::frac(bn0); p[1]=BigNum::frac(bn1);

	BigNum::frac coef = calculateLocalCutCoefficient(F,sigma,p);

	cout << coef << endl;
}

void test_bcl_cut(){
	vector<vbool> mat(5,vbool(5,false));
	mat[0][1] = mat[1][0] = mat[1][2] = mat[2][1] = mat[2][3] = mat[3][2] = mat[3][4] = mat[4][3] = true;
	Graph F(mat);

	mat = vector<vbool>(2,vbool(2,false));
	mat[0][1]=mat[1][0] = true;
	Graph sigma(mat); // a labeled edge

	vector<BigNum::frac> p(4); p[0]=p[1]=p[2]=p[3]=BigNum::frac(bn1,bn2); // p is (1/2,1/2,1/2,1/2)

	BigNum::frac coef = calculateLocalCutCoefficient(F,sigma,p);

	cout << coef << endl;
}

void test_get_coefficient(){
	vector<vbool> mat(3,vbool(3,false));
	mat[0][1] = mat[1][0] = true;
	morphism f(0);
	Flag<Graph> F_prime(Graph(mat),f); // co-cherry

	mat = vector<vbool>(5,vbool(5,false));
	mat[0][1] = mat[1][0] = mat[0][2] = mat[2][0] = mat[1][2] = mat[2][1] = mat[0][4] = mat[4][0] = true;
	Graph F(mat);

	cout << getCoefficient(F,F_prime) << endl;
}

void test_get_coefficient_in_product(){ // answer is 1/90
	vector<vbool> mat(4,vbool(4,false));
	morphism f(2); f[0]=0; f[1]=1;

	mat[0][1] = mat[1][0] = mat[0][2] = mat[2][0] = true;
	Flag<Graph> F1(Graph(mat),f);

	mat[2][3] = mat[3][2] = true;
	Flag<Graph> F2(Graph(mat),f);

	mat = vector<vbool>(6,vbool(6,false));
	mat[0][1] = mat[1][0] = true;
	mat[0][5] = mat[5][0] = true;
	mat[1][4] = mat[4][1] = true;
	mat[2][5] = mat[5][2] = true;
	mat[3][4] = mat[4][3] = true;
	mat[4][5] = mat[5][4] = true;

	Graph F(mat);

	BigNum::frac coef = getCoefficientInProduct (F, F1, F2);

	cout << coef << endl;
}

void test_sdp_mantel(){
	SDPToSDPA sdp(0);
	sdp.addMatrixVariable(2);

	double b = 2.0/3;
	vector<double> c = {};
	vector<MatrixCoeff> d;
	d.push_back({0,1,1,1.0/3});
	d.push_back({0,1,2,1.0/3});

	sdp.addConstraint(b,c,d);

	b = 1.0/3;
	c = {};
	d = {};
	d.push_back({0,1,2,1.0/3});
	d.push_back({0,2,2,1.0/3});

	sdp.addConstraint(b,c,d);

	b = 0;
	c = {};
	d = {};
	d.push_back({0,2,2,1.0});

	sdp.addConstraint(b,c,d);


	sdp.writeToFile("mantel.dat-s");

	double opt_val = solve_sdpa_file("mantel.dat-s");

	if (opt_val < 0) {
		cerr << "Optimal value is negative\n";
	}
	else {
		cout << "Optimal value is: " << opt_val << "\n";
	}
}

void test_sdp_efps(){
	SDPToSDPA sdp(1);
	sdp.addMatrixVariable(2);

	double b = 2.0/3;
	vector<double> c = {0};
	vector<MatrixCoeff> d;
	d.push_back({0,1,1,1.0/3});
	d.push_back({0,1,2,1.0/3});

	sdp.addConstraint(b,c,d);

	b = 1.0/3;
	c = {1.0/3};
	d = {};
	d.push_back({0,1,2,1.0/3});
	d.push_back({0,2,2,1.0/3});

	sdp.addConstraint(b,c,d);

	b = 0;
	c = {0};
	d = {};
	d.push_back({0,2,2,1.0});

	sdp.addConstraint(b,c,d);


	sdp.writeToFile("efps.dat-s");

	double opt_val = solve_sdpa_file("efps.dat-s");

	if (opt_val < 0) {
		cerr << "Optimal value is negative: " << opt_val << "\n";
	}
	else {
		cout << "Optimal value is: " << opt_val << "\n";
	}
}

int old_main(){
	test_sdp_efps();

	return 0;
}