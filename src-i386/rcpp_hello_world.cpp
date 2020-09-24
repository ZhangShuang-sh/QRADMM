#include <Rcpp.h>
#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <math.h>
#include <Rmath.h>
#include <stdio.h>
#include <stdlib.h>
//#include <Eigen/Dense>  

//using Eigen::MatrixXd;  
//using namespace Eigen;
//using namespace Eigen::internal;
//using namespace Eigen::Architecture;

using namespace Rcpp;
using namespace std;
// Generate data toy with dimension of 'dimnum' and sample account equal
// to 'samplenum'.
// In this module, R library 'MASS' is needed.
RcppExport SEXP generate_data(SEXP dimnum, SEXP samplenum, SEXP fun1, SEXP fun2, SEXP fun3) {
	Rcpp::Function mvrnorm("mvrnorm");
	Rcpp::Function rnorm("rnorm");
	Rcpp::Function normcdf("pnorm");
	Rcpp::Function matrixmulti(fun1);
	Rcpp::Function matrixrow(fun2);
	Rcpp::Function matrixtrans(fun3);
	//Function matrixinv(fun2);
	//Function matrixtrans(fun3);
	int samplenumber = Rcpp::as<int>(samplenum);
	int dimensionnumber = Rcpp::as<int>(dimnum);
	Rcpp::NumericVector sigma =
		Rcpp::NumericMatrix(Rcpp::Dimension(dimensionnumber, dimensionnumber));
	Rcpp::NumericVector x =
		Rcpp::NumericMatrix(Rcpp::Dimension(samplenumber, dimensionnumber));
	Rcpp::NumericVector beta =
		Rcpp::NumericMatrix(Rcpp::Dimension(samplenumber, dimensionnumber));
	Rcpp::NumericVector y =
		Rcpp::NumericMatrix(Rcpp::Dimension(samplenumber, 1));
	for (int i = 0; i < dimensionnumber; i++) {
		for (int j = 0; j < samplenumber; j++) {
			beta(j, i) = 0;
		}
	}
	for (int i = 0; i < samplenumber; i++) {
		beta(i, 5) = 1;
		beta(i, 11) = 1;
		beta(i, 14) = 1;
		beta(i, 19) = 1;
		beta(i, 0) = 0.7 * (Rcpp::as<double>(rnorm(1, 0, 1)));
	}
	for (int i = 0; i < dimensionnumber; i++) {
		for (int j = 0; j < dimensionnumber; j++) {
			sigma(i, j) = pow(0.5, abs(i - j));
		}
	}
	x = mvrnorm(samplenumber, rep(0, dimensionnumber), sigma);
	for (int i = 0; i < samplenumber; i++) {
		x(i, 0) = Rcpp::as < double >(normcdf(x(i, 0)));
	}
	for (int i = 0; i < samplenumber; i++) {
		y(i, 0) = Rcpp::as < double >(matrixmulti(matrixrow(x,i), matrixtrans(matrixrow(beta, i))));
	}
	int i = 2;
	return Rcpp::List::create(Rcpp::Named("X") = x, Rcpp::Named("Y") = y, Rcpp::Named("BETA")=beta, Rcpp::Named("CHECK")= matrixmulti(matrixrow(x, i), matrixtrans(matrixrow(beta, i))));
}

// Update beta
RcppExport SEXP Example_of_admm(SEXP fixedpart, SEXP exper_time, SEXP dimnum, SEXP samplenum, SEXP lambda, SEXP gamma, SEXP tau, SEXP z, SEXP r1, SEXP r2, SEXP q1, SEXP q2, SEXP w, SEXP X, SEXP Y, SEXP fun1, SEXP fun2, SEXP fun3, SEXP fun4, SEXP fun5, SEXP fun6)
{
	Rcpp::NumericMatrix Z = Rcpp::clone<Rcpp::NumericMatrix>(z);
	Rcpp::NumericMatrix FIXED = Rcpp::clone<Rcpp::NumericMatrix>(fixedpart);
	Rcpp::NumericMatrix Q1 = Rcpp::clone<Rcpp::NumericMatrix>(q1);
	Rcpp::NumericMatrix Q2 = Rcpp::clone<Rcpp::NumericMatrix>(q2);
	Rcpp::NumericMatrix W = Rcpp::clone<Rcpp::NumericMatrix>(w);
	Rcpp::NumericMatrix x = Rcpp::clone<Rcpp::NumericMatrix>(X);
	Rcpp::NumericMatrix y = Rcpp::clone<Rcpp::NumericMatrix>(Y);
	int dimensionnumber = Rcpp::as<int>(dimnum);
	int samplenumber = Rcpp::as<int>(samplenum);
	Rcpp::NumericVector beta = Rcpp::NumericMatrix(Rcpp::Dimension(dimensionnumber,1));
	//Rcpp::NumericVector y_temp = Rcpp::NumericMatrix(Rcpp::Dimension(1, 1));
	//Rcpp::NumericVector x_temp = Rcpp::NumericMatrix(Rcpp::Dimension(1, dimensionnumber));
	double temp;
	double temp2;
	Rcpp::Function matrixmulti(fun1);
	Rcpp::Function matrixplus(fun2);
	Rcpp::Function matrixminu(fun3);
	Rcpp::Function matrixmultinumber(fun4);
	Rcpp::Function matrixtrans(fun5);
	Rcpp::Function matrixrow(fun6);
	double R1 = Rcpp::as<double>(r1);
	double R2 = Rcpp::as<double>(r2);
	double TAU = Rcpp::as<double>(tau);
	double GAMMA = Rcpp::as<double>(gamma);
	double LAMBDA = Rcpp::as<double>(lambda);
	int times = Rcpp::as<int>(exper_time);
	for (int s = 0; s < times; s++) {
		beta = matrixmulti(FIXED, matrixminu(matrixplus(matrixplus(matrixmulti(matrixmultinumber((R1 / R2), matrixtrans(x)), matrixminu(y, Z)), W), matrixmulti(matrixmultinumber((1 / R2), matrixtrans(x)), Q1)), matrixmultinumber((1 / R2), Q2)));
		for (int i = 0; i < samplenumber; i++) {
			temp = y(i, 0) - Rcpp::as<double>(matrixmulti(matrixrow(x, i), beta)) + Q1(i, 0) / R1;
			if (TAU / R1 < temp) Z(i, 0) = temp - TAU / R1;
			else if (temp < (TAU - 1) / R1) Z(i, 0) = temp - (TAU - 1) / R1;
			else Z(i, 0) = 0;
		}
		for (int j = 0; j < dimensionnumber; j++) {
			temp2 = beta(j,0) + (1 / R2)*Q2(j,0);
			if (GAMMA*LAMBDA < temp2 || GAMMA*LAMBDA < (-1)*temp2) W(j,0) = temp2;
			else if (temp2 <= LAMBDA*GAMMA && temp2 >(1 + (1 / R2))*LAMBDA) W(j, 0) = (R2*(GAMMA - 1)*temp2 - GAMMA*LAMBDA) / (R2*(GAMMA - 1) - 1);
			else if ((-1)*temp2 <= LAMBDA*GAMMA && (-1)*temp2>(1 + (1 / R2))*LAMBDA) W(j, 0) = (R2*(GAMMA - 1)*temp2 + GAMMA*LAMBDA) / (R2*(GAMMA - 1) - 1);
			else if ((temp2 <= ((1 + 1 / R2)*LAMBDA)) && (temp2 > LAMBDA / R2)) W(j, 0) = temp2 - LAMBDA / R2;
			else if (((-1)*temp2 <= ((1 + 1 / R2)*LAMBDA)) && ((-1)*temp2 > LAMBDA / R2)) W(j, 0) = temp2 + LAMBDA / R2;
			else W(j, 0) = 0;
		}
		Q1 = matrixplus(Q1, matrixmultinumber(R1, matrixminu(matrixminu(y, matrixmulti(x, beta)), Z)));
		Q2 = matrixplus(Q2, matrixmultinumber(R2, matrixminu(beta, W)));
	}
	return Rcpp::List::create(Rcpp::Named("beta_result") = beta, Rcpp::Named("z_result") = Z, Rcpp::Named("w_result") = W, Rcpp::Named("q1_result") = Q1, Rcpp::Named("q2_result") = Q2);
}
