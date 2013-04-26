#include <Rcpp.h>
#include "shotgun/common.h"
#include "multitask.h"

SEXP multitask_lasso(SEXP X0, SEXP y0, SEXP lambda0, SEXP eps0)
{
	//convert parameters to Rcpp types
	Rcpp::NumericMatrix A(X0); //using name A to be consistent with Shotgun library
	Rcpp::NumericVector y(y0);
	Rcpp::NumericVector lambda(lambda0); //one-element vector
	Rcpp::NumericVector eps(eps0); //one-element vector
	Rcpp::NumericMatrix result(X0);
	
	shotgun_data data;

	//init data.A_cols and data.A_rows
	int M = A.nrow();
	int N = A.ncol();
	data.A_cols.resize(N);
	data.A_rows.resize(M);
	for (int k = 0; k < M * N; k++) {
		if (A[k] != 0.0) {
			int i = k / N; //row
			int j = k % N; //column
			data.A_cols.at(j).add(i, A[k]);
			data.A_rows.at(i).add(j, A[k]);
		}
	}

	data.y = Rcpp::as< std::vector<valuetype_t> >(y);

	data.nx = N;
	data.ny = M;

	int regpathLength = 0;
	int maxIter = 100;
	int verbose = 1;

	solveLasso(&data, lambda[0], regpathLength, eps[0], maxIter, verbose);
	
	return result;
}
