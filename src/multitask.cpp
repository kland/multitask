#include <Rcpp.h>
#include "shotgun/common.h"
#include "multitask.h"

SEXP multitask_lasso(SEXP X0, SEXP y0, SEXP lambda0, SEXP eps0)
{
	//convert parameters to Rcpp types
	Rcpp::NumericMatrix X(X0);
	Rcpp::NumericVector y(y0);
	Rcpp::NumericVector lambda(lambda0); //one-element vector
	Rcpp::NumericVector eps(eps0); //one-element vector
	
	shotgun_data data; //input to and output from function solveLasso

	//initialize input fields (shotgun calls the design matrix A instead of X)
	int M = X.nrow();
	int N = X.ncol();
	data.A_cols.resize(N);
	data.A_rows.resize(M);
	for (int k = 0; k < M * N; k++) {
		if (X[k] != 0.0) {
			int i = k / N; //row
			int j = k % N; //column
			data.A_cols.at(j).add(i, X[k]);
			data.A_rows.at(i).add(j, X[k]);
		}
	}
	data.y = Rcpp::as< std::vector<valuetype_t> >(y);
	data.nx = N;
	data.ny = M;

	int regpathLength = 0;
	int maxIter = 100;
	int verbose = 0;

	solveLasso(&data, lambda[0], regpathLength, eps[0], maxIter, verbose);
	
	return Rcpp::wrap(data.x);
}
