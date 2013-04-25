#include <Rcpp.h>
#include "shotgun/common.h"
#include "multitask.h"

SEXP multitask(SEXP X0, SEXP y0, SEXP tasks0, SEXP groups0)
{
	/*convert parameters to Rcpp types*/
	Rcpp::NumericMatrix A(X0); /*using name A to be consistent with Shotgun library*/
	Rcpp::NumericVector y(y0);
	Rcpp::CharacterVector tasks(tasks0);
	Rcpp::NumericMatrix groups(groups0);
	Rcpp::NumericMatrix result(X0);
	
	shotgun_data data;

	/*init data.A_cols and data.A_rows*/
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

	double lambda = 1.0;
	int regpathLength = 0;
	double threshold = 1.0e-5;
	int maxIter = 100;
	int verbose = 1;

	solveLasso(&data, lambda, regpathLength, threshold, maxIter, verbose);
	
	return result;
}
