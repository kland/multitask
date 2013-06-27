#include <Rcpp.h>
#include "shotgun/common.h"
#include "multitask.h"

static inline int &elem(Rcpp::IntegerMatrix &A, int i, int j) //returns element at position (i, j) in column-major order matrix A
{
	return A[A.nrow() * j + i];
}


static inline double &elem(Rcpp::NumericMatrix &A, int i, int j) //returns element at position (i, j) in column-major order matrix A
{
	return A[A.nrow() * j + i];
}


SEXP multitask_x_tilde(SEXP X0, SEXP tasks0, SEXP groups0, SEXP d_cur0, SEXP eta_cur0, SEXP K0)
{
	//convert input parameters to Rcpp types or primitive C++ types
	Rcpp::NumericMatrix X(X0);
	Rcpp::CharacterVector tasks(tasks0);
	Rcpp::IntegerMatrix groups(groups0);
	Rcpp::NumericVector d_cur(d_cur0);
	Rcpp::NumericMatrix eta_cur(eta_cur0);
	int K = Rcpp::as<int>(K0);
	assert(K > 0);
	
	int n = X.nrow() / K;
	int p = X.ncol();
	int L = groups.ncol();
	Rcpp::NumericMatrix result(n * K, p * K);

	for (int k = 0; k < K; k++) {	
		for (int j = 0; j < p; j++) {
			//calculate sum for column j
			double sum = 0.0;
			for (int l = 0; l < L; l++) {
				if (elem(groups, j, l)) {
					sum += d_cur[l] * elem(eta_cur, l, k);
				} 
			}
		
			//multiply column j in submatrix k of X with sum
			for (int i = 0; i < n; i++) {
				elem(result, n * k + i, p * k + j) = elem(X, n * k + i, j) * sum;
			}
		}
	}
	return result;
}


SEXP multitask_x_tilde_2(SEXP X0, SEXP tasks0, SEXP groups0, SEXP alpha_new0, SEXP eta_cur0, SEXP K0)
{
	//convert input parameters to Rcpp types or primitive C++ types
	Rcpp::NumericMatrix X(X0);
	Rcpp::CharacterVector tasks(tasks0);
	Rcpp::IntegerMatrix groups(groups0);
	Rcpp::NumericMatrix alpha_new(alpha_new0);
	Rcpp::NumericMatrix eta_cur(eta_cur0);
	int K = Rcpp::as<int>(K0);
	assert(K > 0);
	
	int n = X.nrow() / K;
	int p = X.ncol();
	int L = groups.ncol();
	Rcpp::NumericMatrix result(X.nrow(), L);

	for (int l = 0; l < L; l++) {
		for (int i = 0; i < n * K; i++) {
			int k = i / n;
			double sum = 0.0;
			for (int j = 0; j < p; j++) {
				if (elem(groups, j, l)) {
					sum += elem(X, i, j) * elem(alpha_new, j, k);
				}
			}
			elem(result, i, l) = elem(eta_cur, l, k) * sum;
		}
	}

	return result;
}


SEXP multitask_x_tilde_3(SEXP X0, SEXP tasks0, SEXP groups0, SEXP alpha_new0, SEXP d_new0, SEXP K0, SEXP k0)
{
	//convert input parameters to Rcpp types or primitive C++ types
	Rcpp::NumericMatrix X(X0);
	Rcpp::CharacterVector tasks(tasks0);
	Rcpp::IntegerMatrix groups(groups0);
	Rcpp::NumericMatrix alpha_new(alpha_new0);
	Rcpp::NumericVector d_new(d_new0);
	int K = Rcpp::as<int>(K0);
	assert(K > 0);
	int k = Rcpp::as<int>(k0) - 1; //C++ arrays are zero-based
	assert(k >= 0);

	int n = X.nrow() / K;
	int p = X.ncol();
	int L = groups.ncol();
	Rcpp::NumericMatrix result(n, L);
	
	for (int l = 0; l < L; l++) {
		for (int i = 0; i < n; i++) {
			double sum = 0.0;
			for (int j = 0; j < p; j++) {
				if (elem(groups, j, l)) {
					sum += elem(X, n * k + i, j) * elem(alpha_new, j, k);
				}
			}
			elem(result, i, l) = d_new[l] * sum;
		}
	}
	
	return result;
}


SEXP multitask_beta_new(SEXP groups0, SEXP alpha_new0, SEXP d_new0, SEXP eta_new0, SEXP K0)
{
	//convert input parameters to Rcpp types or primitive C++ types
	Rcpp::IntegerMatrix groups(groups0);
	Rcpp::NumericMatrix alpha_new(alpha_new0);
	Rcpp::NumericVector d_new(d_new0);
	Rcpp::NumericMatrix eta_new(eta_new0);
	int K = Rcpp::as<int>(K0);
	assert(K > 0);
	
	int p = groups.nrow();
	int L = groups.ncol();
	Rcpp::NumericMatrix result(p, K);
	
	for (int k = 0; k < K; k++) {
		for (int j = 0; j < p; j++) {
			double sum = 0.0;
			for (int l = 0; l < L; l++) {
				if (elem(groups, j, l)) {
					sum += d_new[l] * elem(eta_new, l, k);
				}
			}					
			elem(result, j, k) = elem(alpha_new, j, k) * sum;
		}
	}
	return result;
}


SEXP multitask_lasso(SEXP X0, SEXP y0, SEXP lambda0, SEXP model0, SEXP positive0, SEXP eps0)
{
	//convert parameters to Rcpp types
	Rcpp::NumericMatrix X(X0);
	Rcpp::NumericVector y(y0);
	double lambda = Rcpp::as<double>(lambda0);
	double eps = Rcpp::as<double>(eps0);
	int model = Rcpp::as<int>(model0);
	int positive = Rcpp::as<int>(positive0);
	
	shotgun_data data; //input to and output from function solveLasso

	//initialize input fields (shotgun calls the design matrix A instead of X)
	int M = X.nrow();
	int N = X.ncol();
	data.A_cols.resize(N);
	data.A_rows.resize(M);
	for (int k = 0; k < M * N; k++) {
		if (X[k] != 0.0) {
			int i = k % M; //row (R uses column-major order)
			int j = k / M; //column
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
	bool all_zero = false;
	
	if (model == MULTITASK_LINEAR) {
		solveLasso(&data, lambda, positive, regpathLength, eps, maxIter, verbose);
	} else if (model == MULTITASK_LOGISTIC) {
		compute_logreg(&data, lambda, positive, eps, maxIter, verbose, all_zero);
	}else{
		fprintf(stderr, "Unknown method\n");
	}
	
	return Rcpp::wrap(data.x);
}
