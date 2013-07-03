#include <math.h>
#include <stdlib.h>
#include <Rcpp.h>
#include "shotgun/common.h"
#include "multitask.h"

static int sign(double x)
{
	int result;
	
	if (x < 0.0) {
		result = -1;
	} else if (x == 0.0) {
		result = 0;
	} else {
		result = 1;
	}
	return result;
}


static double random(double min, double max) /*return (pseudo) random number between min and max*/
{
	return ((double) rand()) / ((double) RAND_MAX) * (max - min) + min;
}


static inline int &elem(Rcpp::IntegerMatrix &A, int i, int j) //returns element at position (i, j) in column-major order matrix A
{
	return A[A.nrow() * j + i];
}


static inline double &elem(Rcpp::NumericMatrix &A, int i, int j) //returns element at position (i, j) in column-major order matrix A
{
	return A[A.nrow() * j + i];
}


static void print(Rcpp::NumericMatrix A)
{
	for (int i = 0; i < A.nrow(); i++) {
		for (int j = 0; j < A.ncol(); j++) {
			printf("%9f ", elem(A, i, j));
		}
		putchar('\n');
	}
}


static void print(Rcpp::IntegerVector a)
{
	for (int i = 0; i < a.size(); i++) {
		printf("%d ", a[i]);
	}
	putchar('\n');
}


static void print(Rcpp::NumericVector a)
{
	for (int i = 0; i < a.size(); i++) {
		printf("%f ", a[i]);
	}
	putchar('\n');
}


static Rcpp::NumericMatrix x_tilde(Rcpp::NumericMatrix X, 
	int K,
	Rcpp::IntegerMatrix groups, 
	Rcpp::IntegerVector d_cur, 
	Rcpp::NumericMatrix eta_cur)
{
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


static Rcpp::NumericMatrix x_tilde_2(Rcpp::NumericMatrix X, 
	int K,
	Rcpp::IntegerMatrix groups, 
	Rcpp::NumericMatrix alpha_new,
	Rcpp::NumericMatrix eta_cur)
{
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


static Rcpp::NumericMatrix x_tilde_3(Rcpp::NumericMatrix X, 
	int K,
	Rcpp::IntegerMatrix groups, 
	Rcpp::NumericMatrix alpha_new,
	Rcpp::IntegerVector d_new)
{
	assert(K > 0);

	int n = X.nrow() / K;
	int p = X.ncol();
	int L = groups.ncol();
	Rcpp::NumericMatrix result(n * K, L * K);
	
	for (int k = 0; k < K; k++) {
		for (int l = 0; l < L; l++) {
			for (int i = 0; i < n; i++) {
				double sum = 0.0;
				for (int j = 0; j < p; j++) {
					if (elem(groups, j, l)) {
						sum += elem(X, n * k + i, j) * elem(alpha_new, j, k);
					}
				}
				elem(result, n * k + i, L * k + l) = d_new[l] * sum;
			}
		}
	}	
	return result;
}


static Rcpp::NumericMatrix next_beta(int K,
	Rcpp::IntegerMatrix groups, 
	Rcpp::NumericMatrix alpha_new,
	Rcpp::IntegerVector d_new,
	Rcpp::NumericMatrix eta_new)
{
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


static double bic(Rcpp::NumericMatrix X, 
	Rcpp::NumericVector y, 
	Rcpp::NumericMatrix beta_new, 
	double eps, 
	int n)
{
	assert(n > 0);
	
	int p = X.ncol();
	int K = X.nrow() / n;
	
	/*calculate df*/
	int df = 0;
	int N = beta_new.nrow() * beta_new.ncol();
	for (int i = 0; i < N; i++) {
		if (fabs(beta_new[i]) > eps) {
			 df++;
		}
	}

	/*calculate SSe*/
	double SSe = 0.0;
	for (int k = 0; k < K; k++) {
		for (int i = k * n; i < (k + 1) * n; i++) {
			double Xrow_betacol = 0.0;
			for (int j = 0; j < p; j++) {
				Xrow_betacol += elem(X, i, j) * elem(beta_new, j, k);
			}
			SSe += pow(y[i] - Xrow_betacol, 2);
		}
	}
	
	double ll = -y.size() / 2.0 * (log(SSe) - log(y.size()) + log(2.0 * M_PI) + 1);
	double result = -2 * ll + df * log(y.size());

	return result;
}


static Rcpp::NumericVector lasso(Rcpp::NumericMatrix X, 
	Rcpp::NumericVector y, 
	double lambda, 
	int model, 
	bool positive, 
	double eps)
{	
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
		solveLasso(&data, lambda * 2, positive, regpathLength, eps, maxIter, verbose);
	} else if (model == MULTITASK_LOGISTIC) {
		compute_logreg(&data, lambda, positive, eps, maxIter, verbose, all_zero);
	}else{
		fprintf(stderr, "Unknown method\n");
	}
	
	return Rcpp::wrap(data.x);
}


SEXP multitask(SEXP X0, SEXP y0, SEXP K0, SEXP groups0, SEXP lambda0, SEXP model0, SEXP eps0)
{
	Rcpp::NumericVector lasso_result;

	//convert parameters to Rcpp types
	Rcpp::NumericMatrix X(X0);
	Rcpp::NumericVector y(y0);
	int K = Rcpp::as<int>(K0);
	Rcpp::IntegerMatrix groups(groups0);
	double lambda = Rcpp::as<double>(lambda0);
	int model = Rcpp::as<int>(model0);
	double eps = Rcpp::as<double>(eps0);

	assert(K > 0);
	assert(lambda >= 0.0);
	assert(model >= 0);
	assert(model < MULTITASK_MODEL_COUNT);
	assert(eps > 0.0);
	
	int n = y.size() / K;
	int p = groups.nrow();
	int L = groups.ncol();
	
	//initialize start values
	
	Rcpp::NumericMatrix alpha_cur(p, K);
	for (int i = 0; i < p * K; i++) {
		alpha_cur[i] = random(-0.1, 0.1);
	}
	
	Rcpp::NumericMatrix beta_cur = alpha_cur;

	Rcpp::IntegerVector d_cur(L);
	std::fill(d_cur.begin(), d_cur.end(), 1);
	
	Rcpp::NumericMatrix eta_cur(L, K);
	std::fill(eta_cur.begin(), eta_cur.end(), 1.0);
	
	Rcpp::NumericMatrix beta_new;
	
	bool converged = false;

	while (! converged) {
		//update alpha
		Rcpp::NumericMatrix Xtilde = x_tilde(X, K, groups, d_cur, eta_cur);
		lasso_result = lasso(Xtilde, y, lambda, model, false, eps);
		assert(lasso_result.size() == p * K);
		Rcpp::NumericMatrix alpha_new(p, K, lasso_result.begin());
		//printf("C++: alpha_new = \n"); print(alpha_new);

		//update d
		Rcpp::NumericMatrix Xtilde2 = x_tilde_2(X, K, groups, alpha_new, eta_cur);
		lasso_result = lasso(Xtilde2, y, 1.0, model, true, eps);
		assert(lasso_result.size() == L);
		Rcpp::IntegerVector d_new(L);
		for (int i = 0; i < L; i++) {
			d_new[i] = sign(lasso_result[i]);
		}
		//printf("C++: d_new = \n"); print(d_new);
		
		//update eta
		Rcpp::NumericMatrix Xtilde3 = x_tilde_3(X, K, groups, alpha_new, d_new);
		lasso_result = lasso(Xtilde3, y, 1.0, model, true, eps);
		assert(lasso_result.size() == L * K);
		Rcpp::NumericMatrix eta_new(L, K, lasso_result.begin());
		//printf("C++: eta_new = \n"); print(eta_new);

		//update beta
		beta_new = next_beta(K, groups, alpha_new, d_new, eta_new);
		assert(beta_new.nrow() == p);
		assert(beta_new.ncol() == K);
		//printf("C++: beta_new = \n"); print(beta_new);
		
		double max_diff = 0.0;
		for (int i = 0; i < p * K; i++) {
			double diff = fabs(beta_new[i] - beta_cur[i]);
			if (diff > max_diff) {
				max_diff = diff;
			}
		}
		if (max_diff < eps) {
			converged = true;
		}
		
		alpha_cur = alpha_new;
		d_cur = d_new;
		eta_cur = eta_new;
		beta_cur = beta_new;
	}

	Rcpp::List result;
	result["beta"] = beta_cur;
	result["alpha"] = alpha_cur;
	result["d"] = d_cur;
	result["eta"] = eta_cur;
	result["bic"] = bic(X, y, beta_new, eps, n);
	
	return result;
}
