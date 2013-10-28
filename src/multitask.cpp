#include <math.h>
#include <stdlib.h>
#include <Rcpp.h>
#include "shotgun/common.h"
#include "multitask.h"

static int nz(double x,double eps)
{
  int result;
  if (fabs(x) < eps) {
    result = 0;
  } else {
    result = 1;
  }
  return result;
}

static Rcpp::IntegerMatrix nz(Rcpp::NumericMatrix m, double eps)
{
  int nr = m.nrow();
  int nc = m.ncol();
  Rcpp::IntegerMatrix result(nr, nc);
  for(int i=0; i < nr*nc; i++){
    result[i] = nz(m[i],eps);
  }
  return result;
}

static Rcpp::IntegerVector nz(Rcpp::NumericVector v, double eps)
{
  int n = v.size();
  Rcpp::IntegerVector result(n);
  
  for(int i=0; i < n; i++){
    result[i] = nz(v[i],eps);
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
				   Rcpp::IntegerVector nk,
				   Rcpp::IntegerMatrix groups, 
				   Rcpp::NumericVector d_cur, 
				   Rcpp::NumericMatrix eta_cur)
{
  int K = nk.size();
  int n_tot = X.nrow();
  int p = X.ncol();
  int L = groups.ncol();
  Rcpp::NumericMatrix result(n_tot, p * K);

  int idx = 0;
  for (int k = 0; k < K; k++) {
    int n = nk[k];
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
	elem(result, idx + i, p * k + j) = elem(X, idx + i, j) * sum;
      }
    }
    idx += n;
  }
  return result;
}

static Rcpp::NumericMatrix x_tilde_2(Rcpp::NumericMatrix X, 
				     Rcpp::IntegerVector nk,
				     Rcpp::IntegerMatrix groups, 
				     Rcpp::NumericMatrix alpha_new,
				     Rcpp::NumericMatrix eta_cur)
{	
  int K = nk.size();
  int n_tot = X.nrow();
  int p = X.ncol();
  int L = groups.ncol();
  Rcpp::NumericMatrix result(n_tot, L);

  for (int l = 0; l < L; l++) {
    int k = -1;
    int n = 0;    
    for (int i = 0; i < n_tot; i++) {
      if (i == n){
	k +=1;
	n += nk[k];
      }
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
				     Rcpp::IntegerVector nk,
				     Rcpp::IntegerMatrix groups, 
				     Rcpp::NumericMatrix alpha_new,
				     Rcpp::NumericVector d_new)
{
  int K = nk.size();
  int n_tot = X.nrow();
  int p = X.ncol();
  int L = groups.ncol();
  Rcpp::NumericMatrix result(n_tot, L * K);
  
  int idx = 0;
  for (int k = 0; k < K; k++) {
    int n = nk[k];
    for (int l = 0; l < L; l++) {	
      for (int i = 0; i < n; i++) {
	double sum = 0.0;
	for (int j = 0; j < p; j++) {
	  if (elem(groups, j, l)) {
	    sum += elem(X, idx + i, j) * elem(alpha_new, j, k);
	  }
	}
	elem(result, idx + i, L * k + l) = d_new[l] * sum;
      }
    }
    idx += n;
  }	
  return result;
}

static Rcpp::IntegerVector nz_vec(Rcpp::NumericMatrix alpha_new,
				  Rcpp::NumericMatrix eta_new,
				  Rcpp::NumericVector d_new,
				  double eps)
{
  int K = alpha_new.ncol();
  int p = alpha_new.nrow();
  int L = eta_new.nrow();
  
  Rcpp::IntegerVector result(p*K + L*K + L);
 	
  for (int i = 0; i < p*K; i++) {
    result[i] = nz(alpha_new[i],eps);
  }

  for (int i = 0; i < L*K; i++) {
    result[p*K + i] = nz(eta_new[i],eps);
  }

  for (int i = 0; i < L; i++) {
    result[p*K + L*K + i] = nz(d_new[i],eps);
  }

  return result;
}

static Rcpp::IntegerVector nz_vec(Rcpp::NumericMatrix alpha_new,
				  Rcpp::NumericVector d_new,
				  double eps)
{
  int K = alpha_new.ncol();
  int p = alpha_new.nrow();
  int L = d_new.size(); 
  Rcpp::IntegerVector result(p*K + L);
 	
  for (int i = 0; i < p*K; i++) {
    result[i] = nz(alpha_new[i],eps);
  }

  for (int i = 0; i < L; i++) {
    result[p*K + i] = nz(d_new[i],eps);
  }

  return result;
}

static Rcpp::NumericMatrix next_beta(Rcpp::IntegerVector nk,
				     Rcpp::IntegerMatrix groups, 
				     Rcpp::NumericMatrix alpha_new,
				     Rcpp::NumericVector d_new,
				     Rcpp::NumericMatrix eta_new)
{
  int K = nk.size();
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

static int df(Rcpp::NumericMatrix beta_new, double eps)
{
  int result = 0;
  int n = beta_new.nrow() * beta_new.ncol();
  for(int i=0; i < n; i++){
    result += nz(beta_new[i],eps);
  }
  return result;
}


static double bic_linear(Rcpp::NumericMatrix X, 
			 Rcpp::NumericVector y, 
			 Rcpp::NumericMatrix beta_new, 
			 double eps, 
			 Rcpp::IntegerVector nk)
{
  int n_tot = X.nrow();
  int p = X.ncol();
  int K = nk.size();
 

  /*calculate SSe*/
  double SSe = 0.0;  
  int idx = 0;
  
  for (int k = 0; k < K; k++) {
    int n = nk[k];
    for (int i = 0; i < n; i++) {
      double Xrow_betacol = 0.0;
      for (int j = 0; j < p; j++) {
	Xrow_betacol += elem(X, idx+i, j) * elem(beta_new, j, k);
      }
      SSe += pow(y[idx+i] - Xrow_betacol, 2);
    }
    idx += n;
  }
  
  double ll = -n_tot / 2.0 * (log(SSe) - log(n_tot) + log(2.0 * M_PI) + 1);
  double bic = -2 * ll + df(beta_new, eps) * log(n_tot);

  return bic;
}

static double bic_logistic(Rcpp::NumericMatrix X, 
			   Rcpp::NumericVector y, 
			   Rcpp::NumericMatrix beta_new, 
			   double eps, 
			   Rcpp::IntegerVector nk)
{
  int n_tot = X.nrow();
  int p = X.ncol();
  int K = nk.size();
    
  int idx = 0;
  double ll = 0.0;
  for (int k = 0; k < K; k++) {
    int n = nk[k];
    for (int i = 0; i < n; i++) {
      double lp = 0.0;
      for (int j = 0; j < p; j++) {
	lp += elem(X, idx+i, j) * elem(beta_new, j, k);
      }
      ll += y[idx+i] * lp - log(1.0 + exp(lp));
    }
    idx += n;
  }
  
  double bic = -2.0 * ll + df(beta_new, eps) * log(n_tot);
  return bic;
}



static Rcpp::NumericVector lasso(Rcpp::NumericMatrix X, 
				 Rcpp::NumericVector y, 
				 double lambda, 
				 int model, 
				 bool positive, 
				 double eps,
				 int maxiter)
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
  int verbose = 0;
  bool all_zero = false;
  
  if (model == MULTITASK_LINEAR) {
    solveLasso(&data, lambda * 2, positive, regpathLength, eps, maxiter, verbose);
  } else if (model == MULTITASK_LOGISTIC) {
    for (int i = 0; i < data.y.size(); i++) {
      if (data.y[i] == 0.0) {
	data.y[i] = -1.0;
      }
    }	
    compute_logreg(&data, lambda, positive, eps, maxiter, verbose, all_zero);
  }else{
    fprintf(stderr, "Unknown method\n");
  }
  
  return Rcpp::wrap(data.x);
}
static Rcpp::NumericMatrix refit_model(Rcpp::NumericMatrix X,
				       Rcpp::NumericVector y,
				       Rcpp::NumericMatrix beta_new,
				       Rcpp::IntegerVector nk,
				       int model,
				       double eps,
				       int maxiter)
{
  int K = nk.size();
  int p = X.ncol();
  int n;
  Rcpp::NumericMatrix Xtmp(X.nrow(),p);
  Rcpp::NumericMatrix beta_refit(p,K);
  Rcpp::NumericVector lasso_result;
 
  int idx = 0;
  for (int k = 0; k < K; k++) {
    n = nk[k];
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < p; j++) {
	Xtmp(idx+i,j) = X(idx+i,j) * nz(beta_new(j,k),eps);
      }
    }
    idx += n;
  }
   
 
  idx = 0;
  for (int k = 0; k < K; k++) {
    n = nk[k];
    lasso_result = lasso(Xtmp(Rcpp::Range(idx,idx+n-1),Rcpp::_), y[Rcpp::Range(idx,idx+n-1)], 0.0, model, false, eps, maxiter);
    for(int j = 0; j < p; j++){
      beta_refit(j,k) = lasso_result[j];
    }
    idx += n;
  }
  
  return beta_refit;
}


SEXP multitask(SEXP X0, SEXP y0, SEXP nk0, SEXP groups0, SEXP lambda0, SEXP corrfactor0, SEXP model0, SEXP conveps0, SEXP eps0, SEXP maxiter0,SEXP maxitersg0)
{
  Rcpp::NumericVector lasso_result;
  
  //convert parameters to Rcpp types
  Rcpp::NumericMatrix X(X0);
  Rcpp::NumericVector y(y0);
  Rcpp::IntegerVector nk(nk0);
  int K = nk.size();
  Rcpp::IntegerMatrix groups(groups0);
  Rcpp::NumericVector corrfactor(corrfactor0);
  double lambda = Rcpp::as<double>(lambda0);
  int model = Rcpp::as<int>(model0);
  double eps = Rcpp::as<double>(eps0);
  double conveps = Rcpp::as<double>(conveps0);
  int maxiter = Rcpp::as<int>(maxiter0);
  int maxitersg = Rcpp::as<int>(maxitersg0);
 
  assert(K > 0);
  assert(lambda >= 0.0);
  assert(model >= 0);
  assert(model < MULTITASK_MODEL_COUNT);
  assert(eps > 0.0);
  assert(maxiter > 0);
  assert(maxitersg > 0);
 
  int p = groups.nrow();
  int L = groups.ncol();
  
  //initialize start values  
  Rcpp::NumericMatrix alpha_cur(p, K);
  for (int i = 0; i < p * K; i++) {
    alpha_cur[i] = random(-0.1, 0.1);
  }
  
  Rcpp::NumericMatrix beta_cur = alpha_cur;
  
  Rcpp::NumericVector d_cur(L);
  std::fill(d_cur.begin(), d_cur.end(), 1.0);
  
  Rcpp::NumericMatrix eta_cur(L, K);
  std::fill(eta_cur.begin(), eta_cur.end(), 1.0);
	
  Rcpp::IntegerVector nz_cur(p*K + L*K + L);
  std::fill(nz_cur.begin(), nz_cur.end(), 1);

  Rcpp::NumericMatrix beta_new;
  
  bool converged = false;
  int iterations = 0;
    
  do {
    ++iterations;
    if (iterations >= maxiter && maxiter > 0) {
      break;
    }
    //update alpha
    Rcpp::NumericMatrix alpha_new(p, K);
    Rcpp::NumericMatrix Xtilde = x_tilde(X, nk, groups, d_cur, eta_cur);
    int idx = 0; int n; 
    for(int k=0; k < K; k++){
      n = nk[k];
      lasso_result = lasso(Xtilde(Rcpp::Range(idx,idx+n-1), Rcpp::Range(p*k, p*(k+1)-1)), y[Rcpp::Range(idx,idx+n-1)], lambda*corrfactor[k], model, false, eps, maxitersg);
      assert(lasso_result.size() == p);
      for (int j = 0; j < p; j++) {
	alpha_new(j,k) = lasso_result[j];
      }
      idx += n; 
    }
  			
    //update d
    Rcpp::NumericMatrix Xtilde2 = x_tilde_2(X, nk, groups, alpha_new, eta_cur);
    lasso_result = lasso(Xtilde2, y, 1.0, model, true, eps, maxitersg);
    assert(lasso_result.size() == L);
    Rcpp::NumericVector d_new(L);
    for (int i = 0; i < L; i++) {
      d_new[i] = lasso_result[i]/max(lasso_result);
      if(std::isnan(d_new[i])) d_new[i] = 0;
   }
 		
    //update eta
    Rcpp::NumericMatrix eta_new(L, K);
    Rcpp::NumericMatrix Xtilde3 = x_tilde_3(X, nk, groups, alpha_new, d_new);
    idx = 0;
    for(int k=0; k < K; k++){
      n = nk[k];
      lasso_result = lasso(Xtilde3(Rcpp::Range(idx,idx+n-1),Rcpp::Range(L*k, L*(k+1)-1)), y[Rcpp::Range(idx,idx+n-1)], 1.0, model, true, eps, maxitersg);
      assert(lasso_result.size() == L);
      for (int l = 0; l < L; l++) {
	eta_new(l,k) = lasso_result[l];
      }
      idx += n; 
    }
  		
    //update beta
    beta_new = next_beta(nk, groups, alpha_new, d_new, eta_new);
    assert(beta_new.nrow() == p);
    assert(beta_new.ncol() == K);
    
    //check structure convergence
    Rcpp::IntegerVector nz_new = nz_vec(alpha_new,eta_new,d_new,eps);
    
    int nz_diff = 0;
    for (int i = 0; i < nz_new.length(); i++) {
      nz_diff += nz_cur[i] - nz_new[i];
    }

    double max_diff = 0.0;
    for (int i = 0; i < p * K; i++) {
      double diff = fabs(beta_new[i] - beta_cur[i]);
      if (diff > max_diff) {
	max_diff = diff;
      }
    }
    if (max_diff < conveps && nz_diff==0) {
      converged = true;
    }
    
    alpha_cur = alpha_new;
    d_cur = d_new;
    eta_cur = eta_new;
    beta_cur = beta_new;
    nz_cur = nz_new;
  } while (converged != true);
  
  Rcpp::NumericMatrix beta_refit = refit_model(X, y,beta_new, nk, model, eps, maxitersg);

  Rcpp::List result;
  result["beta"] = beta_refit;
  result["d"] = nz(d_cur,eps);
  result["eta"] = nz(eta_cur,eps);
  switch (model) {
  case MULTITASK_LINEAR:
    result["bic"] = bic_linear(X, y, beta_refit, eps, nk);
    break;
  case MULTITASK_LOGISTIC:
    result["bic"] = bic_logistic(X, y, beta_refit, eps, nk);
    break;
  default:
    assert(0);
  }
  result["converged"] = converged;
  
  return result;
}

SEXP grplasso(SEXP X0, SEXP y0, SEXP n0, SEXP groups0, SEXP lambda0, SEXP model0, SEXP conveps0, SEXP eps0, SEXP maxiter0,SEXP maxitersg0)
{
  Rcpp::NumericVector lasso_result;
  
  //convert parameters to Rcpp types
  Rcpp::NumericMatrix X(X0);
  Rcpp::NumericVector y(y0);
  Rcpp::IntegerMatrix groups(groups0);
  double lambda = Rcpp::as<double>(lambda0);
  int model = Rcpp::as<int>(model0);
  double eps = Rcpp::as<double>(eps0);
  double conveps = Rcpp::as<double>(conveps0);
  int maxiter = Rcpp::as<int>(maxiter0);
  int maxitersg = Rcpp::as<int>(maxitersg0);
  Rcpp::IntegerVector n(n0);
  int K = 1;
 
 
  assert(lambda >= 0.0);
  assert(model >= 0);
  assert(model < MULTITASK_MODEL_COUNT);
  assert(eps > 0.0);
  assert(maxiter > 0);
  assert(maxitersg > 0);
 
  int p = groups.nrow();
  int L = groups.ncol();
  
  //initialize start values  
  Rcpp::NumericMatrix alpha_cur(p, K);
  for (int i = 0; i < p * K; i++) {
    alpha_cur[i] = random(-0.1, 0.1);
  }
  
  Rcpp::NumericMatrix beta_cur = alpha_cur;
  
  Rcpp::NumericVector d_cur(L);
  std::fill(d_cur.begin(), d_cur.end(), 1.0);
  
  Rcpp::NumericMatrix eta_cur(L, K);
  std::fill(eta_cur.begin(), eta_cur.end(), 1.0);

  Rcpp::IntegerVector nz_cur(p*K + L);
  std::fill(nz_cur.begin(), nz_cur.end(), 1);

  Rcpp::NumericMatrix beta_new;
 
  bool converged = false;
  int iterations = 0;
  
  do {
    ++iterations;
    if (iterations >= maxiter && maxiter > 0) {
      break;
    }
    //update alpha
    Rcpp::NumericMatrix alpha_new(p, K);
    Rcpp::NumericMatrix Xtilde = x_tilde(X, n, groups, d_cur, eta_cur);
    lasso_result = lasso(Xtilde, y, lambda, model, false, eps, maxitersg);
    assert(lasso_result.size() == p);
    for (int j = 0; j < p; j++) {
      alpha_new(j,1) = lasso_result[j];
    }
    
    //update d
    Rcpp::NumericMatrix Xtilde2 = x_tilde_2(X, n, groups, alpha_new, eta_cur);
    lasso_result = lasso(Xtilde2, y, 1.0, model, true, eps, maxitersg);
    assert(lasso_result.size() == L);
    Rcpp::NumericVector d_new(L);
    for (int i = 0; i < L; i++) {
      d_new[i] = lasso_result[i]/max(lasso_result);
      if(std::isnan(d_new[i])) d_new[i] = 0;
    }
  		
    //update beta
    beta_new = next_beta(n, groups, alpha_new, d_new, eta_cur);
    assert(beta_new.nrow() == p);
    assert(beta_new.ncol() == K);
    
    //check structure convergence
    Rcpp::IntegerVector nz_new = nz_vec(alpha_new,d_new,eps);
  
    int nz_diff = 0;
    for (int i = 0; i < nz_new.length(); i++) {
      nz_diff += nz_cur[i] - nz_new[i];
    }

    double max_diff = 0.0;
    for (int i = 0; i < p * K; i++) {
      double diff = fabs(beta_new[i] - beta_cur[i]);
      if (diff > max_diff) {
	max_diff = diff;
      }
    }
   if (max_diff < conveps && nz_diff==0) {
     converged = true;
    }
    
    alpha_cur = alpha_new;
    d_cur = d_new;
    beta_cur = beta_new;
    nz_cur = nz_new;
  } while (converged != true);
  
  Rcpp::NumericMatrix beta_refit = refit_model(X, y, beta_new, n, model, eps, maxitersg);

  Rcpp::List result;
  result["beta"] = beta_refit;
  result["d"] = nz(d_cur,eps);
  switch (model) {
  case MULTITASK_LINEAR:
    result["bic"] = bic_linear(X, y, beta_new, eps, n);
    break;
  case MULTITASK_LOGISTIC:
    result["bic"] = bic_logistic(X, y, beta_new, eps, n);
    break;
  default:
    assert(0);
  }
  result["converged"] = converged;
  
  return result;
}
