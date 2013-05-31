#ifndef _multitask_MULTITASK_H
#define _multitask_MULTITASK_H

#include <Rcpp.h>

RcppExport
SEXP multitask_x_tilde(SEXP X, SEXP tasks, SEXP groups, SEXP d_cur, SEXP eta_cur, SEXP K, SEXP k);

RcppExport
SEXP multitask_x_tilde_2(SEXP X, SEXP tasks, SEXP groups, SEXP alpha_new, SEXP eta_cur, SEXP K);

RcppExport
SEXP multitask_x_tilde_3(SEXP X, SEXP tasks, SEXP groups, SEXP alpha_new, SEXP d_new, SEXP K, SEXP k);

RcppExport
SEXP multitask_lasso(SEXP X, SEXP y, SEXP lambda, SEXP eps);



#endif
