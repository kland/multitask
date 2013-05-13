#ifndef _multitask_MULTITASK_H
#define _multitask_MULTITASK_H

#include <Rcpp.h>

RcppExport SEXP multitask_x_tilde(SEXP X, SEXP tasks, SEXP groups, SEXP d_cur, SEXP eta_cur, SEXP task_index);

RcppExport SEXP multitask_lasso(SEXP X, SEXP y, SEXP lambda, SEXP eps);

#endif
