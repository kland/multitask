#ifndef _multitask_MULTITASK_H
#define _multitask_MULTITASK_H

#include <Rcpp.h>

RcppExport SEXP multitask(SEXP X, SEXP y, SEXP tasks, SEXP groups, SEXP lambda, SEXP eps);

#endif
