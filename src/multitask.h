#ifndef _multitask_MULTITASK_H
#define _multitask_MULTITASK_H

#include <Rcpp.h>

/*model values*/
#define MULTITASK_LINEAR 0
#define MULTITASK_LOGISTIC 1
#define MULTITASK_MODEL_COUNT 2

RcppExport
SEXP multitask(SEXP X, SEXP y, SEXP K, SEXP groups, SEXP lambda, SEXP model, SEXP eps, SEXP maxiter);

#endif
