#ifndef _multitask_MULTITASK_H
#define _multitask_MULTITASK_H

#include <Rcpp.h>

/*model values*/
#define MULTITASK_LINEAR 0
#define MULTITASK_LOGISTIC 1
#define MULTITASK_MODEL_COUNT 2

RcppExport{
  SEXP multitask(SEXP X, SEXP y, SEXP nk, SEXP groups, SEXP lambda, SEXP corrfactor, SEXP model, SEXP conveps, SEXP eps, SEXP maxiter,SEXP maxitersg);
  SEXP grplasso(SEXP X, SEXP y, SEXP n, SEXP groups, SEXP lambda, SEXP model, SEXP conveps, SEXP eps, SEXP maxiter,SEXP maxitersg);
}

#endif
