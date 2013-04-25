#include <Rcpp.h>
#include "multitask.h"

SEXP multitask(SEXP X0, SEXP y0, SEXP tasks0, SEXP groups0)
{
	/*convert SEXP parameters to Rcpp types*/
	Rcpp::NumericMatrix X(X0);
	Rcpp::NumericVector y(y0);
	Rcpp::CharacterVector tasks(tasks0);
	Rcpp::NumericMatrix groups(groups0);
	Rcpp::NumericMatrix result(X0);
	
	return result;
}
