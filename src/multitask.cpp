#include <Rcpp.h>
#include "multitask.h"

SEXP multitask()
{
	Rcpp::CharacterVector x = Rcpp::CharacterVector::create("foo", "bar");
	Rcpp::NumericVector y = Rcpp::NumericVector::create(0.0, 1.0);
	Rcpp::List z = Rcpp::List::create(x, y);

	return z ;
}
