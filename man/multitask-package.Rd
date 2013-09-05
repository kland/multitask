\name{multitask-package}
\alias{multitask-package}
\alias{multitask}
\docType{package}
\title{
Variable Selection With Group Overlap in a Multitask Setting
}
\description{
More about what it does (maybe more than one line)
~~ A concise (1-5 lines) description of the package ~~
}
\usage{
multitask(X, y, tasks, groups, lambda = NULL, nlambda = 20, model = "linear", standardize = T, eps = 1e-6, maxiter = 500, maxiter.shotgun = 300)	
}
\arguments{
\item{X}{feature matrix of (n * K) x p real numbers where n is the number of replicates, K is the number of tasks (cohorts) and p i s the number of predictors (features)}
\item{y}{response vector of n * K real numbers}
\item{tasks}{indicator vector of n * K characters}
\item{groups}{matrix of p x L integers (booleans) where L is the number of groups (metabolites)}
\item{lambda}{weight vector of positive real numbers}
\item{nlambda}{number of lambda values (only used if lambda is NULL)}
\item{model}{"linear" or "logistic"}
\item{standardize}{boolean value which determines whether the variables should be standardized or not}
\item{eps}{positive real number which sets the convergency threshold}
\item{maxiter}{positive integer which sets the maximum number of iterations}
\item{maxiter.shotgun}{?}
}
\details{
\tabular{ll}{
Package: \tab multitask\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-09-04\cr
License: \tab GPL-3\cr
}
~~ An overview of how to use the package, including the most important ~~
~~ functions ~~
}
\author{
Alexandra Jauhiainen and Karl Landström

Maintainer: Alexandra Jauhiainen <alexandra.jauhiainen@ki.se>, Karl Landström <karl.landstrom@ki.se>
}
\references{
Shotgun C++ package: http://www.select.cs.cmu.edu/projects/shotgun/
}
~~ Optionally other standard keywords, one per line, from file KEYWORDS in ~~
~~ the R documentation directory ~~
\keyword{ package }
\seealso{
~~ Optional links to other man pages, e.g. ~~
~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{
require(multitask)

n <- 2
p <- 2
L <- 2
K <- 2

X <- matrix(c(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0), nrow = n * K, ncol = p)
y <- c(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0)
groups <- matrix(c(1, 0, 1, 1), nrow = p, ncol = L)
tasks <- rep(LETTERS[1:K], each = n)

out <- multitask(X, y, tasks, groups)
}
