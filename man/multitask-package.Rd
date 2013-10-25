\name{multitask-package}
\alias{multitask-package}
\docType{package}
\title{
Variable Selection With Group Overlap in a Multitask Setting
}
\description{
More about what it does (maybe more than one line)
~~ A concise (1-5 lines) description of the package ~~
}
\details{
\tabular{ll}{
Package: \tab multitask\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-10-25\cr
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
\keyword{package}
\keyword{linear regression}
\keyword{logistic regression}
\keyword{sparse models}
\keyword{multi-task problem}
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
