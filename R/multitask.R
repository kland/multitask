lasso <- function (X, y, lambda, eps = 1e-12) {
	.Call("multitask_lasso", X, y, lambda, eps, PACKAGE = "multitask")
}

multitask <- function (X, y, tasks, groups, lambda, eps = 1e-12) {
	#TODO: Implement this function properly
	lasso(X, y, lambda, eps)
}

