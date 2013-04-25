lasso <- function (X, y, tasks, groups, lambda, eps = 1e-12) {
	.Call("multitask_lasso", X, y, tasks, groups, lambda, eps, PACKAGE = "multitask")
}

