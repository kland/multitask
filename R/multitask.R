multitask <- function (X, y, tasks, groups, lambda, eps = 1e-12) {
	.Call("multitask", X, y, tasks, groups, lambda, eps, PACKAGE = "multitask")
}

