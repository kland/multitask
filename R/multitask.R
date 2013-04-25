multitask <- function (X, y, tasks, groups) {
	.Call("multitask", X, y, tasks, groups, PACKAGE = "multitask")
}

