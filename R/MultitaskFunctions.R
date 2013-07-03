multitask <- function(X, y, tasks, groups, lambda, model="linear", eps=1e-12) {

	tasks <- as.factor(tasks)
	K <- length(levels(tasks))    # tasks

	if (model == "linear") {
		model.num <- 0
	} else if(model=="logistic") {
		model.num <- 1
	}
	fit <- .Call("multitask", X, y, K, groups, lambda, model.num, eps, PACKAGE = "multitask")
	fit$tasks <- tasks
	fit$groups <- groups
	fit$lambda <- lambda
	fit
}
