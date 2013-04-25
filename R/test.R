test <- function () {
	n <- 2
	p <- 2
	L <- 2
	K <- 2

	#initialize input parameters

	groups <- matrix(c(1, 0, 1, 1), nrow = p, ncol = L)
	rownames(groups) <- paste('gene', 1:p, sep = ' ')
	colnames(groups) <- paste('group', 1:L, sep = ' ')

	##values for X and y comes from running the script SimStudy-Linear.R once (want the same (test) input for each run).

	X <- matrix(c(0.5498018, 1.3497705, 0.1617772, -0.9073028, 0.9164897, 1.7140090, -0.1020368, -0.4579487), nrow = n * K, ncol = p)
	y <- c(4.3328193, 2.0634363, -1.2631156, 0.7304365)
	tasks <- rep(LETTERS[1:K], each = n)
	lambda <- 40	
	
	#estimate parameters

	out <- lasso(X, y, tasks, groups, lambda)

	print(out)
}
