reload.package <- function () { #reloads the multitask package
	detach("package:multitask", unload = TRUE)
	library.dynam.unload("multitask", system.file(package = "multitask"))
	install.packages(".", repos = NULL, type = "source")
	library("multitask")
}
