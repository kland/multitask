multitask <- function(X, y, tasks, groups, lambda=NULL, nlambda=20, model="linear", standardize=T,eps=1e-6,maxiter=500,maxiter.shotgun=200) {

  tasks <- as.factor(tasks)             # ensure factor format
  K <- length(levels(tasks))            # tasks
  n <- as.numeric(table(tasks))		# replicates
  p <- ncol(X)			        # predictors
  L <- ncol(groups)			# groups

  if (model == "linear") {
    model.num <- 0
  } else if(model=="logistic") {
    model.num <- 1
  }

  if(standardize){
    if(model=="linear")
      y<-as.numeric(unlist(tapply(y,tasks,scale,scale=F,simplify=T)))
    X<-matrix(unlist(apply(X,2,tapply,tasks,scale)),ncol=p)	
  }

  if(is.null(lambda)){
    lams <- c()
    for(k in 1:K){
      task <- levels(tasks)[k]
      if(model=="linear"){
        resids <- y[tasks==task]
      }else if(model=="logistic"){
        resids <- y[tasks==task] - 0.5
      } 
      task <- levels(tasks)[k]
      lams <- c(lams,abs(drop(crossprod(X[tasks==task,], resids))))
    }
    lambda.max <- max(lams)
    lambda.min <- 2*sqrt(max(n)+p)/10
    
    ##lambda.min <- 0.05
    lambda <- seq(lambda.max,lambda.min,length.out=nlambda)
  }

  fit <- NULL; nlam <- length(lambda)
  fit$beta <-  fit$alpha <- array(NA,dim=c(p,K,nlam))
  fit$eta <- array(NA,dim=c(L,K,nlam))
  fit$d <- matrix(NA,nr=L,nc=nlam)
  for(i in 1:nlam){
    lam <- lambda[i]
    temp.fit <- .Call("multitask", X, y, K, groups, lam, model.num, eps, maxiter, maxiter.shotgun, PACKAGE = "multitask")
    fit$beta[,,i] <- temp.fit$beta
    fit$alpha[,,i] <- temp.fit$alpha
    fit$eta[,,i] <- temp.fit$eta
    fit$d[,i] <- temp.fit$d
    fit$bic <- cbind(fit$bic, as.numeric(temp.fit$bic)) 
    fit$converged <- cbind(fit$converged, temp.fit$converged) 

  }
  fit$tasks <- tasks
  fit$groups <- groups
  fit$lambda <- lambda
  fit
}
