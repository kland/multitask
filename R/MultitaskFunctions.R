multitask <- function(X, y, tasks, groups, lambda=NULL, nlambda=20, model="linear", standardize=T,delta=1,conv.eps=1e-3, eps=1e-6,maxiter=100,maxiter.shotgun=300) {

  tasks <- as.factor(tasks)             # ensure factor format
  K <- length(levels(tasks))            # tasks
  p <- ncol(X)			        # predictors
  L <- ncol(groups)			# groups
  nk <- as.numeric(table(tasks)[unique(tasks)])      # number of obs in each task
  nkm <- mean(nk)
    
  if (model == "linear") {
    model.num <- 0
  } else if(model=="logistic") {
    model.num <- 1
    y <- 1*(y==max(y))
  }

  if(standardize){
    if(model=="linear")
      y<-as.numeric(unlist(tapply(y,tasks,scale,scale=F,simplify=T)))
    X<-matrix(unlist(apply(X,2,tapply,tasks,scale)),ncol=p)	
  }

  if(is.null(lambda)){
    lambda <- calcLambda(X,y,tasks,nlambda,model)
  }
  lambda <- sort(lambda,decreasing=T)
     
  fit <- NULL; nlambda <- length(lambda)
  for(i in 1:nlambda){
    corr.factor <- (nkm^delta * nk^(1-delta))/nkm
    temp.fit <- .Call("multitask", X, y, nk, groups, lambda[i], corr.factor,  model.num, conv.eps, eps, maxiter, maxiter.shotgun, PACKAGE = "multitask")
    if(temp.fit$converged){
      fit$converged <- c(fit$converged, temp.fit$converged) 
      fit$beta <- cbind(fit$beta, as.numeric(temp.fit$beta))
      fit$eta <- cbind(fit$eta,as.numeric(temp.fit$eta))
      fit$d <- cbind(fit$d, temp.fit$d)
      fit$bic <- c(fit$bic, as.numeric(temp.fit$bic)) 
      fit$lambda <- c(fit$lambda,lambda[i])
      fit$corr.factor <- cbind(fit$corr.factor,as.numeric(corr.factor))
    }
  }

  if(is.null(fit)){
    warning("No lambdas converged. Try to give larger lambdas or, as a second option, increase the maxiter parameter.")
  }else{   
    fit$beta <- array(fit$beta, dim=c(p,K,ncol(fit$beta)))
    fit$eta <- array(fit$eta, dim=c(L,K,ncol(fit$eta)))
  }

  fit$tasks <- tasks
  fit$groups <- groups
  fit$model <- model
  fit$delta <- delta
  fit
}

grplasso<- function(X, y, groups, lambda=NULL, nlambda=20, model="linear", standardize=T,conv.eps=1e-3, eps=1e-6,maxiter=100,maxiter.shotgun=300) {

  n <- length(y)		        # replicates
  p <- ncol(X)			        # predictors
  L <- ncol(groups)			# groups
  tasks <- as.factor(rep("A",n))        # set tasks, all the same
 
  if (model == "linear") {
    model.num <- 0
  } else if(model=="logistic") {
    model.num <- 1
    y <- 1*(y==max(y))
  }

  if(standardize){
    if(model=="linear")
      y<-as.numeric(unlist(tapply(y,tasks,scale,scale=F,simplify=T)))
    X<-matrix(unlist(apply(X,2,tapply,tasks,scale)),ncol=p)	
  }

   if(is.null(lambda)){
    lambda <- calcLambda(X,y,tasks,nlambda,model)
  }
  lambda <- sort(lambda,decreasing=T)
  
  fit <- NULL; nlambda <- length(lambda)
   for(i in 1:nlambda){
     temp.fit <- .Call("grplasso", X, y, n, groups, lambda[i], model.num, conv.eps, eps, maxiter, maxiter.shotgun, PACKAGE = "multitask")
     if(temp.fit$converged){
       fit$converged <- c(fit$converged, temp.fit$converged) 
       fit$beta <- cbind(fit$beta, as.numeric(temp.fit$beta))
       fit$d <- cbind(fit$d, temp.fit$d)
       fit$bic <- c(fit$bic, as.numeric(temp.fit$bic)) 
       fit$lambda <- c(fit$lambda,lambda[i])
     }
   }

  if(is.null(fit)){
    warning("No lambdas converged. Try to give larger lambdas or, as a second option, increase the maxiter parameter.")
  }
  
  fit$groups <- groups
  fit$model <- model
  fit
}


calcLambda <- function(X,y,tasks,nlambda,model){
    lams <- c()
    n <- as.numeric(table(tasks))	
    K <- length(levels(tasks))
    p <- ncol(X)
    for(k in 1:K){
        task <- levels(tasks)[k]
        if(model=="linear"){
            resids <- y[tasks==task]
    }else if(model=="logistic"){
        resids <- y[tasks==task] - 0.5 #demands 0/1 for y
    } 
        lams <- c(lams,abs(drop(crossprod(X[tasks==task,,drop=F], resids))))
    }
    lambda.max <- max(lams)
    if(model=="linear")
        lambda.min <- 5*sqrt(max(n)+p)/10
    else
        lambda.min <- 2*sqrt(max(n)+p)/10
    lambda <- seq(lambda.max,lambda.min,length.out=nlambda)
    lambda
}
