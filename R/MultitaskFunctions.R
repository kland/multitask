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

  if(model=="logistic")
    y[y==0] <- -1 #format to fit shotgun input
  
  fit <- NULL; nlambda <- length(lambda)
  for(i in 1:nlambda){
    temp.fit <- .Call("multitask", X, y, K, groups, lambda[i], model.num, eps, maxiter, maxiter.shotgun, PACKAGE = "multitask")
    if(!temp.fit$converged){
      break; #skip rest of lambdas, since they are likely to not converge either
    }
    fit$converged <- c(fit$converged, temp.fit$converged) 
    fit$beta <- cbind(fit$beta, as.numeric(temp.fit$beta))
    fit$alpha <- cbind(fit$alpha, as.numeric(temp.fit$alpha))
    fit$eta <- cbind(fit$eta,as.numeric(temp.fit$eta))
    fit$d <- cbind(fit$d, temp.fit$d)
    fit$bic <- c(fit$bic, as.numeric(temp.fit$bic)) 
    fit$lambda <- c(fit$lambda,lambda[i])
  }

  if(is.null(fit)){
    warning("No lambdas converged. Try to give larger lambdas or, as a second option, increase the maxiter parameter.")
  }else{   
    fit$alpha <- array(fit$alpha, dim=c(p,K,ncol(fit$alpha)))
    fit$beta <- array(fit$beta, dim=c(p,K,ncol(fit$beta)))
    fit$eta <- array(fit$eta, dim=c(L,K,ncol(fit$eta)))
  }

  fit$tasks <- tasks
  fit$groups <- groups
  fit
}

calcLambda <- function(X,y,tasks,nlambda,model){
  lams <- c()
  for(k in 1:K){
    task <- levels(tasks)[k]
    if(model=="linear"){
      resids <- y[tasks==task]
    }else if(model=="logistic"){
      resids <- y[tasks==task] - 0.5 #demands 0/1 for y
    } 
    lams <- c(lams,abs(drop(crossprod(X[tasks==task,], resids))))
  }
  lambda.max <- max(lams)
  lambda.min <- 4*sqrt(max(n)+p)/10  
  lambda <- seq(lambda.max,lambda.min,length.out=nlambda)
  lambda
}

