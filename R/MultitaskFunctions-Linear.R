lasso <- function (X, y, lambda, model, positive, eps = 1e-12) {
  if(model=="linear") model <- 1
  else if(model=="logistic") model <- 2
  .Call("multitask_lasso", X, y, lambda, model, as.integer(positive), eps, PACKAGE = "multitask")
}

#lasso2 <- function (X, y, lambda, eps = 1e-12) {
#	.Call("multitask_lasso2", X, y, lambda, eps, PACKAGE = "multitask")
#}


x.tilde <- function (X, tasks, groups, d.cur, eta.cur, K, k) {
	.Call("multitask_x_tilde", X, tasks, groups, d.cur, eta.cur, K, k, PACKAGE = "multitask")
}


#solving the garotte problem 
solveGarotte.linear<-function(y,X,lambda=1,eps=1e-12){
  require(quadprog)
  # quadprog: solving quadratic programming problems of the
  # form min(-d^T b + 1/2 b^T D b) with the constraints A^T b >= b_0.

  n <- length(y) 	# n: number of total replicates	
  K <- ncol(X) 	# K: number of groups
	
  # check for zero or low eigenvalues, extra to check if system is solvable.
  dout<-rep(0,K)
  pick1<-apply(abs(X),2,sum)>eps
  X<-X[,pick1,drop=F]
  if (ncol(X)==0){
    return(dout)
  }
  pick2<-!abs(eigen(t(X) %*% X)$values)<eps
  X<-X[,pick2,drop=F]
	
  # K: new number of groups
  K <- ncol(X)
	
  if(is.null(dim(X)))
    return(dout)
	
  linvec <- t(t(y) %*% X - rep(lambda,K))
  quadmat <- t(X) %*% X
	
  A <- diag(K)
  d0<-rep(0,K)
	
  try(sol<-solve.QP(quadmat, linvec, A, d0)$solution)	
  try(sol[abs(sol)<eps]<-0)
  try(dout[pick1][pick2]<-sol)

  return(dout)
}

multitask.linear<-function(X,y,tasks,groups,lambda,eps=1e-12){

  require(penalized)
  require(glmnet)
		
# initial formatting
  y<-as.numeric(y)
  X<-data.matrix(X)
  tasks<-as.factor(tasks)
	
  n <- as.numeric(table(tasks))	# replicates
  p <- ncol(X)	                # predictors
  K <- length(levels(tasks))    # tasks
  L <- ncol(groups)		# groups

# select random starting points  
  dstart<-runif(L,0.5,5)
  alphastart<-matrix(runif(p*K,-0.1,0.1),nrow=p,ncol=K)
  etastart<-matrix(runif(L*K,0.5,5),nrow=L,ncol=K)
  betastart<-alphastart
		
# variables holding new and current estimates of parameters
  alpha.cur <- alphastart   # size: p x K 
  d.cur <-dstart            # size: L x 1 (vector of length L)
  eta.cur <-etastart        # size: L x K
  beta.cur <- betastart     # size: p x K 


# Here starts the loop that should be moved to C++
#
# Input arguments:
# X,y,tasks,groups,lambda
# alpha.cur,d.cur,eta.cur,beta.cur (i.e. starting values for the parameters, could also be assigned in C++)
#
# Return arguments:
# alpha.cur,d.cur,eta.cur,beta.cur (i.e. the final values for the parameters)
# 

  converged<-FALSE
  while(!converged){
    # 1. update alpha
    alpha.new<-matrix(NA,nrow=p,ncol=K)
    for(k in 1:K){
      task<-levels(tasks)[k]
      # this matrix is of size n x p. Corresponds to equation I (extra notation in paper).
      #Xtilde<-X[tasks==task,] %*% diag(apply(groups %*% diag(d.cur * eta.cur[,k]),1,sum))
      Xtilde <- x.tilde(X, tasks, groups, d.cur, eta.cur, K, k)
      # this is the call to the Lasso solver
      #alpha.fit<-penalized(y[tasks==task],Xtilde,unpenalized = ~0,lambda1=lambda,standardize=F,trace=F)
      #alpha.new[,k]<-coef(alpha.fit,"all")
      #TODO: Find out why result from shotgun lasso and penalized differ
      alpha.new[,k] <- lasso(Xtilde, y[tasks==task], lambda)
    }

#    # alternative, expanded matrix
#    Xtilde <- list()
#    for(k in 1:K){
#      task<-levels(tasks)[k]
#      Xtilde[[k]]<-(X[tasks==task,] %*% diag(apply(groups %*% diag(d.cur * eta.cur[,k]),1,sum)))
#    }
#    Xtilde <- as.matrix(bdiag(Xtilde))
#    #penalized
#    #alpha.fit<-penalized(y,Xtilde,unpenalized = ~0,lambda1=lambda,standardize=F,trace=F)
#    #alpha.new <- matrix(coef(alpha.fit, "all"), nrow = p, ncol = K)
#    #shotgun lasso
#    alpha.new <- matrix(lasso(Xtilde, y, lambda), nrow = p, ncol = K)
    
    # 2. update d
    # Xtilde2 is of size K*n x L after this loop. Corresponds to equation II (extra notation in paper).
    Xtilde2 <- NULL
    for(k in 1:K){
      task<-levels(tasks)[k]
       Xtilde2<-rbind(Xtilde2,((X[tasks==task,] %*% diag(alpha.new[,k])) %*% groups) %*% diag(eta.cur[,k]))
    }
    # this is the call to the quadprog solver
    d.new<-sign(solveGarotte.linear(y,Xtilde2))
 		
    # 3. update eta
    eta.new<-matrix(NA,nrow=L,ncol=K)
    for(k in 1:K){
      task<-levels(tasks)[k]
      # this matrix is of size n x L. Corresponds to equation III (extra notation in paper).
      Xtilde3<-((X[tasks==task,] %*% diag(alpha.new[,k])) %*% groups) %*% diag(d.new)
      eta.new[,k]<-solveGarotte.linear(y[tasks==task],Xtilde3)
    }
 
    # 4. update beta
    beta.new <- alpha.new * (groups %*% (diag(d.new) %*% eta.new))
    
    # check convergence
    if (max(abs(beta.new - beta.cur))<eps){
      converged<-TRUE
    }

    # update current estimates
    alpha.cur<-alpha.new
    d.cur<-d.new
    eta.cur<-eta.new
    beta.cur<-beta.new
  } #end of while loop
  
  fit<-NULL
  fit$beta <- beta.cur
  fit$alpha <- alpha.cur
  fit$d <- d.cur
  fit$eta <- eta.cur
  fit$tasks <- tasks
  fit$groups <- groups
  fit$lambda <- lambda
  fit
}
