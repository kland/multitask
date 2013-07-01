lasso <- function (X, y, lambda, model, positive, eps = 1e-12) {
  if (model == "linear") {
    model <- 0
    lambda <- lambda*2
  } else if(model=="logistic") {
    model <- 1
  }
  .Call("multitask_lasso", X, y, lambda, model, as.integer(positive), eps, PACKAGE = "multitask")
}

x.tilde <- function (X, tasks, groups, d.cur, eta.cur, K) {
  .Call("multitask_x_tilde", X, tasks, groups, d.cur, eta.cur, K, PACKAGE = "multitask")
}

x.tilde.2 <- function (X, tasks, groups, alpha.new, eta.cur, K) {
  .Call("multitask_x_tilde_2", X, tasks, groups, alpha.new, eta.cur, K, PACKAGE = "multitask")
}

x.tilde.3 <- function (X, tasks, groups, alpha.new, d.new, K, k) {
  .Call("multitask_x_tilde_3", X, tasks, groups, alpha.new, d.new, K, k, PACKAGE = "multitask")
}

Beta.new <- function (groups, alpha.new, d.new, eta.new, K) {
  .Call("multitask_beta_new", groups, alpha.new, d.new, eta.new, K, PACKAGE = "multitask")
}

Bic <- function (X, y, beta.new, eps, n) {
  .Call("multitask_bic", X, y, beta.new, eps, n, PACKAGE = "multitask")
}

multitask<-function(X,y,tasks,groups,lambda,model="linear",eps=1e-12){
    
# initial formatting
  y<-as.numeric(y)
  X<-data.matrix(X)
  tasks<-as.factor(tasks)
  
  n <- as.numeric(table(tasks))  # replicates
  p <- ncol(X)                  # predictors
  K <- length(levels(tasks))    # tasks
  L <- ncol(groups)    # groups

# select random starting points  
  dstart<- rep(1,L)
  alphastart<-matrix(runif(p*K,-0.1,0.1),nrow=p,ncol=K)
  etastart<-matrix(1,nrow=L,ncol=K)
  betastart<-alphastart

# variables holding new and current estimates of parameters
  alpha.cur <- alphastart   # size: p x K 
  d.cur <-dstart            # size: L x 1 (vector of length L)
  eta.cur <-etastart        # size: L x K
  beta.cur <- betastart     # size: p x K 
  bic<-NULL

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
    Xtilde <- x.tilde(X, tasks, groups, d.cur, eta.cur, K)
    alpha.new <- matrix(lasso(Xtilde, y, lambda, model="linear", positive=F), nrow = p, ncol = K)
    
    # 2. update d
    Xtilde2 <- x.tilde.2(X, tasks, groups, alpha.new, eta.cur, K)
    d.new <- sign(lasso(Xtilde2,y,lambda=1,model=model,positive=T))
     
    # 3. update eta
    eta.new<-matrix(NA,nrow=L,ncol=K)
    for(k in 1:K){
      task<-levels(tasks)[k]
      # this matrix is of size n x L. Corresponds to equation III (extra notation in paper).
      Xtilde3 <- x.tilde.3(X, tasks, groups, alpha.new, d.new, K, k)
      eta.new[,k]<-lasso(Xtilde3,y[tasks==task],lambda=1,model=model,positive=T)
    }
 
    # 4. update beta
    beta.new <- Beta.new(groups, alpha.new, d.new, eta.new, K)

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

  ## calculate BIC
  bic <- Bic(X, y, beta.new, eps, n[1])
  
  fit<-NULL
  fit$beta <- beta.cur
  fit$alpha <- alpha.cur
  fit$d <- d.cur
  fit$eta <- eta.cur
  fit$tasks <- tasks
  fit$groups <- groups
  fit$lambda <- lambda
  fit$bic<-bic
  fit
}
