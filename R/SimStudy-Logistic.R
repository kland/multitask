simstudy.logistic <- function () {

	# First test case with half of the genes active, and one task killed for one group.
  set.seed(213467)
	#set up parameters (small problem)
	n <- 100	# replicates
	p <- 21		# predictors
	L <- 4		# groups
	K <- 3		# tasks

	#set up groups for genes (predictors)
	groups <- matrix(0,nrow=p,ncol=L)
	rownames(groups)<-paste('gene',1:p,sep=' ');colnames(groups)<-paste('group',1:L,sep=' ')
	offset<-0
	for (g in 1:L){
		groups[(offset+1):min(offset+6,p),g] <- 1
		offset <- offset + 5}

	#generate X (feature) matrix
	X <- array(rnorm(n*p*K), dim=c(n,p,K))

	#generate true regression coefficients
	w <- matrix(0,nrow=p,ncol=K)
	w[groups[,1]+groups[,2]>0,]<-runif(K*sum(groups[,1]+groups[,2]>0),min=0.2,max=2.5)*sample(c(-1,1),size=K*sum(groups[,1]+groups[,2]>0),replace=T )
	w[1:5,3]<-0
	w[9,1]<-0

	colnames(w)<-paste('task ',1:K,sep=' ')
	rownames(w)<-paste('gene',1:p,sep=' ')
	beta.true<-w


	#generate bernoulli responses (logistic regression)
	y<-NULL
	for(k in 1:K){
		linpred <- X[,,k]%*%w[,k] 
		prob <- exp(linpred)/(1 + exp(linpred))
		y<-cbind(y,rbinom(n,1,prob))
	}	

	#reformat data
	y<-as.numeric(y)
	Xtemp<-NULL
	for(k in 1:K){
		Xtemp<-rbind(Xtemp, X[,,k])
	}
	X<-Xtemp
	
	#make task indicator vector
	tasks<-rep(LETTERS[1:K],each=n)

	#run method for a sequence of lambdas
	out <- multitask(X,y,tasks,groups,model="logistic",maxiter=100)
}
