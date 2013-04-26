

#solving the garotte problem 
solveGarotte.logistic<-function(y,X,lambda=1,eps=1e-12){
	require(penalized)
	
	# n: number of replicates
	n <- length(y)
	
	# K: number of groups
	K <- ncol(X)
	
	#check for zero or low eigenvalues
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
		
	try(sol<-penalized(y,X,unpenalized = ~0,lambda1=lambda,standardize=F,trace=F,model="logistic",positive=T))
	try(sol <- coef(sol,"all"))
	
	try(sol[abs(sol)<eps]<-0)
	try(dout[pick1][pick2]<-sol)
	
	return(dout)
}

multitask.logistic<-function(X,y,tasks,groups,steps=10,lambda.range=c(1,0),lambda.val=NULL,standardize=T,verbose=T,plot.it=F,eps=1e-5,warm.start=T){
	require(penalized)
	require(glmnet)
		
	#initial formatting
	y<-as.numeric(y)
	X<-data.matrix(X)
	tasks<-as.factor(tasks)
	
	n <- as.numeric(table(tasks))		# replicates
	p <- ncol(X)						# predictors
	K <- length(levels(tasks))			# tasks
	L <- ncol(groups)					# groups
	if(verbose)
		cat(paste('There are', K ,'tasks,',p ,'predictors, and', L ,'groups.\n'))
	
	# 0.1 standardize
	if(standardize){
		y<-as.numeric(unlist(tapply(y,tasks,scale,scale=F,simplify=T)))
		X<-matrix(unlist(apply(X,2,tapply,tasks,scale)),ncol=p)	
	}
	
	# 0.2 initialize
#	dstart<-rep(1,L)
	dstart<-runif(L,0.5,5)
	alphastart<-NULL
#	for(level in levels(tasks)){
#		gfit<-glmnet(X[tasks==level,],y[tasks==level],standardize=F,family='gaussian')
#		cvfit<-cv.glmnet(X[tasks==level,],y[tasks==level],standardize=F,family='gaussian')
#		alphastart<-cbind(alphastart,as.matrix(predict(gfit,s=cvfit$lambda.min,type='coefficients'))[-1])
#		optfit<-optL1(y[tasks==level],X[tasks==level,],unpenalized = ~0,model='linear',standardize=F,trace=F)
#		#glmnetfit<-cv.glmnet(y[tasks==level],X[tasks==level,],standardize=F,family='gaussian')
#		alphastart<-cbind(alphastart,optfit$fullfit@penalized)
#	}
	alphastart<-matrix(runif(p*K,-2,2),nrow=p,ncol=K)
#	etastart <- matrix(1,nrow=L,ncol=K)
	etastart<-matrix(runif(L*K,0.5,5),nrow=L,ncol=K)
	betastart<-alphastart
	
	# set lambda seq
	init.fit<-penalized(y,X,lambda1=0,unpenalized = ~0,steps=steps,standardize=F,trace=F)
	min.lambda <- Inf
	max.lambda <- -Inf
	for(i in 1:length(init.fit)){
		if (init.fit[[i]]@converged){
			max.lambda<-max(max.lambda,init.fit[[i]]@lambda1)
			min.lambda<-min(init.fit[[i]]@lambda1,min.lambda)
		}
	}
	min.lambda<-max(min.lambda,0.1)
	if(is.null(lambda.val)){
		#lambda.seq<-seq((init.fit[[1]]@lambda1)*lambda.range[1],max(init.fit[[length(init.fit)]]@lambda1*lambda.range[2],1),length.out=steps)
		lambda.range2<-approx(c(1,0),c(max.lambda,min.lambda),xout=lambda.range)$y
		lambda.seq<-seq(lambda.range2[1],lambda.range2[2],length.out=steps)		
	}else{
		lambda.seq=approx(c(1,0),c(max.lambda,min.lambda),xout=lambda.val)$y
	}
	if(verbose)
		cat(paste("Using lambda sequence of length ",length(lambda.seq),".\n",sep=""))
	
	if(verbose){
		cat('Starting values..\n')
		print(dstart)
		print(etastart)
	}
	
	# initialize matrices
	beta.mat<-NULL
	d.mat<-NULL
	alpha.mat<-NULL
	eta.mat<-NULL
	bic<-NULL
	
	alpha.new <- alphastart
	d.new <- dstart
	eta.new <- etastart
	beta.new <- betastart
	

	# estimation for each lambda
	for(lambda in lambda.seq){
		if(warm.start){
			alpha <- as.numeric(alpha.new)
			beta <- as.numeric(beta.new)
			eta<-as.numeric(eta.new)
			d<-d.new
			eta[eta==0]<-1
			d[d==0]<-1
			d.cur<-d 
			eta.cur<-eta.new 
			eta.cur[eta.cur==0]<-1
			beta.cur<-beta.new
			alpha.cur<-alpha.new
		}else{
			alpha <- as.numeric(alphastart)
			d <- dstart
			eta <- etastart
			beta <- as.numeric(betastart)
			beta.cur <- betastart
			alpha.cur <- alphastart
			d.cur <- dstart
			eta.cur<-etastart
			eta<-as.numeric(eta)
		}
		
		converged<-FALSE
		while(!converged){
		# 1. update alpha
			alpha.new<-matrix(NA,nrow=p,ncol=K)
			for(k in 1:K){
				task<-levels(tasks)[k]
				Xtilde<-X[tasks==task,] %*% diag(apply(groups %*% diag(d.cur * eta.cur[,k]),1,sum))
				alpha.fit<-penalized(y[tasks==task],Xtilde,unpenalized = ~0,lambda1=lambda,standardize=F,trace=F,model="logistic")
				alpha.new[,k]<-coef(alpha.fit,"all")
			}
			alpha<-cbind(alpha,as.numeric(alpha.new))
		
		# 2. update d
			Xtilde2<-NULL
			for(k in 1:K){
				task<-levels(tasks)[k]
				Xtilde2<-rbind(Xtilde2,((X[tasks==task,] %*% diag(alpha.new[,k])) %*% groups) %*% diag(eta.cur[,k]))
			}
			d.new<-solveGarotte.logistic(y,Xtilde2)
#d.new[d.new>0]<-1
			d<-cbind(d,d.new)
		
		# 3. update eta
			eta.new<-matrix(NA,nrow=L,ncol=K)
			for(k in 1:K){
				task<-levels(tasks)[k]
				Xtilde2<-((X[tasks==task,] %*% diag(alpha.new[,k])) %*% groups) %*% diag(d.new)
				eta.new[,k]<-solveGarotte.logistic(y[tasks==task],Xtilde2)
			}
#			eta.new[eta.new>0]<-1
			eta<-cbind(eta,as.numeric(eta.new))
		# 4. update beta
			beta.new<- alpha.new * (groups %*% (diag(d.new) %*% eta.new))
			beta<-cbind(beta,as.numeric(beta.new))
		
		# check convergence
			if (max(abs(beta.new - beta.cur))<eps){
#			if (max(abs(c(as.numeric(alpha.new-alpha.cur), d.new-d.cur, as.numeric(eta.new-eta.cur))))<eps){
				if(plot.it){
					layout(matrix(c(1:4),2,2))
					matplot(t(matrix(c(1:ncol(beta)),nrow=nrow(beta),ncol=ncol(beta),byrow=T)),t(beta),type='l',col='green',xlab='iteration',ylab='beta')
					matplot(t(matrix(c(1:ncol(alpha)),nrow=nrow(alpha),ncol=ncol(alpha),byrow=T)),t(alpha),type='l',col='red',xlab='iteration',ylab='alpha')
					matplot(t(matrix(c(1:ncol(d)),nrow=nrow(d),ncol=ncol(d),byrow=T)),t(d),type='l',col='black',xlab='iteration',ylab='d')
					matplot(t(matrix(c(1:ncol(eta)),nrow=nrow(eta),ncol=ncol(eta),byrow=T)),t(eta),type='l',col='blue',xlab='iteration',ylab='eta')
					mtext(paste("lambda (normalized) =", round(lambda/max(lambda.seq),2)), side = 3, line = -3, outer = TRUE,cex=1.2,font=2) 
					}
				
				converged<-TRUE
				alpha.mat<-cbind(alpha.mat,as.numeric(alpha.new))
				beta.mat<-cbind(beta.mat,as.numeric(beta.new))
				d.mat<-cbind(d.mat,d.new)
				eta.mat<-cbind(eta.mat,as.numeric(eta.new))
				#calculate BIC
				df<-sum(abs(beta.new)>eps)
				ll<-NULL
				for(k in 1:K){
					task<-levels(tasks)[k]
					lp <- X[tasks==task,] %*% beta.new[,k]
					ll <- sum(ll,sum(y[tasks==task]*lp - log(1+exp(lp))))
				}
				bic<- c(bic,-2*ll + df*log(sum(n)))
				if(verbose)
					cat(paste('Estimation for lambda =', lambda,'done.\n'))
			}
		# update current estimates
			alpha.cur<-alpha.new
			d.cur<-d.new
			eta.cur<-eta.new
			beta.cur<-beta.new
		}
	}
	fit<-NULL
	fit$lambda<-lambda.seq
	if(is.null(lambda.val)){
		fit$normlambda<-seq(lambda.range[1],lambda.range[2],length.out=steps)
	}else{
		fit$normlambda<-lambda.val
	}
	fit$bic<-bic
	fit$beta<-beta.mat
	fit$alpha<-alpha.mat
	fit$d<-d.mat
	fit$eta<-eta.mat
	fit$tasks<-tasks
	fit
}

lambda.interpolate<-function (lambda, s) {
    if (length(lambda) == 1) {
        nums = length(s)
        left = rep(1, nums)
        right = left
        sfrac = rep(1, nums)
    }
    else {
        s[s > max(lambda)] = max(lambda)
        s[s < min(lambda)] = min(lambda)
        k = length(lambda)
        sfrac <- (lambda[1] - s)/(lambda[1] - lambda[k])
        lambda <- (lambda[1] - lambda)/(lambda[1] - lambda[k])
        coord <- approx(lambda, seq(lambda), sfrac)$y
        left <- floor(coord)
        right <- ceiling(coord)
        sfrac = (sfrac - lambda[right])/(lambda[left] - lambda[right])
        sfrac[left == right] = 1
    }
    list(left = left, right = right, frac = sfrac)
}

predict.multi<-function (fit, newx=NULL, newtasks=NULL,lambda, type = c("fit", "coefficients")){ 
	K<-length(levels(fit$tasks))
	p<-ncol(newx)
	
    type <- match.arg(type)
    if (is.null(newx) & type == "fit") {
        warning("Type=fit with no newx argument; type switched to coefficients")
        type <- "coefficients"}
	if(type=='coefficients')
		return(fit$beta)
			
	newbeta<-fit$beta

	if(missing(lambda)){			
		lambda<-fit$normlambda}else
	{
		lamlist<-lambda.interpolate(fit$normlambda,lambda)
		newbeta<-newbeta[,lamlist$left,drop=FALSE]*lamlist$frac +newbeta[,lamlist$right,drop=FALSE]*(1-lamlist$frac)
	}
	fitted.values<-NULL
	beta.tasks<-rep(levels(fit$tasks),each=p)
	for(k in 1:K){
		task <- levels(fit$tasks)[k]
		fitted.values<-rbind(fitted.values,newx[newtasks==task,] %*% newbeta[beta.tasks==task,])
	}
	return(fitted.values)
}
		

createFolds <- function (y, k = 10, list = TRUE, returnTrain = FALSE) 
{
    if (is.numeric(y)) {
        y <- cut(y, unique(quantile(y, probs = seq(0, 1, length = min(5, length(y))))), include.lowest = TRUE)
    }
    if (k < length(y)) {
        y <- factor(y)
        numInClass <- table(y)
        foldVector <- vector(mode = "integer", length(y))
        for (i in 1:length(numInClass)) {
            seqVector <- rep(1:k, numInClass[i]%/%k)
            if (numInClass[i]%%k > 0) 
			seqVector <- c(seqVector, sample(1:k, numInClass[i]%%k))
            foldVector[which(y == dimnames(numInClass)$y[i])] <- sample(seqVector)
        }
    }
    else foldVector <- seq(along = y)
    if (list) {
        out <- split(seq(along = y), foldVector)
        names(out) <- paste("Fold", gsub(" ", "0", format(seq(along = out))), 
							sep = "")
        if (returnTrain) 
		out <- lapply(out, function(data, y) y[-data], y = seq(along = y))
    }
    else out <- foldVector
    out
}


cv.multi<-function(X,y,tasks,groups,nfolds = 10,steps=10,lambda.range=c(1,0)){
	all.folds <- createFolds(tasks, k=nfolds)
    residmat <- matrix(0, steps, nfolds)
    for (i in seq(nfolds)) {
        omit <- all.folds[[i]]
        fit <- multitask(X[-omit, , drop = FALSE], y[-omit],tasks[-omit],groups,steps=steps,lambda.range=lambda.range,verbose=F)
        pred <- predict.multi(fit, X[omit, , drop = FALSE], tasks[omit],lambda=fit$normlambda,type="fit")
		if (length(omit) == 1) 
			fit <- matrix(fit, nrow = 1)
        residmat[, i] <- apply((y[omit] - pred)^2, 2, mean)
	}
	cv <- apply(residmat, 1, mean)
    cv.error <- sqrt(apply(residmat, 1, var)/nfolds)
	object <- list(lambda=fit$normlambda, cv = cv, cv.error = cv.error)
	return(object)
}
	


#beta.est<-data.matrix(matrix(beta.mat[,which.min(bic)],ncol=K))
#colnames(beta.est)<-paste('task ',1:K,sep=' ')
#rownames(beta.est)<-paste('gene',1:p,sep=' ')
#	beta.true<-w; dimnames(beta.true)<-dimnames(beta.est)
#alpha.est<-data.matrix(matrix(alpha.mat[,which.min(bic)],ncol=K))
#dimnames(alpha.est)<-dimnames(beta.est)
#d.est<-data.matrix(matrix(d.mat[,which.min(bic)],ncol=1))
#rownames(d.est)<-colnames(groups)
#colnames(d.est)<-'Groups'
#eta.est<-data.matrix(matrix(eta.mat[,which.min(bic)],ncol=K))
#colnames(eta.est)<-colnames(beta.est)
#rownames(eta.est)<-colnames(groups)


#alpha <- as.numeric(alpha.new)
#beta <- as.numeric(beta.new)
#eta<-as.numeric(eta.new)
#d<-d.new
#eta[eta==0]<-0.01
#d[d==0]<-0.01

# plotting function
myImagePlot2 <- function(x,color='color',key=T,plot.range=NULL,ylab=T,...){
	min <- min(x)
	max <- max(x)
	yLabels <- rownames(x)
	xLabels <- colnames(x)
	title <-c()
	# check for additional function arguments
	if(key)
		layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(4,1), heights=c(1,1))
	
	# Set color scheme
	if (color=='color'){
		ColorRamp<-rev(heat.colors(256))}
	if (color=='grey'){
		ColorRamp <- grey(seq(1,0,length.out=256))}
	
	if(!is.null(plot.range)){
		min<-plot.range[1];max<-plot.range[2]}
	ColorLevels <- seq(min, max, length=length(ColorRamp))
	
	
	
	# Reverse Y axis
	reverse <- nrow(x) : 1
	yLabels <- yLabels[reverse]
	x <- x[reverse,]
	
	# Data Map
	if(ylab)
		par(mar = c(4,5,2.5,2))
	else
		par(mar = c(4,2,2.5,2))
	image(1:length(xLabels), 1:length(yLabels), t(x), col=ColorRamp, xlab="",
		  ylab="", axes=FALSE, zlim=c(min,max))
	if( !is.null(title) ){
		title(main=title)
	}
	box()
	axis(BELOW<-1, at=1:length(xLabels), labels=xLabels, cex.axis=1,las=2)
	if(ylab){
		axis(LEFT <-2, at=1:length(yLabels), labels=yLabels, las= HORIZONTAL<-1,cex.axis=1)
	}
	axis(TOP <-3, at=c(2,5), labels=c(expression(abs(beta)),expression(abs(hat(beta)-beta))), las= HORIZONTAL<-1,cex.axis=1,tick=F,line=-0.5)
	abline(v=3.5,lwd=2)
	# Color Scale
	if(key){
		par(mar = c(3,2.5,2.5,2))
		image(1, ColorLevels,
			  matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
			  col=ColorRamp,
			  xlab="",ylab="",
			  xaxt="n")
	}
	layout(1)
}




myImagePlotAll <- function(xlst,color='color',key=T,plot.range=NULL,ylab=T,...){
	
	layout(matrix(data=c(1,2,3,7,4,5,6,7), nrow=2, ncol=4,byrow=T), widths=c(4,4,4,1,4,4,4,1),heights=c(1,1,1,0.5,1,1,1,0.5))
#	layout(matrix(data=c(1,2,7,3,4,7,5,6,7), nrow=3, ncol=3,byrow=T), widths=c(4,4,1,4,4,1,4,4,1))
	
	for(i in 1:length(xlst)){
		x<-xlst[[i]]
		min <- min(x)
		max <- max(x)
		yLabels <- rownames(x)
		xLabels <- colnames(x)
		title <-c()
# check for additional function arguments
		
# Set color scheme
		if (color=='color'){
			ColorRamp<-rev(heat.colors(256))}
		if (color=='grey'){
			ColorRamp <- grey(seq(1,0,length.out=256))}
	
		if(!is.null(plot.range)){
			min<-plot.range[1];max<-plot.range[2]}
		ColorLevels <- seq(min, max, length=length(ColorRamp))
	
	
	
# Reverse Y axis
		reverse <- nrow(x) : 1
		yLabels <- yLabels[reverse]
		x <- x[reverse,]
	
# Data Map
		if(ylab)
		par(mar = c(4,5,2.5,2))
		else
		par(mar = c(4,2,2.5,2))
		image(1:length(xLabels), 1:length(yLabels), t(x), col=ColorRamp, xlab="",
		  ylab="", axes=FALSE, zlim=c(min,max))
		if( !is.null(title) ){
			title(main=title)
		}
		box()
		axis(BELOW<-1, at=1:length(xLabels), labels=xLabels, cex.axis=1,las=2)
		if(ylab){
			axis(LEFT <-2, at=1:length(yLabels), labels=yLabels, las= HORIZONTAL<-1,cex.axis=1)
		}
		if(i==3){
			axis(TOP <-3, at=c(10.5,30.5), labels=c(expression(abs(phantom()~beta~phantom())),expression(abs(phantom()~hat(beta)-beta~phantom()))), las= HORIZONTAL<-1,cex.axis=1,tick=F,line=-0.5,cex.axis=1.3)
			abline(v=20.5,lwd=2)
		}else{
			
			axis(TOP <-3, at=c(2,5), labels=c(expression(abs(phantom()~beta~phantom())),expression(abs(phantom()~hat(beta)-beta~phantom()))), las= HORIZONTAL<-1,cex.axis=1,tick=F,line=-0.5,cex.axis=1.3)
			abline(v=3.5,lwd=2)
		}	
## Color Scale
	}
		if(key){
			par(mar = c(3,2.5,2.5,2))
			image(1, ColorLevels,
			  matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
			  col=ColorRamp,
			  xlab="",ylab="",
			  xaxt="n")
		}
		layout(1)
}




myImagePlot <- function(x,color='color',key=T,plot.range=NULL,...){
	min <- min(x)
	max <- max(x)
	yLabels <- rownames(x)
	xLabels <- colnames(x)
	title <-c()
# check for additional function arguments
	if(key)
	layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(4,1), heights=c(1,1))
	
# Set color scheme
	if (color=='color'){
		ColorRamp<-rev(heat.colors(256))}
	if (color=='grey'){
		ColorRamp <- grey(seq(1,0,length.out=256))}
	
	if(!is.null(plot.range)){
		min<-plot.range[1];max<-plot.range[2]}
	ColorLevels <- seq(min, max, length=length(ColorRamp))
	
	
	
# Reverse Y axis
	reverse <- nrow(x) : 1
	yLabels <- yLabels[reverse]
	x <- x[reverse,]
	
# Data Map
	par(mar = c(4,5,2.5,2))
	image(1:length(xLabels), 1:length(yLabels), t(x), col=ColorRamp, xlab="",
		  ylab="", axes=FALSE, zlim=c(min,max))
	if( !is.null(title) ){
		title(main=title)
	}
	box()
	axis(BELOW<-1, at=1:length(xLabels), labels=xLabels, cex.axis=1,las=2)
	axis(LEFT <-2, at=1:length(yLabels), labels=yLabels, las= HORIZONTAL<-1,cex.axis=1)
	# Color Scale
	if(key){
		par(mar = c(3,2.5,2.5,2))
		image(1, ColorLevels,
			  matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
			  col=ColorRamp,
			  xlab="",ylab="",
			  xaxt="n")
	}
	layout(1)
}


myImagePlot4 <- function(x,color='color',key=T,plot.range=NULL,main=NULL,...){
	min <- min(x)
	max <- max(x)
	yLabels <- rownames(x)
	xLabels <- colnames(x)
	
# Set color scheme
	if (color=='color'){
		ColorRamp<-rev(heat.colors(256))}
	if (color=='grey'){
		ColorRamp <- grey(seq(1,0,length.out=256))}
	
	if(!is.null(plot.range)){
		min<-plot.range[1];max<-plot.range[2]}
	ColorLevels <- seq(min, max, length=length(ColorRamp))
	
	
	
# Reverse Y axis
	reverse <- nrow(x) : 1
	yLabels <- yLabels[reverse]
	x <- x[reverse,]
	
# Data Map
	par(mar = c(4,5,2.5,2))
	image(1:length(xLabels), 1:length(yLabels), t(x), col=ColorRamp, xlab="",
		  ylab="", axes=FALSE, zlim=c(min,max))
	if( !is.null(main) ){
		title(main=main,cex.main=1.3)
	}
	box()
	axis(BELOW<-1, at=1:length(xLabels), labels=xLabels, cex.axis=1,las=2)
	axis(LEFT <-2, at=1:length(yLabels), labels=yLabels, las= HORIZONTAL<-1,cex.axis=1)
# Color Scale
	if(key){
		par(mar = c(3,2.5,2.5,2))
		image(1, ColorLevels,
			  matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
			  col=ColorRamp,
			  xlab="",ylab="",
			  xaxt="n")
	}
#	layout(1)
}


