
################################################################
#                                                              #
#                Aalen Additive Mdeol (R/C++)                  #
#                                                              #
################################################################
# Author: Ali Rafei (arafei {AT} umich {DOT} edu)
# Date:   03/02/2018
# The following codes are developed to fit a fully non-paraetric
# Aalen additive Model 
# Method can be either NULL or c("TSL" "JK1")

library(devtools)
library(inline)
library(Rcpp)
library(RcppArmadillo)

sourceCpp("aamC1.cpp")

aam1 <- function(time, censor=NULL, X, weight=NULL, cluster=NULL, strata=NULL, method="no"){
	X <- data.frame(X)
	na <- names(X)
	n <- length(time)
	r <- ncol(X)

	if(is.null(censor)){
		censor <- rep(1, n)
	}
	if(is.null(weight)){
		weight <- rep(1, n)
	}else{
		weight <- weight/sum(weight)
	}
	if(is.null(cluster)){
		cluster <- 1:n
	}
	if(is.null(strata)){
		strata <- rep(1, n)
	}

	# Sorting variables based on time
	so <- order(time)
	time <- time[so]
	censor <- censor[so]
	X <- as.matrix(cbind(rep(1, nrow(X)), X[so, ]))
	weight <- weight[so]
	cluster <- cluster[so]
	strata <- strata[so]

	strata <- as.numeric(as.factor(strata))
	for(h in unique(strata)){
		cluster[strata==h] <- as.numeric(as.factor(cluster[strata==h]))
	}

	outC <- aamC1(time=time, CR=censor, DM=X, WT=weight, ST=strata, CL=cluster, MT=method)
	Ast <- outC$AstX
	Ost <- outC$OstX


	# Extract teh variance of parameters estimates
	V <- matrix(0, n, r+1)
	for(i in 1:n){
		V[i, ] <- diag(Ost[, , i])
	}
	
	time=c(0, time[censor==1])
	cum=t(cbind(rep(0, r+1), Ast[, censor==1]))
	colnames(cum) <- c("Intercept", na)
	V <- rbind(rep(0, r+1), V[censor==1, ])
	colnames(V) <- paste0("V", 0:r)
	list(cum=cbind(time, cum), var.cum=cbind(time, V))
}



plot.aam <- function(object, level=0.05){
	time <- object$cum[, 1]
	B <- object$cum[, -1]
	V <- object$var.cum[, -1]
	r <- ncol(B)
	c.alpha <- qnorm(1 - level/2)
	for (i in 1:r){
		ll <- B[, i] - c.alpha * V[, i]^0.5
		ul <- B[, i] + c.alpha * V[, i]^0.5
		plot(time, B[, i], type="s", lwd=2, main=colnames(B)[i], xlim=c(0, max(time)), ylim=c(min(ll), max(ul)), xlab="time", ylab="cumulative coefficient")
		par(new=T)
		plot(time, ll, type="s", lty=2, main="", xlim=c(0, max(time)), ylim=c(min(ll), max(ul)), xlab="", ylab="")
		par(new=T)
		plot(time, ul, type="s", lty=2, main="", xlim=c(0, max(time)), ylim=c(min(ll), max(ul)), xlab="", ylab="")
		abline(a=0, b=0, col="red")
	}
}

smp.eval <- function(obj, coef, level=0.05){
	c.alpha <- qnorm(1 - level/2)
	Beta <- obj$cum[-1, -1]
	Var <- obj$var.cum[-1, -1]
	ll <- Beta - c.alpha * Var^0.5
	ul <- Beta + c.alpha * Var^0.5
	out <- matrix(0, 3, ncol(Beta))
	out[1, ] <- apply(Beta - coef, 2, function(x)mean(x, na.rm=T))
	out[2, ] <- apply((Beta - coef)^2, 2, function(x)mean(x, na.rm=T))
	out[3, ] <- apply((coef < ul) & (coef > ll), 2, function(x)mean(x, na.rm=T))
	as.numeric(t(out))
}

# A function to draw a stratified two-stage cluster sample with PPS
cpxsmp<- function(strata=NULL, PSU=NULL, PSU.prob=NULL, SSU=NULL, n_psu, ssu_size=NULL){
	if(is.null(PSU))
		strata <- strata
	if(is.null(strata))
		strata <- rep(1, length(PSU))
	if(is.null(SSU))
		SSU <- 1:length(strata)
	if(is.null(n_psu))
		n_psu <- 1
	if(is.null(ssu_size))
		ssu_size <- max(table(strata, PSU))		
	if(is.null(PSU.prob))
		PSU.prob <- rep(1, length(strata))
	ID <- 1:length(strata)
	H <- length(unique(strata))
	smp_cpx <- c()
	for(h in 1:H){
		psu_h <- unique(PSU[strata==h])
		psu_size <- (PSU.prob[strata==h])[!duplicated(PSU[strata==h])]
		psu_id <- psu_h[ppss(psu_size, n_psu)]
		for(j in psu_id){
			tmp <- ID[strata==h & PSU==j]
			smp_cpx <- c(smp_cpx, tmp[sample(1:length(tmp), ssu_size, replace=F)])
		}
	}
	smp_cpx
}
