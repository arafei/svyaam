
rm(list=ls())

# Required packages:
library(MASS)
library(pps)
library(survey)
library(survival)
library(timereg)

setwd("C:\\Users\\arafei\\Desktop\\AAM_analysis\\New")

#idx <- 461:480

#-------------------------------------------------------------------------------
#input:   wtTD is a matrix of (n0.samp) X (dim.parameters), i.e.
#             weighted Taylor deviate 
#         strat is a vector of strata id
#         psu is a vector of psu id
# output:  cs is the variance considering complex sampling
#3/6 revise the program to check each stratum with at least 2 psu's 
#---delete Strata with only one psu for variance estimation
#-------------------------------------------------------------------------------

################################################################
#                                                              #
#                        Aalen Additive Mdeol                  #
#                                                              #
################################################################
# Author: Ali Rafei (arafei {AT} umich {DOT} edu)
# Date:   03/02/2018
# The following codes are developed to fit a fully non-paraetric
# Aalen additive Model 
# Method can be either NULL or c("TSL" "JK1")

source("AAM-Rcpp4.R")


################################################################
#                                                              #
#                        SIMULATION STUDY                      #
#                                                              #
################################################################
# Author: Ali Rafei (arafei {AT} umich {DOT} edu)
# Date:   02/08/2018
# Date:   02/10/2018  Yan Li comments at lines 56, 63, 67, 75


# The following codes are developed to design a simulation study 
#for evaluating the performance of adjustments for complex samples 
# in Aalen's additive regression models

####################################################################

# imput simulation parameters"
# nsim: Number of iterations
# N: population size
# n0: sample size for probability sample
# n1: sample size for nonprobability sample
# JKN: number of replication weights based on probability/nonprobability sample

nSim <- 1000
N <- 1000000
n <- 1000

#model parameters in the population

# generating covariates X_1 and X_2 in the population
cr <- 0.3 #censoring rate
m <- 0 #X1: mean
s <- 0.5 #X1: sd
p <- 0.5 #X2: p(X2=1)
pop <- data.frame(C=rep(NA, N))

# Function to solve H(t)
xtime=function(x){
	u=uniroot(function(t){t+x[2]*log(1+t)+x[3]*(t^2)/2-x[1]},c(0, 50))
u$root
}

# Censoring times
pop$C <- rexp(N, cr)
pop$X1 <- rnorm(N, m, s)
pop$X2 <- rbinom(N, 1, p)

pop$T <- rep(0, N)
tmp <- cbind(rexp(N), pop$X1, pop$X2)
pop$T <- apply(tmp, 1, xtime)
pop$T[pop$T==0] <- min(pop$T[pop$T!=0])

pop$t <- pmin(pop$T, pop$C)
pop$delta <- (pop$T<pop$C)

# generate true population cumulative regression function
pop$b00 <- pop$T
pop$b01 <- log(1+pop$T)
pop$b02 <- (pop$T^2)/2

pop <- pop[order(pop$T), ]
par(mfrow=c(2, 2))
plot(pop$T, pop$b00, type="l", col="red")
plot(pop$T, pop$b01, type="l", col="red")
plot(pop$T, pop$b02, type="l", col="red")


pop$d <- rep(1, N)
aalen_pop <- aalen(formula = Surv(T, d) ~ X1+X2, data = pop, robust=0)
pop$b10 <- rep(NA, N)
pop$b11 <- rep(NA, N)
pop$b12 <- rep(NA, N)
pop[, c("b10", "b11", "b12")] <- aalen_pop$cum[-1, -1]

pdf("pop_plot.pdf", width=20, height=15)
par(mfrow=c(2, 2))
plot(survfit(Surv(t, delta)~1, data=pop), col="blue", conf.int=F, xlab="", ylab="survival", main="S(t)")

plot(pop$T, pop$b10, type="s", col="blue", ylim=c(0, 5), xlim=c(0, 5), xlab="", ylab="Cum. reg. function", main="A0(t)")
par(new=T)
plot(pop$T, pop$b00, type="l", col="red", lty=2, ylim=c(0, 5), xlim=c(0, 5), xlab="", ylab="")
legend(0, 5, c("Fitted CRF", "True CRF"), col=c("blue", "red"), lty=c(1:2), pch=16)

plot(pop$T, pop$b11, type="s", col="blue", ylim=c(0, 2), xlim=c(0, 5), xlab="time", ylab="Cum. reg. function", main="A1(t)")
par(new=T)
plot(pop$T, pop$b01, type="l", col="red", lty=2, ylim=c(0, 2), xlim=c(0, 5), xlab="", ylab="")
legend(0, 2, c("Fitted CRF", "True CRF"), col=c("blue", "red"), lty=c(1:2), pch=16)

plot(pop$T, pop$b12, type="s", col="blue", ylim=c(0, 10), xlim=c(0, 5), xlab="time", ylab="Cum. reg. function", main="A2(t)")
par(new=T)
plot(pop$T, pop$b02, type="l", col="red", lty=2, ylim=c(0, 10), xlim=c(0, 5), xlab="", ylab="")
legend(0, 10, c("Fitted CRF", "True CRF"), col=c("blue", "red"), lty=c(1:2), pch=16)
dev.off()

##########################################################################################

# Read population data
write.csv(pop, "sim_pop.csv", col.names=T, sep=",")
#pop <- read.csv("sim_pop.csv", header=T, sep=",")
#pop$X <- NULL


#--------- Independent sample with unequal probabilities

# non-informative weights
pop$prob0 <- rbeta(N, 1.5, 10)

# White noise variable
e <- rnorm(N, 0, 1)

# informative weights: only depends on the failure time
pop$prob1 <- exp(-1+pop$T)/(1+exp(-1+pop$T))
#pop$prob2 <- exp(-1+pop$t)/(1+exp(-1+pop$t))

# informative weights: depends on both W_1 and W_2
pop$prob2 <- exp(-2-2*pop$X1-pop$X2)/(1+exp(-2-2*pop$X1-pop$X2))

# informative weights: depends on both W_1 and W_2, and failure time 
pop$prob3 <- exp(-2-2*pop$X1-pop$X2+pop$T)/(1+exp(-2-2*pop$X1-pop$X2+pop$T))
#pop$prob5 <- exp(-2-2*pop$X1-pop$X2+pop$t)/(1+exp(-2-2*pop$X1-pop$X2+pop$t))


#--------- Complex sample design

# 1/2-stage clustered design (PSU size: 100) without stratification
#  with equal weights
nstrata <- 2
npsu <- 10000
nssu <- 100

pop <- pop[sample(1:N, N, replace=F), ]
pop$strata0 <- rep(1, N)
pop$strata1 <- rep(1:nstrata, rep(N/nstrata, nstrata))

xx <- 3*pop$T + e
pop <- pop[order(pop$strata1, xx), ]
pop$psu <- rep(1:npsu, rep(N/npsu, npsu))

# with unequal weights
tmp <- as.data.frame(sapply(pop[, c("T", "X1", "X2")], function(x)ave(x, pop$strata1, pop$psu, FUN = mean)))
pop$prob4 <- exp(-1-5*tmp$X1-2*tmp$X2)/(1+exp(-1-5*tmp$X1-2*tmp$X2))
#pop$prob5 <- exp(-4-2*tmp$X1-tmp$X2+0.5*tmp$T)/(1+exp(-4-2*tmp$X1-tmp$X2+0.5*tmp$T))


# stratified cluster sampling design

smp10 <- pop[cpxsmp(strata=pop$strata1, PSU=pop$psu, n_psu=50, ssu_size=10, PSU.prob=pop$prob6), ]

by(pop$T, pop$strata1, mean)
mean(pop$T)
mean(smp10$T)

smp10$wght=1/smp10$prob4
svyd <- svydesign(ids=~psu, strata=~strata1, weights=~wght, data=smp10, nested=T)
svymean(~T, design=svyd, deff=T)
confint(svymean(~T, design=svyd))

idx_smp <- as.data.frame(matrix(NA, n, nSim))
idx_smp0 <- as.data.frame(matrix(NA, n, nSim))
idx_smp1 <- as.data.frame(matrix(NA, n, nSim))
idx_smp2 <- as.data.frame(matrix(NA, n, nSim))
idx_smp3 <- as.data.frame(matrix(NA, n, nSim))
idx_smp4 <- as.data.frame(matrix(NA, n, nSim))
idx_smp5 <- as.data.frame(matrix(NA, n, nSim))

ID <- 1:N
for(i in 1:nSim){
	idx_smp[, i] <- sample(1:N, n, replace=F)
	idx_smp0[, i] <- ppss(pop$prob0, n)
	idx_smp1[, i] <- ppss(pop$prob1, n)
	idx_smp2[, i] <- ppss(pop$prob2, n)
	idx_smp3[, i] <- ppss(pop$prob3, n)
	idx_smp4[, i] <- cpxsmp(strata=pop$strata0, PSU=pop$psu, n_psu=10, ssu_size=100, PSU.prob=pop$prob4)
	idx_smp5[, i] <- cpxsmp(strata=pop$strata1, PSU=pop$psu, n_psu=20, ssu_size=25, PSU.prob=pop$prob4)
}

write.csv(idx_smp, "idx_smp.csv", col.names=T, sep=",")
write.csv(idx_smp0, "idx_smp0.csv", col.names=T, sep=",")
write.csv(idx_smp1, "idx_smp1.csv", col.names=T, sep=",")
write.csv(idx_smp2, "idx_smp2.csv", col.names=T, sep=",")
write.csv(idx_smp3, "idx_smp3.csv", col.names=T, sep=",")
write.csv(idx_smp4, "idx_smp4.csv", col.names=T, sep=",")
write.csv(idx_smp5, "idx_smp5.csv", col.names=T, sep=",")



#idx_smp0 <- read.csv("idx_smp0.csv", header=T, sep=",")
#idx_smp0$X <- NULL
#idx_smp1 <- read.csv("idx_smp1.csv", header=T, sep=",")
#idx_smp1$X <- NULL
#idx_smp2 <- read.csv("idx_smp2.csv", header=T, sep=",")
#idx_smp2$X <- NULL

#idx_smp0 <- idx_smp0[, idx]
#idx_smp1 <- idx_smp1[, idx]
#idx_smp2 <- idx_smp2[, idx]

sim_res1 <- as.data.frame(matrix(NA, nSim, 15*9))
sim_res2 <- as.data.frame(matrix(NA, nSim, 15*9))

for(i in 1:nSim){
	#draw simple random sample
	smp <- pop[idx_smp[, i], ]
	smp0 <- pop[idx_smp0[, i], ]
	smp1 <- pop[idx_smp1[, i], ]
	smp2 <- pop[idx_smp2[, i], ]
	smp3 <- pop[idx_smp3[, i], ]
	smp4 <- pop[idx_smp4[, i], ]
	smp5 <- pop[idx_smp5[, i], ]

	smp <- smp[order(smp$t), ]
	smp0 <- smp0[order(smp0$t), ]
	smp1 <- smp1[order(smp1$t), ]
	smp2 <- smp2[order(smp2$t), ]
	smp3 <- smp3[order(smp3$t), ]
	smp4 <- smp4[order(smp4$t), ]
	smp5 <- smp5[order(smp5$t), ]	

	#Calculate sampling weights
	smp0$wght <- 1/smp$prob0
	smp1$wght <- 1/smp1$prob1
	smp2$wght <- 1/smp2$prob2
	smp3$wght <- 1/smp3$prob3
	smp4$wght <- 1/smp4$prob4
	smp5$wght <- 1/smp5$prob4

	#Fit the Aalen additive models
	aalen <- aam1(smp$t, smp$delta, smp[, c("X1", "X2")], method="no")	
	aalen0 <- aam1(smp0$t, smp0$delta, smp0[, c("X1", "X2")], method="no")
	aalen0_adj <- aam1(smp0$t, smp0$delta, smp0[, c("X1", "X2")],  weight=smp0$wght, method="no")
	aalen1 <- aam1(smp1$t, smp1$delta, smp1[, c("X1", "X2")], method="no")
	aalen1_adj <- aam1(smp1$t, smp1$delta, smp1[, c("X1", "X2")], weight=smp1$wght, method="no")
	aalen2 <- aam1(smp2$t, smp2$delta, smp2[, c("X1", "X2")], weight=smp2$wght, method="no")
	aalen2_adj <- aam1(smp2$t, smp2$delta, smp2[, c("X1", "X2")], method="no")
	aalen3 <- aam1(smp3$t, smp3$delta, smp3[, c("X1", "X2")], method="no")
	aalen3_adj <- aam1(smp3$t, smp3$delta, smp3[, c("X1", "X2")], weight=smp3$wght, method="no")
	aalen4 <- aam1(smp4$t, smp4$delta, smp4[, c("X1", "X2")], method="no")
	aalen4_TSL <- aam1(smp4$t, smp4$delta, smp4[, c("X1", "X2")], weight=smp4$wght, strata=smp4$strata0, cluster=smp4$psu, method="TSL")
	aalen4_JK1 <- aam1(smp4$t, smp4$delta, smp4[, c("X1", "X2")], weight=smp4$wght, strata=smp4$strata0, cluster=smp4$psu, method="JK1")
	aalen5 <- aam1(smp5$t, smp5$delta, smp5[, c("X1", "X2")], method="no")
	aalen5_TSL <- aam1(smp5$t, smp5$delta, smp5[, c("X1", "X2")], weight=smp5$wght, strata=smp5$strata1, cluster=smp5$psu, method="TSL")
	aalen5_JK1 <- aam1(smp5$t, smp5$delta, smp5[, c("X1", "X2")], weight=smp5$wght, strata=smp5$strata1, cluster=smp5$psu, method="JK1")
	
	sim_res1[i, 1:9] <- smp.eval(aalen, smp[smp$delta==1, c("b00", "b01", "b02")])
	sim_res2[i, 1:9] <- smp.eval(aalen, smp[smp$delta==1, c("b10", "b11", "b12")])
	sim_res1[i, 10:18] <- smp.eval(aalen0, smp0[smp0$delta==1, c("b00", "b01", "b02")])
	sim_res2[i, 10:18] <- smp.eval(aalen0, smp0[smp0$delta==1, c("b10", "b11", "b12")])
	sim_res1[i, 19:27] <- smp.eval(aalen0_adj, smp0[smp0$delta==1, c("b00", "b01", "b02")])
	sim_res2[i, 19:27] <- smp.eval(aalen0_adj, smp0[smp0$delta==1, c("b10", "b11", "b12")])
	sim_res1[i, 28:36] <- smp.eval(aalen1, smp1[smp1$delta==1, c("b00", "b01", "b02")])
	sim_res2[i, 28:36] <- smp.eval(aalen1, smp1[smp1$delta==1, c("b10", "b11", "b12")])
	sim_res1[i, 37:45] <- smp.eval(aalen1_adj, smp1[smp1$delta==1, c("b00", "b01", "b02")])
	sim_res2[i, 37:45] <- smp.eval(aalen1_adj, smp1[smp1$delta==1, c("b10", "b11", "b12")])
	sim_res1[i, 46:54] <- smp.eval(aalen2, smp2[smp2$delta==1, c("b00", "b01", "b02")])
	sim_res2[i, 46:54] <- smp.eval(aalen2, smp2[smp2$delta==1, c("b10", "b11", "b12")])
	sim_res1[i, 55:63] <- smp.eval(aalen2_adj, smp2[smp2$delta==1, c("b00", "b01", "b02")])
	sim_res2[i, 55:63] <- smp.eval(aalen2_adj, smp2[smp2$delta==1, c("b10", "b11", "b12")])
	sim_res1[i, 64:72] <- smp.eval(aalen3, smp3[smp3$delta==1, c("b00", "b01", "b02")])
	sim_res2[i, 64:72] <- smp.eval(aalen3, smp3[smp3$delta==1, c("b10", "b11", "b12")])
	sim_res1[i, 73:81] <- smp.eval(aalen3_adj, smp3[smp3$delta==1, c("b00", "b01", "b02")])
	sim_res2[i, 73:81] <- smp.eval(aalen3_adj, smp3[smp3$delta==1, c("b10", "b11", "b12")])
	sim_res1[i, 82:90] <- smp.eval(aalen4, smp4[smp4$delta==1, c("b00", "b01", "b02")])
	sim_res2[i, 82:90] <- smp.eval(aalen4, smp4[smp4$delta==1, c("b10", "b11", "b12")])
	sim_res1[i, 91:99] <- smp.eval(aalen4_TSL, smp4[smp4$delta==1, c("b00", "b01", "b02")])
	sim_res2[i, 91:99] <- smp.eval(aalen4_TSL, smp4[smp4$delta==1, c("b10", "b11", "b12")])
	sim_res1[i, 100:108] <- smp.eval(aalen4_JK1, smp4[smp4$delta==1, c("b00", "b01", "b02")])
	sim_res2[i, 100:108] <- smp.eval(aalen4_JK1, smp4[smp4$delta==1, c("b10", "b11", "b12")])
	sim_res1[i, 109:117] <- smp.eval(aalen5, smp5[smp5$delta==1, c("b00", "b01", "b02")])
	sim_res2[i, 109:117] <- smp.eval(aalen5, smp5[smp5$delta==1, c("b10", "b11", "b12")])
	sim_res1[i, 118:126] <- smp.eval(aalen5_TSL, smp5[smp5$delta==1, c("b00", "b01", "b02")])
	sim_res2[i, 118:126] <- smp.eval(aalen5_TSL, smp5[smp5$delta==1, c("b10", "b11", "b12")])
	sim_res1[i, 127:135] <- smp.eval(aalen5_JK1, smp5[smp5$delta==1, c("b00", "b01", "b02")])
	sim_res2[i, 127:135] <- smp.eval(aalen5_JK1, smp5[smp5$delta==1, c("b10", "b11", "b12")])

	write.csv(sim_res1, "sim_res1.csv", col.names=T, sep=",")
	write.csv(sim_res2, "sim_res2.csv", col.names=T, sep=",")
}

tb1 <- matrix(apply(sim_res1, 2, function(x)mean(x, na.rm=T)), , 9, byrow=T)
#td1 <- data.frame(Measure=rep(c("Bias", "rMSE", "Cov"), 45), tb1)
tb2 <- matrix(apply(sim_res2, 2, function(x)mean(x, na.rm=T)), , 9, byrow=T)
#td2 <- data.frame(Measure=rep(c("Bias", "rMSE", "Cov"), 45), tb2)

	write.csv(tb1, "tb1_final1.csv", col.names=T, sep=",")
	write.csv(tb2, "tb2_final1.csv", col.names=T, sep=",")

