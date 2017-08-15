
## ALL minrms scores
library(data.table)
W <-fread("Documents/protein_interaction/zz_all_minrms_scores/zz_distances_for_query_10000")
W <- as.matrix(W)
hist(W)
plot(ecdf(W))

idx <- sample(5822045, 1000)
W1 <- W[idx]
hist(W1)
plot(ecdf(W1))



#####################################
#####################################

## X is 5000 minrms score for 5000 randomely picked pair of triplets
X <- read.table("Documents/protein_interaction/z_random_pairs_rms_10MIL_vals")
X <- as.matrix(X)
dim(X)
hist(X, breaks=50)
hist(X, breaks=100, xlim = c(0,0.8))
plot(ecdf(X))
plot(ecdf(X), xlim=c(0.01,2), ylim=c(0.01,1))
which(X == min(X))




## Y is MIN minrms score of 1000 randomely picked triplet
Y <- read.table("Documents/protein_interaction/zzz_values")
Y <- as.matrix(Y)
dim(Y)
hist(Y, breaks = 50)
hist(Y, breaks = 1000, xlim = c(0,0.05))
plot(ecdf(Y))
plot(ecdf(Y), xlim = c(0.045,0.1))
grid(nx = NULL, ny = NULL, col = "darkgray", lty = "dotted")



#####################################
###### Fitting a Distribution #######
#####################################
# Distributions available: "beta", "cauchy", "chi-squared", "exponential", "f", "gamma", "geometric", 
# "log-normal", "lognormal", "logistic", "negative binomial", "normal", "Poisson", "t" and "weibull". 

library(fitdistrplus)
Z <- as.vector(Y)

## use 'method = c("mle", "mme", "qme", "mge")' to specify method by which to estimate parameters
fit.norm <- fitdist(Z, "norm")
plot(fit.norm)
summary(fit.norm)

fit.weibull <- fitdist(Z, "weibull")
plot(fit.weibull)
summary(fit.weibull)
fit.weibull$estimate




#################################################
#### Generalized Extereme Value Distribution ####
#################################################


library(extRemes)
library(RCurl)


Z <- as.vector(Y)
# parameter 'method' can be used to specify the method by which parameters are estimated, maximum likelihood by deafult

# Gumbel
fit0 <- fevd(Z, type="Gumbel")
fit0
plot(fit0)
plot(fit0, "trace")
return.level(fit0)

# Generalized Extreme Value distribution
fit1 <- fevd(Z, type="GEV")
fit1
plot(fit1)
plot(fit1, "trace") 
return.level(fit1)







########################################
### calculating Weibull parameters #####
########################################


D <- 5822045

S <- sort(X)

Fr1 <- 0.2
r1 <- S[2000]
Fr2 <- 0.6
r2 <- S[6000]

alpha <- ( log(-D*log(Fr1)) - log(-D*log(Fr2)) ) / (log(r1) - log(r2))
sigma <- exp( (alpha*log(r1) - log(-D*log(Fr1))) / (alpha) )

alpha
sigma




###############################################################
###############################################################
###############################################################
###############################################################


library(data.table)
X <-fread("Documents/protein_interaction/z_random_pairs_rms_500MIL_values")
X <- as.matrix(X)
dim(X)
hist(X, breaks = 1000, xlab = "Min RMS Score", main = "")
hist(X, breaks = 1000, xlim = c(0,5), xlab = "Min RMS Score", main = "")
hist(X, breaks = 2000, xlim = c(0,1), xlab = "Min RMS Score", main = "")
hist(X, breaks = 2000, xlim = c(0,0.5), xlab = "Min RMS Score", main = "")

Y <- X[X<.8]
Y <- as.matrix(Y)
dim(Y)
hist(Y, breaks = 200, xlab = "Min RMS Score", main = "")


## CDF
plot(ecdf(X))


D <- 5822045

alpha <- ( log(-D*log(1-Fr1)) - log(-D*log(1-Fr2)) ) / (log(r1) - log(r2))
sigma <- exp( (alpha*log(r1) - log(-D*log(1-Fr1))) / (alpha) )
alpha
sigma


##############################
##############################
##############################

Z <- fread("Documents/protein_interaction/z_random_pairs_rms_500MIL_values")
Z <- as.matrix(Z)


idx <- sample(5822045, 100000)
X <- Z[idx]
X <- as.matrix(X)


######### best hits excluding threshold: 9 #########
Y <- fread("Documents/protein_interaction/sample_and_all_best_rms/10k/z_best_rms_10K_sample_copy_values_threshold_9")
Y <- as.matrix(Y)
hist(Y, breaks = 100)
plot(ecdf(Y))
###########################

D<- 1000
  
cdf <- ecdf(X)
plot(X,(1-cdf(X))^D,ylim=c(0,1),type="l");
par(new=T);
plot(X,1-cdf(X),ylim=c(0,1),type="l");
par(new=T);
plot(X,dweibull(r,shape=alpha,scale=sigma),type="l",col=2);
par(new=T);
hist(Y, breaks = 50,xlim=c(0,2))



S <- X[order(X)]
S <- as.matrix(S)

Fr1 <- 0.3
r1 <- S[2000]
Fr2 <- 0.4
r2 <- S[6000]


D <- 5822045

alpha <- ( log(-D*log(1-Fr1)) - log(-D*log(1-Fr2)) ) / (log(r1) - log(r2))
sigma <- exp( (alpha*log(r1) - log(-D*log(1-Fr1))) / (alpha) )
alpha
sigma




                                                                                                                                                                       
                                                                                                                                                                                 




