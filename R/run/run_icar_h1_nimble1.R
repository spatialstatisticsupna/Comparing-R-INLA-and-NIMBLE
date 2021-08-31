
# ------             Spanish breast cancer data - NIMBLE 1             ------ #
# -------   Spatio-temporal models: ICAR - RW1 - Type I, II, III, IV  ------- #
# --------------                H1 hyperpriors                 -------------- #


# Load NIMBLE package
library(nimble)

# Set working directory
setwd("")



#################################
##  Load and prepare the data  ##
#################################

# Load the data
load("BreastCancer_data.Rdata")
head(Data)

N <- length(unique(Data$Area))
t <- length(unique(Data$Year))


# Spatial random effects: specify adjacency and weights vectors
num <- diag(Rs)
adj <- c()

for (i in 1:dim(Rs)[1]){
  neighbours <- unname(which(Rs[i, ]==-1))
  adj <- c(adj, neighbours)
}
weights <- rep(1, length(adj))


# Temporal random effects: specify adjacency and weights vectors
weights.rw1 <- c()
adj.rw1 <- c()
num.rw1 <- c()

for(i in 1:1) {
  weights.rw1[i] <- 1; adj.rw1[i] <- i+1; num.rw1[i] <- 1
}
for(i in 2:(t-1)) {
  weights.rw1[2+(i-2)*2] <- 1; adj.rw1[2+(i-2)*2] <- i-1
  weights.rw1[3+(i-2)*2] <- 1; adj.rw1[3+(i-2)*2] <- i+1; num.rw1[i] <- 2
}
for(i in t:t) {
  weights.rw1[(i-2)*2 + 2] <- 1; adj.rw1[(i-2)*2 + 2] <- i-1; num.rw1[i] <- 1
}


# Counts and expected cases in matrix form
expected <-  matrix(Data$Expected, nrow = 50)
y <- matrix(Data$Counts, nrow = 50)




#######################
##  Load the models  ##
#######################

source("icar_models/icar_h1_nimble1.R")





##################################
##  n.chains, n.iter, n.burnin  ##
##################################

num.chains <- 3
num.iter <- 30000
num.burnin <- 5000
num.thin <- 75




############################
##  ICAR model - Type I   ##
############################
code <- icar.type1.h1


## 1. Define data ##
data <- list(y = y)


## 2. Define constants ##
constants <- list(N = N, t = t, E = expected, 
                  adj = adj, num = num, weights = weights, L = length(adj),
                  adj.rw1 = adj.rw1, num.rw1 = num.rw1, 
                  weights.rw1 = weights.rw1, L.rw1 = length(adj.rw1))



## 3. Define initial values ##
inits <- list(alpha0 = rnorm(1,0,0.1), sd.u = runif(1,0,1),
              sd.temp = runif(1,0,1), sd.t1 = runif(1,0,1),
              spat.u = rnorm(N), temp = rnorm(t),
              stType1 = matrix(rnorm(N*t), nrow=N, ncol=t))



## 4. Run the model ##
model <- nimbleModel(code, constants = constants, data = data, inits = inits)
cmodel <- compileNimble(model)

conf <- configureMCMC(model, monitors = c('alpha0', 'sd.u', 'sd.temp', 'sd.t1',
                                          'tau.u', 'tau.temp', 'tau.t1', 
                                          'var.u', 'var.temp', 'var.t1',
                                          'spat.u', 'temp', 'stType1',
                                          'Etemp', 'Espat', 'Est', 'RR',
                                          'Dev', 'sumDev', 'lambda'))

MCMC <- buildMCMC(conf, enableWAIC = TRUE)
cMCMC <- compileNimble(MCMC, project=cmodel, resetFunctions = TRUE)


## 5. Obtain the samples ##
result.nimble <- runMCMC(cMCMC, niter = num.iter, nburnin = num.burnin, 
                         nchains = num.chains, thin = num.thin, 
                         samplesAsCodaMCMC = TRUE, summary=TRUE, WAIC = TRUE)

## 6. Save results ##
save(result.nimble, file = "icar_typeI_model_h1_nimble1.Rdata")




#############################
##  ICAR model - Type II   ##
#############################
code <- icar.type2.h1


## 1. Define data ##
data <- list(y = y)


## 2. Define constants ##
constants <- list(N = N, t = t, E = expected, 
                  adj = adj, num = num, weights = weights, L = length(adj),
                  adj.rw1 = adj.rw1, num.rw1 = num.rw1, 
                  weights.rw1 = weights.rw1, L.rw1 = length(adj.rw1),
                  adj.t2 = adj.rw1, num.t2 = num.rw1, weights.t2 = weights.rw1,
                  L.t2 = length(adj.rw1))



## 3. Define initial values ##
inits <- list(alpha0 = rnorm(1,0,0.1), sd.u = runif(1,0,1),
              sd.temp = runif(1,0,1), sd.t2 = runif(1,0,1),
              spat.u = rnorm(N), temp = rnorm(t), 
              stType2 = matrix(rnorm(N*t), nrow=N, ncol=t))



## 4. Run the model ##
model <- nimbleModel(code, constants = constants, data = data, inits = inits)
cmodel <- compileNimble(model)

conf <- configureMCMC(model, monitors = c('alpha0', 'sd.u', 'sd.temp', 'sd.t2',
                                          'tau.u', 'tau.temp', 'tau.t2', 
                                          'var.u', 'var.temp', 'var.t2',
                                          'spat.u', 'temp', 'stType2',
                                          'Etemp', 'Espat', 'Est', 'RR',
                                          'Dev', 'sumDev', 'lambda'))

MCMC <- buildMCMC(conf, enableWAIC = TRUE)
cMCMC <- compileNimble(MCMC, project=cmodel, resetFunctions = TRUE)


## 5. Obtain the samples ##
result.nimble <- runMCMC(cMCMC, niter = num.iter, nburnin = num.burnin, 
                         nchains = num.chains, thin = num.thin, 
                         samplesAsCodaMCMC = TRUE, summary=TRUE, WAIC = TRUE)

## 6. Save results ##
save(result.nimble, file = "icar_typeII_model_h1_nimble1.Rdata")




##############################
##  ICAR model - Type III   ##
##############################
code <- icar.type3.h1


## 1. Define data ##
data <- list(y = y)


## 2. Define constants ##
constants <- list(N = N, t = t, E = expected, 
                  adj = adj, num = num, weights = weights, L = length(adj), 
                  adj.rw1 = adj.rw1, num.rw1 = num.rw1, weights.rw1 = weights.rw1, 
                  L.rw1 = length(adj.rw1),
                  adj.t3 = adj, num.t3 = num, weights.t3 = weights,
                  L.t3 = length(adj))


## 3. Define initial values ##
inits <- list(alpha0 = rnorm(1,0,0.1), sd.u = runif(1,0,1),
              sd.temp = runif(1,0,1), sd.t3 = runif(1,0,1),
              spat.u = rnorm(N), temp = rnorm(t), 
              stType3 = matrix(rnorm(N*t), nrow=N, ncol=t))



## 4. Run the model ##
model <- nimbleModel(code, constants = constants, data = data, inits = inits)
cmodel <- compileNimble(model)

conf <- configureMCMC(model, monitors = c('alpha0', 'sd.u', 'sd.temp', 'sd.t3',
                                          'tau.u', 'tau.temp', 'tau.t3', 
                                          'var.u', 'var.temp', 'var.t3',
                                          'spat.u', 'temp', 'stType3',
                                          'Etemp', 'Espat', 'Est', 'RR',
                                          'Dev', 'sumDev', 'lambda'))

MCMC <- buildMCMC(conf, enableWAIC = TRUE)
cMCMC <- compileNimble(MCMC, project=cmodel, resetFunctions = TRUE)


## 5. Obtain the samples ##
result.nimble <- runMCMC(cMCMC, niter = num.iter, nburnin = num.burnin, 
                         nchains = num.chains, thin = num.thin, 
                         samplesAsCodaMCMC = TRUE, summary=TRUE, WAIC = TRUE)
 
## 6. Save results ##
save(result.nimble, file = "icar_typeIII_model_h1_nimble1.Rdata")




#############################
##  ICAR model - Type IV   ##
#############################
code <- icar.type4.h1


## 0. Temporal structure matrix for a RW1 prior ##
Dm <- diff(diag(t),differences=1)
Rt <- t(Dm)%*%Dm

cov.Rt <- MASS::ginv(Rt)    # Generalized Inverse of Rt
chol.cov.Rt <- chol(cov.Rt)



## 1. Define data ##
data <- list(y = y)


## 2. Define constants ##
constants <- list(N = N, t = t, E = expected, L = length(adj), L.rw1 = length(adj.rw1),
                  adj = adj, num = num, weights = weights,
                  adj.rw1 = adj.rw1, num.rw1 = num.rw1, weights.rw1 = weights.rw1,
                  Achol = chol.cov.Rt, nt= N*t)


## 3. Define initial values ##
inits <- list(alpha0 = rnorm(1,0,0.1), sd.u = runif(1,0,1),
              sd.temp = runif(1,0,1), sd.t4 = runif(1,0,1),
              spat.u = rnorm(N), temp = rnorm(t),
              stType4 = matrix(rnorm(N*t), nrow=N, ncol=t))



## 4. Run the model ##
model <- nimbleModel(code, constants = constants, data = data, inits = inits)
cmodel <- compileNimble(model)

conf <- configureMCMC(model, monitors = c('alpha0', 'sd.u', 'sd.temp', 'sd.t4',
                                          'tau.u', 'tau.temp', 'tau.t4', 
                                          'var.u', 'var.temp', 'var.t4',
                                          'spat.u', 'temp', 'stType4',
                                          'Etemp', 'Espat', 'Est', 'Ealpha0.temp', 
                                          'RR', 'Dev', 'sumDev', 'lambda'))

MCMC <- buildMCMC(conf, enableWAIC = TRUE)
cMCMC <- compileNimble(MCMC, project=cmodel, resetFunctions = TRUE)


## 5. Obtain the samples ##
result.nimble <- runMCMC(cMCMC, niter = num.iter, nburnin = num.burnin, 
                         nchains = num.chains, thin = num.thin, 
                         samplesAsCodaMCMC = TRUE, summary=TRUE, WAIC = TRUE)

## 6. Save results ##
save(result.nimble, file = "icar_typeIV_model_h1_nimble1.Rdata")






