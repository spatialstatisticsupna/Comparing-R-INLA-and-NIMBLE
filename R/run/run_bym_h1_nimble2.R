
# ------              Spanish breast cancer data - NIMBLE 2            ------ #
# -------    Spatio-temporal models: BYM - RW1 - Type I, II, III, IV  ------- #
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

source("bym_models/bym_h1_nimble2.R")




##################################
##  n.chains, n.iter, n.burnin  ##
##################################

num.chains <- 3
num.iter <- 30000
num.burnin <- 5000
num.thin <- 75




###########################
##  BYM model - Type I   ##
###########################
code <- bym.type1.h1


## 1. Define data ##
data <- list(y = y)


## 2. Define constants ##
constants <- list(N = N, t = t, E = expected, 
                  adj = adj, num = num, weights = weights, L = length(adj),
                  temp.zero = 0)


## 3. Define initial values ##
inits <- list(alpha0 = rnorm(1,0,0.1), sd.u = runif(1,0,1), sd.v = runif(1,0,1),
              sd.temp = runif(1,0,1), sd.t1 = runif(1,0,1),
              spat.u = rnorm(50), temp = c(0, rnorm(20)),
              stType1 = matrix(rnorm(N*t), nrow=N, ncol=t))



## 4. Run the model ##

model <- nimbleModel(code, constants = constants, data = data, inits = inits)
cmodel <- compileNimble(model)

conf <- configureMCMC(model, monitors = c('alpha0', 'sd.u', 'sd.v', 'sd.temp', 'sd.t1',
                                          'tau.u', 'tau.v', 'tau.temp', 'tau.t1', 
                                          'var.u', 'var.v', 'var.temp', 'var.t1',
                                          'spat.u', 'spat.v', 'temp', 'stType1',
                                          'Eu', 'Ev', 'Etemp', 'Espat', 'Est', 'RR',
                                          'Dev', 'sumDev', 'lambda'))

MCMC <- buildMCMC(conf, enableWAIC = TRUE)
cMCMC <- compileNimble(MCMC, project=cmodel, resetFunctions = TRUE)


## 5. Obtain the samples ##
result.nimble <- runMCMC(cMCMC, niter = num.iter, nburnin = num.burnin, 
                         nchains = num.chains, thin = num.thin, 
                         samplesAsCodaMCMC = TRUE, summary=TRUE, WAIC = TRUE)

## 6. Save results ##
save(result.nimble, file = "bym_typeI_model_h1_nimble2.Rdata")





############################
##  BYM model - Type II   ##
############################
code <- bym.type2.h1



## 1. Define data ##
data <- list(y = y)


## 2. Define constants ##
constants <- list(N = N, t = t, E = expected, 
                  adj = adj, num = num, weights = weights, L = length(adj),
                  adj.t2 = adj.rw1, num.t2 = num.rw1, weights.t2 = weights.rw1,
                  L.t2 = length(adj.rw1), temp.zero = 0)


## 3. Define initial values ##
inits <- list(alpha0 = rnorm(1,0,0.1), sd.u = runif(1,0,1), sd.v = runif(1,0,1),
              sd.temp = runif(1,0,1), sd.t2 = runif(1,0,1),
              spat.u = rnorm(50), temp = c(0, rnorm(20)),
              stType2 = matrix(rnorm(N*t), nrow=N, ncol=t))



## 4. Run the model ##
model <- nimbleModel(code, constants = constants, data = data, inits = inits)
cmodel <- compileNimble(model)

conf <- configureMCMC(model, monitors = c('alpha0', 'sd.u', 'sd.v', 'sd.temp', 'sd.t2',
                                          'tau.u', 'tau.v', 'tau.temp', 'tau.t2', 
                                          'var.u', 'var.v', 'var.temp', 'var.t2',
                                          'spat.u', 'spat.v', 'temp', 'stType2',
                                          'Eu', 'Ev', 'Etemp', 'Espat', 'Est', 'RR',
                                          'Dev', 'sumDev', 'lambda'))

MCMC <- buildMCMC(conf, enableWAIC = TRUE)
cMCMC <- compileNimble(MCMC, project=cmodel, resetFunctions = TRUE)


## 5. Obtain the samples ##
result.nimble <- runMCMC(cMCMC, niter = num.iter, nburnin = num.burnin, 
                         nchains = num.chains, thin = num.thin, 
                         samplesAsCodaMCMC = TRUE, summary=TRUE, WAIC = TRUE)

### 6. Save results ###
save(result.nimble, file = "bym_typeII_model_h1_nimble2.Rdata")





#############################
##  BYM model - Type III   ##
#############################
code <- bym.type3.h1



## 1. Define data ##
data <- list(y = y)


## 2. Define constants ##
constants <- list(N = N, t = t, E = expected, 
                  adj = adj, num = num, weights = weights, L = length(adj), 
                  adj.t3 = adj, num.t3 = num, weights.t3 = weights,
                  L.t3 = length(adj), temp.zero = 0)



## 3. Define initial values ##
inits <- list(alpha0 = rnorm(1,0,0.1), sd.u = runif(1,0,1), sd.v = runif(1,0,1),
              sd.temp = runif(1,0,1), sd.t3 = runif(1,0,1),
              spat.u = rnorm(50), temp = c(0, rnorm(20)),
              stType3 = matrix(rnorm(N*t), nrow=N, ncol=t))



## 4. Run the model ##

model <- nimbleModel(code, constants = constants, data = data, inits = inits)
cmodel <- compileNimble(model)

conf <- configureMCMC(model, monitors = c('alpha0', 'sd.u', 'sd.v', 'sd.temp', 'sd.t3',
                                          'tau.u', 'tau.v', 'tau.temp', 'tau.t3', 
                                          'var.u', 'var.v', 'var.temp', 'var.t3',
                                          'spat.u', 'spat.v', 'temp', 'stType3',
                                          'Eu', 'Ev', 'Etemp', 'Espat', 'Est', 'RR',
                                          'Dev', 'sumDev', 'lambda'))


MCMC <- buildMCMC(conf, enableWAIC = TRUE)
cMCMC <- compileNimble(MCMC, project=cmodel, resetFunctions = TRUE)



## 5. Obtain the samples ##
result.nimble <- runMCMC(cMCMC, niter = num.iter, nburnin = num.burnin, 
                         nchains = num.chains, thin = num.thin, 
                         samplesAsCodaMCMC = TRUE, summary=TRUE, WAIC = TRUE)
 

## 6. Save results ##
save(result.nimble, file = "bym_typeIII_model_h1_nimble2.Rdata")




############################
##  BYM model - Type IV   ##
############################
code <- bym.type4.h1


## 0. Temporal structure matrix for a RW1 prior ###
Dm <- diff(diag(t),differences=1)
Rt <- t(Dm)%*%Dm

cov.Rt <- MASS::ginv(Rt)    # Generalized Inverse of Rt
chol.cov.Rt <- chol(cov.Rt)


## 1. Define data ##
data <- list(y = y)


## 2. Define constants ##
constants <- list(N = N, t = t, E = expected, L = length(adj),
                  adj = adj, num = num, weights = weights,
                  Achol = chol.cov.Rt, nt= N*t, temp.zero = 0)


## 3. Define initial values ##
inits <- list(alpha0 = rnorm(1,0,0.1), sd.u = runif(1,0,1), sd.v = runif(1,0,1),
              sd.temp = runif(1,0,1), sd.t4 = runif(1,0,1),
              spat.u = rnorm(50), temp = c(0, rnorm(20)),
              stType4 = matrix(rnorm(N*t), nrow=N, ncol=t))



## 4. Run the model ##
model <- nimbleModel(code, constants = constants, data = data, inits = inits)

cmodel <- compileNimble(model)

conf <- configureMCMC(model, monitors = c('alpha0', 'sd.u', 'sd.v', 'sd.temp', 'sd.t4',
                                          'tau.u', 'tau.v', 'tau.temp', 'tau.t4', 
                                          'var.u', 'var.v', 'var.temp', 'var.t4',
                                          'spat.u', 'spat.v', 'temp', 'stType4',
                                          'Eu', 'Ev', 'Etemp', 'Espat', 'Est', 'RR',
                                          'Dev', 'sumDev', 'lambda'))

MCMC <- buildMCMC(conf, enableWAIC = TRUE)
cMCMC <- compileNimble(MCMC, project=cmodel, resetFunctions = TRUE)

## 5. Obtain the samples ##
result.nimble <- runMCMC(cMCMC, niter = num.iter, nburnin = num.burnin, 
                         nchains = num.chains, thin = num.thin, 
                         samplesAsCodaMCMC = TRUE, summary=TRUE, WAIC = TRUE)

## 6. Save results ##
save(result.nimble, file = "bym_typeIV_model_h1_nimble2.Rdata")






