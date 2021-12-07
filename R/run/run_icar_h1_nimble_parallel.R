
# ------        Spanish breast cancer data - NIMBLE 1      ------ #
# -------   Spatio-temporal models: ICAR - RW1 - Type IV  ------- #
# -------                 H1 hyperpriors                  ------- #


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




#############################
##  ICAR model - Type IV   ##
#############################

## 0. Temporal structure matrix for a RW1 prior ##
Dm <- diff(diag(t),differences=1)
Rt <- t(Dm)%*%Dm

cov.Rt <- MASS::ginv(Rt)    # Generalized Inverse of Rt
chol.cov.Rt <- chol(cov.Rt)



## 1. Define data ##
data <- list(y = y)


## 2. Define constants ##

# Nimble 1
constants <- list(N = N, t = t, E = expected, L = length(adj), L.rw1 = length(adj.rw1),
                  adj = adj, num = num, weights = weights,
                  adj.rw1 = adj.rw1, num.rw1 = num.rw1, weights.rw1 = weights.rw1,
                  Achol = chol.cov.Rt, nt= N*t)


## 3. Define initial values ##

# Nimble 1
inits <- list(alpha0 = rnorm(1,0,0.1), sd.u = runif(1,0,1),
              sd.temp = runif(1,0,1), sd.t4 = runif(1,0,1),
              spat.u = rnorm(N), temp = rnorm(t),
              stType4 = matrix(rnorm(N*t), nrow=N, ncol=t))


## 4. Variables to retrive ##
param <- c('alpha0', 'sd.u', 'sd.temp', 'sd.t4',
           'tau.u', 'tau.temp', 'tau.t4', 
           'var.u', 'var.temp', 'var.t4',
           'spat.u', 'temp', 'stType4',
           'Etemp', 'Espat', 'Est', 'RR',
           'Dev', 'sumDev', 'mu')




## RUN IN PARALLEL ##

library(parallel)
parallelize <- makeCluster(3)   # Run 3 chains, each one in a different core


icar.type4.parallel <- function(seed, data, constants, inits, param) {
  library(nimble)
  source("./icar_models/icar_h1_nimble1.R")
  code <- icar.type4.h1
  num.iter <- 30000
  num.burnin <- 5000
  num.thin <- 75

  model <- nimbleModel(code, constants = constants, data = data, 
                            inits = inits)
  cmodel <- compileNimble(model)
  conf <- configureMCMC(model, monitors = param)
  MCMC <- buildMCMC(conf, enableWAIC = TRUE)
  cMCMC <- compileNimble(MCMC, project=cmodel)
  
  result.nimble <- runMCMC(cMCMC, niter = num.iter, nburnin = num.burnin, 
                           thin = num.thin, samplesAsCodaMCMC = TRUE, 
                           summary=TRUE, WAIC = TRUE)
  
  return(result.nimble)
}


result.chains <- parLapply(cl = parallelize, X = 1:3, 
                          fun = icar.type4.parallel, 
                          # seed = set.seed(123),
                          data = data,
                          constants=constants,
                          inits=inits,
                          param=param)

stopCluster(parallelize)

save(result.chains, file = "icar_type4_model_h1_nimble1.Rdata")



