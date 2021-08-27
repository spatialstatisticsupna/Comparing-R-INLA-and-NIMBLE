
# --------------     Spanish breast cancer data - NIMBLE 2     -------------- #
# -------   Spatio-temporal models: ICAR - RW1 - Type I, II, III, IV  ------- #
# --------------                H1 hyperpriors                 -------------- #



############################
##  ICAR - RW1 - Type I   ##
############################
icar.type1.h1 <- nimbleCode({
  # Priors
  alpha0 ~ dflat()
  sd.u ~ dunif(0, 100)
  tau.u <- 1/sd.u^2
  var.u <- sd.u^2
  sd.temp ~ dunif(0, 100)
  tau.temp <- 1/sd.temp^2
  var.temp <- sd.temp^2
  sd.t1 ~ dunif(0, 100)
  tau.t1 <- 1/sd.t1^2
  var.t1 <- sd.t1^2
  
  # Spatial random effects
  spat.u[1:N] ~ dcar_normal(adj[1:L], weights[1:L], num[1:N], 1, zero_mean = 1)
  
  # Temporal random effects
  temp[1]  <- temp.zero
  for (k in 2:t){
    temp[k] ~ dnorm(temp[k-1], 1)
  }
  
  for(i in 1:N){ 
    for (k in 1:t){
      # Spatio-temporal random effects (type I)
      stType1[i,k] ~ dnorm(0, sd=sd.t1)
      y[i,k] ~ dpois(lambda[i,k])
      log(lambda[i,k]) <- log(E[i,k]) + alpha0 + sd.u*spat.u[i] + sd.temp*temp[k] + stType1.centered[i,k]
      RR[i,k] <- exp(alpha0 + sd.u*spat.u[i] + sd.temp*temp[k] + stType1.centered[i,k])
      
      # Spatio-temporal trends
      Est[i,k] <- exp(stType1.centered[i,k])

      # Deviance
      Dev[i,k] <- -2*(-lambda[i,k] + y[i,k]*log(lambda[i,k])-lfactorial(y[i,k]))
    }
    
    # Spatial pattern
    Espat[i] <- exp(sd.u*spat.u[i])
  }
  
  stType1.mean <- sum(stType1[1:N, 1:t])/(N*t)
  
  for(i in 1:N){ 
    for (k in 1:t){
      stType1.centered[i,k] <- stType1[i, k] - stType1.mean
    }}
  
  # Temporal trend
  for(k in 1:t){
    Etemp[k] <- exp(sd.temp*temp[k])
  }
  
  # Posterior mean deviance
  sumDev <- sum(Dev[1:N, 1:t])
})




#############################
##  ICAR - RW1 - Type II   ##
#############################
icar.type2.h1 <- nimbleCode({
  # Priors
  alpha0 ~ dflat()
  sd.u ~ dunif(0, 100)
  tau.u <- 1/sd.u^2
  var.u <- sd.u^2
  sd.temp ~ dunif(0, 100)
  tau.temp <- 1/sd.temp^2
  var.temp <- sd.temp^2
  sd.t2 ~ dunif(0, 100)
  tau.t2 <- 1/sd.t2^2
  var.t2 <- sd.t2^2
  
  # Spatial random effects
  spat.u[1:N] ~ dcar_normal(adj[1:L], weights[1:L], num[1:N], 1, zero_mean = 1)
  
  # Temporal random effects
  temp[1]  <- temp.zero
  for (k in 2:t){
    temp[k] ~ dnorm(temp[k-1], 1)
  }
  
  # Spatio-temporal random effects (type II)
  for (i in 1:N){
    stType2[i, 1:t] ~ dcar_normal(adj.t2[1:L.t2], weights.t2[1:L.t2], num.t2[1:t], 1, zero_mean = 1) 
  }
  
  for(i in 1:N){ 
    for (k in 1:t){
      y[i,k] ~ dpois(lambda[i,k])
      log(lambda[i,k]) <- log(E[i,k]) + alpha0 + sd.u*spat.u[i] + sd.temp*temp[k] + sd.t2*stType2[i,k]
      RR[i,k] <- exp(alpha0 + sd.u*spat.u[i] + sd.temp*temp[k] + sd.t2*stType2[i,k])
      
      # Spatio-temporal trends
      Est[i,k] <- exp(sd.t2*stType2[i,k])

      # Deviance
      Dev[i,k] <- -2*(-lambda[i,k] + y[i,k]*log(lambda[i,k])-lfactorial(y[i,k]))
    }
    # Spatial pattern
    Espat[i] <- exp(sd.u*spat.u[i])
  }
  
  # Temporal trend
  for(k in 1:t){
    Etemp[k] <- exp(sd.temp*temp[k])
  }
  
  # Posterior mean deviance
  sumDev <- sum(Dev[1:N, 1:t])
})





#############################
##  ICAR - RW1 - Type III  ##
#############################
icar.type3.h1 <- nimbleCode({
  # Priors
  alpha0 ~ dflat()
  sd.u ~ dunif(0, 100)
  tau.u <- 1/sd.u^2
  var.u <- sd.u^2
  sd.temp ~ dunif(0, 100)
  tau.temp <- 1/sd.temp^2
  var.temp <- sd.temp^2
  sd.t3 ~ dunif(0, 100)
  tau.t3 <- 1/sd.t3^2
  var.t3 <- sd.t3^2
  
  # Spatial random effects
  spat.u[1:N] ~ dcar_normal(adj[1:L], weights[1:L], num[1:N], 1, zero_mean = 1)
  
  # Temporal random effects
  temp[1]  <- temp.zero
  for (k in 2:t){
    temp[k] ~ dnorm(temp[k-1], 1)
  }
  
  # Spatio-temporal random effects (type III)
  for (k in 1:t){
    stType3[1:N, k] ~ dcar_normal(adj.t3[1:L.t3], weights.t3[1:L.t3], num.t3[1:N], 1, zero_mean = 1) 
  }
  
  for(i in 1:N){ 
    for (k in 1:t){
      y[i,k] ~ dpois(lambda[i,k])
      log(lambda[i,k]) <- log(E[i,k]) + alpha0 + sd.u*spat.u[i] + sd.temp*temp[k] + sd.t3*stType3[i,k]
      RR[i,k] <- exp(alpha0 + sd.u*spat.u[i] + sd.temp*temp[k] + sd.t3*stType3[i,k])
      
      # Spatio-temporal trends
      Est[i,k] <- exp(sd.t3*stType3[i,k])
      
      # Deviance
      Dev[i,k] <- -2*(-lambda[i,k] + y[i,k]*log(lambda[i,k])-lfactorial(y[i,k]))
    }
    # Spatial pattern
    Espat[i] <- exp(sd.u*spat.u[i])
  }
  
  # Temporal trend
  for(k in 1:t){
    Etemp[k] <- exp(sd.temp*temp[k])
  }
  
  # Posterior mean deviance
  sumDev <- sum(Dev[1:N, 1:t])
})




#############################
##  ICAR - RW1 - Type IV   ##
#############################
icar.type4.h1 <- nimbleCode({
  # Priors
  alpha0 ~ dflat()
  sd.u ~ dunif(0, 100)
  tau.u <- 1/sd.u^2
  var.u <- sd.u^2
  sd.temp ~ dunif(0, 100)
  tau.temp <- 1/sd.temp^2
  var.temp <- sd.temp^2
  sd.t4 ~ dunif(0, 100)
  tau.t4 <- 1/sd.t4^2
  var.t4 <- sd.t4^2
  
  # Spatial random effects
  spat.u[1:N] ~ dcar_normal(adj[1:L], weights[1:L], num[1:N], 1, zero_mean = 1)
  
  # Temporal random effects
  temp[1]  <- temp.zero
  for (k in 2:t){
    temp[k] ~ dnorm(temp[k-1], 1)
  }
  
  # Spatio-temporal random effects (type IV)
  for (k in 1:t){
    u.t4[1:N, k] ~ dcar_normal(adj[1:L], weights[1:L], num[1:N], 1, zero_mean = 1)
  }
  
  for (i in 1:N){
    stType4[i, 1:t] <- u.t4[i, 1:t] %*% Achol[1:t, 1:t]
  }
  
  for(i in 1:N) {
    for (k in 1:t){ 
      y[i,k] ~ dpois(lambda[i,k])
      log(lambda[i,k]) <- log(E[i,k]) + alpha0 + sd.u*spat.u[i] + sd.temp*temp[k] + sd.t4*stType4[i,k]
      RR[i,k] <- exp(alpha0 + sd.u*spat.u[i] + sd.temp*temp[k] + sd.t4*stType4[i,k])
      
      # Spatio-temporal trend
      Est[i,k] <- exp(sd.t4*stType4[i,k])
      
      # Deviance
      Dev[i,k] <- -2*(-lambda[i,k] + y[i,k]*log(lambda[i,k])-lfactorial(y[i,k]))
    }
    # Spatial pattern
    Espat[i] <- exp(sd.u*spat.u[i])
  }
  
  # Temporal trend
  for(k in 1:t){
    Etemp[k] <- exp(sd.temp*temp[k])
  }
  
  # Posterior mean deviance
  sumDev <- sum(Dev[1:N, 1:t])
})




