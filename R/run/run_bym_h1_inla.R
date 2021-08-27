
# --------------         Spanish breast cancer data - INLA     -------------- #
# -------   Spatio-temporal models: BYM - RW1 - Type I, II, III, IV   ------- #
# --------------                 H1 hyperpriors                -------------- #


# Load the R package
library(INLA)

# Set working directory
setwd("")



#################################
##  Load and prepare the data  ##
#################################

load("BreastCancer_data.Rdata")
head(Data)

N <- length(unique(Data$Area))
t <- length(unique(Data$Year))


# Prepare the data
Data.INLA <- data.frame(O=Data$Counts, E=Data$Expected,
                        ID.area=Data$Area, ID.year=Data$Year, 
                        ID.area.year=seq(1,N*t))


# Spatial neighborhood matrix: Rs

# RW1 precision matrix
Dm <- diff(diag(t),differences=1)
Rt <- t(Dm)%*%Dm




################################################
##  H1: Define the Uniform(0,inf) hyperprior  ##
################################################

sdunif="expression:
  logdens=-log_precision/2;
  return(logdens)"




#######################
##  Load the models  ##
#######################

source("bym_models/bym_h1_inla.R")





##########################
##  BYM model - Type I  ##
##########################

model <- bym.type1.h1



# Run the model
result.inla <- inla(model, family="poisson", data=Data.INLA, E=E,
                                  control.predictor=list(compute=TRUE, cdf=c(log(1))),
                                  control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                                  control.inla=list(strategy="laplace"))

# Save results
save(t.inla, result.inla, file = "bym_typeI_model_h1_inla.Rdata")





#############################################
##  BYM model - Type II/Type III/Type IV   ##
#############################################

model <- bym.types.h1


# Constraints Type II
R <- kronecker(Rt, diag(N))
rank.def <- N
A.constr <- kronecker(matrix(1,1,t), diag(N))
A.constr <- A.constr[-1,]
e <- rep(0,N-1)


# Constraints Type III
R <- kronecker(diag(t), Rs)
rank.def <- t
A.constr <- kronecker(diag(t),matrix(1,1,N))
A.constr <- A.constr[-1,]
e <- rep(0,t-1)


# Constraints Type IV
R <- kronecker(Rt, Rs)
rank.def <- N+t-1
A1 <- kronecker(matrix(1,1,t),diag(N))
A2 <- kronecker(diag(t),matrix(1,1,N))
A.constr <- rbind(A1[-1,], A2[-1,])
e <- rep(0, N+t-2)


# Run the model
result.inla <- inla(model, family="poisson", data=Data.INLA, E=E,
                                  control.predictor=list(compute=TRUE, cdf=c(log(1))),
                                  control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                                  control.inla=list(strategy="laplace"))


# Save results
save(result.inla, file = "bym_typeIV_model_h1_inla.Rdata")   # Type IV



