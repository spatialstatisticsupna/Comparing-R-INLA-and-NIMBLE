
# --------------         Spanish breast cancer data - INLA     -------------- #
# -------   Spatio-temporal models: ICAR - RW1 - Type I, II, III, IV  ------- #
# --------------               H1 hyperpriors                  -------------- #


# Load INLA package
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

# H1 hyperpriors
source("icar_models/icar_h1_inla.R")
  



############################
##  ICAR model - Type I   ##
############################
model <- icar.type1.h1

# Run the model
typeI.inla <- inla(model, family="poisson", data=Data.INLA, E=E,
                    control.predictor=list(compute=TRUE, cdf=c(log(1))),
                    control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                    control.inla=list(strategy="laplace"))

# Save results
save(typeI.inla, file = "icar_typeI_model_h1_inla.Rdata")




##############################################
##  ICAR model - Type II/Type III/Type IV   ##
##############################################
model <- icar.types.h1   # This model permits to fit Type II, III and IV interactions


# The following code defines the constraints for each type of space-time interactions.
# Execute the corresponding code for your model.

# Type II
R <- kronecker(Rt, diag(N))
rank.def <- N
A.constr <- kronecker(matrix(1,1,t), diag(N))
A.constr <- A.constr[-1,]
e <- rep(0,N-1)


# Type III
R <- kronecker(diag(t), Rs)
rank.def <- t
A.constr <- kronecker(diag(t),matrix(1,1,N))
A.constr <- A.constr[-1,]
e <- rep(0,t-1)


# Type IV
R <- kronecker(Rt, Rs)
rank.def <- N+t-1
A1 <- kronecker(matrix(1,1,t),diag(N))
A2 <- kronecker(diag(t),matrix(1,1,N))
A.constr <- rbind(A1[-1,], A2[-1,])
e <- rep(0, N+t-2)


# Define linear combinations in INLA to obtain the posterior distribution of
# the intercept + temporal random effects
lc <- inla.make.lincombs("(Intercept)"=matrix(1,t,1), ID.year=diag(t))



# Run the model
typeIV.inla <- inla(model, family="poisson", data=Data.INLA, E=E, lincomb=lc,
                    control.predictor=list(compute=TRUE, cdf=c(log(1))),
                    control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                    control.inla=list(strategy="laplace"))

# Save results
save(typeIV.inla, file = "icar_typeIV_model_h1_inla.Rdata")   # Type IV



