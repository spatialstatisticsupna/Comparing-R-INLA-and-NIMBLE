
# --------------         Spanish breast cancer data - INLA     -------------- #
# -------   Spatio-temporal models: BYM - RW1 - Type I, II, III, IV   ------- #
# --------------                 H1 hyperpriors                -------------- #


rm(list=ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))


# Load INLA package
library(INLA)


# Define the space-time interaction
st.interac <- "TypeI"
# st.interac <- "TypeII"
# st.interac <- "TypeIII"
# st.interac <- "TypeIV"



#################################
##  Load and prepare the data  ##
#################################

load("../BreastCancer_data.Rdata")
head(Data)

N <- length(unique(Data$Area))
t <- length(unique(Data$Year))


# Prepare the data
Data.INLA <- data.frame(O=Data$Counts, E=Data$Expected,
                        ID.area=Data$Area, ID.year=Data$Year, 
                        ID.area.year=seq(1,N*t))


# Spatial neighborhood matrix: Rs
Rs <- as(Rs, "Matrix")


# RW1 precision matrix
Dm <- diff(diag(t),differences=1)
Rt <- as(t(Dm)%*%Dm, "Matrix")



################################################
##  H1: Define the Uniform(0,inf) hyperprior  ##
################################################

sdunif="expression:
  logdens=-log_precision/2;
  return(logdens)"



#######################
##  Load the models  ##
#######################

source("../bym_models/bym_h1_inla.R")



##########################
##  BYM model - Type I  ##
##########################

if(st.interac=="TypeI") {
  model <- bym.type1.h1
  
  # Run the model
  typeI.inla <- inla(model, family="poisson", data=Data.INLA, E=E,
                     control.predictor=list(compute=TRUE, cdf=c(log(1))),
                     control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                     control.inla=list(strategy="laplace"))
  
  # Save results
  save(typeI.inla, file = "bym_typeI_model_h1_inla.Rdata")
}



############################################
##  BYM model - Type II/Type III/Type IV  ##
############################################

# The following code defines the constraints for each type of space-time interactions.
# Execute the corresponding code for your model.

if(st.interac=="TypeII") {
  R <- kronecker(Rt, diag(N))
  rank.def <- N
  A.constr <- kronecker(matrix(1,1,t), diag(N))
  A.constr <- A.constr[-1,]
  e <- rep(0,N-1)
}


if(st.interac=="TypeIII") {
  R <- kronecker(diag(t), Rs)
  rank.def <- t
  A.constr <- kronecker(diag(t),matrix(1,1,N))
  A.constr <- A.constr[-1,]
  e <- rep(0,t-1)
}


if(st.interac=="TypeIV") {
  R <- kronecker(Rt, Rs)
  rank.def <- N+t-1
  A1 <- kronecker(matrix(1,1,t),diag(N))
  A2 <- kronecker(diag(t),matrix(1,1,N))
  A.constr <- rbind(A1[-1,], A2[-1,])
  e <- rep(0, N+t-2)
}


# Define linear combinations in INLA to obtain the posterior distribution of
# the intercept + temporal random effects
lc <- inla.make.lincombs("(Intercept)"=matrix(1,t,1), ID.year=diag(t))



if(st.interac %in% c("TypeII", "TypeIII", "TypeIV")) {
  model <- bym.types.h1
  
  # Run the model
  types.inla <- inla(model, family="poisson", data=Data.INLA, E=E, lincomb=lc,
                      control.predictor=list(compute=TRUE, cdf=c(log(1))),
                      control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                      control.inla=list(strategy="laplace"))
}


# Save results
if(st.interac=="TypeII") {
  save(types.inla, file = "bym_typeII_model_h1_inla.Rdata")
}

if(st.interac=="TypeIII") {
  save(types.inla, file = "bym_typeIII_model_h1_inla.Rdata")
}

if(st.interac=="TypeIV") {
  save(types.inla, file = "bym_typeIV_model_h1_inla.Rdata")
}




