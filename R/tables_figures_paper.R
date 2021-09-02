
# -------- Tables and figures of Spanish Breast Cancer data analysis -------- #
# --------            ICAR spatial prior, H1 hyperpriors             -------- #



# Before running this code it is necessary to fit the spatio-temporal models 
# with Type I, Type II, Type III and Type IV interactions with R-INLA, Nimble 1      
# and Nimble 2.


rm(list=ls())
# Load packages
library(INLA)
library(RColorBrewer)
library(tmap)

# Set the working directory
setwd("")




#####################################
##  Load the data and the results  ##
#####################################

# Load .Rdata file that containts all the spatio-temporal models fitted.
# The name of each model should be typex.method where
# x=I,II,III or IV and method=inla,nimble1 or nimble2

load("")

Models.INLA <- list(typeI.inla, typeII.inla, typeIII.inla, typeIV.inla)
Models.NIMBLE1 <- list(typeI.nimble1, typeII.nimble1, typeIII.nimble1, typeIV.nimble1)
Models.NIMBLE2 <- list(typeI.nimble2, typeII.nimble2, typeIII.nimble2, typeIV.nimble2)


# Load Spanish breast cancer data
load("BreastCancer_data.Rdata")
t <- 21
n <- 50
nt <- n*t




###############################
##  Model selection criteria  #
###############################

## INLA ##
Table2.INLA <- data.frame(mean.deviance=round(unlist(lapply(Models.INLA, function(x) x$dic$mean.deviance)), 4),
                          p.eff=round(unlist(lapply(Models.INLA, function(x) x$dic$p.eff)), 4),
                          DIC=round(unlist(lapply(Models.INLA, function(x) x$dic$dic)), 4),
                          WAIC=round(unlist(lapply(Models.INLA, function(x) x$waic$waic)), 4),
                          row.names=c("Type I","Type II","Type III","Type IV"))
print(Table2.INLA)



## NIMBLE ##
DIC.nimble <- function(summary, result){
  dev <- summary[grep('sumDev', rownames(summary), value=TRUE), -c(2)]
  lambda <- summary[grep('lambda', rownames(summary), value=TRUE), -c(2)]
  
  dev2 <- 0
  for (i in 1:nt){
    dev2 <- dev2 -2*(-lambda[i, 1] + Data$Counts[i]*log(lambda[i, 1])-lfactorial(Data$Counts[i]))
  }
  pD <- dev[1] - dev2
  data.frame(mean.dev=dev[1],
             pD=dev[1] - dev2,     # Effective number of parameters
             DIC=dev[1]+pD,        
             WAIC=result$WAIC[1])       
}


# Nimble 1
Table2.NIMBLE1 <- lapply(Models.NIMBLE1, function(x) DIC.nimble(x$summary$all.chains, x))
Table2.NIMBLE1 <- data.frame(matrix(unlist(Table2.NIMBLE1), 4, 4, byrow=T),
                             row.names = c("Type I","Type II","Type III","Type IV"))
colnames(Table2.NIMBLE1) <- c("mean.deviance", "p.eff", "DIC", "WAIC")

print(Table2.NIMBLE1)
  

# Nimble 2
Table2.NIMBLE2 <- lapply(Models.NIMBLE2, function(x) DIC.nimble(x$summary$all.chains, x))
Table2.NIMBLE2 <- data.frame(matrix(unlist(Table2.NIMBLE2), 4, 4, byrow=T),
                             row.names = c("Type I","Type II","Type III","Type IV"))
colnames(Table2.NIMBLE2) <- c("mean.deviance", "p.eff", "DIC", "WAIC")

print(Table2.NIMBLE2)




# From now on all the tables and figures are for the model with Type IV 
# spatio-temporal interaction.

###############################################
##  Posterior means and standard deviations  ##
###############################################

## INLA ##

model <- typeIV.inla

marg.sd.spat <- inla.tmarginal(function(x) sqrt(1/x), model$marginals.hyperpar$`Precision for ID.area`)
marg.sd.temp <- inla.tmarginal(function(x) sqrt(1/x), model$marginals.hyperpar$`Precision for ID.year`)
marg.sd.st <- inla.tmarginal(function(x) sqrt(1/x), model$marginals.hyperpar$`Precision for ID.area.year`)

summary.INLA <- function(marginal){
  data.frame(mean=inla.emarginal(function(x) x, marginal),
             sd=sqrt(inla.emarginal(function(x) x^2, marginal)-(inla.emarginal(function(x) x, marginal))^2),
             median=inla.qmarginal(c(0.5), marginal),
             CI.L=inla.qmarginal(c(0.025), marginal),
             CI.U=inla.qmarginal(c(0.975), marginal))
}


Table3.INLA <- rbind(model$summary.fixed[, 1:2],
                     summary.INLA(marg.sd.spat)[, 1:2],
                     summary.INLA(marg.sd.temp)[, 1:2],
                     summary.INLA(marg.sd.st)[, 1:2])

Table3.INLA <- data.frame(round(Table3.INLA, 4), 
                          row.names = c("intercept", "sd.spat", "sd.temp", "sd.st"))

print(Table3.INLA)



## NIMBLE ##

summary.NIMBLE <- function(model, name){
  summary <- model$summary$all.chains[grep(name, rownames(model$summary$all.chains), value=TRUE), ]
  data.frame(mean=summary[1],
             sd=summary[3],
             median=summary[2],
             CI.L=summary[4], 
             CI.U=summary[5])
}


# NIMBLE 1
alpha0.nimble1 <- typeIV.nimble1$summary$all.chains[grep("alpha0", rownames(typeIV.nimble1$summary$all.chains), value=TRUE), ]

Table3.NIMBLE1 <- rbind(alpha0.nimble1[22, c(1,3)],
                        summary.NIMBLE(typeIV.nimble1, "sd.u")[, 1:2],
                        summary.NIMBLE(typeIV.nimble1, "sd.temp")[, 1:2],
                        summary.NIMBLE(typeIV.nimble1, "sd.t4")[, 1:2])

Table3.NIMBLE1 <- data.frame(round(Table3.NIMBLE1, 4), 
                          row.names = c("intercept", "sd.spat", "sd.temp", "sd.st"))

print(Table3.NIMBLE1)


# NIMBLE 2
alpha0.nimble2 <- typeIV.nimble2$summary$all.chains[grep("alpha0", rownames(typeIV.nimble2$summary$all.chains), value=TRUE), ]

Table3.NIMBLE2 <- rbind(alpha0.nimble2[22, c(1,3)],
                        summary.NIMBLE(typeIV.nimble2, "sd.u")[, 1:2],
                        summary.NIMBLE(typeIV.nimble2, "sd.temp")[, 1:2],
                        summary.NIMBLE(typeIV.nimble2, "sd.t4")[, 1:2])

Table3.NIMBLE2 <- data.frame(round(Table3.NIMBLE2, 4), 
                             row.names = c("intercept", "sd.spat", "sd.temp", "sd.st"))

print(Table3.NIMBLE2)




#################################################
## Spatial patterns (posterior mean estimates) ##
#################################################

Carto_ESP$spat.inla <- unlist(lapply(typeIV.inla$marginals.random$ID.area, function(x) inla.emarginal(exp,x)))
Carto_ESP$spat.nimble1 <- typeIV.nimble1$summary$all.chains[grep('Espat', rownames(typeIV.nimble1$summary$all.chains), value=TRUE), 1]
Carto_ESP$spat.nimble2 <- typeIV.nimble2$summary$all.chains[grep('Espat', rownames(typeIV.nimble2$summary$all.chains), value=TRUE), 1]


color.pal <- brewer.pal(6,"Oranges")[1:6]
values <- c(0.77,0.83,0.91,1,1.1,1.2,1.3)


Spat.pattern <- tm_shape(Carto_ESP) +
  tm_polygons(col=c("spat.inla", "spat.nimble1", "spat.nimble2"), palette=color.pal, title="", 
              legend.show=T, legend.reverse=T, style="fixed", breaks=values, interval.closure="left") +
  tm_layout(main.title="", main.title.position="center",
            legend.outside=T, legend.outside.position="right", legend.frame=F, legend.outside.size=0.2,
            panel.labels=c("R-INLA", "Nimble 1", "Nimble 2"),
            panel.label.size = 1.7,
            outer.margins=c(0.02,0.01,0.02,0.01), legend.text.size=1.4) + 
  tm_facets(nrow=1, ncol=3)

tmap_mode("plot")
print(Spat.pattern)




#############################
## Temporal trends (gamma) ##
#############################

temp.inla <- unlist(lapply(typeIV.inla$marginals.random$ID.year, function(x) inla.emarginal(exp,x)))
temp.nimble1 <- typeIV.nimble1$summary$all.chains[grep('Etemp', rownames(typeIV.nimble1$summary$all.chains), value=TRUE), ]
temp.nimble2 <- typeIV.nimble2$summary$all.chains[grep('Etemp', rownames(typeIV.nimble2$summary$all.chains), value=TRUE), ]


t.from <- 1990
t.to <- 2010
x <- 1:t


# INLA
aux.inla <- lapply(typeIV.inla$marginals.random$ID.year, function(x) inla.tmarginal(exp,x))
q1.inla <- unlist(lapply(aux.inla, function(x) inla.qmarginal(0.025,x)))
q2.inla <- unlist(lapply(aux.inla, function(x) inla.qmarginal(0.975,x)))


# NIMBLE 1
q1.nimble1 <- temp.nimble1[, 4]
q2.nimble1 <- temp.nimble1[, 5]


# NIMBLE 2
q1.nimble2 <- temp.nimble2[, 4]
q2.nimble2 <- temp.nimble2[, 5]



par(mfrow=c(1,3), mar=c(5.1, 4.5, 4.1, 2.1))
# INLA
plot(range(x),c(0.7,1.3),type="n",xlab="year",ylab=expression(exp(gamma[t])), xaxt="n", main="R-INLA", 
     cex.lab=1.3, cex.axis=1.4, cex.main=1.8, cex.sub=1.5)
axis(1, at=seq(1,t,4), labels=seq(t.from,t.to,4), las=0)
X.Vec <- c(x, tail(x, 1), rev(x), x[1])
Y.Vec <- c(q1.inla, tail(q2.inla, 1), rev(q2.inla), q1.inla[1])
polygon(X.Vec, Y.Vec, col = "grey", border = NA)
lines(temp.inla)
abline(h=1,lty=2)
# NIMBLE 1
plot(range(x),c(0.7,1.3),type="n",xlab="year",ylab=expression(exp(gamma[t])), xaxt="n", main="Nimble 1",
     cex.lab=1.3, cex.axis=1.2, cex.main=1.8, cex.sub=1.5)
axis(1, at=seq(1,t,4), labels=seq(t.from,t.to,4), las=0)
X.Vec <- c(x, tail(x, 1), rev(x), x[1])
Y.Vec <- c(q1.nimble1, tail(q2.nimble1, 1), rev(q2.nimble1), q1.nimble1[1])
polygon(X.Vec, Y.Vec, col = "grey", border = NA)
lines(temp.nimble1[, 1])
abline(h=1,lty=2)
# NIMBLE 2
plot(range(x),c(0.7,1.3),type="n",xlab="year",ylab=expression(exp(gamma[t])), xaxt="n", main="Nimble 2",
     cex.lab=1.3, cex.axis=1.2, cex.main=1.8, cex.sub=1.5)
axis(1, at=seq(1,t,4), labels=seq(t.from,t.to,4), las=0)
X.Vec <- c(x, tail(x, 1), rev(x), x[1])
Y.Vec <- c(q1.nimble2, tail(q2.nimble2, 1), rev(q2.nimble2), q1.nimble2[1])
polygon(X.Vec, Y.Vec, col = "grey", border = NA)
lines(temp.nimble2[, 1])
abline(h=1,lty=2)

dev.off()




######################################
## Temporal trends (alpha0 + gamma) ##
######################################

alpha0.temp.inla <- unlist(lapply(typeIV.inla$marginals.lincomb.derived, function(x) inla.emarginal(exp,x)))
alpha0.temp.nimb1 <- typeIV.nimble1$summary$all.chains[grep('Ealpha0.temp', rownames(typeIV.nimble1$summary$all.chains), value=TRUE), ]
alpha0.temp.nimb2 <- typeIV.nimble2$summary$all.chains[grep('Ealpha0.temp', rownames(typeIV.nimble2$summary$all.chains), value=TRUE), ]


t.from <- 1990
t.to <- 2010
x <- 1:t


# INLA
aux.inla.2 <- lapply(typeIV.inla$marginals.lincomb.derived, function(x) inla.tmarginal(exp,x))
q1.inla.2 <- unlist(lapply(aux.inla.2, function(x) inla.qmarginal(0.025,x)))
q2.inla.2 <- unlist(lapply(aux.inla.2, function(x) inla.qmarginal(0.975,x)))


# NIMBLE 1
q1.nimble1.2 <- alpha0.temp.nimb1[, 4]
q2.nimble1.2 <- alpha0.temp.nimb1[, 5]


# NIMBLE 2
q1.nimble2.2 <- alpha0.temp.nimb2[, 4]
q2.nimble2.2 <- alpha0.temp.nimb1[, 5]



par(mfrow=c(1,3), mar=c(5.1, 4.5, 4.1, 2.1))
# INLA
plot(range(x),c(0.7,1.3),type="n",xlab="year",ylab=expression(exp(alpha [0]+gamma[t])), xaxt="n", main="R-INLA", 
     cex.lab=1.3, cex.axis=1.4, cex.main=1.8, cex.sub=1.5)
axis(1, at=seq(1,t,4), labels=seq(t.from,t.to,4), las=0)
X.Vec <- c(x, tail(x, 1), rev(x), x[1])
Y.Vec <- c(q1.inla.2, tail(q2.inla.2, 1), rev(q2.inla.2), q1.inla.2[1])
polygon(X.Vec, Y.Vec, col = "grey", border = NA)
lines(alpha0.temp.inla)
abline(h=1,lty=2)
# NIMBLE 1
plot(range(x),c(0.7,1.3),type="n",xlab="year",ylab=expression(exp(alpha [0]+gamma[t])), xaxt="n", main="Nimble 1",
     cex.lab=1.3, cex.axis=1.2, cex.main=1.8, cex.sub=1.5)
axis(1, at=seq(1,t,4), labels=seq(t.from,t.to,4), las=0)
X.Vec <- c(x, tail(x, 1), rev(x), x[1])
Y.Vec <- c(q1.nimble1.2, tail(q2.nimble1.2, 1), rev(q2.nimble1.2), q1.nimble1.2[1])
polygon(X.Vec, Y.Vec, col = "grey", border = NA)
lines(alpha0.temp.nimb1[, 1])
abline(h=1,lty=2)
# NIMBLE 2
plot(range(x),c(0.7,1.3),type="n",xlab="year",ylab=expression(exp(alpha [0]+gamma[t])), xaxt="n", main="Nimble 2",
     cex.lab=1.3, cex.axis=1.2, cex.main=1.8, cex.sub=1.5)
axis(1, at=seq(1,t,4), labels=seq(t.from,t.to,4), las=0)
X.Vec <- c(x, tail(x, 1), rev(x), x[1])
Y.Vec <- c(q1.nimble2.2, tail(q2.nimble2.2, 1), rev(q2.nimble2.2), q1.nimble2.2[1])
polygon(X.Vec, Y.Vec, col = "grey", border = NA)
lines(alpha0.temp.nimb2[, 1])
abline(h=1,lty=2)

dev.off()




######################################################################
## Dispersion plot (posterior mean estimates of the relative risks) ##
######################################################################

RR.inla <- typeIV.inla$summary.fitted.values[, 1]
RR.nimble1 <- typeIV.nimble1$summary$all.chains[grep('RR', rownames(typeIV.nimble1$summary$all.chains), value=TRUE), 1]
RR.nimble2 <- typeIV.nimble2$summary$all.chains[grep('RR', rownames(typeIV.nimble2$summary$all.chains), value=TRUE), 1]


par(mfrow=c(1,2))
plot(RR.nimble1, RR.inla, xlab="Nimble 1", ylab="R-INLA",
     cex.axis=1, cex.lab=1.2)
abline(0,1)

plot(RR.nimble1, RR.nimble2, xlab="Nimble 1", ylab="Nimble 2",
     cex.axis=1, cex.lab=1.2)
abline(0,1)

dev.off()
