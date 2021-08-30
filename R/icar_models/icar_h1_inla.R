
# --------------         Spanish breast cancer data - INLA     -------------- #
# -------   Spatio-temporal models: ICAR - RW1 - Type I, II, III, IV  ------- #
# --------------                H1 hyperpriors                 -------------- #



###########################
##  ICAR - RW1 - Type I  ##
###########################
icar.type1.h1 <- O ~ 1 + f(ID.area, model="besag", graph=Rs, scale.model=FALSE, 
                           constr=TRUE, hyper=list(prec=list(prior=sdunif))) + 
                         f(ID.year, model="rw1", constr=TRUE, 
                           hyper=list(prec=list(prior=sdunif))) +
                         f(ID.area.year, model="iid", constr=TRUE, 
                           hyper=list(prec=list(prior=sdunif)))




#############################################
##  ICAR - RW1 - Type II/Type III/Type IV  ##
#############################################
icar.types.h1 <- O ~ 1 + f(ID.area, model="besag", graph=Rs, scale.model=FALSE, 
                           constr=TRUE, hyper=list(prec=list(prior=sdunif))) + 
                         f(ID.year, model="rw1", constr=TRUE,
                           hyper=list(prec=list(prior=sdunif))) +
                         f(ID.area.year, model="generic0", Cmatrix=R, 
                           rankdef=rank.def, constr=TRUE,
                           extraconstr=list(A=A.constr, e=e), 
                           hyper=list(prec=list(prior=sdunif)))

