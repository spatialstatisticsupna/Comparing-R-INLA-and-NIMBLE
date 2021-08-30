
# --------------         Spanish breast cancer data - INLA     -------------- #
# -------   Spatio-temporal models: BYM - RW1 - Type I, II, III, IV   ------- #
# --------------                H1 hyperpriors                 -------------- #


##########################
##  BYM - RW1 - Type I  ##
##########################
bym.type1.h1 <- O ~ 1 + f(ID.area, model="bym", graph=Rs, scale.model=FALSE, 
                          constr=TRUE, hyper=list(prec.unstruct=list(prior=sdunif),
                          prec.spatial=list(prior=sdunif))) + 
                         f(ID.year, model="rw1", constr=TRUE, 
                            hyper=list(prec=list(prior=sdunif))) +
                         f(ID.area.year, model="iid", constr=TRUE, 
                            hyper=list(prec=list(prior=sdunif)))



#############################################
##  BYM - RW1 - Type II/Type III/Type IV   ##
#############################################
bym.types.h1 <- O ~ 1 + f(ID.area, model="bym", graph=Rs, scale.model=FALSE, 
                          constr=TRUE, hyper=list(prec.unstruct=list(prior=sdunif),
                          prec.spatial=list(prior=sdunif))) + 
                        f(ID.year, model="rw1", constr=TRUE,
                          hyper=list(prec=list(prior=sdunif))) +
                        f(ID.area.year, model="generic0", Cmatrix=R, 
                          rankdef=rank.def, constr=TRUE,
                          extraconstr=list(A=A.constr, e=e), 
                          hyper=list(prec=list(prior=sdunif)))

