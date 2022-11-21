
source("simulation-study.R") 


##############################
# Situation of widehat(pi)=pi
##############################
 


simula(Nrep=1000,ene=100,mu0_D=2,mu_D=4,sigma_D=1,mu0_H=0.5,mu_H=1,sigma_H=sqrt(8/3), sigmaX=1/3, tipo.mis_H=2,tipo.mis_D=2,estim_p_H="preal",estim_p_D="preal",regre_D="lineal",regre_H="lineal")


simula(Nrep=1000,ene=100,mu0_D=2,mu_D=4,sigma_D=1,mu0_H=0.5,mu_H=1,sigma_H=sqrt(8/3), sigmaX=1/3, tipo.mis_H=3,tipo.mis_D=3,estim_p_H="preal",estim_p_D="preal",regre_D="lineal",regre_H="lineal")


simula(Nrep=1000,ene=100,mu0_D=2,mu_D=4,sigma_D=1,mu0_H=0.5,mu_H=1,sigma_H=sqrt(8/3), sigmaX=1/3, tipo.mis_H=3,tipo.mis_D=2,estim_p_H="preal",estim_p_D="preal",regre_D="lineal",regre_H="lineal")


simula(Nrep=1000,ene=100,mu0_D=2,mu_D=4,sigma_D=1,mu0_H=0.5,mu_H=1,sigma_H=sqrt(8/3), sigmaX=1/3, tipo.mis_H=2,tipo.mis_D=3,estim_p_H="preal",estim_p_D="preal",regre_D="lineal",regre_H="lineal")



##############################
# Situation of widehat(pi) 
# using a logistic model
##############################
 
simula(Nrep=1000,ene=100,mu0_D=2,mu_D=4,sigma_D=1,mu0_H=0.5,mu_H=1,sigma_H=sqrt(8/3), sigmaX=1/3, tipo.mis_H=2,tipo.mis_D=2,estim_p_H="plogit",estim_p_D="plogit",regre_D="lineal",regre_H="lineal")


simula(Nrep=1000,ene=100,mu0_D=2,mu_D=4,sigma_D=1,mu0_H=0.5,mu_H=1,sigma_H=sqrt(8/3), sigmaX=1/3, tipo.mis_H=3,tipo.mis_D=3,estim_p_H="plogit",estim_p_D="plogit",regre_D="lineal",regre_H="lineal")


simula(Nrep=1000,ene=100,mu0_D=2,mu_D=4,sigma_D=1,mu0_H=0.5,mu_H=1,sigma_H=sqrt(8/3), sigmaX=1/3, tipo.mis_H=3,tipo.mis_D=2,estim_p_H="plogit",estim_p_D="plogit",regre_D="lineal",regre_H="lineal")


simula(Nrep=1000,ene=100,mu0_D=2,mu_D=4,sigma_D=1,mu0_H=0.5,mu_H=1,sigma_H=sqrt(8/3), sigmaX=1/3, tipo.mis_H=2,tipo.mis_D=3,estim_p_H="plogit",estim_p_D="plogit",regre_D="lineal",regre_H="lineal")



##############################
# Situation of widehat(pi)
# estimated using MCAR
##############################
 
simula(Nrep=1000,ene=100,mu0_D=2,mu_D=4,sigma_D=1,mu0_H=0.5,mu_H=1,sigma_H=sqrt(8/3), sigmaX=1/3, tipo.mis_H=2,tipo.mis_D=2,estim_p_H="pcte",estim_p_D="pcte",regre_D="lineal",regre_H="lineal")


simula(Nrep=1000,ene=100,mu0_D=2,mu_D=4,sigma_D=1,mu0_H=0.5,mu_H=1,sigma_H=sqrt(8/3), sigmaX=1/3, tipo.mis_H=3,tipo.mis_D=3,estim_p_H="pcte",estim_p_D="pcte",regre_D="lineal",regre_H="lineal")


simula(Nrep=1000,ene=100,mu0_D=2,mu_D=4,sigma_D=1,mu0_H=0.5,mu_H=1,sigma_H=sqrt(8/3), sigmaX=1/3, tipo.mis_H=3,tipo.mis_D=2,estim_p_H="pcte",estim_p_D="pcte",regre_D="lineal",regre_H="lineal")


simula(Nrep=1000,ene=100,mu0_D=2,mu_D=4,sigma_D=1,mu0_H=0.5,mu_H=1,sigma_H=sqrt(8/3), sigmaX=1/3, tipo.mis_H=2,tipo.mis_D=3,estim_p_H="pcte",estim_p_D="pcte",regre_D="lineal",regre_H="lineal")


############################################################
# Diseased population \widehat(pi) using a logistic model
# Healthy population of \widehat(pi) estimated using MCAR
############################################################

simula(Nrep=1000,ene=100,mu0_D=2,mu_D=4,sigma_D=1,mu0_H=0.5,mu_H=1,sigma_H=sqrt(8/3), sigmaX=1/3, tipo.mis_H=2,tipo.mis_D=2,estim_p_H="pcte",estim_p_D="plogit",regre_D="lineal",regre_H="lineal")


simula(Nrep=1000,ene=100,mu0_D=2,mu_D=4,sigma_D=1,mu0_H=0.5,mu_H=1,sigma_H=sqrt(8/3), sigmaX=1/3, tipo.mis_H=3,tipo.mis_D=3,estim_p_H="pcte",estim_p_D="plogit",regre_D="lineal",regre_H="lineal")


simula(Nrep=1000,ene=100,mu0_D=2,mu_D=4,sigma_D=1,mu0_H=0.5,mu_H=1,sigma_H=sqrt(8/3), sigmaX=1/3, tipo.mis_H=3,tipo.mis_D=2,estim_p_H="pcte",estim_p_D="plogit",regre_D="lineal",regre_H="lineal")


simula(Nrep=1000,ene=100,mu0_D=2,mu_D=4,sigma_D=1,mu0_H=0.5,mu_H=1,sigma_H=sqrt(8/3), sigmaX=1/3, tipo.mis_H=2,tipo.mis_D=3,estim_p_H="pcte",estim_p_D="plogit",regre_D="lineal",regre_H="lineal")



############################################################
# \widehat(pi) using true propensity
# Misspecification of the regression function
############################################################ 


simula(Nrep=1000,ene=100,mu0_D=2,mu_D=4,sigma_D=1,mu0_H=0.5,mu_H=1,sigma_H=sqrt(8/3), sigmaX=1/3, tipo.mis_H=2,tipo.mis_D=2,estim_p_H="preal",estim_p_D="preal",regre_D="nointer",regre_H="nointer")

simula(Nrep=1000,ene=100,mu0_D=2,mu_D=4,sigma_D=1,mu0_H=0.5,mu_H=1,sigma_H=sqrt(8/3), sigmaX=1/3, tipo.mis_H=3,tipo.mis_D=3,estim_p_H="preal",estim_p_D="preal",regre_D="nointer",regre_H="nointer")

simula(Nrep=1000,ene=100,mu0_D=2,mu_D=4,sigma_D=1,mu0_H=0.5,mu_H=1,sigma_H=sqrt(8/3), sigmaX=1/3, tipo.mis_H=3,tipo.mis_D=2,estim_p_H="preal",estim_p_D="preal",regre_D="nointer",regre_H="nointer")

simula(Nrep=1000,ene=100,mu0_D=2,mu_D=4,sigma_D=1,mu0_H=0.5,mu_H=1,sigma_H=sqrt(8/3), sigmaX=1/3, tipo.mis_H=2,tipo.mis_D=3,estim_p_H="preal",estim_p_D="preal",regre_D="nointer",regre_H="nointer")


############################################################
# \widehat(pi) using logistic model
# Misspecification of the regression function
############################################################ 

 
simula(Nrep=1000,ene=100,mu0_D=2,mu_D=4,sigma_D=1,mu0_H=0.5,mu_H=1,sigma_H=sqrt(8/3), sigmaX=1/3, tipo.mis_H=2,tipo.mis_D=2,estim_p_H="plogit",estim_p_D="plogit",regre_D="nointer",regre_H="nointer")

simula(Nrep=1000,ene=100,mu0_D=2,mu_D=4,sigma_D=1,mu0_H=0.5,mu_H=1,sigma_H=sqrt(8/3), sigmaX=1/3, tipo.mis_H=3,tipo.mis_D=3,estim_p_H="plogit",estim_p_D="plogit",regre_D="nointer",regre_H="nointer")

simula(Nrep=1000,ene=100,mu0_D=2,mu_D=4,sigma_D=1,mu0_H=0.5,mu_H=1,sigma_H=sqrt(8/3), sigmaX=1/3, tipo.mis_H=3,tipo.mis_D=2,estim_p_H="plogit",estim_p_D="plogit",regre_D="nointer",regre_H="nointer")

simula(Nrep=1000,ene=100,mu0_D=2,mu_D=4,sigma_D=1,mu0_H=0.5,mu_H=1,sigma_H=sqrt(8/3), sigmaX=1/3, tipo.mis_H=2,tipo.mis_D=3,estim_p_H="plogit",estim_p_D="plogit",regre_D="nointer",regre_H="nointer")



