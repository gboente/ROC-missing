rm(list=ls())


#################################################
# Auxiliary Functions
#################################################

source("functions_FY.R")

#################################################
#Required Library
#################################################

library(np)

#################################################
#Load Data
#################################################

load("autobomile.RData")
attach(datos)

#################################################
# Define Variables
# riesgo: CLASSIFICATION VARIABLE
# delta: missing indicator
# X: Covariates using to predict drop-outs
# V: Factor used to predict drop-outs
# Z: Covariates used to fit a regression model to the biomarker
# Y: Biomarker
#################################################

riesgo=rep(0,length(V1))
riesgo[V1 <=0]<- 1

delta=as.numeric(!is.na(V2))

X=cbind(V12,V17)
summary(X)
 
V= V7

Z=cbind(V12,V13,V21)
summary(Z)

#################################################
# Define Definitive Variables
#################################################

Y_D=V2[riesgo==0]
Y_H=V2[riesgo==1]
delta_D=delta[riesgo==0]
delta_H=delta[riesgo==1]
x_D=X[riesgo==0,]
x_H=X[riesgo==1,]

V_D= V[riesgo==0]
V_H= V[riesgo==1]


z_D=Z[riesgo==0,]
z_H=Z[riesgo==1,]

#################################################
# Consider  a Logistic Model
#################################################
a=glm(delta_H~x_H[,1]+x_H[,2]+V_H,family="binomial") 

py_H=unname(a$fitted.values)


b=glm(delta_D~x_D[,1]+x_D[,2]+V_D,family="binomial")
py_D=unname(b$fitted.values)

summary(a)
summary(b)
 

#################################################
# Fit a Regression Linear Model
# Using Data at Hand 
#################################################

y_D_1=Y_D[delta_D==1]
y_H_1=Y_H[delta_H==1]
z_D_1=z_D[delta_D==1,]
z_H_1=z_H[delta_H==1,]

#################################################
# Prepare Data to Use predict.lm
#################################################

aux=z_D_1[,1]
aux2=z_D_1[,2]
aux3=z_D_1[,3]
regresion_Y_D <- lm(y_D_1~aux+aux2+aux3)
m.est_D= predict.lm(regresion_Y_D,newdata=data.frame(aux=z_D[,1],aux2=z_D[,2],aux3=z_D[,3]))

#################################################
# Healthy Sample
#################################################

aux=z_H_1[,1]
aux2=z_H_1[,2]
aux3=z_H_1[,3]
regresion_Y_H <- lm(y_H_1~aux+aux2+aux3)
m.est_H= predict.lm(regresion_Y_H,newdata=data.frame(aux=z_H[,1],aux2=z_H[,2],aux3=z_H[,3]))

#################################################
#Summary Measures
#################################################

summary(regresion_Y_D)
summary(regresion_Y_H)
summary(m.est_D)
summary(m.est_H)

#################################################
# Grid of p Values
#################################################

grilla_p<- seq(0.01,1-0.01,0.01)

#################################################
# ROC estimators
#################################################

roc_KER_grilla=rep(NA,length(grilla_p))
roc_IPW_grilla=rep(NA,length(grilla_p))
roc_MUL_grilla=rep(NA,length(grilla_p))


for(i in 1:length(grilla_p)){

  roc_KER_grilla[i]<-roc_KER(y_D=Y_D,y_H=Y_H,delta_D,delta_H,pyx_D=py_D,pyx_H=py_H, p=grilla_p[i])
  roc_IPW_grilla[i]<-roc_IPW(y_D=Y_D,y_H=Y_H,delta_D,delta_H,pyx_D=py_D,pyx_H=py_H, p=grilla_p[i])
  roc_MUL_grilla[i]<-roc_MUL(y_D=Y_D,y_H=Y_H,delta_D,delta_H,pyx_D=py_D,pyx_H=py_H,  m.est_D,m.est_H,p=grilla_p[i])

}

pyx_D=rep(1,length(py_D))
pyx_H=rep(1,length(py_H))
 

#################################################
# ROC Computed With Data at Hand
#################################################
 


pyx_D_1=rep(1,length(y_D_1))
pyx_H_1=rep(1,length(y_H_1))

delta_D_1=rep(1,length(y_D_1))
delta_H_1=rep(1,length(y_H_1))

roc_simp_grilla=rep(NA,length(grilla_p))

for(i in 1:length(grilla_p)){
  roc_simp_grilla[i]<-roc_IPW(y_D=y_D_1,y_H=y_H_1,delta_D_1,delta_H_1,pyx_D=pyx_D_1,pyx_H=pyx_H_1, p=grilla_p[i])
}

#################################################
# AUC Estimators
#################################################

mean(roc_KER_grilla)
mean(roc_IPW_grilla)
mean(roc_MUL_grilla)
mean(roc_simp_grilla) 

#################################################
# Plots
#################################################

res=data.frame(grilla_p,roc_IPW_grilla,roc_KER_grilla,roc_MUL_grilla,roc_simp_grilla)
library(ggplot2)
cols <- c("IPW"="red","KER"="blue","MUL"="green","SIMP"="black")

ggplot(data=res, aes(x=grilla_p))+
  geom_line(aes(y=roc_IPW_grilla,color="IPW"))+
  geom_line(aes(y=roc_KER_grilla,color="KER"))+
  geom_line(aes(y=roc_MUL_grilla,color="MUL"))+
  geom_line(aes(y=roc_simp_grilla,color="SIMP"),lty=2,lwd=1)+
  scale_colour_manual(name="Estimations",values=cols) +
  ylab("True Positive Rate") + xlab("False Positive Rate") +
  geom_line(aes(y=grilla_p,color="black"))


pdf("ROC.pdf", bg='transparent')
par(mar=c(4.5,4.5,3,3))

plot(grilla_p,roc_IPW_grilla, type="n", xlim=c(0,1), ylim=c(0,1), xlab="p", ylab=expression(hat(ROC)),cex.lab=1.3,cex=1.2)
lines(grilla_p,roc_simp_grilla, col="black",lwd=3,lty=2)
lines(grilla_p,roc_IPW_grilla, col="red",lwd=3,lty=1)
lines(grilla_p,roc_KER_grilla, col="blue",lwd=3,lty=1)
lines(grilla_p,roc_MUL_grilla, col="gray40",lwd=3,lty=1)
lines(grilla_p,grilla_p,type="l", col="gray50",lty=1,lwd=1)
dev.off()


pdf("ROC_IPW.pdf", bg='transparent')
par(mar=c(4.5,4.5,3,3))
plot(grilla_p,roc_IPW_grilla, type="n", xlim=c(0,1), ylim=c(0,1), xlab="p", ylab=expression(hat(ROC)),cex.lab=1.3,cex=1.2)
lines(grilla_p,roc_simp_grilla, col="black",lwd=3,lty=2)
lines(grilla_p,roc_IPW_grilla, col="gray40",lwd=3,lty=1)
lines(grilla_p,grilla_p,type="l", col="gray50",lty=1,lwd=1)
dev.off()

pdf("ROC_KER.pdf", bg='transparent')
par(mar=c(4.5,4.5,3,3))
plot(grilla_p,roc_IPW_grilla, type="n", xlim=c(0,1), ylim=c(0,1), xlab="p", ylab=expression(hat(ROC)),cex.lab=1.3,cex=1.2)
lines(grilla_p,roc_simp_grilla, col="black",lwd=3,lty=2)
lines(grilla_p,roc_KER_grilla, col="gray40",lwd=3,lty=1)
lines(grilla_p,grilla_p,type="l", col="gray50",lty=1,lwd=1)
dev.off()


pdf("ROC_CONV.pdf", bg='transparent')
par(mar=c(4.5,4.5,3,3))
plot(grilla_p,roc_IPW_grilla, type="n", xlim=c(0,1), ylim=c(0,1), xlab="p", ylab=expression(hat(ROC)),cex.lab=1.3,cex=1.2)
lines(grilla_p,roc_simp_grilla, col="black",lwd=3,lty=2)
lines(grilla_p,roc_MUL_grilla, col="gray40",lwd=3,lty=1)
lines(grilla_p,grilla_p,type="l", col="gray50",lty=1,lwd=1)
dev.off()

