source("functions_FY.R")

##################################
# Summary measures
##################################

DEV<- function(rocas,rocashat){
	logrocas<- ifelse(rocas>= 10^(-100),log(rocas),log(10^(-100)))       
	logrocas.m1<- ifelse(1-rocas>= 10^(-100),log(1-rocas),log(10^(-100)))       
	DEV<- -rocashat*logrocas-(1-rocashat)*logrocas.m1
	return(DEV)
}


DIF2<- function(rocas,rocashat){
	rocas<- ifelse(rocas>= 10^(-100),rocas,10^(-100))
	rocashat<- ifelse(rocashat>= 10^(-100),rocashat,10^(-100)) 
	DIF2<- (rocashat/rocas-1)**2
	return(DIF2)
}

 

##################################
# Drop-out Model
##################################

pobsY<-function(v,tipo.mis){
                t=v-0.5
                if(tipo.mis==0){pobs<-1}
		if(tipo.mis==1){pobs<-0.7}
	 	if(tipo.mis==2){		
				pobs<-1/(1+exp(2*t))
		}
	
		if(tipo.mis==3){		
				pobs<-exp(t/2)/(1+exp(t/2))
		}

		return(pobs)
}


##################################
# True ROC Curve
##################################

roc_real<-function(mu0_D,mu_D, mu0_H, mu_H, sigma_D, sigma_H,sigmaX, p){
	FHinv=qnorm(1-p)*sqrt(sigma_H^2+mu_H^2*sigmaX^2)+mu0_H
	sigma=sqrt(sigma_D^2+mu_D^2*sigmaX^2)
	resultado=1-pnorm((FHinv-mu0_D)/sigma)
	return(resultado)

}


##################################
# Data Generation
##################################

genero<-function(jota, ene=100,mu0_D=2,mu_D=4,sigma_D=2,mu0_H=0.5,mu_H=1,sigma_H=1.5, sigmaX=1/3, tipo.mis_D=2,tipo.mis_H=2){
		set.seed(123+jota)

		##########################################
		# Regression Model Y_D= mu(X) + sigma*eps
		##########################################
		 
  		x_D<-rnorm(ene,0,1)*sigmaX                    #Covariate x Diseased
  		x_H<-rnorm(ene,0,1)*sigmaX                    #Covariate x Healthy

  		eps_D <- rnorm(ene,0,1)                       #Errors Diseased  
  		eps_H <- rnorm(ene,0,1)                       #Errors Healthy

  		Y_D <- mu0_D + mu_D*x_D  + sigma_D*eps_D      #Biomarker Diseased
  		Y_H <- mu0_H + mu_H*x_H  + sigma_H*eps_H      #Biomarker Healthy

		##########################################
		# Generate Missing data 
		##########################################

		delta_D=delta_H=  py_D=py_H=rep(1,ene)

		u1=runif(ene)
		u2=runif(ene)

		for (i in 1: ene){
        		py_D[i]=pobsY(x_D[i],tipo.mis=tipo.mis_D)
			py_H[i]=pobsY(x_H[i],tipo.mis=tipo.mis_H)
		
			if (u1[i]>py_D[i]){delta_D[i]=0}
      
			if (u2[i]>py_H[i]){delta_H[i]=0}
		}
  
		return(list(Y_D=Y_D,Y_H=Y_H, x_D=x_D,x_H=x_H, delta_D=delta_D,delta_H=delta_H, py_TRUE_D=py_D, py_TRUE_H=py_H))

}



##################################
# Simulation function
##################################

simula<-function(Nrep=1000,ene=100,mu0_D=2,mu_D=4,sigma_D=2,mu0_H=0.5,mu_H=1,sigma_H=1.5, sigmaX=1/3, tipo.mis_H=2,tipo.mis_D=2,estim_p_H="preal",estim_p_D="preal",regre_D="lineal",regre_H="lineal"){

	################################################
	# Names of the files where to save the results
	################################################
	
 	ARCHIVOCOEF<-paste("regre_CL_n",ene,"_H_",estim_p_H,"_D_",estim_p_D, regre_H, regre_D,"_mis_H_",  tipo.mis_H,"_mis_D_",  tipo.mis_D,  ".txt",sep="") 

 	ARCHIVOIPW<-paste("IPW_n",ene,"_H_",estim_p_H,"_D_",estim_p_D, regre_H, regre_D,"_mis_H_",  tipo.mis_H,"_mis_D_",  tipo.mis_D, ".txt",sep="") 

 	ARCHIVOMUL<-paste("MUL_n",ene,"_H_",estim_p_H,"_D_",estim_p_D,regre_H, regre_D, "_mis_H_",  tipo.mis_H,"_mis_D_",  tipo.mis_D,".txt",sep="") 

 
 	ARCHIVOKER<-paste("KER_n",ene,"_H_",estim_p_H,"_D_",estim_p_D, "_mis_H_",  tipo.mis_H,"_mis_D_",  tipo.mis_D, "_dim1.txt",sep="") 

 
 	ARCHIVOROCKER<-paste("ROC_KER_n",ene,"_H_",estim_p_H,"_D_",estim_p_D,"_mis_H_",  tipo.mis_H,"_mis_D_",  tipo.mis_D, "_dim1.txt",sep="") 


 	ARCHIVOROCIPW<-paste("ROC_IPW_n",ene,"_H_",estim_p_H,"_D_",estim_p_D,regre_H, regre_D,"_mis_H_",  tipo.mis_H,"_mis_D_",  tipo.mis_D, ".txt",sep="") 

 	ARCHIVOROCMUL<-paste("ROC_MUL_n",ene,"_H_",estim_p_H,"_D_",estim_p_D,regre_H, regre_D,"_mis_H_",  tipo.mis_H,"_mis_D_",  tipo.mis_D, ".txt",sep="") 

   	grilla_p<- seq(0.01,1-0.01,0.01)
 	tiempo1=Sys.time()

	################################################
	# Start the replications
	################################################

  	for(j in 1:Nrep){

   		print(j)

   		
	##########################################
	# Generate data
	##########################################
		
		datos=genero(j, ene=ene,mu0_D,mu_D,sigma_D,mu0_H,mu_H,sigma_H, sigmaX, tipo.mis_D,tipo.mis_H)

		Y_D=datos$Y_D
		Y_H=datos$Y_H

		grilla_H=sort(Y_H)
		grilla_D=sort(Y_D)
 
		x_D=datos$x_D
		x_H=datos$x_H

 		delta_D=datos$delta_D
		delta_H=datos$delta_H

		ptrue_D=datos$py_TRUE_D 
		ptrue_H=datos$py_TRUE_H 


		Y_D[delta_D==0]<- NA

		Y_H[delta_H==0]<- NA

		y_D_1=Y_D[delta_D==1]
		y_H_1=Y_H[delta_H==1]

		x_D_1=x_D[delta_D==1]
		x_H_1=x_H[delta_H==1]
		
		################################################
		# Drop-out estimator for the 
		# IPW and Smoothed method
		################################################

		################################################
		# True propensity 
		################################################

		if(estim_p_H=="preal"){py_H=ptrue_H}
		if(estim_p_D=="preal"){py_D=ptrue_D}

		################################################
		# Logistic fit
		################################################
	
		if(estim_p_D=="plogit"){
		a=glm(delta_D~x_D,family="binomial")
		py_D=unname(a$fitted.values)
		}
		if(estim_p_H=="plogit"){
		a=glm(delta_H~x_H,family="binomial")
		py_H=unname(a$fitted.values)
		}

		################################################
		# MCAR model
		################################################

		if(estim_p_D=="pcte"){
		a=mean(delta_D)
		py_D=rep(a,ene)
		}
		if(estim_p_H=="pcte"){
		a=mean(delta_H)
		py_H=rep(a,ene)
		}


		######################################
		# Regression estimator
		# for the convolution based method
		######################################

		######################################
		# Estimated as y= a + b*x + epsilon
		######################################
 
		if(regre_D=="lineal"){
			regresion_Y_D <- lm(y_D_1~x_D_1)
			m.est_D= regresion_Y_D$coef[1]+ regresion_Y_D$coef[2]*x_D
			vec_D=regresion_Y_D$coef
		}
 		
		if(regre_H=="lineal"){
			regresion_Y_H <- lm(y_H_1~x_H_1)
	  		m.est_H= regresion_Y_H$coef[1]+ regresion_Y_H$coef[2]*x_H
			vec_H=regresion_Y_H$coef 
			}

		######################################
		# Misspecified Regression Model
		######################################

		######################################
		# Estimated as y= b*x + epsilon
		# that is without intercept
		######################################

		if(regre_D=="nointer"){
			regresion_Y_D <- lm(y_D_1~x_D_1-1)
			m.est_D= regresion_Y_D$coef[1]*x_D
			vec_D=regresion_Y_D$coef
		}
 		
		if(regre_H=="nointer"){
			regresion_Y_H <- lm(y_H_1~x_H_1-1)
	  		m.est_H= regresion_Y_H$coef[1]*x_H
			vec_H=regresion_Y_H$coef 
			}

 		################################################
		# Save the estimated regression parameters
		################################################
		
		vec.cl<-unname(c(j,vec_D,vec_H))

		lvec=length(vec.cl)
	
		write(vec.cl,file=ARCHIVOCOEF,ncolumns=lvec,append=T)
   
		################################################
		# Compute the ROC curve estimators
		################################################


  	  	roc_real_grilla=roc_IPW_grilla<-roc_MUL_grilla<-roc_KER_grilla<-rep(NA,length(grilla_p))
  
      		for(i in 1:length(grilla_p)){
			roc_real_grilla[i]<- roc_real(mu0_D,mu_D, mu0_H, mu_H, sigma_D, sigma_H,sigmaX, p=grilla_p[i]) 
        		roc_IPW_grilla[i]<- roc_IPW(y_D=Y_D,y_H=Y_H,delta_D,delta_H,pyx_D=py_D,pyx_H=py_H, p=grilla_p[i]) 
			roc_MUL_grilla[i]<- roc_MUL(y_D=Y_D,y_H=Y_H,delta_D,delta_H,pyx_D=py_D,pyx_H=py_H,m.est_D,m.est_H, p=grilla_p[i]) 
         		roc_KER_grilla[i]<- roc_KER(y_D=Y_D,y_H=Y_H,delta_D,delta_H,pyx_D=py_D,pyx_H=py_H, p=grilla_p[i]) 

       		}
     
 

		##############################################
		# Save the IPW estimator of the ROC
		##############################################

		roc_cl_grilla=roc_IPW_grilla
		vec.auc<- c(j, roc_cl_grilla)
   
		lvec=length(vec.auc)
	
     		write(t(vec.auc),file=ARCHIVOROCIPW,ncolumns=lvec,append=T)


		##############################################
		# Save the Convolution estimator of the ROC
		##############################################

		roc_cl_grilla=roc_MUL_grilla
		vec.auc<- c(j, roc_cl_grilla)
   
		lvec=length(vec.auc)
	
     		write(t(vec.auc),file=ARCHIVOROCMUL,ncolumns=lvec,append=T)


		###############################################
		# Save the Smoother  estimator of the ROC
		###############################################

		roc_cl_grilla=roc_KER_grilla
		vec.auc<- c(j, roc_cl_grilla)
   
		lvec=length(vec.auc)
	
     		write(t(vec.auc),file=ARCHIVOROCKER,ncolumns=lvec,append=T)

		###########################################
		# SUMMARY MEASURES
		########################################### 

		AUC_real=mean(roc_real_grilla)

		####################################################
		# Summary for the IPW estimator of the ROC and AUC
		####################################################
		
		roc_cl_grilla=roc_IPW_grilla
		MSE.cl <- mean((roc_real_grilla-roc_cl_grilla)*(roc_real_grilla-roc_cl_grilla))
      		DEV.cl <- mean(DEV(roc_real_grilla,roc_cl_grilla))
      		DIF2.cl  <- mean(DIF2(roc_real_grilla,roc_cl_grilla))
      		MAXDIF.cl  <-   max(abs(roc_real_grilla-roc_cl_grilla))
        	AUC<-mean(roc_cl_grilla)
		vec.auc<- c(j, MSE.cl ,DEV.cl , DIF2.cl  , MAXDIF.cl , AUC,AUC_real)
   
		lvec=length(vec.auc)
	
     		write(t(vec.auc),file=ARCHIVOIPW,ncolumns=lvec,append=T)

		############################################################
		# Summary for the Convolution estimator of the ROC and AUC
		############################################################

		roc_cl_grilla=roc_MUL_grilla
	 
		MSE.cl <- mean((roc_real_grilla-roc_cl_grilla)*(roc_real_grilla-roc_cl_grilla))
      		DEV.cl <- mean(DEV(roc_real_grilla,roc_cl_grilla))
      		DIF2.cl  <- mean(DIF2(roc_real_grilla,roc_cl_grilla))
      		MAXDIF.cl  <-   max(abs(roc_real_grilla-roc_cl_grilla))
        	AUC<-mean(roc_cl_grilla)
		vec.auc<- c(j, MSE.cl ,DEV.cl , DIF2.cl  , MAXDIF.cl , AUC,AUC_real)
   
   
		lvec=length(vec.auc)
	
     		write(t(vec.auc),file=ARCHIVOMUL,ncolumns=lvec,append=T)


		########################################################
		# Summary for the Smoothed estimator of the ROC and AUC
		########################################################
 		
		roc_cl_grilla=roc_KER_grilla
		MSE.cl <- mean((roc_real_grilla-roc_cl_grilla)*(roc_real_grilla-roc_cl_grilla))
      		DEV.cl <- mean(DEV(roc_real_grilla,roc_cl_grilla))
      		DIF2.cl  <- mean(DIF2(roc_real_grilla,roc_cl_grilla))
      		MAXDIF.cl  <-   max(abs(roc_real_grilla-roc_cl_grilla))
        	AUC<-mean(roc_cl_grilla)
		vec.auc<- c(j, MSE.cl ,DEV.cl , DIF2.cl  , MAXDIF.cl , AUC,AUC_real)
   
		lvec=length(vec.auc)
	
     		write(t(vec.auc),file=ARCHIVOKER,ncolumns=lvec,append=T)
     

 		}
			tiempo2=Sys.time()


		print(tiempo2-tiempo1)

}

  
 

 
 
		 

 

 
