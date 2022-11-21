
#################################################
# Quantile Computation 
#################################################

weighted.fractile<-function (y, w, p){
    w <- w/sum(w)
    a <- 1 - p
    b <- p
    ox <- order(y)
    y <- y[ox]
    w <- w[ox]
    k <- 1
    low <- cumsum(c(0, w))
    up <- sum(w) - low
    df <- a * low - b * up
    repeat {
        if (df[k] < 0) 
            k <- k + 1
        else if (df[k] == 0) 
            return((w[k] * y[k] + w[k - 1] * y[k - 1])/(w[k] + 
                w[k - 1]))
        else return(y[k - 1])
    }
}




#################################################
# IPW Estimator of the Distribution Function
#################################################

FY.IPW<-function(y,delta,pyx,puntos){

	nsamp=length(y)

 	tau1=rep(0,nsamp)


	for (i in 1: nsamp){
    
     		tau1[i]=  delta[i]/pyx[i]   
      
     	}

	tau= tau1/sum(tau1)
 
 

	lp=length(puntos)

	FY.IPW=rep(NA,lp)

	for (ele in 1:lp){
		FY.IPW[ele]=sum(tau*1*(y<=puntos[ele]),na.rm=TRUE)

	}
 	return(FY.IPW)

}

#################################################
# ROC(p) IPW Estimator
#################################################

roc_IPW<- function(y_D,y_H,delta_D,delta_H,pyx_D,pyx_H, p) 
{
  	n_H=length(y_H)

 	tau1=rep(0,n_H)


	for (i in 1: n_H){
    
     		tau1[i]=  delta_H[i]/pyx_H[i]   
      
     	}

	tau= tau1/sum(tau1)
 
	punto<- unname(weighted.fractile(y_H,tau,1-p))
	roc_IPW <- 1-FY.IPW(y=y_D,delta=delta_D,pyx=pyx_D,puntos=punto)
  	
	return(roc_IPW)
}


#######################################################
# Convolution Estimator of the Distribution Function
#######################################################

FY.MUL<-function(y,delta,pyx,m.est,puntos){

	nsamp=length(y)

 

	tau= rep(1,nsamp)/nsamp
	kapa= delta/sum(delta)

 
	
	pes=ysomb= matrix(rep(0,nsamp*nsamp),nsamp,nsamp)


	for (i in 1: nsamp){
		for (j in 1: nsamp){
			pes[i,j]<-kapa[i]*tau[j]
			ysomb[i,j]<-m.est[j]+(y[i]-m.est[i])
 
		}}

	lp=length(puntos)

	FY.mul=rep(NA,lp)

	for (ele in 1:lp){
		FY.mul[ele]=sum(pes*1*(ysomb<=puntos[ele]),na.rm=TRUE)

	}
 	return(FY.mul)

}

#################################################
# ROC(p) Convolution Based Estimator
#################################################

roc_MUL<- function(y_D,y_H,delta_D,delta_H,pyx_D,pyx_H, m.est_D,m.est_H, p) 
{
  	n_H=length(y_H)

 	tau= rep(1,n_H)/n_H
	kapa= delta_H/sum(delta_H)
 
	pes=ysomb=matrix(rep(0,n_H*n_H),n_H,n_H)


	for (i in 1: n_H){
		for (j in 1: n_H){
			pes[i,j]<-kapa[i]*tau[j]
			ysomb[i,j]<-m.est_H[j]+(y_H[i]-m.est_H[i])
		}}

	pes.vec=as.vector(pes)
	ysomb.vec=as.vector(ysomb)
 

	ysomb.vec[pes.vec==0]=NA
 
 
	punto<- unname(weighted.fractile(y=ysomb.vec, w=pes.vec,p=1-p))
	roc_MUL <- 1-FY.MUL(y=y_D,delta=delta_D,pyx=pyx_D,m.est=m.est_D,puntos=punto)
  	
	return(roc_MUL)
}

#################################################
# Integrated Kernel
#################################################

nucleo=function(x){
  a=0.75*(x-x^3/3+2/3)
  b=a*(abs(x)<1)+1*(1<=x)
  return(b)
}


#################################################
# Kernel Estimator
#################################################

FY.KER<-function(y,delta,pyx,ache,puntos){
  
  nsamp=length(y)
  
  tau1=rep(0,nsamp)
  
  
  for (i in 1: nsamp){
    
    tau1[i]=  delta[i]/pyx[i]   
    
  }
  
  tau= tau1/sum(tau1)
  
  arg=(puntos-y)/ache
  FY.KER =sum(tau*nucleo(arg),na.rm=TRUE)
  
  
  return(FY.KER)
  
}


#################################################
# ROC(p)   Smoothed Estimator 
#################################################

roc_KER<- function(y_D,y_H,delta_D,delta_H,pyx_D,pyx_H, p) 
{
  n_H=length(y_H)
  n_D=length(y_D)
  
  
  cn=1+1.8*n_D^(-1/5)
  ache= cn*sqrt(5*p*(1-p))/sqrt(2*n_H)
  
#################################################
# Pseudo Data
#################################################
  Z=rep(NA,n_D)
  
  for (i in 1:n_D){
    punto=y_D[i]
    Z[i]=1-FY.IPW(y=y_H,delta=delta_H,pyx=pyx_H,puntos= punto) 
  }
  
  Z[delta_D==0]<- NA
  
  tau1=rep(0,n_H)
  
  
  for (i in 1: n_H){
    
    tau1[i]=  delta_H[i]/pyx_H[i]   
    
  }
  
  tau= tau1/sum(tau1)
  
  
  roc_KER <- FY.KER(y=Z,delta=delta_D,pyx=pyx_D,ache=ache,puntos=p)
  
  return(roc_KER)
}


