#============================================================================
#                             Infinite FHMM   
#============================================================================

library(Rcpp)
library(mhsmm)
library(MCMCpack)
library(armspp)
library(actuar)
library(TruncatedDistributions)

DMS <- TRUE
if(DMS){
  path <- "/home/maoudek/Rsim/Article2/"
} 
# else {
#   path <- "C:/Users/DellPC/Dropbox/Mylène/MCMC in HMM/"
# }

setwd(path)
#----------------------------------------------------------------------------
#----------------- F O N C T I O N S    M C M C -----------------------------
#----------------------------------------------------------------------------

g<-function(vector){
  Sortie<-matrix(NA,2,2)
  for(i in c(0,1)) for(j in c(0,1)) Sortie[i+1,j+1]<-sum((vector[-1]==j)*(vector[-length(vector)]==i))
  if(vector[1]==0) Sortie[1,1]<-Sortie[1,1]+1
  if(vector[1]==1) Sortie[1,2]<-Sortie[1,2]+1
  colnames(Sortie)<-c(0,1)
  rownames(Sortie)<-c(0,1)
  return(Sortie)
}

# Probabilites lissees non stationnaires for S
f.p_st<-function(ech_r,k,S,m0,P){
  n<-length(ech_r)
  K<-ncol(S)
  prob<-matrix(NA,2,n)  # Les proba de C_t sachant x_{1:t}
  
  if(K==1) {
    s_t<-rep(1,n)
  }else if(K>2){
    s_t<-apply(t((m0[-k])^t(S[,-k])),1,function(x) prod(x))
  }else{
    s_t<-(m0[-k])^(S[,-k])
  }
  
  V_t<-s_t[1]*c(1,m0[k])
  temp<-log(P[1,])+dnorm(ech_r[1],mean=0,sd=sqrt(V_t),log = TRUE)
  w<-c(1/(1+exp(temp[2]-temp[1])),1/(1+exp(temp[1]-temp[2])))
  prob[,1]<-w
  
  for(i in 2:n){
    V_t<-s_t[i]*c(1,m0[k])
    temp<-log(w%*%P)+dnorm(ech_r[i],mean=0,sd=sqrt(V_t),log = TRUE)
    w<-c(1/(1+exp(temp[2]-temp[1])),1/(1+exp(temp[1]-temp[2])))
    prob[,i]<-w
  }
  
  return(prob)
}

# Vraissemblance.
f.loglik<-function(ech_r,K,S,m0){
  if(K>1){
    V_t<-apply(t(m0^t(S)),1,function(x) prod(x))
  }else{
    V_t<-m0^S
  }
  
  sum(dnorm(ech_r,mean=0,sd=sqrt(V_t),log=TRUE))
}

vit <- 0.6
jumpmin <- 10^(-12)
jumpmax <- 10
jump_m0 <- 2
target <- 0.234
acc_m0 <- 0

f.p_K<-function(x,alpha,N)  return( alpha*sum(((1-x)^(1:N))/(1:N)) + (alpha-1)*log(x) + N*log(1-x) )

#----------------------------------------------------------------------------
#------------------------------ D A T A  ------------------------------------
#----------------------------------------------------------------------------

####### DATA
set.seed(19812)

tru.K<-10 #5, 7, 10
nnstep<-100
sortie<-matrix(0,nnstep,min(tru.K+6,13)+5)
colnames(sortie)<-c(paste("K = ",max(1,(tru.K-6)):(tru.K+6),sep=""),"tru.K","m0","tru.m0","loglik","tru.loglik")
sstep<-1
n<-5000 #2500, 5000, 10000
tru.p<-seq(0.05,0.95,length.out = tru.K)
tru.m0<-10 #3, 10

filename<-paste0("MCMC_K",tru.K,"_m0",tru.m0,"_n",n)

for(sstep in 1:nnstep){
  
  tru.P<-matrix(c(tru.p[1],1-tru.p[1],1-tru.p[1],tru.p[1]),2,2,T)
  MC_sim<- sim.mc(c(1,0), tru.P, n+1)-1
  S<-MC_sim[-1]
  
  if(tru.K>1) {
    for(k in 2:tru.K){
      tru.P<-matrix(c(tru.p[k],1-tru.p[k],1-tru.p[k],tru.p[k]),2,2,T)
      MC_sim<- sim.mc(c(1,0), tru.P, n+1)-1
      S<-cbind(S,MC_sim[-1])
    }
    s_t<-rowSums(S)
    s_k<-colSums(S)
    V_t<-apply(t(tru.m0^t(S)),1,function(x) prod(x))
  }else{
    s_t<-S
    V_t<-tru.m0^S
  }
  
  r_t<-rnorm(n,mean = 0,sd=sqrt(V_t))
  
  tru.lik<-f.loglik(ech_r=r_t,K=tru.K,S=S,m0=tru.m0)
  
  ####### CONSTANTES
  nstep<-10000
  
  ####### INITIALISATION
  all_p<-list()
  all_m0<-list()
  all_K<-numeric(nstep)
  all_alpha<-numeric(nstep)
  all_lik<-numeric(nstep)
  K<-K0<-1
  
  S<-matrix(rbinom(n*K,size=1,prob=0.5),n,K)
  p<-c(rep(0.5,K),0)
  m0<-rep(1,K)
  
  lambda_post<-lambda0<-2.0
  eta_post<-eta0<-2.0
  
  a0<-1.0
  b0<-1.0
  a_post<-a0+K
  b_post<-b0+sum(1/(1:n))
  
  step<-1
  debut<-1
  alpha<-5
  
  for(step in 1:nstep){
    
    # alpha<-rgamma(1,a_post,b_post)
    
    if(K>1){
      nij<-g(S[,1])
      a<-nij[1,2]+nij[2,1]
      b<-nij[1,1]+nij[2,2]
      p[1]<-rtbeta(1,alpha=a,beta=b+1,a=p[2],b=1)
      if(K>2) for(k in 2:(K-1)) {
        nij<-g(S[,k])
        a<-nij[1,2]+nij[2,1]
        b<-nij[1,1]+nij[2,2]
        p[k]<-rtbeta(1,alpha=a,beta=b+1,a=p[k+1],b=p[k-1])
      }
      nij<-g(S[,K])
      a<-nij[1,2]+nij[2,1]
      b<-nij[1,1]+nij[2,2]
      p[K]<-rtbeta(1,alpha=a,beta=b+1,a=0,b=p[K-1])
      s_t<-rowSums(S)
    }else{
      nij<-g(S)
      a<-nij[1,2]+nij[2,1]
      b<-nij[1,1]+nij[2,2]
      p[1]<-rtbeta(1,alpha=a,beta=b+1,a=0,b=1)
      s_t<-S
    }
    
    p[K+1]<-arms(100,log_pdf=f.p_K,lower=0,upper=p[K],
                 metropolis=TRUE,arguments = list(a=alpha,N=n), include_n_evaluations = FALSE)[100]
    
    mu<-runif(1,0,p[K])
    
    while(p[K+1]>mu){
      p[K+2]<-arms(100,log_pdf=f.p_K,lower=0,upper=p[K+1],
                   metropolis=TRUE, arguments = list(a=alpha,N=n), include_n_evaluations = FALSE)[100]
      lambda_post[K+1]<-lambda0
      eta_post[K+1]<-eta0
      m0[K+1]<-rinvgamma(1,shape=lambda_post[K+1],scale=eta_post[K+1])
      K<-length(p)-1
    }
    S<-cbind(S,matrix(0,n,max(0,K-ncol(S))))
    
    for(k in 1:K){
      P<-matrix(c(1-p[k],p[k],p[k],1-p[k]),2,2,T)
      proba<-f.p_st(ech_r=r_t,k=k,S=S,m0=m0,P=P)
      S[n,k]<-sample(c(0,1),1,prob=proba[,n])
      
      for(i in (n-1):1){
        prob<-(proba[,i]*P[,S[i+1,k]+1])/sum(proba[,i]*P[,S[i+1,k]+1])
        S[i,k]<-sample(c(0,1),1,prob=prob)
      }
    }
    
    s_k<- colSums(S)
    if((sum(s_k<1)!=0)){
      # if((length(which(s_k<1))>0)) print(s_k)
      p<-p[-which((s_k<1))]
      lambda_post<-lambda_post[-which((s_k<1))]
      eta_post<-eta_post[-which((s_k<1))]
      m0<-m0[-which((s_k<1))]
      S<-as.matrix(S[,-which((s_k<1))])
      K<-length(p)-1
    }
    
    lambda_post<-0.5*colSums(S)+lambda0
    eta_post<-rep(eta0,K)
    for(k in 1:K){
      if(K>2){
        temp<-apply(t((m0[-k])^t(S[,-k])),1,function(x) prod(x))
      }else if(K==2){
        temp<-(m0[-k])^(S[,-k])
      }else{
        temp<-1
      }
      
      eta_post[k]<-eta0+(0.5)*sum((S[,k]*(r_t^2))/temp)
    }
    for(k in 1:K){
      m0[k]<-rinvgamma(1,shape=lambda_post[k],scale=eta_post[k])
    }
    
    lik<-f.loglik(ech_r=r_t,K=K,S=S,m0=m0)
    
    all_m0[[step]]<-m0
    all_K[step]<-K
    all_p[[step]]<-p
    all_alpha[step]<-alpha
    all_lik[step]<-lik
    
    # if(step%%20==0) {
    #   print(paste("=== : ",round(100*step/nstep,2) , "%" ,"==== "))
    #   print(paste("K :",K))
    #   print(paste("p :",round(p,5)))
    #   print(paste("m0 :",round(m0,5)))
    #   print(s_k)
    #   print(paste("lik : ",round(lik,5), ", TRUE : ", round(tru.lik,5)))
    # }
    # 
    # if(step%%100==0) {
    #   par(mfrow=c(2,1))
    #   print(plot(all_K[debut:(step)],type="l"))
    #   print(lines(rep(tru.K,step-debut+1),col="red"))
    #   temp<-all_lik[debut:(step)]
    #   temp[(temp>tru.lik+2500)]<-tru.lik+2500
    #   temp[(temp<tru.lik-2500)]<-tru.lik-2500
    #   print(plot(temp,type="l"),ylim=c(tru.lik-2500,tru.lik+2500))
    #   print(lines(rep(tru.lik,step-debut+1),col="red"))
    # }
    
  }
  
  # par(mfrow=c(2,1))
  # print(plot(all_K[1:(nstep)],type="l"))
  # print(lines(rep(tru.K,nstep),col="red"))
  # temp<-all_lik[1:(nstep)]
  # temp[(temp>tru.lik+2500)]<-tru.lik+2500
  # temp[(temp<tru.lik-2500)]<-tru.lik-2500
  # print(plot(temp,type="l"),ylim=c(tru.lik-2500,tru.lik+2500))
  # print(lines(rep(tru.lik,nstep),col="red"))
  
  if(nstep!=step) nstep<-step-1
  debut<-floor(nstep*0.75)+1
  
  # Y<-all_p[debut:nstep]
  # auxi <- max(sapply(Y , length))
  # res <- sapply(Y , function(u) c(u,rep(NA,auxi-length(u))))
  # res <- t(res)
  # round(colMeans(res),5)
  # round(1-tru.p,5)
  # for(k in 1:K){
  #   print(plot(res[,k],type='l'))
  #   print(lines(rep(1-tru.p[k],length(res[,k])),type="l",col="red"))
  # }
  #
  # plot(all_sigma[debut:nstep],type='l')
  # lines(rep(tru.sigma,(step-1)),col="red")
  # 
  # plot(all_m0[debut:nstep],type='l')
  # lines(rep(tru.m0,(step-1)),col="red")
  
  Y <- all_m0[debut:nstep]
  auxi <- max(sapply(Y , length))
  res <- sapply(Y , function(u) c(u,rep(NA,auxi-length(u))))
  res <- apply(t(res),1,function(x) (prod(x,na.rm=TRUE))^(1/length(x[!is.na(x)])))
  sortie[sstep,"m0"] <- mean(res)
  
  res <- table(all_K[debut:nstep])
  sortie[sstep,colnames(sortie) %in% paste("K = ", names(res),sep="")] <- as.numeric(res[paste("K = ", names(res),sep="") %in% colnames(sortie)])
  sortie[sstep,"tru.m0"] <- tru.m0
  sortie[sstep,"tru.K"] <- tru.K
  sortie[sstep,"loglik"] <- mean(all_lik[debut:nstep])
  sortie[sstep,"tru.loglik"] <- tru.lik
  
  write.csv(sortie,file=paste0(filename,".csv"))
  
}