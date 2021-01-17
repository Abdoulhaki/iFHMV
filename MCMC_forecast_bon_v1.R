#============================================================================
#                             G R A P H I C S    
#============================================================================

library(Rcpp)
library(mhsmm)
library(MCMCpack)
library(matrixcalc)
library(armspp)
library(actuar)
library(TruncatedNormal)
library(TruncatedDistributions)

DMS <- TRUE
if(DMS){
  sourceCpp("/home/maoudek/Rsim/Article2/MCMC_function.cpp")
  path.Data <- "/home/maoudek/Rsim/Article1/"
  path.Work <- "/home/maoudek/Rsim/Article2/Forecast/"
} else{
  sourceCpp('C:/Users/DellPC/Dropbox/MCMC_function.cpp')
  path.Work <- "C:/Users/DellPC/Dropbox/Abdoul/These/Article1/LaTeX/Graphs"
  path.Data <- "C:/Users/DellPC/Dropbox/Abdoul/These/Article1/Data"
}


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
f.p_st<-function(ech,k,S,m0,P){
  n<-length(ech)
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
  temp<-log(P[1,])+dnorm(ech[1],mean=0,sd=sqrt(V_t),log = TRUE)
  w<-c(1/(1+exp(temp[2]-temp[1])),1/(1+exp(temp[1]-temp[2])))
  prob[,1]<-w
  
  for(i in 2:n){
    V_t<-s_t[i]*c(1,m0[k])
    temp<-log(w%*%P)+dnorm(ech[i],mean=0,sd=sqrt(V_t),log = TRUE)
    w<-c(1/(1+exp(temp[2]-temp[1])),1/(1+exp(temp[1]-temp[2])))
    prob[,i]<-w
  }
  
  return(prob)
}

# Vraissemblance.
f.loglik<-function(ech,K,S,m0,r_t){
  if(K>1){
    V_t<-apply(t(m0^t(S)),1,function(x) prod(x))
  }else{
    V_t<-m0^S
  }
  
  sum(dnorm(ech,mean=0,sd=sqrt(V_t),log=TRUE))
}

vit <- 0.6
jumpmin <- 10^(-12)
jumpmax <- 10
jump_m0 <- 2
target <- 0.234
acc_m0 <- 0

f.p_K<-function(x,alpha,N)  return( alpha*sum(((1-x)^(1:N))/(1:N)) + (alpha-1)*log(x) + N*log(1-x) )

#----------------------------------------------------------------------------
#----------- G I B B S   M A I N   F O N C T I O N S   ----------------------
#----------------------------------------------------------------------------

MCMCGibbs_est<-function(ech,r_t,nstep=10000,alpha=5,seed=19812,ovars=NULL,show.progress=TRUE){
  set.seed(seed)
  n<-length(ech)
  vars<-c("n",paste("K = ",1:20,sep=""),"N+","m0","loglik_min","loglik_q1","loglik_q2",
            "loglik_mean","loglik_q3","loglik_max","AIC_min","AIC_q1","AIC_q2","AIC_mean",
            "AIC_q3","AIC_max","BIC_min","BIC_q1","BIC_q2","BIC_mean","BIC_q3","BIC_max")
  sortie<-matrix(0,1,length(vars))
  colnames(sortie)<-vars

  ####### INITIALISATION
  all_p<-list()
  all_m0<-list()
  all_K<-numeric(nstep)
  all_lik<-numeric(nstep)

  if(is.null(ovars)){
    K<-K0<-1
    S<-matrix(rbinom(n*K,size=1,prob=0.5),n,K)
    p<-c(rep(0.5,K),0)
    m0<-rep(1,K)
  }else{
    p<-as.vector(ovars$p)
    m0<-as.vector(ovars$m0)
    K<-length(m0)

    P<-matrix(c(1-p[1],p[1],p[1],1-p[1]),2,2,T)
    MC_sim<- sim.mc(c(1,0), P, n+1)-1
    S<-MC_sim[-1]

    if(K>1) for(k in 2:K){
        P<-matrix(c(1-p[k],p[k],p[k],1-p[k]),2,2,T)
        MC_sim<- sim.mc(c(1,0), P, n+1)-1
        S<-cbind(S,MC_sim[-1])
    }

    s_k<- colSums(S)
    if((sum(s_k<1)!=0)){
      p<-p[-which((s_k<1))]
      m0<-m0[-which((s_k<1))]
      S<-as.matrix(S[,-which((s_k<1))])
      K<-length(p)-1
    }
  }

  lambda0<-2.0
  eta0<-2.0
  lambda_post<-rep(lambda0,K)
  eta_post<-rep(eta0,K)

  a0<-1.0
  b0<-1.0
  a_post<-a0+K
  b_post<-b0+sum(1/(1:n))

  alpha<-alpha
  debut<-1

  for(step in 1:nstep){

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
      proba<-f.p_st(ech=ech,k=k,S=S,m0=m0,P=P)
      S[n,k]<-sample(c(0,1),1,prob=proba[,n])

      for(i in (n-1):1){
        prob<-(proba[,i]*P[,S[i+1,k]+1])/sum(proba[,i]*P[,S[i+1,k]+1])
        S[i,k]<-sample(c(0,1),1,prob=prob)
      }
    }

    s_k<- colSums(S)
    if((sum(s_k<1)!=0)){
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

      eta_post[k]<-eta0+(0.5)*sum((S[,k]*(ech^2))/temp)
    }
    for(k in 1:K){
      m0[k]<-rinvgamma(1,shape=lambda_post[k],scale=eta_post[k])
    }

    lik<-f.loglik(ech=ech,r_t=r_t,K=K,S=S,m0=m0)

    all_m0[[step]]<-m0
    all_K[step]<-K
    all_p[[step]]<-p
    all_lik[step]<-lik

    if(show.progress){
      if(step%%(nstep/500)==0) {
        print(paste("=== : ",round(100*step/nstep,2) , "%" ,"==== "))
        print(paste("K :",K))
        print(paste("p :",round(p,5)))
        print(paste("m0 :",round(m0,5)))
        print(s_k)
        print(paste("lik : ",round(lik,5)))
      }

      if(step%%(nstep/100)==0) {
        par(mfrow=c(2,1))
        print(plot(all_K[debut:(step)],type="l"))
        print(plot(all_lik[debut:(step)],type="l"))
      }

      if(step%%(nstep/2)==0) {
        debut<-step-1
      }
    }

  }

  ovars_new<-list(p=p,m0=m0)

  debut<-floor(nstep*0.75)+1

  Y <- all_m0
  auxi <- max(sapply(Y , length))
  res <- sapply(Y , function(u) c(u,rep(NA,auxi-length(u))))
  res <- apply(t(res),1,function(x) (prod(x,na.rm=TRUE))^(1/length(x[!is.na(x)])))
  sortie[1,"m0"] <- mean(res)

  res <- table(all_K[debut:nstep])
  sortie[1,colnames(sortie) %in% paste("K = ", names(res),sep="")] <- as.numeric(res[paste("K = ", names(res),sep="") %in% colnames(sortie)])
  sortie[1,"n"] <- n
  sortie[1,"N+"] <- 2*max(as.numeric(names(table(all_K[debut:nstep]))))
  sortie[1,"loglik_min"] <- quantile((all_lik[debut:nstep]))[1]
  sortie[1,"loglik_q1"] <- quantile((all_lik[debut:nstep]))[2]
  sortie[1,"loglik_q2"] <- quantile((all_lik[debut:nstep]))[3]
  sortie[1,"loglik_mean"] <- mean(all_lik[debut:nstep])
  sortie[1,"loglik_q3"] <- quantile((all_lik[debut:nstep]))[4]
  sortie[1,"loglik_max"] <- quantile((all_lik[debut:nstep]))[5]

  sortie[1,"AIC_min"] <- as.numeric(sortie[1,"loglik_min"]) - as.numeric(sortie[1,"N+"])
  sortie[1,"AIC_q1"] <- as.numeric(sortie[1,"loglik_q1"]) - as.numeric(sortie[1,"N+"])
  sortie[1,"AIC_q2"] <- as.numeric(sortie[1,"loglik_q2"]) - as.numeric(sortie[1,"N+"])
  sortie[1,"AIC_mean"] <- as.numeric(sortie[1,"loglik_mean"]) - as.numeric(sortie[1,"N+"])
  sortie[1,"AIC_q3"] <- as.numeric(sortie[1,"loglik_q3"]) - as.numeric(sortie[1,"N+"])
  sortie[1,"AIC_max"] <- as.numeric(sortie[1,"loglik_max"]) - as.numeric(sortie[1,"N+"])

  sortie[1,"BIC_min"] <- as.numeric(sortie[1,"loglik_min"]) - 0.5*as.numeric(sortie[1,"N+"])*log(n)
  sortie[1,"BIC_q1"] <- as.numeric(sortie[1,"loglik_q1"]) - 0.5*as.numeric(sortie[1,"N+"])*log(n)
  sortie[1,"BIC_q2"] <- as.numeric(sortie[1,"loglik_q2"]) - 0.5*as.numeric(sortie[1,"N+"])*log(n)
  sortie[1,"BIC_mean"] <- as.numeric(sortie[1,"loglik_mean"]) - 0.5*as.numeric(sortie[1,"N+"])*log(n)
  sortie[1,"BIC_q3"] <- as.numeric(sortie[1,"loglik_q3"]) - 0.5*as.numeric(sortie[1,"N+"])*log(n)
  sortie[1,"BIC_max"] <- as.numeric(sortie[1,"loglik_max"]) - 0.5*as.numeric(sortie[1,"N+"])*log(n)

  return(list(all_m0=all_m0,all_K=all_K,all_lik=all_lik,sortie=sortie,ovars=ovars_new))
}

MCMCGibbs_for<-function(ech,r_t,nstep=100,alpha=5,H=100,seed=19812,ovars=NULL){
  set.seed(seed)
  n<-length(ech)
  
  ####### INITIALISATION
  all_density_prev<-numeric(nstep)
  all_lik<-numeric(nstep)
  H<-H
  rt2_sim<-matrix(0,nstep,H)
  
  if(is.null(ovars)){
    K<-K0<-1
    S<-matrix(rbinom(n*K,size=1,prob=0.5),n,K)
    p<-c(rep(0.5,K),0)
    m0<-rep(1,K)
  }else{
    p<-as.vector(ovars$p)
    m0<-as.vector(ovars$m0)
    K<-length(m0)
    
    P<-matrix(c(1-p[1],p[1],p[1],1-p[1]),2,2,T)
    MC_sim<- sim.mc(c(1,0), P, n+1)-1
    S<-MC_sim[-1]
    
    if(K>1) for(k in 2:K){
      P<-matrix(c(1-p[k],p[k],p[k],1-p[k]),2,2,T)
      MC_sim<- sim.mc(c(1,0), P, n+1)-1
      S<-cbind(S,MC_sim[-1])
    }
    
    s_k<- colSums(S)
    if((sum(s_k<1)!=0)){
      p<-p[-which((s_k<1))]
      m0<-m0[-which((s_k<1))]
      S<-as.matrix(S[,-which((s_k<1))])
      K<-length(p)-1
    }
  }
  
  lambda0<-2.0
  eta0<-2.0
  lambda_post<-rep(lambda0,K)
  eta_post<-rep(eta0,K)
  
  a0<-1.0
  b0<-1.0
  a_post<-a0+K
  b_post<-b0+sum(1/(1:n))
  
  alpha<-alpha
  debut<-1
  
  for(step in 1:nstep){
    
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
      proba<-f.p_st(ech=ech,k=k,S=S,m0=m0,P=P)
      S[n,k]<-sample(c(0,1),1,prob=proba[,n])
      
      for(i in (n-1):1){
        prob<-(proba[,i]*P[,S[i+1,k]+1])/sum(proba[,i]*P[,S[i+1,k]+1])
        S[i,k]<-sample(c(0,1),1,prob=prob)
      }
    }
    
    s_k<- colSums(S)
    if((sum(s_k<1)!=0)){
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
      
      eta_post[k]<-eta0+(0.5)*sum((S[,k]*(ech^2))/temp)
    }
    for(k in 1:K){
      m0[k]<-rinvgamma(1,shape=lambda_post[k],scale=eta_post[k])
    }
    
    lik<-f.loglik(ech=ech,r_t=r_t,K=K,S=S,m0=m0)
    
    temp<-rt2(ech=ech,m0=m0,p=p,H=H,r_t=r_t)
    
    rt2_sim[step,]<-temp$rt_sim
    all_lik[step]<-lik
    all_density_prev[step]<-temp$dens_prev
  }
  
  return(list(all_lik=all_lik,all_density_prev=all_density_prev,rt2_sim=rt2_sim))
}

#----------------------------------------------------------------------------
#----------------------  D A T A    L O A D I N G  --------------------------
#----------------------------------------------------------------------------

# #-------------------------------------------------------------------------------
# #Data input
# #-------------------------------------------------------------------------------
setwd(path.Data)

data_source <- "OxfordMan" #"OxfordMan"

index_set<-c(".AEX", ".AORD", ".BFX", ".BSESN", ".BVLG", ".BVSP", ".DJI", ".FCHI", ".FTMIB", ".FTSE", ".GDAXI",
             ".GSPTSE", ".HSI", ".IBEX", ".IXIC", ".KS11", ".KSE", ".MXX", ".N225", ".NSEI", ".OMXC20", ".OMXHPI",
             ".OMXSPI", ".OSEAX", ".RUT", ".SMSI", ".SPX", ".SSEC", ".SSMI", ".STI", ".STOXX50E")

use.kernel <- TRUE #realized kernel instead of realized variance

start.date <- as.Date("2000-01-03")
end.date   <- as.Date("2020-06-06")

setwd(path.Data)
source("realized_functions.r") #load functions for realized volatility

donnees<-list()
for(i in 1:length(index_set)){
  index<-index_set[i]
  source("realized_data_traitement.r") #load and organize realized volatility data
  
  RV<-cbind(r=100*Data$RV$r)
  names(RV)<-as.Date(rownames(Data$RV))
  
  assign(index_set[i],RV)
  donnees[[i]]<-get(index_set[i])
}

index <- ".STOXX50E"
donne<-donnees[[which(index_set==index)]]
# donne<-donne[,"r"]-mean(donne[,"r"])
# n<-length(r_t)

#----------------------------------------------------------------------------
#----------------------  F O R E C A S T I N G ------------------------------
#----------------------------------------------------------------------------
setwd(path.Work)
# setwd("C:/Users/DellPC/Dropbox/Mylène/MCMC in HMM/")
  
t0<-63*12 #(9 ans)
H_range<-c(1,5,10,25,50,75,100)

date_debut<-as.Date(names(donne)[length(donne)-t0])
date_fin<-as.Date(names(donne)[length(donne)])

filename <- paste("Forecast3_iFHMV_", index,"_", date_debut, "_", date_fin, sep="")

model<-expand.grid(date = names(donne)[(length(donne)-t0):(length(donne))])
vars<-c(paste0("var_R_for",H_range),paste0("var_R_tru",H_range),'predict_loglik','N+','loglik', 'AIC', 'BIC',"m0")
model_add <- matrix(0, nrow=nrow(model), ncol=length(vars))
colnames(model_add) <- vars
model <- cbind(model, model_add)

for_R_var <- rep(0,length(H_range))
R_var <- rep(0,length(H_range))
opt<-NULL
ovars<-NULL
est_nstep<-300
for_nstep<-25
debut<-1

for(t in 0:(t0-1)){
# for(t in 718:(t0-1)){
  ech<-donne[1:(length(donne)-t0+t),]
  ech_r<-donne[(length(donne)-t0+t+1):length(donne)]
  
  update_date<-seq(0,t0,by=63)
  
  if(t %in% update_date){
    sim<-MCMCGibbs_est(ech,r_t=ech_r[1],nstep=est_nstep,alpha=5,seed=19812,ovars=ovars,show.progress=FALSE)
    ovars<-sim$ovars
    model[(t+1):(t+63),"m0"] <- round(sim$sortie[1,"m0"],5)
    model[(t+1):(t+63),"N+"] <- round(sim$sortie[1,"N+"],5)
  }
  
  sim<-MCMCGibbs_for(ech,r_t=ech_r[1],nstep=for_nstep,alpha=5,H=100,seed=19812,ovars=ovars)
    
  model[t+1,"loglik"]<- round(mean(sim$all_lik[debut:for_nstep]),5)
  model[t+1,'AIC'] <- round(model[t+1,"loglik"] - model[t+1,"N+"],5)
  model[t+1,'BIC'] <- round(model[t+1,"loglik"] - 0.5*model[t+1,"N+"]*log(n),5)
  
  model[t+1,"predict_loglik"]<-round(mean(sim$all_density_prev[debut:for_nstep]),5)
  rt2_sim<-colMeans(sim$rt2_sim[debut:for_nstep,])
  for(k in 1:length(H_range)) for_R_var[k] <- sum(rt2_sim[1:H_range[k]])
  names(for_R_var) <- paste0("var_R_for",H_range)
  model[(t+1),colnames(model) %in% names(for_R_var)]<-for_R_var
  
  for(k in 1:length(H_range)) R_var[k] <- sum((ech_r[1:H_range[k]])^2)
  names(R_var) <- paste0("var_R_tru",H_range)
  model[(t+1),colnames(model) %in% names(R_var)]<-R_var
  
  write.csv(model,file=paste0(filename,".csv"),row.names = FALSE)
  print(paste("===",round(100*(t+1)/nrow(model),2) , "%" ,"===="))
}

