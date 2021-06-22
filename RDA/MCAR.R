library(maps)
#Import California map
ca.county = map("county","california", fill=TRUE, plot=FALSE)
library(readr)
library(spdep)
library(maptools)
library(classInt)
library(RColorBrewer)
library(tidyr)
library(MASS)
library(Matrix)
library(mapproj)
library(rjags)
library(R2jags)
library(lattice)
library(mvtnorm)
library(matrixStats)
library(fields)
library(boot)
library(blockmatrix)
library(ggmcmc)
library(mcmc)
library(magic)
library(msos)
library(AICcmodavg)
library(coda)
library(invgamma)
library(mcmcse)
library(LaplacesDemon)
library(gtools)
library(ggpubr)

#Import covariates
rate_5y <- read.csv("age_adjusted.csv")
covariates <- read.csv("covariates.csv")
race <- read.csv("race.csv")
sex <- read.csv("sex.csv")
insurance <- read.csv("insurance.csv")
smoking <- read.csv("smoking.csv")
smoking$smoking <- as.numeric(substr(smoking$Cigarette.Smoking.Rate., 1,4))

#Import incidence data for 4 cancers in California
rate_CA = rate_5y[substr(rate_5y$State_county,1,2) == "CA",]

rate_lung = rate_CA[rate_CA$Site_recode_ICD_O_3_WHO_2008=="Lung and Bronchus",]
rate_lung = rate_lung[order(readr::parse_number(as.character(rate_lung$State_county))),]
rate_lung$E_count = (sum(rate_lung$Count) / sum(rate_lung$Population)) * rate_lung$Population
rate_lung$standard_ratio = rate_lung$Count / rate_lung$E_count

rate_esophagus = rate_CA[rate_CA$Site_recode_ICD_O_3_WHO_2008=="Esophagus",]
rate_esophagus = rate_esophagus[order(readr::parse_number(as.character(rate_esophagus$State_county))),]
rate_esophagus$E_count = (sum(rate_esophagus$Count) / sum(rate_esophagus$Population)) * rate_esophagus$Population
rate_esophagus$standard_ratio = rate_esophagus$Count / rate_esophagus$E_count

rate_larynx = rate_CA[rate_CA$Site_recode_ICD_O_3_WHO_2008=="Larynx",]
rate_larynx = rate_larynx[order(readr::parse_number(as.character(rate_larynx$State_county))),]
rate_larynx$E_count = (sum(rate_larynx$Count) / sum(rate_larynx$Population)) * rate_larynx$Population
rate_larynx$standard_ratio = rate_larynx$Count / rate_larynx$E_count

rate_colrect = rate_CA[rate_CA$Site_recode_ICD_O_3_WHO_2008=="Colon and Rectum",]
rate_colrect = rate_colrect[order(readr::parse_number(as.character(rate_colrect$State_county))),]
rate_colrect$E_count = (sum(rate_colrect$Count) / sum(rate_colrect$Population)) * rate_colrect$Population
rate_colrect$standard_ratio = rate_colrect$Count / rate_colrect$E_count

county.ID <- sapply(strsplit(ca.county$names, ","), function(x) x[2])
ca.poly = map2SpatialPolygons(ca.county, IDs=county.ID)
ca.poly$rate_lung = rate_lung$standard_ratio
ca.poly$rate_esophagus = rate_esophagus$standard_ratio
ca.poly$rate_larynx = rate_larynx$standard_ratio
ca.poly$rate_colrect = rate_colrect$standard_ratio
ca.poly$smoking = smoking$smoking

ca.coords = coordinates(ca.poly)
################## Data construction and adjacency matrix ###########

## Data

county_attribute = covariates[substr(covariates$State_county,1,2) == "CA",]
county_attribute$state = extract_numeric(county_attribute$State_county)
county_attribute1 = data.frame(county_attribute[order(county_attribute$state),])
county_attribute1$V_Persons_age_18_ACS_2012_2016 = as.numeric(county_attribute1$V_Persons_age_18_ACS_2012_2016)/100
county_attribute1$V_Persons_age_65_ACS_2012_2016 = as.numeric(county_attribute1$V_Persons_age_65_ACS_2012_2016)/100
county_attribute1$VHighschooleducationACS2012201 = as.numeric(county_attribute1$VHighschooleducationACS2012201)/100
county_attribute1$VFamiliesbelowpovertyACS201220 = as.numeric(county_attribute1$VFamiliesbelowpovertyACS201220)/100
county_attribute1$V_Unemployed_ACS_2012_2016 = as.numeric(county_attribute1$V_Unemployed_ACS_2012_2016)/100

race1 = race[substr(race$State_county,1,2) == "CA"&race$Race_recode_White_Black_Other=="Black",]
sex1 = sex[substr(sex$State_county,1,2) == "CA"&sex$Sex=="Male",]
insurance1 = insurance[substr(insurance$State_county,1,2) == "CA"&insurance$Insurance_Recode_2007=="Uninsured",]

rate_lung1 = cbind(rate_lung, smoking$smoking, county_attribute1[,2:6], race1$Row_Percent, sex1$Row_Percent,insurance1$Row_Percent)
colnames(rate_lung1) = c("county", "site", "rate", "count", "population", "E_count", "standard_ratio", "smoking", "young","old", "highschool", "poverty", "unemployed", "black", "male", "uninsured")
rate_esophagus1 = cbind(rate_esophagus, smoking$smoking, county_attribute1[,2:6], race1$Row_Percent, sex1$Row_Percent,insurance1$Row_Percent)
colnames(rate_esophagus1) = c("county", "site", "rate", "count", "population", "E_count", "standard_ratio", "smoking", "young","old", "highschool", "poverty", "unemployed", "black", "male", "uninsured")
rate_larynx1 = cbind(rate_larynx, smoking$smoking, county_attribute1[,2:6], race1$Row_Percent, sex1$Row_Percent,insurance1$Row_Percent)
colnames(rate_larynx1) = c("county", "site", "rate", "count", "population", "E_count", "standard_ratio", "smoking", "young","old", "highschool", "poverty", "unemployed", "black", "male", "uninsured")
rate_colrect1 = cbind(rate_colrect, smoking$smoking, county_attribute1[,2:6], race1$Row_Percent, sex1$Row_Percent,insurance1$Row_Percent)
colnames(rate_colrect1) = c("county", "site", "rate", "count", "population", "E_count", "standard_ratio", "smoking", "young","old", "highschool", "poverty", "unemployed", "black", "male", "uninsured")

## Adjacency matrix and neighbor info
ca.neighbors = poly2nb(ca.poly)
n=length(ca.neighbors)

Adj=sapply(ca.neighbors,function(x,n) {v=rep(0,n);v[x]=1;v},n)
colnames(Adj)=county.ID
ca.coord = coordinates(ca.poly)
ca.latrange=round(quantile(ca.coord[,2],c(0.25,0.75)))
ca.albersproj=mapproject(ca.coord[,1],ca.coord[,2],projection = "albers",param=ca.latrange)

# Neighboring counties information
neighborvec0 = NULL
neighbor_list0 = NULL
neighbor_name = NULL
for(i in 1:(n-1)){
  for(j in (i+1):n){
    if(Adj[i,j] == 1){
      neighborvec0 = c(neighborvec0, paste(i, ",", j, sep=""))
      neighbor_list0 = rbind(neighbor_list0, c(i,j))
      neighbor_name = c(neighbor_name, paste(colnames(Adj)[i], ", ", colnames(Adj)[j], sep=""))
    }
  }
}

## Specify the order and reorder the map with new neighbor info
perm=order(ca.albersproj$x-ca.albersproj$y)
Adj_new=Adj[perm,perm]

n=nrow(Adj_new)
ni=rowSums(Adj_new)
maxn=max(ni)
neimat=matrix(0,n,maxn)
neighbors=lapply(1:n,function(x) which(Adj_new[x,]==1))
#N(i): 2:n
dneighbors=sapply(2:n,function(i) intersect(neighbors[[i]],1:(i-1)))
#n<i: 2:n
dni=sapply(dneighbors,length)
original_perm = 1:58
index2=c(1,which(dni==0)+1)

final_perm=c(original_perm[perm][index2],
             original_perm[perm][-index2])
final_perm[order(final_perm)]

Minc = Adj[final_perm,final_perm]
n=nrow(Minc)
ni=rowSums(Minc)
maxn=max(ni)
neimat=matrix(0,n,maxn)
neighbors=lapply(1:n,function(x) which(Minc[x,]==1))
#N(i): 2:n
dneighbors=sapply(2:n,function(i) intersect(neighbors[[i]],1:(i-1)))
#n<i: 2:n

dni=sapply(dneighbors,length)
nmax=max(dni)
cni=cumsum(dni)
dneimat=sapply(dneighbors, function(nei,nmax,n) c(nei,rep(n+1,nmax+1-length(nei))),nmax,n)
udnei=unlist(dneighbors)

ni_wo = sapply(neighbors,length)
cni_wo = cumsum(ni_wo)
udnei_wo = unlist(neighbors)
cn = c(0, cni)
ns = dni

#DAGAR
region = seq(1:n)
cn = c(0, cni)
ns = dni
index = list()
for(i in 1:(n-2)){
  index[[i]] = region[-(udnei[(cn[i+1] + 1):(cn[i+1] + ns[i+1])])]
}
index1 = unlist(index)

#CAR
D=matrix(0,n,n)
for(i in 1:n) D[i,i]=ni[i]

### Functions for MCMC updates

# compute stick-breaking weights
makeprobs<-function(v){
  m<-length(v) 
  p <-matrix(0,m,m)
  for(j in 2:m){
    p[1:(j-1),j]<-1
  }
  probs<-exp(log(v)+log(1-v)%*%p)
  probs[m]<-1-sum(probs[1:(m-1)])
  probs
}

# compute u weights from F_r and prob 
makeu<- function(F_r,probs){
  
  m1=length(F_r)
  m2=length(probs)
  u=rep(0,m1)
  for (k in 1:m1){
    for (l in 1:m2){
      if (sum(probs[1:(l-1)])<F_r[k] && F_r[k]<sum(probs[1:l])){
        u[k]=l
      }
      if (u[k]==0){
        u[k]=1
      }
    }
  }
  return(u)
}

# Compute presicion matrices of DAGAR 
Dinv_new <- function(Rho, n, cn, ns, udnei,q){
  Tau <- list()
  B <- list()
  invD <- list()
  for(i in 1:q){
    Tau[[i]] <- diag(n)
    B[[i]] <- matrix(0, n, n)
    for(j in 3:n){
      Tau[[i]][j,j] <- (1 + (ns[j-1] - 1) * Rho[i]^2) / (1 - Rho[i]^2)
      for(k in (udnei[(cn[j-1] + 1):(cn[j-1] + ns[j-1])])){
        B[[i]][j,k] <- Rho[i] / (1 + (ns[j-1] - 1) * Rho[i]^2)
      }
    }
    invD[[i]] <- t(diag(n) - B[[i]]) %*% Tau[[i]] %*% (diag(n) - B[[i]])
  }
  return(invD)
}

# Jacobian matrix for updating A
Jacob_A <- function(A){
  prod = 1
  for(i in 1:nrow(A)){
    prod = prod*A[i,i]^(nrow(A)-i+1)
  }
  return(2^nrow(A)*prod)
}

# Data in model
Y1 = rate_lung1$count[final_perm]
Y2 = rate_esophagus1$count[final_perm]
Y3 = rate_larynx1$count[final_perm]
Y4 = rate_colrect1$count[final_perm]

E1 = rate_lung1$E_count[final_perm]
E2 = rate_esophagus1$E_count[final_perm]
E3 = rate_larynx1$E_count[final_perm]
E4 = rate_colrect1$E_count[final_perm]

X1 = as.matrix(cbind(1,rate_lung1[,8:16]))[final_perm,]
X2 = as.matrix(cbind(1,rate_esophagus1[,8:16]))[final_perm,]
X3 = as.matrix(cbind(1,rate_larynx1[,8:16]))[final_perm,]
X4 = as.matrix(cbind(1,rate_colrect1[,8:16]))[final_perm,]

Y = c(Y1,Y2,Y3,Y4)
E = c(E1, E2, E3, E4)
X = as.matrix(bdiag(bdiag(X1[,1], X2[,1]), bdiag(X3[,1],X4[,1])))

### MCAR model
# Metropolis within Gibbs Sampler for MCMC updating
ARDP_joint_diff_car<-function(y=Y, x1=as.matrix(X1[,1]), x2=as.matrix(X2[,1]), x3 = as.matrix(X3[,1]), x4 = as.matrix(X4[,1]), 
                              X = X, E = E, Minc, alpha=1, q = 4, n.atoms=15, runs=10000, burn=1000){
  #y:       data
  #x:       covariates
  #n.atoms: number of atoms in the mixture dist.
  #theta:      the theta's (iid from baseline) in the model
  #alpha:     v~beta(1,alpha)
  #u:       the index indicator of spatial random effect 
  #rho:     sptatial autocorrelation parameter in DAGAR
  #Minc:       0-1 adjacency matrix 
  #V_r:     covariance matrix of joint MDAGAR
  #Q:       presicion matrix of DAGAR
  #r:       random effects following DAGAR
  #F_r:     Marginal CDF of r
  #taued:   presicion in Baseline of DP for disease d   
  #taus:    precision for theta
  
  
  
  
  nq<-length(y)
  n <- nq / q
  p<-ncol(X)
  y1 <- y[1:n]
  y2 <- y[(n+1):(2*n)]
  y3 <- y[(2*n+1):(3*n)]
  y4 <- y[(3*n+1):(4*n)]
  
  E1 <- E[1:n]
  E2 <- E[(n+1):(2*n)]
  E3 <- E[(2*n+1):(3*n)]
  E4 <- E[(3*n+1):(4*n)]
  
  sigmasq_beta = 10000
  keepbeta=matrix(0,runs,p)
  keepphi=matrix(0,runs,nq)
  keeptheta=matrix(0,runs,n.atoms)
  keepA=matrix(0,runs,10)
  keepu=matrix(0,runs,nq)
  keeprho1=keeprho2=keeprho3=keeprho4=keeptaus = rep(0,runs)
  keepv=matrix(0,runs,n.atoms)
  keepr=matrix(0,runs,nq)
  keepFr=matrix(0,runs,nq)
  
  
  #initial values
  
  theta=rep(0,n.atoms)
  
  beta1<-rep(0,ncol(x1))
  beta2<-rep(0,ncol(x2))
  beta3<-rep(0,ncol(x3))
  beta4<-rep(0,ncol(x4))
  beta = c(beta1, beta2, beta3, beta4)
  
  taus = 1
  
  c=2
  d=0.1
  d2=1
  v<-rep(.1,n.atoms)
  probs=makeprobs(v)
  
  #A matrix
  rho=c(0.5, 0.5, 0.5, 0.5)
  eta = log(rho/(1-rho))
  A = matrix(0, q, q)
  for(i in 1:q){
    for(j in 1:i){
      if(j == i){
        A[i, j] = exp(rnorm(1))
      }else{
        A[i, j] = rnorm(1)
      }
    }
  }
  
  invQ1 = solve(D - rho[1] * Minc)
  invQ2 = solve(D - rho[2] * Minc)
  invQ3 = solve(D - rho[3] * Minc)
  invQ4 = solve(D - rho[4] * Minc)
  
  kprod = kronecker(A, diag(n))
  invQ = as.matrix(bdiag(bdiag(invQ1, invQ2), bdiag(invQ3, invQ4)))
  
  Vr = as.matrix(forceSymmetric(kprod %*% invQ %*% t(kprod)))
  r = rmvnorm(1, rep(0, nq), Vr)
  
  
  F_r=pnorm(r,0,sqrt(diag(Vr)))
  u=makeu(F_r,probs)
  phi <- theta[u]
  
  nu = 2
  R = 0.1 * diag(q)
  
  
  acceptr=acceptv=acceptrho=acceptA=0
  acceptbeta1 = acceptbeta2=acceptbeta3=acceptbeta4 = 0
  accepttheta = 0
  
  count<-afterburn<-0
  burn = burn + 1
  
  for(iter in 1:runs){
    
    if(iter %% 100 == 0){
      print(iter)
      print(acceptA/(iter-1)) 
      print(acceptr/nq/(iter-1))
      print(acceptrho/(iter-1)) 
      print(acceptv/n.atoms/(iter-1))
      print(accepttheta/n.atoms/(iter-1))
      print(acceptbeta1/(iter-1))
      print(acceptbeta2/(iter-1))
      print(acceptbeta3/(iter-1))
      print(acceptbeta4/(iter-1))
    }
    ######################
    ###   update beta  ###
    ###################### 
    # update beta (intercept only model)
    
    pro_beta1=rnorm(1,beta1,0.02)
    MHrate1=sum(-E1*exp(pro_beta1*x1+theta[u][1:n])+y1*(pro_beta1*x1+theta[u][1:n]))-
      sum(-E1*exp(beta1*x1+theta[u][1:n])+y1*(beta1*x1+theta[u][1:n]))  
    
    if(runif(1,0,1)<exp(MHrate1)){
      beta1<-pro_beta1
      acceptbeta1=acceptbeta1+1
    } 
    
    pro_beta2=rnorm(1,beta2,0.05)
    MHrate2=sum(-E2*exp(pro_beta2*x2+theta[u][(n+1):(2*n)])+y2*(pro_beta2*x2+theta[u][(n+1):(2*n)]))-
      sum(-E2*exp(beta2*x2+theta[u][(n+1):(2*n)])+y2*(beta2*x2+theta[u][(n+1):(2*n)]))  
    
    if(runif(1,0,1)<exp(MHrate2)){
      beta2<-pro_beta2
      acceptbeta2=acceptbeta2+1
    } 
    
    pro_beta3=rnorm(1,beta3,0.05)
    MHrate3=sum(-E3*exp(pro_beta3*x3+theta[u][(2*n+1):(3*n)])+y3*(pro_beta3*x3+theta[u][(2*n+1):(3*n)]))-
      sum(-E3*exp(beta3*x3+theta[u][(2*n+1):(3*n)])+y3*(beta3*x3+theta[u][(2*n+1):(3*n)]))  
    
    if(runif(1,0,1)<exp(MHrate3)){
      beta3<-pro_beta3
      acceptbeta3=acceptbeta3+1
    } 
    
    pro_beta4=rnorm(1,beta4,0.02)
    MHrate4=sum(-E4*exp(pro_beta4*x4+theta[u][(3*n+1):(4*n)])+y4*(pro_beta4*x4+theta[u][(3*n+1):(4*n)]))-
      sum(-E4*exp(beta4*x4+theta[u][(3*n+1):(4*n)])+y4*(beta4*x4+theta[u][(3*n+1):(4*n)]))  
    
    if(runif(1,0,1)<exp(MHrate4)){
      beta4<-pro_beta4
      acceptbeta4=acceptbeta4+1
    } 
    
    beta <- c(beta1, beta2, beta3, beta4)
    
    
    #########################
    ###   update theta    ###
    ######################### 
    
    u1 = u[1:n]
    u2 = u[(n+1):(2*n)]
    u3 = u[(2*n+1):(3*n)]
    u4 = u[(3*n+1):(4*n)]
    
    
    
    for (j in 1:n.atoms){
      
      count[j]=sum(u==j)
      
      if (count[j]==0){
        theta[j]=rnorm(1,0,sqrt(1/taus))
      }else{
        pro_theta=rnorm(1,theta[j],0.06)
        
        MHrate<-sum(-(E[u==j])*exp(X[u==j,] %*% beta+pro_theta)+(y[u==j])*(X[u==j,] %*% beta+pro_theta))-taus/2*pro_theta^2-
          sum(-(E[u==j])*exp(X[u==j,] %*% beta+theta[j])+(y[u==j])*(X[u==j,] %*% beta+theta[j]))+taus/2*theta[j]^2
        if(runif(1,0,1)<exp(MHrate)){
          theta[j]<-pro_theta
          accepttheta=accepttheta+1
        } 
      }
    }
    
    
    
    ######################
    ###   update r     ###
    ###################### 
    
    for (k in 1:nq){
      pro_r=r;pro_Fr=F_r;pro_u=u
      pro_r[k]=rnorm(1,r[k],1.7)
      pro_Fr[k]=pnorm(pro_r[k],0,sqrt(Vr[k,k]))			
      pro_u[k]=makeu(pro_Fr[k],probs)
      
      MH=dmvnorm(pro_r,mean=rep(0, nq),sigma=Vr,log=T)+dpois(y[k],E[k]*exp(as.numeric(X[k,]%*%beta)+theta[pro_u[k]]),log=T)-
        dmvnorm(r,mean=rep(0, nq),sigma=Vr,log=T)-dpois(y[k],E[k]*exp(as.numeric(X[k,]%*%beta)+theta[u[k]]),log=T)
      if(runif(1,0,1)<exp(MH)){
        r[k]=pro_r[k]
        F_r[k]=pro_Fr[k]
        u[k]=pro_u[k]
        acceptr=acceptr+1
      }
    }
    
    
    ######################
    ###   update   v   ###
    ###################### 
    
    for (j in 1:n.atoms){

      pro_v=v
      pro_v[j]<-rnorm(1,v[j],0.06)
      while (pro_v[j]<0 || pro_v[j]>1) {pro_v[j]<-rnorm(1,v[j],0.06)}
  
      pro_probs=makeprobs(pro_v)
      pro_u=makeu(F_r,pro_probs)
      
      MH=log(dbeta(pro_v[j],1,alpha))+sum(dpois(y,E*exp(X%*%beta+theta[pro_u]),log=T))-
        log(dbeta(v[j],1,alpha))-sum(dpois(y,E*exp(X%*%beta+theta[u]),log=T))
      
      
      if(runif(1,0,1)<exp(MH)){
        v[j]=pro_v[j]
        probs=pro_probs
        u=pro_u
        acceptv=acceptv+1	
      }
      # }
    }
    
    
    ######################
    ###   update taus  ###
    ###################### 
    
    
    taus=rgamma(1,shape=1/2*n.atoms+c,rate=sum(theta^2)/2+d2)
    
    
    ######################
    ###   update rho  ###
    ###################### 
    
    
    pro_eta1 = rnorm(1, eta[1], 1.1)
    pro_eta2 = rnorm(1, eta[2], 1.1)
    pro_eta3 = rnorm(1, eta[3], 1.1)
    pro_eta4 = rnorm(1, eta[4], 1.1)
    pro_eta = c(pro_eta1, pro_eta2, pro_eta3, pro_eta4)
    pro_rho = exp(pro_eta)/(1+exp(pro_eta))
    
    pro_invQ1 = solve(D - pro_rho[1] * Minc)
    pro_invQ2 = solve(D - pro_rho[2] * Minc)
    pro_invQ3 = solve(D - pro_rho[3] * Minc)
    pro_invQ4 = solve(D - pro_rho[4] * Minc)
    
    kprod = kronecker(A, diag(n))
    pro_invQ = as.matrix(bdiag(bdiag(pro_invQ1, pro_invQ2), bdiag(pro_invQ3, pro_invQ4)))
    pro_Vr = kprod %*% pro_invQ %*% t(kprod)
    
    MH = dmvnorm(r,mean=rep(0, nq),sigma=as.matrix(forceSymmetric(pro_Vr)),log=T) + 
      log(pro_rho[1]) + log(1-pro_rho[1]) + log(pro_rho[2]) + log(1-pro_rho[2]) +
      log(pro_rho[3]) + log(1-pro_rho[3]) + log(pro_rho[4]) + log(1-pro_rho[4]) - 
      dmvnorm(r,mean=rep(0, nq),sigma=Vr,log=T) - log(rho[1]) - log(1-rho[1]) - log(rho[2]) - log(1-rho[2]) -
      log(rho[3]) - log(1-rho[3]) - log(rho[4]) - log(1-rho[4])
    
    if(runif(1,0,1)<exp(MH)){
      eta = pro_eta
      rho = pro_rho
      Vr = as.matrix(forceSymmetric(pro_Vr))
      acceptrho=acceptrho+1	
    }
    
    ######################
    ###   update A  ###
    ###################### 

    pro_A = matrix(0, q, q)
    for(i in 1:q){
      for(j in 1:i){
        if(j == i){
          pro_A[i, j] = exp(rnorm(1, log(A[i, j]), 0.05))
        }else{
          pro_A[i, j] = rnorm(1, A[i, j], 0.05)
        }
      }
    }
    
    invQ1 = solve(D - rho[1] * Minc)
    invQ2 = solve(D - rho[2] * Minc)
    invQ3 = solve(D - rho[3] * Minc)
    invQ4 = solve(D - rho[4] * Minc)
    
    Sigma =  A %*% t(A)
    
    pro_kprod = kronecker(pro_A, diag(n))
    invQ = as.matrix(bdiag(bdiag(invQ1, invQ2), bdiag(invQ3, invQ4)))
    pro_Vr = pro_kprod %*% as.matrix(forceSymmetric(invQ)) %*% t(pro_kprod)
    pro_Sigma =  pro_A %*% t(pro_A)
    
    lpA = -(nu+4)/2*logdet(Sigma) - 1/2*sum(diag(nu*R%*%solve(Sigma))) + log(Jacob_A(A)) + sum(log(diag(A)))
    pro_lpA = -(nu+4)/2*logdet(pro_Sigma) - 1/2*sum(diag(nu*R%*%solve(pro_Sigma))) + log(Jacob_A(pro_A)) + sum(log(diag(pro_A)))
    MH = dmvnorm(r,mean=rep(0, nq),sigma=as.matrix(forceSymmetric(pro_Vr)),log=T) + pro_lpA - 
      dmvnorm(r,mean=rep(0, nq),sigma=Vr,log=T) - lpA
    
    if(runif(1,0,1)<exp(MH)){
      A = pro_A
      Vr = as.matrix(forceSymmetric(pro_Vr))
      acceptA=acceptA+1	
    }
    
    
    ######################
    ### record samples ###
    ###################### 
    
    keeptheta[iter,] = theta
    keepphi[iter,] = theta[u]
    
    keepA[iter,] = as.vector(A)[-c(5,9,10,13:15)]
    
    keepbeta[iter,]=beta
    keeptaus[iter]=taus
    keeprho1[iter]=rho[1]
    keeprho2[iter]=rho[2]
    keeprho3[iter]=rho[3]
    keeprho4[iter]=rho[4]
    keepu[iter,]=u
    keepv[iter,]=v
    keepr[iter,]=r
    keepFr[iter,]=F_r
    
    #	cat("iteration = ", i, "acceptance rate of r = ", acceptr/(n*i),"acceptance rate of v = ",acceptv/(n.atoms*i), "\n")
    
  }
  
  list(beta=keepbeta[burn:runs,], A = keepA[burn:runs,],phi=keepphi[burn:runs,],theta=keeptheta[burn:runs,],u=keepu[burn:runs,],
       v=keepv[burn:runs,],r=keepr[burn:runs,],Fr=keepFr[burn:runs,], taus=keeptaus[burn:runs],
       rho1=keeprho1[burn:runs], rho2=keeprho2[burn:runs], rho3=keeprho3[burn:runs], rho4=keeprho4[burn:runs])
}
set.seed(123)
mcmc_samples=ARDP_joint_diff_car(y=Y, x1=as.matrix(X1[,1]), x2=as.matrix(X2[,1]), x3 = as.matrix(X3[,1]), x4 = as.matrix(X4[,1]), 
                                 X = X, E=E, Minc, alpha=1, q = 4, n.atoms=15, runs=30000, burn=20000)

sample.mcmc = cbind(mcmc_samples$taus, mcmc_samples$rho1, mcmc_samples$rho2, mcmc_samples$rho3, mcmc_samples$rho4)

phis <- mcmc_samples$phi
phis1 <- phis[,1:58][,order(final_perm)]
phis2 <- phis[,59:116][,order(final_perm)]
phis3 <- phis[,117:174][,order(final_perm)]
phis4 <- phis[,175:232][,order(final_perm)]
phis_origin <- cbind(phis1, phis2, phis3, phis4)

####### difference boundaries for each cancer individually #######
vij_samples1 <- t(apply(phis1, 1, function(x){
  diff_b1 <- NULL
  for(i in 1:nrow(neighbor_list0)){
    if(x[neighbor_list0[i,1]] != x[neighbor_list0[i,2]]){
      diff_b1 <- c(diff_b1, 1)
    }else{
      diff_b1 <- c(diff_b1, 0)
    }
  }
  return(diff_b1)
}))



vij_samples2 <- t(apply(phis2, 1, function(x){
  diff_b1 <- NULL
  for(i in 1:nrow(neighbor_list0)){
    if(x[neighbor_list0[i,1]] != x[neighbor_list0[i,2]]){
      diff_b1 <- c(diff_b1, 1)
    }else{
      diff_b1 <- c(diff_b1, 0)
    }
  }
  return(diff_b1)
}))

vij_samples3 <- t(apply(phis3, 1, function(x){
  diff_b1 <- NULL
  for(i in 1:nrow(neighbor_list0)){
    if(x[neighbor_list0[i,1]] != x[neighbor_list0[i,2]]){
      diff_b1 <- c(diff_b1, 1)
    }else{
      diff_b1 <- c(diff_b1, 0)
    }
  }
  return(diff_b1)
}))



vij_samples4 <- t(apply(phis4, 1, function(x){
  diff_b1 <- NULL
  for(i in 1:nrow(neighbor_list0)){
    if(x[neighbor_list0[i,1]] != x[neighbor_list0[i,2]]){
      diff_b1 <- c(diff_b1, 1)
    }else{
      diff_b1 <- c(diff_b1, 0)
    }
  }
  return(diff_b1)
}))


# probabilities
pvij1 <- apply(vij_samples1, 2, mean)
pvij2 <- apply(vij_samples2, 2, mean)
pvij3 <- apply(vij_samples3, 2, mean)
pvij4 <- apply(vij_samples4, 2, mean)

names(pvij1) <- neighbor_name
names(pvij2) <- neighbor_name
names(pvij3) <- neighbor_name
names(pvij4) <- neighbor_name

# Estimated FDR curves
T_edge1 = seq(0, 135, 1)[-1]
threshold1 = sort(pvij1, decreasing = TRUE)[T_edge1]
threshold2 = sort(pvij2, decreasing = TRUE)[T_edge1]
threshold3 = sort(pvij3, decreasing = TRUE)[T_edge1]
threshold4 = sort(pvij4, decreasing = TRUE)[T_edge1]
FDR_est1 <- rep(0, length(T_edge1))
FDR_est2 <- rep(0, length(T_edge1))
FDR_est3 <- rep(0, length(T_edge1))
FDR_est4 <- rep(0, length(T_edge1))


for(i in 1:length(threshold1)){
  th1 <- threshold1[i]
  est_diff1 <- as.numeric(pvij1 >= th1)
  FDR_est1[i] <- sum((1-pvij1) * est_diff1)  / sum(est_diff1)
  
  th2 <- threshold2[i]
  est_diff2 <- as.numeric(pvij2 >= th2)
  FDR_est2[i] <- sum((1-pvij2) * est_diff2)  / sum(est_diff2)
  
  th3 <- threshold3[i]
  est_diff3 <- as.numeric(pvij3 >= th3)
  FDR_est3[i] <- sum((1-pvij3) * est_diff3)  / sum(est_diff3)
  
  th4 <- threshold4[i]
  est_diff4 <- as.numeric(pvij4 >= th4)
  FDR_est4[i] <- sum((1-pvij4) * est_diff4)  / sum(est_diff4)
}

plot(FDR_est1,type="l", lty = 1, xlab="Number of edges selected", ylab = "Estimated FDR", ylim=c(0,0.5))
lines(FDR_est2,type="l", lty = 2)
lines(FDR_est3,type="l", lty = 3)
lines(FDR_est4,type="l", lty = 4)
legend("topleft", legend=c("Lung", "Esophageal", "Layrnx", "Colorectal"),
       lty=1:4, cex=0.7)

# Select a threshold for FDR and identify difference boundaries under the threshold
alpha_n = 0.025
T1 = sum(FDR_est1<=alpha_n)
T2 = sum(FDR_est2<=alpha_n)
T3 = sum(FDR_est3<=alpha_n)
T4 = sum(FDR_est4<=alpha_n)

est_diff1 <- as.numeric(pvij1 >= threshold1[T1])
est_diff1_1 <- as.numeric(pvij1 >= threshold1[75])
name_diff1_car <- names(pvij1[order(pvij1,decreasing = T)][1:75])

est_diff2 <- as.numeric(pvij2 >= threshold2[T2])
name_diff2_car <- names(pvij2[order(pvij2,decreasing = T)][1:75])


est_diff3 <- as.numeric(pvij3 >= threshold3[T3])
name_diff3_car <- names(pvij3[order(pvij3,decreasing = T)][1:75])


est_diff4 <- as.numeric(pvij4 >= threshold4[T4])
name_diff4_car <- names(pvij4[order(pvij4,decreasing = T)][1:75])

# Table for top 75 pairs of neighbors
name_diff_car = cbind(name_diff1_car, name_diff2_car, name_diff3_car, name_diff4_car)
colnames(name_diff_car) = c("Lung", "Esophageal", "Layrnx", "Coloretal")


# Plot boundaries
if (!require(gpclib)) install.packages("gpclib", type="source")
gpclibPermit()
p = list()
p_est=list()

library(plyr)
ca.poly@data$id <- rownames(ca.poly@data)
ca.poly.f <- fortify(ca.poly, region = "id")
ca.poly.df <- join(ca.poly.f, ca.poly@data, by = "id")

# Correct boundary data for plots
path=list()
for(i in 1: nrow(neighbor_list0)){
  r1 <- ca.poly.df[ca.poly.df$id %in% neighbor_list0[i,1],1:2]
  r2 <- ca.poly.df[ca.poly.df$id %in% neighbor_list0[i,2],1:2]
  edges <- generics::intersect(r1, r2)
  path[[i]] = edges
}
path[[2]][nrow(path[[2]])+1,] <- path[[2]][1,]
path[[2]] <- path[[2]][-1,]
path[[11]][nrow(path[[11]])+1,] <- path[[11]][1,]
path[[11]] <- path[[11]][-1,]
path[[19]][nrow(path[[19]])+1,] <- path[[19]][1,]
path[[19]] <- path[[19]][-1,]
path[[25]][nrow(path[[25]])+1,] <- path[[25]][1,]
path[[25]] <- path[[25]][-1,]
path[[27]][nrow(path[[27]])+1,] <- path[[27]][1,]
path[[27]] <- path[[27]][-1,]
path[[31]][nrow(path[[31]])+1,] <- path[[31]][1,]
path[[31]] <- path[[31]][-1,]
path[[39]][nrow(path[[39]])+1,] <- path[[39]][1,]
path[[39]] <- path[[39]][-1,]
path[[43]][nrow(path[[43]])+1,] <- path[[43]][1,]
path[[43]] <- path[[43]][-1,]
path[[46]][nrow(path[[46]])+1,] <- path[[46]][1,]
path[[46]] <- path[[46]][-1,]
path[[64]][nrow(path[[64]])+1,] <- path[[64]][1,]
path[[64]] <- path[[64]][-1,]
path[[75]][nrow(path[[75]])+1,] <- path[[75]][1,]
path[[75]] <- path[[75]][-1,]
path[[76]][nrow(path[[76]])+1,] <- path[[76]][1,]
path[[76]] <- path[[76]][-1,]
path[[81]][nrow(path[[81]])+1,] <- path[[81]][1,]
path[[81]] <- path[[81]][-1,]
path[[100]][nrow(path[[100]])+1,] <- path[[100]][1,]
path[[100]] <- path[[100]][-1,]
path[[125]][nrow(path[[125]])+1,] <- path[[125]][1,]
path[[125]] <- path[[125]][-1,]
path[[130]][nrow(path[[130]])+1,] <- path[[130]][1,]
path[[130]] <- path[[130]][-1,]
path[[133]][nrow(path[[133]])+1,] <- path[[133]][1,]
path[[133]] <- path[[133]][-1,]

#saveRDS(path, "path.rds")

edge_plot1 <- ggplot() +
  geom_polygon(data = ca.poly.df,  aes(long, lat, group = group),
               color = "black",
               fill = "white"
  )
for(i in which(est_diff1 == 1)){
  edge_plot1 = edge_plot1 + geom_path(aes_string(x = path[[i]][,1], y = path[[i]][,2]), color = "dodgerblue") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(size=16, face="bold",hjust = 0.5)) + 
    ggtitle(paste("Lung (T = ", T1, ")", sep=""))
}

for(i in which(est_diff1_1 == 1)){
  edge_plot1 = edge_plot1 + geom_path(aes_string(x = path[[i]][,1], y = path[[i]][,2]), color = "red")
}
edge_plot1

edge_plot2 <- ggplot() +
  geom_polygon(data = ca.poly.df,  aes(long, lat, group = group),
               color = "black",
               fill = "white"
  )
for(i in which(est_diff2 == 1)){
  edge_plot2 = edge_plot2 + geom_path(aes_string(x = path[[i]][,1], y = path[[i]][,2]), color = "red") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(size=16, face="bold",hjust = 0.5)) + 
    ggtitle(paste("Esophageal (T = ", T2, ")", sep=""))
}
edge_plot2


edge_plot3 <- ggplot() +
  geom_polygon(data = ca.poly.df,  aes(long, lat, group = group),
               color = "black",
               fill = "white"
  )
for(i in which(est_diff3 == 1)){
  edge_plot3 = edge_plot3 + geom_path(aes_string(x = path[[i]][,1], y = path[[i]][,2]), color = "red") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(size=16, face="bold",hjust = 0.5)) + 
    ggtitle(paste("Larynx (T = ", T3, ")", sep=""))
}
edge_plot3


edge_plot4 <- ggplot() +
  geom_polygon(data = ca.poly.df,  aes(long, lat, group = group),
               color = "black",
               fill = "white"
  )
for(i in which(est_diff4 == 1)){
  edge_plot4 = edge_plot4 + geom_path(aes_string(x = path[[i]][,1], y = path[[i]][,2]), color = "red") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(size=16, face="bold",hjust = 0.5)) + 
    ggtitle(paste("Colorectal (T = ", T4, ")", sep=""))
}
edge_plot4

#pdf("est_diff_cancer_car_pois.pdf", height = 10, width = 12)
ggarrange(edge_plot1, edge_plot2, edge_plot3, edge_plot4, nrow = 2, ncol = 2)
#dev.off()

# model fitting: D score
Y_rep1 <- matrix(0, nrow = n, ncol = 10000)
Y_rep2 <- matrix(0, nrow = n, ncol = 10000)
Y_rep3 <- matrix(0, nrow = n, ncol = 10000)
Y_rep4 <- matrix(0, nrow = n, ncol = 10000)
for(i in 1:10000){
  
  beta1 <- mcmc_samples$beta[i,1]
  beta2 <- mcmc_samples$beta[i,2]
  beta3 <- mcmc_samples$beta[i,3]
  beta4 <- mcmc_samples$beta[i,4]
  
  W1_mcmc = mcmc_samples$phi[i,1:58]
  W2_mcmc = mcmc_samples$phi[i,59:116]
  W3_mcmc = mcmc_samples$phi[i,117:174]
  W4_mcmc = mcmc_samples$phi[i,175:232]
  
  Y_rep1[,i] <- rpois(n, E1 * exp(X1[,1] * beta1 +  W1_mcmc)) / E1
  Y_rep2[,i] <- rpois(n, E2 * exp(X2[,1] * beta2 +  W2_mcmc)) / E2
  Y_rep3[,i] <- rpois(n, E3 * exp(X3[,1] * beta3 +  W3_mcmc)) / E3
  Y_rep4[,i] <- rpois(n, E4 * exp(X4[,1] * beta4 +  W4_mcmc)) / E4
  
}
mu_rep1 = rowMeans(Y_rep1)
mu_rep2 = rowMeans(Y_rep2)
mu_rep3 = rowMeans(Y_rep3)
mu_rep4 = rowMeans(Y_rep4)

var_rep1 = rowVars(Y_rep1)
var_rep2 = rowVars(Y_rep2)
var_rep3 = rowVars(Y_rep3)
var_rep4 = rowVars(Y_rep4)

G.latent1 = sum((Y1/E1 - mu_rep1)^2)
P.latent1 = sum(var_rep1)
D.latent1 = G.latent1 + P.latent1

G.latent2 = sum((Y2/E2 - mu_rep2)^2)
P.latent2 = sum(var_rep2)
D.latent2 = G.latent2 + P.latent2

G.latent3 = sum((Y3/E3 - mu_rep3)^2)
P.latent3 = sum(var_rep3)
D.latent3 = G.latent3 + P.latent3

G.latent4 = sum((Y4/E4 - mu_rep4)^2)
P.latent4 = sum(var_rep4)
D.latent4 = G.latent4 + P.latent4

D_CAR = c(D.latent1, D.latent2, D.latent3, D.latent4)

# You can also apply other analysis included in MDAGAR.R here for more details.


