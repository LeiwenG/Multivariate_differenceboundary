library(maps)
ca.county = map("county","california", fill=TRUE, plot=FALSE)
library(spdep)
library(maptools)
library(mapproj)
library(stringr)
library(classInt)
library(RColorBrewer)
library(rjags)
library(R2jags)
library(dplyr)
library(caret)
library(Matrix)

county.ID <- sapply(strsplit(ca.county$names, ","), function(x) x[2])
ca.poly = map2SpatialPolygons(ca.county, IDs=county.ID)
ca.coords = coordinates(ca.poly)

n_county <- length(county.ID)


###################################
### Assign true mean on the map ###
###################################

# discrete random effects for disease 1

phi_true1 <- readRDS("phi_true1_25.rds")

# discrete random effects for disease 2

phi_true2 <- readRDS("phi_true2_25.rds")


## Adjacency matrix
ca.neighbors = poly2nb(ca.poly)
n=length(ca.neighbors)

Adj=sapply(ca.neighbors,function(x,n) {v=rep(0,n);v[x]=1;v},n)
colnames(Adj)=county_id

num_edge <- sum(Adj)/2

neighborvec0 = NULL
neighbor_list0 = NULL
for(i in 1:(n_county-1)){
  for(j in (i+1):n_county){
    if(Adj[i,j] == 1){
      neighborvec0 = c(neighborvec0, paste(i, ",", j, sep=""))
      neighbor_list0 = rbind(neighbor_list0, c(i,j))
    }
  }
}

## Reorder the map
ca.latrange=round(quantile(ca.coords[,2],c(0.25,0.75)))
ca.albersproj=mapproject(ca.coords[,1],ca.coords[,2],projection = "albers",param=ca.latrange)

perm=order(ca.albersproj$x-ca.albersproj$y)
colnames(Adj)[perm]

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

region = seq(1:n)
index = list()
for(i in 1:(n-2)){
  index[[i]] = region[-(udnei[(cn[i+1] + 1):(cn[i+1] + ns[i+1])])]
}
index1 = unlist(index)
mns = max(dni) + 1

#######true difference boundary########

true_diff1 <- NULL
true_diff2 <- NULL
true_diff3 <- NULL
true_diff4 <- NULL
for(i in 1:nrow(neighbor_list0)){
  if(phi_true1[neighbor_list0[i,1]] != phi_true1[neighbor_list0[i,2]]){
    true_diff1 <- c(true_diff1, 1)
  }else{
    true_diff1 <- c(true_diff1, 0)
  }
  if(phi_true2[neighbor_list0[i,1]] != phi_true2[neighbor_list0[i,2]]){
    true_diff2 <- c(true_diff2, 1)
  }else{
    true_diff2 <- c(true_diff2, 0)
  }
  if(phi_true1[neighbor_list0[i,1]] != phi_true2[neighbor_list0[i,2]]){
    true_diff3 <- c(true_diff3, 1)
  }else{
    true_diff3 <- c(true_diff3, 0)
  }
  if(phi_true2[neighbor_list0[i,1]] != phi_true1[neighbor_list0[i,2]]){
    true_diff4 <- c(true_diff4, 1)
  }else{
    true_diff4 <- c(true_diff4, 0)
  }
}

sum(true_diff1)
sum(true_diff2)
sum(true_diff3)
sum(true_diff4)

true_diffm <- rbind(true_diff1, true_diff2, true_diff3, true_diff4)



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

### MDAGAR model
# Metropolis within Gibbs Sampler for MCMC updating 
ARDP_joint_diff<-function(y=Y1, x1=Xo1, x2=Xo2, X = Xb1, Minc, alpha=1, q = 2, n.atoms=15, runs=10000, burn=1000){
  #y:       data
  #x:       covariates
  #n.atoms: number of atoms in the mixture dist.
  #theta:      the theta's (iid from baseline) in the model
  #alpha:     v~beta(1,alpha)
  #u:       the index indicator of spatial random effect 
  #rho:     sptatial autocorrelation parameter in DAGAR
  #Minc:    0-1 adjacency matrix 
  #V_r:     covariance matrix of joint MDAGAR
  #Q:       presicion matrix of independent DAGAR
  #r:       random effects following DAGAR
  #F_r:     Marginal CDF of r
  #taued:   presicion in Baseline of DP for disease d   
  #taus:    precision for theta
  
  
  
  
  nq<-length(y)
  n <- nq / q
  p<-ncol(X)
  y1 <- y[1:n]
  y2 <- y[(n+1):(2*n)]
  
  sigmasq_beta = 10000
  keepbeta=matrix(0,runs,p)
  keepphi=matrix(0,runs,nq)
  keeptheta=matrix(0,runs,n.atoms)
  keepu=matrix(0,runs,nq)
  keeprho1=keeprho2=keeptaue1=keeptaue2 = keeptaus = rep(0,runs)
  keepv=matrix(0,runs,n.atoms)
  keepr=matrix(0,runs,nq)
  keepFr=matrix(0,runs,nq)
  
  
  #initial values
  
  theta=rep(0,n.atoms)
  beta<-rep(0,p)
  taue1=1
  taue2=1
  taus = 1
  
  c=2
  d=0.1
  d2 = 1
  v<-rep(.5,n.atoms)
  probs=makeprobs(v)
  
  #A matrix
  rho=c(0.1, 0.1)
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
  
  Q = Dinv_new(Rho = rho, n, cn, ns, udnei,q=2)
  invQ1 = solve(Q[[1]])
  invQ2 = solve(Q[[2]])
  f1 = rmvnorm(1, rep(0, n), invQ1)
  f2 = rmvnorm(1, rep(0, n), invQ2)
  f = c(f1, f2)
  kprod = kronecker(A, diag(n))
  invQ = as.matrix(bdiag(invQ1, invQ2))
  r = kprod %*% f
  
  Vr = as.matrix(forceSymmetric(kprod %*% invQ %*% t(kprod)))
  #r = rmvnorm(1, rep(0, nq), Vr)
  
  
  F_r=pnorm(r,0,sqrt(diag(Vr)))
  u=makeu(F_r,probs)
  phi <- theta[u]
  
  nu = 2
  R = 0.1 * diag(q)
  
  
  acceptr=acceptv=acceptrho=acceptA=0
  
  count<-afterburn<-0
  burn = burn + 1

  for(iter in 1:runs){
    #print(iter)
    if(iter %% 100 == 0){
      print(iter)
      print(acceptA/(iter-1)) 
      print(acceptr/116/(iter-1))
      print(acceptrho/(iter-1)) 
      print(acceptv/n.atoms/(iter-1))
    }
    ######################
    ###   update beta  ###
    ###################### 
    # update beta (intercept only model)
    
    M1 = solve(taue1 * t(x1) %*% x1 + 1/sigmasq_beta * diag(ncol(x1)))
    m1 = taue1 *  t(x1) %*% (y1 - theta[u][1:n])
    mu1 = M1 %*% m1
    
    M2 = solve(taue2 * t(x2) %*% x2 + 1/sigmasq_beta * diag(ncol(x2)))
    m2 = taue2 *  t(x2) %*% (y2 - theta[u][(n+1):(2*n)])
    mu2 = M2 %*% m2
    
    beta1 <- rmvnorm(1, mu1, M1)
    beta2 <- rmvnorm(1, mu2, M2)
    
    beta <- c(beta1, beta2)
    
    
    #########################
    ###   update theta    ###
    ######################### 
    
    u1 = u[1:n]
    u2 = u[(n+1):(2*n)]
    
    for (j in 1:n.atoms){
      M = 1/(taue1 * sum(u1 == j) + taue2 * sum(u2 == j) + taus)
      m = taue1 * sum(y1[u1 == j] - as.vector(x1[u1 == j, ] %*% t(beta1))) + taue2 * sum(y2[u2 == j] - as.vector(x2[u2 == j, ] %*% t(beta2)))
      theta[j] <- rnorm(1, M*m, sqrt(M))
    }
    
    
    
    ######################
    ###   update r     ###
    ###################### 
    
    
    tauvec = c(rep(taue1, n), rep(taue2, n))
    
    for (k in 1:nq){
      pro_r=r;pro_Fr=F_r;pro_u=u
      pro_r[k]=rnorm(1,r[k],1.8)
      pro_Fr[k]=pnorm(pro_r[k],0,sqrt(Vr[k,k]))			
      pro_u[k]=makeu(pro_Fr[k],probs)
      
      MH=dmvnorm(t(pro_r),mean=rep(0, nq),sigma=Vr,log=T)+dnorm(y[k],as.numeric(X[k,]%*%beta)+theta[pro_u[k]], sqrt(1/tauvec[k]),log=T)-
        dmvnorm(t(r),mean=rep(0, nq),sigma=Vr,log=T)-dnorm(y[k],as.numeric(X[k,]%*%beta)+theta[u[k]], sqrt(1/tauvec[k]),log=T) 
      
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
      pro_v[j]<-rnorm(1,v[j],0.45)
      if(pro_v[j] > 0 & pro_v[j] < 1){
        pro_probs=makeprobs(pro_v)
        pro_u=makeu(F_r,pro_probs)
        
        MH=log(dbeta(pro_v[j],1,alpha))+dmvnorm(t(as.matrix(y1)),mean=as.vector(x1%*%t(beta1))+theta[pro_u[1:n]],sigma=1/taue1*diag(n),log=T)+
          dmvnorm(t(as.matrix(y2)),mean=as.vector(x2%*%t(beta2))+theta[pro_u[(n+1):(2*n)]],sigma=1/taue2*diag(n),log=T)-
          log(dbeta(v[j],1,alpha))-dmvnorm(t(as.matrix(y1)),mean=as.vector(x1%*%t(beta1))+theta[u[1:n]],sigma=1/taue1*diag(n),log=T) -
          dmvnorm(t(as.matrix(y2)),mean=as.vector(x2%*%t(beta2))+theta[u[(n+1):(2*n)]],sigma=1/taue2*diag(n),log=T)
        
        
        if(runif(1,0,1)<exp(MH)){
          v[j]=pro_v[j]
          probs=pro_probs
          u=pro_u
          acceptv=acceptv+1	
        }
      }
    }
    
    ######################
    ###   update taued  ###
    ###################### 
    
    taue1=rgamma(1,shape=n/2+c,rate=t(y1-as.vector(x1%*%t(beta1))-theta[u[1:n]])%*%(y1-as.vector(x1%*%t(beta1))-theta[u[1:n]])/2+d)
    taue2=rgamma(1,shape=n/2+c,rate=t(y2-as.vector(x2%*%t(beta2))-theta[u[(n+1):(2*n)]])%*%(y2-as.vector(x2%*%t(beta2))-theta[u[(n+1):(2*n)]])/2+d)
    
    
    ######################
    ###   update taus  ###
    ###################### 
    
    
    taus=rgamma(1,shape=1/2*n.atoms+c,rate=sum(theta^2)/2+d)
    
    
    ######################
    ###   update rho  ###
    ###################### 
    
    #0.8
    pro_eta1 = rnorm(1, eta[1], 0.8)
    pro_eta2 = rnorm(1, eta[2], 0.8)
    pro_eta = c(pro_eta1, pro_eta2)
    pro_rho = exp(pro_eta)/(1+exp(pro_eta))
    pro_Q = Dinv_new(Rho = pro_rho, n, cn, ns, udnei,q=2)
    pro_invQ1 = solve(pro_Q[[1]])
    pro_invQ2 = solve(pro_Q[[2]])

    kprod = kronecker(A, diag(n))
    pro_invQ = as.matrix(bdiag(pro_invQ1, pro_invQ2))
    pro_Vr = kprod %*% pro_invQ %*% t(kprod)
    
    MH = dmvnorm(t(r),mean=rep(0, nq),sigma=as.matrix(forceSymmetric(pro_Vr)),log=T) + log(pro_rho[1]) + log(1-pro_rho[1]) + log(pro_rho[2]) + log(1-pro_rho[2]) - 
      dmvnorm(t(r),mean=rep(0, nq),sigma=Vr,log=T) - log(rho[1]) - log(1-rho[1]) - log(rho[2]) - log(1-rho[2])
    
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
          pro_A[i, j] = exp(rnorm(1, log(A[i, j]), 0.1))
          #pro_A[i, j] = rlnorm(1, log(A[i, j]), 0.1)
        }else{
          pro_A[i, j] = rnorm(1, A[i, j], 0.1)
        }
      }
    }
    
    Q = Dinv_new(Rho = rho, n, cn, ns, udnei,q=2)
    invQ1 = solve(Q[[1]])
    invQ2 = solve(Q[[2]])
    Sigma =  A %*% t(A)
    #f1 = rmvnorm(1, rep(0, n), invQ1)
    #f2 = rmvnorm(1, rep(0, n), invQ2)
    #f = c(f1, f2)
    pro_kprod = kronecker(pro_A, diag(n))
    invQ = as.matrix(bdiag(invQ1, invQ2))
    pro_Vr = pro_kprod %*% as.matrix(forceSymmetric(invQ)) %*% t(pro_kprod)
    pro_Sigma =  pro_A %*% t(pro_A)
    #pro_Vr = as.matrix(forceSymmetric(kronecker(pro_Sigma, solve(Q))))
    
    lpA = -(nu+4)/2*logdet(Sigma) - 1/2*sum(diag(nu*R%*%solve(Sigma))) + log(Jacob_A(A)) + sum(log(diag(A)))
    pro_lpA = -(nu+4)/2*logdet(pro_Sigma) - 1/2*sum(diag(nu*R%*%solve(pro_Sigma))) + log(Jacob_A(pro_A)) + sum(log(diag(pro_A)))
    MH = dmvnorm(t(r),mean=rep(0, nq),sigma=as.matrix(forceSymmetric(pro_Vr)),log=T) + pro_lpA - 
      dmvnorm(t(r),mean=rep(0, nq),sigma=Vr,log=T) - lpA
    
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
    
    
    keepbeta[iter,]=beta
    keeptaus[iter]=taus
    keeprho1[iter]=rho[1]
    keeprho2[iter]=rho[2]
    keeptaue1[iter]=taue1
    keeptaue2[iter]=taue2
    keepu[iter,]=u
    keepv[iter,]=v
    keepr[iter,]=r
    keepFr[iter,]=F_r
    
    #	cat("iteration = ", i, "acceptance rate of r = ", acceptr/(n*i),"acceptance rate of v = ",acceptv/(n.atoms*i), "\n")
    
  }
  
  list(beta=keepbeta[burn:runs,],phi=keepphi[burn:runs,],theta=keeptheta[burn:runs,],u=keepu[burn:runs,],
       v=keepv[burn:runs,],r=keepr[burn:runs,],Fr=keepFr[burn:runs,],taue1=keeptaue1[burn:runs],
       taue2=keeptaue2[burn:runs],taus=keeptaus[burn:runs],rho1=keeprho1[burn:runs], rho2=keeprho2[burn:runs])
}

T_edge = c(60, 65, 70, 75, 80, 85, 90, 95, 100, 105) 
spec1 <- matrix(0, 30, length(T_edge))
sens1 <- matrix(0, 30, length(T_edge))
spec2 <- matrix(0, 30, length(T_edge))
sens2 <- matrix(0, 30, length(T_edge))
spec12 <- matrix(0, 30, length(T_edge))
sens12 <- matrix(0, 30, length(T_edge))
spec21 <- matrix(0, 30, length(T_edge))
sens21 <- matrix(0, 30, length(T_edge))
mcmc_samples <- list()

# Simulate 30 datasets
for(seed in 1:30){
  print(seed)

  set.seed(seed)
  X1=cbind(1, rnorm(n_county))
  X2=cbind(1, rnorm(n_county))
  beta1=c(2,5)
  beta2=c(1,6)
  #beta = c(beta1, beta2)
  sigma_esq1 = 0.1
  sigma_esq2 = 0.1
  y1 = X1 %*% beta1 + (phi_true1 - mean(phi_true1))  + sqrt(sigma_esq1) * rnorm(n_county)
  y2 = X2 %*% beta2 + (phi_true2 - mean(phi_true2))  + sqrt(sigma_esq2) * rnorm(n_county)

  Xo1 <- X1[final_perm,]
  yo1 <- y1[final_perm]
  Xo2 <- X2[final_perm,]
  yo2 <- y2[final_perm]
  
  Y1 = c(yo1, yo2)
  Xb1 = as.matrix(bdiag(Xo1, Xo2))
  
  mcmc_samples[[seed]]=ARDP_joint_diff(y=Y1, x1=Xo1, x2=Xo2, X = Xb1, Minc, alpha=1, q = 2,
                                       n.atoms=15, runs=25000, burn=15000)
  
  phis <-  mcmc_samples[[seed]]$phi
  phis1 <- phis[,1:58][,order(final_perm)]
  phis2 <- phis[,59:116][,order(final_perm)]
  phis_origin <- cbind(phis1, phis2)
  
  # estimate difference boundaries
  vij_samples1 <- t(apply(phis1, 1, function(x){
    #x <- phis[1,]
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
    #x <- phis[1,]
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
  
  vij_samples12 <- t(apply(phis_origin, 1, function(x){
    #x <- phis[1,]
    diff_b1 <- NULL
    for(i in 1:nrow(neighbor_list0)){
      if(x[neighbor_list0[i,1]] != x[neighbor_list0[i,2]+58]){
        diff_b1 <- c(diff_b1, 1)
      }else{
        diff_b1 <- c(diff_b1, 0)
      }
    }
    return(diff_b1)
  }))
  
  vij_samples21 <- t(apply(phis_origin, 1, function(x){
    #x <- phis[1,]
    diff_b1 <- NULL
    for(i in 1:nrow(neighbor_list0)){
      if(x[neighbor_list0[i,1]+58] != x[neighbor_list0[i,2]]){
        diff_b1 <- c(diff_b1, 1)
      }else{
        diff_b1 <- c(diff_b1, 0)
      }
    }
    return(diff_b1)
  }))
  
  
  pvij1 <- apply(vij_samples1, 2, mean)
  pvij2 <- apply(vij_samples2, 2, mean)
  pvij12 <- apply(vij_samples12, 2, mean)
  pvij21 <- apply(vij_samples21, 2, mean)
  pvijm <- rbind(pvij1, pvij2, pvij12, pvij21)
  
  threshold1 = sort(pvij1, decreasing = TRUE)[T_edge]
  threshold2 = sort(pvij2, decreasing = TRUE)[T_edge]
  threshold12 = sort(pvij12, decreasing = TRUE)[T_edge]
  threshold21 = sort(pvij21, decreasing = TRUE)[T_edge]
  thresholdm <- rbind(threshold1, threshold2, threshold12, threshold21)
  
  # calculate sensitivity and specificity
  for(i in 1:length(T_edge)){
    est_diff1 <- factor(as.numeric(pvijm[1,] >= threshold1[i]), levels = c(0,1))
    conf_matrix1 <-table(est_diff1,true_diff1)
    spec1[seed,i] <- sensitivity(conf_matrix1)
    sens1[seed,i] <- specificity(conf_matrix1)
    
    est_diff2 <- factor(as.numeric(pvijm[2,] >= threshold2[i]), levels = c(0,1))
    conf_matrix2 <-table(est_diff2,true_diff2)
    spec2[seed,i] <- sensitivity(conf_matrix2)
    sens2[seed,i] <- specificity(conf_matrix2)
    
    est_diff12 <- factor(as.numeric(pvijm[3,] >= threshold12[i]), levels = c(0,1))
    conf_matrix12 <-table(est_diff12,true_diff3)
    spec12[seed,i] <- sensitivity(conf_matrix12)
    sens12[seed,i] <- specificity(conf_matrix12)
    
    est_diff21 <- factor(as.numeric(pvijm[4,] >= threshold21[i]), levels = c(0,1))
    conf_matrix21<-table(est_diff21,true_diff4)
    spec21[seed,i] <- sensitivity(conf_matrix21)
    sens21[seed,i] <- specificity(conf_matrix21)
  }
  
  save.image("simBARDP_dagarc25.RData")
}

