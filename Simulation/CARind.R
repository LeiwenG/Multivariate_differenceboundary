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

#DAGAR
region = seq(1:n)
index = list()
for(i in 1:(n-2)){
  index[[i]] = region[-(udnei[(cn[i+1] + 1):(cn[i+1] + ns[i+1])])]
}
index1 = unlist(index)
mns = max(dni) + 1

#CAR
D=matrix(0,n_county,n_county)
for(i in 1:n) D[i,i]=ni[i]

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


## CAR
D=matrix(0,n_county,n_county)
for(i in 1:n) D[i,i]=ni[i]

# CAR_ind model
sink("MCAR_ARDP_ind.txt")
cat("
    model
    {
    
    for (i in 1:k)
    {
    mu[i] <- X[i,] %*% beta + phi[i]
    Y[i] ~ dnorm(mu[i], taue1)
    phi[i] <- inprod(theta, u[i,])
    Fr[i] <- pnorm(r[i], 0, 1/invQ1[i,i])
    for(h in 2:H){
    u[i, h] <- ifelse(Fr[i] > sum_pi[h-1] && Fr[i] < sum_pi[h], 1, 0)
    }
    u[i,1] <- ifelse(sum(u[i, 2:H])==0, 1, 0)
    }
    
    for (i in (k+1):(2*k))
    {
    mu[i] <- X[i,] %*% beta + phi[i]
    Y[i] ~ dnorm(mu[i], taue2)
    phi[i] <- inprod(theta, u[i,])
    Fr[i] <- pnorm(r[i], 0, 1/invQ2[i-k,i-k])
    for(h in 2:H){
    u[i, h] <- ifelse(Fr[i] > sum_pi[h-1] && Fr[i] < sum_pi[h], 1, 0)
    }
    u[i,1] <- ifelse(sum(u[i, 2:H])==0, 1, 0)
    }
    
    pi[1] <- V[1]
    prod[1] <- 1 - pi[1]
    for (h in 2:H){
    prod[h] <- prod[h-1] * (1-V[h])
    pi[h] <- V[h] * prod[h-1]
    }
    
    for (h in 1:H){
    theta[h] ~ dnorm(0,taus)
    V[h] ~ dbeta(1,alpha)
    sum_pi[h] <- sum(pi[1:h])
    }
    
    Q1 <- tau1 * (D - rho1 * Minc)
    Q2 <- tau2 * (D - rho2 * Minc)
    r1[1:k] ~ dmnorm(rep(0, k), Q1)
    r2[1:k] ~ dmnorm(rep(0, k), Q2)
    r[1:k] <- r1[1:k]
    r[(k+1):(2*k)] <- r2[1:k]
    
    invQ1 <- inverse(Q1)
    invQ2 <- inverse(Q2)
    
    rho1 ~ dunif(0, 0.98)
    rho2 ~ dunif(0, 0.98)
    tau1 ~ dgamma(2, 0.1)
    tau2 ~ dgamma(2, 0.1)
    taue1 ~ dgamma(2, 0.1)
    taue2 ~ dgamma(2, 0.1)
    vare1 <- 1/taue1
    vare2 <- 1/taue2
    taus ~ dgamma(2,0.1)
    beta[1:ncolumn] ~ dmnorm(rep(0,ncolumn), (0.0001*I));
    }
    ", fill = TRUE)
sink()

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

# Generate 50 datasets
for(seed in 1:50){
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
  #y = X %*% beta + phi_true  + sqrt(sigma_e2) * rnorm(n_county)
  Xo1 <- X1[final_perm,]
  yo1 <- y1[final_perm]
  Xo2 <- X2[final_perm,]
  yo2 <- y2[final_perm]
  
  Y1 = c(yo1, yo2)
  Xb1 = as.matrix(bdiag(Xo1, Xo2))
  
  model.data1 <- list(k = n_county,  I = diag(ncol(Xb1)), X = Xb1, Y = Y1,
                      alpha = 1, ncolumn = ncol(Xb1),  H = 15, D = D, Minc = Minc)
  model.inits <- rep(list(list(rho1 = 0.1, rho2 = 0.1, tau1 = 1, tau2 = 1, taue1 = 1, taue2 = 1, 
                               taus = 1, beta = rep(0, ncol(Xb1)))),2)
  model.param <- c("beta", "rho1", "rho2", "tau1", "tau2", "vare1", "vare2", 
                   "taus", "phi", "u")
  
  run.CAR1 <- jags(model.data1, model.inits, model.param, "MCAR_ARDP_ind.txt",
                   n.chains = 2, n.iter = 25000,n.burnin = 20000, n.thin = 1)
  mcmc_samples[[seed]] <- run.CAR1$BUGSoutput$sims.matrix
  
  phis <-  mcmc_samples[[seed]][,6:121]
  phis1 <- phis[,1:58][,order(final_perm)]
  phis2 <- phis[,59:116][,order(final_perm)]
  phis_origin <- cbind(phis1, phis2)
  table(apply(phis2, 1, function(x) length(unique(x))))
  table(apply(phis1, 1, function(x) length(unique(x))))
  
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
  
  # calculate sensitivity and specificity
  for(i in 1:length(T_edge)){
    est_diff1 <- factor(as.numeric(pvijm[1,] >= threshold1[i]), levels = c(0, 1))
    conf_matrix1 <-table(est_diff1,true_diff1)
    spec1[seed,i] <- sensitivity(conf_matrix1)
    sens1[seed,i] <- specificity(conf_matrix1)
    
    est_diff2 <- factor(as.numeric(pvijm[2,] >= threshold2[i]), levels = c(0, 1))
    conf_matrix2 <-table(est_diff2,true_diff2)
    spec2[seed,i] <- sensitivity(conf_matrix2)
    sens2[seed,i] <- specificity(conf_matrix2)
    
    est_diff12 <- factor(as.numeric(pvijm[3,] >= threshold12[i]), levels = c(0, 1))
    conf_matrix12 <-table(est_diff12,true_diff3)
    spec12[seed,i] <- sensitivity(conf_matrix12)
    sens12[seed,i] <- specificity(conf_matrix12)
    
    est_diff21 <- factor(as.numeric(pvijm[4,] >= threshold21[i]), levels = c(0, 1))
    conf_matrix21<-table(est_diff21,true_diff4)
    spec21[seed,i] <- sensitivity(conf_matrix21)
    sens21[seed,i] <- specificity(conf_matrix21)
  }
  
  save.image("simBARDP_car_indc25.RData")
}

