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

## CAR
D=matrix(0,n,n)
for(i in 1:n) D[i,i]=ni[i]

# Data in model
Y1 = rate_lung1$count[final_perm]
Y2 = rate_esophagus1$count[final_perm]
Y3 = rate_larynx1$count[final_perm]
Y4 = rate_colrect1$count[final_perm]

E1 = rate_lung1$E_count[final_perm]
E2 = rate_esophagus1$E_count[final_perm]
E3 = rate_larynx1$E_count[final_perm]
E4 = rate_colrect1$E_count[final_perm]

X1 = as.matrix(cbind(1,rate_lung1[,8:14]))[final_perm,]
X2 = as.matrix(cbind(1,rate_esophagus1[,8:14]))[final_perm,]
X3 = as.matrix(cbind(1,rate_larynx1[,8:14]))[final_perm,]
X4 = as.matrix(cbind(1,rate_colrect1[,8:14]))[final_perm,]

Y = c(Y1,Y2,Y3,Y4)
E = c(E1, E2, E3, E4)
X = as.matrix(bdiag(bdiag(X1[,1], X2[,1]), bdiag(X3[,1],X4[,1])))

# Set up model in JAGS
sink("MCAR_ARDP_ind.txt")
cat("
    model
    {
    
    for (i in 1:k)
    {
    log(mu[i]) <- log(E[i]) + X[i,] %*% beta + phi[i]
    Y[i] ~ dpois(mu[i])
    phi[i] <- inprod(theta, u[i,])
    Fr[i] <- pnorm(r[i], 0, 1/invQ1[i,i])
    for(h in 2:H){
    u[i, h] <- ifelse(Fr[i] > sum_pi[h-1] && Fr[i] < sum_pi[h], 1, 0)
    }
    u[i,1] <- ifelse(sum(u[i, 2:H])==0, 1, 0)
    }
    
    for (i in (k+1):(2*k))
    {
    log(mu[i]) <- log(E[i]) + X[i,] %*% beta + phi[i]
    Y[i] ~ dpois(mu[i])
    phi[i] <- inprod(theta, u[i,])
    Fr[i] <- pnorm(r[i], 0, 1/invQ2[i-k,i-k])
    for(h in 2:H){
    u[i, h] <- ifelse(Fr[i] > sum_pi[h-1] && Fr[i] < sum_pi[h], 1, 0)
    }
    u[i,1] <- ifelse(sum(u[i, 2:H])==0, 1, 0)
    }
    
    for (i in (2*k+1):(3*k))
    {
    log(mu[i]) <- log(E[i]) + X[i,] %*% beta + phi[i]
    Y[i] ~ dpois(mu[i])
    phi[i] <- inprod(theta, u[i,])
    Fr[i] <- pnorm(r[i], 0, 1/invQ3[i-2*k,i-2*k])
    for(h in 2:H){
    u[i, h] <- ifelse(Fr[i] > sum_pi[h-1] && Fr[i] < sum_pi[h], 1, 0)
    }
    u[i,1] <- ifelse(sum(u[i, 2:H])==0, 1, 0)
    }
    
    for (i in (3*k+1):(4*k))
    {
    log(mu[i]) <- log(E[i]) + X[i,] %*% beta + phi[i]
    Y[i] ~ dpois(mu[i])
    phi[i] <- inprod(theta, u[i,])
    Fr[i] <- pnorm(r[i], 0, 1/invQ4[i-3*k,i-3*k])
    for(h in 2:H){
    u[i, h] <- ifelse(Fr[i] > sum_pi[h-1] && Fr[i] < sum_pi[h], 1, 0)
    }
    u[i,1] <- ifelse(sum(u[i, 2:H])==0, 1, 0)
    }
    
    
    pi[1] <- V[1]
    prod[1] <- 1 - pi[1]
    for (h in 2:(H-1)){
    prod[h] <- prod[h-1] * (1-V[h])
    pi[h] <- V[h] * prod[h-1]
    }
    pi[H] <- 1-sum(pi[1:(H-1)])
    
    for (h in 1:H){
    theta[h] ~ dnorm(0,taus)
    V[h] ~ dbeta(1,alpha)
    sum_pi[h] <- sum(pi[1:h])
    }
    
    Q1 <- tau1 * (D - rho1 * Minc)
    Q2 <- tau2 * (D - rho2 * Minc)
    Q3 <- tau3 * (D - rho3 * Minc)
    Q4 <- tau4 * (D - rho4 * Minc)
    
    r1[1:k] ~ dmnorm(rep(0, k), Q1)
    r2[1:k] ~ dmnorm(rep(0, k), Q2)
    r3[1:k] ~ dmnorm(rep(0, k), Q3)
    r4[1:k] ~ dmnorm(rep(0, k), Q4)
    r[1:k] <- r1[1:k]
    r[(k+1):(2*k)] <- r2[1:k]
    r[(2*k+1):(3*k)] <- r3[1:k]
    r[(3*k+1):(4*k)] <- r4[1:k]
    
    invQ1 <- inverse(Q1)
    invQ2 <- inverse(Q2)
    invQ3 <- inverse(Q3)
    invQ4 <- inverse(Q4)
    
    rho1 ~ dunif(0, 0.98)
    rho2 ~ dunif(0, 0.98)
    rho3 ~ dunif(0, 0.98)
    rho4 ~ dunif(0, 0.98)
    tau1 ~ dgamma(2, 0.1)
    tau2 ~ dgamma(2, 0.1)
    tau3 ~ dgamma(2, 0.1)
    tau4 ~ dgamma(2, 0.1)
    
    taus ~ dgamma(2,1)
    beta[1:ncolumn] ~ dmnorm(rep(0,ncolumn), (0.0001*I));
    }
    ", fill = TRUE)
sink()

model.data1 <- list(k = n,  I = diag(ncol(X)), X = X, Y = Y, E = E, 
                    alpha = 1, ncolumn = ncol(X),  H = 15,  D = D, Minc = Minc)
model.inits <- rep(list(list(rho1 = 0.1, rho2 = 0.1, rho3 = 0.1, rho4 = 0.1, tau1 = 1, tau2 = 1, tau3 = 1, tau4 = 1,  
                             taus = 1, beta = rep(0, ncol(X)))),2)
model.param <- c("beta", "rho1", "rho2", "rho3", "rho4", "tau1", "tau2", "tau3", "tau4", 
                 "taus", "phi", "u")
set.seed(123)
mcmc_samples <- jags(model.data1, model.inits, model.param, "MCAR_ARDP_ind.txt",
                 n.chains = 2, n.iter = 25000,n.burnin = 20000, n.thin = 1)

sample.mcmc = mcmc_samples$BUGSoutput$sims.matrix[,c(246, 238:245)]
phis <- mcmc_samples$BUGSoutput$sims.matrix[,6:237]

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
alpha_n = 0.1
T1 = sum(FDR_est1<=alpha_n)
T2 = sum(FDR_est2<=alpha_n)
T3 = sum(FDR_est3<=alpha_n)
T4 = sum(FDR_est4<=alpha_n)

est_diff1 <- as.numeric(pvij1 >= threshold1[T1])
est_diff2 <- as.numeric(pvij2 >= threshold2[T2])
est_diff3 <- as.numeric(pvij3 >= threshold3[T3])
est_diff4 <- as.numeric(pvij4 >= threshold4[T4])

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

#pdf("est_diff_cancer_car_ind_pois.pdf", height = 10, width = 12)
ggarrange(edge_plot1, edge_plot2, edge_plot3, edge_plot4, nrow = 2, ncol = 2)
#dev.off()

########### Shared difference boundary for each pair of cancers ##########
vij_samplesc12 <- t(apply(cbind(phis1,phis2), 1, function(x){
  diff_bc1 <- NULL
  for(i in 1:nrow(neighbor_list0)){
    if(x[neighbor_list0[i,1]] != x[neighbor_list0[i,2]] & x[neighbor_list0[i,1]+58] != x[neighbor_list0[i,2]+58]){
      diff_bc1 <- c(diff_bc1, 1)
    }else{
      diff_bc1 <- c(diff_bc1, 0)
    }
  }
  return(diff_bc1)
}))

vij_samplesc13 <- t(apply(cbind(phis1,phis3), 1, function(x){
  diff_bc1 <- NULL
  for(i in 1:nrow(neighbor_list0)){
    if(x[neighbor_list0[i,1]] != x[neighbor_list0[i,2]] & x[neighbor_list0[i,1]+58] != x[neighbor_list0[i,2]+58]){
      diff_bc1 <- c(diff_bc1, 1)
    }else{
      diff_bc1 <- c(diff_bc1, 0)
    }
  }
  return(diff_bc1)
}))

vij_samplesc14 <- t(apply(cbind(phis1,phis4), 1, function(x){
  diff_bc1 <- NULL
  for(i in 1:nrow(neighbor_list0)){
    if(x[neighbor_list0[i,1]] != x[neighbor_list0[i,2]] & x[neighbor_list0[i,1]+58] != x[neighbor_list0[i,2]+58]){
      diff_bc1 <- c(diff_bc1, 1)
    }else{
      diff_bc1 <- c(diff_bc1, 0)
    }
  }
  return(diff_bc1)
}))

vij_samplesc23 <- t(apply(cbind(phis2,phis3), 1, function(x){
  diff_bc1 <- NULL
  for(i in 1:nrow(neighbor_list0)){
    if(x[neighbor_list0[i,1]] != x[neighbor_list0[i,2]] & x[neighbor_list0[i,1]+58] != x[neighbor_list0[i,2]+58]){
      diff_bc1 <- c(diff_bc1, 1)
    }else{
      diff_bc1 <- c(diff_bc1, 0)
    }
  }
  return(diff_bc1)
}))

vij_samplesc24 <- t(apply(cbind(phis2,phis4), 1, function(x){
  diff_bc1 <- NULL
  for(i in 1:nrow(neighbor_list0)){
    if(x[neighbor_list0[i,1]] != x[neighbor_list0[i,2]] & x[neighbor_list0[i,1]+58] != x[neighbor_list0[i,2]+58]){
      diff_bc1 <- c(diff_bc1, 1)
    }else{
      diff_bc1 <- c(diff_bc1, 0)
    }
  }
  return(diff_bc1)
}))

vij_samplesc34 <- t(apply(cbind(phis3,phis4), 1, function(x){
  diff_bc1 <- NULL
  for(i in 1:nrow(neighbor_list0)){
    if(x[neighbor_list0[i,1]] != x[neighbor_list0[i,2]] & x[neighbor_list0[i,1]+58] != x[neighbor_list0[i,2]+58]){
      diff_bc1 <- c(diff_bc1, 1)
    }else{
      diff_bc1 <- c(diff_bc1, 0)
    }
  }
  return(diff_bc1)
}))
# probabilities
pvijc12 <- apply(vij_samplesc12, 2, mean)
pvijc13 <- apply(vij_samplesc13, 2, mean)
pvijc14 <- apply(vij_samplesc14, 2, mean)
pvijc23 <- apply(vij_samplesc23, 2, mean)
pvijc24 <- apply(vij_samplesc24, 2, mean)
pvijc34 <- apply(vij_samplesc34, 2, mean)


# Estimated FDR curves
thresholdc12 = sort(pvijc12, decreasing = TRUE)[T_edge1]
thresholdc13 = sort(pvijc13, decreasing = TRUE)[T_edge1]
thresholdc14 = sort(pvijc14, decreasing = TRUE)[T_edge1]
thresholdc23 = sort(pvijc23, decreasing = TRUE)[T_edge1]
thresholdc24 = sort(pvijc24, decreasing = TRUE)[T_edge1]
thresholdc34 = sort(pvijc34, decreasing = TRUE)[T_edge1]

FDR_estc12 <- rep(0, length(T_edge1))
FDR_estc13 <- rep(0, length(T_edge1))
FDR_estc14 <- rep(0, length(T_edge1))
FDR_estc23 <- rep(0, length(T_edge1))
FDR_estc24 <- rep(0, length(T_edge1))
FDR_estc34 <- rep(0, length(T_edge1))


for(i in 1:length(threshold1)){
  thc12 <- thresholdc12[i]
  est_diffc12 <- as.numeric(pvijc12 >= thc12)
  FDR_estc12[i] <- sum((1-pvijc12) * est_diffc12)  / sum(est_diffc12)
  
  thc13 <- thresholdc13[i]
  est_diffc13 <- as.numeric(pvijc13 >= thc13)
  FDR_estc13[i] <- sum((1-pvijc13) * est_diffc13)  / sum(est_diffc13)
  
  thc14 <- thresholdc14[i]
  est_diffc14 <- as.numeric(pvijc14 >= thc14)
  FDR_estc14[i] <- sum((1-pvijc14) * est_diffc14)  / sum(est_diffc14)
  
  thc23 <- thresholdc23[i]
  est_diffc23 <- as.numeric(pvijc23 >= thc23)
  FDR_estc23[i] <- sum((1-pvijc23) * est_diffc23)  / sum(est_diffc23)
  
  thc24 <- thresholdc24[i]
  est_diffc24 <- as.numeric(pvijc24 >= thc24)
  FDR_estc24[i] <- sum((1-pvijc24) * est_diffc24)  / sum(est_diffc24)
  
  thc34 <- thresholdc34[i]
  est_diffc34 <- as.numeric(pvijc34 >= thc34)
  FDR_estc34[i] <- sum((1-pvijc34) * est_diffc34)  / sum(est_diffc34)
  
}
##### Use the same threshold as above and identify shared boundaries ######
alpha_nc = 0.1
T12c = sum(FDR_estc12<=alpha_nc)
T13c = sum(FDR_estc13<=alpha_nc)
T14c = sum(FDR_estc14<=alpha_nc)
T23c = sum(FDR_estc23<=alpha_nc)
T24c = sum(FDR_estc24<=alpha_nc)
T34c = sum(FDR_estc34<=alpha_nc)
Tallc = sum(FDR_estcall<=alpha_nc)


est_diffc12 <- as.numeric(pvijc12 >= thresholdc12[T12c])
neighbor_list_diffc12 <- neighbor_list0[est_diffc12 == 1, ]

est_diffc13 <- as.numeric(pvijc13 >= thresholdc13[T13c])
neighbor_list_diffc13 <- neighbor_list0[est_diffc13 == 1, ]

est_diffc14 <- as.numeric(pvijc14 >= thresholdc14[T14c])
neighbor_list_diffc14 <- neighbor_list0[est_diffc14 == 1, ]

est_diffc23 <- as.numeric(pvijc23 >= thresholdc23[T23c])
neighbor_list_diffc23 <- neighbor_list0[est_diffc23 == 1, ]

est_diffc24 <- as.numeric(pvijc24 >= thresholdc24[T24c])
neighbor_list_diffc24 <- neighbor_list0[est_diffc24 == 1, ]

est_diffc34 <- as.numeric(pvijc34 >= thresholdc34[T34c])
neighbor_list_diffc34 <- neighbor_list0[est_diffc34 == 1, ]

edge_plotc12 <- ggplot() +
  geom_polygon(data = ca.poly.df,  aes(long, lat, group = group),
               color = "black",
               fill = "white"
  )
for(i in which(est_diffc12==1)){
  edge_plotc12 = edge_plotc12 + geom_path(aes_string(x = path[[i]][,1], y = path[[i]][,2]), color = "red") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(size=14, face="bold",hjust = 0.5)) + 
    ggtitle(paste("Lung, Esophageal (T = ", T12c, ")", sep=""))
}
edge_plotc12


edge_plotc13 <- ggplot() +
  geom_polygon(data = ca.poly.df,  aes(long, lat, group = group),
               color = "black",
               fill = "white"
  )
for(i in which(est_diffc13==1)){
  edge_plotc13 = edge_plotc13 + geom_path(aes_string(x = path[[i]][,1], y = path[[i]][,2]), color = "red") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(size=14, face="bold",hjust = 0.5)) + 
    ggtitle(paste("Lung, Layrnx (T = ", T13c, ")", sep=""))
}
edge_plotc13


edge_plotc14 <- ggplot() +
  geom_polygon(data = ca.poly.df,  aes(long, lat, group = group),
               color = "black",
               fill = "white"
  )
for(i in which(est_diffc14==1)){
  edge_plotc14 = edge_plotc14 + geom_path(aes_string(x = path[[i]][,1], y = path[[i]][,2]), color = "red") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(size=14, face="bold",hjust = 0.5)) + 
    ggtitle(paste("Lung, Colorectal (T = ", T14c, ")", sep=""))
}
edge_plotc14

edge_plotc23 <- ggplot() +
  geom_polygon(data = ca.poly.df,  aes(long, lat, group = group),
               color = "black",
               fill = "white"
  )
for(i in which(est_diffc23==1)){
  edge_plotc23 = edge_plotc23 + geom_path(aes_string(x = path[[i]][,1], y = path[[i]][,2]), color = "red") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(size=14, face="bold",hjust = 0.5)) + 
    ggtitle(paste("Esophageal, Layrnx (T = ", T23c, ")", sep=""))
}
edge_plotc23


edge_plotc24 <- ggplot() +
  geom_polygon(data = ca.poly.df,  aes(long, lat, group = group),
               color = "black",
               fill = "white"
  )
for(i in which(est_diffc24==1)){
  edge_plotc24 = edge_plotc24 + geom_path(aes_string(x = path[[i]][,1], y = path[[i]][,2]), color = "red") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(size=14, face="bold",hjust = 0.5)) + 
    ggtitle(paste("Esophageal, Colorectal (T = ", T24c, ")", sep=""))
}
edge_plotc24


edge_plotc34 <- ggplot() +
  geom_polygon(data = ca.poly.df,  aes(long, lat, group = group),
               color = "black",
               fill = "white"
  )
for(i in which(est_diffc34==1)){
  edge_plotc34 = edge_plotc34 + geom_path(aes_string(x = path[[i]][,1], y = path[[i]][,2]), color = "red") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(size=14, face="bold",hjust = 0.5)) + 
    ggtitle(paste("Larynx, Colorectal (T = ", T34c, ")", sep=""))
}
edge_plotc34

#pdf("est_share_car_ind_pois.pdf", height = 10, width = 15)
ggarrange(edge_plotc12, edge_plotc13, edge_plotc14, 
          edge_plotc23, edge_plotc24, edge_plotc34, nrow = 2, ncol = 3)
#dev.off()

# model fitting: D score
Y_rep1 <- matrix(0, nrow = n, ncol = 10000)
Y_rep2 <- matrix(0, nrow = n, ncol = 10000)
Y_rep3 <- matrix(0, nrow = n, ncol = 10000)
Y_rep4 <- matrix(0, nrow = n, ncol = 10000)
for(i in 1:10000){
  
  beta1 <- mcmc_samples$BUGSoutput$sims.matrix[i,1]
  beta2 <- mcmc_samples$BUGSoutput$sims.matrix[i,2]
  beta3 <- mcmc_samples$BUGSoutput$sims.matrix[i,3]
  beta4 <- mcmc_samples$BUGSoutput$sims.matrix[i,4]
  
  W1_mcmc = phis[i,1:58]
  W2_mcmc = phis[i,59:116]
  W3_mcmc = phis[i,117:174]
  W4_mcmc = phis[i,175:232]
  
  #set.seed(seed)
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

D_CARind = c(D.latent1, D.latent2, D.latent3, D.latent4)
