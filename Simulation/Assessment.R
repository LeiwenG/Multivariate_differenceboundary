library(dplyr)
library(caret)
library(Matrix)
library(LaplacesDemon)
library(matrixStats)
library(monomvn)

###### MDAGAR #####
load("simBARDP_dagarc25.RData")

# "disease difference" in the same region
disease_diff <- as.numeric(phi_true1 != phi_true2)
sum(disease_diff)

# estimate for "disease difference"
T_edged = c(10, 15, 20, 22, 25, 30)
specd <- matrix(0, 50, length(T_edged))
sensd <- matrix(0, 50, length(T_edged))
for(seed in 1:50){
  print(seed)

  phis <-  mcmc_samples[[seed]]$phi
  phis1 <- phis[,1:58][,order(final_perm)]
  phis2 <- phis[,59:116][,order(final_perm)]
  phis_origin <- cbind(phis1, phis2)
  
  disease_samples <- t(apply(phis_origin, 1, function(x){
    #x <- phis_origin[1,]
    diff_d <- NULL
    for(i in 1:ncol(phis1)){
      if(x[i] != x[i+ncol(phis1)]){
        diff_d <- c(diff_d, 1)
      }else{
        diff_d <- c(diff_d, 0)
      }
    }
    return(diff_d)
  }))
  
  pd <- apply(disease_samples, 2, mean)
  threshold_d = sort(pd, decreasing = TRUE)[T_edged]
  
  for(i in 1:length(T_edged)){
    est_diffd <- factor(as.numeric(pd >= threshold_d[i]), levels = c(0,1))
    conf_matrixd <-table(est_diffd, disease_diff)
    specd[seed,i] <- sensitivity(conf_matrixd)
    sensd[seed,i] <- specificity(conf_matrixd)
  }
}

table_d <- cbind(colMeans(specd), colMeans(sensd))

# specificity and specificity for difference boundaries
spec1_mean_dagar = colMeans(spec1)
sens1_mean_dagar = colMeans(sens1)

spec2_mean_dagar = colMeans(spec2)
sens2_mean_dagar = colMeans(sens2)

spec12_mean_dagar = colMeans(spec12)
sens12_mean_dagar = colMeans(sens12)

spec21_mean_dagar = colMeans(spec21)
sens21_mean_dagar = colMeans(sens21)

table = cbind(spec1_mean_dagar, sens1_mean_dagar, spec2_mean_dagar, sens2_mean_dagar,
              spec12_mean_dagar, sens12_mean_dagar, spec21_mean_dagar, sens21_mean_dagar)

# model fitting
WAIC <- rep(0,50)
D_dagar <- rep(0,50)
KL <- matrix(0, nrow = 50, ncol = 10000)
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
  m1 <- X1 %*% beta1 + (phi_true1 - mean(phi_true1))
  m2 <- X2 %*% beta2 + (phi_true2 - mean(phi_true2))
  y1 = m1 + sqrt(sigma_esq1) * rnorm(n_county)
  y2 = m2 + sqrt(sigma_esq2) * rnorm(n_county)

  Xo1 <- X1[final_perm,]
  yo1 <- y1[final_perm]
  Xo2 <- X2[final_perm,]
  yo2 <- y2[final_perm]
  
  Y1 = c(yo1, yo2)
  Xb1 = as.matrix(bdiag(Xo1, Xo2))
  
  Y2 = c(yo2, yo1)
  Xb2 = as.matrix(bdiag(Xo2, Xo1))
  sample.mcmc1 <- mcmc_samples[[seed]]
  sample.mcmc <- cbind(sample.mcmc1$beta, sample.mcmc1$phi, sample.mcmc1$taue1, sample.mcmc1$taue2)
  
  #WAIC
  PL_single <- matrix(0, nrow = 2*n, ncol = 10000)
  for(i in 1:10000){
    theta <- sample.mcmc[i,]
    
    beta1 <- theta[1:2]
    beta2 <- theta[3:4]
    sigmasq1 <- 1/theta[121]
    sigmasq2 <- 1/theta[122]
    
    W1_mcmc = theta[5:62]
    W2_mcmc = theta[63:120]
    
    for(j in 1:n){
      v1 <- sigmasq1
      PL_single[j,i] <- as.numeric(-1/2*log(2*pi) - 1/2*log(v1) - 1/(2*v1)*(yo1[j]-Xo1[j,]%*%as.vector(beta1) - W1_mcmc[j])^2)
      v2 <- sigmasq2
      PL_single[(j+n),i] <- as.numeric(-1/2*log(2*pi) - 1/2*log(v2) - 1/(2*v2)*(yo2[j]-Xo2[j,]%*%as.vector(beta2) - W2_mcmc[j])^2)
    }
  }
  WAIC[seed] <- WAIC(PL_single)$WAIC
  
  #D
  Y_rep1 <- matrix(0, nrow = n, ncol = 10000)
  Y_rep2 <- matrix(0, nrow = n, ncol = 10000)
  for(i in 1:10000){
    theta <- sample.mcmc[i,]
    beta1 <- theta[1:2]
    beta2 <- theta[3:4]
    sigmasq1 <- 1/theta[121]
    sigmasq2 <- 1/theta[122]
    
    W1_mcmc = theta[5:62]
    W2_mcmc = theta[63:120]
    
    set.seed(seed)
    z1 <- rnorm(n, 0, 1)
    z2 <- rnorm(n, 0, 1)
    Y_rep1[,i] <- Xo1 %*% as.vector(beta1) +  W1_mcmc + sqrt(as.numeric(sigmasq1)) * z1
    Y_rep2[,i] <- Xo2 %*% as.vector(beta2) +  W2_mcmc + sqrt(as.numeric(sigmasq2)) * z2
    
  }
  mu_rep1 = rowMeans(Y_rep1)
  mu_rep2 = rowMeans(Y_rep2)
  
  var_rep1 = rowVars(Y_rep1)
  var_rep2 = rowVars(Y_rep2)
  
  G.latent1 = sum((yo1 - mu_rep1)^2)
  P.latent1 = sum(var_rep1)
  D.latent1 = G.latent1 + P.latent1
  
  G.latent2 = sum((yo2 - mu_rep2)^2)
  P.latent2 = sum(var_rep2)
  D.latent2 = G.latent2 + P.latent2
  
  D_dagar[seed] = D.latent1 + D.latent2
 
  #KL
  KLi = NULL
  V = as.matrix(bdiag(sigma_esq1*diag(n), sigma_esq2*diag(n)))
  for(i in 1:10000){
    theta <- sample.mcmc[i,]
    
    beta1 <- theta[1:2]
    beta2 <- theta[3:4]
    sigmasq1 <- 1/theta[121]
    sigmasq2 <- 1/theta[122]
    
    W1_mcmc = theta[5:62]
    W2_mcmc = theta[63:120]
    
    m_est1 <- Xo1 %*% as.vector(beta1) +  W1_mcmc 
    m_est2 <- Xo2 %*% as.vector(beta2) +  W2_mcmc 
    V_est = as.matrix(bdiag(sigmasq1*diag(n), sigmasq2*diag(n)))
    KLi[i] = kl.norm(mu1 = c(m1[final_perm], m2[final_perm]), S1= V, mu2 = c(m_est1, m_est2), S2= V_est)
  }
  
  KL[seed,] = KLi
}

saveRDS(WAIC,"WAIC_ARDP_Jdagarc25.rds")
saveRDS(D_dagar,"D_ARDP_Jdagarc25.rds")
saveRDS(KL,"KL_ARDP_Jdagarc25.rds")

###### MCAR #####
load("simBARDP_carc25_1.RData")

# "disease difference" in the same region
disease_diff <- as.numeric(phi_true1 != phi_true2)
sum(disease_diff)

# estimate for "disease difference"
T_edged = c(10, 15, 20, 22, 25, 30)
specd <- matrix(0, 30, length(T_edged))
sensd <- matrix(0, 30, length(T_edged))
for(seed in 1:50){
  print(seed)
  
  phis <-  mcmc_samples[[seed]]$phi
  phis1 <- phis[,1:58][,order(final_perm)]
  phis2 <- phis[,59:116][,order(final_perm)]
  phis_origin <- cbind(phis1, phis2)

  disease_samples <- t(apply(phis_origin, 1, function(x){
    #x <- phis_origin[1,]
    diff_d <- NULL
    for(i in 1:ncol(phis1)){
      if(x[i] != x[i+ncol(phis1)]){
        diff_d <- c(diff_d, 1)
      }else{
        diff_d <- c(diff_d, 0)
      }
    }
    return(diff_d)
  }))
  
  pd <- apply(disease_samples, 2, mean)
  threshold_d = sort(pd, decreasing = TRUE)[T_edged]
  
  for(i in 1:length(T_edged)){
    est_diffd <- factor(as.numeric(pd >= threshold_d[i]), levels = c(0,1))
    conf_matrixd <-table(est_diffd, disease_diff)
    specd[seed,i] <- sensitivity(conf_matrixd)
    sensd[seed,i] <- specificity(conf_matrixd)
  }
}

table_d <- cbind(colMeans(specd), colMeans(sensd))

# specificity and specificity for difference boundaries
spec1_mean_car = colMeans(spec1)
sens1_mean_car = colMeans(sens1)

spec2_mean_car = colMeans(spec2)
sens2_mean_car = colMeans(sens2)

spec12_mean_car = colMeans(spec12)
sens12_mean_car = colMeans(sens12)

spec21_mean_car = colMeans(spec21)
sens21_mean_car = colMeans(sens21)

table = cbind(spec1_mean_car, sens1_mean_car, spec2_mean_car, sens2_mean_car,
              spec12_mean_car, sens12_mean_car, spec21_mean_car, sens21_mean_car)

# model fitting
WAIC_car <- rep(0,50)
D_car <- rep(0,50)
Dn_car <- rep(0,50)
KL_car <- matrix(0, nrow = 50, ncol = 10000)
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
  m1 <- X1 %*% beta1 + (phi_true1 - mean(phi_true1))
  m2 <- X2 %*% beta2 + (phi_true2 - mean(phi_true2))
  y1 = m1 + sqrt(sigma_esq1) * rnorm(n_county)
  y2 = m2 + sqrt(sigma_esq2) * rnorm(n_county)

  Xo1 <- X1[final_perm,]
  yo1 <- y1[final_perm]
  Xo2 <- X2[final_perm,]
  yo2 <- y2[final_perm]
  
  
  Y1 = c(yo1, yo2)
  Xb1 = as.matrix(bdiag(Xo1, Xo2))
  
  Y2 = c(yo2, yo1)
  Xb2 = as.matrix(bdiag(Xo2, Xo1))
  sample.mcmc1 <- mcmc_samples[[seed]]
  sample.mcmc <- cbind(sample.mcmc1$beta, sample.mcmc1$phi, sample.mcmc1$taue1, sample.mcmc1$taue2)
  #WAIC
  PL_single <- matrix(0, nrow = 2*n, ncol = 10000)
  for(i in 1:10000){
    theta <- sample.mcmc[i,]
    
    beta1 <- theta[1:2]
    beta2 <- theta[3:4]
    sigmasq1 <- 1/theta[121]
    sigmasq2 <- 1/theta[122]
    
    W1_mcmc = theta[5:62]
    W2_mcmc = theta[63:120]
    
    for(j in 1:n){
      v1 <- sigmasq1
      PL_single[j,i] <- as.numeric(-1/2*log(2*pi) - 1/2*log(v1) - 1/(2*v1)*(yo1[j]-Xo1[j,]%*%as.vector(beta1) - W1_mcmc[j])^2)
      v2 <- sigmasq2
      PL_single[(j+n),i] <- as.numeric(-1/2*log(2*pi) - 1/2*log(v2) - 1/(2*v2)*(yo2[j]-Xo2[j,]%*%as.vector(beta2) - W2_mcmc[j])^2)
    }
  }
  WAIC_car[seed] <- WAIC(PL_single)$WAIC
  
  #D
  Y_rep1 <- matrix(0, nrow = n, ncol = 10000)
  Y_rep2 <- matrix(0, nrow = n, ncol = 10000)
  for(i in 1:10000){
    theta <- sample.mcmc[i,]
    
    beta1 <- theta[1:2]
    beta2 <- theta[3:4]
    sigmasq1 <- 1/theta[121]
    sigmasq2 <- 1/theta[122]
    
    W1_mcmc = theta[5:62]
    W2_mcmc = theta[63:120]
    
    set.seed(seed)
    z1 <- rnorm(n, 0, 1)
    z2 <- rnorm(n, 0, 1)
    Y_rep1[,i] <- Xo1 %*% as.vector(beta1) +  W1_mcmc + sqrt(as.numeric(sigmasq1)) * z1
    Y_rep2[,i] <- Xo2 %*% as.vector(beta2) +  W2_mcmc + sqrt(as.numeric(sigmasq2)) * z2
    
  }
  mu_rep1 = rowMeans(Y_rep1)
  mu_rep2 = rowMeans(Y_rep2)
  
  var_rep1 = rowVars(Y_rep1)
  var_rep2 = rowVars(Y_rep2)
  
  G.latent1 = sum((yo1 - mu_rep1)^2)
  P.latent1 = sum(var_rep1)
  D.latent1 = G.latent1 + P.latent1
  
  G.latent2 = sum((yo2 - mu_rep2)^2)
  P.latent2 = sum(var_rep2)
  D.latent2 = G.latent2 + P.latent2
  
  D_car[seed] = D.latent1 + D.latent2
  
  #KL
  KLi = NULL
  V = as.matrix(bdiag(sigma_esq1*diag(n), sigma_esq2*diag(n)))
  for(i in 1:10000){
    theta <- sample.mcmc[i,]
    
    beta1 <- theta[1:2]
    beta2 <- theta[3:4]
    sigmasq1 <- 1/theta[121]
    sigmasq2 <- 1/theta[122]
    
    W1_mcmc = theta[5:62]
    W2_mcmc = theta[63:120]
    
    m_est1 <- Xo1 %*% as.vector(beta1) +  W1_mcmc 
    m_est2 <- Xo2 %*% as.vector(beta2) +  W2_mcmc 
    V_est = as.matrix(bdiag(sigmasq1*diag(n), sigmasq2*diag(n)))
    KLi[i] = kl.norm(mu1 = c(m1[final_perm], m2[final_perm]), S1= V, mu2 = c(m_est1, m_est2), S2= V_est)
  }
  
  KL_car[seed,] = KLi
  
}

saveRDS(WAIC_car,"WAIC_ARDP_Jcarc25_1.rds")
saveRDS(D_car,"D_ARDP_Jcarc25_1.rds")
saveRDS(KL_car,"KL_ARDP_Jcarc25_1.rds")


###### DAGAR_ind #####
load("simBARDP_dagar_indc25_v1.RData")

# "disease difference" in the same region
disease_diff <- as.numeric(phi_true1 != phi_true2)
sum(disease_diff)

# estimate for "disease difference"
T_edged = c(10, 15, 20, 22, 25, 30)
specd <- matrix(0, 50, length(T_edged))
sensd <- matrix(0, 50, length(T_edged))
for(seed in 1:50){
  print(seed)
  #seed = 1
  phis <-  mcmc_samples[[seed]][,6:121]
  phis1 <- phis[,1:58][,order(final_perm)]
  phis2 <- phis[,59:116][,order(final_perm)]
  phis_origin <- cbind(phis1, phis2)
  table(apply(phis2, 1, function(x) length(unique(x))))
  table(apply(phis1, 1, function(x) length(unique(x))))
  
  disease_samples <- t(apply(phis_origin, 1, function(x){
    #x <- phis_origin[1,]
    diff_d <- NULL
    for(i in 1:ncol(phis1)){
      if(x[i] != x[i+ncol(phis1)]){
        diff_d <- c(diff_d, 1)
      }else{
        diff_d <- c(diff_d, 0)
      }
    }
    return(diff_d)
  }))
  
  pd <- apply(disease_samples, 2, mean)
  threshold_d = sort(pd, decreasing = TRUE)[T_edged]
  
  for(i in 1:length(T_edged)){
    est_diffd <- as.numeric(pd >= threshold_d[i])
    conf_matrixd <-table(est_diffd, disease_diff)
    specd[seed,i] <- sensitivity(conf_matrixd)
    sensd[seed,i] <- specificity(conf_matrixd)
  }
}

table_d <- cbind(colMeans(specd), colMeans(sensd))


# specificity and specificity for difference boundaries
spec1_mean_dagar = colMeans(spec1)
sens1_mean_dagar = colMeans(sens1)

spec2_mean_dagar = colMeans(spec2)
sens2_mean_dagar = colMeans(sens2)

spec12_mean_dagar = colMeans(spec12)
sens12_mean_dagar = colMeans(sens12)

spec21_mean_dagar = colMeans(spec21)
sens21_mean_dagar = colMeans(sens21)

table = cbind(spec1_mean_dagar, sens1_mean_dagar, spec2_mean_dagar, sens2_mean_dagar,
              spec12_mean_dagar, sens12_mean_dagar, spec21_mean_dagar, sens21_mean_dagar)

# model fitting
WAIC <- rep(0,50)
D_dagar <- rep(0,50)
KL <- matrix(0, nrow = 50, ncol = 10000)
for(seed in 1:50){
  print(seed)
 
  set.seed(seed)
  X1=cbind(1, rnorm(n_county))
  X2=cbind(1, rnorm(n_county))
  beta1=c(2,5)
  beta2=c(1,6)
  sigma_esq1 = 0.1
  sigma_esq2 = 0.1
  m1 <- X1 %*% beta1 + (phi_true1 - mean(phi_true1))
  m2 <- X2 %*% beta2 + (phi_true2 - mean(phi_true2))
  y1 = m1 + sqrt(sigma_esq1) * rnorm(n_county)
  y2 = m2 + sqrt(sigma_esq2) * rnorm(n_county)
  Xo1 <- X1[final_perm,]
  yo1 <- y1[final_perm]
  Xo2 <- X2[final_perm,]
  yo2 <- y2[final_perm]
  
  
  Y1 = c(yo1, yo2)
  Xb1 = as.matrix(bdiag(Xo1, Xo2))
  
  sample.mcmc <- mcmc_samples[[seed]]
  #WAIC
  PL_single <- matrix(0, nrow = 2*n, ncol = 10000)
  for(i in 1:10000){
    theta <- sample.mcmc[i,]
    
    beta1 <- theta[1:2]
    beta2 <- theta[3:4]
    sigmasq1 <- theta[1867]
    sigmasq2 <- theta[1868]
    
    W1_mcmc = theta[6:63]
    W2_mcmc = theta[64:121]
    
    for(j in 1:n){
      v1 <- sigmasq1
      PL_single[j,i] <- as.numeric(-1/2*log(2*pi) - 1/2*log(v1) - 1/(2*v1)*(yo1[j]-Xo1[j,]%*%as.vector(beta1) - W1_mcmc[j])^2)
      v2 <- sigmasq2
      PL_single[(j+n),i] <- as.numeric(-1/2*log(2*pi) - 1/2*log(v2) - 1/(2*v2)*(yo2[j]-Xo2[j,]%*%as.vector(beta2) - W2_mcmc[j])^2)
    }
  }
  WAIC[seed] <- WAIC(PL_single)$WAIC
}
#D
Y_rep1 <- matrix(0, nrow = n, ncol = 10000)
Y_rep2 <- matrix(0, nrow = n, ncol = 10000)
for(i in 1:10000){
  theta <- sample.mcmc[i,]
  beta1 <- theta[1:2]
  beta2 <- theta[3:4]
  sigmasq1 <- theta[1867]
  sigmasq2 <- theta[1868]
  
  W1_mcmc = theta[6:63]
  W2_mcmc = theta[64:121]
  
  set.seed(seed)
  z1 <- rnorm(n, 0, 1)
  z2 <- rnorm(n, 0, 1)
  Y_rep1[,i] <- Xo1 %*% as.vector(beta1) +  W1_mcmc + sqrt(as.numeric(sigmasq1)) * z1
  Y_rep2[,i] <- Xo2 %*% as.vector(beta2) +  W2_mcmc + sqrt(as.numeric(sigmasq2)) * z2
  
}
mu_rep1 = rowMeans(Y_rep1)
mu_rep2 = rowMeans(Y_rep2)

var_rep1 = rowVars(Y_rep1)
var_rep2 = rowVars(Y_rep2)

G.latent1 = sum((yo1 - mu_rep1)^2)
P.latent1 = sum(var_rep1)
D.latent1 = G.latent1 + P.latent1

G.latent2 = sum((yo2 - mu_rep2)^2)
P.latent2 = sum(var_rep2)
D.latent2 = G.latent2 + P.latent2

D_dagar[seed] = D.latent1 + D.latent2

#KL
KLi = NULL
V = as.matrix(bdiag(sigma_esq1*diag(n), sigma_esq2*diag(n)))
for(i in 1:10000){
  theta <- sample.mcmc[i,]
  
  beta1 <- theta[1:2]
  beta2 <- theta[3:4]
  sigmasq1 <- theta[1867]
  sigmasq2 <- theta[1868]
  
  W1_mcmc = theta[6:63]
  W2_mcmc = theta[64:121]
  
  m_est1 <- Xo1 %*% as.vector(beta1) +  W1_mcmc 
  m_est2 <- Xo2 %*% as.vector(beta2) +  W2_mcmc 
  V_est = as.matrix(bdiag(sigmasq1*diag(n), sigmasq2*diag(n)))
  KLi[i] = kl.norm(mu1 = c(m1[final_perm], m2[final_perm]), S1= V, mu2 = c(m_est1, m_est2), S2= V_est)
}

KL[seed,] = KLi

saveRDS(WAIC,"WAIC_ARDP_dagar_indc25_v1.rds")
saveRDS(D_dagar,"D_ARDP_dagar_indc25_v1.rds")
saveRDS(KL,"KL_ARDP_dagar_indc25_v1.rds")


###### CAR_ind #####
load("simBARDP_car_indc25.RData")

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

disease_diff <- as.numeric(phi_true1 != phi_true2)
sum(disease_diff)

true_diffm <- rbind(true_diff1, true_diff2, true_diff3, true_diff4)

# disease difference
T_edged = c(10, 15, 20,22, 25, 30)
specd <- matrix(0, 50, length(T_edged))
sensd <- matrix(0, 50, length(T_edged))
for(seed in 1:50){
  print(seed)
  #seed = 1
  phis <-  mcmc_samples[[seed]][,6:121]
  phis1 <- phis[,1:58][,order(final_perm)]
  phis2 <- phis[,59:116][,order(final_perm)]
  phis_origin <- cbind(phis1, phis2)
  table(apply(phis2, 1, function(x) length(unique(x))))
  table(apply(phis1, 1, function(x) length(unique(x))))
  
  disease_samples <- t(apply(phis_origin, 1, function(x){
    diff_d <- NULL
    for(i in 1:ncol(phis1)){
      if(x[i] != x[i+ncol(phis1)]){
        diff_d <- c(diff_d, 1)
      }else{
        diff_d <- c(diff_d, 0)
      }
    }
    return(diff_d)
  }))
  
  pd <- apply(disease_samples, 2, mean)
  threshold_d = sort(pd, decreasing = TRUE)[T_edged]
  
  for(i in 1:length(T_edged)){
    est_diffd <- as.numeric(pd >= threshold_d[i])
    conf_matrixd <-table(est_diffd, disease_diff)
    specd[seed,i] <- sensitivity(conf_matrixd)
    sensd[seed,i] <- specificity(conf_matrixd)
  }
}

table_d <- cbind(colMeans(specd), colMeans(sensd))


# specificity and specificity for difference boundaries
spec1_mean_car = colMeans(spec1)
sens1_mean_car = colMeans(sens1)

spec2_mean_car = colMeans(spec2)
sens2_mean_car = colMeans(sens2)

spec12_mean_car = colMeans(spec12)
sens12_mean_car = colMeans(sens12)

spec21_mean_car = colMeans(spec21)
sens21_mean_car = colMeans(sens21)

table = cbind(spec1_mean_car, sens1_mean_car, spec2_mean_car, sens2_mean_car,
              spec12_mean_car, sens12_mean_car, spec21_mean_car, sens21_mean_car)

# model fitting
WAIC_car <- rep(0,50)
D_car <- rep(0,50)
KL_car <- matrix(0, nrow = 50, ncol = 10000)
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
  m1 <- X1 %*% beta1 + (phi_true1 - mean(phi_true1))
  m2 <- X2 %*% beta2 + (phi_true2 - mean(phi_true2))
  y1 = m1 + sqrt(sigma_esq1) * rnorm(n_county)
  y2 = m2 + sqrt(sigma_esq2) * rnorm(n_county)
  #y = X %*% beta + phi_true  + sqrt(sigma_e2) * rnorm(n_county)
  Xo1 <- X1[final_perm,]
  yo1 <- y1[final_perm]
  Xo2 <- X2[final_perm,]
  yo2 <- y2[final_perm]
  
  
  Y1 = c(yo1, yo2)
  Xb1 = as.matrix(bdiag(Xo1, Xo2))
  
  sample.mcmc <- mcmc_samples[[seed]]
  #WAIC
  PL_single <- matrix(0, nrow = 2*n, ncol = 10000)
  for(i in 1:10000){
    theta <- sample.mcmc[i,]
    
    beta1 <- theta[1:2]
    beta2 <- theta[3:4]
    sigmasq1 <- theta[1867]
    sigmasq2 <- theta[1868]
    
    W1_mcmc = theta[6:63]
    W2_mcmc = theta[64:121]
    
    for(j in 1:n){
      v1 <- sigmasq1
      PL_single[j,i] <- as.numeric(-1/2*log(2*pi) - 1/2*log(v1) - 1/(2*v1)*(yo1[j]-Xo1[j,]%*%as.vector(beta1) - W1_mcmc[j])^2)
      v2 <- sigmasq2
      PL_single[(j+n),i] <- as.numeric(-1/2*log(2*pi) - 1/2*log(v2) - 1/(2*v2)*(yo2[j]-Xo2[j,]%*%as.vector(beta2) - W2_mcmc[j])^2)
    }
  }
  WAIC_car[seed] <- WAIC(PL_single)$WAIC
  
  #D
  Y_rep1 <- matrix(0, nrow = n, ncol = 10000)
  Y_rep2 <- matrix(0, nrow = n, ncol = 10000)
  for(i in 1:10000){
    theta <- sample.mcmc[i,]
    
    beta1 <- theta[1:2]
    beta2 <- theta[3:4]
    sigmasq1 <- theta[1867]
    sigmasq2 <- theta[1868]
    
    W1_mcmc = theta[6:63]
    W2_mcmc = theta[64:121]
    
    set.seed(seed)
    z1 <- rnorm(n, 0, 1)
    z2 <- rnorm(n, 0, 1)
    Y_rep1[,i] <- Xo1 %*% as.vector(beta1) +  W1_mcmc + sqrt(as.numeric(sigmasq1)) * z1
    Y_rep2[,i] <- Xo2 %*% as.vector(beta2) +  W2_mcmc + sqrt(as.numeric(sigmasq2)) * z2
    
  }
  mu_rep1 = rowMeans(Y_rep1)
  mu_rep2 = rowMeans(Y_rep2)
  
  var_rep1 = rowVars(Y_rep1)
  var_rep2 = rowVars(Y_rep2)
  
  G.latent1 = sum((yo1 - mu_rep1)^2)
  P.latent1 = sum(var_rep1)
  D.latent1 = G.latent1 + P.latent1
  
  G.latent2 = sum((yo2 - mu_rep2)^2)
  P.latent2 = sum(var_rep2)
  D.latent2 = G.latent2 + P.latent2
  
  D_car[seed] = D.latent1 + D.latent2
  
  #KL
  KLi = NULL
  V = as.matrix(bdiag(sigma_esq1*diag(n), sigma_esq2*diag(n)))
  for(i in 1:10000){
    theta <- sample.mcmc[i,]
    
    beta1 <- theta[1:2]
    beta2 <- theta[3:4]
    sigmasq1 <- theta[1867]
    sigmasq2 <- theta[1868]
    
    W1_mcmc = theta[6:63]
    W2_mcmc = theta[64:121]
    
    m_est1 <- Xo1 %*% as.vector(beta1) +  W1_mcmc 
    m_est2 <- Xo2 %*% as.vector(beta2) +  W2_mcmc 
    V_est = as.matrix(bdiag(sigmasq1*diag(n), sigmasq2*diag(n)))
    KLi[i] = kl.norm(mu1 = c(m1[final_perm], m2[final_perm]), S1= V, mu2 = c(m_est1, m_est2), S2= V_est)
  }
  
  KL_car[seed,] = KLi
}

saveRDS(WAIC_car,"WAIC_ARDP_car_indc25.rds")
saveRDS(D_car,"D_ARDP_car_indc25.rds")
saveRDS(KL_car,"KL_ARDP_car_indc25.rds")
