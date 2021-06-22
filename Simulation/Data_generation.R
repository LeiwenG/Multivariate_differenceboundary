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
library(mvtnorm)
library(msos)
library(hesim)
library(LaplacesDemon)
library(blockmatrix)

county.ID <- sapply(strsplit(ca.county$names, ","), function(x) x[2])
ca.poly = map2SpatialPolygons(ca.county, IDs=county.ID)
ca.coords = coordinates(ca.poly)
n_county <- length(county.ID)


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
projmat=cbind(ca.albersproj$x,ca.albersproj$y)

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
projmat=projmat[final_perm,]
dmat=as.matrix(dist(projmat))
dmat=dmat/mean(dmat[which(Minc==1)])
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

### Specify Covariance matrix for spatial components
q=2
rho = c(0.2, 0.8)

G1 = rho[1]^dmat
G2 = rho[2]^dmat

a11 = 1
a12 = 1
a22 = 1
A = matrix(c(a11, a12, 0, a22), 2, 2)
kprod = kronecker(A, diag(n))
invQ = as.matrix(bdiag(G1, G2))
Vr = as.matrix(kprod %*% invQ %*% t(kprod))
Cov = rep(0, n)
for(i in 1:n){
  Cov[i] <- Vr[i, i+n]/(sqrt(Vr[i,i])*sqrt(Vr[n+i, n+i]))
}

alpha = 1
sigmas_sq = 4

### Functions to generate discrete random effects

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

K=15
phi_true1 = list()
phi_true2 = list()
nq = n*q

### Generate discrete random effects
seed=25
set.seed(seed)
r = rmvnorm(1, rep(0, nq), Vr)
v = rbeta(K,1,alpha)
while(v[1] > 0.8){v = rbeta(K,1,alpha)}
probvec = makeprobs(v)
F_r=pnorm(r,0,sqrt(diag(Vr)))
u=makeu(F_r,probvec)

thetavec = rnorm(K, 0, sqrt(sigmas_sq))
phivec = thetavec[u]

phi_true1 = phivec[1:n][order(final_perm)]
phi_true2 = phivec[(n+1):(2*n)][order(final_perm)]

saveRDS(phi_true1, "phi_true1_25.rds")
saveRDS(phi_true2, "phi_true2_25.rds")

ca.poly$phi_true1 <- phi_true1
ca.poly$phi_true2 <- phi_true2

values <- sort(unique(phi_true1))
color.pallete = brewer.pal(length(values),"Blues")
col1 <- rep(0, length(phi_true1))
for(i in 1:length(values)){
  col1[phi_true1 == values[i]] <- color.pallete[i]
}

col2 <- rep(0, length(phi_true2))
for(i in 1:length(values)){
  col2[phi_true2 == values[i]] <- color.pallete[i]
}


pdf("sim_boundary2.pdf", height = 5, width = 10)
par(mfrow=c(1,2))
plot(ca.poly, col = col1)
leg.txt1 = c("-2.67", "-1.73", "-0.98","0.42","0.77")
#legend("bottomleft", title="", legend=leg.txt1, cex=1.25, bty="n", horiz = FALSE, 
#      fill = color.pallete)
plot(ca.poly, col = col2)
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("center", legend=leg.txt1, xpd = TRUE, cex=1.25, bty="n", horiz = FALSE, 
       fill = color.pallete)
dev.off()



