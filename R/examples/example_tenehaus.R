library(tidyverse)
# IMPORTANT: blocks MUST be a named list
source("R/tenenhaus_cv.R")

pp <- c(10,10,10);
p <- sum(pp)
s  <- c(1:5);
r <- 1;
max_iter <- 15000;
p<-sum(pp)

n=30

set.seed(2025)
p <- sum(pp)
pp1=pp[1];pp2=pp[2];pp3=pp[3]
u1 <- matrix(0, pp[1],r);
u2 <- matrix(0, pp[2],r);
u3 <- matrix(0, pp[3],r);

# Generate Sigma and Sigma0
Sigma = diag(sum(pp));
T1 = toeplitz(0.5^(0:(pp[1]-1)));
T2 = toeplitz(0.7^(0:(pp[2]-1)));
T3 = toeplitz(0.9^(0:(pp[3]-1)));
u1[s,1:r] <- matrix( rnorm(length(s)*r,mean=0,sd=1), length(s), r)
u1 <- u1 %*%(pracma::sqrtm(t(u1[s,1:r]) %*% T1[s,s] %*% u1[s,1:r])$Binv);
u2[s,1:r] = matrix(rnorm(length(s)*r,mean=0,sd=1), length(s), r) 
u2 <- u2 %*% (pracma::sqrtm(t(u2[s,1:r]) %*% T2[s,s] %*% u2[s,1:r])$Binv);
u3[s,1:r] = matrix( rnorm(length(s)*r,mean=0,sd=1), length(s), r) 
u3 <- u3 %*%(pracma::sqrtm(t(u3[s,1:r]) %*% T3[s,s] %*% u3[s,1:r])$Binv);

idx1 = 1:pp[1];
idx2 = (pp[1]+1):(pp[1]+pp[2]);
idx3 = (pp[1]+pp[2]+1):(pp[1]+pp[2]+pp[3]);
Sigma[idx1, idx1] = T1;
Sigma[idx2, idx2] = T2;
Sigma[idx3, idx3] = T3;


SigmaD = Sigma;

Sigma[idx1, idx2] <-T1 %*% u1 %*% t(u2) %*% T2;
Sigma[idx1, idx3] <-T1 %*% u1 %*% t(u3) %*% T3;
Sigma[idx2, idx3] <-T2 %*% u2 %*% t(u3) %*% T3;
Sigma = Sigma + t(Sigma) - SigmaD+diag(p);


Mask = matrix(0, sum(pp),sum(pp));
Mask[idx1,idx1] <- matrix(1,pp[1],pp[1]);
Mask[idx2,idx2] <- matrix(1,pp[2],pp[2]);
Mask[idx3,idx3] <- matrix(1,pp[3],pp[3]);
Sigma0 = Sigma * Mask;



heatmap(Sigma, Rowv=NA, Colv=NA)
heatmap(Sigma0, Rowv=NA, Colv=NA)
# Generate data for SGCA
X <-mvrnorm(n, rep(0, p) , Sigma)
S <- t(X)%*% X / n
sigma0hat <- S * Mask

C_0 <- solve(Sigma0) %*% Sigma %*% solve(Sigma0)
can_vec=pracma::sqrtm(Sigma0)$Binv %*% svd(pracma::sqrtm(Sigma0)$B %*% C_0 %*%pracma::sqrtm(Sigma0)$B)$u
a=can_vec[,1]

blocks=list()
blocks[[1]]=X[,idx1]
blocks[[2]]=X[,idx2]
blocks[[3]]=X[,idx3]
# colnames(blocks[[1]]) <- sapply(1:ncol(blocks[[1]]), function(x){paste0("V1_", x)})    
# colnames(blocks[[2]]) <- sapply(1:ncol(blocks[[2]]), function(x){paste0("V2_", x)}) 
# colnames(blocks[[3]]) <- sapply(1:ncol(blocks[[3]]), function(x){paste0("V3_", x)})     
lambda_values <- c(0, 0.001, 0.1, 0.25, 0.5, 0.75, 1)

cv_unsup <- rgcca_unsupervised_cv_tau(
  blocks = blocks,
  lambda_values = lambda_values,
  kfold = 5,
  n_cores = 4
)

cv_unsup$best_tau
fit_best <- cv_unsup$fit_full

