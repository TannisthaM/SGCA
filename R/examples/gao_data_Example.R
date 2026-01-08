library(MASS)
library(stats)
library(pracma)
library(tidyverse)
library(ggplot2)
library(Matrix)
library(tidyr)
library(geigen)
library(dplyr)
library(expm)
library(foreach)
library(doParallel)

pp <- c(10,10,10);
p <- sum(pp)
s  <- c(1:5);
r <- 1;
k <- 20;
eta <- 0.001;
max_iter <- 15000;
p<-sum(pp)
rho<-1
lambda <- 0.01;

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

# Generate data for SGCA
X <-mvrnorm(n, rep(0, p) , Sigma)
S <- t(X)%*% X / n
sigma0hat <- S * Mask

C_0 <- solve(Sigma0) %*% Sigma %*% solve(Sigma0)
can_vec=pracma::sqrtm(Sigma0)$Binv %*% svd(pracma::sqrtm(Sigma0)$B %*% C_0 %*%pracma::sqrtm(Sigma0)$B)$u
a=can_vec[,1]
scale <- a %*% pracma::sqrtm(diag(r)+t(a) %*% Sigma %*% a/lambda)$B;

#########################SGCA WITHOUT CV##################################

ag <- sgca_init(A=S, B=sigma0hat, rho=0.5*sqrt(log(p)/n),K=r ,nu=1,trace=FALSE)
ainit <- init_process(ag$Pi, r)
final <- sgca_tgd(A=S, B=sigma0hat,r,ainit,k,lambda = 0.01, eta=0.001,convergence=1e-6,maxiter=15000, plot = TRUE)

####################SGCA WITH CV OVER LAMBDA###############################

lambda_values <- c(0, 1e-5,1e-4,1e-3,1e-2, 0.1, 1,10,100,1000,1e+4,1e+5)

cv <- sgca_cv_init_and_final(
  X = X, pp = pp, r = r, k = k,
  lambda_grid = lambda_values,
  rho_scale = 0.5,   # fixed
  nfold = 5, ncores = 5,
  eta = 0.001
)

cv$init_mean
cv$best_lambda

U_init_hat  <- cv$U_full_init
U_final_hat <- cv$U_full_final















