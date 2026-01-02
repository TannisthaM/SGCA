###SIMULATION EXAMPLE################################################3
library(ggplot)
source("R/admm_sgca.R")
source("R/parallel_cv.R")
source("R/utils.R")
source("R/sgca_cv.R")

pp <- c(10,10,10);
p <- sum(pp)
s  <- c(1:3);
r <- 1;
rho<-1

n=100000

set.seed(2023)
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
Sigma = Sigma + t(Sigma) - SigmaD;


Mask = matrix(0, sum(pp),sum(pp));
Mask[idx1,idx1] <- matrix(1,pp[1],pp[1]);
Mask[idx2,idx2] <- matrix(1,pp[2],pp[2]);
Mask[idx3,idx3] <- matrix(1,pp[3],pp[3]);
Sigma0 = Sigma * Mask;

X <-mvrnorm(n, rep(0, p) , Sigma)
S <- t(X)%*% X / n
sigma0hat <- S * Mask
C_0 <- solve(Sigma0) %*% Sigma %*% solve(Sigma0)

fit_admm_oracle <- admm_sgca(S, sigma0hat, 0, r,
                      rho = 1,
                      p_list = pp,
                      penalty = c("l1"),
                      penalize = c("all"),   # used for l1 and l21_rows
                      
                      max_iter = 500000,
                      abs_tol = 1e-18, rel_tol = 1e-17,
                      adapt_rho = FALSE, mu = 10, tau_incr = 2, tau_decr = 2,
                      sparsity_threshold = 1e-4,
                      verbose = TRUE)

fit_admm_GT <- admm_sgca(Sigma0, Sigma, 0, r,
                             rho = 1,
                             p_list = pp,
                             penalty = c("l1"),
                             penalize = c("all"),   # used for l1 and l21_rows
                             
                             max_iter = 500000,
                             abs_tol = 1e-18, rel_tol = 1e-17,
                             adapt_rho = FALSE, mu = 10, tau_incr = 2, tau_decr = 2,
                             sparsity_threshold = 1e-4,
                             verbose = TRUE)


fit <- cv_admm_sgca_C(
  X, p_list, lambdas = 10^seq(-3, 0, length.out = 12), r = 1, K = 5,
  penalty = "l1", penalize = "all", loss_part = "all", relative_loss = TRUE,
  rho = 1, max_iter = 200, adapt_rho = TRUE,   # used only if admm_sgca supports them
  parallel = TRUE, nb_cores = 5, blas_threads = 1,
  verbose = TRUE
)




