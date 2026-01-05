library(tidyverse)
# IMPORTANT: blocks MUST be a named list
source("R/tenenhaus_cv.R")
source("R/multicca_cv.R")
source("R/gao_cv_functions.R")
source("R/sgcar.R")
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
a_star=can_vec[,1]

blocks=list()
blocks[[1]]=X[,idx1]
blocks[[2]]=X[,idx2]
blocks[[3]]=X[,idx3]
names(blocks) <- paste0("block", seq_along(blocks))

block_names <- names(blocks)
lambda_values <- c(0.000001, 0.001, 0.1, 0.25, 0.5, 0.75, 1, 10, 100)

cv_unsup <- rgcca_unsupervised_cv_tau(
  blocks = blocks,
  lambda_values = lambda_values,
  kfold = 5,
  n_cores = 4
)

cv_unsup$best_tau
cv_unsup$fold_scores
colnames(blocks[[1]]) <- sapply(1:ncol(blocks[[1]]), function(x){paste0("block1_V", x)})
colnames(blocks[[2]]) <- sapply(1:ncol(blocks[[2]]), function(x){paste0("block2_V", x)})
colnames(blocks[[3]]) <- sapply(1:ncol(blocks[[3]]), function(x){paste0("block3_V", x)})

fit_best <- cv_unsup$fit_full
a_tenehaus = rgcca_transform(fit_best, blocks_test = blocks)
print(paste0("Error for Tenehaus is: ", sqrt(mean((c(fit_best$a[[1]], fit_best$a[[2]],
                                                    fit_best$a[[3]])- a_star)^2) )))


cv_out <- MultiCCA_unsup_cv(
  xlist = blocks,
  lambda_values = lambda_values,   # common penalty for all blocks
  type = "standard",
  ncomponents = 1,
  niter = 25,
  nfold = 5,
  workers = 4
)

cv_out$best_penalty
fit_best <- cv_out$fit_full
# weights:
fit_best$ws
print(paste0("Error for MultiCCA is: ", sqrt(mean((c(fit_best$ws[[1]], fit_best$ws[[2]],
                                                     fit_best$ws[[3]])- a_star)^2) )))


#########################SGCA WITHOUT CV##################################
k=10

ag <- sgca_init(A=S, B=sigma0hat, rho=0.5*sqrt(log(p)/n),K=r ,nu=1,trace=FALSE)
ainit <- init_process(ag$Pi, r)
final <- sgca_tgd(A=S, B=sigma0hat, r,ainit,k,lambda = 0.001, real_astar = a_star,
                  eta=0.001,convergence=1e-6,maxiter=15000, plot = TRUE)
print(paste0("Error for Gao SGCA is: ", sqrt(mean((final- a_star)^2) )))
####################SGCA WITH CV OVER LAMBDA###############################

lambda_values <- c(1e-5,1e-4,1e-3,1e-2, 0.1, 1,10,100,1000,1e+4,1e+5)

X, pp, r, k,
lambda_grid,
rho_scale = 1,
nfold = 5,
eta = 1e-3,
ridge_B = 1e-6,
ridge_norm = 1e-8,
nu = 1,
epsilon = 5e-3,
maxiter_admm = 1000,
convergence = 1e-6,
maxiter_tgd = 15000,
parallel=TRUE,
renorm_by_sigma0 = TRUE,  # TRUE: normalize on B_val; FALSE: normalize on B_train
center = TRUE,
ncores = max(1, parallel::detectCores() - 1),
seed = 2023


cv <- gao_gca_cv_init_and_final(
  X = X, pp = pp, r = r, k = k,
  lambda_grid = lambda_values,
  rho_scale = 0.5,   # fixed
  nfold = 5, ncores = 5,
  parallel=TRUE,
  eta = 0.001
)

cv$init_mean
cv$best_lambda

U_init_hat  <- cv$U_full_init
U_final_hat <- cv$U_full_final
print(paste0("Error for ADMM is: ", sqrt(mean((cv$U_full_final- a_star)^2) )))

fit_admm <- admm_sgca(Sigma, Sigma0, 0.0001, r,
                      rho = 1,
                      penalty = "l1")
print(paste0("Error for ADMM is: ", sqrt(mean((fit_admm$U- a_star)^2) )))


fit_admm_cv <- sgcar_cv_folds(
    X,
    p_list=pp,
    lambda_values,
    r,
    K = 10,
    folds = NULL,
    seed = NULL,
    penalty = "l1",
    loss_weight = NULL,
    relative_loss = FALSE,
    lambda_order =  "increasing",
    warm_start =  "CZU",
    penalize = "all",
    # solver controls
    rho = 1,
    # parallel controls
    parallel = FALSE,
    verbose = TRUE
)
print(paste0("Error for ADMM is: ", sqrt(mean((fit_admm_cv$fit_min$U- a_star)^2) )))

apply(fit_admm_cv$cvloss, 1, mean)

