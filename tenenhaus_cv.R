library(MASS)
library(stats)
library(pracma)
library(tidyverse)
library(ggplot2)
library(Matrix)
library(tidyr)
library(geigen)
library(RGCCA)
library(dplyr)
library(expm)
library(foreach)
library(doParallel)

rgcca_holdout_score <- function(Y_test, connection = NULL, scheme = "factorial", bias = TRUE, comp = 1) {
  J <- length(Y_test)
  if (is.null(connection)) connection <- 1 - diag(J)
  
  g <- if (is.function(scheme)) scheme else switch(
    scheme,
    horst    = function(x) x,
    centroid = function(x) abs(x),
    factorial= function(x) x^2,
    stop("Unknown scheme: ", scheme)
  )
  
  cov_b <- function(u, v) {
    u <- u[, comp]; v <- v[, comp]
    u <- u - mean(u); v <- v - mean(v)
    denom <- if (bias) length(u) else (length(u) - 1)
    sum(u * v) / denom
  }
  
  score <- 0
  for (j in 1:(J - 1)) for (k in (j + 1):J) {
    if (connection[j, k] != 0) {
      score <- score + connection[j, k] * g(cov_b(Y_test[[j]], Y_test[[k]]))
    }
  }
  score
}

rgcca_unsupervised_cv_tau <- function(
    blocks, lambda_values,
    connection = NULL,
    scheme = "factorial",
    ncomp = 1,
    kfold = 5,
    n_cores = max(1, parallel::detectCores() - 1),
    seed = 1,
    scale = TRUE,
    scale_block = TRUE,
    bias = TRUE
) {
  # ---- ENSURE BLOCK NAMES ----
  if (is.null(names(blocks)) || any(names(blocks) == "")) {
    names(blocks) <- paste0("block", seq_along(blocks))
  }
  block_names <- names(blocks)
  
  J <- length(blocks)
  n <- nrow(blocks[[1]])
  if (is.null(connection)) connection <- 1 - diag(J)
  
  set.seed(seed)
  fold_id <- sample(rep(seq_len(kfold), length.out = n))
  
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  future::plan(future::multisession, workers = n_cores)
  
  fold_scores <- future_lapply(seq_len(kfold), function(f) {
    tr <- which(fold_id != f)
    te <- which(fold_id == f)
    
    train_blocks <- lapply(blocks, function(B) B[tr, , drop = FALSE])
    test_blocks  <- lapply(blocks, function(B) B[te, , drop = FALSE])
    
    # ---- PRESERVE NAMES (this fixes your error) ----
    names(train_blocks) <- block_names
    names(test_blocks)  <- block_names
    
    sapply(lambda_values, function(tau_scalar) {
      tau_vec <- rep(tau_scalar, J)
      
      fit <- tryCatch(
        rgcca(
          blocks = train_blocks,
          connection = connection,
          method = "rgcca",
          tau = tau_vec,
          ncomp = ncomp,
          scheme = scheme,
          scale = scale,
          scale_block = scale_block,
          bias = bias,
          verbose = FALSE
        ),
        error = function(e) NULL
      )
      if (is.null(fit)) return(NA_real_)
      
      trans <- tryCatch(rgcca_transform(fit, blocks_test = test_blocks),
                        error = function(e) NULL)
      if (is.null(trans)) return(NA_real_)
      
      # rgcca_transform() may return either list-of-Y or an object containing $Y
      Y_test <- if (!is.null(trans$Y)) trans$Y else trans
      
      rgcca_holdout_score(Y_test, connection = connection, scheme = scheme, bias = bias, comp = 1)
    })
  }, future.seed = TRUE)
  
  score_mat <- do.call(rbind, fold_scores)  # kfold x length(lambda_values)
  cv_mean <- colMeans(score_mat, na.rm = TRUE)
  cv_sd   <- apply(score_mat, 2, sd, na.rm = TRUE)
  
  best_idx <- which.max(replace(cv_mean, is.na(cv_mean), -Inf))
  best_tau <- lambda_values[best_idx]
  
  fit_full <- rgcca(
    blocks = blocks,
    connection = connection,
    method = "rgcca",
    tau = rep(best_tau, J),
    ncomp = ncomp,
    scheme = scheme,
    scale = scale,
    scale_block = scale_block,
    bias = bias,
    verbose = TRUE
  )
  
  list(
    lambda_values = lambda_values,
    fold_scores = score_mat,
    cv_mean = cv_mean,
    cv_sd = cv_sd,
    best_tau = best_tau,
    fit_full = fit_full
  )
}



# IMPORTANT: blocks MUST be a named list

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
                           
lambda_values <- c(0, 0.1, 0.25, 0.5, 0.75, 1)

cv_unsup <- rgcca_unsupervised_cv_tau(
  blocks = blocks,
  lambda_values = lambda_values,
  kfold = 5,
  n_cores = 4
)

cv_unsup$best_tau
fit_best <- cv_unsup$fit_full

