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
library(PMA)
library(future.apply)

# --- Split a single matrix X into a list of blocks, given block sizes pp ---
make_xlist_from_X <- function(X, pp, block_names = paste0("X", seq_along(pp))) {
  stopifnot(sum(pp) == ncol(X))
  idx <- split(seq_len(ncol(X)), rep(seq_along(pp), times = pp))
  xlist <- lapply(idx, function(cols) X[, cols, drop = FALSE])
  names(xlist) <- block_names
  xlist
}

# --- Standardize using TRAIN statistics (avoid leakage) ---
standardize_train_test <- function(Xtr, Xte) {
  mu <- colMeans(Xtr)
  sdv <- apply(Xtr, 2, sd)
  sdv[!is.finite(sdv) | sdv == 0] <- 1
  Xtr_s <- sweep(sweep(Xtr, 2, mu, "-"), 2, sdv, "/")
  Xte_s <- sweep(sweep(Xte, 2, mu, "-"), 2, sdv, "/")
  list(Xtr = Xtr_s, Xte = Xte_s)
}

# --- Held-out score: sum_{comp} sum_{i<j} cor(Xi w_i, Xj w_j) ---
multicca_score <- function(xlist_test, ws, ncomponents = 1) {
  K <- length(xlist_test)
  score <- 0
  for (comp in seq_len(ncomponents)) {
    for (i in 2:K) {
      yi <- drop(xlist_test[[i]] %*% ws[[i]][, comp])
      for (j in 1:(i - 1)) {
        yj <- drop(xlist_test[[j]] %*% ws[[j]][, comp])
        cval <- suppressWarnings(cor(yi, yj))
        if (!is.finite(cval)) cval <- 0
        score <- score + cval
      }
    }
  }
  score
}

# --- Main CV function ---
MultiCCA_unsup_cv <- function(
    xlist,
    lambda_values,          # vector OR KxL matrix of candidate penalties
    type = "standard",      # or a length-K vector
    ncomponents = 1,
    niter = 25,
    nfold = 5,
    seed = 1,
    standardize = TRUE,
    parallel = TRUE,
    workers = max(1, parallel::detectCores() - 1),
    trace = FALSE
) {
  K <- length(xlist)
  n <- nrow(xlist[[1]])
  stopifnot(all(vapply(xlist, nrow, integer(1)) == n))
  
  # Build penalty matrix: K x L (each column = one candidate penalty vector)
  if (is.null(dim(lambda_values))) {
    L <- length(lambda_values)
    penalty_mat <- matrix(rep(lambda_values, each = K), nrow = K, ncol = L)
  } else {
    penalty_mat <- as.matrix(lambda_values)
    stopifnot(nrow(penalty_mat) == K)
    L <- ncol(penalty_mat)
  }
  
  set.seed(seed)
  fold_id <- sample(rep(seq_len(nfold), length.out = n))
  
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  if (parallel) future::plan(future::multisession, workers = workers) else future::plan(future::sequential)
  
  fold_scores <- future_lapply(seq_len(nfold), function(f) {
    tr <- which(fold_id != f)
    te <- which(fold_id == f)
    
    # Split train/test by block
    xtr <- lapply(xlist, function(Xk) Xk[tr, , drop = FALSE])
    xte <- lapply(xlist, function(Xk) Xk[te, , drop = FALSE])
    
    # Standardize using train stats, and reuse SVD init within this fold
    if (standardize) {
      scaled <- Map(standardize_train_test, xtr, xte)
      xtr_s <- lapply(scaled, `[[`, "Xtr")
      xte_s <- lapply(scaled, `[[`, "Xte")
    } else {
      xtr_s <- xtr
      xte_s <- xte
    }
    
    # Precompute SVD-based init ws for this fold (saves time across L penalties)
    ws_init <- lapply(xtr_s, function(Xk) {
      v <- svd(Xk)$v
      matrix(v[, seq_len(ncomponents), drop = FALSE], ncol = ncomponents)
    })
    
    scores <- rep(NA_real_, L)
    for (ell in seq_len(L)) {
      pen_vec <- penalty_mat[, ell]
      
      fit <- tryCatch(
        MultiCCA(
          xlist = xtr_s,
          penalty = pen_vec,
          ws = ws_init,
          niter = niter,
          type = type,
          ncomponents = ncomponents,
          trace = trace,
          standardize = FALSE  # we already standardized (if requested)
        ),
        error = function(e) NULL
      )
      
      if (!is.null(fit)) {
        scores[ell] <- multicca_score(xte_s, fit$ws, ncomponents = ncomponents)
      }
    }
    scores
  }, future.seed = TRUE)
  
  score_mat <- do.call(rbind, fold_scores)   # nfold x L
  cv_mean <- colMeans(score_mat, na.rm = TRUE)
  cv_sd   <- apply(score_mat, 2, sd, na.rm = TRUE)
  
  best_idx <- which.max(replace(cv_mean, is.na(cv_mean), -Inf))
  best_penalty <- penalty_mat[, best_idx]
  
  # Refit on full data with best penalty
  fit_full <- MultiCCA(
    xlist = xlist,
    penalty = best_penalty,
    ws = NULL,
    niter = niter,
    type = type,
    ncomponents = ncomponents,
    trace = trace,
    standardize = standardize
  )
  
  list(
    penalty_mat = penalty_mat,
    fold_scores = score_mat,
    cv_mean = cv_mean,
    cv_sd = cv_sd,
    best_idx = best_idx,
    best_penalty = best_penalty,
    fit_full = fit_full
  )
}

# ---- Example usage ----
# Suppose you have X and block sizes pp:
# xlist <- make_xlist_from_X(X, pp)

# Candidate penalty values (for type="standard", must be between 1 and sqrt(p_k))
lambda_values <- seq(1.2, 5, length.out = 10)

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

