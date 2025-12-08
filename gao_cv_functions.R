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
library(future.apply)

sqrtm_sym <- function(M, ridge = 0){
  M <- (M + t(M))/2
  ev <- eigen(M, symmetric = TRUE)
  d  <- pmax(ev$values, ridge)
  ev$vectors %*% diag(sqrt(d), length(d)) %*% t(ev$vectors)
}

invsqrt_sym <- function(M, ridge = 1e-8){
  M <- (M + t(M))/2
  ev <- eigen(M, symmetric = TRUE)
  d  <- pmax(ev$values, ridge)
  ev$vectors %*% diag(1/sqrt(d), length(d)) %*% t(ev$vectors)
}

## ---- init_process that always returns a p x r matrix (NO sqrtm needed) ----
init_process <- function(Pi, r){
  s <- svd(Pi)
  if (r == 1){
    s$u[, 1, drop = FALSE] * sqrt(s$d[1])
  } else {
    s$u[, 1:r, drop = FALSE] %*% diag(sqrt(s$d[1:r]), r, r)
  }
}

## ---- hard threshold robust to vectors ----
hard <- function(U, k){
  if (is.null(dim(U))) U <- matrix(U, ncol = 1)
  p <- nrow(U); r_local <- ncol(U)
  k <- min(k, p)
  
  if (r_local > 1){
    row_energy <- rowSums(U^2)
    keep <- order(row_energy, decreasing = TRUE)[seq_len(k)]
    U[-keep, ] <- 0
  } else {
    absu <- abs(U[, 1])
    keep <- order(absu, decreasing = TRUE)[seq_len(k)]
    U[-keep, 1] <- 0
  }
  U
}

## ---- safer sgca_init (symmetrize + tau guard) ----
sgca_init_fixed <- function(A,B,rho,K,nu=1,epsilon=5e-3,maxiter=1000,trace=FALSE){
  A <- (A + t(A))/2
  B <- (B + t(B))/2
  p <- nrow(B)
  
  evB <- eigen(B, symmetric = TRUE)
  vals <- pmax(evB$values, 0)
  sqB  <- evB$vectors %*% diag(sqrt(vals), p, p) %*% t(evB$vectors)
  
  tau <- 4 * nu * (max(vals)^2)
  if (!is.finite(tau) || tau <= 0) tau <- 1
  
  criteria <- Inf
  iter <- 0
  H <- Pi <- oldPi <- diag(1, p)
  Gamma <- matrix(0, p, p)
  
  while(criteria > epsilon && iter < maxiter){
    for (j in 1:20){
      Pi <- updatePi(B, sqB, A, H, Gamma, nu, rho, Pi, tau)
    }
    H <- updateH(sqB, Gamma, nu, Pi, K)
    Gamma <- Gamma + (sqB %*% Pi %*% sqB - H) * nu
    criteria <- sqrt(sum((Pi - oldPi)^2))
    oldPi <- Pi
    iter <- iter + 1
    if (trace) cat("iter:", iter, "crit:", criteria, "\n")
  }
  list(Pi=Pi,H=H,Gamma=Gamma,iteration=iter,convergence=criteria)
}

## ---- sgca_tgd that NEVER calls pracma::sqrtm (uses eigen + ridge) ----
sgca_tgd_safe <- function(A, B, r, init, k, lambda = 0.01, eta = 0.01,
                          convergence = 1e-3, maxiter = 10000,
                          ridge_norm = 1e-8){
  
  A <- (A + t(A))/2
  B <- (B + t(B))/2
  
  ut <- hard(init, k)
  ut <- ut %*% sqrtm_sym(diag(r) + (t(ut) %*% A %*% ut) / lambda, ridge = 0)
  
  criteria <- Inf
  iter <- 0
  
  while(criteria > convergence && iter < maxiter){
    G <- t(ut) %*% B %*% ut - diag(r)
    grad <- -A %*% ut + lambda * B %*% ut %*% G
    vt <- ut - eta * grad
    vt <- hard(vt, k)
    
    if (!all(is.finite(vt))) {
      return(matrix(NA_real_, nrow = nrow(A), ncol = r))
    }
    
    criteria <- sqrt(sum((ut - vt)^2))
    ut <- vt
    iter <- iter + 1
  }
  
  ut %*% invsqrt_sym(t(ut) %*% B %*% ut, ridge = ridge_norm)
}


make_mask_pp <- function(pp){
  p <- sum(pp)
  Mask <- matrix(0, p, p)
  start <- 1
  for (g in seq_along(pp)){
    idx <- start:(start + pp[g] - 1)
    Mask[idx, idx] <- 1
    start <- start + pp[g]
  }
  Mask
}

trace_obj <- function(S, U){
  SU <- S %*% U
  sum(U * SU)
}

## ---- 5-fold CV for lambda (parallel across folds) ----
sgca_cv_lambda <- function(X, pp, r, k, lambda_grid,
                           nfold = 5,
                           eta = 1e-3,
                           rho_scale = 0.5,
                           ridge_B = 1e-6,
                           ridge_norm = 1e-8,
                           nu = 1,
                           epsilon = 5e-3,
                           maxiter_admm = 1000,
                           convergence = 1e-6,
                           maxiter_tgd = 15000,
                           renorm_by_sigma0 = TRUE,
                           center = TRUE,
                           ncores = max(1, parallel::detectCores() - 1),
                           seed = 2023){
  
  X <- as.matrix(X)
  n <- nrow(X); p <- ncol(X)
  stopifnot(sum(pp) == p)
  
  Mask <- make_mask_pp(pp)
  set.seed(seed)
  fold_id <- sample(rep(1:nfold, length.out = n))
  
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  future::plan(future::multisession, workers = ncores)
  
  fold_scores <- future_lapply(1:nfold, function(f){
    
    tr <- which(fold_id != f)
    va <- which(fold_id == f)
    
    Xtr <- X[tr, , drop = FALSE]
    Xva <- X[va, , drop = FALSE]
    
    if (center){
      mu <- colMeans(Xtr)
      Xtr <- sweep(Xtr, 2, mu, "-")
      Xva <- sweep(Xva, 2, mu, "-")
    }
    
    ntr <- nrow(Xtr); nva <- nrow(Xva)
    
    S_tr <- crossprod(Xtr) / ntr; S_tr <- (S_tr + t(S_tr))/2
    S_va <- crossprod(Xva) / nva; S_va <- (S_va + t(S_va))/2
    
    B_tr <- S_tr * Mask
    B_va <- S_va * Mask
    
    # tiny ridge to avoid fold singularities
    if (ridge_B > 0){
      B_tr <- B_tr + ridge_B * diag(p)
      B_va <- B_va + ridge_B * diag(p)
    }
    B_tr <- (B_tr + t(B_tr))/2
    B_va <- (B_va + t(B_va))/2
    
    rho <- rho_scale * sqrt(log(p) / ntr)
    
    ag <- tryCatch(
      sgca_init_fixed(A = S_tr, B = B_tr, rho = rho, K = r,
                      nu = nu, epsilon = epsilon, maxiter = maxiter_admm),
      error = function(e) NULL
    )
    if (is.null(ag)) return(rep(NA_real_, length(lambda_grid)))
    
    ainit <- init_process(ag$Pi, r)
    
    vapply(lambda_grid, function(lam){
      tryCatch({
        U_tr <- sgca_tgd_safe(A = S_tr, B = B_tr, r = r, init = ainit, k = k,
                              lambda = lam, eta = eta,
                              convergence = convergence, maxiter = maxiter_tgd,
                              ridge_norm = ridge_norm)
        if (anyNA(U_tr)) return(NA_real_)
        
        U_use <- if (renorm_by_sigma0){
          U_tr %*% invsqrt_sym(t(U_tr) %*% B_va %*% U_tr, ridge = ridge_norm)
        } else U_tr
        
        trace_obj(S_va, U_use)
      }, error = function(e) NA_real_)
    }, numeric(1))
    
  }, future.seed = TRUE)
  
  score_mat <- do.call(rbind, fold_scores)
  cv_mean <- colMeans(score_mat, na.rm = TRUE)
  cv_sd   <- apply(score_mat, 2, sd, na.rm = TRUE)
  
  best_idx <- which.max(replace(cv_mean, is.na(cv_mean), -Inf))
  best_lambda <- lambda_grid[best_idx]
  
  # Refit full
  Xfull <- if (center) scale(X, center = TRUE, scale = FALSE) else X
  nfull <- nrow(Xfull)
  
  S_full <- crossprod(Xfull) / nfull; S_full <- (S_full + t(S_full))/2
  B_full <- S_full * Mask
  if (ridge_B > 0) B_full <- B_full + ridge_B * diag(p)
  B_full <- (B_full + t(B_full))/2
  
  rho_full <- rho_scale * sqrt(log(p) / nfull)
  
  ag_full <- sgca_init_fixed(A = S_full, B = B_full, rho = rho_full, K = r,
                             nu = nu, epsilon = epsilon, maxiter = maxiter_admm)
  
  ainit_full <- init_process(ag_full$Pi, r)
  
  U_full <- sgca_tgd_safe(A = S_full, B = B_full, r = r, init = ainit_full, k = k,
                          lambda = best_lambda, eta = eta,
                          convergence = convergence, maxiter = maxiter_tgd,
                          ridge_norm = ridge_norm)
  
  list(best_lambda = best_lambda,
       lambda_grid = lambda_grid,
       cv_mean = cv_mean,
       cv_sd = cv_sd,
       fold_scores = score_mat,
       U_full = U_full)
}

## 5-fold CV:
## - init-only score (single number per fold; rho fixed by rho_scale)
## - init+TGD score (vector over lambda_grid per fold)
## Returns BOTH full-data canonical vectors: init-only and init+TGD.
sgca_cv_init_and_final <- function(
    X, pp, r, k,
    lambda_grid,
    rho_scale = 1,          # <-- FIXED (no CV over rho)
    nfold = 5,
    eta = 1e-3,
    ridge_B = 1e-6,
    ridge_norm = 1e-8,
    nu = 1,
    epsilon = 5e-3,
    maxiter_admm = 1000,
    convergence = 1e-6,
    maxiter_tgd = 15000,
    renorm_by_sigma0 = TRUE,  # TRUE: normalize on B_val; FALSE: normalize on B_train
    center = TRUE,
    ncores = max(1, parallel::detectCores() - 1),
    seed = 2023
){
  
  X <- as.matrix(X)
  n <- nrow(X); p <- ncol(X)
  stopifnot(sum(pp) == p)
  
  Mask <- make_mask_pp(pp)
  
  set.seed(seed)
  fold_id <- sample(rep(1:nfold, length.out = n))
  
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  future::plan(future::multisession, workers = ncores)
  
  lambda_grid <- as.numeric(lambda_grid)
  
  fold_res <- future_lapply(1:nfold, function(f){
    
    tr <- which(fold_id != f)
    va <- which(fold_id == f)
    
    Xtr <- X[tr, , drop = FALSE]
    Xva <- X[va, , drop = FALSE]
    
    if (center){
      mu <- colMeans(Xtr)
      Xtr <- sweep(Xtr, 2, mu, "-")
      Xva <- sweep(Xva, 2, mu, "-")
    }
    
    ntr <- nrow(Xtr); nva <- nrow(Xva)
    
    S_tr <- crossprod(Xtr) / ntr; S_tr <- (S_tr + t(S_tr))/2
    S_va <- crossprod(Xva) / nva; S_va <- (S_va + t(S_va))/2
    
    B_tr <- S_tr * Mask
    B_va <- S_va * Mask
    
    if (ridge_B > 0){
      B_tr <- B_tr + ridge_B * diag(p)
      B_va <- B_va + ridge_B * diag(p)
    }
    B_tr <- (B_tr + t(B_tr))/2
    B_va <- (B_va + t(B_va))/2
    
    # ---- FIXED rho (no tuning) ----
    rho <- rho_scale * 1
    
    ag <- tryCatch(
      sgca_init_fixed(A = S_tr, B = B_tr, rho = rho, K = r,
                      nu = nu, epsilon = epsilon, maxiter = maxiter_admm),
      error = function(e) NULL
    )
    if (is.null(ag)){
      return(list(init_score = NA_real_,
                  final_scores = rep(NA_real_, length(lambda_grid))))
    }
    
    ainit <- init_process(ag$Pi, r)
    
    # ---- init-only canonical vectors + held-out scoring ----
    U0 <- hard(ainit, k)
    U0_use <- if (renorm_by_sigma0){
      U0 %*% invsqrt_sym(t(U0) %*% B_va %*% U0, ridge = ridge_norm)
    } else {
      U0 %*% invsqrt_sym(t(U0) %*% B_tr %*% U0, ridge = ridge_norm)
    }
    init_score <- trace_obj(S_va, U0_use)
    
    # ---- init+TGD scores over lambda ----
    final_scores <- vapply(lambda_grid, function(lam){
      tryCatch({
        U_tr <- sgca_tgd_safe(A = S_tr, B = B_tr, r = r, init = ainit, k = k,
                              lambda = lam, eta = eta,
                              convergence = convergence, maxiter = maxiter_tgd,
                              ridge_norm = ridge_norm)
        if (anyNA(U_tr)) return(NA_real_)
        
        U_use <- if (renorm_by_sigma0){
          U_tr %*% invsqrt_sym(t(U_tr) %*% B_va %*% U_tr, ridge = ridge_norm)
        } else U_tr
        
        trace_obj(S_va, U_use)
      }, error = function(e) NA_real_)
    }, numeric(1))
    
    list(init_score = init_score, final_scores = final_scores)
    
  }, future.seed = TRUE)
  
  ## ---- aggregate CV results ----
  init_vec  <- vapply(fold_res, `[[`, numeric(1), "init_score")
  final_mat <- do.call(rbind, lapply(fold_res, `[[`, "final_scores"))  # nfold x nlambda
  
  init_mean <- mean(init_vec, na.rm = TRUE)
  init_sd   <- sd(init_vec, na.rm = TRUE)
  
  final_mean <- colMeans(final_mat, na.rm = TRUE)
  final_sd   <- apply(final_mat, 2, sd, na.rm = TRUE)
  
  best_lambda_idx <- which.max(replace(final_mean, is.na(final_mean), -Inf))
  best_lambda <- lambda_grid[best_lambda_idx]
  
  ## ---- refit on FULL data (both init-only and init+TGD) ----
  Xfull <- if (center) scale(X, center = TRUE, scale = FALSE) else X
  nfull <- nrow(Xfull)
  
  S_full <- crossprod(Xfull) / nfull; S_full <- (S_full + t(S_full))/2
  B_full <- S_full * Mask
  if (ridge_B > 0) B_full <- B_full + ridge_B * diag(p)
  B_full <- (B_full + t(B_full))/2
  
  rho_full <- rho_scale * sqrt(log(p) / nfull)
  
  ag_full <- sgca_init_fixed(A = S_full, B = B_full, rho = rho_full, K = r,
                             nu = nu, epsilon = epsilon, maxiter = maxiter_admm)
  ainit_full <- init_process(ag_full$Pi, r)
  
  # init-only on full
  U_full_init <- hard(ainit_full, k)
  U_full_init <- U_full_init %*% invsqrt_sym(t(U_full_init) %*% B_full %*% U_full_init,
                                             ridge = ridge_norm)
  
  # init + TGD on full with best lambda
  U_full_final <- sgca_tgd_safe(A = S_full, B = B_full, r = r, init = ainit_full, k = k,
                                lambda = best_lambda, eta = eta,
                                convergence = convergence, maxiter = maxiter_tgd,
                                ridge_norm = ridge_norm)
  
  list(
    rho_scale = rho_scale,
    
    # init-only CV
    init_fold_scores = init_vec,
    init_mean = init_mean,
    init_sd = init_sd,
    
    # final CV (lambda)
    lambda_grid = lambda_grid,
    final_fold_scores = final_mat,
    final_mean = final_mean,
    final_sd = final_sd,
    best_lambda = best_lambda,
    
    # full-data canonical vectors
    U_full_init = U_full_init,
    U_full_final = U_full_final
  )
}



