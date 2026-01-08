library(MASS)
library(stats)
library(pracma)
library(tidyverse)
library(Matrix)
library(tidyr)

matmul <- function(A, B) {
  if (requireNamespace("SMUT", quietly = TRUE)) {
    SMUT::eigenMapMatMult(A, B)
  } else {
    A %*% B
  }
}

soft_threshold <- function(X, T) {
  # T can be scalar or same shape as X
  sign(X) * pmax(abs(X) - T, 0)
}

sym_inv_sqrt <- function(S, eps = 1e-10) {
  S <- (S + t(S))/2
  ee <- eigen(S, symmetric = TRUE)
  vals <- pmax(ee$values, eps)
  V <- ee$vectors
  matmul(matmul(V, diag(1/sqrt(vals), nrow = length(vals))), t(V))
}

top_eigs_sym <- function(A, r) {
  A <- (A + t(A))/2
  if (requireNamespace("RSpectra", quietly = TRUE) && r < nrow(A)) {
    out <- RSpectra::eigs_sym(A, k = r, which = "LM")
    list(values = Re(out$values), vectors = Re(out$vectors))
  } else {
    ev <- eigen(A, symmetric = TRUE)
    list(values = ev$values[seq_len(r)], vectors = ev$vectors[, seq_len(r), drop = FALSE])
  }
}

.block_indices <- function(plist) {
  edges <- c(0, cumsum(plist))
  lapply(seq_along(plist), function(i) (edges[i] + 1):edges[i + 1])
}


# Allow passing S plus optional S0 (otherwise derive S0 by copying diagonal blocks)
.build_sigma_pair <- function(S, plist, S0 = NULL) {
  S <- (S + t(S)) / 2
  ptot <- nrow(S)
  stopifnot(sum(plist) == ptot)
  
  if (is.null(S0)) {
    S0 <- matrix(0, ptot, ptot)
    idxs <- .block_indices(plist)
    for (idx in idxs) S0[idx, idx] <- S[idx, idx]
  } else {
    S0 <- (S0 + t(S0)) / 2
  }
  list(Sigma = S, Sigma0 = S0)
}

## row/group proxes only used if you select penalties other than "l1"
.prox_l21_rows <- function(X, tau, mask, row_weights) {
  Z <- X
  p <- nrow(X)
  for (i in seq_len(p)) {
    v <- X[i, ] * mask[i, ]
    nrm <- sqrt(sum(v^2))
    if (nrm > 0) {
      shrink <- max(1 - (tau * row_weights[i]) / nrm, 0)
      v <- v * shrink
    }
    Z[i, ] <- v + X[i, ] * (1 - mask[i, ])
  }
  (Z + t(Z))/2
}

.prox_l21_groups <- function(X, lambda_over_rho, groups_l21 = list(), group_weights = NULL) {
  if (length(groups_l21) == 0) return(X)
  if (is.null(group_weights)) group_weights <- rep(1, length(groups_l21))
  Z <- X
  for (g in seq_along(groups_l21)) {
    idx <- groups_l21[[g]]
    if (is.matrix(idx) && ncol(idx) == 2) {
      vals <- X[cbind(idx[,1], idx[,2])]
    } else {
      vals <- X[idx]
    }
    nrm <- sqrt(sum(vals^2))
    shrink <- if (nrm > 0) max(1 - lambda_over_rho * group_weights[g] / nrm, 0) else 0
    if (is.matrix(idx) && ncol(idx) == 2) {
      Z[cbind(idx[,1], idx[,2])] <- vals * shrink
    } else {
      Z[idx] <- vals * shrink
    }
  }
  (Z + t(Z))/2
}




#' Canonical Correlation Analysis via Reduced Rank Regression (RRR)
#'
#' Estimates canonical directions using various RRR solvers and penalties.
#'
#' @param Sigma Matrix of predictors.
#' @param Sigma0 Matrix of responses.

#' @return A list with elements:
#' \itemize{
#'   \item U: Canonical direction matrix for X (p x r)
#'   \item V: Canonical direction matrix for Y (q x r)
#'   \item cor: Canonical covariances
#'   \item loss: The prediction error 1/n * || XU - YV ||^2
#' }

#' @export
admm_sgca <- function(Sigma, Sigma0, lambda, r,
                      rho = 1,
                      p_list = NULL,
                      penalty = c("l1","l21_rows","l21_groups"),
                      penalize = c("offdiag","all", "block"),   # used for l1 and l21_rows
                      weight = NULL,                   # l1 weights (matrix)
                      row_weights = NULL,              # l21_rows per-row weights (vector)
                      groups_l21 = NULL,               # list of 2-col index matrices (l21_groups)
                      group_weights = NULL,            # optional weights per group
                      symmetrize_z = TRUE,
                      max_iter = 4000,
                      abs_tol = 1e-4, rel_tol = 1e-3,
                      adapt_rho = FALSE, mu = 10, tau_incr = 2, tau_decr = 2,
                      sparsity_threshold = 1e-4,
                      verbose = TRUE) {
  
  penalty  <- match.arg(penalty)
  penalize <- match.arg(penalize)
  p_all <- nrow(Sigma0)
  
  # EVD of Σ0
  ev0 <- eigen((Sigma0 + t(Sigma0))/2, symmetric = TRUE)
  U0 <- ev0$vectors
  lam2 <- pmax(ev0$values, 0)
  Lam2 <- diag(lam2, nrow = length(lam2))
  
  # constants for the C update in Σ0-basis
  Sigma_tilde <- matmul(t(U0), matmul(Sigma, U0))
  const_rhs   <- matmul(matmul(Lam2, Sigma_tilde), Lam2)
  denom <- rho + outer(lam2^2, lam2^2, `*`)
  
  # init
  C <- matrix(0, p_all, p_all)
  Z <- C
  U <- C  # scaled dual
  
  # masks & weights for l1 / l21_rows
  mask <- matrix(1, p_all, p_all)
  if (penalize == "offdiag"){
    diag(mask) <- 0
  }
  if (penalize == "block"){
    count = 0
    for (u in p_list){
      mask[(count+1):(count+u), (count+1):(count+u)] <-0
      count = count + u
    }
  }
  

  if (is.null(weight))      weight      <- matrix(1, p_all, p_all)
  if (is.null(row_weights)) row_weights <- rep(1, p_all)
  
  for (iter in seq_len(max_iter)) {
    Z_prev <- Z
    
    # ----- C update (closed form) -----
    rhs_tilde <- rho * matmul(t(U0), matmul(Z - U, U0)) + const_rhs
    
    
    
    C_tilde <- rhs_tilde / denom
    C <- matmul(U0, matmul(C_tilde, t(U0)))
    
    # ----- Z update (prox) -----
    Z_tilde <- C + U
    if (penalty == "l1") {
      if (all(weight == 1) && all(mask == 1)) {
        Z <- soft_threshold(Z_tilde, lambda / rho)
      } else {
        Z <- Z_tilde
        idx <- (mask == 1)
        Z[idx] <- soft_threshold(Z_tilde[idx], (lambda / rho) * weight[idx])
      }
    } else if (penalty == "l21_rows") {
      Z <- .prox_l21_rows(Z_tilde, tau = (lambda / rho), mask = mask, row_weights = row_weights)
    } else if (penalty == "l21_groups") {
      Z <- .prox_l21_groups(Z_tilde, lambda_over_rho = (lambda / rho),
                            groups_l21 = groups_l21 %||% list(),
                            group_weights = group_weights)
    }
    
    if (symmetrize_z) Z <- (Z + t(Z))/2  # keep things well-behaved for the next C-step
    
    # ----- dual update (scaled) -----
    U <- U + (C - Z)
    
    # ----- diagnostics -----
    r_norm <- norm(C - Z, "F")
    s_norm <- rho * norm(Z - Z_prev, "F")
    eps_pri  <- sqrt(p_all * p_all) * abs_tol + rel_tol * max(norm(C, "F"), norm(Z, "F"))
    eps_dual <- sqrt(p_all * p_all) * abs_tol + rel_tol * rho * norm(U, "F")
    
    if (verbose && iter %% 50 == 0) {
      cat(sprintf("iter %5d  r=%.3e  s=%.3e  eps_pri=%.3e  eps_dual=%.3e  rho=%.3g\n",
                  iter, r_norm, s_norm, eps_pri, eps_dual, rho))
    }
    if (r_norm <= eps_pri && s_norm <= eps_dual) break
    
    # ----- adaptive rho -----
    if (adapt_rho) {
      if (r_norm > mu * s_norm) {
        rho <- rho * tau_incr; U <- U / tau_incr
        denom <- rho + outer(lam2^2, lam2^2, `*`)
      } else if (s_norm > mu * r_norm) {
        rho <- rho / tau_decr; U <- U * tau_decr
        denom <- rho + outer(lam2^2, lam2^2, `*`)
      }
    }
  }
  
  # Step 2: eigenvectors of Σ0^{1/2} C Σ0^{1/2}
  Sigma0_sqrt <- matmul(U0, matmul(diag(sqrt(lam2), nrow = length(lam2)), t(U0)))
  target <- matmul(Sigma0_sqrt, matmul(C, Sigma0_sqrt))
  eU <- top_eigs_sym(target, r)
  U_svd <- eU$vectors
  
  # Step 3: normalization to enforce U^T Σ0 U = I
  B <- matmul(t(U_svd), matmul(Sigma0, U_svd))
  U_canon <- matmul(U_svd, sym_inv_sqrt(B))
  
  list(C = C, 
       U = U_canon,
       cor = sqrt(eU$values),
       C_sparsity = mean(abs(C) < sparsity_threshold),
       U_sparsity = mean(abs(U_canon) < sparsity_threshold))
}
