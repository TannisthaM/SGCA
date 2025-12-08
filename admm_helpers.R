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

## --- helpers
matmul <- function(A, B) {
  if (requireNamespace("SMUT", quietly = TRUE)) {
    SMUT::eigenMapMatMult(A, B)
  } else {
    A %*% B
  }
}

## Your Procrustes subspace distance
subdistance <- function(A, B){
  svdresult = svd(t(A) %*% B)
  U = svdresult$u
  V = svdresult$v
  O = U %*% t(V)
  l = norm(A %*% O - B, type = "F")
  return(l)
}

soft_threshold <- function(X, T) {
  sign(X) * pmax(abs(X) - T, 0)
}

sym_inv_sqrt <- function(S, eps = 1e-10) {
  S <- (S + t(S))/2
  ee <- eigen(S, symmetric = TRUE)
  vals <- pmax(ee$values, eps)
  V <- ee$vectors
  matmul(matmul(V, diag(1/sqrt(vals), nrow = length(vals))), t(V))
}

## NEW: symmetric square-root (needed to whiten by Σ0^{1/2})
sym_sqrt <- function(S, eps = 1e-10) {
  S <- (S + t(S))/2
  ee <- eigen(S, symmetric = TRUE)
  vals <- pmax(ee$values, eps)
  V <- ee$vectors
  matmul(matmul(V, diag(sqrt(vals), nrow = length(vals))), t(V))
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

.block_diag_cov <- function(X, p_list, use = "pairwise.complete.obs") {
  p <- ncol(X); stopifnot(sum(p_list) == p)
  Sigma0 <- matrix(0, p, p)
  start <- 1
  for (sz in p_list) {
    idx <- start:(start + sz - 1)
    Sigma0[idx, idx] <- stats::cov(X[, idx, drop = FALSE], use = use)
    start <- start + sz
  }
  (Sigma0 + t(Sigma0)) / 2
}

.make_folds <- function(n, K, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  idx <- sample.int(n)
  fold_id <- cut(seq_len(n), breaks = K, labels = FALSE)
  out <- integer(n); out[idx] <- fold_id
  out
}

.mask_from_part <- function(p, part = c("all","offdiag","block_off"), p_list = NULL) {
  part <- match.arg(part)
  M <- matrix(TRUE, p, p)
  if (part == "offdiag") {
    diag(M) <- FALSE
  } else if (part == "block_off") {
    if (is.null(p_list)) stop("p_list is required when part='block_off'.")
    start <- 1
    for (sz in p_list) {
      idx <- start:(start + sz - 1)
      M[idx, idx] <- FALSE
      start <- start + sz
    }
  }
  M
}

.cv_loss_C <- function(C, Sigma_val, Sigma0_val,
                       part = c("all","offdiag","block_off"),
                       p_list = NULL, weight = NULL,
                       relative = TRUE) {
  p <- nrow(Sigma_val)
  part <- match.arg(part)
  A <- (Sigma_val + t(Sigma_val)) / 2
  B <- Sigma0_val %*% ((C + t(C))/2) %*% Sigma0_val
  R <- A - B
  M <- .mask_from_part(p, part, p_list)
  W <- matrix(0, p, p); W[M] <- 1
  if (!is.null(weight)) {
    if (length(weight) == 1) weight <- matrix(weight, p, p)
    stopifnot(all(dim(weight) == c(p, p)))
    W <- W * weight
  }
  num <- sum((R * R) * W)
  if (!relative) return(num)
  den <- sum((A * A) * W)
  if (den <= .Machine$double.eps) return(NA_real_)
  num / den
}

## (Optional) row/group proxes only used if you select penalties other than "l1"
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

## --- canonical directions from (Σ, Σ0): U = Σ0^{-1/2} V_r
.U_from_sigma_pair <- function(Sigma, Sigma0, r, eps = 1e-10) {
  S0is <- sym_inv_sqrt(Sigma0, eps = eps)
  M <- matmul(S0is, matmul((Sigma + t(Sigma)) / 2, S0is))
  M <- (M + t(M)) / 2
  ev <- top_eigs_sym(M, r)
  list(U = matmul(S0is, ev$vectors), values = ev$values)  # U is Σ0-orthonormal
}

## --- Loss by U using your Procrustes subspace distance in Σ0 metric
## metric = "subdist" (raw Frobenius ∈ [0, sqrt(2r)]) or "subdist_rel" (0–1 normalized)
.cv_loss_U <- function(U_hat, Sigma_val, Sigma0_val, r,
                       metric = c("subdist","subdist_rel","proj","angles")) {
  metric <- match.arg(metric)
  p <- nrow(Sigma0_val)
  r <- min(r, ncol(U_hat), p)
  
  # Reference canonical basis on validation data
  U_ref <- .U_from_sigma_pair(Sigma_val, Sigma0_val, r)$U     # Σ0_val-orthonormal
  
  # Orthonormalize U_hat in Σ0_val metric
  U_hat_r <- U_hat[, seq_len(r), drop = FALSE]
  G_hat <- matmul(t(U_hat_r), matmul(Sigma0_val, U_hat_r))
  U_hat_orth <- matmul(U_hat_r, sym_inv_sqrt((G_hat + t(G_hat))/2))  # Σ0_val-orthonormal
  
  if (metric %in% c("proj","angles")) {
    S <- svd(matmul(t(U_ref), matmul(Sigma0_val, U_hat_orth)))$d
    S <- pmin(pmax(Re(S), 0), 1)
    if (metric == "proj") return(1 - mean(S^2))
    return(mean(acos(S)^2))
  }
  
  # Procrustes in Σ0_val metric: whiten to Euclidean orthonormal columns
  #S0_sqrt <- sym_sqrt(Sigma0_val)
  #A <- matmul(S0_sqrt, U_ref)
  #B <- matmul(S0_sqrt, U_hat_orth)
  #l <- subdistance(A, B)                  # ∈ [0, sqrt(2r)]
  l <- subdistance(U_hat_orth,U_ref) 
  if (metric == "subdist") return(l)
  (l^2) / (2 * r)                         # ∈ [0,1]
}

## --- helper: run code with base::browser suppressed (no change to admm_sgca)
.with_browser_suppressed <- function(expr) {
  # Mask base::browser() with a no-op in .GlobalEnv, then restore.
  existed <- exists("browser", envir = .GlobalEnv, inherits = FALSE)
  if (existed) old <- get("browser", envir = .GlobalEnv, inherits = FALSE)
  assign("browser", function(...) invisible(NULL), envir = .GlobalEnv)
  on.exit({
    if (existed) assign("browser", old, envir = .GlobalEnv)
    else rm("browser", envir = .GlobalEnv)
  }, add = TRUE)
  eval(substitute(expr), envir = parent.frame())
}         











