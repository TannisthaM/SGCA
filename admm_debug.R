library(MASS)
library(stats)
library(pracma)
library(tidyverse)
library(ggplot2)
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
                      verbose = TRUE,
                      C0 = NULL) {                     # <-- NEW: reference matrix for distance
  
  penalty  <- match.arg(penalty)
  penalize <- match.arg(penalize)
  p_all <- nrow(Sigma0)
  
  # ----- NEW: set/validate C0 and preallocate histories -----
  if (is.null(C0)) {
    C0 <- matrix(0, p_all, p_all)
  } else if (!all(dim(C0) == c(p_all, p_all))) {
    stop("C0 must be a p x p matrix with p = nrow(Sigma0).")
  }
  dist_C_C0   <- rep(NA_real_, max_iter)  # ‖C - C0‖_F history
  r_norm_hist <- rep(NA_real_, max_iter)  # ‖C - Z‖_F history
  
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
      mask[(count+1):(count+u), (count+1):(count+u)] <- 0
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
    debug <- TRUE
    debug_every <- 1L   # break every iteration (or set 50, 100, etc.)
    if (debug && (iter %% debug_every == 0 ||
                  any(!is.finite(C)) || any(!is.finite(Z)) ||
                  r_norm > 1e6 || is.na(r_norm) || is.na(s_norm))) {
      browser()  # inspect C, Z, U, denom, etc.
    }
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
    
    # ----- NEW: record histories -----
    r_norm_hist[iter] <- r_norm
    dist_C_C0[iter]   <- norm(C - C0, "F")
    
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
  
  n_iter <- iter  # number of iterations actually executed
  
  list(
    C = C,Z=Z,
    U = U_canon,
    C_sparsity = mean(abs(C) < sparsity_threshold),
    U_sparsity = mean(abs(U_canon) < sparsity_threshold),
    # NEW: histories (trimmed to n_iter)
    r_norm_hist = r_norm_hist[seq_len(n_iter)],   # per-iter ‖C - Z‖_F
    dist_C_C0   = dist_C_C0[seq_len(n_iter)],     # per-iter ‖C - C0‖_F
    n_iter = n_iter
  )
}


###SIMULATION EXAMPLE################################################3

pp <- c(10,10,10);
p <- sum(pp)
s  <- c(1:3);
r <- 1;
rho<-1

n=100000

set.seed(2023)
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

###############checking if sigma and sigma0 are pd###########################33

Sigma_svd=svd(Sigma)
Sigma_svd$d
Sigma0_svd=svd(Sigma0)
Sigma0_svd$d

#############################################################################3

X <-mvrnorm(n, rep(0, p) , Sigma)
S <- t(X)%*% X / n
sigma0hat <- S * Mask
C_0 <- solve(Sigma0) %*% Sigma %*% solve(Sigma0)

fit_admm_oracle <- admm_sgca(S, sigma0hat, 0, r,
                             rho = 1,
                             p_list = pp,
                             penalty = c("l1"),
                             penalize = c("all"),   # used for l1 and l21_rows
                             
                             max_iter = 1,
                             abs_tol = 1e-18, rel_tol = 1e-17,
                             adapt_rho = FALSE, mu = 10, tau_incr = 2, tau_decr = 2,
                             sparsity_threshold = 1e-4,
                             verbose = TRUE)

fit_admm_GT <- admm_sgca(Sigma, Sigma0, 0, r,
                         rho = 1,
                         p_list = pp,
                         penalty = c("l1"),
                         penalize = c("all"),   # used for l1 and l21_rows
                         
                         max_iter = 1,
                         abs_tol = 1e-18, rel_tol = 1e-17,
                         adapt_rho = FALSE, mu = 10, tau_incr = 2, tau_decr = 2,
                         sparsity_threshold = 1e-4,
                         verbose = TRUE)

