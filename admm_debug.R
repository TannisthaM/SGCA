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
