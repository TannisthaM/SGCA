library(MASS)
library(stats)
library(pracma)
library(tidyverse)
library(ggplot2)
library(Matrix)
library(tidyr)


soft_threshold <- function(X, tau) {
  sign(X) * pmax(abs(X) - tau, 0)
}

admm_l1_matrix <- function(Sigma0, Sigma, lambda, rho = 1.0, max_iters = 3000,
                           abs_tol = 1e-5, rel_tol = 1e-4,
                           ensure_symmetric = TRUE, verbose = FALSE,
                           use_sparse = FALSE,
                           C0 = NULL,                     # <- ground-truth C (optional, but recommended)
                           print_every = 50) {            # <- how often to print progress if verbose
  # Ensure matrices
  Sigma0 <- as.matrix(Sigma0)
  Sigma  <- as.matrix(Sigma)
  n <- nrow(Sigma0)
  n2 <- n * n
  
  # Optional ground-truth handling
  if (!is.null(C0)) {
    C0 <- as.matrix(C0)
    if (!all(dim(C0) == c(n, n))) {
      stop("C0 must be an n x n matrix with the same n as Sigma0.")
    }
    C0_F <- norm(C0, "F")
    if (!is.finite(C0_F) || C0_F == 0) C0_F <- 1.0
  } else {
    C0_F <- 1.0
  }
  
  # Precompute pieces used every iteration
  S0sq     <- Sigma0 %*% Sigma0               # Σ0^2
  S0_S_S0  <- Sigma0 %*% Sigma %*% Sigma0     # Σ0 Σ Σ0
  
  # ---- Build linear operator on vec(C):  [ (S0sq^T ⊗ S0sq) + rho*I ] vec(C) = vec(RHS) ----
  if (use_sparse) {
    S0sq_sp <- Matrix::Matrix(S0sq, sparse = TRUE)
    K <- Matrix::kronecker(Matrix::t(S0sq_sp), S0sq_sp) + rho * Matrix::Diagonal(n2)
    Kfac <- try(Matrix::Cholesky(K, LDL = FALSE), silent = TRUE)  # SPD factorization
    solveK <- if (!inherits(Kfac, "try-error")) {
      function(b) as.numeric(Matrix::solve(Kfac, b))
    } else {
      function(b) as.numeric(Matrix::solve(K, b))
    }
  } else {
    K <- kronecker(t(S0sq), S0sq) + rho * diag(n2)
    Kchol <- try(chol(K), silent = TRUE)  # SPD Cholesky
    solveK <- if (!inherits(Kchol, "try-error")) {
      function(b) as.numeric(backsolve(Kchol, forwardsolve(t(Kchol), b)))
    } else {
      function(b) as.numeric(solve(K, b))
    }
  }
  
  # ---- ADMM state ----
  C <- matrix(0, n, n)
  Z <- matrix(0, n, n)
  U <- matrix(0, n, n)
  
  history <- vector("list", max_iters)
  
  for (iter in seq_len(max_iters)) {
    # ---- C-step: solve Σ0^2 C Σ0^2 + ρ C = ρ(Z-U) + Σ0 Σ Σ0 via Kronecker ----
    RHS  <- rho * (Z - U) + S0_S_S0
    vecC <- solveK(as.vector(RHS))
    C    <- matrix(vecC, n, n)
    if (ensure_symmetric) C <- 0.5 * (C + t(C))
    
    # ---- Z-step: soft-threshold ----
    Z_prev <- Z
    Z <- soft_threshold(C + U, lambda / rho)
    if (ensure_symmetric) Z <- 0.5 * (Z + t(Z))
    
    # ---- U-step (scaled dual) ----
    U <- U + (C - Z)
    
    # ---- Diagnostics & stopping (SAME as admm_sgca) ----
    r_norm <- norm(C - Z, "F")
    s_norm <- rho * norm(Z - Z_prev, "F")
    eps_pri  <- sqrt(n * n) * abs_tol + rel_tol * max(norm(C, "F"), norm(Z, "F"))
    eps_dual <- sqrt(n * n) * abs_tol + rel_tol * rho * norm(U, "F")
    
    # Objective (optional monitor)
    obj_val <- 0.5 * norm(Sigma0 %*% C %*% Sigma0 - Sigma, "F")^2 + lambda * sum(abs(Z))
    
    # ---- C-loss versus ground truth (each iteration) ----
    if (!is.null(C0)) {
      c_loss     <- norm(C - C0, "F")
      c_rel_loss <- c_loss / C0_F
    } else {
      c_loss <- NA_real_
      c_rel_loss <- NA_real_
    }
    
    # Progress print
    if (verbose && (iter %% print_every == 0)) {
      cat(sprintf("iter %5d  obj=%.6e  r=%.3e  s=%.3e  eps_pri=%.3e  eps_dual=%.3e  rho=%.3g  C_loss=%.3e  C_rel=%.3e\n",
                  iter, obj_val, r_norm, s_norm, eps_pri, eps_dual, rho, c_loss, c_rel_loss))
    }
    
    history[[iter]] <- list(
      iter = iter,
      obj = obj_val,
      r = r_norm,
      s = s_norm,
      eps_pri = eps_pri,
      eps_dual = eps_dual,
      rho = rho,
      C_loss = c_loss,
      C_rel_loss = c_rel_loss
    )
    
    # Stop if converged
    if (r_norm <= eps_pri && s_norm <= eps_dual) {
      history <- history[1:iter]
      if (verbose) cat(sprintf("ADMM finished in %d iterations.\n", iter))
      break
    }
    
    # If reached max_iters without break, announce it
    if (iter == max_iters && verbose) {
      cat(sprintf("ADMM reached max_iters = %d.\n", max_iters))
    }
  }
  
  iters <- length(history)
  
  # Return traces as simple numeric vectors as well (convenient for plotting)
  c_loss_trace     <- vapply(history, function(h) h$C_loss,     numeric(1))
  c_rel_loss_trace <- vapply(history, function(h) h$C_rel_loss, numeric(1))
  
  list(
    C = C, Z = Z, U = U,
    history = history,              # per-iteration records incl. C_loss
    iters = iters,                  # number of iterations actually run
    c_loss_trace = c_loss_trace,    # numeric vector: ||C - C0||_F per iter (NA if no C0)
    c_rel_loss_trace = c_rel_loss_trace
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

fit_GT <- admm_l1_matrix(Sigma0, Sigma, 0, rho = 1.0, max_iters = 10000,
                      abs_tol = 1e-18, rel_tol = 1e-17,C0=C_0,
                      ensure_symmetric = TRUE, verbose = TRUE,print_every = 50)

fit_oracle <- admm_l1_matrix(sigma0hat, S, 0, rho = 1.0, max_iters = 10000,
                      abs_tol = 1e-18, rel_tol = 1e-17,C0=C_0,
                      ensure_symmetric = TRUE, verbose = TRUE,print_every = 50)

#############C_hat heatmap GT###############################################

C_mat = data.frame(fit_GT$C)
colnames(C_mat) = 1:ncol(C_mat)
C_mat["x"] = 1:nrow(C_mat)
C_mat_long = pivot_longer(C_mat, cols = -c("x"))

# Ensure `name` is treated as a factor with levels in ascending order
C_mat_long$name <- factor(C_mat_long$name, levels = as.character(sort(as.numeric(unique(C_mat_long$name)))))

# Create the plot
m <- max(abs(C_mat_long$value), na.rm = TRUE)  # keep white exactly at 0
ggplot(C_mat_long, aes(x = x, y = name, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue",
                       midpoint = 0, limits = c(-m, m)) +
  labs(
    x = "X-axis (Rows)",
    y = "Y-axis (Columns)",
    fill = "Value",
    title = "Heatmap of C_hat Matrix,lambda=0,n=100000,iter=10000,GT"
  ) +
  theme_minimal()

##############CVXR within ADMM function###############################################

admm_l1_matrix <- function(Sigma0, Sigma, lambda, rho = 1.0, max_iters = 3000,
                           abs_tol = 1e-5, rel_tol = 1e-4,
                           ensure_symmetric = TRUE, verbose = FALSE,
                           solver = c("SCS","OSQP","CVXOPT")) {
  solver <- match.arg(solver)
  
  # Ensure plain matrices
  Sigma0 <- as.matrix(Sigma0)
  Sigma  <- as.matrix(Sigma)
  n <- nrow(Sigma0)
  
  # Precompute fixed pieces
  S0sq      <- Sigma0 %*% Sigma0                  # Σ0^2
  S0_S_S0   <- Sigma0 %*% Sigma %*% Sigma0        # Σ0 Σ Σ0
  
  # ADMM state
  C <- matrix(0, n, n)
  Z <- matrix(0, n, n)
  U <- matrix(0, n, n)
  
  # ----- CVXR model for the C-step (built once) -----
  Cvar <- if (ensure_symmetric) CVXR::Variable(n, n, symmetric = TRUE) else CVXR::Variable(n, n)
  RHSpar <- CVXR::Parameter(n, n, name = "RHS")  # iteration-dependent RHS
  
  # Solve the linear matrix equation exactly:
  constraints <- list(S0sq %*% Cvar %*% S0sq + rho * Cvar == RHSpar)
  prob <- CVXR::Problem(CVXR::Minimize(0), constraints)
  
  history <- vector("list", max_iters)
  
  for (k in seq_len(max_iters)) {
    # ---- C-step: solve Σ0^2 C Σ0^2 + ρ C = ρ(Z-U) + Σ0 Σ Σ0 via CVXR ----
    RHS_k <- rho * (Z - U) + S0_S_S0
    res <- try(
      CVXR::solve(prob, solver = solver, warm_start = TRUE, verbose = FALSE,
                  params = list(RHS = RHS_k)),
      silent = TRUE
    )
    if (inherits(res, "try-error") || is.null(CVXR::value(Cvar))) {
      stop(sprintf("C-subproblem failed to solve at iter %d.", k))
    }
    C <- CVXR::value(Cvar)
    if (ensure_symmetric) C <- 0.5 * (C + t(C))
    
    # ---- Z-step: soft-threshold ----
    Z_prev <- Z
    Z <- soft_threshold(C + U, lambda / rho)
    if (ensure_symmetric) Z <- 0.5 * (Z + t(Z))
    
    # ---- U-step (scaled dual) ----
    U <- U + (C - Z)
    
    # ---- Residuals & stopping (now matches admm_sgca) ----
    r_norm <- norm(C - Z, "F")
    s_norm <- rho * norm(Z - Z_prev, "F")
    eps_pri  <- sqrt(n * n) * abs_tol + rel_tol * max(norm(C, "F"), norm(Z, "F"))
    eps_dual <- sqrt(n * n) * abs_tol + rel_tol * rho * norm(U, "F")
    
    # (Optional) objective value for monitoring
    obj_val <- 0.5 * norm(Sigma0 %*% C %*% Sigma0 - Sigma, "F")^2 + lambda * sum(abs(Z))
    
    if (verbose) {
      cat(sprintf("iter %4d  obj=%.6e  r=%.3e  s=%.3e  eps_pri=%.3e  eps_dual=%.3e  rho=%.3g\n",
                  k, obj_val, r_norm, s_norm, eps_pri, eps_dual, rho))
    }
    history[[k]] <- list(iter = k, obj = obj_val, r = r_norm, s = s_norm,
                         eps_pri = eps_pri, eps_dual = eps_dual, rho = rho)
    
    if (r_norm <= eps_pri && s_norm <= eps_dual) {
      history <- history[1:k]
      break
    }
  }
  
  list(C = C, Z = Z, U = U, history = history, iters = length(history))
}














