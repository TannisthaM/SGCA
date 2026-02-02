library(SGCAR)
library(MASS)
library(Matrix)
library(future.apply)
library(PMA)
library(pracma)
library(RGCCA)
library(geigen)
library(CVXR)

subdistance <- function(A, B){
  svdresult = svd(t(A) %*% B)
  U = svdresult$u
  V = svdresult$v
  O = U %*% t(V)
  l = norm(A %*% O - B, type = c("F"))
  return(l)
}

library(CVXR)

suppressWarnings(library(CVXR, warn.conflicts = FALSE))

solve_C_supported_idx123 <- function(S, sigma0hat, pp1, pp2, pp3,
                                     solver = NULL, ridge = 0, verbose = FALSE) {
  S <- as.matrix(S)
  sigma0hat <- as.matrix(sigma0hat)
  p <- nrow(S)
  stopifnot(ncol(S) == p, all(dim(sigma0hat) == c(p, p)))
  
  # Indices
  idx1 <- 1:pp1
  idx2 <- (pp1 + 1):(pp1 + pp2)
  idx3 <- (pp1 + pp2 + 1):(pp1 + pp2 + pp3)
  idx_keep <- c(idx1, idx2, idx3)
  q <- length(idx_keep)
  
  if (max(idx_keep) > p) stop("pp1+pp2+pp3 must be <= p.")
  
  A <- sigma0hat[, idx_keep, drop = FALSE]   # p x q
  B <- sigma0hat[idx_keep, , drop = FALSE]   # q x p
  
  # Decision variable: only the active submatrix, symmetric
  Ck <- Variable(q, q, symmetric = TRUE)
  
  resid <- S - A %*% Ck %*% B
  obj <- sum_squares(resid)                  
  
  if (ridge > 0) obj <- obj + ridge * sum_squares(Ck)
  
  prob <- Problem(Minimize(obj))
  
  # Choose solver
  if (is.null(solver)) {
    solver <- if (requireNamespace("osqp", quietly = TRUE)) "OSQP" else "SCS"
  }
  
  res <- solve(prob, solver = solver, verbose = verbose)
  
  # Embed back into full p x p symmetric matrix with zeros elsewhere
  C_hat <- matrix(0, p, p)
  C_hat[idx_keep, idx_keep] <- res$getValue(Ck)
  C_hat <- (C_hat + t(C_hat)) / 2
  
  list(
    C_hat = C_hat,
    C_keep = res$getValue(Ck),
    idx1 = idx1, idx2 = idx2, idx3 = idx3, idx_keep = idx_keep,
    objective = res$value,
    status = res$status,
    solver = solver
  )
}



pp=c(10,10,10);   N=c(30, 1000, 10000, 50000, 100000)

s  <- 1:5
r  <- 1
k  <- 20
eta <- 0.001
max_iter <- 15000
rho <- 1
lambda <- 0.01

lambda_values <- c(0, 1e-7,1e-6,1e-5,1e-4,1e-3,1e-2, 0.1, 1,10,100,1000,1e+4,1e+5,1e+6,1e+7)
tau_values    <- c(0.000001, 0.001, 0.1, 0.25, 0.5, 0.75, 1)

sgca_ncores   <- 4
rgcca_ncores  <- 4
mcca_workers  <- 4
admm_ncores   <- 4

p <- sum(pp)
p_list <- pp

sgca2 <- list(); rgcca2 <- list(); sgcca2 <- list(); multicca <- list()
admm  <- list(); oracle <- list(); oracle2 <- list(); naive <- list()

naive$time <- c(); naive$loss <- c()
sgca2$time.final <- c(); sgca2$loss <- list()
rgcca2$time <- c(); rgcca2$loss <- c()
sgcca2$time <- c(); sgcca2$loss <- c()
multicca$time <- c(); multicca$loss <- c()
admm$time <- c(); admm$loss <- c()
oracle$time <- c();  oracle$loss <- c()
oracle2$time <- c(); oracle2$loss <- c()
# pp is your p_list, X is n x p
idxs <- .block_indices(pp)

# Ledoit–Wolf shrinkage per block (each returns a covariance matrix)
Sigma0_blocks_LW <- lapply(idxs, function(id) {
  cvCovEst::linearShrinkLWEst(dat = X[, id, drop = FALSE])
})

# Reassemble block diagonal Sigma0

C0_Cmat <- list()
cv_fit  <- list()

for (j in 1:5) {
  n=N[j]
  
  set.seed(2023)
  
  pp1 <- pp[1]; pp2 <- pp[2]; pp3 <- pp[3]
  p <- sum(pp)
  u1 <- matrix(0, pp1, r)
  u2 <- matrix(0, pp2, r)
  u3 <- matrix(0, pp3, r)
  
  idx1 <- 1:pp1
  idx2 <- (pp1+1):(pp1+pp2)
  idx3 <- (pp1+pp2+1):(pp1+pp2+pp3)
  
  Sigma <- diag(sum(pp))
  lambda = 0.9
  T1 <- toeplitz(0.5^(0:(pp1-1)))
  T2 <- toeplitz(0.7^(0:(pp2-1)))
  T3 <- toeplitz(0.9^(0:(pp3-1)))
  
  u1[s,1:r] <- matrix(rnorm(length(s)*r, mean=0, sd=1), length(s), r)
  u1 <- u1 %*% (pracma::sqrtm(t(u1[,1:r]) %*% T1 %*% u1[,1:r])$Binv)
  
  u2[s,1:r] <- matrix(rnorm(length(s)*r, mean=0, sd=1), length(s), r)
  u2 <- u2 %*% (pracma::sqrtm(t(u2[,1:r]) %*% T2 %*% u2[,1:r])$Binv)
  
  u3[s,1:r] <- matrix(rnorm(length(s)*r, mean=0, sd=1), length(s), r)
  u3 <- u3 %*% (pracma::sqrtm(t(u3[,1:r]) %*% T3 %*% u3[,1:r])$Binv)
  
  u_tot = rbind(u1, u2, u3)
  
  
  idx1 <- 1:pp1
  idx2 <- (pp1+1):(pp1+pp2)
  idx3 <- (pp1+pp2+1):(pp1+pp2+pp3)
  
  Sigma[idx1, idx1] <- T1
  Sigma[idx2, idx2] <- T2
  Sigma[idx3, idx3] <- T3
  
  Sigma[idx1, idx2] <-  lambda * T1 %*% u1 %*% t(u2) %*% T2
  Sigma[idx1, idx3] <- lambda * T1 %*% u1 %*% t(u3) %*% T3
  Sigma[idx2, idx3] <- lambda * T2 %*% u2 %*% t(u3) %*% T3
  Sigma[idx2, idx1] <-  lambda * T2 %*% u2 %*% t(u1) %*% T1
  Sigma[idx3, idx1] <- lambda * T3 %*% u3 %*% t(u1) %*% T1
  Sigma[idx3, idx2] <- lambda * T3 %*% u3 %*% t(u2) %*% T2
  
  Sigma <- (Sigma + t(Sigma))/2  # numerical symmetry
  
  # ---- Mask and Sigma0 (unchanged names)
  Mask <- matrix(0, p, p)
  Mask[idx1, idx1] <- 1
  Mask[idx2, idx2] <- 1
  Mask[idx3, idx3] <- 1
  Sigma0 <- Sigma * Mask
  
  
  X <- MASS::mvrnorm(n, rep(0, p), Sigma)
  S <- crossprod(X) / n
  S <- (S + t(S))/2
  sigma0hat <- S * Mask
  
  C_0 <- solve(Sigma0) %*% Sigma %*% solve(Sigma0)
  can_vec <- pracma::sqrtm(Sigma0)$Binv %*%
    svd(pracma::sqrtm(Sigma0)$B %*% Sigma %*% pracma::sqrtm(Sigma0)$B)$u
  a <- can_vec[, 1]  # keep as vector like your code
  result = geigen(Sigma,Sigma0)
  evalues <- result$values
  evectors <-result$vectors
  a = evectors[, ncol(evectors)]
  
  tryCatch({
    start_time <- proc.time()
    
    cv <- SGCAR:::gao_gca_cv_init_and_final(
      X = X, pp = pp, r = r, k = k,
      lambda_grid = lambda_values,
      rho_scale = 0.5,
      nfold = 5,
      ncores = sgca_ncores,     # capped by SLURM
      parallel = TRUE,
      eta = 0.001
    )
    
    end_time <- proc.time()
    final.time <- (end_time - start_time)
    
    ainit  <- cv$U_full_init
    afinal <- cv$U_full_final
    
    sgca2$time.final[j] <- final.time[3]
    sgca2$loss[[j]] <- c(subdistance(ainit, a), subdistance(afinal, a))
  }, error = function(e) {
    cat("Error in SGCA iteration", j, "n=", n, ":", conditionMessage(e), "\n")
    sgca2$time.final[j] <- NA
    sgca2$loss[[j]] <- NA
  })
  
  blocks <- list()
  blocks[[1]] <- X[, idx1]
  blocks[[2]] <- X[, idx2]
  blocks[[3]] <- X[, idx3]
  
  # ---- RGCCA CV (method="rgcca")
  tryCatch({
    start_time <- proc.time()
    
    cv_unsup <- SGCAR:::rgcca_unsupervised_cv_tau(
      blocks = blocks,
      lambda_values = tau_values,
      kfold = 5,
      n_cores = rgcca_ncores,   # capped by SLURM
      method = "rgcca",
      seed = seed_rep
    )
    
    end_time <- proc.time()
    output <- cv_unsup$fit_full
    time_taken <- end_time - start_time
    
    U_SVD <- as.matrix(c(output$astar$block1[,1], output$astar$block2[,1], output$astar$block3[,1]), ncol=1)
    
    rgcca2$time[j] <- time_taken[3]
    rgcca2$loss[j] <- subdistance(U_SVD[,1:r], a)
  }, error = function(e) {
    cat("Error in RGCCA iteration", j, "n=", n, ":", conditionMessage(e), "\n")
    rgcca2$time[j] <- NA
    rgcca2$loss[j] <- NA
  })
  
  # ---- SGCCA CV (method="sgcca")
  tryCatch({
    start_time <- proc.time()
    
    cv_unsup <- SGCAR:::rgcca_unsupervised_cv_tau(
      blocks = blocks,
      lambda_values = tau_values,
      kfold = 5,
      n_cores = rgcca_ncores,   # capped by SLURM
      method = "sgcca",
      seed = seed_rep
    )
    
    end_time <- proc.time()
    output <- cv_unsup$fit_full
    time_taken <- end_time - start_time
    
    U_SVD <- as.matrix(c(output$astar$block1[,1], output$astar$block2[,1], output$astar$block3[,1]), ncol=1)
    
    sgcca2$time[j] <- time_taken[3]
    sgcca2$loss[j] <- subdistance(U_SVD[,1:r], a)
  }, error = function(e) {
    cat("Error in SGCCA iteration", j, "n=", n, ":", conditionMessage(e), "\n")
    sgcca2$time[j] <- NA
    sgcca2$loss[j] <- NA
  })
  
  # ---- MultiCCA CV
  tryCatch({
    start_time <- proc.time()
    
    cv_out <- SGCAR:::MultiCCA_unsup_cv(
      xlist = blocks,
      lambda_values = tau_values,
      type = "standard",
      ncomponents = 1,
      niter = 25,
      nfold = 5,
      workers = mcca_workers   # capped by SLURM
    )
    
    end_time <- proc.time()
    fit_best <- cv_out$fit_full
    time_taken <- end_time - start_time
    
    U_SVD <- as.matrix(c(fit_best$ws[[1]], fit_best$ws[[2]], fit_best$ws[[3]]), ncol=1)
    
    multicca$time[j] <- time_taken[3]
    multicca$loss[j] <- subdistance(U_SVD[,1:r], a)
  }, error = function(e) {
    cat("Error in MultiCCA iteration", j, "- skipping\n")
    multicca$time[j] <- NA
    multicca$loss[j] <- NA
  })
  
  # ---- ADMM CV (sgcar_cv_folds)
  tryCatch({
    Sigma0hat_LW <- as.matrix(Matrix::bdiag(Sigma0_blocks_LW))
    Sigma0hat_LW <- (Sigma0hat_LW + t(Sigma0hat_LW)) / 2
    temp=c()
    for(i in 1:16){
      fit_admm_GT <- .admm_sgca_run(
        prep = SGCAR:::.admm_sgca_prepare(
          Sigma  = S,
          Sigma0 = Sigma0hat_LW,
          rho = NA,
          p_list = p_list,
          penalty = "l1",
          penalize = "all",
          symmetrize_z = TRUE
        ),
        lambda = lambda_values[i],
        warm_start =  "none",
        r = 1,
        verbose=TRUE,
        adapt_rho = TRUE,
        compute_canon = TRUE,
        abs_tol = 1e-10,
        rel_tol = 1e-10,
        max_iter = 100000
      )
      temp[i]=min(sqrt(mean((fit_admm_GT$U[,1] - a)^2)),sqrt(mean((fit_admm_GT$U[,1] + a)^2)))
    }
    start_time <- proc.time()
    
    # Build args exactly like your call, but optionally pass cores if supported
    fit_admm_GT= .admm_sgca_run(
      prep = SGCAR:::.admm_sgca_prepare(
        Sigma  = S,
        Sigma0 = linearShrinkLWEst(X)*Mask,
        rho = NA,
        p_list = p_list,
        penalty = "l1",
        penalize = "all",
        symmetrize_z = TRUE
      ),
      lambda = lambda_values[which.min(temp)],
      warm_start =  "none",
      r = 1,
      verbose=TRUE,
      adapt_rho = TRUE,
      compute_canon = TRUE,
      abs_tol = 1e-10,
      rel_tol = 1e-10,
      max_iter = 100000
    )
    
    end_time <- proc.time()
    admm$cv[j] <- (end_time - start_time)[3]
    admm$loss[j] <- min(sqrt(mean((fit_admm_GT$U[,1] - a)^2)),
                        sqrt(mean((fit_admm_GT$U[,1] + a)^2)))
    
    
  }, error = function(e) {
    cat("Error in ADMM iteration", j, "n=", n, ":", conditionMessage(e), "\n")
    admm$time[j] <- NA
    admm$loss[j] <- NA
  })
  
  # ---- Oracle / Naive (as in your code)
  tryCatch({
    start_time <- proc.time()
    newid <- c(s, pp[1]+s, pp[1]+pp[2]+s)
    
    fit=solve_C_supported_idx123(S,sigma0hat,pp1,pp2,pp3)
    Sigma0_sqrt <- pracma::sqrtm(sigma0hat)$B
    target <- matmul(Sigma0_sqrt, matmul(fit$C_hat, Sigma0_sqrt))
    eU <- top_eigs_sym(target, r)
    U_svd <- eU$vectors
    # Not over ---- need to compute the correct normalization
    if (r == 1){
      U_canon <- matmul(fit$C_hat, Sigma0_sqrt) %*% U_svd  * 1/eU$values
    }else{
      U_canon <- matmul(fit$C_hat, Sigma0_sqrt) %*% U_svd %*% diag(1/eU$values)
    }
    end_time <- proc.time()
    oracle2$time[j] <- (end_time - start_time)[3]
    oracle2$loss[j] <- min(sqrt(mean((U_canon[,1] - a)^2)),
                           sqrt(mean((U_canon[,1] + a)^2)))
    
    start_time <- proc.time()
    sqrt_m_sigma <- pracma::sqrtm(sigma0hat)
    C_oracle <- solve(sigma0hat) %*% S %*% solve(sigma0hat)
    target_mat <- sqrt_m_sigma$B %*% C_oracle %*% sqrt_m_sigma$B
    SVD <- svd(target_mat)
    lam_inv <- diag(1 / SVD$d)
    U_canon <- C_oracle %*% sqrt_m_sigma$B %*% SVD$u %*% lam_inv
    
    end_time <- proc.time()
    naive$time[j] <- (end_time - start_time)[3]
    naive$loss[j] <- subdistance(U_canon[,1], a)
    
    
    
    ###########loss for r=1 ######################################
    ###########loss for r=1 ######################################
    
    fit_admm_GT <- .admm_sgca_run(
      prep = SGCAR:::.admm_sgca_prepare(
        Sigma  = Sigma,
        Sigma0 = Sigma0,
        rho = NA,
        p_list = p_list,
        penalty = "l1",
        penalize = "all",
        symmetrize_z = TRUE
      ),
      lambda = 0,
      warm_start =  "none",
      r = 1,
      verbose=TRUE,
      adapt_rho = TRUE,
      compute_canon = TRUE,
      abs_tol = 1e-10,
      rel_tol = 1e-10,
      max_iter = 100000
    )
    oracle$loss[j] <- min(sqrt(mean((fit_admm_GT$U[,1] - a)^2)),
                          sqrt(mean((fit_admm_GT$U[,1] + a)^2)))
  }, error = function(e) {
    cat("Error in Oracle/Naive iteration", j, "- skipping\n")
    oracle2$time[j] <- NA; oracle2$loss[j] <- NA
    oracle$time[j]  <- NA; oracle$loss[j]  <- NA
    naive$time[j]   <- NA; naive$loss[j]   <- NA
  })
  
  cat("Done j=", j, " n=", n, "\n")
}

expt=list(sgca2=sgca2,multicca=multicca,admm=admm,oracle=oracle,oracle2=oracle2,naive=naive)

library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)

## N values
N <- c(30, 1000, 10000, 50000, 100000)

## ---------------------------------------------
## Extract SGCA2 loss components (Gao Fantope / Final)
## ---------------------------------------------
gao_fantope <- vapply(expt$sgca2$loss, function(x) x[1], numeric(1))
gao_final   <- vapply(expt$sgca2$loss, function(x) x[2], numeric(1))

## =========================
## TIME PLOT
## =========================
time_df <- tibble(
  N = N,
  `SGCA2 final`   = expt$sgca2$time.final,
  `Multicca time` = expt$multicca$time,
  `ADMM cv`       = expt$admm$cv,
  `Oracle1`       = expt$oracle$time,
  `Oracle2`       = expt$oracle2$time,
  `Naive`         = expt$naive$time
)

time_long <- time_df %>%
  pivot_longer(
    cols = -N,
    names_to  = "method",
    values_to = "time"
  )

time_plot <- ggplot(time_long, aes(x = N, y = time, color = method)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  scale_x_log10(breaks = N) +
  labs(
    x = "n (log scale)",
    y = "Time",
    color = "Method",
    title = "Computation time vs n"
  ) +
  theme_minimal()

print(time_plot)

## =========================
## LOSS PLOT
## =========================
loss_df <- tibble(
  N = N,
  `Gao Fantope` = gao_fantope,
  `Gao Final`   = gao_final,
  `Multicca`    = expt$multicca$loss,
  `ADMM`        = expt$admm$loss,
  `Oracle1`     = expt$oracle$loss,
  `Oracle2`     = expt$oracle2$loss,
  `Naive`       = expt$naive$loss
)

loss_long <- loss_df %>%
  pivot_longer(
    cols = -N,
    names_to  = "method",
    values_to = "loss"
  )

loss_plot <- ggplot(loss_long, aes(x = N, y = loss, color = method)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  scale_x_log10(breaks = N) +
  labs(
    x = "n (log scale)",
    y = "Loss",
    color = "Method",
    title = "Loss vs n"
  ) +
  theme_minimal()

print(loss_plot)


library(ggplot2); C=fit_admm_GT$C

plot_C_heatmap <- function(C, title = "Heatmap of C") {
  stopifnot(is.matrix(C), nrow(C) == ncol(C))
  
  df <- expand.grid(i = seq_len(nrow(C)), j = seq_len(ncol(C)))
  df$value <- as.vector(C)
  
  ggplot(df, aes(x = j, y = i, fill = value)) +
    geom_tile() +
    scale_y_reverse() +                  # puts (1,1) in top-left like a matrix
    coord_equal() +
    scale_fill_gradient2(midpoint = 0) + # centered at 0 (no manual colors)
    labs(title = title, x = "Column", y = "Row", fill = "Value") +
    theme_minimal() +
    theme(panel.grid = element_blank())
}

plot_C_heatmap(C, "C (30x30)")

