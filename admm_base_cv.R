## --- deps
library(MASS)
library(stats)
library(pracma)
library(tidyverse)
library(ggplot2)
library(Matrix)
library(tidyr)
library(geigen)
library(dplyr)
library(expm)
library(foreach)
library(doParallel)

##########################base ADMM Function###################################################################

admm_sgca <- function(Sigma, Sigma0, lambda, r,
                      rho = 1,
                      p_list = NULL,
                      penalty = c("l1","l21_rows","l21_groups"),
                      penalize = c("offdiag","all", "block"),
                      weight = NULL,
                      row_weights = NULL,
                      groups_l21 = NULL,
                      group_weights = NULL,
                      symmetrize_z = TRUE,
                      max_iter = 4000,
                      abs_tol = 1e-4, rel_tol = 1e-3,
                      adapt_rho = FALSE, mu = 10, tau_incr = 2, tau_decr = 2,
                      sparsity_threshold = 1e-4,
                      verbose = TRUE,
                      C0 = NULL) {
  penalty  <- match.arg(penalty)
  penalize <- match.arg(penalize)
  p_all <- nrow(Sigma0)
  
  if (is.null(C0)) {
    C0 <- matrix(0, p_all, p_all)
  } else if (!all(dim(C0) == c(p_all, p_all))) {
    stop("C0 must be a p x p matrix with p = nrow(Sigma0).")
  }
  dist_C_C0   <- rep(NA_real_, max_iter)
  r_norm_hist <- rep(NA_real_, max_iter)
  
  ev0 <- eigen((Sigma0 + t(Sigma0))/2, symmetric = TRUE)
  U0 <- ev0$vectors
  lam2 <- pmax(ev0$values, 0)
  Lam2 <- diag(lam2, nrow = length(lam2))
  
  Sigma_tilde <- matmul(t(U0), matmul(Sigma, U0))
  const_rhs   <- matmul(matmul(Lam2, Sigma_tilde), Lam2)
  
  denom <- rho + outer(lam2^2, lam2^2, `*`)
  
  C <- matrix(0, p_all, p_all)
  Z <- C
  U <- C
  
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
    
    rhs_tilde <- rho * matmul(t(U0), matmul(Z - U, U0)) + const_rhs
    C_tilde <- rhs_tilde / denom
    C <- matmul(U0, matmul(C_tilde, t(U0)))
    
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
    
    if (symmetrize_z) Z <- (Z + t(Z))/2
    
    U <- U + (C - Z)
    
    r_norm <- norm(C - Z, "F")
    s_norm <- rho * norm(Z - Z_prev, "F")
    eps_pri  <- sqrt(p_all * p_all) * abs_tol + rel_tol * max(norm(C, "F"), norm(Z, "F"))
    eps_dual <- sqrt(p_all * p_all) * abs_tol + rel_tol * rho * norm(U, "F")
    
    r_norm_hist[iter] <- r_norm
    dist_C_C0[iter]   <- norm(C - C0, "F")
    
    if (verbose && iter %% 50 == 0) {
      cat(sprintf("iter %5d  r=%.3e  s=%.3e  eps_pri=%.3e  eps_dual=%.3e  rho=%.3g\n",
                  iter, r_norm, s_norm, eps_pri, eps_dual, rho))
    }
    if (r_norm <= eps_pri && s_norm <= eps_dual) break
    
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
  
  ## NOTE: the following section is left AS IS in your code.
  Sigma0_sqrt <- matmul(U0, matmul(diag(sqrt(lam2), nrow = length(lam2)), t(U0)))
  SVD=svd(Sigma0_sqrt%*%C%*%Sigma0_sqrt)$u
  lam_inv=diag(1/svd(Sigma0_sqrt%*%C%*%Sigma0_sqrt)$d)
  U_canon=C%*%Sigma0_sqrt%*%SVD%*%lam_inv
  
  n_iter <- iter
  list(
    C = C, Z = Z,
    U = U_canon,
    C_sparsity = mean(abs(C) < sparsity_threshold),
    U_sparsity = mean(abs(U_canon) < sparsity_threshold),
    r_norm_hist = r_norm_hist[seq_len(n_iter)],
    dist_C_C0   = dist_C_C0[seq_len(n_iter)],
    n_iter = n_iter
  )
}



#############################ADMM CV wrt C#####################################################################

cv_admm_sgca_C <- function(
    X, p_list, lambdas, r,
    K = 5, folds = NULL, seed = NULL,
    penalty  = c("l1","l21_rows","l21_groups"),
    penalize = c("offdiag","all","block"),
    loss_part = c("auto","all","offdiag","block_off"),
    loss_weight = NULL,
    relative_loss = TRUE,
    suppress_browser = TRUE,
    verbose = TRUE,
    
    # foreach/doParallel controls
    parallel     = TRUE,
    nb_cores     = NULL,     # if NULL -> min(K, detectCores()-1)
    blas_threads = 1,        # threads per worker (1 avoids oversubscription)
    rng_seed     = NULL,     # optional: reproducible per-(fold,lambda) fits
    ...
) {
  penalty   <- match.arg(penalty)
  penalize  <- match.arg(penalize)
  loss_part <- match.arg(loss_part)
  
  n <- nrow(X); p <- ncol(X)
  stopifnot(sum(p_list) == p)
  
  # ---------- Folds ----------
  if (is.null(folds)) {
    folds <- .make_folds(n, K, seed = seed)
  } else {
    stopifnot(length(folds) == n)
    K <- length(unique(folds))
  }
  
  # ---------- Effective loss part ----------
  part_eff <- if (loss_part == "auto") {
    if (penalize == "offdiag") "offdiag"
    else if (penalize == "block") "block_off"
    else "all"
  } else loss_part
  
  L <- length(lambdas)
  if (verbose) {
    cat(sprintf("CV by C: %d folds × %d lambdas; n=%d, p=%d\n", K, L, n, p))
    cat(sprintf("  penalty=%s, penalize=%s, validation loss on: %s (relative=%s)\n",
                penalty, penalize, part_eff, relative_loss))
  }
  
  # ---------- Helper: block-diagonal-from-S ----------
  bd_from_S <- function(S, p_list) {
    ptot <- sum(p_list)
    S0 <- matrix(0, ptot, ptot)
    offs <- c(0L, cumsum(p_list))
    for (j in seq_along(p_list)) {
      id <- (offs[j] + 1L):offs[j+1L]
      S0[id, id] <- S[id, id]
    }
    S0
  }
  
  # ---------- Containers ----------
  cv_mat <- matrix(NA_real_, nrow = L, ncol = K,
                   dimnames = list(paste0("lambda=", signif(lambdas, 4)),
                                   paste0("Fold", seq_len(K))))
  
  # ---------- Full-data second-moment & S0 (no centering) ----------
  n_full  <- nrow(X)
  S_full  <- crossprod(X) / n_full
  S_full  <- 0.5 * (S_full + t(S_full))
  S0_full <- bd_from_S(S_full, p_list)
  
  # ---------- Precompute per-fold second moments via down-dating ----------
  if (verbose) cat("Precomputing per-fold second-moment matrices (down-dated)...\n")
  fold_stats <- vector("list", K)
  for (k in seq_len(K)) {
    idx_val <- which(folds == k)
    Xval    <- X[idx_val, , drop = FALSE]
    n_val   <- nrow(Xval)
    n_tr    <- n_full - n_val
    
    M_val <- crossprod(Xval)
    S_va  <- M_val / n_val
    S_va  <- 0.5 * (S_va + t(S_va))
    
    S_tr  <- (n_full * S_full - M_val) / n_tr
    S_tr  <- 0.5 * (S_tr + t(S_tr))
    
    S0_tr <- bd_from_S(S_tr, p_list)
    S0_va <- bd_from_S(S_va, p_list)
    
    fold_stats[[k]] <- list(S_tr = S_tr, S0_tr = S0_tr, S_va = S_va, S0_va = S0_va)
  }
  
  # ---------- Filter ... to supported args ----------
  admm_fun <- admm_sgca
  dots <- list(...)
  admm_formals <- tryCatch(setdiff(names(formals(admm_fun)), "..."),
                           error = function(e) character())
  valid_dots <- if (length(admm_formals)) dots[names(dots) %in% admm_formals] else list()
  
  # Optional browser suppressor (sequential + refit only)
  wb_fun <- tryCatch(get(".with_browser_suppressed", mode = "function"),
                     error = function(e) NULL)
  if (is.null(wb_fun)) wb_fun <- function(x) x
  
  .call_admm <- function(S_in, S0_in, lam) {
    base_args <- list(Sigma = S_in, Sigma0 = S0_in,
                      lambda = lam, r = r,
                      p_list = p_list,
                      penalty = penalty, penalize = penalize,
                      verbose = FALSE)
    if (suppress_browser) wb_fun(do.call(admm_fun, c(base_args, valid_dots)))
    else                  do.call(admm_fun, c(base_args, valid_dots))
  }
  
  # ============================
  # foreach/doParallel over folds
  # ============================
  if (isTRUE(parallel)) {
    if (!requireNamespace("doParallel", quietly = TRUE) ||
        !requireNamespace("foreach",    quietly = TRUE)) {
      warning("foreach/doParallel not available; running sequentially.")
      parallel <- FALSE
    }
  }
  
  if (isTRUE(parallel)) {
    nphys <- tryCatch(parallel::detectCores(logical = TRUE), error = function(e) 1L)
    cores <- if (is.null(nb_cores)) max(1L, min(K, nphys - 1L))
    else                   max(1L, min(K, as.integer(nb_cores)))
    
    if (verbose) cat(sprintf("Running CV in parallel over %d folds (doParallel)\n", cores))
    
    cl <- parallel::makeCluster(cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    doParallel::registerDoParallel(cl)
    
    # ---- IMPORTANT: export ALL dependencies that workers need ----
    # (This is the main fix.)
    need_global <- c(
      "admm_sgca", "matmul", "soft_threshold",
      ".cv_loss_C", ".mask_from_part",
      ".prox_l21_rows", ".prox_l21_groups",
      "%||%"
    )
    # only export those that exist (avoids errors if you never defined %||% etc.)
    need_global <- need_global[vapply(need_global, exists, logical(1), envir = .GlobalEnv, inherits = TRUE)]
    parallel::clusterExport(cl, varlist = need_global, envir = .GlobalEnv)
    
    # export local objects needed inside workers
    parallel::clusterExport(
      cl,
      varlist = c("lambdas","L","r","p_list","part_eff","loss_weight","relative_loss",
                  "penalty","penalize","valid_dots","rng_seed"),
      envir = environment()
    )
    
    # Thread controls (optional)
    parallel::clusterCall(
      cl,
      function(bt) {
        if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
          RhpcBLASctl::blas_set_num_threads(as.integer(bt))
          RhpcBLASctl::omp_set_num_threads(as.integer(bt))
        }
        Sys.setenv(OMP_NUM_THREADS = as.character(bt),
                   MKL_NUM_THREADS = as.character(bt),
                   OPENBLAS_NUM_THREADS = as.character(bt))
        NULL
      },
      as.integer(blas_threads)
    )
    
    # ---- key change: iterate fs = fold_stats so workers receive ONLY their fold ----
    res_list <- foreach::foreach(
      fid = seq_len(K),
      fs  = fold_stats,
      .inorder = TRUE,
      .errorhandling = "pass"
    ) %dopar% {
      if (!is.null(rng_seed)) set.seed(as.integer(rng_seed + 1000L * fid))
      
      losses <- rep(NA_real_, length(lambdas))
      for (li in seq_along(lambdas)) {
        lam <- lambdas[li]
        if (!is.null(rng_seed)) set.seed(as.integer(rng_seed + 1000L * fid + li))
        
        losses[li] <- tryCatch({
          base_args <- list(Sigma = fs$S_tr, Sigma0 = fs$S0_tr,
                            lambda = lam, r = r, p_list = p_list,
                            penalty = penalty, penalize = penalize,
                            verbose = FALSE)
          
          fit <- do.call(admm_sgca, c(base_args, valid_dots))
          .cv_loss_C(fit$C, fs$S_va, fs$S0_va,
                     part = part_eff, p_list = p_list,
                     weight = loss_weight, relative = relative_loss)
        }, error = function(e) NA_real_)
      }
      losses
    }
    
    # if any fold returned an error object, show it (helps debugging)
    err_idx <- which(vapply(res_list, inherits, logical(1), "error"))
    if (length(err_idx) && verbose) {
      message("Parallel fold(s) errored. Example error:\n  ",
              conditionMessage(res_list[[err_idx[1]]]))
    }
    
    # Combine into L×K
    for (k in seq_len(K)) {
      col <- res_list[[k]]
      if (is.numeric(col) && length(col) == L) cv_mat[, k] <- col else cv_mat[, k] <- NA_real_
    }
    
  } else {
    if (verbose) cat("Running CV sequentially\n")
    for (k in seq_len(K)) {
      if (verbose) cat(sprintf("  Fold %d/%d\n", k, K))
      fs <- fold_stats[[k]]
      if (!is.null(rng_seed)) set.seed(as.integer(rng_seed + 1000L * k))
      
      for (li in seq_len(L)) {
        lam <- lambdas[li]
        if (!is.null(rng_seed)) set.seed(as.integer(rng_seed + 1000L * k + li))
        
        loss <- tryCatch({
          fit <- .call_admm(fs$S_tr, fs$S0_tr, lam)
          .cv_loss_C(fit$C, fs$S_va, fs$S0_va,
                     part = part_eff, p_list = p_list,
                     weight = loss_weight, relative = relative_loss)
        }, error = function(e) NA_real_)
        
        cv_mat[li, k] <- loss
        if (verbose && (li %% max(1, floor(L/5)) == 0)) {
          cat(sprintf("    λ=%g  loss=%.5f\n", lam, loss))
        }
      }
    }
  }
  
  # ---------- Refill missing (fold, λ) cells sequentially ----------
  na_idx <- which(!is.finite(cv_mat), arr.ind = TRUE)
  if (nrow(na_idx) > 0L) {
    if (verbose) cat(sprintf("Refilling %d missing (fold, λ) cells sequentially...\n", nrow(na_idx)))
    for (ii in seq_len(nrow(na_idx))) {
      li <- na_idx[ii, 1]; k <- na_idx[ii, 2]
      fs <- fold_stats[[k]]
      lam <- lambdas[li]
      if (!is.null(rng_seed)) set.seed(as.integer(rng_seed + 1000L * k + li))
      loss <- tryCatch({
        fit <- .call_admm(fs$S_tr, fs$S0_tr, lam)
        .cv_loss_C(fit$C, fs$S_va, fs$S0_va,
                   part = part_eff, p_list = p_list,
                   weight = loss_weight, relative = relative_loss)
      }, error = function(e) NA_real_)
      cv_mat[li, k] <- loss
    }
  }
  
  # ---------- Choose λ among those evaluated on ALL folds ----------
  row_ok <- rowSums(is.finite(cv_mat)) == K
  if (!any(row_ok)) {
    stop("No lambda was successfully evaluated on all folds. Check solver errors or relax settings.")
  }
  
  cvm_sel  <- rowMeans(cv_mat[row_ok, , drop = FALSE], na.rm = FALSE)
  best_rel_idx <- which.min(cvm_sel)
  lambda_min   <- lambdas[row_ok][best_rel_idx]
  
  if (verbose) cat(sprintf("Refitting on full data at lambda_min = %g\n", lambda_min))
  fit_min <- .call_admm(S_full, S0_full, lambda_min)
  
  # report overall means (with na.rm=TRUE)
  full_cvm  <- rowMeans(cv_mat, na.rm = TRUE)
  full_cvsd <- apply(cv_mat, 1, function(x) {
    m <- sum(is.finite(x)); if (m <= 1) return(NA_real_)
    stats::sd(x, na.rm = TRUE) / sqrt(m)
  })
  
  out <- list(
    lambdas = lambdas,
    cvloss  = cv_mat,
    cvm     = full_cvm,
    cvsd    = full_cvsd,
    lambda_min = lambda_min,
    fit_min = fit_min,
    folds = folds,
    K = K,
    r = r,
    penalty = penalty,
    penalize = penalize,
    loss_part = part_eff,
    relative_loss = relative_loss
  )
  class(out) <- "cv_admm_sgca_C"
  out
}








