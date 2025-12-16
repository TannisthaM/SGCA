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

########################## base ADMM Function ###################################################################

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

############################# ADMM CV wrt C (parallel over lambdas) #############################################

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
    nb_cores     = NULL,     # if NULL -> min(L, detectCores()-1)
    blas_threads = 1,        # threads per worker (1 avoids oversubscription)
    rng_seed     = NULL,     # optional: reproducible per-(fold,lambda) fits
    ...
) {
  penalty   <- match.arg(penalty)
  penalize  <- match.arg(penalize)
  loss_part <- match.arg(loss_part)
  
  .log_line <- function(...) if (isTRUE(verbose)) cat(..., "\n")
  
  n <- nrow(X); p <- ncol(X)
  stopifnot(sum(p_list) == p)
  
  .log_line(sprintf("[cv] start | pid=%s | parallel requested=%s", Sys.getpid(), parallel))
  
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
    .log_line(sprintf("[cv] config | K=%d folds | L=%d lambdas | n=%d | p=%d", K, L, n, p))
    .log_line(sprintf("[cv] config | penalty=%s | penalize=%s | loss_part=%s | relative_loss=%s",
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
  
  # ---------- Full-data second-moment & S0 ----------
  n_full  <- nrow(X)
  S_full  <- crossprod(X) / n_full
  S_full  <- 0.5 * (S_full + t(S_full))
  S0_full <- bd_from_S(S_full, p_list)
  
  # ---------- Precompute per-fold moments ----------
  if (verbose) .log_line("[cv] precompute | building per-fold second moments (down-dated)...")
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
  if (verbose) .log_line("[cv] precompute | done")
  
  # ---------- Filter ... to supported args ----------
  admm_fun <- admm_sgca
  dots <- list(...)
  admm_formals <- tryCatch(setdiff(names(formals(admm_fun)), "..."),
                           error = function(e) character())
  valid_dots <- if (length(admm_formals)) dots[names(dots) %in% admm_formals] else list()
  
  # Optional browser suppressor
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
  
  # ---------- sequential runner ----------
  .run_sequential_cv <- function(reason = NULL) {
    if (!is.null(reason)) message(sprintf("[cv][SEQ] Running sequentially because: %s", reason))
    else .log_line("[cv][SEQ] Running CV sequentially")
    
    for (k in seq_len(K)) {
      if (verbose) .log_line(sprintf("[cv][SEQ] Fold %d/%d", k, K))
      fs <- fold_stats[[k]]
      for (li in seq_len(L)) {
        lam <- lambdas[li]
        loss <- tryCatch({
          fit <- .call_admm(fs$S_tr, fs$S0_tr, lam)
          .cv_loss_C(fit$C, fs$S_va, fs$S0_va,
                     part = part_eff, p_list = p_list,
                     weight = loss_weight, relative = relative_loss)
        }, error = function(e) NA_real_)
        cv_mat[li, k] <- loss
      }
    }
    invisible(NULL)
  }
  
  # ============================
  # foreach/doParallel over lambdas (FIXED EXPORT OF DOT-NAMES)
  # ============================
  parallel_ok_to_try <- isTRUE(parallel) &&
    requireNamespace("doParallel", quietly = TRUE) &&
    requireNamespace("foreach", quietly = TRUE) &&
    (L > 1L)
  
  did_parallel <- FALSE
  parallel_failed_reason <- NULL
  
  if (parallel_ok_to_try) {
    cl <- NULL
    
    ok <- tryCatch({
      nphys <- tryCatch(parallel::detectCores(logical = TRUE), error = function(e) 1L)
      cores <- if (is.null(nb_cores)) max(1L, min(L, nphys - 1L))
      else                   max(1L, min(L, as.integer(nb_cores)))
      
      .log_line(sprintf("[cv][PAR] Attempting parallel over lambdas | detected_cores=%d | using_cores=%d", nphys, cores))
      if (cores <= 1L) stop("only 1 core available/requested (cores<=1)")
      
      cl <- parallel::makeCluster(cores)
      doParallel::registerDoParallel(cl)
      
      .log_line(sprintf("[cv][PAR] foreach backend: %s | workers: %d",
                        foreach::getDoParName(), foreach::getDoParWorkers()))
      .log_line(sprintf("[cv][PAR] master pid=%s", Sys.getpid()))
      
      # load a few common deps on workers
      parallel::clusterEvalQ(cl, {
        suppressPackageStartupMessages({
          library(stats)
          library(Matrix)
          library(pracma)
          # these are only to support %||% if your code uses it:
          if (requireNamespace("rlang", quietly = TRUE)) library(rlang)
          if (requireNamespace("purrr", quietly = TRUE)) library(purrr)
        })
        NULL
      })
      
      # ---- export from wherever these symbols exist on master (and INCLUDE dot-names) ----
      sym_required <- c("admm_sgca", ".cv_loss_C")
      sym_optional <- c("soft_threshold", "matmul",
                        ".mask_from_part", ".prox_l21_rows", ".prox_l21_groups", "%||%")
      sym_to_ship  <- unique(c(sym_required, sym_optional))
      
      export_env <- new.env(parent = emptyenv())
      found <- logical(length(sym_to_ship)); names(found) <- sym_to_ship
      
      for (nm in sym_to_ship) {
        if (exists(nm, inherits = TRUE)) {
          assign(nm, get(nm, inherits = TRUE), envir = export_env)
          found[nm] <- TRUE
        } else {
          found[nm] <- FALSE
        }
      }
      
      missing_required <- sym_required[!found[sym_required]]
      if (length(missing_required) > 0L) {
        stop(sprintf("Required symbol(s) not found on master: %s",
                     paste(missing_required, collapse = ", ")))
      }
      
      # ***** FIX IS RIGHT HERE: all.names=TRUE so dot-names are included *****
      shipped <- ls(export_env, all.names = TRUE)
      parallel::clusterExport(cl, varlist = shipped, envir = export_env)
      
      .log_line(sprintf("[cv][PAR] exported symbols to workers: %s", paste(shipped, collapse = ", ")))
      
      # export local objects workers need
      parallel::clusterExport(
        cl,
        varlist = c("fold_stats","K",
                    "lambdas","L","r","p_list","part_eff","loss_weight","relative_loss",
                    "penalty","penalize","valid_dots","rng_seed"),
        envir = environment()
      )
      .log_line("[cv][PAR] exported local objects to workers")
      
      # verify required symbol exists on each worker
      chk <- parallel::clusterCall(cl, function() exists(".cv_loss_C", inherits = TRUE))
      .log_line(sprintf("[cv][PAR] worker check: .cv_loss_C exists on %d/%d workers", sum(unlist(chk)), length(chk)))
      if (!all(unlist(chk))) stop("At least one worker does not have .cv_loss_C after export; aborting parallel run.")
      
      # threads per worker
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
      .log_line(sprintf("[cv][PAR] set BLAS/OMP threads per worker to %d", as.integer(blas_threads)))
      
      .log_line("[cv][PAR] launching foreach %dopar% over lambdas ...")
      
      res_list <- foreach::foreach(
        li  = seq_len(L),
        lam = lambdas,
        .inorder = TRUE,
        .errorhandling = "pass"
      ) %dopar% {
        
        pid  <- Sys.getpid()
        host <- tryCatch(Sys.info()[["nodename"]], error = function(e) NA_character_)
        
        losses_k <- rep(NA_real_, K)
        first_err <- NULL
        
        for (k in seq_len(K)) {
          fs <- fold_stats[[k]]
          losses_k[k] <- tryCatch({
            base_args <- list(Sigma = fs$S_tr, Sigma0 = fs$S0_tr,
                              lambda = lam, r = r, p_list = p_list,
                              penalty = penalty, penalize = penalize,
                              verbose = FALSE)
            
            fit <- do.call(admm_sgca, c(base_args, valid_dots))
            
            .cv_loss_C(fit$C, fs$S_va, fs$S0_va,
                       part = part_eff, p_list = p_list,
                       weight = loss_weight, relative = relative_loss)
          }, error = function(e) {
            if (is.null(first_err)) first_err <<- conditionMessage(e)
            NA_real_
          })
        }
        
        list(li = li, lambda = lam, pid = pid, host = host, losses = losses_k, err = first_err)
      }
      
      # summarize pids used
      err_obj <- vapply(res_list, inherits, logical(1), "error")
      ok_idx <- which(!err_obj & vapply(res_list, function(x) is.list(x) && !is.null(x$losses), logical(1)))
      if (length(ok_idx) == 0L) stop("all parallel tasks failed (no valid results)")
      
      pids_used <- unique(vapply(res_list[ok_idx], function(x) x$pid, integer(1)))
      hosts_used <- unique(vapply(res_list[ok_idx], function(x) as.character(x$host), character(1)))
      
      .log_line(sprintf("[cv][PAR] completed | valid lambda results=%d/%d", length(ok_idx), L))
      .log_line(sprintf("[cv][PAR] unique worker PIDs observed: %s", paste(pids_used, collapse = ", ")))
      if (length(hosts_used)) .log_line(sprintf("[cv][PAR] worker hosts observed: %s", paste(hosts_used, collapse = ", ")))
      
      # fill cv_mat
      for (li0 in seq_len(L)) {
        x <- res_list[[li0]]
        if (is.list(x) && !is.null(x$losses) && length(x$losses) == K) {
          cv_mat[li0, ] <- x$losses
        } else {
          cv_mat[li0, ] <- NA_real_
        }
      }
      
      # if still NA, show an example error
      any_fold_err <- vapply(res_list, function(x) is.list(x) && !is.null(x$err) && nzchar(x$err), logical(1))
      if (any(any_fold_err)) {
        idx <- which(any_fold_err)[1]
        .log_line(sprintf("[cv][PAR] Example fold-level error from worker pid=%s, lambda=%g: %s",
                          res_list[[idx]]$pid, res_list[[idx]]$lambda, res_list[[idx]]$err))
      }
      
      did_parallel <- TRUE
      TRUE
      
    }, error = function(e) {
      parallel_failed_reason <<- conditionMessage(e)
      FALSE
    }, finally = {
      if (!is.null(cl)) {
        try(parallel::stopCluster(cl), silent = TRUE)
        try(foreach::registerDoSEQ(), silent = TRUE)
      }
    })
    
    if (!isTRUE(ok)) {
      message(sprintf("[cv][PAR] Parallel attempt FAILED -> falling back to sequential. Reason: %s",
                      parallel_failed_reason))
      did_parallel <- FALSE
      .run_sequential_cv(reason = parallel_failed_reason)
    } else {
      .log_line("[cv][PAR] Parallel attempt succeeded (cv_mat filled from parallel results).")
    }
    
  } else {
    .run_sequential_cv(reason = "parallel disabled / unavailable")
  }
  
  # ---------- NA check and refill sequentially ----------
  n_na <- sum(!is.finite(cv_mat))
  .log_line(sprintf("[cv] post-run | non-finite cv_mat cells = %d (out of %d)", n_na, length(cv_mat)))
  
  na_idx <- which(!is.finite(cv_mat), arr.ind = TRUE)
  if (nrow(na_idx) > 0L) {
    message(sprintf("[cv][REFILL] Refilling %d missing (fold, λ) cells SEQUENTIALLY.", nrow(na_idx)))
    for (ii in seq_len(nrow(na_idx))) {
      li <- na_idx[ii, 1]; k <- na_idx[ii, 2]
      fs <- fold_stats[[k]]
      lam <- lambdas[li]
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
  
  .log_line(sprintf("[cv] refit | lambda_min = %g (refitting on full data)", lambda_min))
  fit_min <- .call_admm(S_full, S0_full, lambda_min)
  
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
    relative_loss = relative_loss,
    did_parallel = did_parallel,
    parallel_failed_reason = parallel_failed_reason
  )
  class(out) <- "cv_admm_sgca_C"
  out
}




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
    nb_cores     = NULL,     # if NULL -> min(L*K, detectCores()-1)
    blas_threads = 1,        # threads per worker (1 avoids oversubscription)
    rng_seed     = NULL,     # optional: reproducible per-(fold,lambda) fits
    ...
) {
  penalty   <- match.arg(penalty)
  penalize  <- match.arg(penalize)
  loss_part <- match.arg(loss_part)
  
  .log_line <- function(...) if (isTRUE(verbose)) cat(..., "\n")
  
  n <- nrow(X); p <- ncol(X)
  stopifnot(sum(p_list) == p)
  
  .log_line(sprintf("[cv] start | pid=%s | parallel requested=%s", Sys.getpid(), parallel))
  
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
    .log_line(sprintf("[cv] config | K=%d folds | L=%d lambdas | n=%d | p=%d", K, L, n, p))
    .log_line(sprintf("[cv] config | penalty=%s | penalize=%s | loss_part=%s | relative_loss=%s",
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
  if (verbose) .log_line("[cv] precompute | building per-fold second moments (down-dated)...")
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
  if (verbose) .log_line("[cv] precompute | done")
  
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
  
  # ---------- sequential runner (fallback / refill) ----------
  .run_sequential_cv <- function(reason = NULL) {
    if (!is.null(reason)) message(sprintf("[cv][SEQ] Running sequentially because: %s", reason))
    else .log_line("[cv][SEQ] Running CV sequentially")
    
    for (k in seq_len(K)) {
      if (verbose) .log_line(sprintf("[cv][SEQ] Fold %d/%d", k, K))
      fs <- fold_stats[[k]]
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
      }
    }
    invisible(NULL)
  }
  
  # ============================
  # PARALLEL over (lambda, fold) PAIRS  (L*K tasks)
  # ============================
  parallel_ok_to_try <- isTRUE(parallel) &&
    requireNamespace("doParallel", quietly = TRUE) &&
    requireNamespace("foreach", quietly = TRUE) &&
    (L > 1L || K > 1L)
  
  did_parallel <- FALSE
  parallel_failed_reason <- NULL
  
  if (parallel_ok_to_try) {
    cl <- NULL
    
    ok <- tryCatch({
      
      nphys <- tryCatch(parallel::detectCores(logical = TRUE), error = function(e) 1L)
      n_tasks <- L * K
      
      cores <- if (is.null(nb_cores)) max(1L, min(n_tasks, nphys - 1L))
      else                   max(1L, min(n_tasks, as.integer(nb_cores)))
      
      .log_line(sprintf("[cv][PAR] mode=grid (lambda×fold) | tasks=%d (=L*K) | detected_cores=%d | using_cores=%d",
                        n_tasks, nphys, cores))
      if (cores <= 1L) stop("only 1 core available/requested (cores<=1)")
      
      cl <- parallel::makeCluster(cores)
      doParallel::registerDoParallel(cl)
      
      .log_line(sprintf("[cv][PAR] foreach backend: %s | workers: %d",
                        foreach::getDoParName(), foreach::getDoParWorkers()))
      .log_line(sprintf("[cv][PAR] master pid=%s", Sys.getpid()))
      
      # Load lightweight deps on workers
      parallel::clusterEvalQ(cl, {
        suppressPackageStartupMessages({
          library(stats)
          library(Matrix)
          library(pracma)  # matmul
          if (requireNamespace("rlang", quietly = TRUE)) library(rlang)
          if (requireNamespace("purrr", quietly = TRUE)) library(purrr)
        })
        NULL
      })
      
      # ---- Robust export (INCLUDING dot-names) ----
      sym_required <- c("admm_sgca", ".cv_loss_C")
      sym_optional <- c("soft_threshold", "matmul",
                        ".mask_from_part", ".prox_l21_rows", ".prox_l21_groups",
                        ".with_browser_suppressed", "%||%")
      sym_to_ship  <- unique(c(sym_required, sym_optional))
      
      export_env <- new.env(parent = emptyenv())
      found <- logical(length(sym_to_ship)); names(found) <- sym_to_ship
      
      for (nm in sym_to_ship) {
        if (exists(nm, inherits = TRUE)) {
          assign(nm, get(nm, inherits = TRUE), envir = export_env)
          found[nm] <- TRUE
        } else {
          found[nm] <- FALSE
        }
      }
      
      missing_required <- sym_required[!found[sym_required]]
      if (length(missing_required) > 0L) {
        stop(sprintf("Required symbol(s) not found on master: %s",
                     paste(missing_required, collapse = ", ")))
      }
      
      shipped <- ls(export_env, all.names = TRUE)   # <-- CRITICAL for dot-functions
      parallel::clusterExport(cl, varlist = shipped, envir = export_env)
      .log_line(sprintf("[cv][PAR] exported symbols to workers: %s", paste(shipped, collapse = ", ")))
      
      # Export local objects needed in workers
      parallel::clusterExport(
        cl,
        varlist = c("fold_stats","K","L","lambdas","r","p_list",
                    "part_eff","loss_weight","relative_loss",
                    "penalty","penalize","valid_dots","rng_seed"),
        envir = environment()
      )
      .log_line("[cv][PAR] exported local objects to workers")
      
      # Verify required symbol exists on each worker
      chk <- parallel::clusterCall(cl, function() exists(".cv_loss_C", inherits = TRUE))
      .log_line(sprintf("[cv][PAR] worker check: .cv_loss_C exists on %d/%d workers", sum(unlist(chk)), length(chk)))
      if (!all(unlist(chk))) stop("At least one worker missing .cv_loss_C after export")
      
      # Thread controls
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
      .log_line(sprintf("[cv][PAR] set BLAS/OMP threads per worker to %d", as.integer(blas_threads)))
      
      .log_line("[cv][PAR] launching foreach %dopar% over (lambda,fold) grid ...")
      
      # ---- One task per (lambda, fold) pair ----
      # Use ii -> (li,k) mapping to avoid shipping a grid object.
      res_list <- foreach::foreach(
        ii = seq_len(n_tasks),
        .inorder = FALSE,
        .errorhandling = "pass",
        .options.snow = list(preschedule = FALSE, chunkSize = 1)  # better load-balancing for uneven runtimes
      ) %dopar% {
        
        # map ii -> (li,k) with ordering: ii = (li-1)*K + k
        li <- ((ii - 1L) %/% K) + 1L
        k  <- ((ii - 1L) %%  K) + 1L
        
        lam <- lambdas[li]
        fs  <- fold_stats[[k]]
        
        if (!is.null(rng_seed)) set.seed(as.integer(rng_seed + 1000L * k + li))
        
        pid  <- Sys.getpid()
        host <- tryCatch(Sys.info()[["nodename"]], error = function(e) NA_character_)
        err  <- NULL
        
        loss <- tryCatch({
          base_args <- list(
            Sigma   = fs$S_tr,
            Sigma0  = fs$S0_tr,
            lambda  = lam,
            r       = r,
            p_list  = p_list,
            penalty = penalty,
            penalize= penalize,
            verbose = FALSE
          )
          
          fit <- do.call(admm_sgca, c(base_args, valid_dots))
          
          .cv_loss_C(
            fit$C, fs$S_va, fs$S0_va,
            part     = part_eff,
            p_list   = p_list,
            weight   = loss_weight,
            relative = relative_loss
          )
        }, error = function(e) {
          err <<- conditionMessage(e)
          NA_real_
        })
        
        list(li = li, k = k, lambda = lam, pid = pid, host = host, loss = loss, err = err)
      }
      
      # Handle foreach-level errors
      err_obj <- vapply(res_list, inherits, logical(1), "error")
      if (any(err_obj)) {
        .log_line(sprintf("[cv][PAR] WARNING: %d/%d tasks returned a foreach error object.",
                          sum(err_obj), length(res_list)))
        .log_line(sprintf("[cv][PAR] Example foreach error: %s",
                          conditionMessage(res_list[[which(err_obj)[1]]])))
      }
      
      # Fill cv_mat
      for (x in res_list) {
        if (is.list(x) && !is.null(x$li) && !is.null(x$k) && !is.null(x$loss)) {
          cv_mat[x$li, x$k] <- x$loss
        }
      }
      
      # Report how parallel it really was
      ok_items <- vapply(res_list, function(x) is.list(x) && !is.null(x$pid), logical(1))
      pids_used <- unique(vapply(res_list[ok_items], function(x) x$pid, integer(1)))
      hosts_used <- unique(vapply(res_list[ok_items], function(x) as.character(x$host), character(1)))
      
      .log_line(sprintf("[cv][PAR] completed | unique worker PIDs observed: %s", paste(pids_used, collapse = ", ")))
      if (length(hosts_used)) .log_line(sprintf("[cv][PAR] worker hosts observed: %s", paste(hosts_used, collapse = ", ")))
      
      # If still many NA, print one example fold-level error
      any_task_err <- vapply(res_list, function(x) is.list(x) && !is.null(x$err) && nzchar(x$err), logical(1))
      if (any(any_task_err)) {
        idx <- which(any_task_err)[1]
        .log_line(sprintf("[cv][PAR] Example task error | lambda=%g fold=%d pid=%s: %s",
                          res_list[[idx]]$lambda, res_list[[idx]]$k, res_list[[idx]]$pid, res_list[[idx]]$err))
      }
      
      did_parallel <- TRUE
      TRUE
      
    }, error = function(e) {
      parallel_failed_reason <<- conditionMessage(e)
      FALSE
    }, finally = {
      if (!is.null(cl)) {
        try(parallel::stopCluster(cl), silent = TRUE)
        try(foreach::registerDoSEQ(), silent = TRUE)
      }
    })
    
    if (!isTRUE(ok)) {
      message(sprintf("[cv][PAR] Parallel attempt FAILED -> falling back to sequential. Reason: %s",
                      parallel_failed_reason))
      did_parallel <- FALSE
      .run_sequential_cv(reason = parallel_failed_reason)
    } else {
      .log_line("[cv][PAR] Parallel attempt succeeded (cv_mat filled from parallel results).")
    }
    
  } else {
    .run_sequential_cv(reason = "parallel disabled / unavailable")
  }
  
  # ---------- NA check and refill sequentially ----------
  n_na <- sum(!is.finite(cv_mat))
  .log_line(sprintf("[cv] post-run | non-finite cv_mat cells = %d (out of %d)", n_na, length(cv_mat)))
  
  na_idx <- which(!is.finite(cv_mat), arr.ind = TRUE)
  if (nrow(na_idx) > 0L) {
    message(sprintf("[cv][REFILL] Refilling %d missing (fold, λ) cells SEQUENTIALLY.", nrow(na_idx)))
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
  
  .log_line(sprintf("[cv] refit | lambda_min = %g (refitting on full data)", lambda_min))
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
    relative_loss = relative_loss,
    did_parallel = did_parallel,
    parallel_failed_reason = parallel_failed_reason
  )
  class(out) <- "cv_admm_sgca_C"
  out
}






