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
  S_full  <- 0.5 * (S_full + t(S_full))   # symmetrize
  S0_full <- bd_from_S(S_full, p_list)
  
  # ---------- Precompute per-fold second moments via down-dating ----------
  if (verbose) cat("Precomputing per-fold second-moment matrices (down-dated)...\n")
  fold_stats <- vector("list", K)
  for (k in seq_len(K)) {
    idx_val <- which(folds == k)
    Xval    <- X[idx_val, , drop = FALSE]
    n_val   <- nrow(Xval)
    n_tr    <- n_full - n_val
    
    M_val <- crossprod(Xval)          # X_val' X_val
    S_va  <- M_val / n_val
    S_va  <- 0.5 * (S_va + t(S_va))
    
    S_tr  <- (n_full * S_full - M_val) / n_tr
    S_tr  <- 0.5 * (S_tr + t(S_tr))
    
    S0_tr <- bd_from_S(S_tr, p_list)
    S0_va <- bd_from_S(S_va, p_list)
    
    fold_stats[[k]] <- list(
      S_tr  = S_tr,
      S0_tr = S0_tr,
      S_va  = S_va,
      S0_va = S0_va
    )
  }
  
  # ---------- Freeze solver & filter ... to supported args (with aliases) ----------
  admm_fun <- admm_sgca
  dots <- list(...)
  admm_formals <- tryCatch(setdiff(names(formals(admm_fun)), "..."),
                           error = function(e) character())
  valid_dots <- if (length(admm_formals)) dots[names(dots) %in% admm_formals] else list()
  # Map common aliases -> solver names if present
  alias_map <- c("max_iter" = "maxit", "rho" = "rho0", "adapt_rho" = "adapt")
  for (from in names(alias_map)) {
    to <- alias_map[[from]]
    if (!is.null(dots[[from]]) && (to %in% admm_formals) && !(from %in% admm_formals)) {
      valid_dots[[to]] <- dots[[from]]
    }
  }
  ignored_dots <- setdiff(names(dots), union(names(valid_dots), names(alias_map)))
  if (length(ignored_dots) && isTRUE(verbose)) {
    message("Ignoring unsupported admm_sgca() args: ", paste(ignored_dots, collapse = ", "))
  }
  
  # Optional browser suppressor for sequential calls
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
  
  # Identify the package that defines admm_sgca (for compiled code on workers)
  admm_pkg <- tryCatch({
    nm <- environmentName(environment(admm_fun))
    nm <- sub("^(package|namespace):", "", nm)
    if (identical(nm, "R_GlobalEnv")) NULL else nm
  }, error = function(e) NULL)
  
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
    
    # Initialize workers (pass args explicitly — no missing globals)
    parallel::clusterCall(
      cl,
      function(bt, admm_pkg) {
        if (!is.null(admm_pkg)) {
          try(suppressMessages(requireNamespace(admm_pkg, quietly = TRUE)), silent = TRUE)
        }
        if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
          RhpcBLASctl::blas_set_num_threads(as.integer(bt))
        }
        Sys.setenv(OMP_NUM_THREADS = as.character(bt),
                   MKL_NUM_THREADS = as.character(bt))
        NULL
      },
      as.integer(blas_threads),
      admm_pkg
    )
    
    # Bring %dopar% into scope
    `%dopar%` <- foreach::`%dopar%`
    
    # Export exact objects referenced in the body
    res_list <- foreach::foreach(
      fid = seq_len(K),
      .inorder = TRUE,
      .errorhandling = "pass",
      .export = c("fold_stats","lambdas","L","p_list","part_eff","loss_weight",
                  "relative_loss",".cv_loss_C","penalty","penalize","r",
                  "admm_fun","valid_dots","rng_seed","wb_fun")
    ) %dopar% {
      fs <- fold_stats[[fid]]
      losses <- rep(NA_real_, length(lambdas))
      
      if (!is.null(rng_seed)) set.seed(as.integer(rng_seed + 1000L * fid))
      
      for (li in seq_along(lambdas)) {
        lam <- lambdas[[li]]
        if (!is.null(rng_seed)) set.seed(as.integer(rng_seed + 1000L * fid + li))
        losses[li] <- tryCatch({
          base_args <- list(Sigma = fs$S_tr, Sigma0 = fs$S0_tr,
                            lambda = lam, r = r, p_list = p_list,
                            penalty = penalty, penalize = penalize,
                            verbose = FALSE)
          fit <- do.call(admm_fun, c(base_args, valid_dots))
          C_hat <- fit$C
          .cv_loss_C(C_hat, fs$S_va, fs$S0_va,
                     part = part_eff, p_list = p_list,
                     weight = loss_weight, relative = relative_loss)
        }, error = function(e) NA_real_)
      }
      losses
    }
    
    # Combine into L×K matrix (guard if a fold returned an error object)
    for (k in seq_len(K)) {
      col <- res_list[[k]]
      if (is.numeric(col) && length(col) == L) cv_mat[, k] <- col
      else cv_mat[, k] <- NA_real_
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
          fit  <- .call_admm(fs$S_tr, fs$S0_tr, lam)
          C_hat <- fit$C
          .cv_loss_C(C_hat, fs$S_va, fs$S0_va,
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
  
  # ---------- Refill missing (fold, λ) cells sequentially (fair CV) ----------
  na_idx <- which(!is.finite(cv_mat), arr.ind = TRUE)
  if (nrow(na_idx) > 0L) {
    if (verbose) cat(sprintf("Refilling %d missing (fold, λ) cells sequentially...\n", nrow(na_idx)))
    for (ii in seq_len(nrow(na_idx))) {
      li <- na_idx[ii, 1]; k <- na_idx[ii, 2]
      fs <- fold_stats[[k]]
      lam <- lambdas[li]
      if (!is.null(rng_seed)) set.seed(as.integer(rng_seed + 1000L * k + li))
      loss <- tryCatch({
        fit  <- .call_admm(fs$S_tr, fs$S0_tr, lam)
        C_hat <- fit$C
        .cv_loss_C(C_hat, fs$S_va, fs$S0_va,
                   part = part_eff, p_list = p_list,
                   weight = loss_weight, relative = relative_loss)
      }, error = function(e) NA_real_)
      cv_mat[li, k] <- loss
    }
  }
  
  # ---------- Aggregate CV: choose among λ with ALL K folds ----------
  row_ok <- rowSums(is.finite(cv_mat)) == K
  if (!any(row_ok)) {
    stop("No lambda was successfully evaluated on all folds. ",
         "Check solver errors or relax settings.")
  }
  
  cvm_sel  <- rowMeans(cv_mat[row_ok, , drop = FALSE], na.rm = FALSE)
  cvsd_sel <- apply(cv_mat[row_ok, , drop = FALSE], 1, function(x) {
    m <- length(x); if (m <= 1) return(NA_real_)
    stats::sd(x, na.rm = FALSE) / sqrt(m)
  })
  
  best_rel_idx <- which.min(cvm_sel)
  lambda_min   <- lambdas[row_ok][best_rel_idx]
  
  if (verbose) cat(sprintf("Refitting on full data at lambda_min = %g\n", lambda_min))
  fit_min <- .call_admm(S_full, S0_full, lambda_min)
  
  # For reporting, also compute cvm/cvsd over all rows (with na.rm=TRUE)
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


fit <- cv_admm_sgca_C(
  X, p_list, lambdas = 10^seq(-3, 0, length.out = 12), r = 1, K = 5,
  penalty = "l1", penalize = "all", loss_part = "all", relative_loss = TRUE,
  rho = 1, max_iter = 200, adapt_rho = TRUE,   # used only if admm_sgca supports them
  parallel = TRUE, nb_cores = 5, blas_threads = 1,
  verbose = TRUE
)
