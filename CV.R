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

subdistance <- function(A, B){
  svdresult = svd(t(A) %*% B);
  U = svdresult$u
  V = svdresult$v
  O = U %*% t(V);
  l = norm(A %*% O-B, type = c('F'));
  return(l)
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

# (Removed .choose_lambda_1se)

.mask_from_part <- function(p, part = c("all","offdiag","block_off"), p_list = NULL) {
  part <- match.arg(part)
  M <- matrix(TRUE, p, p)
  if (part == "offdiag") {
    diag(M) <- FALSE
  } else if (part == "block_off") {
    if (is.null(p_list)) stop("p_list is required when part='block_off'.")
    M[,] <- TRUE
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
  
  # Match training objective: ||Sigma_val - Sigma0_val * C * Sigma0_val||_F^2
  A <- (Sigma_val + t(Sigma_val)) / 2
  B <- Sigma0_val %*% ((C + t(C))/2) %*% Sigma0_val   # symmetrize C for robustness
  R <- A - B
  
  # Optional masking/weighting
  mask_from_part <- function(p, part, p_list) {
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
  M <- mask_from_part(p, part, p_list)
  W <- matrix(0, p, p); W[M] <- 1
  if (!is.null(weight)) {
    if (length(weight) == 1) weight <- matrix(weight, p, p)
    stopifnot(all(dim(weight) == c(p, p)))
    W <- W * weight
  }
  
  num <- sum((R * R) * W)
  if (!relative) return(num)
  
  # Use ||Sigma_val||_F^2 (under the same mask/weights) as the scale
  den <- sum((A * A) * W)
  if (den <= .Machine$double.eps) return(NA_real_)
  num / den
}

##################ADMM function############################################
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
  SVD=svd(Sigma0_sqrt%*%C%*%Sigma0_sqrt)$u
  lam_inv=diag(1/svd(Sigma0_sqrt%*%C%*%Sigma0_sqrt)$d)
  U_canon=C%*%Sigma0_sqrt%*%SVD%*%lam_inv
  
  
  
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



# ---------- main: CV scoring by C ----------
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


pp <- c(10,10,10);
p <- sum(pp)
s  <- c(1:3);
r <- 1;

p<-sum(pp)
rho<-1

lambda_values <- c(0, 1e-5,1e-4,1e-3,1e-2, 0.1, 1,10,100,1000,1e+4,1e+5)

p_list <- c(10,10,10)
p <- sum(p_list)
C0_Cmat=list()
cv_fit=list()
admm=list(); oracle=list();oracle2=list();naive=list();oracle_c_loss=c()
oracle$time=c()
oracle$loss=c()
naive$time=c()
naive$loss=c()
oracle2$time=c()
oracle2$loss=c()
admm_cv_time=c()
N=c(30,1000,10000,50000,100000)
for(j in 1:5){
  n=N[j]
  #for(i in 1:10){
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

  # Generate data for SGCA
  X <-mvrnorm(n, rep(0, p) , Sigma)
  S <- t(X)%*% X / n
  sigma0hat <- S * Mask
  
  C_0 <- solve(Sigma0) %*% Sigma %*% solve(Sigma0)
  can_vec=pracma::sqrtm(Sigma0)$Binv %*% svd(pracma::sqrtm(Sigma0)$B %*% C_0 %*%pracma::sqrtm(Sigma0)$B)$u
  a=can_vec[,1]
  C0_Cmat[[j]] <- C_0
  result <- tryCatch(
    {
      
      
      start_time <- proc.time()
      cv_fit[[j]]<-cv_admm_sgca_C(
         X, p_list, lambda_values, r = 1, K = 5,
         penalty = "l1", penalize = "all",
         loss_part = "all",         # <-- matches the stated objective
         relative_loss = TRUE,
         rho = 1, max_iter = 4000, adapt_rho = TRUE, verbose = FALSE
       )
      end_time <- proc.time()
      admm$time[j] <- end_time-start_time  # Store NA for failed iterations
      admm$loss[j] <- subdistance(cv_fit[[j]]$fit_min$U[,1], a)
      
      output  # Return the output to be used in the result variable
    },
    error = function(e) {
      cat("Error in iteration", j, "- skipping\n")
      admm$time[j] <- NA  # Store NA for failed iterations
      admm$loss[j] <- NA
      return(NULL)
    }
    
  )
  
  result <- tryCatch(
    {
      #result <- rgcca(blocks,method="rgcca")
      start_time <- proc.time()
      newid=c(s,pp[1]+s,pp[1]+pp[2]+s)
      
      sqrt_m_sigma = pracma::sqrtm(sigma0hat[newid,newid])
      result = svd(sqrt_m_sigma$Binv %*% S[newid,newid] %*% sqrt_m_sigma$Binv )
      evalues <- result$d
      evectors <-result$u
      u_oracle <- (sqrt_m_sigma$Binv %*% as.matrix(evectors))[,1]
      C_oracle = sqrt_m_sigma$Binv %*% as.matrix(evectors) %*% diag(evalues) %*% t(as.matrix(evectors)) %*% sqrt_m_sigma$Binv
      end_time <- proc.time()
      oracle2$time[j]=end_time - start_time
      oracle2$loss[j]=subdistance(u_oracle,a[newid])
      oracle_c_loss[j]=norm(C_oracle -C_0, "F")
      output  # Return the output to be used in the result variable
    },
    error = function(e) {
      cat("Error in iteration", j, "- skipping\n")
      oracle2$time[j] <- NA  # Store NA for failed iterations
      oracle2$loss[j] <- NA
      return(NULL)
    }
    
  )
  print((j))
  
}



