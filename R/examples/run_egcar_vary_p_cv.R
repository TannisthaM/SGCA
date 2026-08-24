#!/usr/bin/env Rscript

# =============================================================================
# EGCAR simulation with three views, pairwise ADMM, theoretical penalty,
# five-fold held-out-Rayleigh cross-validation, and an oracle-support estimator.
#
# Population model (Section 6.1 of the manuscript):
#   Sigma_kk[i,j] = alpha^|i-j|,
#   U_k = G_k (G_k^T Sigma_kk G_k)^(-1/2),
#   Sigma_kl = signal_lambda Sigma_kk U_k U_l^T Sigma_ll,
#   C0_kl = Sigma_kk^(-1) Sigma_kl Sigma_ll^(-1)
#         = signal_lambda U_k U_l^T.
#
# Sparse pairwise estimator:
#   argmin_C 0.5 || n^(-1) X_k C X_l^T - I_n ||_F^2 + rho ||C||_1.
#
# Methods:
#   1. Oracle: true active rows are supplied; values are estimated from data.
#   2. EGCAR-theory: rho = sqrt(log(p_total)/n).
#   3. EGCAR-CV: rho chosen from 10^(-5),...,10^4 by five-fold validation
#      using the held-out generalized Rayleigh score in equation (40).
#
# The row threshold is fixed at tau = 1e-4 for every EGCAR fit.
# For very large p, candidate-row screening is fitted inside each training fold.
# =============================================================================

options(stringsAsFactors = FALSE, warn = 1)

EXPERIMENT <- "vary_p_cv"
SCRIPT_NAME <- "run_egcar_vary_p_cv.R"
X_LABEL <- "Number of covariates per block p"
TAU <- 1e-4
RHO_GRID <- 10^(-5:4)
RANKS <- c(1L, 2L, 5L)
N_REPS <- 10L

make_config_table <- function() {
  rows <- list()
  id <- 1L
  for (rank in RANKS) {
    for (p in c(30L, 50L, 100L, 500L, 1000L, 5000L)) {
      rows[[id]] <- data.frame(
        config_id = id,
        panel = "n150",
        n = 150L,
        p_per_block = p,
        rank = rank,
        stringsAsFactors = FALSE
      )
      id <- id + 1L
    }
  }
  do.call(rbind, rows)
}
CONFIGS <- make_config_table()

get_env_num <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(default)
  ans <- suppressWarnings(as.numeric(value))
  if (!is.finite(ans)) stop(sprintf("%s must be numeric", name))
  ans
}
get_env_int <- function(name, default) as.integer(round(get_env_num(name, default)))

CFG <- list(
  K = 3L,
  support_size = get_env_int("EGCAR_SUPPORT_SIZE", 5L),
  signal_lambda = get_env_num("EGCAR_SIGNAL_LAMBDA", 0.80),
  alpha = get_env_num("EGCAR_TOEPLITZ_ALPHA", 0.30),
  theory_rho_constant = get_env_num("EGCAR_THEORY_RHO_CONSTANT", 1.0),
  cv_folds = get_env_int("EGCAR_CV_FOLDS", 5L),
  mu = get_env_num("EGCAR_ADMM_MU", 1.0),
  max_iter = get_env_int("EGCAR_ADMM_MAX_ITER", 2000L),
  abs_tol = get_env_num("EGCAR_ADMM_ABS_TOL", 1e-5),
  rel_tol = get_env_num("EGCAR_ADMM_REL_TOL", 1e-4),
  covariance_ridge = get_env_num("EGCAR_COVARIANCE_RIDGE", 1e-4),
  oracle_ridge = get_env_num("EGCAR_ORACLE_RIDGE", 1e-8),
  score_floor = get_env_num("EGCAR_SCORE_FLOOR", 1e-10),
  screen_cap = get_env_int("EGCAR_SCREEN_CAP", 200L),
  screen_support_multiplier = get_env_int("EGCAR_SCREEN_SUPPORT_MULTIPLIER", 8L),
  screen_sqrt_multiplier = get_env_num("EGCAR_SCREEN_SQRT_MULTIPLIER", 2.5),
  heatmap_max_per_block = get_env_int("EGCAR_HEATMAP_MAX_PER_BLOCK", 30L),
  base_seed = get_env_int("EGCAR_BASE_SEED", 912407L)
)

if (CFG$K != 3L) stop("This simulation is written for exactly three blocks")
if (CFG$support_size < max(RANKS)) stop("EGCAR_SUPPORT_SIZE must be at least 5")
if (!(CFG$signal_lambda > 0 && CFG$signal_lambda < 1)) {
  stop("EGCAR_SIGNAL_LAMBDA must lie in (0,1)")
}
if (abs(CFG$alpha) >= 1) stop("EGCAR_TOEPLITZ_ALPHA must have absolute value < 1")
if (CFG$mu <= 0) stop("EGCAR_ADMM_MU must be positive")
if (CFG$cv_folds < 2L) stop("EGCAR_CV_FOLDS must be at least 2")

usage <- function() {
  cat(sprintf("Worker: Rscript %s <config_id> <rep_id> <outdir>\n", SCRIPT_NAME))
  cat(sprintf("Aggregate: Rscript %s aggregate 0 <outdir>\n", SCRIPT_NAME))
}

available_cores <- function() {
  z <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1")))
  if (!is.finite(z) || z < 1L) z <- 1L
  z
}

make_dirs <- function(outdir) {
  for (d in c("metrics", "cv", "snapshots", "plots", "config")) {
    dir.create(file.path(outdir, d), recursive = TRUE, showWarnings = FALSE)
  }
}

atomic_csv <- function(x, path, row.names = FALSE) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp.", Sys.getpid(), ".", sample.int(1e9, 1L))
  write.csv(x, tmp, row.names = row.names)
  if (!file.rename(tmp, path)) {
    unlink(tmp)
    stop("Could not atomically write ", path)
  }
}

atomic_rds <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp.", Sys.getpid(), ".", sample.int(1e9, 1L))
  saveRDS(x, tmp, compress = TRUE)
  if (!file.rename(tmp, path)) {
    unlink(tmp)
    stop("Could not atomically write ", path)
  }
}

fnorm <- function(A) sqrt(sum(A * A))
symm <- function(A) (A + t(A)) / 2
soft_threshold <- function(A, threshold) sign(A) * pmax(abs(A) - threshold, 0)

inverse_sqrt <- function(A, floor = 1e-12) {
  ee <- eigen(symm(A), symmetric = TRUE)
  if (min(ee$values) <= floor) stop("Matrix is not positive definite in inverse_sqrt")
  ee$vectors %*% (diag(1 / sqrt(ee$values), nrow = length(ee$values)) %*% t(ee$vectors))
}

matrix_roots <- function(A, ridge = 0) {
  d <- nrow(A)
  if (d == 0L) return(list(sqrt = matrix(0, 0, 0), invsqrt = matrix(0, 0, 0)))
  scale <- max(1, mean(diag(A)))
  Areg <- symm(A) + diag(ridge * scale, d)
  ee <- eigen(Areg, symmetric = TRUE)
  vals <- pmax(ee$values, .Machine$double.eps)
  list(
    sqrt = ee$vectors %*% (diag(sqrt(vals), d) %*% t(ee$vectors)),
    invsqrt = ee$vectors %*% (diag(1 / sqrt(vals), d) %*% t(ee$vectors))
  )
}

toeplitz_submatrix <- function(rows, cols, alpha) {
  if (!length(rows) || !length(cols)) return(matrix(0, length(rows), length(cols)))
  alpha ^ abs(outer(rows, cols, "-"))
}

# ------------------------ population and data generation ----------------------

make_population <- function(p_vec, rank, rep_id) {
  # The seed does not depend on n or p. Thus the active signal is identical
  # across the dimension/sample-size curve for a fixed rank and replicate.
  population_seed <- CFG$base_seed + 100000L * rank + rep_id
  set.seed(population_seed)

  supports <- U_active <- SigmaU <- vector("list", CFG$K)
  support <- seq_len(CFG$support_size)
  Sigma_ss <- toeplitz_submatrix(support, support, CFG$alpha)

  for (k in seq_len(CFG$K)) {
    G <- matrix(rnorm(CFG$support_size * rank), CFG$support_size, rank)
    U <- G %*% inverse_sqrt(crossprod(G, Sigma_ss %*% G))
    check <- crossprod(U, Sigma_ss %*% U)
    if (max(abs(check - diag(rank))) > 1e-7) stop("Population U normalization failed")

    supports[[k]] <- support
    U_active[[k]] <- U
    SigmaU[[k]] <- toeplitz_submatrix(seq_len(p_vec[k]), support, CFG$alpha) %*% U
  }

  C0 <- list()
  for (k in 1:(CFG$K - 1L)) {
    for (ell in (k + 1L):CFG$K) {
      key <- sprintf("%d_%d", k, ell)
      C0[[key]] <- list(
        rows = supports[[k]], cols = supports[[ell]],
        value = CFG$signal_lambda * U_active[[k]] %*% t(U_active[[ell]]),
        k = k, ell = ell
      )
    }
  }

  list(
    p_vec = as.integer(p_vec), rank = rank, supports = supports,
    U_active = U_active, SigmaU = SigmaU, C0 = C0,
    population_seed = population_seed
  )
}

ar1_gaussian <- function(n, p, alpha) {
  E <- matrix(rnorm(n * p), n, p)
  if (p > 1L && alpha != 0) {
    innovation_sd <- sqrt(1 - alpha^2)
    for (j in 2:p) E[, j] <- alpha * E[, j - 1L] + innovation_sd * E[, j]
  }
  E
}

sample_views <- function(pop, n, config_id, rep_id) {
  data_seed <- CFG$base_seed + 10000000L + 10000L * rep_id + config_id
  set.seed(data_seed)
  z <- matrix(rnorm(n * pop$rank), n, pop$rank)
  adjustment <- sqrt(1 - CFG$signal_lambda) - 1
  X <- vector("list", CFG$K)

  for (k in seq_len(CFG$K)) {
    Y <- ar1_gaussian(n, pop$p_vec[k], CFG$alpha)
    factor_scores <- Y[, pop$supports[[k]], drop = FALSE] %*% pop$U_active[[k]]
    X[[k]] <- Y +
      (adjustment * factor_scores + sqrt(CFG$signal_lambda) * z) %*%
      t(pop$SigmaU[[k]])
  }

  list(X = X, data_seed = data_seed)
}

center_full <- function(X) {
  lapply(X, function(A) sweep(A, 2L, colMeans(A), "-"))
}

center_train_validation <- function(X, train_idx, validation_idx) {
  means <- lapply(X, function(A) colMeans(A[train_idx, , drop = FALSE]))
  train <- validation <- vector("list", length(X))
  for (k in seq_along(X)) {
    train[[k]] <- sweep(X[[k]][train_idx, , drop = FALSE], 2L, means[[k]], "-")
    validation[[k]] <- sweep(X[[k]][validation_idx, , drop = FALSE], 2L, means[[k]], "-")
  }
  list(train = train, validation = validation)
}

# ------------------------- candidate-row screening ----------------------------

screen_rows <- function(X, true_supports = NULL) {
  start <- proc.time()[3]
  n <- nrow(X[[1]])
  p_vec <- vapply(X, ncol, integer(1))

  budgets <- vapply(p_vec, function(p) {
    rank_budget <- max(1L, n - 1L)
    root_budget <- as.integer(ceiling(CFG$screen_sqrt_multiplier * sqrt(p)))
    support_budget <- CFG$support_size * CFG$screen_support_multiplier
    min(p, max(support_budget, min(CFG$screen_cap, max(rank_budget, root_budget))))
  }, integer(1))

  if (all(budgets >= p_vec)) {
    sets <- lapply(p_vec, seq_len)
    scores <- lapply(p_vec, function(p) rep(NA_real_, p))
  } else {
    grams <- lapply(X, tcrossprod)
    pair_scores <- vector("list", CFG$K)
    aggregate_scores <- lapply(p_vec, numeric)
    for (k in seq_len(CFG$K)) pair_scores[[k]] <- vector("list", CFG$K)

    for (k in seq_len(CFG$K)) {
      for (ell in setdiff(seq_len(CFG$K), k)) {
        GX <- grams[[ell]] %*% X[[k]]
        score <- pmax(colSums(X[[k]] * GX) / n^2, 0)
        pair_scores[[k]][[ell]] <- score
        aggregate_scores[[k]] <- aggregate_scores[[k]] + score
      }
    }

    sets <- vector("list", CFG$K)
    for (k in seq_len(CFG$K)) {
      m <- budgets[k]
      if (m >= p_vec[k]) {
        sets[[k]] <- seq_len(p_vec[k])
        next
      }
      other_views <- setdiff(seq_len(CFG$K), k)
      per_pair <- max(1L, floor(m / length(other_views)))
      candidate <- integer(0)
      for (ell in other_views) {
        candidate <- union(
          candidate,
          head(order(pair_scores[[k]][[ell]], decreasing = TRUE), per_pair)
        )
      }
      if (length(candidate) > m) {
        ord <- order(aggregate_scores[[k]][candidate], decreasing = TRUE)
        candidate <- candidate[ord[seq_len(m)]]
      } else if (length(candidate) < m) {
        ord <- order(aggregate_scores[[k]], decreasing = TRUE)
        filler <- ord[!(ord %in% candidate)]
        candidate <- c(candidate, head(filler, m - length(candidate)))
      }
      sets[[k]] <- sort(candidate)
    }
    scores <- aggregate_scores
  }

  recall <- NA_real_
  if (!is.null(true_supports)) {
    total <- sum(vapply(true_supports, length, integer(1)))
    retained <- sum(vapply(seq_len(CFG$K), function(k) {
      sum(true_supports[[k]] %in% sets[[k]])
    }, integer(1)))
    recall <- retained / total
  }

  list(
    sets = sets, scores = scores, budgets = budgets, recall = recall,
    seconds = unname(proc.time()[3] - start)
  )
}

# ------------------------------ ADMM solver -----------------------------------

covariance_eigensystem <- function(X) {
  n <- nrow(X)
  p <- ncol(X)
  if (p > n) {
    # Thin SVD implementation of Section 4.2: Sigmahat = V diag(d^2/n) V^T.
    sv <- svd(X, nu = 0L, nv = min(n, p))
    values <- pmax(sv$d^2 / n, 0)
    tolerance <- max(values, 0) * max(n, p) * .Machine$double.eps
    keep <- which(values > tolerance)
    if (!length(keep)) {
      return(list(vectors = matrix(0, p, 0), values = numeric(0), thin = TRUE))
    }
    return(list(
      vectors = sv$v[, keep, drop = FALSE],
      values = values[keep],
      thin = length(keep) < p
    ))
  }

  S <- crossprod(X) / n
  ee <- eigen(symm(S), symmetric = TRUE)
  list(vectors = ee$vectors, values = pmax(ee$values, 0), thin = FALSE)
}

make_pair_stats <- function(Xk, Xell) {
  n <- nrow(Xk)
  Skell <- crossprod(Xk, Xell) / n
  ek <- covariance_eigensystem(Xk)
  eell <- covariance_eigensystem(Xell)
  list(
    n = n, pk = ncol(Xk), pell = ncol(Xell),
    Skell = Skell,
    Qk = ek$vectors, Qell = eell$vectors,
    dk = ek$values, dell = eell$values,
    use_thin = ek$thin || eell$thin
  )
}

pair_admm <- function(stats, rho, warm_Z = NULL) {
  start <- proc.time()[3]
  pk <- stats$pk
  pell <- stats$pell

  # KKT shortcut: C=0 is optimal whenever ||Sigmahat_kl||_max <= rho.
  if (max(abs(stats$Skell)) <= rho) {
    return(list(
      C = matrix(0, pk, pell), iterations = 0L, converged = TRUE,
      primal = 0, dual = 0, seconds = unname(proc.time()[3] - start)
    ))
  }

  covariance_products <- outer(stats$dk, stats$dell, "*")
  denom <- covariance_products + CFG$mu
  if (is.null(warm_Z) || !all(dim(warm_Z) == c(pk, pell))) {
    Z <- matrix(0, pk, pell)
  } else {
    Z <- warm_Z
  }
  H <- matrix(0, pk, pell)
  C <- Z
  converged <- FALSE
  primal <- dual <- NA_real_

  for (iter in seq_len(CFG$max_iter)) {
    G <- stats$Skell + CFG$mu * (Z - H)
    if (isTRUE(stats$use_thin)) {
      # Exact low-rank-plus-ridge solution:
      # C = G/mu - V_k [{d_k d_l A}/{mu(mu+d_k d_l)}] V_l^T,
      # A = V_k^T G V_l. This avoids completing V_k,V_l to p x p bases.
      if (!length(stats$dk) || !length(stats$dell)) {
        C <- G / CFG$mu
      } else {
        A <- crossprod(stats$Qk, G %*% stats$Qell)
        correction <- (covariance_products * A) / (CFG$mu * denom)
        C <- G / CFG$mu - stats$Qk %*% correction %*% t(stats$Qell)
      }
    } else {
      transformed <- crossprod(stats$Qk, G %*% stats$Qell) / denom
      C <- stats$Qk %*% transformed %*% t(stats$Qell)
    }

    Z_old <- Z
    Z <- soft_threshold(C + H, rho / CFG$mu)
    H <- H + C - Z

    primal <- fnorm(C - Z)
    dual <- CFG$mu * fnorm(Z - Z_old)
    eps_primal <- sqrt(pk * pell) * CFG$abs_tol +
      CFG$rel_tol * max(fnorm(C), fnorm(Z))
    eps_dual <- sqrt(pk * pell) * CFG$abs_tol +
      CFG$rel_tol * CFG$mu * fnorm(H)

    if (primal <= eps_primal && dual <= eps_dual) {
      converged <- TRUE
      break
    }
  }

  list(
    C = Z, iterations = as.integer(iter), converged = converged,
    primal = primal, dual = dual,
    seconds = unname(proc.time()[3] - start)
  )
}

oracle_pair <- function(Xk, Xell) {
  start <- proc.time()[3]
  n <- nrow(Xk)
  Sk <- crossprod(Xk) / n
  Sell <- crossprod(Xell) / n
  Skell <- crossprod(Xk, Xell) / n
  ridge_k <- CFG$oracle_ridge * max(1, mean(diag(Sk)))
  ridge_ell <- CFG$oracle_ridge * max(1, mean(diag(Sell)))
  left <- solve(Sk + diag(ridge_k, nrow(Sk)), Skell)
  C <- t(solve(Sell + diag(ridge_ell, nrow(Sell)), t(left)))
  list(C = C, seconds = unname(proc.time()[3] - start))
}

block_repr <- function(rows, cols, value, k, ell) {
  list(
    rows = as.integer(rows), cols = as.integer(cols), value = value,
    k = as.integer(k), ell = as.integer(ell)
  )
}

extract_block <- function(R, rows, cols) {
  ans <- matrix(0, length(rows), length(cols))
  if (!length(rows) || !length(cols)) return(ans)
  ir <- match(rows, R$rows)
  jc <- match(cols, R$cols)
  rr <- which(!is.na(ir))
  cc <- which(!is.na(jc))
  if (length(rr) && length(cc)) {
    ans[rr, cc] <- R$value[ir[rr], jc[cc], drop = FALSE]
  }
  ans
}

row_norms_from_blocks <- function(blocks, p_vec) {
  accum <- lapply(p_vec, numeric)
  for (R in blocks) {
    accum[[R$k]][R$rows] <- accum[[R$k]][R$rows] + rowSums(R$value^2)
    accum[[R$ell]][R$cols] <- accum[[R$ell]][R$cols] + colSums(R$value^2)
  }
  lapply(accum, sqrt)
}

# -------------------------- localized GCA loading -----------------------------

top_algebraic_eigen <- function(A, rank) {
  d <- nrow(A)
  if (d < rank) stop("Matrix dimension is smaller than requested rank")

  use_rspectra <- d > max(250L, 5L * rank) &&
    requireNamespace("RSpectra", quietly = TRUE) && rank < d

  if (use_rspectra) {
    ans <- tryCatch(
      RSpectra::eigs_sym(A, k = rank, which = "LA", opts = list(tol = 1e-8, maxitr = 10000L)),
      error = function(e) NULL
    )
    if (!is.null(ans)) {
      ord <- order(ans$values, decreasing = TRUE)
      return(list(values = ans$values[ord], vectors = ans$vectors[, ord, drop = FALSE]))
    }
  }

  ee <- eigen(symm(A), symmetric = TRUE)
  list(values = ee$values[seq_len(rank)], vectors = ee$vectors[, seq_len(rank), drop = FALSE])
}

fit_global_loading <- function(X, blocks, rank, forced_sets = NULL) {
  p_vec <- vapply(X, ncol, integer(1))
  if (is.null(forced_sets)) {
    norms <- row_norms_from_blocks(blocks, p_vec)
    sets <- lapply(norms, function(z) which(z > TAU))
  } else {
    sets <- lapply(forced_sets, as.integer)
    norms <- row_norms_from_blocks(blocks, p_vec)
  }

  sizes <- vapply(sets, length, integer(1))
  total_selected <- sum(sizes)
  if (total_selected < rank || total_selected == 0L) {
    return(list(valid = FALSE, sets = sets, L = NULL, row_norms = norms,
                eigenvalues = rep(NA_real_, rank), selected_total = total_selected))
  }

  roots <- vector("list", CFG$K)
  for (k in seq_len(CFG$K)) {
    if (!sizes[k]) {
      roots[[k]] <- matrix_roots(matrix(0, 0, 0), CFG$covariance_ridge)
    } else {
      S <- crossprod(X[[k]][, sets[[k]], drop = FALSE]) / nrow(X[[k]])
      roots[[k]] <- matrix_roots(S, CFG$covariance_ridge)
    }
  }

  offsets <- c(0L, cumsum(sizes))
  Rmat <- matrix(0, total_selected, total_selected)
  for (key in names(blocks)) {
    B <- blocks[[key]]
    k <- B$k
    ell <- B$ell
    if (!sizes[k] || !sizes[ell]) next
    Ckl <- extract_block(B, sets[[k]], sets[[ell]])
    Rkl <- roots[[k]]$sqrt %*% Ckl %*% roots[[ell]]$sqrt
    rr <- (offsets[k] + 1L):offsets[k + 1L]
    cc <- (offsets[ell] + 1L):offsets[ell + 1L]
    Rmat[rr, cc] <- Rkl
    Rmat[cc, rr] <- t(Rkl)
  }

  ee <- top_algebraic_eigen(symm(Rmat), rank)
  L <- vector("list", CFG$K)
  for (k in seq_len(CFG$K)) {
    if (!sizes[k]) {
      L[[k]] <- matrix(0, 0, rank)
    } else {
      rr <- (offsets[k] + 1L):offsets[k + 1L]
      L[[k]] <- roots[[k]]$invsqrt %*% ee$vectors[rr, , drop = FALSE]
    }
  }

  list(
    valid = TRUE, sets = sets, L = L, row_norms = norms,
    eigenvalues = ee$values, selected_total = total_selected
  )
}

heldout_rayleigh_score <- function(loading, X_validation) {
  if (is.null(loading) || !isTRUE(loading$valid)) return(NA_real_)
  largest_block <- which.max(vapply(loading$L, nrow, integer(1)))
  rank <- ncol(loading$L[[largest_block]])
  n_val <- nrow(X_validation[[1]])
  Y <- vector("list", CFG$K)

  for (k in seq_len(CFG$K)) {
    if (!length(loading$sets[[k]])) {
      Y[[k]] <- matrix(0, n_val, rank)
    } else {
      Y[[k]] <- X_validation[[k]][, loading$sets[[k]], drop = FALSE] %*% loading$L[[k]]
    }
  }

  Ysum <- Reduce("+", Y)
  A <- crossprod(Ysum) / n_val
  Q <- Reduce("+", lapply(Y, function(Z) crossprod(Z) / n_val))
  q_eig <- eigen(symm(Q), symmetric = TRUE, only.values = TRUE)$values
  if (min(q_eig) <= CFG$score_floor || any(!is.finite(q_eig))) return(NA_real_)

  score <- sum(diag(solve(Q, A)))
  if (is.finite(score)) score else NA_real_
}

population_subspace_distance <- function(loading, pop) {
  if (is.null(loading) || !isTRUE(loading$valid)) return(NA_real_)
  rank <- pop$rank
  A <- matrix(0, rank, rank)
  cross <- matrix(0, rank, rank)

  for (k in seq_len(CFG$K)) {
    S_hat <- loading$sets[[k]]
    if (!length(S_hat)) next
    Lk <- loading$L[[k]]
    Sigma_hat_hat <- toeplitz_submatrix(S_hat, S_hat, CFG$alpha)
    Sigma_hat_true <- toeplitz_submatrix(S_hat, pop$supports[[k]], CFG$alpha)
    A <- A + crossprod(Lk, Sigma_hat_hat %*% Lk)
    cross <- cross + crossprod(Lk, Sigma_hat_true %*% pop$U_active[[k]]) / sqrt(CFG$K)
  }

  A <- symm(A)
  if (min(eigen(A, symmetric = TRUE, only.values = TRUE)$values) <= 1e-12) return(NA_real_)
  normalized_cross <- inverse_sqrt(A) %*% cross
  singular_values <- svd(normalized_cross, nu = 0, nv = 0)$d
  singular_values <- pmin(1, pmax(0, singular_values))
  sqrt(sum(1 - singular_values^2))
}

# ---------------------------- method fitting ----------------------------------

pair_index <- data.frame(
  key = c("1_2", "1_3", "2_3"),
  k = c(1L, 1L, 2L), ell = c(2L, 3L, 3L),
  stringsAsFactors = FALSE
)

fit_pair_blocks <- function(X, candidate_sets, rho, cores = 1L, warm_blocks = NULL) {
  Xs <- lapply(seq_len(CFG$K), function(k) X[[k]][, candidate_sets[[k]], drop = FALSE])
  stats <- lapply(seq_len(nrow(pair_index)), function(i) {
    make_pair_stats(Xs[[pair_index$k[i]]], Xs[[pair_index$ell[i]]])
  })

  fit_one <- function(i) {
    key <- pair_index$key[i]
    warm <- if (!is.null(warm_blocks) && !is.null(warm_blocks[[key]])) warm_blocks[[key]]$value else NULL
    fit <- pair_admm(stats[[i]], rho, warm_Z = warm)
    list(
      key = key, k = pair_index$k[i], ell = pair_index$ell[i],
      fit = fit,
      block = block_repr(
        candidate_sets[[pair_index$k[i]]],
        candidate_sets[[pair_index$ell[i]]],
        fit$C, pair_index$k[i], pair_index$ell[i]
      )
    )
  }

  ncores <- min(3L, max(1L, as.integer(cores)))
  if (.Platform$OS.type == "unix" && ncores > 1L) {
    fits <- parallel::mclapply(seq_len(nrow(pair_index)), fit_one,
                              mc.cores = ncores, mc.set.seed = FALSE)
  } else {
    fits <- lapply(seq_len(nrow(pair_index)), fit_one)
  }

  blocks <- setNames(lapply(fits, `[[`, "block"), vapply(fits, `[[`, character(1), "key"))
  list(
    blocks = blocks,
    iterations = sum(vapply(fits, function(z) z$fit$iterations, integer(1))),
    all_converged = all(vapply(fits, function(z) z$fit$converged, logical(1))),
    cumulative_pair_seconds = sum(vapply(fits, function(z) z$fit$seconds, numeric(1)))
  )
}

fit_egcar_final <- function(X_centered, pop, rho, rank, cores = 1L) {
  start <- proc.time()[3]
  screen <- screen_rows(X_centered, pop$supports)
  pair_fit <- fit_pair_blocks(X_centered, screen$sets, rho, cores = cores)
  loading <- fit_global_loading(X_centered, pair_fit$blocks, rank)
  list(
    blocks = pair_fit$blocks, loading = loading, rho = rho,
    screen = screen, iterations = pair_fit$iterations,
    converged = pair_fit$all_converged,
    cumulative_pair_seconds = pair_fit$cumulative_pair_seconds,
    final_fit_seconds = unname(proc.time()[3] - start)
  )
}

fit_oracle <- function(X_centered, pop, rank, cores = 1L) {
  start <- proc.time()[3]

  fit_one <- function(i) {
    k <- pair_index$k[i]
    ell <- pair_index$ell[i]
    Xk <- X_centered[[k]][, pop$supports[[k]], drop = FALSE]
    Xell <- X_centered[[ell]][, pop$supports[[ell]], drop = FALSE]
    fit <- oracle_pair(Xk, Xell)
    list(
      key = pair_index$key[i], fit = fit,
      block = block_repr(pop$supports[[k]], pop$supports[[ell]], fit$C, k, ell)
    )
  }

  ncores <- min(3L, max(1L, as.integer(cores)))
  if (.Platform$OS.type == "unix" && ncores > 1L) {
    fits <- parallel::mclapply(seq_len(nrow(pair_index)), fit_one,
                              mc.cores = ncores, mc.set.seed = FALSE)
  } else {
    fits <- lapply(seq_len(nrow(pair_index)), fit_one)
  }

  blocks <- setNames(lapply(fits, `[[`, "block"), vapply(fits, `[[`, character(1), "key"))
  loading <- fit_global_loading(X_centered, blocks, rank, forced_sets = pop$supports)
  list(
    blocks = blocks, loading = loading, rho = NA_real_,
    screen = list(recall = 1, sets = pop$supports, seconds = 0),
    iterations = 0L, converged = TRUE,
    cumulative_pair_seconds = sum(vapply(fits, function(z) z$fit$seconds, numeric(1))),
    final_fit_seconds = unname(proc.time()[3] - start)
  )
}

# ------------------------- cross-validation -----------------------------------

make_folds <- function(n, folds, seed) {
  set.seed(seed)
  perm <- sample.int(n)
  ids <- rep(seq_len(folds), length.out = n)
  out <- integer(n)
  out[perm] <- ids
  out
}

fit_one_cv_fold <- function(fold_id, fold_ids, X_raw, pop, rank) {
  validation_idx <- which(fold_ids == fold_id)
  train_idx <- which(fold_ids != fold_id)
  split <- center_train_validation(X_raw, train_idx, validation_idx)
  screen <- screen_rows(split$train, pop$supports)

  # Pair statistics are computed once per fold and reused for all penalties.
  Xs <- lapply(seq_len(CFG$K), function(k) {
    split$train[[k]][, screen$sets[[k]], drop = FALSE]
  })
  stats <- lapply(seq_len(nrow(pair_index)), function(i) {
    make_pair_stats(Xs[[pair_index$k[i]]], Xs[[pair_index$ell[i]]])
  })

  states <- setNames(vector("list", nrow(pair_index)), pair_index$key)
  rows <- vector("list", length(RHO_GRID))
  rho_order <- sort(RHO_GRID, decreasing = TRUE)

  for (j in seq_along(rho_order)) {
    rho <- rho_order[j]
    blocks <- list()
    total_iterations <- 0L
    all_converged <- TRUE
    cumulative_pair_seconds <- 0

    for (i in seq_len(nrow(pair_index))) {
      key <- pair_index$key[i]
      warm <- states[[key]]
      fit <- pair_admm(stats[[i]], rho, warm_Z = warm)
      states[[key]] <- fit$C
      blocks[[key]] <- block_repr(
        screen$sets[[pair_index$k[i]]],
        screen$sets[[pair_index$ell[i]]],
        fit$C, pair_index$k[i], pair_index$ell[i]
      )
      total_iterations <- total_iterations + fit$iterations
      all_converged <- all_converged && fit$converged
      cumulative_pair_seconds <- cumulative_pair_seconds + fit$seconds
    }

    loading <- fit_global_loading(split$train, blocks, rank)
    score <- heldout_rayleigh_score(loading, split$validation)
    rows[[j]] <- data.frame(
      fold = fold_id, rho = rho, validation_score = score,
      selected_rows = loading$selected_total,
      screening_recall = screen$recall,
      screening_seconds = screen$seconds,
      cumulative_pair_seconds = cumulative_pair_seconds,
      admm_iterations = total_iterations,
      all_converged = all_converged,
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, rows)
}

select_rho_cv <- function(X_raw, pop, rank, cores = 1L) {
  start <- proc.time()[3]
  n <- nrow(X_raw[[1]])
  folds <- min(CFG$cv_folds, n)
  fold_seed <- pop$population_seed + 7000000L + n
  fold_ids <- make_folds(n, folds, fold_seed)

  fold_worker <- function(f) fit_one_cv_fold(f, fold_ids, X_raw, pop, rank)
  ncores <- min(folds, max(1L, as.integer(cores)))
  if (.Platform$OS.type == "unix" && ncores > 1L) {
    pieces <- parallel::mclapply(seq_len(folds), fold_worker,
                                mc.cores = ncores, mc.set.seed = FALSE)
  } else {
    pieces <- lapply(seq_len(folds), fold_worker)
  }
  detail <- do.call(rbind, pieces)

  summary_rows <- lapply(RHO_GRID, function(rho) {
    z <- detail[detail$rho == rho, , drop = FALSE]
    finite <- is.finite(z$validation_score)
    data.frame(
      rho = rho,
      mean_validation_score = if (any(finite)) mean(z$validation_score[finite]) else NA_real_,
      se_validation_score = if (sum(finite) > 1L) sd(z$validation_score[finite]) / sqrt(sum(finite)) else 0,
      valid_folds = sum(finite),
      mean_selected_rows = if (any(finite)) mean(z$selected_rows[finite]) else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  summary <- do.call(rbind, summary_rows)

  eligible <- summary[summary$valid_folds == folds & is.finite(summary$mean_validation_score), , drop = FALSE]
  if (!nrow(eligible)) {
    eligible <- summary[summary$valid_folds >= ceiling(folds / 2) & is.finite(summary$mean_validation_score), , drop = FALSE]
  }
  if (!nrow(eligible)) {
    warning("No CV penalty produced a valid Rayleigh score; using the smallest rho")
    selected_rho <- min(RHO_GRID)
  } else {
    best <- max(eligible$mean_validation_score)
    # Exact-score ties are resolved toward the larger, sparser penalty.
    selected_rho <- max(eligible$rho[abs(eligible$mean_validation_score - best) <= 1e-12])
  }

  list(
    rho = selected_rho, detail = detail, summary = summary,
    folds = fold_ids, tuning_seconds = unname(proc.time()[3] - start)
  )
}

# ------------------------------- metrics --------------------------------------

compare_blocks_to_truth <- function(blocks, C0) {
  diff2 <- true2 <- est2 <- inner <- 0
  for (key in names(C0)) {
    T <- C0[[key]]
    E <- blocks[[key]]
    rows <- sort(unique(c(T$rows, E$rows)))
    cols <- sort(unique(c(T$cols, E$cols)))
    Tmat <- extract_block(T, rows, cols)
    Emat <- extract_block(E, rows, cols)
    D <- Emat - Tmat
    diff2 <- diff2 + sum(D^2)
    true2 <- true2 + sum(Tmat^2)
    est2 <- est2 + sum(Emat^2)
    inner <- inner + sum(Emat * Tmat)
  }
  list(
    c_euclidean_distance = sqrt(2 * diff2),
    c_squared_loss = 2 * diff2,
    c_relative_error = sqrt(diff2 / true2),
    c_cosine_similarity = if (est2 > 0 && true2 > 0) inner / sqrt(est2 * true2) else 0,
    c0_frobenius = sqrt(2 * true2)
  )
}

support_metrics <- function(fit, pop, oracle = FALSE) {
  selected <- if (oracle) pop$supports else fit$loading$sets
  tp <- fp <- fn <- 0L
  for (k in seq_len(CFG$K)) {
    tp <- tp + sum(selected[[k]] %in% pop$supports[[k]])
    fp <- fp + sum(!(selected[[k]] %in% pop$supports[[k]]))
    fn <- fn + sum(!(pop$supports[[k]] %in% selected[[k]]))
  }
  precision <- if (tp + fp > 0) tp / (tp + fp) else 0
  recall <- if (tp + fn > 0) tp / (tp + fn) else 1
  f1 <- if (precision + recall > 0) 2 * precision * recall / (precision + recall) else 0
  c(row_precision = precision, row_recall = recall, row_f1 = f1)
}

method_metric_row <- function(method, fit, pop, method_wall_seconds,
                              tuning_seconds = 0, oracle = FALSE) {
  cm <- compare_blocks_to_truth(fit$blocks, pop$C0)
  sm <- support_metrics(fit, pop, oracle = oracle)
  data.frame(
    method = method,
    rho = fit$rho,
    tau = TAU,
    c_euclidean_distance = cm$c_euclidean_distance,
    c_squared_loss = cm$c_squared_loss,
    c_relative_error = cm$c_relative_error,
    c_cosine_similarity = cm$c_cosine_similarity,
    c0_frobenius = cm$c0_frobenius,
    gca_subspace_distance = population_subspace_distance(fit$loading, pop),
    row_precision = unname(sm["row_precision"]),
    row_recall = unname(sm["row_recall"]),
    row_f1 = unname(sm["row_f1"]),
    selected_rows = fit$loading$selected_total,
    screening_recall = fit$screen$recall,
    admm_iterations = fit$iterations,
    admm_converged = fit$converged,
    cumulative_pair_seconds = fit$cumulative_pair_seconds,
    final_fit_seconds = fit$final_fit_seconds,
    tuning_seconds = tuning_seconds,
    method_wall_seconds = method_wall_seconds,
    stringsAsFactors = FALSE
  )
}

# -------------------------------- worker --------------------------------------

run_worker <- function(config_id, rep_id, outdir) {
  if (config_id < 1L || config_id > nrow(CONFIGS)) stop("config_id is out of range")
  if (rep_id < 1L || rep_id > N_REPS) stop("rep_id is out of range")
  make_dirs(outdir)

  cfg <- CONFIGS[config_id, , drop = FALSE]
  n <- cfg$n
  p <- cfg$p_per_block
  rank <- cfg$rank
  p_vec <- rep(as.integer(p), CFG$K)
  p_total <- sum(p_vec)
  cores <- available_cores()

  cat(sprintf(
    "Experiment=%s config=%d/%d rep=%d n=%d p/block=%d rank=%d cores=%d\n",
    EXPERIMENT, config_id, nrow(CONFIGS), rep_id, n, p, rank, cores
  ))

  pop <- make_population(p_vec, rank, rep_id)
  dat <- sample_views(pop, n, config_id, rep_id)
  X_centered <- center_full(dat$X)
  theory_rate <- sqrt(log(p_total) / n)
  theory_rho <- CFG$theory_rho_constant * theory_rate

  theory_start <- proc.time()[3]
  theory_fit <- fit_egcar_final(X_centered, pop, theory_rho, rank, cores)
  theory_wall <- unname(proc.time()[3] - theory_start)

  cv_start <- proc.time()[3]
  cv <- select_rho_cv(dat$X, pop, rank, cores)
  cv_fit <- fit_egcar_final(X_centered, pop, cv$rho, rank, cores)
  cv_wall <- unname(proc.time()[3] - cv_start)

  oracle_start <- proc.time()[3]
  oracle_fit <- fit_oracle(X_centered, pop, rank, cores)
  oracle_wall <- unname(proc.time()[3] - oracle_start)

  common <- data.frame(
    experiment = EXPERIMENT,
    config_id = config_id,
    panel = cfg$panel,
    n = n,
    p_per_block = p,
    p_total = p_total,
    rank = rank,
    rep_id = rep_id,
    K = CFG$K,
    support_per_block = CFG$support_size,
    signal_lambda = CFG$signal_lambda,
    toeplitz_alpha = CFG$alpha,
    theory_rate = theory_rate,
    theory_rho_constant = CFG$theory_rho_constant,
    cv_folds = CFG$cv_folds,
    cv_grid = paste(format(RHO_GRID, scientific = TRUE), collapse = ";"),
    population_seed = pop$population_seed,
    data_seed = dat$data_seed,
    stringsAsFactors = FALSE
  )

  metrics <- rbind(
    cbind(common, method_metric_row(
      "Oracle", oracle_fit, pop, oracle_wall, tuning_seconds = 0, oracle = TRUE
    )),
    cbind(common, method_metric_row(
      "EGCAR-theory", theory_fit, pop, theory_wall, tuning_seconds = 0, oracle = FALSE
    )),
    cbind(common, method_metric_row(
      "EGCAR-CV", cv_fit, pop, cv_wall, tuning_seconds = cv$tuning_seconds, oracle = FALSE
    ))
  )

  tag <- sprintf("c%03d_rep%02d", config_id, rep_id)
  cv_detail <- cbind(
    data.frame(config_id = config_id, rep_id = rep_id, n = n,
               p_per_block = p, rank = rank, stringsAsFactors = FALSE),
    cv$detail
  )
  cv_summary <- cbind(
    data.frame(config_id = config_id, rep_id = rep_id, n = n,
               p_per_block = p, rank = rank, selected_rho = cv$rho,
               stringsAsFactors = FALSE),
    cv$summary
  )
  atomic_csv(cv_detail, file.path(outdir, "cv", paste0("cv_detail_", tag, ".csv")))
  atomic_csv(cv_summary, file.path(outdir, "cv", paste0("cv_summary_", tag, ".csv")))

  atomic_csv(data.frame(
    config_id = config_id, rep_id = rep_id, n = n, p_per_block = p, rank = rank,
    theory_rho = theory_rho, selected_cv_rho = cv$rho, tau = TAU,
    theory_screen_sizes = paste(vapply(theory_fit$screen$sets, length, integer(1)), collapse = ";"),
    cv_screen_sizes = paste(vapply(cv_fit$screen$sets, length, integer(1)), collapse = ";"),
    stringsAsFactors = FALSE
  ), file.path(outdir, "config", paste0("task_config_", tag, ".csv")))

  if (rep_id == 1L) {
    atomic_rds(list(
      experiment = EXPERIMENT, config = cfg, rep_id = rep_id,
      p_vec = p_vec, rank = rank, C0 = pop$C0, supports = pop$supports,
      theory = theory_fit, cv = cv_fit, oracle = oracle_fit,
      selected_cv_rho = cv$rho, theory_rho = theory_rho
    ), file.path(outdir, "snapshots", paste0("snapshot_", tag, ".rds")))
  }

  # Completion marker: written last so the shell script can count finished tasks.
  atomic_csv(metrics, file.path(outdir, "metrics", paste0("metrics_", tag, ".csv")))
  cat(sprintf(
    "Done: Oracle distance=%.4g, theory=%.4g, CV=%.4g, selected rho=%.3g\n",
    metrics$c_euclidean_distance[metrics$method == "Oracle"],
    metrics$c_euclidean_distance[metrics$method == "EGCAR-theory"],
    metrics$c_euclidean_distance[metrics$method == "EGCAR-CV"],
    cv$rho
  ))
}

# --------------------------- aggregation and plotting -------------------------

summarize_metrics <- function(M) {
  grouping <- c("experiment", "panel", "n", "p_per_block", "p_total", "rank", "method")
  metrics <- c(
    "c_euclidean_distance", "c_squared_loss", "c_relative_error",
    "c_cosine_similarity", "gca_subspace_distance",
    "row_precision", "row_recall", "row_f1", "selected_rows",
    "screening_recall", "final_fit_seconds", "tuning_seconds",
    "method_wall_seconds", "theory_rate", "rho"
  )
  key <- do.call(interaction, c(M[grouping], list(drop = TRUE, lex.order = TRUE)))
  pieces <- split(M, key)
  ans <- lapply(pieces, function(d) {
    out <- d[1, grouping, drop = FALSE]
    out$repetitions <- length(unique(d$rep_id))
    for (v in metrics) {
      x <- d[[v]]
      x <- x[is.finite(x)]
      out[[paste0(v, "_mean")]] <- if (length(x)) mean(x) else NA_real_
      out[[paste0(v, "_sd")]] <- if (length(x) > 1L) sd(x) else 0
      out[[paste0(v, "_se")]] <- if (length(x) > 1L) sd(x) / sqrt(length(x)) else 0
    }
    out
  })
  S <- do.call(rbind, ans)
  rownames(S) <- NULL
  S[order(S$panel, S$rank, S$n, S$p_per_block, S$method), , drop = FALSE]
}

method_style <- function(methods) {
  all_methods <- c("Oracle", "EGCAR-theory", "EGCAR-CV")
  idx <- match(methods, all_methods)
  list(col = c(1, 2, 4)[idx], pch = c(16, 17, 15)[idx], lty = c(1, 2, 3)[idx])
}

plot_mean_se <- function(S, xvar, metric, ylab, path,
                         logx = TRUE, logy = FALSE, add_theory_rate = FALSE,
                         main = "") {
  methods <- c("Oracle", "EGCAR-theory", "EGCAR-CV")
  methods <- methods[methods %in% unique(S$method)]
  style <- method_style(methods)

  yname <- paste0(metric, "_mean")
  sename <- paste0(metric, "_se")
  yy <- S[[yname]]
  ee <- S[[sename]]
  extra <- numeric(0)
  if (add_theory_rate) extra <- unique(S$theory_rate_mean)
  ylim <- range(c(yy - ee, yy + ee, extra), finite = TRUE)
  if (!all(is.finite(ylim)) || diff(ylim) == 0) ylim <- c(0, max(1, yy, na.rm = TRUE))
  if (logy) {
    positive <- c(yy - ee, yy + ee, extra)
    positive <- positive[is.finite(positive) & positive > 0]
    if (!length(positive)) positive <- c(1e-8, 1)
    ylim <- range(positive)
    ylim <- c(max(min(positive) * 0.8, 1e-12), max(positive) * 1.25)
  } else {
    pad <- 0.06 * diff(ylim)
    if (!is.finite(pad) || pad == 0) pad <- 0.1
    ylim <- ylim + c(-pad, pad)
  }

  pdf(path, width = 7.4, height = 5.4)
  on.exit(dev.off(), add = TRUE)
  plot(range(S[[xvar]], finite = TRUE), ylim, type = "n",
       log = paste0(if (logx) "x" else "", if (logy) "y" else ""),
       xlab = X_LABEL, ylab = ylab, main = main)
  grid()

  for (i in seq_along(methods)) {
    d <- S[S$method == methods[i], , drop = FALSE]
    d <- d[order(d[[xvar]]), , drop = FALSE]
    x <- d[[xvar]]
    y <- d[[yname]]
    se <- d[[sename]]
    lines(x, y, type = "b", lwd = 1.7, col = style$col[i],
          pch = style$pch[i], lty = style$lty[i])
    lower <- y - se
    upper <- y + se
    if (logy) lower <- pmax(lower, .Machine$double.eps)
    arrows(x, lower, x, upper, angle = 90, code = 3, length = 0.035,
           col = style$col[i])
  }

  legend_labels <- methods
  legend_col <- style$col
  legend_pch <- style$pch
  legend_lty <- style$lty
  if (add_theory_rate) {
    rate_data <- unique(S[, c(xvar, "theory_rate_mean"), drop = FALSE])
    rate_data <- rate_data[order(rate_data[[xvar]]), , drop = FALSE]
    lines(rate_data[[xvar]], rate_data$theory_rate_mean,
          lwd = 2, lty = 4, col = 6)
    legend_labels <- c(legend_labels, "sqrt(log(p_total)/n) reference")
    legend_col <- c(legend_col, 6)
    legend_pch <- c(legend_pch, NA)
    legend_lty <- c(legend_lty, 4)
  }
  legend("topright", legend_labels, col = legend_col, pch = legend_pch,
         lty = legend_lty, lwd = 1.7, bty = "n", cex = 0.9)
}

assemble_display_matrix <- function(blocks, sets) {
  sizes <- vapply(sets, length, integer(1))
  offsets <- c(0L, cumsum(sizes))
  A <- matrix(0, sum(sizes), sum(sizes))
  for (B in blocks) {
    k <- B$k
    ell <- B$ell
    if (!sizes[k] || !sizes[ell]) next
    rr <- (offsets[k] + 1L):offsets[k + 1L]
    cc <- (offsets[ell] + 1L):offsets[ell + 1L]
    M <- extract_block(B, sets[[k]], sets[[ell]])
    A[rr, cc] <- M
    A[cc, rr] <- t(M)
  }
  list(A = A, boundaries = cumsum(sizes))
}

display_sets_for_snapshot <- function(snap) {
  p_vec <- snap$p_vec
  if (max(p_vec) <= 150L) return(lapply(p_vec, seq_len))

  out <- vector("list", CFG$K)
  for (k in seq_len(CFG$K)) {
    truth <- snap$supports[[k]]
    nt <- snap$theory$loading$row_norms[[k]]
    nc <- snap$cv$loading$row_norms[[k]]
    top_t <- head(order(nt, decreasing = TRUE), CFG$heatmap_max_per_block)
    top_c <- head(order(nc, decreasing = TRUE), CFG$heatmap_max_per_block)
    candidate <- unique(c(truth, top_t, top_c))
    score <- pmax(nt[candidate], nc[candidate])
    ord <- order(score, decreasing = TRUE)
    keep <- head(candidate[ord], CFG$heatmap_max_per_block)
    out[[k]] <- sort(unique(c(truth, keep)))
  }
  out
}

matrix_heatmap <- function(A, title, zlim, boundaries = NULL) {
  palette <- colorRampPalette(c("navy", "white", "firebrick3"))(101)
  image(seq_len(nrow(A)), seq_len(ncol(A)), A,
        col = palette, zlim = zlim, axes = FALSE,
        xlab = "", ylab = "", main = title, useRaster = TRUE, asp = 1)
  box()
  if (!is.null(boundaries) && length(boundaries) > 1L) {
    cuts <- boundaries[-length(boundaries)] + 0.5
    abline(v = cuts, h = cuts, lwd = 1.1)
  }
}

make_heatmap_pdf <- function(outdir) {
  files <- list.files(file.path(outdir, "snapshots"), "^snapshot_.*\\.rds$", full.names = TRUE)
  if (!length(files)) return(invisible(NULL))
  files <- sort(files)
  pdf(file.path(outdir, "plots", "C0_and_reconstructed_C_heatmaps_rep1.pdf"),
      width = 13.5, height = 4.5, onefile = TRUE)
  for (f in files) {
    snap <- readRDS(f)
    sets <- display_sets_for_snapshot(snap)
    matrices <- list(
      "Ground-truth C0" = assemble_display_matrix(snap$C0, sets),
      "Oracle" = assemble_display_matrix(snap$oracle$blocks, sets),
      "EGCAR-theory" = assemble_display_matrix(snap$theory$blocks, sets),
      "EGCAR-CV" = assemble_display_matrix(snap$cv$blocks, sets)
    )
    zmax <- max(abs(unlist(lapply(matrices, function(z) z$A))))
    if (!is.finite(zmax) || zmax == 0) zmax <- 1
    par(mfrow = c(1, 4), mar = c(1, 1, 3, 1), oma = c(0, 0, 2, 0))
    for (nm in names(matrices)) {
      matrix_heatmap(matrices[[nm]]$A, nm, c(-zmax, zmax), matrices[[nm]]$boundaries)
    }
    mtext(sprintf("n=%d, p/block=%d, rank=%d, CV rho=%.3g",
                  snap$config$n, snap$config$p_per_block, snap$rank,
                  snap$selected_cv_rho), outer = TRUE)
  }
  dev.off()
}

aggregate_results <- function(outdir) {
  make_dirs(outdir)
  files <- list.files(file.path(outdir, "metrics"), "^metrics_c[0-9]+_rep[0-9]+\\.csv$", full.names = TRUE)
  if (!length(files)) stop("No worker metric files were found")
  M <- do.call(rbind, lapply(files, read.csv, check.names = FALSE))
  rownames(M) <- NULL
  M <- M[order(M$config_id, M$rep_id, M$method), , drop = FALSE]
  atomic_csv(M, file.path(outdir, "all_metrics.csv"))

  S <- summarize_metrics(M)
  atomic_csv(S, file.path(outdir, "summary_metrics_mean_over_10_repetitions.csv"))

  cv_files <- list.files(file.path(outdir, "cv"), "^cv_summary_.*\\.csv$", full.names = TRUE)
  if (length(cv_files)) {
    CV <- do.call(rbind, lapply(cv_files, read.csv, check.names = FALSE))
    atomic_csv(CV, file.path(outdir, "all_cv_penalty_summaries.csv"))
  }

  for (rank in RANKS) {
    d <- S[S$rank == rank, , drop = FALSE]
    prefix <- file.path(outdir, "plots", sprintf("n150_rank%d", rank))
    title_suffix <- sprintf("n=150, rank=%d; means over 10 repetitions", rank)
    plot_mean_se(d, "p_per_block", "c_euclidean_distance",
                 "Euclidean/Frobenius distance ||C-hat - C0||_F",
                 paste0(prefix, "_C_euclidean_distance.pdf"),
                 logx = TRUE, main = title_suffix)
    plot_mean_se(d, "p_per_block", "c_squared_loss",
                 "Squared Frobenius loss ||C-hat - C0||_F^2",
                 paste0(prefix, "_C_squared_loss.pdf"),
                 logx = TRUE, logy = TRUE, main = title_suffix)
    plot_mean_se(d, "p_per_block", "c_relative_error",
                 "Relative C-matrix error",
                 paste0(prefix, "_relative_C_error_with_theory_rate.pdf"),
                 logx = TRUE, logy = TRUE, add_theory_rate = TRUE,
                 main = title_suffix)
    plot_mean_se(d, "p_per_block", "gca_subspace_distance",
                 "GCA loading-subspace distance ||sin Theta||_F",
                 paste0(prefix, "_GCA_subspace_distance.pdf"),
                 logx = TRUE, main = title_suffix)
    plot_mean_se(d, "p_per_block", "method_wall_seconds",
                 "Wall-clock seconds, including CV tuning",
                 paste0(prefix, "_total_wall_time.pdf"),
                 logx = TRUE, logy = TRUE, main = title_suffix)
    plot_mean_se(d, "p_per_block", "final_fit_seconds",
                 "Final-fit wall-clock seconds, excluding CV tuning",
                 paste0(prefix, "_final_fit_time.pdf"),
                 logx = TRUE, logy = TRUE, main = title_suffix)
  }

  make_heatmap_pdf(outdir)
  cat("Aggregation complete. Results are in: ", outdir, "\n", sep = "")
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  usage()
  quit(status = 2L)
}

if (tolower(args[1]) == "aggregate") {
  aggregate_results(args[3])
} else {
  config_id <- suppressWarnings(as.integer(args[1]))
  rep_id <- suppressWarnings(as.integer(args[2]))
  if (!is.finite(config_id) || !is.finite(rep_id)) {
    usage()
    quit(status = 2L)
  }
  run_worker(config_id, rep_id, args[3])
}
