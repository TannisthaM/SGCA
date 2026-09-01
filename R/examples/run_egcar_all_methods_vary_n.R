#!/usr/bin/env Rscript

# =============================================================================
# EGCAR and multiview benchmark simulation: increasing sample size
#
# Methods:
#   Oracle1          Population Sigma/Sigma0, exact pairwise C construction,
#                    followed by the population spectral loading step.
#   Oracle2          True row supports, but pairwise C values and loadings are
#                    estimated from the sample.
#   EGCAR-theory     Pairwise l1-ADMM with rho=sqrt(log(p_total)/n).
#   EGCAR-CV         A single common rho for all three pairs, selected by
#                    five-fold held-out generalized-Rayleigh CV over 1e-5,...,1e4.
#   SGCA             Gao-Ma implementation used by the uploaded SGCAR script.
#   RGCCA / SGCCA    Unsupervised CV wrappers used by the uploaded SGCAR script;
#                    a direct RGCCA-package fallback supplies r components if an
#                    installed wrapper returns only one component.
#
# Population model: Section 6.1 of the EGCAR manuscript.
#   Sigma_kk[i,j] = alpha^|i-j|
#   U_k = G_k (G_k' Sigma_kk G_k)^(-1/2)
#   Sigma_kl = signal_lambda Sigma_kk U_k U_l' Sigma_ll
#   C0_kl = signal_lambda U_k U_l'
#
# Important: C0 has ZERO diagonal blocks.  This is C*=Sigma0^{-1} Omega
# Sigma0^{-1}, not Sigma0^{-1} Sigma Sigma0^{-1}.
# =============================================================================

options(stringsAsFactors = FALSE, warn = 1)
Sys.setenv(
  OMP_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1"
)

ARGS <- commandArgs(trailingOnly = TRUE)
MODE <- if (length(ARGS)) ARGS[[1]] else "help"

# ------------------------------- libraries -----------------------------------

R_LIB_USER <- path.expand(Sys.getenv("R_LIBS_USER", "~/Rlibs"))
dir.create(R_LIB_USER, recursive = TRUE, showWarnings = FALSE)
.libPaths(unique(c(R_LIB_USER, .libPaths())))

CRAN_PACKAGES <- c(
  "MASS", "Matrix", "future", "future.apply", "PMA", "pracma",
  "RGCCA", "RSpectra", "ggplot2", "remotes"
)

SGCAR_INTERNALS <- c(
  "gao_gca_cv_init_and_final",
  "rgcca_unsupervised_cv_tau"
)

install_required_packages <- function() {
  missing <- CRAN_PACKAGES[!vapply(
    CRAN_PACKAGES, requireNamespace, logical(1), quietly = TRUE
  )]

  if (length(missing)) {
    message("Installing missing CRAN packages into ", R_LIB_USER, ": ",
            paste(missing, collapse = ", "))
    install.packages(
      missing,
      repos = "https://cloud.r-project.org",
      lib = R_LIB_USER,
      dependencies = c("Depends", "Imports", "LinkingTo"),
      Ncpus = {
        requested <- suppressWarnings(as.integer(Sys.getenv("EGCAR_INSTALL_CORES", "2")))
        if (!is.finite(requested) || requested < 1L) requested <- 2L
        min(requested, 8L)
      }
    )
  }

  if (!requireNamespace("SGCAR", quietly = TRUE)) {
    local_source <- Sys.getenv("SGCAR_LOCAL_SOURCE", "")
    github_repo <- Sys.getenv("SGCAR_GITHUB_REPO", "")

    if (nzchar(local_source) && (file.exists(local_source) || dir.exists(local_source))) {
      message("Installing SGCAR from local source: ", local_source)
      if (dir.exists(local_source)) {
        remotes::install_local(
          local_source, lib = R_LIB_USER,
          dependencies = NA, upgrade = "never"
        )
      } else {
        install.packages(
          local_source, repos = NULL, type = "source", lib = R_LIB_USER
        )
      }
    } else if (nzchar(github_repo)) {
      message("Installing SGCAR from GitHub repository: ", github_repo)
      remotes::install_github(
        github_repo, lib = R_LIB_USER,
        dependencies = NA, upgrade = "never"
      )
    } else {
      stop(
        "The custom package SGCAR is not installed.  The uploaded benchmark ",
        "script also requires it.  Set SGCAR_LOCAL_SOURCE to a package ",
        "directory/tar.gz, or SGCAR_GITHUB_REPO to owner/repository, then ",
        "rerun install_packages."
      )
    }
  }

  invisible(check_required_packages(stop_on_missing = TRUE))
}

check_required_packages <- function(stop_on_missing = TRUE) {
  packages <- c(CRAN_PACKAGES, "SGCAR")
  rows <- lapply(packages, function(pkg) {
    ok <- requireNamespace(pkg, quietly = TRUE)
    data.frame(
      package = pkg,
      installed = ok,
      version = if (ok) as.character(utils::packageVersion(pkg)) else NA_character_,
      library = if (ok) find.package(pkg) else NA_character_,
      stringsAsFactors = FALSE
    )
  })
  tab <- do.call(rbind, rows)
  print(tab, row.names = FALSE)

  missing <- tab$package[!tab$installed]
  missing_internal <- character(0)
  if (requireNamespace("SGCAR", quietly = TRUE)) {
    ns <- asNamespace("SGCAR")
    missing_internal <- SGCAR_INTERNALS[
      !vapply(SGCAR_INTERNALS, exists, logical(1), envir = ns, inherits = FALSE)
    ]
    if (length(missing_internal)) {
      message("Missing required SGCAR internal functions: ",
              paste(missing_internal, collapse = ", "))
    }
  }

  if (stop_on_missing && (length(missing) || length(missing_internal))) {
    stop(
      "Package preflight failed. Run this script with mode install_packages, ",
      "then rerun check_packages."
    )
  }
  invisible(tab)
}

if (identical(MODE, "install_packages")) {
  install_required_packages()
  quit(save = "no", status = 0)
}
if (identical(MODE, "check_packages")) {
  check_required_packages(stop_on_missing = TRUE)
  quit(save = "no", status = 0)
}

# Worker and aggregation modes need these packages.  The custom SGCAR package
# is checked explicitly because this program calls two functions used in the
# uploaded benchmark program through SGCAR:::.
check_required_packages(stop_on_missing = MODE %in% c("worker", "aggregate"))

suppressPackageStartupMessages({
  library(Matrix)
  library(ggplot2)
})

# ------------------------------- configuration -------------------------------

get_env_num <- function(name, default) {
  x <- Sys.getenv(name, "")
  if (!nzchar(x)) return(default)
  out <- suppressWarnings(as.numeric(x))
  if (!is.finite(out)) stop("Environment variable ", name, " is not numeric")
  out
}
get_env_int <- function(name, default) {
  as.integer(round(get_env_num(name, default)))
}

CFG <- list(
  K = 3L,
  support_size = get_env_int("EGCAR_SUPPORT_SIZE", 5L),
  alpha = get_env_num("EGCAR_TOEPLITZ_ALPHA", 0.30),
  signal_lambda = get_env_num("EGCAR_SIGNAL_LAMBDA", 0.80),
  row_tau = get_env_num("EGCAR_ROW_TAU", 1e-4),
  theory_rho_constant = get_env_num("EGCAR_THEORY_RHO_CONSTANT", 1.0),
  cv_folds = get_env_int("EGCAR_CV_FOLDS", 5L),
  admm_mu = get_env_num("EGCAR_ADMM_MU", 1.0),
  admm_max_iter = get_env_int("EGCAR_ADMM_MAX_ITER", 500L),
  admm_abs_tol = get_env_num("EGCAR_ADMM_ABS_TOL", 1e-4),
  admm_rel_tol = get_env_num("EGCAR_ADMM_REL_TOL", 1e-3),
  covariance_ridge = get_env_num("EGCAR_COVARIANCE_RIDGE", 1e-6),
  oracle_ridge = get_env_num("EGCAR_ORACLE_RIDGE", 1e-8),
  c_reconstruction_ridge = get_env_num("EGCAR_C_RECON_RIDGE", 1e-8),
  screen_full_max_p = get_env_int("EGCAR_SCREEN_FULL_MAX_P", 200L),
  screen_cap = get_env_int("EGCAR_SCREEN_CAP", 250L),
  screen_support_multiplier = get_env_int("EGCAR_SCREEN_SUPPORT_MULTIPLIER", 8L),
  screen_sqrt_multiplier = get_env_num("EGCAR_SCREEN_SQRT_MULTIPLIER", 2.5),
  sgca_k = get_env_int("EGCAR_SGCA_K", 20L),
  sgca_eta = get_env_num("EGCAR_SGCA_ETA", 0.001),
  sgca_rho_scale = get_env_num("EGCAR_SGCA_RHO_SCALE", 0.5),
  full_matrix_save_max = get_env_int("EGCAR_FULL_MATRIX_SAVE_MAX", 1000L),
  base_seed = get_env_int("EGCAR_BASE_SEED", 20260901L)
)

if (CFG$support_size < 5L) {
  warning("The requested ranks include r=5; support_size has been raised to 5.")
  CFG$support_size <- 5L
}
if (!(CFG$signal_lambda > 0 && CFG$signal_lambda < 1)) {
  stop("EGCAR_SIGNAL_LAMBDA must lie strictly between 0 and 1")
}
if (!(abs(CFG$alpha) < 1)) stop("EGCAR_TOEPLITZ_ALPHA must satisfy |alpha|<1")

RHO_GRID <- 10^seq(-5, 4)
SGCA_LAMBDA_GRID <- c(0, RHO_GRID, 1e5)
RGCCA_TUNING_GRID <- c(1e-6, 0.001, 0.1, 0.25, 0.5, 0.75, 1)
METHOD_LEVELS <- c(
  "Oracle1-population", "Oracle2-support", "EGCAR-theory", "EGCAR-CV",
  "SGCA", "RGCCA", "SGCCA"
)
N_REPS <- 10L
EXPERIMENT <- "vary_n"

build_configs <- function() {
  A <- expand.grid(
    rank = c(1L, 2L, 5L),
    n = c(200L, 250L, 300L, 500L, 1000L, 10000L),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  A$p_per_block <- 100L
  A$panel <- "p1=p2=p3=100"

  B <- expand.grid(
    rank = c(1L, 2L, 5L),
    n = c(250L, 300L, 450L, 1000L, 5000L, 10000L),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  B$p_per_block <- 150L
  B$panel <- "p1=p2=p3=150"

  out <- rbind(A, B)
  out <- out[order(out$p_per_block, out$rank, out$n),
             c("panel", "p_per_block", "n", "rank")]
  rownames(out) <- NULL
  out
}

CONFIGS <- build_configs()
CONFIGS$config_id <- seq_len(nrow(CONFIGS))
EXPECTED_TASKS <- nrow(CONFIGS) * N_REPS

available_cores <- function() {
  from_slurm <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "")))
  if (is.finite(from_slurm) && from_slurm > 0L) return(max(1L, from_slurm - 1L))
  detected <- parallel::detectCores(logical = FALSE)
  if (!is.finite(detected)) detected <- 1L
  max(1L, detected - 1L)
}

# -------------------------------- utilities ----------------------------------

fnorm <- function(A) sqrt(sum(A * A))
symm <- function(A) (A + t(A)) / 2
soft_threshold <- function(A, threshold) sign(A) * pmax(abs(A) - threshold, 0)

atomic_csv <- function(x, path, row.names = FALSE) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp.", Sys.getpid(), ".", sample.int(1e9, 1L))
  utils::write.csv(x, tmp, row.names = row.names)
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

make_output_dirs <- function(outdir) {
  for (d in c("metrics", "fits", "cv", "plots", "logs")) {
    dir.create(file.path(outdir, d), recursive = TRUE, showWarnings = FALSE)
  }
}

inverse_sqrt <- function(A, floor = 1e-12) {
  ee <- eigen(symm(A), symmetric = TRUE)
  scale <- max(1, max(abs(ee$values)))
  if (min(ee$values) <= floor * scale) {
    stop("Matrix is rank deficient in inverse_sqrt")
  }
  ee$vectors %*%
    (diag(1 / sqrt(ee$values), nrow = length(ee$values)) %*% t(ee$vectors))
}

matrix_roots <- function(A, ridge = 0) {
  d <- nrow(A)
  if (d == 0L) {
    return(list(sqrt = matrix(0, 0, 0), invsqrt = matrix(0, 0, 0)))
  }
  scale <- max(1, mean(diag(A)))
  Areg <- symm(A) + diag(ridge * scale, d)
  ee <- eigen(Areg, symmetric = TRUE)
  vals <- pmax(ee$values, .Machine$double.eps * max(1, max(abs(ee$values))))
  list(
    sqrt = ee$vectors %*% (diag(sqrt(vals), d) %*% t(ee$vectors)),
    invsqrt = ee$vectors %*% (diag(1 / sqrt(vals), d) %*% t(ee$vectors))
  )
}

toeplitz_submatrix <- function(rows, cols, alpha = CFG$alpha) {
  if (!length(rows) || !length(cols)) return(matrix(0, length(rows), length(cols)))
  alpha ^ abs(outer(rows, cols, "-"))
}

pair_index <- data.frame(
  key = c("1_2", "1_3", "2_3"),
  k = c(1L, 1L, 2L),
  ell = c(2L, 3L, 3L),
  stringsAsFactors = FALSE
)

# -------------------------- block representations ----------------------------

make_dense_block <- function(rows, cols, value, k, ell) {
  list(
    type = "dense", rows = as.integer(rows), cols = as.integer(cols),
    value = as.matrix(value), k = as.integer(k), ell = as.integer(ell)
  )
}

make_lowrank_block <- function(rows, cols, A, M, B, k, ell) {
  list(
    type = "lowrank", rows = as.integer(rows), cols = as.integer(cols),
    A = as.matrix(A), M = as.matrix(M), B = as.matrix(B),
    k = as.integer(k), ell = as.integer(ell)
  )
}

block_sqnorm <- function(B) {
  if (identical(B$type, "dense")) return(sum(B$value^2))
  GA <- crossprod(B$A)
  GB <- crossprod(B$B)
  val <- sum(B$M * (GA %*% B$M %*% GB))
  max(0, as.numeric(val))
}

block_inner_with_truth <- function(B, T) {
  stopifnot(identical(T$type, "dense"))
  if (identical(B$type, "dense")) {
    rows <- intersect(B$rows, T$rows)
    cols <- intersect(B$cols, T$cols)
    if (!length(rows) || !length(cols)) return(0)
    Eb <- B$value[match(rows, B$rows), match(cols, B$cols), drop = FALSE]
    Tb <- T$value[match(rows, T$rows), match(cols, T$cols), drop = FALSE]
    return(sum(Eb * Tb))
  }

  ir <- match(T$rows, B$rows)
  jc <- match(T$cols, B$cols)
  keep_r <- which(!is.na(ir))
  keep_c <- which(!is.na(jc))
  if (!length(keep_r) || !length(keep_c)) return(0)
  Eon <- B$A[ir[keep_r], , drop = FALSE] %*% B$M %*%
    t(B$B[jc[keep_c], , drop = FALSE])
  Ton <- T$value[keep_r, keep_c, drop = FALSE]
  sum(Eon * Ton)
}

block_diff2_from_truth <- function(B, T) {
  val <- block_sqnorm(B) + sum(T$value^2) - 2 * block_inner_with_truth(B, T)
  max(0, as.numeric(val))
}

compare_C_blocks <- function(blocks, C0) {
  diff2 <- true2 <- est2 <- inner <- 0
  for (key in names(C0)) {
    T <- C0[[key]]
    B <- blocks[[key]]
    diff2 <- diff2 + block_diff2_from_truth(B, T)
    true2 <- true2 + sum(T$value^2)
    est2 <- est2 + block_sqnorm(B)
    inner <- inner + block_inner_with_truth(B, T)
  }
  list(
    c_frobenius_distance = sqrt(2 * diff2),
    c_squared_loss = 2 * diff2,
    c_relative_error = if (true2 > 0) sqrt(diff2 / true2) else NA_real_,
    c_cosine_similarity = if (est2 > 0 && true2 > 0) inner / sqrt(est2 * true2) else 0,
    c0_frobenius = sqrt(2 * true2)
  )
}

materialize_block <- function(B, p_k, p_ell) {
  out <- matrix(0, p_k, p_ell)
  if (!length(B$rows) || !length(B$cols)) return(out)
  if (identical(B$type, "dense")) {
    out[B$rows, B$cols] <- B$value
  } else {
    out[B$rows, B$cols] <- B$A %*% B$M %*% t(B$B)
  }
  out
}

materialize_C <- function(blocks, p_vec) {
  offsets <- c(0L, cumsum(p_vec))
  p_total <- sum(p_vec)
  out <- matrix(0, p_total, p_total)
  for (key in names(blocks)) {
    B <- blocks[[key]]
    rr <- offsets[B$k] + seq_len(p_vec[B$k])
    cc <- offsets[B$ell] + seq_len(p_vec[B$ell])
    M <- materialize_block(B, p_vec[B$k], p_vec[B$ell])
    out[rr, cc] <- M
    out[cc, rr] <- t(M)
  }
  out
}

materialize_population_covariances <- function(pop) {
  p_vec <- pop$p_vec
  offsets <- c(0L, cumsum(p_vec))
  p_total <- sum(p_vec)
  Sigma0 <- matrix(0, p_total, p_total)
  for (k in seq_len(CFG$K)) {
    ii <- (offsets[k] + 1L):offsets[k + 1L]
    Sigma0[ii, ii] <- toeplitz_submatrix(seq_len(p_vec[k]), seq_len(p_vec[k]))
  }
  Sigma <- Sigma0
  for (i in seq_len(nrow(pair_index))) {
    k <- pair_index$k[i]
    ell <- pair_index$ell[i]
    ii <- (offsets[k] + 1L):offsets[k + 1L]
    jj <- (offsets[ell] + 1L):offsets[ell + 1L]
    Sigma_kl <- CFG$signal_lambda * pop$SigmaU[[k]] %*% t(pop$SigmaU[[ell]])
    Sigma[ii, jj] <- Sigma_kl
    Sigma[jj, ii] <- t(Sigma_kl)
  }
  list(Sigma = Sigma, Sigma0 = Sigma0)
}

# ------------------------ population and sampling -----------------------------

make_population <- function(p_vec, rank, rep_id) {
  if (rank > CFG$support_size) {
    stop("rank cannot exceed support_size under this simulation")
  }
  population_seed <- CFG$base_seed + 100000L * rank + rep_id
  set.seed(population_seed)

  support <- seq_len(CFG$support_size)
  Sigma_ss <- toeplitz_submatrix(support, support)
  supports <- U_active <- SigmaU <- vector("list", CFG$K)

  for (k in seq_len(CFG$K)) {
    ok <- FALSE
    for (attempt in seq_len(100L)) {
      G <- matrix(rnorm(CFG$support_size * rank), CFG$support_size, rank)
      gram <- crossprod(G, Sigma_ss %*% G)
      ev <- eigen(symm(gram), symmetric = TRUE, only.values = TRUE)$values
      if (min(ev) > 1e-10 * max(ev)) {
        ok <- TRUE
        break
      }
    }
    if (!ok) stop("Could not generate a well-conditioned active loading matrix")

    U <- G %*% inverse_sqrt(gram)
    check <- crossprod(U, Sigma_ss %*% U)
    if (max(abs(check - diag(rank))) > 1e-7) {
      stop("Population loading normalization failed")
    }

    supports[[k]] <- support
    U_active[[k]] <- U
    SigmaU[[k]] <- toeplitz_submatrix(seq_len(p_vec[k]), support) %*% U
  }

  C0 <- list()
  for (i in seq_len(nrow(pair_index))) {
    k <- pair_index$k[i]
    ell <- pair_index$ell[i]
    key <- pair_index$key[i]
    C0[[key]] <- make_dense_block(
      supports[[k]], supports[[ell]],
      CFG$signal_lambda * U_active[[k]] %*% t(U_active[[ell]]),
      k, ell
    )
  }

  L0 <- list(
    valid = TRUE,
    sets = supports,
    L = lapply(U_active, function(U) U / sqrt(CFG$K)),
    eigenvalues = rep((CFG$K - 1) * CFG$signal_lambda, rank),
    selected_total = CFG$K * CFG$support_size
  )

  list(
    p_vec = as.integer(p_vec), rank = rank, supports = supports,
    U_active = U_active, SigmaU = SigmaU, Sigma_ss = Sigma_ss,
    C0 = C0, L0 = L0, population_seed = population_seed
  )
}

ar1_gaussian <- function(n, p, alpha = CFG$alpha) {
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
  Z <- matrix(rnorm(n * pop$rank), n, pop$rank)
  adjustment <- sqrt(1 - CFG$signal_lambda) - 1
  X <- vector("list", CFG$K)

  for (k in seq_len(CFG$K)) {
    Y <- ar1_gaussian(n, pop$p_vec[k])
    YU <- Y[, pop$supports[[k]], drop = FALSE] %*% pop$U_active[[k]]
    X[[k]] <- Y +
      (adjustment * YU + sqrt(CFG$signal_lambda) * Z) %*%
      t(pop$SigmaU[[k]])
  }
  names(X) <- paste0("block", seq_len(CFG$K))
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
    validation[[k]] <- sweep(
      X[[k]][validation_idx, , drop = FALSE], 2L, means[[k]], "-"
    )
  }
  list(train = train, validation = validation)
}

# -------------------------- common screening ---------------------------------

screen_rows <- function(X, true_supports = NULL) {
  start <- proc.time()[3]
  n <- nrow(X[[1]])
  p_vec <- vapply(X, ncol, integer(1))

  budgets <- vapply(p_vec, function(p) {
    if (p <= CFG$screen_full_max_p) return(as.integer(p))
    rank_budget <- max(1L, n - 1L)
    root_budget <- as.integer(ceiling(CFG$screen_sqrt_multiplier * sqrt(p)))
    support_budget <- CFG$support_size * CFG$screen_support_multiplier
    as.integer(min(
      p,
      max(support_budget, min(CFG$screen_cap, max(rank_budget, root_budget)))
    ))
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

# ------------------------------- pairwise ADMM -------------------------------

covariance_eigensystem <- function(X) {
  n <- nrow(X)
  p <- ncol(X)
  if (p > n) {
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
  ek <- covariance_eigensystem(Xk)
  eell <- covariance_eigensystem(Xell)
  list(
    n = n, pk = ncol(Xk), pell = ncol(Xell),
    Skell = crossprod(Xk, Xell) / n,
    Qk = ek$vectors, Qell = eell$vectors,
    dk = ek$values, dell = eell$values,
    use_thin = ek$thin || eell$thin
  )
}

pair_admm <- function(stats, rho, warm_Z = NULL) {
  start <- proc.time()[3]
  pk <- stats$pk
  pell <- stats$pell

  if (max(abs(stats$Skell)) <= rho) {
    return(list(
      C = matrix(0, pk, pell), iterations = 0L, converged = TRUE,
      primal = 0, dual = 0, seconds = unname(proc.time()[3] - start)
    ))
  }

  covariance_products <- outer(stats$dk, stats$dell, "*")
  denom <- covariance_products + CFG$admm_mu
  if (is.null(warm_Z) || !all(dim(warm_Z) == c(pk, pell))) {
    Z <- matrix(0, pk, pell)
  } else {
    Z <- warm_Z
  }
  H <- matrix(0, pk, pell)
  C <- Z
  converged <- FALSE
  primal <- dual <- NA_real_

  for (iter in seq_len(CFG$admm_max_iter)) {
    G <- stats$Skell + CFG$admm_mu * (Z - H)
    if (isTRUE(stats$use_thin)) {
      if (!length(stats$dk) || !length(stats$dell)) {
        C <- G / CFG$admm_mu
      } else {
        A <- crossprod(stats$Qk, G %*% stats$Qell)
        correction <- (covariance_products * A) /
          (CFG$admm_mu * denom)
        C <- G / CFG$admm_mu -
          stats$Qk %*% correction %*% t(stats$Qell)
      }
    } else {
      transformed <- crossprod(stats$Qk, G %*% stats$Qell) / denom
      C <- stats$Qk %*% transformed %*% t(stats$Qell)
    }

    Z_old <- Z
    Z <- soft_threshold(C + H, rho / CFG$admm_mu)
    H <- H + C - Z

    primal <- fnorm(C - Z)
    dual <- CFG$admm_mu * fnorm(Z - Z_old)
    eps_primal <- sqrt(pk * pell) * CFG$admm_abs_tol +
      CFG$admm_rel_tol * max(fnorm(C), fnorm(Z))
    eps_dual <- sqrt(pk * pell) * CFG$admm_abs_tol +
      CFG$admm_rel_tol * CFG$admm_mu * fnorm(H)

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

fit_pair_blocks <- function(X, candidate_sets, rho, cores = 1L, warm_blocks = NULL) {
  Xs <- lapply(seq_len(CFG$K), function(k) {
    X[[k]][, candidate_sets[[k]], drop = FALSE]
  })
  stats <- lapply(seq_len(nrow(pair_index)), function(i) {
    make_pair_stats(Xs[[pair_index$k[i]]], Xs[[pair_index$ell[i]]])
  })

  fit_one <- function(i) {
    key <- pair_index$key[i]
    warm <- NULL
    if (!is.null(warm_blocks) && !is.null(warm_blocks[[key]]) &&
        identical(warm_blocks[[key]]$type, "dense")) {
      warm <- warm_blocks[[key]]$value
    }
    fit <- pair_admm(stats[[i]], rho, warm_Z = warm)
    list(
      key = key, fit = fit,
      block = make_dense_block(
        candidate_sets[[pair_index$k[i]]],
        candidate_sets[[pair_index$ell[i]]],
        fit$C, pair_index$k[i], pair_index$ell[i]
      )
    )
  }

  ncores <- min(3L, max(1L, as.integer(cores)))
  if (.Platform$OS.type == "unix" && ncores > 1L) {
    fits <- parallel::mclapply(
      seq_len(nrow(pair_index)), fit_one,
      mc.cores = ncores, mc.set.seed = FALSE
    )
  } else {
    fits <- lapply(seq_len(nrow(pair_index)), fit_one)
  }

  blocks <- setNames(
    lapply(fits, `[[`, "block"),
    vapply(fits, `[[`, character(1), "key")
  )
  list(
    blocks = blocks,
    iterations = sum(vapply(fits, function(z) z$fit$iterations, integer(1))),
    converged = all(vapply(fits, function(z) z$fit$converged, logical(1))),
    cumulative_pair_seconds = sum(vapply(fits, function(z) z$fit$seconds, numeric(1)))
  )
}

# ------------------------ row norms and GCA loadings --------------------------

row_norms_from_blocks <- function(blocks, p_vec) {
  accum <- lapply(p_vec, numeric)
  for (B in blocks) {
    if (identical(B$type, "dense")) {
      accum[[B$k]][B$rows] <- accum[[B$k]][B$rows] + rowSums(B$value^2)
      accum[[B$ell]][B$cols] <- accum[[B$ell]][B$cols] + colSums(B$value^2)
    } else {
      # Exact row/column norms of A M B' without materializing the block.
      GB <- crossprod(B$B)
      GA <- crossprod(B$A)
      AM <- B$A %*% B$M
      BMt <- B$B %*% t(B$M)
      accum[[B$k]][B$rows] <- accum[[B$k]][B$rows] +
        rowSums((AM %*% GB) * AM)
      accum[[B$ell]][B$cols] <- accum[[B$ell]][B$cols] +
        rowSums((BMt %*% GA) * BMt)
    }
  }
  lapply(accum, function(z) sqrt(pmax(z, 0)))
}

top_algebraic_eigen <- function(A, rank) {
  d <- nrow(A)
  if (d < rank) stop("Matrix dimension is smaller than requested rank")
  use_rspectra <- d > max(250L, 5L * rank) &&
    requireNamespace("RSpectra", quietly = TRUE) && rank < d
  if (use_rspectra) {
    ans <- tryCatch(
      RSpectra::eigs_sym(
        A, k = rank, which = "LA",
        opts = list(tol = 1e-8, maxitr = 10000L)
      ),
      error = function(e) NULL
    )
    if (!is.null(ans)) {
      ord <- order(ans$values, decreasing = TRUE)
      return(list(values = ans$values[ord], vectors = ans$vectors[, ord, drop = FALSE]))
    }
  }
  ee <- eigen(symm(A), symmetric = TRUE)
  list(
    values = ee$values[seq_len(rank)],
    vectors = ee$vectors[, seq_len(rank), drop = FALSE]
  )
}

extract_dense_subblock <- function(B, rows, cols) {
  if (!identical(B$type, "dense")) stop("Expected a dense block")
  out <- matrix(0, length(rows), length(cols))
  ir <- match(rows, B$rows)
  jc <- match(cols, B$cols)
  rr <- which(!is.na(ir))
  cc <- which(!is.na(jc))
  if (length(rr) && length(cc)) {
    out[rr, cc] <- B$value[ir[rr], jc[cc], drop = FALSE]
  }
  out
}

fit_global_loading <- function(X, blocks, rank, forced_sets = NULL) {
  p_vec <- vapply(X, ncol, integer(1))
  norms <- row_norms_from_blocks(blocks, p_vec)
  sets <- if (is.null(forced_sets)) {
    lapply(norms, function(z) which(z > CFG$row_tau))
  } else {
    lapply(forced_sets, as.integer)
  }

  sizes <- vapply(sets, length, integer(1))
  total_selected <- sum(sizes)
  if (total_selected < rank || total_selected == 0L) {
    return(list(
      valid = FALSE, sets = sets, L = NULL, row_norms = norms,
      eigenvalues = rep(NA_real_, rank), selected_total = total_selected
    ))
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
    Ckl <- if (identical(B$type, "dense")) {
      extract_dense_subblock(B, sets[[k]], sets[[ell]])
    } else {
      ir <- match(sets[[k]], B$rows)
      jc <- match(sets[[ell]], B$cols)
      out <- matrix(0, sizes[k], sizes[ell])
      rr <- which(!is.na(ir)); cc <- which(!is.na(jc))
      if (length(rr) && length(cc)) {
        out[rr, cc] <- B$A[ir[rr], , drop = FALSE] %*% B$M %*%
          t(B$B[jc[cc], , drop = FALSE])
      }
      out
    }
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

fit_population_loading <- function(pop, blocks, rank) {
  sets <- pop$supports
  sizes <- vapply(sets, length, integer(1))
  roots <- replicate(CFG$K, matrix_roots(pop$Sigma_ss, 0), simplify = FALSE)
  offsets <- c(0L, cumsum(sizes))
  Rmat <- matrix(0, sum(sizes), sum(sizes))

  for (key in names(blocks)) {
    B <- blocks[[key]]
    Ckl <- extract_dense_subblock(B, sets[[B$k]], sets[[B$ell]])
    Rkl <- roots[[B$k]]$sqrt %*% Ckl %*% roots[[B$ell]]$sqrt
    rr <- (offsets[B$k] + 1L):offsets[B$k + 1L]
    cc <- (offsets[B$ell] + 1L):offsets[B$ell + 1L]
    Rmat[rr, cc] <- Rkl
    Rmat[cc, rr] <- t(Rkl)
  }

  ee <- top_algebraic_eigen(symm(Rmat), rank)
  L <- vector("list", CFG$K)
  for (k in seq_len(CFG$K)) {
    rr <- (offsets[k] + 1L):offsets[k + 1L]
    L[[k]] <- roots[[k]]$invsqrt %*% ee$vectors[rr, , drop = FALSE]
  }
  list(
    valid = TRUE, sets = sets, L = L, eigenvalues = ee$values,
    selected_total = sum(sizes), row_norms = row_norms_from_blocks(blocks, pop$p_vec)
  )
}

heldout_rayleigh_score <- function(loading, X_validation) {
  if (is.null(loading) || !isTRUE(loading$valid)) return(NA_real_)
  rank <- ncol(loading$L[[which.max(vapply(loading$L, nrow, integer(1))) ]])
  n_val <- nrow(X_validation[[1]])
  Y <- vector("list", CFG$K)
  for (k in seq_len(CFG$K)) {
    if (!length(loading$sets[[k]])) {
      Y[[k]] <- matrix(0, n_val, rank)
    } else {
      Y[[k]] <- X_validation[[k]][, loading$sets[[k]], drop = FALSE] %*%
        loading$L[[k]]
    }
  }
  Ysum <- Reduce("+", Y)
  A <- crossprod(Ysum) / n_val
  Q <- Reduce("+", lapply(Y, function(Z) crossprod(Z) / n_val))
  roots <- tryCatch(matrix_roots(Q, CFG$covariance_ridge), error = function(e) NULL)
  if (is.null(roots)) return(NA_real_)
  score <- sum(diag(roots$invsqrt %*% A %*% roots$invsqrt))
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
    Sigma_hat_hat <- toeplitz_submatrix(S_hat, S_hat)
    Sigma_hat_true <- toeplitz_submatrix(S_hat, pop$supports[[k]])
    A <- A + crossprod(Lk, Sigma_hat_hat %*% Lk)
    cross <- cross +
      crossprod(Lk, Sigma_hat_true %*% pop$U_active[[k]]) / sqrt(CFG$K)
  }

  roots <- tryCatch(matrix_roots(A, 0), error = function(e) NULL)
  if (is.null(roots)) return(NA_real_)
  singular_values <- svd(roots$invsqrt %*% cross, nu = 0, nv = 0)$d
  singular_values <- pmin(1, pmax(0, singular_values))
  sqrt(sum(1 - singular_values^2))
}

euclidean_loading_procrustes <- function(loading, pop) {
  if (is.null(loading) || !isTRUE(loading$valid)) return(NA_real_)
  rank <- pop$rank
  GA <- matrix(0, rank, rank)
  GB <- matrix(0, rank, rank)
  cross <- matrix(0, rank, rank)

  for (k in seq_len(CFG$K)) {
    Lk <- loading$L[[k]]
    GA <- GA + crossprod(Lk)
    Bfull <- pop$U_active[[k]] / sqrt(CFG$K)
    GB <- GB + crossprod(Bfull)
    common <- intersect(loading$sets[[k]], pop$supports[[k]])
    if (length(common)) {
      Arows <- Lk[match(common, loading$sets[[k]]), , drop = FALSE]
      Brows <- Bfull[
        match(common, pop$supports[[k]]), , drop = FALSE
      ]
      cross <- cross + crossprod(Arows, Brows)
    }
  }

  ra <- tryCatch(matrix_roots(GA, 0), error = function(e) NULL)
  rb <- tryCatch(matrix_roots(GB, 0), error = function(e) NULL)
  if (is.null(ra) || is.null(rb)) return(NA_real_)
  s <- svd(ra$invsqrt %*% cross %*% rb$invsqrt, nu = 0, nv = 0)$d
  s <- pmin(1, pmax(0, s))
  sqrt(sum(1 - s^2))
}

# --------------------------- EGCAR estimators --------------------------------

fit_egcar_final <- function(X, pop, rho, rank, cores, fixed_screen = NULL) {
  start <- proc.time()[3]
  screen <- if (is.null(fixed_screen)) screen_rows(X, pop$supports) else fixed_screen
  pair_fit <- fit_pair_blocks(X, screen$sets, rho, cores = cores)
  c_seconds <- unname(proc.time()[3] - start)
  load_start <- proc.time()[3]
  loading <- fit_global_loading(X, pair_fit$blocks, rank)
  loading_step_seconds <- unname(proc.time()[3] - load_start)
  list(
    valid = isTRUE(loading$valid), blocks = pair_fit$blocks, loading = loading,
    rho = rho, selected_parameter = rho, screen = screen,
    iterations = pair_fit$iterations, converged = pair_fit$converged,
    cumulative_pair_seconds = pair_fit$cumulative_pair_seconds,
    c_seconds = c_seconds,
    loading_step_seconds = loading_step_seconds,
    canonical_seconds = c_seconds + loading_step_seconds,
    total_seconds = c_seconds + loading_step_seconds,
    backend = "pairwise-ADMM", error_message = ""
  )
}

oracle_pair_sample <- function(Xk, Xell) {
  n <- nrow(Xk)
  Sk <- crossprod(Xk) / n
  Sell <- crossprod(Xell) / n
  Skell <- crossprod(Xk, Xell) / n
  rk <- CFG$oracle_ridge * max(1, mean(diag(Sk)))
  rl <- CFG$oracle_ridge * max(1, mean(diag(Sell)))
  left <- solve(Sk + diag(rk, nrow(Sk)), Skell)
  t(solve(Sell + diag(rl, nrow(Sell)), t(left)))
}

fit_oracle2 <- function(X, pop, rank, cores) {
  start <- proc.time()[3]
  fit_one <- function(i) {
    k <- pair_index$k[i]; ell <- pair_index$ell[i]
    C <- oracle_pair_sample(
      X[[k]][, pop$supports[[k]], drop = FALSE],
      X[[ell]][, pop$supports[[ell]], drop = FALSE]
    )
    make_dense_block(pop$supports[[k]], pop$supports[[ell]], C, k, ell)
  }
  ncores <- min(3L, max(1L, cores))
  if (.Platform$OS.type == "unix" && ncores > 1L) {
    vals <- parallel::mclapply(
      seq_len(nrow(pair_index)), fit_one,
      mc.cores = ncores, mc.set.seed = FALSE
    )
  } else {
    vals <- lapply(seq_len(nrow(pair_index)), fit_one)
  }
  blocks <- setNames(vals, pair_index$key)
  c_seconds <- unname(proc.time()[3] - start)
  load_start <- proc.time()[3]
  loading <- fit_global_loading(X, blocks, rank, forced_sets = pop$supports)
  load_seconds <- unname(proc.time()[3] - load_start)
  list(
    valid = isTRUE(loading$valid), blocks = blocks, loading = loading,
    rho = NA_real_, selected_parameter = NA_real_,
    screen = list(sets = pop$supports, recall = 1, seconds = 0),
    iterations = 0L, converged = TRUE, cumulative_pair_seconds = c_seconds,
    c_seconds = c_seconds, loading_step_seconds = load_seconds,
    canonical_seconds = c_seconds + load_seconds,
    total_seconds = c_seconds + load_seconds,
    backend = "sample-oracle-support", error_message = ""
  )
}

fit_oracle1 <- function(pop, rank) {
  start <- proc.time()[3]
  blocks <- list()
  for (i in seq_len(nrow(pair_index))) {
    k <- pair_index$k[i]; ell <- pair_index$ell[i]; key <- pair_index$key[i]
    Sigma_kl_ss <- CFG$signal_lambda * pop$Sigma_ss %*%
      pop$U_active[[k]] %*% t(pop$U_active[[ell]]) %*% pop$Sigma_ss
    left <- solve(pop$Sigma_ss, Sigma_kl_ss)
    C <- t(solve(pop$Sigma_ss, t(left)))
    blocks[[key]] <- make_dense_block(
      pop$supports[[k]], pop$supports[[ell]], C, k, ell
    )
  }
  c_seconds <- unname(proc.time()[3] - start)
  load_start <- proc.time()[3]
  loading <- fit_population_loading(pop, blocks, rank)
  load_seconds <- unname(proc.time()[3] - load_start)
  list(
    valid = TRUE, blocks = blocks, loading = loading,
    rho = NA_real_, selected_parameter = NA_real_,
    screen = list(sets = pop$supports, recall = 1, seconds = 0),
    iterations = 0L, converged = TRUE, cumulative_pair_seconds = c_seconds,
    c_seconds = c_seconds, loading_step_seconds = load_seconds,
    canonical_seconds = c_seconds + load_seconds,
    total_seconds = c_seconds + load_seconds,
    backend = "population-Sigma", error_message = ""
  )
}

# -------------------------- EGCAR cross-validation ----------------------------

make_folds <- function(n, folds, seed) {
  set.seed(seed)
  perm <- sample.int(n)
  ids <- rep(seq_len(folds), length.out = n)
  out <- integer(n)
  out[perm] <- ids
  out
}

fit_one_egcar_cv_fold <- function(fold_id, fold_ids, X_raw, pop, rank) {
  val_idx <- which(fold_ids == fold_id)
  train_idx <- which(fold_ids != fold_id)
  split <- center_train_validation(X_raw, train_idx, val_idx)
  screen <- screen_rows(split$train, pop$supports)
  Xs <- lapply(seq_len(CFG$K), function(k) {
    split$train[[k]][, screen$sets[[k]], drop = FALSE]
  })
  stats <- lapply(seq_len(nrow(pair_index)), function(i) {
    make_pair_stats(Xs[[pair_index$k[i]]], Xs[[pair_index$ell[i]]])
  })

  states <- setNames(vector("list", nrow(pair_index)), pair_index$key)
  rho_order <- sort(RHO_GRID, decreasing = TRUE)
  rows <- vector("list", length(rho_order))

  for (j in seq_along(rho_order)) {
    rho <- rho_order[j]
    blocks <- list()
    total_iter <- 0L
    all_converged <- TRUE
    cumulative_seconds <- 0

    for (i in seq_len(nrow(pair_index))) {
      key <- pair_index$key[i]
      fit <- pair_admm(stats[[i]], rho, warm_Z = states[[key]])
      states[[key]] <- fit$C
      blocks[[key]] <- make_dense_block(
        screen$sets[[pair_index$k[i]]],
        screen$sets[[pair_index$ell[i]]],
        fit$C, pair_index$k[i], pair_index$ell[i]
      )
      total_iter <- total_iter + fit$iterations
      all_converged <- all_converged && fit$converged
      cumulative_seconds <- cumulative_seconds + fit$seconds
    }

    loading <- fit_global_loading(split$train, blocks, rank)
    score <- heldout_rayleigh_score(loading, split$validation)
    rows[[j]] <- data.frame(
      fold = fold_id, rho = rho, validation_score = score,
      selected_rows = loading$selected_total,
      screening_recall = screen$recall,
      admm_iterations = total_iter,
      all_converged = all_converged,
      cumulative_pair_seconds = cumulative_seconds,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

select_egcar_rho_cv <- function(X_raw, pop, rank, cores) {
  start <- proc.time()[3]
  folds <- min(CFG$cv_folds, nrow(X_raw[[1]]))
  fold_seed <- pop$population_seed + 7000000L + nrow(X_raw[[1]])
  fold_ids <- make_folds(nrow(X_raw[[1]]), folds, fold_seed)
  worker <- function(f) fit_one_egcar_cv_fold(f, fold_ids, X_raw, pop, rank)
  ncores <- min(folds, max(1L, cores))

  if (.Platform$OS.type == "unix" && ncores > 1L) {
    pieces <- parallel::mclapply(
      seq_len(folds), worker, mc.cores = ncores, mc.set.seed = FALSE
    )
  } else {
    pieces <- lapply(seq_len(folds), worker)
  }
  detail <- do.call(rbind, pieces)

  summary <- do.call(rbind, lapply(RHO_GRID, function(rho) {
    z <- detail[detail$rho == rho, , drop = FALSE]
    finite <- is.finite(z$validation_score)
    data.frame(
      rho = rho,
      mean_validation_score = if (any(finite)) mean(z$validation_score[finite]) else NA_real_,
      se_validation_score = if (sum(finite) > 1L) {
        stats::sd(z$validation_score[finite]) / sqrt(sum(finite))
      } else 0,
      valid_folds = sum(finite),
      mean_selected_rows = if (any(finite)) mean(z$selected_rows[finite]) else NA_real_,
      stringsAsFactors = FALSE
    )
  }))

  eligible <- summary[
    summary$valid_folds == folds & is.finite(summary$mean_validation_score),
    , drop = FALSE
  ]
  if (!nrow(eligible)) {
    eligible <- summary[
      summary$valid_folds >= ceiling(folds / 2) &
        is.finite(summary$mean_validation_score),
      , drop = FALSE
    ]
  }
  if (!nrow(eligible)) {
    selected <- min(RHO_GRID)
  } else {
    best <- max(eligible$mean_validation_score)
    selected <- max(eligible$rho[
      abs(eligible$mean_validation_score - best) <= 1e-12
    ])
  }

  list(
    rho = selected, detail = detail, summary = summary, fold_ids = fold_ids,
    tuning_seconds = unname(proc.time()[3] - start)
  )
}

# ------------------------ benchmark loading utilities -------------------------

extract_rgcca_weight_list <- function(fit, K = CFG$K) {
  candidates <- list(fit$astar, fit$a, fit$weights, fit$weight)
  for (obj in candidates) {
    if (is.list(obj) && length(obj) >= K) {
      mats <- lapply(seq_len(K), function(k) as.matrix(obj[[k]]))
      if (all(vapply(mats, nrow, integer(1)) > 0L)) return(mats)
    }
  }
  stop("Could not locate block weight matrices in the RGCCA/SGCCA fit")
}

canonicalize_raw_loading <- function(X, candidate_sets, raw_L, rank) {
  if (length(raw_L) != CFG$K) stop("raw_L must contain one matrix per block")
  L <- vector("list", CFG$K)
  for (k in seq_len(CFG$K)) {
    M <- as.matrix(raw_L[[k]])
    if (nrow(M) != length(candidate_sets[[k]])) {
      stop("Loading row count does not match candidate-set size in block ", k)
    }
    if (ncol(M) < rank) stop("Benchmark returned fewer than r components")
    L[[k]] <- M[, seq_len(rank), drop = FALSE]
    L[[k]][!is.finite(L[[k]])] <- 0
  }

  n <- nrow(X[[1]])
  Y <- lapply(seq_len(CFG$K), function(k) {
    X[[k]][, candidate_sets[[k]], drop = FALSE] %*% L[[k]]
  })
  Q <- Reduce("+", lapply(Y, function(Z) crossprod(Z) / n))
  eeq <- eigen(symm(Q), symmetric = TRUE)
  tol <- max(1, max(abs(eeq$values))) * 1e-10
  if (sum(eeq$values > tol) < rank) {
    stop("Benchmark loading matrix is rank deficient under the sample Sigma0 metric")
  }
  Qinvhalf <- eeq$vectors %*%
    (diag(1 / sqrt(pmax(eeq$values, tol)), rank) %*% t(eeq$vectors))
  L <- lapply(L, function(A) A %*% Qinvhalf)

  Y <- lapply(seq_len(CFG$K), function(k) {
    X[[k]][, candidate_sets[[k]], drop = FALSE] %*% L[[k]]
  })
  B <- matrix(0, rank, rank)
  for (i in seq_len(nrow(pair_index))) {
    k <- pair_index$k[i]; ell <- pair_index$ell[i]
    Rkl <- crossprod(Y[[k]], Y[[ell]]) / n
    B <- B + Rkl + t(Rkl)
  }
  ee <- eigen(symm(B), symmetric = TRUE)
  rotation <- ee$vectors[, seq_len(rank), drop = FALSE]
  L <- lapply(L, function(A) A %*% rotation)

  # Remove exactly zero rows while retaining original variable indices.
  sets <- vector("list", CFG$K)
  Ltrim <- vector("list", CFG$K)
  for (k in seq_len(CFG$K)) {
    keep <- which(rowSums(L[[k]]^2) > 1e-20)
    sets[[k]] <- candidate_sets[[k]][keep]
    Ltrim[[k]] <- L[[k]][keep, , drop = FALSE]
  }

  list(
    valid = sum(vapply(sets, length, integer(1))) >= rank,
    sets = sets, L = Ltrim,
    eigenvalues = ee$values[seq_len(rank)],
    selected_total = sum(vapply(sets, length, integer(1))),
    row_norms = NULL
  )
}

reconstruct_C_in_loading_spans <- function(X, loading) {
  if (is.null(loading) || !isTRUE(loading$valid)) {
    stop("Cannot reconstruct C from an invalid loading")
  }
  start <- proc.time()[3]
  n <- nrow(X[[1]])
  blocks <- list()
  for (i in seq_len(nrow(pair_index))) {
    k <- pair_index$k[i]; ell <- pair_index$ell[i]; key <- pair_index$key[i]
    Lk <- loading$L[[k]]; Ll <- loading$L[[ell]]
    Yk <- X[[k]][, loading$sets[[k]], drop = FALSE] %*% Lk
    Yl <- X[[ell]][, loading$sets[[ell]], drop = FALSE] %*% Ll
    Qk <- crossprod(Yk) / n
    Ql <- crossprod(Yl) / n
    Rkl <- crossprod(Yk, Yl) / n
    rk <- CFG$c_reconstruction_ridge * max(1, mean(diag(Qk)))
    rl <- CFG$c_reconstruction_ridge * max(1, mean(diag(Ql)))
    left <- solve(Qk + diag(rk, nrow(Qk)), Rkl)
    M <- t(solve(Ql + diag(rl, nrow(Ql)), t(left)))
    blocks[[key]] <- make_lowrank_block(
      loading$sets[[k]], loading$sets[[ell]], Lk, M, Ll, k, ell
    )
  }
  list(blocks = blocks, seconds = unname(proc.time()[3] - start))
}

call_with_supported_args <- function(fun, args) {
  fml <- setdiff(names(formals(fun)), "...")
  do.call(fun, args[names(args) %in% fml])
}

extract_scalar_parameter <- function(x, patterns = c("best", "lambda", "tau", "sparsity")) {
  if (!is.list(x)) return(NA_real_)
  nm <- names(x)
  if (is.null(nm)) return(NA_real_)
  for (pat in patterns) {
    hits <- grep(pat, nm, ignore.case = TRUE)
    for (h in hits) {
      val <- x[[h]]
      if (is.numeric(val) && length(val) == 1L && is.finite(val)) return(as.numeric(val))
    }
  }
  NA_real_
}

fit_sgca_benchmark <- function(X, pop, rank, cores, screen, seed) {
  start <- proc.time()[3]
  Xs <- lapply(seq_len(CFG$K), function(k) {
    X[[k]][, screen$sets[[k]], drop = FALSE]
  })
  pp <- vapply(Xs, ncol, integer(1))
  Xstack <- do.call(cbind, Xs)
  set.seed(seed)

  fun <- getFromNamespace("gao_gca_cv_init_and_final", "SGCAR")
  ans <- call_with_supported_args(fun, list(
    X = Xstack,
    pp = pp,
    r = rank,
    k = min(CFG$sgca_k, sum(pp)),
    lambda_grid = SGCA_LAMBDA_GRID,
    rho_scale = CFG$sgca_rho_scale,
    nfold = CFG$cv_folds,
    ncores = min(CFG$cv_folds, cores),
    parallel = cores > 1L,
    eta = CFG$sgca_eta
  ))

  raw <- ans$U_full_final
  if (is.null(raw)) raw <- ans$U_final
  if (is.null(raw)) stop("SGCA result does not contain U_full_final or U_final")
  raw <- as.matrix(raw)
  if (nrow(raw) != sum(pp) && ncol(raw) == sum(pp)) raw <- t(raw)
  if (nrow(raw) != sum(pp)) stop("Unexpected SGCA loading dimension")

  raw_list <- vector("list", CFG$K)
  offsets <- c(0L, cumsum(pp))
  for (k in seq_len(CFG$K)) {
    raw_list[[k]] <- raw[(offsets[k] + 1L):offsets[k + 1L], , drop = FALSE]
  }
  loading <- canonicalize_raw_loading(X, screen$sets, raw_list, rank)
  canonical_seconds <- unname(proc.time()[3] - start)
  rec <- reconstruct_C_in_loading_spans(X, loading)
  c_seconds <- canonical_seconds + rec$seconds

  list(
    valid = TRUE, blocks = rec$blocks, loading = loading,
    rho = NA_real_, selected_parameter = extract_scalar_parameter(ans),
    screen = screen, iterations = NA_integer_, converged = TRUE,
    cumulative_pair_seconds = NA_real_, c_seconds = c_seconds,
    loading_step_seconds = canonical_seconds,
    canonical_seconds = canonical_seconds,
    total_seconds = c_seconds,
    tuning_seconds = canonical_seconds,
    backend = "SGCAR:::gao_gca_cv_init_and_final", error_message = ""
  )
}

rgcca_fit_once <- function(blocks, method, parameter, rank) {
  K <- length(blocks)
  common <- list(
    blocks = blocks,
    connection = 1 - diag(K),
    ncomp = rep(rank, K),
    scheme = "factorial",
    scale = FALSE,
    scale_block = FALSE,
    init = "svd",
    verbose = FALSE,
    quiet = TRUE,
    comp_orth = TRUE,
    method = method
  )
  if (identical(method, "rgcca")) {
    common$tau <- rep(parameter, K)
  } else {
    mins <- vapply(blocks, function(B) 1 / sqrt(ncol(B)), numeric(1))
    common$sparsity <- pmax(rep(parameter, K), mins)
  }
  call_with_supported_args(RGCCA::rgcca, common)
}

direct_rgcca_unsupervised_cv <- function(blocks, method, rank, cores, seed) {
  n <- nrow(blocks[[1]])
  fold_ids <- make_folds(n, min(CFG$cv_folds, n), seed)
  folds <- sort(unique(fold_ids))

  fold_worker <- function(f) {
    train_idx <- which(fold_ids != f)
    val_idx <- which(fold_ids == f)
    means <- lapply(blocks, function(B) colMeans(B[train_idx, , drop = FALSE]))
    tr <- lapply(seq_along(blocks), function(k) {
      sweep(blocks[[k]][train_idx, , drop = FALSE], 2L, means[[k]], "-")
    })
    va <- lapply(seq_along(blocks), function(k) {
      sweep(blocks[[k]][val_idx, , drop = FALSE], 2L, means[[k]], "-")
    })
    rows <- vector("list", length(RGCCA_TUNING_GRID))
    local_sets <- lapply(tr, function(B) seq_len(ncol(B)))

    for (j in seq_along(RGCCA_TUNING_GRID)) {
      par <- RGCCA_TUNING_GRID[j]
      score <- tryCatch({
        fit <- rgcca_fit_once(tr, method, par, rank)
        raw <- extract_rgcca_weight_list(fit)
        loading <- canonicalize_raw_loading(tr, local_sets, raw, rank)
        heldout_rayleigh_score(loading, va)
      }, error = function(e) NA_real_)
      rows[[j]] <- data.frame(fold = f, parameter = par, score = score)
    }
    do.call(rbind, rows)
  }

  ncores <- min(length(folds), max(1L, cores))
  if (.Platform$OS.type == "unix" && ncores > 1L) {
    pieces <- parallel::mclapply(
      folds, fold_worker, mc.cores = ncores, mc.set.seed = FALSE
    )
  } else {
    pieces <- lapply(folds, fold_worker)
  }
  detail <- do.call(rbind, pieces)
  summary <- do.call(rbind, lapply(RGCCA_TUNING_GRID, function(par) {
    z <- detail[detail$parameter == par, , drop = FALSE]
    finite <- is.finite(z$score)
    data.frame(
      parameter = par,
      mean_score = if (any(finite)) mean(z$score[finite]) else NA_real_,
      valid_folds = sum(finite)
    )
  }))
  eligible <- summary[summary$valid_folds == length(folds) & is.finite(summary$mean_score), ]
  if (!nrow(eligible)) eligible <- summary[is.finite(summary$mean_score), ]
  if (!nrow(eligible)) stop("No valid RGCCA/SGCCA CV parameter")
  best <- eligible$parameter[which.max(eligible$mean_score)]
  fit_full <- rgcca_fit_once(blocks, method, best, rank)
  list(
    fit_full = fit_full, selected_parameter = best,
    detail = detail, summary = summary, backend = "direct-RGCCA-fallback"
  )
}

fit_rgcca_benchmark <- function(X, pop, rank, cores, screen, seed,
                                method = c("rgcca", "sgcca")) {
  method <- match.arg(method)
  start <- proc.time()[3]
  blocks_screened <- lapply(seq_len(CFG$K), function(k) {
    X[[k]][, screen$sets[[k]], drop = FALSE]
  })

  wrapper <- getFromNamespace("rgcca_unsupervised_cv_tau", "SGCAR")
  wrapper_args <- list(
    blocks = blocks_screened,
    lambda_values = RGCCA_TUNING_GRID,
    kfold = CFG$cv_folds,
    n_cores = min(CFG$cv_folds, cores),
    method = method,
    seed = seed,
    ncomp = rank,
    ncomponents = rank,
    r = rank,
    rank = rank
  )

  wrapped <- tryCatch(
    call_with_supported_args(wrapper, wrapper_args),
    error = function(e) e
  )

  backend <- "SGCAR:::rgcca_unsupervised_cv_tau"
  cv_detail <- cv_summary <- NULL
  selected_parameter <- NA_real_
  fit_full <- NULL

  if (!inherits(wrapped, "error")) {
    fit_full <- wrapped$fit_full
    selected_parameter <- extract_scalar_parameter(wrapped)
    raw_try <- tryCatch(extract_rgcca_weight_list(fit_full), error = function(e) NULL)
    enough <- !is.null(raw_try) && all(vapply(raw_try, ncol, integer(1)) >= rank)
    if (!enough) fit_full <- NULL
  }

  if (is.null(fit_full)) {
    direct <- direct_rgcca_unsupervised_cv(
      blocks_screened, method, rank, cores,
      seed = seed + if (method == "rgcca") 100L else 200L
    )
    fit_full <- direct$fit_full
    selected_parameter <- direct$selected_parameter
    cv_detail <- direct$detail
    cv_summary <- direct$summary
    backend <- direct$backend
  }

  raw <- extract_rgcca_weight_list(fit_full)
  loading <- canonicalize_raw_loading(X, screen$sets, raw, rank)
  canonical_seconds <- unname(proc.time()[3] - start)
  rec <- reconstruct_C_in_loading_spans(X, loading)
  c_seconds <- canonical_seconds + rec$seconds

  list(
    valid = TRUE, blocks = rec$blocks, loading = loading,
    rho = NA_real_, selected_parameter = selected_parameter,
    screen = screen, iterations = NA_integer_, converged = TRUE,
    cumulative_pair_seconds = NA_real_, c_seconds = c_seconds,
    loading_step_seconds = canonical_seconds,
    canonical_seconds = canonical_seconds,
    total_seconds = c_seconds, tuning_seconds = canonical_seconds,
    backend = backend, error_message = "",
    benchmark_cv_detail = cv_detail, benchmark_cv_summary = cv_summary
  )
}

# ------------------------------- metrics --------------------------------------

empty_failed_fit <- function(message, screen = NULL) {
  if (is.null(screen)) screen <- list(sets = NULL, recall = NA_real_, seconds = 0)
  list(
    valid = FALSE, blocks = NULL, loading = NULL, rho = NA_real_,
    selected_parameter = NA_real_, screen = screen, iterations = NA_integer_,
    converged = FALSE, cumulative_pair_seconds = NA_real_,
    c_seconds = NA_real_, loading_step_seconds = NA_real_,
    canonical_seconds = NA_real_, total_seconds = NA_real_,
    tuning_seconds = NA_real_, backend = NA_character_, error_message = message
  )
}

safe_fit <- function(label, expression, screen = NULL) {
  tryCatch(
    force(expression),
    error = function(e) {
      message(label, " failed: ", conditionMessage(e))
      empty_failed_fit(conditionMessage(e), screen)
    }
  )
}

row_support_metrics <- function(fit, pop) {
  if (!isTRUE(fit$valid) || is.null(fit$loading)) {
    return(c(row_precision = NA, row_recall = NA, row_f1 = NA, selected_rows = NA))
  }
  tp <- fp <- fn <- 0L
  for (k in seq_len(CFG$K)) {
    selected <- fit$loading$sets[[k]]
    truth <- pop$supports[[k]]
    tp <- tp + sum(selected %in% truth)
    fp <- fp + sum(!(selected %in% truth))
    fn <- fn + sum(!(truth %in% selected))
  }
  precision <- if (tp + fp > 0) tp / (tp + fp) else 0
  recall <- if (tp + fn > 0) tp / (tp + fn) else 1
  f1 <- if (precision + recall > 0) 2 * precision * recall / (precision + recall) else 0
  c(
    row_precision = precision, row_recall = recall, row_f1 = f1,
    selected_rows = sum(vapply(fit$loading$sets, length, integer(1)))
  )
}

method_metric_row <- function(method, fit, pop, common, tuning_seconds = NA_real_) {
  if (isTRUE(fit$valid) && !is.null(fit$blocks)) {
    cm <- compare_C_blocks(fit$blocks, pop$C0)
    cdist <- population_subspace_distance(fit$loading, pop)
    edist <- euclidean_loading_procrustes(fit$loading, pop)
  } else {
    cm <- list(
      c_frobenius_distance = NA_real_, c_squared_loss = NA_real_,
      c_relative_error = NA_real_, c_cosine_similarity = NA_real_,
      c0_frobenius = sqrt(2 * sum(vapply(pop$C0, function(B) sum(B$value^2), numeric(1))))
    )
    cdist <- edist <- NA_real_
  }
  sm <- row_support_metrics(fit, pop)

  data.frame(
    common,
    method = method,
    status = if (isTRUE(fit$valid)) "ok" else "failed",
    error_message = fit$error_message %||% "",
    backend = fit$backend %||% NA_character_,
    rho = fit$rho %||% NA_real_,
    selected_parameter = fit$selected_parameter %||% NA_real_,
    tau = CFG$row_tau,
    c_frobenius_distance = cm$c_frobenius_distance,
    c_squared_loss = cm$c_squared_loss,
    c_relative_error = cm$c_relative_error,
    c_cosine_similarity = cm$c_cosine_similarity,
    c0_frobenius = cm$c0_frobenius,
    canonical_subspace_distance = cdist,
    canonical_euclidean_subspace_distance = edist,
    row_precision = unname(sm["row_precision"]),
    row_recall = unname(sm["row_recall"]),
    row_f1 = unname(sm["row_f1"]),
    selected_rows = unname(sm["selected_rows"]),
    screening_recall = fit$screen$recall %||% NA_real_,
    screening_seconds = fit$screen$seconds %||% NA_real_,
    admm_iterations = fit$iterations %||% NA_integer_,
    admm_converged = fit$converged %||% NA,
    cumulative_pair_seconds = fit$cumulative_pair_seconds %||% NA_real_,
    c_seconds = fit$c_seconds %||% NA_real_,
    canonical_seconds = fit$canonical_seconds %||% NA_real_,
    total_seconds = fit$total_seconds %||% NA_real_,
    tuning_seconds = if (is.finite(tuning_seconds)) tuning_seconds else
      (fit$tuning_seconds %||% NA_real_),
    stringsAsFactors = FALSE
  )
}

`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x

# -------------------------------- worker --------------------------------------

run_worker <- function(config_id, rep_id, outdir) {
  if (config_id < 1L || config_id > nrow(CONFIGS)) stop("config_id out of range")
  if (rep_id < 1L || rep_id > N_REPS) stop("rep_id out of range")
  make_output_dirs(outdir)

  cfg <- CONFIGS[config_id, , drop = FALSE]
  n <- as.integer(cfg$n)
  p_block <- as.integer(cfg$p_per_block)
  rank <- as.integer(cfg$rank)
  p_vec <- rep(p_block, CFG$K)
  p_total <- sum(p_vec)
  cores <- available_cores()

  cat(sprintf(
    "%s: config %d/%d, rep %d, n=%d, p/block=%d, r=%d, workers=%d\n",
    EXPERIMENT, config_id, nrow(CONFIGS), rep_id, n, p_block, rank, cores
  ))

  pop <- make_population(p_vec, rank, rep_id)
  dat <- sample_views(pop, n, config_id, rep_id)
  X <- center_full(dat$X)
  common_screen <- screen_rows(X, pop$supports)
  theory_rate <- sqrt(log(p_total) / n)
  theory_rho <- CFG$theory_rho_constant * theory_rate
  seed_method <- CFG$base_seed + 20000000L + 10000L * config_id + rep_id

  oracle1 <- safe_fit("Oracle1", fit_oracle1(pop, rank))
  oracle2 <- safe_fit("Oracle2", fit_oracle2(X, pop, rank, cores))
  egcar_theory <- safe_fit(
    "EGCAR-theory",
    fit_egcar_final(X, pop, theory_rho, rank, cores, fixed_screen = common_screen),
    common_screen
  )

  egcar_cv_info <- safe_fit("EGCAR-CV tuning", {
    cv <- select_egcar_rho_cv(dat$X, pop, rank, cores)
    final <- fit_egcar_final(X, pop, cv$rho, rank, cores, fixed_screen = common_screen)
    final$tuning_seconds <- cv$tuning_seconds
    final$c_seconds <- cv$tuning_seconds + final$c_seconds
    final$canonical_seconds <- cv$tuning_seconds + final$canonical_seconds
    final$total_seconds <- final$canonical_seconds
    final$cv_detail <- cv$detail
    final$cv_summary <- cv$summary
    final
  }, common_screen)

  sgca <- safe_fit(
    "SGCA",
    fit_sgca_benchmark(X, pop, rank, cores, common_screen, seed_method + 10L),
    common_screen
  )
  rgcca <- safe_fit(
    "RGCCA",
    fit_rgcca_benchmark(
      X, pop, rank, cores, common_screen, seed_method + 20L, method = "rgcca"
    ),
    common_screen
  )
  sgcca <- safe_fit(
    "SGCCA",
    fit_rgcca_benchmark(
      X, pop, rank, cores, common_screen, seed_method + 30L, method = "sgcca"
    ),
    common_screen
  )

  fits <- list(
    `Oracle1-population` = oracle1,
    `Oracle2-support` = oracle2,
    `EGCAR-theory` = egcar_theory,
    `EGCAR-CV` = egcar_cv_info,
    SGCA = sgca,
    RGCCA = rgcca,
    SGCCA = sgcca
  )

  common <- data.frame(
    experiment = EXPERIMENT,
    config_id = config_id,
    panel = as.character(cfg$panel),
    n = n,
    p_per_block = p_block,
    p_total = p_total,
    rank = rank,
    rep_id = rep_id,
    K = CFG$K,
    support_per_block = CFG$support_size,
    signal_lambda = CFG$signal_lambda,
    toeplitz_alpha = CFG$alpha,
    theory_rate = theory_rate,
    theory_rho = theory_rho,
    cv_grid = paste(format(RHO_GRID, scientific = TRUE), collapse = ";"),
    population_seed = pop$population_seed,
    data_seed = dat$data_seed,
    screen_sizes = paste(vapply(common_screen$sets, length, integer(1)), collapse = ";"),
    common_screen_recall = common_screen$recall,
    stringsAsFactors = FALSE
  )

  metrics <- do.call(rbind, lapply(names(fits), function(method) {
    method_metric_row(method, fits[[method]], pop, common)
  }))

  tag <- sprintf("cfg%03d_rep%02d", config_id, rep_id)
  atomic_csv(metrics, file.path(outdir, "metrics", paste0("metrics_", tag, ".csv")))

  if (!is.null(egcar_cv_info$cv_detail)) {
    d <- egcar_cv_info$cv_detail
    d$config_id <- config_id; d$rep_id <- rep_id; d$n <- n
    d$p_per_block <- p_block; d$rank <- rank
    atomic_csv(d, file.path(outdir, "cv", paste0("egcar_cv_detail_", tag, ".csv")))
  }
  if (!is.null(egcar_cv_info$cv_summary)) {
    s <- egcar_cv_info$cv_summary
    s$config_id <- config_id; s$rep_id <- rep_id; s$n <- n
    s$p_per_block <- p_block; s$rank <- rank
    atomic_csv(s, file.path(outdir, "cv", paste0("egcar_cv_summary_", tag, ".csv")))
  }

  snapshot <- list(
    config = cfg,
    rep_id = rep_id,
    population = list(
      p_vec = pop$p_vec, rank = pop$rank, supports = pop$supports,
      U_active = pop$U_active, Sigma_ss = pop$Sigma_ss,
      C0 = pop$C0, L0 = pop$L0
    ),
    screen = common_screen,
    fits = fits,
    metrics = metrics,
    package_versions = vapply(
      c(CRAN_PACKAGES, "SGCAR"),
      function(pkg) as.character(utils::packageVersion(pkg)), character(1)
    )
  )

  if (p_total <= CFG$full_matrix_save_max) {
    snapshot$full_matrices <- lapply(fits, function(fit) {
      if (!isTRUE(fit$valid) || is.null(fit$blocks)) return(NULL)
      list(
        C = materialize_C(fit$blocks, p_vec),
        L = {
          out <- matrix(0, p_total, rank)
          offsets <- c(0L, cumsum(p_vec))
          for (k in seq_len(CFG$K)) {
            rr <- offsets[k] + fit$loading$sets[[k]]
            out[rr, ] <- fit$loading$L[[k]]
          }
          out
        }
      )
    })
    snapshot$full_matrices$C0 <- materialize_C(pop$C0, p_vec)
    population_covariances <- materialize_population_covariances(pop)
    snapshot$full_matrices$Sigma <- population_covariances$Sigma
    snapshot$full_matrices$Sigma0 <- population_covariances$Sigma0
    snapshot$full_matrices$L0 <- {
      out <- matrix(0, p_total, rank)
      offsets <- c(0L, cumsum(p_vec))
      for (k in seq_len(CFG$K)) {
        rr <- offsets[k] + pop$supports[[k]]
        out[rr, ] <- pop$U_active[[k]] / sqrt(CFG$K)
      }
      out
    }
  }

  atomic_rds(snapshot, file.path(outdir, "fits", paste0("fit_", tag, ".rds")))
  cat("Saved ", tag, "\n", sep = "")
}

# ------------------------------- aggregation ---------------------------------

mean_se <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(c(mean = NA_real_, se = NA_real_, n = 0))
  c(
    mean = mean(x),
    se = if (length(x) > 1L) stats::sd(x) / sqrt(length(x)) else 0,
    n = length(x)
  )
}

summarize_metrics <- function(M) {
  numeric_metrics <- c(
    "c_frobenius_distance", "c_squared_loss", "c_relative_error",
    "c_cosine_similarity", "canonical_subspace_distance",
    "canonical_euclidean_subspace_distance", "row_precision", "row_recall",
    "row_f1", "selected_rows", "screening_recall", "c_seconds",
    "canonical_seconds", "total_seconds", "tuning_seconds"
  )
  keys <- unique(M[c(
    "experiment", "panel", "n", "p_per_block", "p_total", "rank", "method", "theory_rate", "theory_rho"
  )])
  rows <- vector("list", nrow(keys))
  for (i in seq_len(nrow(keys))) {
    keep <- rep(TRUE, nrow(M))
    for (nm in names(keys)) keep <- keep & M[[nm]] == keys[[nm]][i]
    z <- M[keep, , drop = FALSE]
    row <- keys[i, , drop = FALSE]
    for (metric in numeric_metrics) {
      ms <- mean_se(z[[metric]])
      row[[paste0(metric, "_mean")]] <- unname(ms["mean"])
      row[[paste0(metric, "_se")]] <- unname(ms["se"])
      row[[paste0(metric, "_n")]] <- unname(ms["n"])
    }
    row$successful_repetitions <- sum(z$status == "ok")
    rows[[i]] <- row
  }
  do.call(rbind, rows)
}

plot_metric <- function(S, rank_value, xvar, metric, ylab, path,
                        facet = TRUE, log_x = TRUE, log_y = TRUE,
                        add_rate = FALSE) {
  D <- S[S$rank == rank_value, , drop = FALSE]
  mean_col <- paste0(metric, "_mean")
  se_col <- paste0(metric, "_se")
  D$mean_value <- D[[mean_col]]
  D$se_value <- D[[se_col]]
  D$method <- factor(D$method, levels = METHOD_LEVELS)
  D$plot_value <- if (log_y) pmax(D$mean_value, 1e-12) else D$mean_value
  D$lower <- if (log_y) pmax(D$mean_value - D$se_value, 1e-12) else
    D$mean_value - D$se_value
  D$upper <- D$mean_value + D$se_value

  p <- ggplot(D, aes_string(x = xvar, y = "plot_value", color = "method")) +
    geom_line(linewidth = 0.7, na.rm = TRUE) +
    geom_point(size = 2, na.rm = TRUE) +
    geom_errorbar(
      aes(ymin = lower, ymax = upper), width = 0, linewidth = 0.35, na.rm = TRUE
    ) +
    labs(
      x = if (xvar == "n") "Sample size n" else "Features per block",
      y = ylab, color = "Method",
      title = paste0(ylab, ", r = ", rank_value),
      caption = if (log_y) "Values at numerical zero are displayed at 1e-12." else NULL
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom", legend.box = "vertical")

  if (log_x) p <- p + scale_x_log10()
  if (log_y) p <- p + scale_y_log10()
  if (facet && length(unique(D$panel)) > 1L) {
    p <- p + facet_wrap(~panel, scales = "free_x")
  }

  if (add_rate) {
    rate <- unique(D[c("panel", "n", "p_per_block", "p_total", "theory_rate")])
    names(rate)[names(rate) == "theory_rate"] <- "plot_value"
    p <- p +
      geom_line(
        data = rate,
        aes_string(x = xvar, y = "plot_value", group = "panel"),
        inherit.aes = FALSE, linetype = 2, linewidth = 0.9, color = "black"
      ) +
      labs(subtitle = "Dashed black curve: sqrt(log(p_total)/n)")
  }

  ggsave(path, p, width = 8.5, height = 5.8)
  ggsave(sub("\\.pdf$", ".png", path), p, width = 8.5, height = 5.8, dpi = 180)
  p
}

aggregate_results <- function(outdir) {
  make_output_dirs(outdir)
  files <- list.files(
    file.path(outdir, "metrics"), pattern = "^metrics_.*\\.csv$", full.names = TRUE
  )
  if (!length(files)) stop("No metric files found in ", file.path(outdir, "metrics"))
  M <- do.call(rbind, lapply(files, utils::read.csv, check.names = FALSE))
  M$method <- as.character(M$method)
  atomic_csv(M, file.path(outdir, "all_metrics.csv"))
  S <- summarize_metrics(M)
  atomic_csv(S, file.path(outdir, "summary_mean_over_10_repetitions.csv"))

  xvar <- if (EXPERIMENT == "vary_n") "n" else "p_per_block"
  all_plots <- list()
  for (r in c(1L, 2L, 5L)) {
    plist <- list(
      plot_metric(
        S, r, xvar, "c_frobenius_distance",
        "Frobenius distance ||C-hat - C0||F",
        file.path(outdir, "plots", sprintf("C_frobenius_rank_%d.pdf", r))
      ),
      plot_metric(
        S, r, xvar, "c_squared_loss",
        "Squared Frobenius loss ||C-hat - C0||F^2",
        file.path(outdir, "plots", sprintf("C_squared_loss_rank_%d.pdf", r))
      ),
      plot_metric(
        S, r, xvar, "c_relative_error",
        "Relative C error",
        file.path(outdir, "plots", sprintf("C_relative_with_rate_rank_%d.pdf", r)),
        add_rate = TRUE
      ),
      plot_metric(
        S, r, xvar, "c_seconds",
        "Time through C construction (seconds)",
        file.path(outdir, "plots", sprintf("C_time_rank_%d.pdf", r))
      ),
      plot_metric(
        S, r, xvar, "canonical_subspace_distance",
        "Sigma0 sine-theta loading-subspace distance",
        file.path(outdir, "plots", sprintf("canonical_error_rank_%d.pdf", r))
      ),
      plot_metric(
        S, r, xvar, "canonical_euclidean_subspace_distance",
        "Euclidean sine-theta loading-subspace distance",
        file.path(outdir, "plots", sprintf("canonical_euclidean_error_rank_%d.pdf", r))
      ),
      plot_metric(
        S, r, xvar, "canonical_seconds",
        "Time through canonical directions (seconds)",
        file.path(outdir, "plots", sprintf("canonical_time_rank_%d.pdf", r))
      )
    )
    all_plots[[as.character(r)]] <- plist
    pdf(file.path(outdir, "plots", sprintf("all_metrics_rank_%d.pdf", r)),
        width = 8.5, height = 5.8, onefile = TRUE)
    for (p in plist) print(p)
    dev.off()
  }

  status <- aggregate(
    rep_id ~ method + panel + n + p_per_block + rank,
    data = M[M$status == "ok", , drop = FALSE], FUN = length
  )
  names(status)[names(status) == "rep_id"] <- "successful_repetitions"
  atomic_csv(status, file.path(outdir, "successful_repetitions.csv"))

  writeLines(
    c(
      paste("Experiment:", EXPERIMENT),
      paste("Metric files found:", length(files), "of expected", EXPECTED_TASKS),
      paste("Generated:", format(Sys.time())),
      paste("R:", R.version.string),
      paste("Library paths:", paste(.libPaths(), collapse = " | "))
    ),
    file.path(outdir, "AGGREGATION_COMPLETE.txt")
  )
  cat("Aggregation complete: ", outdir, "\n", sep = "")
}

# -------------------------------- dispatch ------------------------------------

usage <- function() {
  cat(
    "Usage:\n",
    "  Rscript SCRIPT.R check_packages\n",
    "  Rscript SCRIPT.R install_packages\n",
    "  Rscript SCRIPT.R worker <config_id> <rep_id> <out_dir>\n",
    "  Rscript SCRIPT.R aggregate 0 0 <out_dir>\n",
    sep = ""
  )
}

if (identical(MODE, "worker")) {
  if (length(ARGS) < 4L) { usage(); stop("Missing worker arguments") }
  run_worker(as.integer(ARGS[[2]]), as.integer(ARGS[[3]]), ARGS[[4]])
} else if (identical(MODE, "aggregate")) {
  if (length(ARGS) < 4L) { usage(); stop("Missing aggregation output directory") }
  aggregate_results(ARGS[[4]])
} else {
  usage()
}
