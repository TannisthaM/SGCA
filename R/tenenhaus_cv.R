library(MASS)
library(stats)
library(pracma)
library(tidyverse)
library(ggplot2)
library(Matrix)
library(tidyr)
library(geigen)
library(RGCCA)
library(dplyr)
library(expm)
library(future.apply)
library(foreach)
library(doParallel)

rgcca_holdout_score <- function(Y_test, connection = NULL, scheme = "factorial", bias = TRUE, comp = 1) {
  J <- length(Y_test)
  if (is.null(connection)) connection <- 1 - diag(J)
  
  g <- if (is.function(scheme)) scheme else switch(
    scheme,
    horst    = function(x) x,
    centroid = function(x) abs(x),
    factorial= function(x) x^2,
    stop("Unknown scheme: ", scheme)
  )
  
  cov_b <- function(u, v) {
    u <- u[, comp]; v <- v[, comp]
    u <- u - mean(u); v <- v - mean(v)
    denom <- if (bias) length(u) else (length(u) - 1)
    sum(u * v) / denom
  }
  
  score <- 0
  for (j in 1:(J - 1)) for (k in (j + 1):J) {
    if (connection[j, k] != 0) {
      score <- score + connection[j, k] * g(cov_b(Y_test[[j]], Y_test[[k]]))
    }
  }
  score
}


rgcca_unsupervised_cv_tau <- function(
    blocks, lambda_values,
    connection = NULL,
    scheme = "factorial",
    ncomp = 1,
    kfold = 5,
    n_cores = max(1, parallel::detectCores() - 1),
    seed = 1,
    scale = TRUE,
    scale_block = TRUE,
    bias = TRUE,
    method = c("rgcca", "sgcca")   # NEW
) {
  method <- match.arg(method)      # NEW
  
  # ---- ENSURE BLOCK NAMES ----
  if (is.null(names(blocks)) || any(names(blocks) == "")) {
    names(blocks) <- paste0("block", seq_along(blocks))
  }
  block_names <- names(blocks)
  
  # Force matrices + assign unique colnames per block
  blocks <- lapply(seq_along(blocks), function(j) {
    nm <- names(blocks)[j]
    B  <- as.matrix(blocks[[j]])
    
    if (is.null(colnames(B))) {
      colnames(B) <- paste0(nm, "_V", seq_len(ncol(B)))
    } else {
      # still make them unique + prefixed to avoid cross-block duplicates
      colnames(B) <- paste0(nm, "_", make.unique(colnames(B)))
    }
    
    B
  })
  
  J <- length(blocks)
  n <- nrow(blocks[[1]])
  if (is.null(connection)) connection <- 1 - diag(J)
  
  set.seed(seed)
  fold_id <- sample(rep(seq_len(kfold), length.out = n))
  
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  future::plan(future::multisession, workers = n_cores)
  
  fold_scores <- future.apply::future_lapply(
    seq_len(kfold),
    function(f) {
      tr <- which(fold_id != f)
      te <- which(fold_id == f)
      
      train_blocks <- lapply(blocks, function(B) B[tr, , drop = FALSE])
      test_blocks  <- lapply(blocks, function(B) B[te, , drop = FALSE])
      
      names(train_blocks) <- block_names
      names(test_blocks)  <- block_names
      
      sapply(lambda_values, function(tau_scalar) {
        tau_vec <- rep(tau_scalar, J)
        
        fit <- tryCatch(
          RGCCA::rgcca(
            blocks = train_blocks,
            connection = connection,
            method = method,              # CHANGED (was "rgcca")
            tau = tau_vec,
            ncomp = ncomp,
            scheme = scheme,
            scale = scale,
            scale_block = scale_block,
            bias = bias,
            verbose = FALSE
          ),
          error = function(e) {
            message("rgcca() failed: fold=", f, " tau=", tau_scalar,
                    " method=", method, " :: ", conditionMessage(e))
            NULL
          }
        )
        if (is.null(fit)) return(NA_real_)
        
        trans <- tryCatch(
          RGCCA::rgcca_transform(fit, blocks_test = test_blocks),
          error = function(e) {
            message("rgcca_transform() failed: fold=", f, " tau=", tau_scalar,
                    " method=", method, " :: ", conditionMessage(e))
            NULL
          }
        )
        if (is.null(trans)) return(NA_real_)
        
        rgcca_holdout_score(trans, connection = connection, scheme = scheme, bias = bias, comp = 1)
      })
    },
    future.seed = TRUE,
    future.packages = c("RGCCA")   # crucial for multisession workers
  )
  
  score_mat <- do.call(rbind, fold_scores)  # kfold x length(lambda_values)
  cv_mean <- colMeans(score_mat, na.rm = TRUE)
  cv_sd   <- apply(score_mat, 2, sd, na.rm = TRUE)
  
  best_idx <- which.max(replace(cv_mean, is.na(cv_mean), -Inf))
  best_tau <- lambda_values[best_idx]
  
  fit_full <- RGCCA::rgcca(
    blocks = blocks,
    connection = connection,
    method = method,               # CHANGED (was "rgcca")
    tau = rep(best_tau, J),
    ncomp = ncomp,
    scheme = scheme,
    scale = scale,
    scale_block = scale_block,
    bias = bias,
    verbose = TRUE
  )
  
  list(
    lambda_values = lambda_values,
    fold_scores = score_mat,
    cv_mean = cv_mean,
    cv_sd = cv_sd,
    best_tau = best_tau,
    fit_full = fit_full,
    method = method               # NEW (handy for debugging)
  )
}


