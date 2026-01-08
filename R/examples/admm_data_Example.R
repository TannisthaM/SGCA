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


###############################Load/copy-paste content from DGP.R######################################################################

lambda_values <- c(0, 1e-5,1e-4,1e-3,1e-2, 0.1, 1,10,100,1000,1e+4,1e+5)

cv_fit <- cv_admm_sgca_C(
          X, p_list, lambdas = lambda_values, r = 1, K = 5,
          penalty = "l1", penalize = "all", loss_part = "all", relative_loss = TRUE,
          rho = 1, max_iter = 200, adapt_rho = TRUE,   # used only if admm_sgca supports them
          parallel = TRUE, nb_cores = 1, blas_threads = 1,
          verbose = TRUE
        )


fit <- admm_sgca(S, sigma0hat, cv_fit$lambda_min, r,
                         rho = 1,
                         p_list = pp,
                         penalty = c("l1"),
                         penalize = c("all"),   # used for l1 and l21_rows
                         
                         max_iter = 200,
                         abs_tol = 1e-18, rel_tol = 1e-17,
                         adapt_rho = FALSE, mu = 10, tau_incr = 2, tau_decr = 2,
                         sparsity_threshold = 1e-4,
                         verbose = TRUE)

