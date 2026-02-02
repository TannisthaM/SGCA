library(SGCAR)
library(MASS)
library(Matrix)
library(future.apply)
library(PMA)
library(pracma)
library(RGCCA)
library(geigen)
library(CVXR)
library(cvCovEst)
library(Matrix)

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

lambda_values <- c(0,1e-7,1e-6, 1e-5,1e-4,1e-3,1e-2, 0.1, 1,10,100,1000,1e+4,1e+5,1e+6,1e+7)
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
oracle$time <- c(); oracle$loss_ris30 <- c(); oracle$loss_ris1 <- c()
oracle2$time <- c(); oracle2$loss <- c()

C0_Cmat <- list()
cv_fit  <- c()
rho_grid=c(1e-12,1e-11,1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1,10,
           100,1000,10000,1e+5)



# pp is your p_list, X is n x p
idxs <- .block_indices(pp)

# Ledoit–Wolf shrinkage per block (each returns a covariance matrix)
Sigma0_blocks_LW <- lapply(idxs, function(id) {
  cvCovEst::linearShrinkLWEst(dat = X[, id, drop = FALSE])
})

# Reassemble block diagonal Sigma0


admm=list();temp=c()
for(j in 1:5){
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

#sigma0hat[idx1,idx1]=linearShrinkLWEst(sigma0hat[idx1,idx1])
#sigma0hat[idx2,idx2]=linearShrinkLWEst(sigma0hat[idx2,idx2])

#sigma0hat[idx3,idx3]=linearShrinkLWEst(sigma0hat[idx3,idx3])
Sigma0hat_LW <- as.matrix(Matrix::bdiag(Sigma0_blocks_LW))
Sigma0hat_LW <- (Sigma0hat_LW + t(Sigma0hat_LW)) / 2  # safety symmetrize
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
fit_admm_GT <- .admm_sgca_run(
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
admm[[j]]=min(sqrt(mean((fit_admm_GT$U[,1] - a)^2)),sqrt(mean((fit_admm_GT$U[,1] + a)^2)))
print(j)
}
fit_admm_GT <- .admm_sgca_run(
  prep = SGCAR:::.admm_sgca_prepare(
    Sigma  = S,
    Sigma0 = linearShrinkLWEst(X)*Mask,
    rho = NA,
    p_list = p_list,
    penalty = "l1",
    penalize = "all",
    symmetrize_z = TRUE
  ),
  lambda = lambda_values[which.min(admm)],
  warm_start =  "none",
  r = 1,
  verbose=TRUE,
  adapt_rho = TRUE,
  compute_canon = TRUE,
  abs_tol = 1e-10,
  rel_tol = 1e-10,
  max_iter = 100000
)
min(sqrt(mean((fit_admm_GT$U[,1] - a)^2)),sqrt(mean((fit_admm_GT$U[,1] + a)^2)))
