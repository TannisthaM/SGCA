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

gca_to_cca <-
  function(a_estimate, S, pp){
    p1 = pp[1];
    p2 = pp[2];
    p = p1 + p2;
    sigmaxhat = S[1:p1,1:p1];
    sigmayhat = S[(p1+1):p,(p1+1):p];
    u_estimate = a_estimate[1:p1,] %*% pracma::sqrtm(t(a_estimate[1:p1,]) %*% sigmaxhat %*% a_estimate[1:p1,])$Binv;
    v_estimate = a_estimate[(p1+1):p,] %*% pracma::sqrtm(t(a_estimate[(p1+1):p,]) %*% sigmayhat %*% a_estimate[(p1+1):p,])$Binv;
    l = list("u" = u_estimate, "v" = v_estimate)
    return(l)
  }


Soft <-
  function(a,b){
    if(b<0) stop("Can soft-threshold by a nonnegative quantity only.")
    return(sign(a)*pmax(0,abs(a)-b))
  }

updatePi <-
  function(B,sqB,A,H,Gamma,nu,rho,Pi,tau){
    C <- Pi + 1/tau*A-nu/tau*B%*%Pi%*%B+nu/tau*sqB%*%(H-Gamma)%*%sqB
    D <- rho/tau
    return(Soft(C,D))
  }

updateH <-
  function(sqB,Gamma,nu,Pi,K){
    
    temp <- 1/nu * Gamma + sqB%*%Pi%*%sqB
    temp <- (temp+t(temp))/2
    svdtemp <- eigen(temp)
    d <- svdtemp$values
    p <- length(d)
    if(sum(pmin(1,pmax(d,0)))<=K){
      dfinal <- pmin(1,pmax(d,0))
      return(svdtemp$vectors%*%diag(dfinal)%*%t(svdtemp$vectors))
    }
    fr <- function(x){
      sum(pmin(1,pmax(d-x,0)))
    }
    # Vincent Vu Fantope Projection
    knots <- unique(c((d-1),d))
    knots <- sort(knots,decreasing=TRUE)
    temp <- which(sapply(knots,fr)<=K)
    lentemp <- tail(temp,1)
    a=knots[lentemp]
    b=knots[lentemp+1]
    fa <- sum(pmin(pmax(d-a,0),1))
    fb <- sum(pmin(pmax(d-b,0),1))
    theta <- a+ (b-a)*(K-fa)/(fb-fa)
    dfinal <- pmin(1,pmax(d-theta,0))
    res <- svdtemp$vectors%*%diag(dfinal)%*%t(svdtemp$vectors)
    return(res)
  }

hard <-
  function(U, k){
    if(r>1){
      truncate.value <- sort(apply(U, 1, FUN = function(x) sum(x^2)),decreasing=TRUE)[k]
      U[which(apply(U, 1, FUN = function(x) sum(x^2))<truncate.value)] <- 0
    }else{
      truncate.value <- sort(abs(U),decreasing=TRUE)[k]
      U[which(abs(U)<truncate.value)] <- 0
    }
    return(U)
  }



# function to convert Pi from Fantope to input of TGD
# Inputs:
# =======
# Pi:         Output of sgca_init$Pi
# r:          Latent dimension

# Outputs:
# ========
# ainit:   Initialization for the generalized eigenspace
init_process <-
  function(Pi, r){
    ainit = svd(Pi)
    uest <- ainit$u
    dest <- diag(ainit$d)
    if (r == 1){
      ainit <- uest[,1] * sqrt(dest[1:r,1:r])
    } else
      ainit <- uest[,1:r] %*% pracma::sqrtm(dest[1:r,1:r])$B
    return(ainit)
  }



# Function for Intialization via generalized fantope projection
# Inputs:
# =======
# A, B:       Pair of matrix to calculate generalized eigenspace (sigma, sigma0)
# nu:         Parameter of ADMM, default set to 1
# K:          nuclear norm constraint, equal to r
# rho:     penalty parameter on the l_1 norm of the solution, scaled by
#             sqrt(log(max(p1,p2))/n)
# epsilon:    tolerance level for convergence in ADMM
# maxiter:    maximum number of iterations in ADMM
# trace:      if set to True will print all iterations 

# Outputs:
# ========
# $Pi:     optimum of the convex program

sgca_init <-
  function(A,B,rho,K,nu=1,epsilon=5e-3,maxiter=1000,trace=FALSE){
    p <- nrow(B)
    eigenB <- eigen(B)
    sqB <- eigenB$vectors%*%sqrt(diag(pmax(eigenB$values,0)))%*%t(eigenB$vectors)	
    tau <- 4*nu*eigenB$values[1]^2	
    criteria <- 1e10
    i <- 1
    # Initialize parameters
    H <- Pi <- oldPi <-  diag(1,p,p)
    Gamma <- matrix(0,p,p)
    # While loop for the iterations
    while(criteria > epsilon && i <= maxiter){
      for (i in 1:20){
        Pi <- updatePi(B,sqB,A,H,Gamma,nu,rho,Pi,tau)
      }
      #Pi <- updatePi(B,sqB,A,H,Gamma,nu,lambda,Pi,tau)
      
      H <- updateH(sqB,Gamma,nu,Pi,K)
      Gamma <- Gamma + (sqB%*%Pi%*%sqB-H) * nu	
      criteria <- sqrt(sum((Pi-oldPi)^2))
      oldPi <- Pi
      i <- i+1
      if(trace==TRUE)
      {
        print(i)
        print(criteria)
      }
    }
    return(list(Pi=Pi,H=H,Gamma=Gamma,iteration=i,convergence=criteria))
    
  }

# Function for Thredholded Gradient Descent
# Inputs:
# =======
# A, B:         Pair of matrix to calculate generalized eigenspace (sigma, sigma0)
# r:            Latent dimension
# init:         Initialize estimator of generlized eigenspace, obtained from sgca_int. Need to be k-sparse
# lambda:       penalty parameter of Lagrangian function f(L), default to be 0.01
# convergence:  tolerance level for convergence in TGD
# maxiter:      maximum number of iterations in TGD
# plot:         if set to True will plot intermediate iterations, need to specify scale variable V (as shown in paper)
#               default set to False

# Outputs:
# ========
# final_estimate:   final estimation of leading r sparse generalized eigenspace

sgca_tgd <-
  function(A, B, r, init, k, lambda = 0.01, eta=0.01, convergence=1e-3, maxiter=10000, plot = FALSE){
    #perform hard thresholding
    init <- hard(init, k)
    u <- init
    criteria <- 1e10
    iter <- 0
    error <- rep(0, maxiter)
    # renormalization 
    ut <- init %*% pracma::sqrtm(diag(r)+t(u) %*% A %*% u/lambda)$B;
    
    while(criteria > convergence && iter <= maxiter){
      #perform gradient descent
      grad <- -A %*% ut + lambda * B %*% ut %*% (t(ut) %*% B %*% ut- diag(r));
      vt <- ut - eta * grad
      
      
      # Perform truncation
      vt <- hard(vt, k)
      
      criteria <- sqrt(sum((ut-vt)^2))
      
      ut <- vt
      iter <- iter+1
      if (plot){
        error[iter] <- subdistance(vt, scale)
      }
    }
    if (plot){
      plot(error[1:iter], type='l',  main="Matrix distance and iterations", 
           xlab="Number of iterations", ylab="Matrix distance",)
    }
    final_estimate <- ut %*% pracma::sqrtm(t(ut) %*% B %*% ut)$Binv
    return(final_estimate)
  }

