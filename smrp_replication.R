# ------------------ Estimation of VAR(1) Parameters ------------------

# Estimate Recursion Matrix A
#
# Load corpcor library for pseudoinverse function
library("corpcor")
estimate_a <- function(s) {
  n <- ncol(s)
  #print(n)
  m <- nrow(s)
  #print(m)
  mat1 <- matrix(0,n,n)
  mat2 <- matrix(0,n,n)
  # Calculate the two summation matrices
  for(t in 2:m) {
    t1_iter = s[t,] %*% t(s[t-1,])
    mat1 <- mat1 + t1_iter
    t2_iter = s[t-1,] %*% t(s[t-1,]) 
    mat2 <- mat2 + t2_iter
  }
  # Take the Moore-Penrose pseudoinverse of the second matrix
  mpi <- pseudoinverse(mat2)
  # Multiply the matrices
  a <- mat1 %*% mpi
  # cat("\n")
  # print(nrow(a))
  # print(ncol(a))
  return(a)
}

# Estimate Sigma
#
estimate_sigma <- function(s,a) {
  n <- ncol(s)
  m <- nrow(s)
  coeff <- 1/(n*(m-1))
  summation <- 0
  for(t in 2:m) {
    iter <- ((s[t,] - (a %*% s[t-1,]))^2)
    summation <- summation + iter
  }
  sigma <- sqrt(coeff*summation)
  return(sigma)
}

# Estimate Covariance Matrix K of Noise
#
estimate_k <- function(s,a) {
  n <- ncol(s)
  m <- nrow(s)
  summation <- 0
  coeff <- 1/(m-1)
  for(t in 2:m) {
    # print(nrow(t(s[t-1,])))
    # print(ncol(t(s[t-1,])))
    iter <- (s[t,] - (a %*% s[t-1,]))
    summation <- summation + (iter %*% t(iter))
  }
  k = coeff * summation
  return(k)
}

# Estimate stationary covariance matrix G 
#
# For now, we assume Delta = 0.1, iteration # = 1000. 
# TODO: figure out how to choose these constants
#
estimate_g <- function(s,k,a) {
  # Take G(0) as sample covariance matrix
  g <- cov(s)
  delta <- 0.1
  iter_len <- 1000
  for(i in 1:iter_len) {
    g <- g - delta * (g - (a %*% g %*% t(a)) - k)
  }
  return(g)
}

# -------------------- Portfolio Selection Methods --------------------

# Generalized Eigenvalue Problem
# ------------------------------
# AGA(T)x = (lambda)Gx,
# where lambda is the mean reversion parameter
#
# solve(a,b,...) solves the equation a %*% x = b for x
#
# eigen(x) computes the eigenvalues and eigenvectors of x
# Returns a list with components:
#   values  - a vector containing the eigenvalues of x in decreasing order
#   vectors - a matrix whose columns contain the eigenvectors of x 
#
# estimate_eigen <- function(a,g) {
#   aga <- a %*% g %*% t(a)
#   s <- solve(g,aga)
#   eigen_ret <- eigen(s)
#   cat("\n")
#   # print("eigen ret")
#   # print(eigen_ret)
#   eigenvalue <- eigen_ret$values[1]
#   eigenvector <- eigen_ret$vectors[1]
#   c(eigen_ret$value[1],eigen_ret$vectors[1])
# }
estimate_eigen <- function(a,b) {
  s <- solve(b,a)
  eigen_ret <- eigen(s)
  eigenvalue <- eigen_ret$values[1]
  eigenvector <- eigen_ret$vectors[1]
  c(eigenvalue,eigenvector)
}


# -------------- Constrained Optimization Problem --------------
# 
#         x(opt) =  arg max     (x(T)AGA(T)x)/(x(T)Gx
#               x(E)R^n,card(x)<=L
#
# We explore 4 methods of solving this problem - 
# Exhaustive Search, Greedy Method, Truncation Method and Simulated Annealing

# Exhaustive Search Method
# -------------------------
# Brute force method of constructing all L-dimensional submatrices of G 
# and AGA(T) and solving each corresponding eigenvalue problem
#
exhaustive_search <- function(a,g,l) {
  combinations <- combn(1:ncol(a),l)
  iter_len = ncol(combinations)
  opt <- 0
  for(i in 1:iter_len) {
    # cat("\n")
    # print(a[combinations[,i],combinations[,i]])
    curr_a <- a[combinations[,i],combinations[,i]]
    curr_g <- g[combinations[,i],combinations[,i]]
    curr_aga <- curr_a %*% curr_g %*% t(curr_a)
    curr_opt <- estimate_eigen(curr_aga,curr_g)
    # cat("\n")
    # print("Curr_opt[1]")
    # print(curr_opt[1])
    # cat("\n")
    # print("opt[1]")
    # print(opt[1])
    if(curr_opt[1]>opt[1]){
      opt <- curr_opt
    }
  }
  opt
}

# Greedy Method
# -------------------------
#
greedy <- function(a,g,l) {
  n <- ncol(a)
  opt <- 0
  x <- c()
  if(l>n){
    ret(-1)
  }
  #while(length(x)<l){
    for(i in 1:n){
      curr_a <- a[,combinations[,i]]
      curr_g <- g[,combinations[,i]]
      curr_aga <- curr_a %*% curr_g %*% t(curr_a)
      curr_opt <- estimate_eigen(curr_aga,curr_g)
      if(curr_opt[1]>opt[1]){
        opt <- curr_opt
      }
    }
  #}
}

calculate.eigen <- function(a.est, g.est, cols) {
  AGAt <- as.matrix(a.est) %*% g.est %*% t(a.est)
  sol <- eigen(solve(g.est,AGAt))
  x <- rep(0, ncol(a.est))
  x[cols] <- sol$vectors[,1]
}

truncated.search <- function(a.est, g.est, L) {
  x_opt <- calculate.eigen(a.est, g.est, ncol(a.est))
  calculate.eigen(a.est, g.est, order(x_opt,decreasing=TRUE)[1:L])
}

library(readxl)
data <- read_excel("final_train.xlsx")
data <- as.matrix(data)
  
a <- estimate_a(data)
# cat("\n")
# print(a)

sigma <- estimate_sigma(data,a)
cat("\n")
print(sigma)

k <- estimate_k(data,a)
# cat("\n")
# print(k)

g <- estimate_g(data,k,a)
# cat("\n")
# print(g)

# exh <- exhaustive_search(a,g,4)
# cat("\n")
# print(exh)
# 
# trunc <- truncated.search(a,g,4)
# cat("\n")
# print(trunc)
# 
# s <- 0
# for(i in 1:length(trunc)) {
#   s <- s + trunc[i]
# }
# cat("\n")
# print(s)

# s <- 0
# for(i in 1:length(exh)) {
#   s <- s + exhaustive[i]
# }
# cat("\n")
# print(s)

