# Import ------------------------------------------------------------------

import <- function(){
  packages = c('clrdag', 'car', 'svMisc', 'MASS', 'readxl')
  lapply(packages, require, character.only = TRUE)
  setwd("~/Desktop/HW3")
}

load_data <- function(){
  sheets <- c("1. cd3cd28", "2. cd3cd28icam2", "3. cd3cd28+aktinhib","4. cd3cd28+g0076",
              "5. cd3cd28+psitect","6. cd3cd28+u0126","7. cd3cd28+ly", "8. pma", "9. b2camp")
  
  l <- vector("list", 9)
  for(i in 1:9) l[[i]] <- data.frame(read_excel("cytometry-data.xlsx", sheet = sheets[i]))
  
  return(l)
}

# Estimates ---------------------------------------------------------------

sigma_sum <- function(A, X, i, j, idx){
  
  # this function computes the internal sum of the sigma estimate formula
  
  s = 0
  for(k in idx[idx != j]) s = s + X[i, k]*A[j, k]
  return(s)
  
}

sigma_square <- function(A, X, j, idx){
  
  # this function computes the sum of the squared terms in the
  # estimated sigma formula
  
  ss = 0
  for(i in 1:dim(X)[1]) ss = ss + (X[i, j] - sigma_sum(A, X, i, j, idx))**2
  return(ss)
  
}

sigma_est <-function(A, X, idx){
  
  # this function estimates the sigma of a given matrix
  
  c = 1/(dim(X)[1]*length(idx))
  sss = 0
  for(j in 1:length(idx)) sss = sss + sigma_square(A, X, j, idx)
  return(c*sss)
  
}

lik_est <- function(A, X, idx, sig){
  
  # this function estimates the logLikelihood
  # of an adjacency matrix given the data
  
  c = 1/(2*sig)
  a = (dim(X)[1]/2)*log(sig)
  s = 0
  for(j in 1:length(idx)) s = s + c*sigma_square(A, X, j, idx) + a
  return(-s)
  
}

# Linkage Test ------------------------------------------------------------

generate_matrix <- function(p){
  
  # this function generates a random adjacency
  # matrix given a sparsity parameter
  
  sparsity = 1/p
  A <- matrix(rbinom(p*p, 1, sparsity)*sign(runif(p*p, -1, 1)), p, p)
  A[upper.tri(A, diag=TRUE)] <- 0
  
  return(A)
  
}

generate_data <- function(A, n, p){
  
  # given an adjacency matrix this function generates
  # a data matrix with a gaussian stochastic error
  
  Sigma = solve(diag(p) - A)
  X <- matrix(rnorm(n*p), n, p) %*% t(Sigma)
  return(X)
  
}

split_data <- function(X){
  
  # This function splits the data in train and test set
  # given a split proportion
  
  n <- dim(X)[1]
  idx <- sample(1:n, n/2)
  train <- X[idx,]
  test <- X[-idx,]
  
  return(list(train=train, test=test))
  
}

generate_edge <- function(idx){
  
  # this function generates a random edge, diagonal excluded
  
  a = sample(idx, 1)
  b = sample(idx[idx != a], 1)
  return(list(a = a, b = b))
  
}

Un <- function(tr, te, idx, mu, D){
  
  # estimate the Adjacency matrix
  out_tr <- MLEdag(X=tr,D=D,tau=0.35,mu=mu,rho=1.2,trace_obj=FALSE)
  out_te <- MLEdag(X=te, D=D,tau=0.35,mu=mu,rho=1.2,trace_obj=FALSE)
  A.0 = out_tr$A.H0; A = out_te$A.H1
  # estimate sigma for constrained and unconstrained
  sig_alt = sigma_est(A, te, idx)
  sig_null = sigma_est(A.0, tr, idx)
  # compute the test statistic
  x = lik_est(A, tr, idx, sig_alt) - lik_est(A.0, tr, idx, sig_null)
  return(x)
  
}

linkage_test <- function(X, D, idx, alpha, mu){
  
  # split the data
  split = split_data(X)
  tr = split$train; te = split$test
  
  # compute the test statistics U_n and W_n
  Un_base = Un(tr, te, idx, mu, D)
  Un_swap = Un(te, tr, idx, mu, D)
  Wn = (exp(Un_base) + exp(Un_swap))/2
  
  # test decision
  lik_ratio = Un_base > -log(alpha)
  lik_ratio_swap = Wn > 1/alpha
  
  return(list(lrt=lik_ratio, crossfit=lik_ratio_swap))
}

compute_size_power <- function(M, D, A, n, alpha, mu){
  
  # this function computes the size of a linkage test
  # running M tests and computing the size as the number of rejections
  # over the number of simulations under the null hypothesis
  
  p <- dim(D)[1]
  idx <- c(1:p)
  size_lrt = rep(NA, M)
  size_crossfit = rep(NA, M)
  
  for(m in 1:M){
    X = generate_data(A, n, p)
    exit = linkage_test(X, D, idx, alpha, mu)
    size_lrt[m] = exit$lrt
    size_crossfit[m] = exit$crossfit
    progress(m, M)
  }
  
  cat('\nLRT: ', sum(size_lrt)/M, '\n')
  cat('Crossfit: ', sum(size_crossfit)/M)
}

# Pathway Test ------------------------------------------------------------

generate_path <- function(p, len){
  
  # this function generates a matrix D with a random path of a given length
  
  D = matrix(0, p, p)
  v <- sample(c(1:(p-len)), 1)
  for(i in v:(v+len-1)) D[i, i+1] = 1
  return(D)
  
}

lik_est_mod <- function(A, X, idx, i){
  
  # this is a modified version of the logLikelihood that computes the logL 
  # on a modified adjacency matrix where we put a zero for one of the edges
  # belonging to the path. Vectorizing the function we can pass directly 
  # the vector of cells of the path instead of a single index
  
  A[i] = 0 # modifying the adjacency matrix
  sig = sigma_est(A, X, idx)
  c = 1/(2*sig)
  a = (dim(X)[1]/2)*log(sig)
  s = 0
  for(j in 1:length(idx)) s = s + c*sigma_square(A, X, j, idx) + a
  return(-s)
  
}

lik_est_mod = Vectorize(lik_est_mod, vectorize.args = 'i') # vectorizing over the index

find_denominator <- function(A, D, X, idx){
  
  # this function finds the denominator of the Un test statistc
  
  path = which(D != 0)
  return(max(lik_est_mod(A, X, idx, path))) # find the max along all the logL
  
}

Un_path <- function(tr, te, idx, mu, D){
  
  # computes the test statistic for the pathway test
  
  out_tr <- MLEdag(X=tr,D=D,tau=0.35,mu=mu,rho=1.2,trace_obj=FALSE)
  out_te <- MLEdag(X=te, D=D,tau=0.35,mu=mu,rho=1.2,trace_obj=FALSE)
  A.0 = out_tr$A.H0; A = out_te$A.H1
  sig_alt = sigma_est(A, te, idx)
  denominator = find_denominator(A.0, D, tr, idx)
  x = lik_est(A, tr, idx, sig_alt) - denominator
  
  return(x)
  
}

generate_hypo <- function(p, len){
  
  # this function generates a null hypothesis for the pathway test
  # by setting in A to 0 some cell that is a nonzero in D
  
  A <- generate_matrix(p)
  D <- generate_path(p, len)
  path <- which(D != 0) # find the path cells
  
  k <- sample(c(1:(len-1)), 1) # select a random number of edges < path length
  idx <- sample(path, k) # pick k random cells to equal to zero
  A[idx] = 0
  
  return(list(D=D, A=A))
  
}

generate_alt_hypo <- function(p, len){
  
  # this function generates an alternative hypothesis for the pathway test
  # by setting in A to a certain strength every cell that is a nonzero in D
  
  A <- generate_matrix(p)
  D <- generate_path(p, len)
  path <- which(D != 0) # find the path cells
  A[path] = runif(len, -1, 1)
  
  return(list(D=D, A=A))
  
}

pathway_test <- function(X, D, idx, alpha, mu){
  
  # this function implements a universal linkage test
  
  # split the data
  split = split_data(X)
  tr = split$train; te = split$test
  
  # compute the test statistics Un and Wn
  Un_base = Un_path(tr, te, idx, mu, D)
  Un_swap = Un_path(te, tr, idx, mu, D)
  Wn = (exp(Un_base) + exp(Un_swap))/2
  
  # test decision
  lik_ratio = Un_base > -log(alpha)
  lik_ratio_swap = Wn > 1/alpha
  
  return(list(lrt=lik_ratio, crossfit=lik_ratio_swap))
  
}

compute_size_power_path <- function(M, D, A, n, alpha, mu){
  
  # this function computes the size of a linkage test
  # running M tests and computing the size as the number of rejections
  # over the number of simulations under the null hypothesis
  
  p <- dim(D)[1]
  idx <- c(1:p)
  size_lrt = rep(NA, M)
  size_crossfit = rep(NA, M)
  
  for(m in 1:M){
    X = generate_data(A, n, p)
    exit = pathway_test(X, D, idx, alpha, mu)
    size_lrt[m] = exit$lrt
    size_crossfit[m] = exit$crossfit
    progress(m, M)
  }
  
  cat('LRT: ', sum(size_lrt)/M, '\n')
  cat('Crossfit: ', sum(size_crossfit)/M)
}

# Cytometry ---------------------------------------------------------------

check_normality <- function(data, col, scale = F){
  
  qqPlot(data[,col], ylab = col)
  if(isFALSE(scale)){
    hist(data[,col], col = 'tomato', border = 'white', prob = T, xlab = col, main = 'boxcox')
    curve(dnorm(x), col = 'firebrick', lwd = 1.5, add = T)
  }
  else{
    hist(scale(data[,col]), xlim = c(-4, 4), col = 'tomato', border = 'white', prob = T, xlab = col, main = 'scale')
    curve(dnorm(x), col = 'firebrick', lwd = 1.5, add = T)
  }
}

normalizer <- function(x){
  
  out = boxcox(x~1, lambda = seq(-2,2,0.001), plotit = F)
  lam = out$x[which.max(out$y)]
  
  if(lam==0) x = log(x)
  else x = (x^lam-1)/lam
  
  return(scale(x))
  
}





