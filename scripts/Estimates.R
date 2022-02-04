# Sigma and Likelihood estimation -----------------------------------------

# Estimate sigma ----------------------------------------------------------

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

# Estimate Likelihood -----------------------------------------------------

lik_est <- function(A, X, idx, sig){
  
  # this function estimates the logLikelihood
  # of an adjacency matrix given the data
  
  c = 1/(2*sig)
  a = (dim(X)[1]/2)*log(sig)
  s = 0
  for(j in 1:length(idx)) s = s + c*sigma_square(A, X, j, idx) + a
  return(-s)
  
}