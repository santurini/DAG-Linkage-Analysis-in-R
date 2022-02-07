# Utilities --------------------------------------------------------------

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

# Test statistic ----------------------------------------------------------

Un_path <- function(tr, te, idx, mu, D){
  
  # computes the test statistic for the pathway test
  
  out <- MLEdag(X=tr,D=D,tau=0.35,mu=mu,rho=1.2,trace_obj=FALSE)
  A.0 = out$A.H0; A = out$A.H1
  sig_alt = sigma_est(A, te, idx)
  denominator = find_denominator(A.0, D, tr, idx)
  x = lik_est(A, tr, idx, sig_alt) - denominator
  return(x)
  
}

# Hypothesis ------------------------------------------------------------

generate_hypo_path <- function(p, n, len){
  
  # this function generates a null hypothesis for the linkage test
  # by setting in A to 0 every cell that is a nonzero in D
  
  A <- generate_matrix(p)
  D <- generate_path(p, len)
  path <- which(D != 0) # find the path cells
  
  k <- sample(c(1:(len-1)), 1) # select a random number of edges < path length
  idx <- sample(path, k) # pick k random cells to equal to zero
  A[idx] = 0
  
  X = generate_data(A, n, p)
  
  return(list(D=D, X=X))
  
}

generate_alt_hypo_path <- function(p, n, len){
  
  # this function generates a null hypothesis for the linkage test
  # by setting in A to 0 every cell that is a nonzero in D
  
  A <- generate_matrix(p)
  D <- generate_path(p, len)
  path <- which(D != 0) # find the path cells
  A[path] = runif(len, -1, 1)
  
  X = generate_data(A, n, p)
  
  return(list(D=D, X=X))
  
}

# Pathway test ------------------------------------------------------------

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

# Size --------------------------------------------------------------------

compute_size_path <- function(M, p, n, alpha, len, mu){
  
  # this function computes the size of a linkage test
  # running M tests and computing the size as the number of rejections
  # over the number of simulations under the null hypothesis
  
  idx = c(1:p)
  hypo = generate_hypo_path(p, n, len)
  D = hypo$D; X = hypo$X
  
  size_lrt = rep(NA, M)
  size_crossfit = rep(NA, M)
  
  for(m in 1:M){
    exit <- pathway_test(X, D, idx, alpha, mu)
    size_lrt[m] = exit$lrt
    size_crossfit[m] = exit$crossfit
    progress(m, M)
  }
  
  cat('LRT size: ', sum(size_lrt)/M, '\n')
  cat('Crossfit size: ', sum(size_crossfit)/M)
}

compute_size_path(100, 10, 150, 0.05, 3, 1)

# Power -------------------------------------------------------------------

compute_power_path <- function(M, p, n, alpha, len, mu){
  
  # this function computes the power of a linkage test running M tests
  # and computing the power as one minus the number of acceptances
  # over the number of simulations under the alternative hypothesis
  
  idx = c(1:p)
  hypo = generate_alt_hypo_path(p, n, len)
  D = hypo$D; X = hypo$X
  
  power_lrt = rep(NA, M)
  power_crossfit = rep(NA, M)
  
  for(m in 1:M){
    exit <- pathway_test(X, D, idx, alpha, mu)
    power_lrt[m] = exit$lrt
    power_crossfit[m] = exit$crossfit
    progress(m, M)
  }
  
  cat('LRT power: ', 1 - (M-sum(power_lrt))/M, '\n')
  cat('Crossfit power: ', 1 - (M-sum(power_crossfit))/M)
}

compute_power_path(100, 10, 150, 0.05, 3, 1)


