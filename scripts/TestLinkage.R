# Utilities ---------------------------------------------------------------

generate_matrix <- function(p){
  
  # this function generates a random adjacency
  # matrix given a sparsity parameter
  
  sparsity = 2/p
  A <- matrix(rbinom(p*p, 1, sparsity)*sign(runif(p*p, -1, 1))*runif(p*p, 0.7, 1), p, p)
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

split_data <- function(X, prop = 0.8){
  
  # This function splits the data in train and test set
  # given a split proportion
  
  split <- sample(c(rep(0, prop * nrow(X)), rep(1, (1-prop) * nrow(X))))
  train <- X[split == 0, ]
  test <- X[split == 1, ]
  return(list(train=train, test=test))
  
}

# Hypothesis --------------------------------------------------------------

generate_edge <- function(idx){
  
  # this function generates a random edge, diagonal excluded
  
  a = sample(idx, 1)
  b = sample(idx[idx != a], 1)
  return(list(a = a, b = b))
  
}

generate_hypo <- function(p, idx, num_edges){
  
  # this function generates a null hypothesis for the linkage test
  # by setting in A to 0 every cell that is a nonzero in D
  
  A <- generate_matrix(p)
  D <- matrix(0, p, p)
  for(i in 1:num_edges){
    edge = generate_edge(idx)
    D[edge$a, edge$b] = 1
    A[edge$a, edge$b] = 0
  }
  return(list(D=D, A=A))
}

generate_alt_hypo <- function(p, idx, num_edges, type = 'random'){
  
  # this function generates an alternative hypothesis for the linkage test
  # the parameter type chooses the type of alternative hypothesis
  # type: random (default), distant, close (w.r.t the null hypothesis)
  
  A <- generate_matrix(p)
  D <- matrix(0, p, p)
  
  if( type == 'random'){
    edges = matrix(NA, num_edges, 2)
    for(i in 1:num_edges){
      edge = generate_edge(idx)
      D[edge$a, edge$b] = 1
      A[edge$a, edge$b] = 0
      edges[i, 1] = edge$a; edges[i, 2] = edge$b
    }
    # select one of the random edges and put it equal to 1
    row_idx = sample(c(1:num_edges), 1)
    A[edges[row_idx, 1], edges[row_idx, 2]] = 1
  }
  
  if( type == 'close'){
    for(i in 1:num_edges){
      edge = generate_edge(idx)
      D[edge$a, edge$b] = 1
      A[edge$a, edge$b] = 0.05  # to make it close to H0 we set a very small value
    }
  } 
  
  if( type == 'distant'){
    for(i in 1:num_edges){
      edge = generate_edge(idx)
      D[edge$a, edge$b] = 1
      A[edge$a, edge$b] = sign(runif(1, -1, 1))  # the farest value to 0 is {-1, 1}
    }
  }
  return(list(D=D, A=A))
}

# Test statistics ---------------------------------------------------------

Un <- function(train, test, idx, D){
  # estimate the Adjacency matrix
  out <- MLEdag(X=train,D=D,tau=0.35,mu=1,rho=1.2,trace_obj=FALSE)
  A.0 = out$A.H0; A = out$A.H1
  # estimate sigma for constrained and unconstrained
  sig_alt = sigma_est(A, test, idx)
  sig_null = sigma_est(A.0, train, idx)
  # compute the test statistic
  x = lik_est(A, train, idx, sig_alt) - lik_est(A.0, train, idx, sig_null)
  return(x)
} 

# Universal test function ----------------------------------------------------------

linkage_test <- function(p, n, alpha, num_edges, alternative = NULL){
  
  # this function implements a universal linkage test
  
  idx = c(1:p)
  
  # generate hypothesis (null or alternative)
  if(is.null(alternative)){
    hypo = generate_hypo(p, idx, num_edges)
    D = hypo$D
    A = hypo$A
  }
  else{
    hypo = generate_alt_hypo(p, idx, num_edges, alternative)
    D = hypo$D
    A = hypo$A
  }
  
  # generate the data
  X = generate_data(A, n, p)
  
  # split the data
  split_train_test = split_data(X, 0.5)
  train = split_train_test$train; test = split_train_test$test
  
  # compute the test statistics Un and Wn
  Un_base = Un(train, test, idx, D); Un_swap = Un(test, train, idx, D)
  Wn = (exp(Un_base) + exp(Un_swap))/2
  
  # test decision
  lik_ratio = Un_base > -log(alpha)
  lik_ratio_swap = Wn > 1/alpha
  
  return(list(lrt=lik_ratio, crossfit=lik_ratio_swap))
}

# Size --------------------------------------------------------------------

compute_size <- function(M, p, n, alpha, num_edges){
  
  # this function computes the size of a linkage test
  # running M tests and computing the size as the number of rejections
  # over the number of simulations under the null hypothesis
  
  size_lrt = rep(NA, M)
  size_crossfit = rep(NA, M)
  for(m in 1:M){
    test <- linkage_test(p, n, alpha, num_edges)
    size_lrt[m] = test$lrt
    size_crossfit[m] = test$crossfit
  }
  
  cat('LRT size: ', sum(size_lrt)/M, '\n')
  cat('Crossfit size: ', sum(size_crossfit)/M)
}

compute_size(50, 10, 150, alpha, 1)

# Power -------------------------------------------------------------------

compute_power <- function(M, p, n, alpha, num_edges, alternative = 'random'){
  
  # this function computes the power of a linkage test running M tests
  # and computing the power as one minus the number of acceptances
  # over the number of simulations under the alternative hypothesis
  
  power_lrt = rep(NA, M)
  power_crossfit = rep(NA, M)
  for(m in 1:M){
    test <- linkage_test(p, n, alpha, num_edges, alternative)
    power_lrt[m] = test$lrt
    power_crossfit[m] = test$crossfit
  }
  
  cat('LRT power: ', 1 - (M-sum(power_lrt))/M, '\n')
  cat('Crossfit power: ', 1 - (M-sum(power_crossfit))/M)
}

compute_power(50, 10, 100, alpha, 1)

