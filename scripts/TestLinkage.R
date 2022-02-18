# Parameters --------------------------------------------------------------

p <- 25
alpha <- 0.05 # alpha = {0.01, 0.05}
mu <- 1
M <- 1000 # M = {1000, 2000} 
n <- 500 # n = {200, 300, 500}
idx <- c(1:p)

# Size --------------------------------------------------------------------

###############################################################################
############################## Single-Edge ####################################
###############################################################################

# generating the matrix D and A compatible with the null hypothesis
D <- matrix(0, p, p); D[2,1] = 1
A <- generate_matrix(p); A[2,1] = 0 # compatible with H0

# computing the size as the number of times we reject H0 knowing H0 is true.
# the test returns true if we reject H0

compute_size_power(M, D, A, n, alpha, mu)

###############################################################################
############################ More edges (hub) #################################
###############################################################################

# generating the matrix D and A compatible with the null hypothesis
D <- matrix(0, p, p)
D[2:16,c(p-1,p)] = 1

# compatible with H0
A      <- matrix(0, p, p)     
A[, 1] <- sign( runif( p, min = -1, max = 1 ) )
A[1,1] <- 0

# computing the size as the number of times we reject H0 knowing H0 is true.

compute_size_power(1000, D, A, 1000, 0.01, mu)

# Power -------------------------------------------------------------------

###############################################################################
############################## Single-Edge ####################################
###############################################################################

# generating the matrix D and A compatible with the alternative hypothesis
D <- matrix(0, p, p); D[2,1] = 1
A <- generate_matrix(p); A[2,1] = 1 # compatible with H1

# computing the power as the number of times we correctly reject H0 knowing H1 is true.

compute_size_power(M, D, A, 200, 0.01, mu)

###############################################################################
############################ More edges (hub) #################################
###############################################################################

# generating the matrix D and A compatible with the alternative hypothesis
D <- matrix(0, p, p)
D[2:16,c(p-1,p)] = 1

# four different hypothesis compatible with H1
A <- generate_matrix(p)
#A[2:16,c(p-1,p)] = 0.1 # nearest
A[2:16,c(p-1,p)] = 0.3
#A[2:16,c(p-1,p)] = 0.5
#A[2:16,c(p-1,p)] = 1 # the most distant

# computing the power as the number of times we correctly reject H0 knowing H1 is true.

compute_size_power(M, D, A, n, alpha, mu)

