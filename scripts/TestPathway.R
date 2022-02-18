# Parameters --------------------------------------------------------------

p <- 15
len <- 3
alpha <- 0.05
mu <- 1
M <- 200 # M = {200, 1000, 2000}
n = 500 # n = {500, 1000, 1500}
idx = c(1:p)

# Size -------------------------------------------------------------------

# generating the matrix D and A compatible with the null hypothesis
hypo = generate_hypo(p, len)
D = hypo$D; A = hypo$A

# computing the size as the number of times we reject H0 knowing H0 is true.
# the test returns true if we reject H0

compute_size_power_path(M, D, A, n, alpha, mu)

# Power -------------------------------------------------------------------

# generating the matrix D and A compatible with the alternative hypothesis
hypo = generate_alt_hypo(p, len)
D = hypo$D; A = hypo$A

# computing the power as the number of times we correctly reject H0 knowing H1 is true.

compute_size_power_path(M, D, A, n, alpha, mu)
