# Collecting the data -----------------------------------------------------

import()
l = load_data()
data = rbind(l[[1]],l[[2]],l[[3]],l[[4]],l[[5]],l[[6]],l[[7]],l[[8]],l[[9]])

# Checking (pseudo) normality ------------------------------------------------------

par(mfrow = c(1, 2))
check_normality(data, 'praf')

# Pre-processing -----------------------------------------------------------

# Box-Cox

scaled = apply(data, MARGIN = 2, FUN = normalizer)

par(mfrow=c(2, 2))
check_normality(scaled, 'praf')
check_normality(data, 'praf', scale = T)

# Simulation -------------------------------------------------------------

colnames(scaled)
p = dim(scaled)[2]
idx = c(1:p)
D = matrix(0, p, p)

D[5, 7] = 1 # PIP3 --> AKT, missing (H0)
D[6, 7] = 1 # ERK --> AKT, reported (H1) 
D[9, 8] = 1 # PKC --> PKA, reported (H1) 
D[9, 11] = 1 # PKC --> PJNK, expected

M = 500
mu = 1
alpha = 0.05
lrt = rep(NA, M)
crf = rep(NA, M)

for(m in 1:M){
  
  out = linkage_test(scaled, D, idx, alpha, mu)
  lrt[m] = out$lrt # true se rifiuto H0
  crf[m] = out$crossfit
  progress(m, M)
  
}

mle = MLEdag(X=scaled,D=D,tau=0.35,mu=mu,rho=1.2,trace_obj=FALSE)
mle$pval < alpha
sum(lrt)/M
sum(crf)/M


# Path simulation --------------------------------------------------------

D = matrix(0, p, p)
# PKC --> RAF --> MEK --> ERK --> AKT
D[9, 1] = 1; D[1, 2] = 1; D[2, 6] = 1; D[6, 7] = 1

lrt = rep(NA, M)
crf = rep(NA, M)
mu = 1
M = 100
for(m in 1:M){
  
  out = pathway_test(scaled, D, idx, alpha, mu)
  lrt[m] = out$lrt # true se rifiuto H0
  crf[m] = out$crossfit
  progress(m, M)
  
}

mle = MLEdag(X=scaled,D=D,tau=0.35,mu=mu,rho=1.2,trace_obj=FALSE)
mle$pval < alpha
sum(lrt)/M
sum(crf)/M

