# Import ------------------------------------------------------------------

packages = c('clrdag', 'car', 'svMisc', 'MASS', 'readxl')
lapply(packages, require, character.only = TRUE)

# Parameters --------------------------------------------------------------

M = 200
mu = 1
alpha = 0.05

# Collecting the data -----------------------------------------------------

sheets <- c("1. cd3cd28", "2. cd3cd28icam2", "3. cd3cd28+aktinhib","4. cd3cd28+g0076",
            "5. cd3cd28+psitect","6. cd3cd28+u0126","7. cd3cd28+ly", "8. pma", "9. b2camp")

l <- vector("list", 9)
for(i in 1:9) l[[i]] <- data.frame(read_excel("cytometry-data.xlsx", sheet = sheets[i]))

data = rbind(l[[1]],l[[2]],l[[3]],l[[4]],l[[5]],l[[6]],l[[7]],l[[8]],l[[9]])

apply(data9,2,mean) #zero mean col cavolo
boxplot.matrix(as.matrix(data9), use.cols = TRUE)

# Checking (pseudo) normality ------------------------------------------------------

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

check_normality(data9, 'praf')

# Pre-processing -----------------------------------------------------------

# Box-Cox

normalizer <- function(x){
  
  out = boxcox(x~1, lambda = seq(-2,2,0.001), plotit = F)
  lam = out$x[which.max(out$y)]
  
  if(lam==0) x = log(x)
  else x = (x^lam-1)/lam
  
  return(scale(x))
  
}

scaled = apply(data9, MARGIN = 2, FUN = normalizer)

par(mfrow=c(2, 2))
check_normality(scaled, 'praf')
check_normality(data9, 'praf', scale = T)

# Simulescion -------------------------------------------------------------

colnames(scaled)
p = dim(scaled)[2]
idx = c(1:p)
D = matrix(0, p, p)
D[5, 7] = 1

mle = rep(NA, M)
nostra = rep(NA, M)

for(m in 1:M){
  #out = MLEdag(X=scaled,D=D,tau=0.35,mu=mu,rho=1.2,trace_obj=FALSE)
  out1 = linkage_test(scaled, D, idx, alpha, mu)
  nostra[m] = out1$lrt
  #mle[m] = out$pval < 0.05
  progress(m, M)
}

sum(mle)
sum(nostra)

# stamo a testa presenza di arco che non deve esserci quindi tutti false stonks
