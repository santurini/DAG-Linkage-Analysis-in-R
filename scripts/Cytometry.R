
# Collecting the data -----------------------------------------------------

library("car")
data = rbind(data1,data2,data3,data4,data5,data6,data7,data8,data9)

apply(data9,2,mean) #zero mean col cavolo
boxplot.matrix(as.matrix(data9), use.cols = TRUE)

# Checking normality ------------------------------------------------------

par(mfcol=c(1,1))
#praf
qqnorm(data9$praf, pch = 1, frame = FALSE)
qqline(data9$praf, col = "steelblue", lwd = 2)
qqPlot(data9$praf)
#pmek
qqnorm(data9$pmek, pch = 1, frame = FALSE)
qqline(data9$pmek, col = "steelblue", lwd = 2)
qqPlot(data9$pmek)
#plcg"
qqnorm(data9$plcg, pch = 1, frame = FALSE)
qqline(data9$plcg, col = "steelblue", lwd = 2)
qqPlot(data9$plcg)
#"PIP2"
qqnorm(data9$PIP2, pch = 1, frame = FALSE)
qqline(data9$PIP2, col = "steelblue", lwd = 2)
qqPlot(data9$PIP2)
#"PIP3"
qqnorm(data9$PIP3, pch = 1, frame = FALSE)
qqline(data9$PIP3, col = "steelblue", lwd = 2)
qqPlot(data9$PIP3)
#"p44/42"
qqnorm(data9$`p44/42`, pch = 1, frame = FALSE)
qqline(data9$`p44/42`, col = "steelblue", lwd = 2)
qqPlot(data9$`p44/42`)
#"pakts473" 
qqnorm(data9$pakts473, pch = 1, frame = FALSE)
qqline(data9$pakts473, col = "steelblue", lwd = 2)
qqPlot(data9$pakts473)
#PKA"
qqnorm(data9$PKA, pch = 1, frame = FALSE)
qqline(data9$PKA, col = "steelblue", lwd = 2)
qqPlot(data9$PKA)
#PKC
qqnorm(data9$PKC, pch = 1, frame = FALSE)
qqline(data9$PKC, col = "steelblue", lwd = 2)
qqPlot(data9$PKC)
#"P38"
qqnorm(data9$P38, pch = 1, frame = FALSE)
qqline(data9$P38, col = "steelblue", lwd = 2)
qqPlot(data9$P38)
#pjnk"
qqnorm(data9$pjnk, pch = 1, frame = FALSE)
qqline(data9$pjnk, col = "steelblue", lwd = 2)
qqPlot(data9$pjnk)

# Pre-processing -----------------------------------------------------------

# Scaling
scaled = scale(data9)
round(apply(scaled,2,mean), 2)
boxplot.matrix(as.matrix(scaled), use.cols = TRUE)

# Box-Cox
library(MASS)
normalizer <- function(x){
  out = boxcox(x~1, lambda = seq(-2,2,0.001), plotit = F)
  lam = out$x[which.max(out$y)]
  if(lam==0) x = log(x)
  else x = (x^lam-1)/lam
  return(x - mean(x))
}

scaled = apply(data9, MARGIN = 2, FUN = normalizer)


# Linkage test ------------------------------------------------------

linkage_test_ad <- function(X, D, alpha){
  
  idx = c(1:dim(X)[2])
  
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


colnames(scaled)
p = dim(scaled)[2]
n = dim(scaled)
D = matrix(0, p, p)
D[5, 7] = 1

# p-value > 0.05 not refuse H0
test = seq(1, 5, 1)

MLEdag(scaled, D=D, tau=0.35,mu=1,rho=1.2,trace_obj=FALSE)
linkage_test_ad(scaled, D, alpha = 0.05)
