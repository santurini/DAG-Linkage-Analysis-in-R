values = c(0.00001,0.0001,0.001,0.01,0.1, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80,0.85, 0.90, 0.95, 1.0, 1.1,1.3,1.5,1.7,1.9,2.1,3,5,10)

cr = rep(1, length(values))

lr = c(0.75,0.74,0.73,0.75,0.745,0.82,0.77,0.775,0.715,0.75,0.755,0.765,0.8,0.75,0.755,0.72,0.71,0.695,0.76,0.77,0.78,0.795,0.75,0.715,0.735,0.775,0.73,0.74,0.745,0.765,0.725,0.74)

par(mfcol=c(1,2))
plot(values,lr, xlim=c(0,1),ylim=c(0.7,0.9),type="b",col="blue",xlab="mu values",ylab="Split LRT")
plot(values,lr, xlim=c(0,5),ylim=c(0.7,0.9),type="b",col="blue",xlab="mu values",ylab="Split LRT")

library(ggplot2)
library(cowplot)

# Basic scatter plot
g1=ggplot(data=NULL, aes(x=values, y=lr)) + geom_jitter(size=2.5, shape=18,color="blue") +xlim(0,1) +  xlab("mu values")+  ylab("Split LRt") 

g2=ggplot(data=NULL, aes(x=values, y=lr)) + geom_jitter(size=2, shape=18,color="blue") +xlim(0,10) +  xlab("mu values")+  ylab("Split LRt") 
plot_grid(g1, g2, labels = "AUTO")
