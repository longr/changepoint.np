#data <- c(rnorm(100,2,1), rnorm(100,5,1), rnorm(100,3,1))
#penalty="MBIC"
#pen.value=0
#method="NPPELT"
#test.stat="empirical_distribution"
#class=TRUE
#minseglen=1
#pen.value = 10
#source("sumstattest.R")
#dyn.load("~/Dropbox/changepoint.np/src/sumstattest.so")
# source("NPPELT.R")
# dyn.load("~/Dropbox/changepoint.np/src/NPPELT.so")
# set.seed(1)
# data <- rnorm(100)
# #data
# #sumstat <- nonparametric.ed.sumstat(data, K = 4)
# #sumstat
# #sumstatTest(sumstat,nquantiles = 4)
#
# nonparametric.ed.sumstat = function(data,K=nquantiles){ # This now takes into account the integral transformation
#   ##USE K points in integral
#   n <- length(data)
#   if(K>n) K=n
#   Q <- matrix(0,K,n+1)
#   x=sort(data)
#   yK= -1 + (2*(1:K)/K-1/K)
#   c=-log(2*n-1)
#   pK=(1+exp(c*yK))^-1
#   for (i in 1:K){
#     j=as.integer((n-1)*pK[i] + 1)
#     Q[i,-1] <- cumsum(data<x[j])+0.5*cumsum(data==x[j])
#   }
#   return(Q)
# }
# sumstat <- nonparametric.ed.sumstat(data, K = 4)
# NPPELT(sumstat,pen=2*log(length(data)), cost_func = "nonparametric.ed", minseglen = 1, nquantile = 4)
