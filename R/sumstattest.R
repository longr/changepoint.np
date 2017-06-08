
nonparametric.ed.sumstat = function(data,K=10){ # This now takes into account the integral transformation
  ##USE K points in integral
  n <- length(data)
  if(K>n) K=n
  Q <- matrix(0,K,n+1)
  x=sort(data)
  yK= -1 + (2*(1:K)/K-1/K)
  c=-log(2*n-1)
  pK=(1+exp(c*yK))^-1
  for (i in 1:K){
    j=as.integer((n-1)*pK[i] + 1)
    Q[i,-1] <- cumsum(data<x[j])+0.5*cumsum(data==x[j])
  }
  return(Q)
}

sumstatTest = function(sumstat,nquantiles = 4){
  storage.mode(sumstat) = 'double'
  storage.mode(nquantiles)='integer'
  sumstatout=rep(0,nquantiles) # sets up null vector for changepoint answer
  storage.mode(sumstatout)='double'
  n=dim(sumstat)[2]-1
  lastchangelike=rep(0,n)
  storage.mode(lastchangelike) = 'double'
  lastchangecpts=rep(0,n)
  storage.mode(lastchangecpts) = 'integer'
  numchangecpts=rep(0,n)
  storage.mode(numchangecpts) = 'integer'
  sumlastchangelike = 0;
  cptsout = rep(0,n)
  storage.mode(cptsout) = 'integer'
  pen = 20
  cost_func = "nonparametric.ed";
  answer=.C('sumstattest',cost_func, sumstat,nquantiles, as.integer(n), as.double(pen), cptsout, lastchangelike, lastchangecpts, numchangecpts)

  return(list(answer[[2]],answer[[6]]))
}
