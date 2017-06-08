############# Wrapper for chngpnt
### internal function to call the robust fpop algorithm using
### some standard parametrization of the loss function
### for now this include the L1 loss (Laplace), Huber loss
### and outlier loss (Note that to be coherent with changepoint the loss
### corresponds to the test.stat parameter)

fpop_intern <- structure(function(x,test.stat="Outlier",pen.value,lthreshold=NA,class=TRUE){
### x, A vector of double : the signal to be segmented
### test.stat, The assumed test statisticis or distribution. It can be either
### L1, Huber or Outlier
### lthreshold, threshold parameters for the Huber and Outlier test.stat.
### By default for Huber we take 1.345
### By default for Outlier we take 3

    
################################################################################
## provide an lthreshold value if lthreshold value if it is not provided
################################################################################
if(is.na(lthreshold)){
  lthreshold <- switch(test.stat,
                       L1 = NA,
                       Huber   = 1.345,
                       Outlier = 3)
}

################################################################################
## call core fpop function
################################################################################
tmp.res <- switch(test.stat,
           ### L1
           L1 = Rob_seg(x, lambda=pen.value, lthreshold=0, lslope=-1),

           ### Huber
           Huber   = Rob_seg(x, lambda=pen.value, lthreshold=lthreshold,
                             lslope=-2*lthreshold),

           ### Outlier
           Outlier = Rob_seg(x, lambda=pen.value, lthreshold=lthreshold)
         )

################################################################################
## add some missing argument for changepoints
################################################################################

tmp.res$method <- "FPOP"
tmp.res$test.stat <- test.stat
tmp.res$minseglen <- 1
## (Note for the outlier loss the min length depends on the penalty and threshold)

################################################################################
## format the object for changepoint
################################################################################

##
link.param <- matrix(ncol=3, byrow=T,
data=c("signal",     "data",         "keep",
  "n",          "n",            "discard",
  "lambda",     "pen.value",    "keep",
  "lthreshold", "lthreshold",   "keep",
  "rthreshold", "rthreshold",   "keep",
  "lslope",     "lslope",       "keep",
  "rslope",     "rslope",       "keep",
  "min",        "min",          "discard",
  "max",        "max",          "discard",
  "path",       "path",         "discard",
  "cost",       "cost",         "keep",
  "mean",       "mean.path",    "discard",
  "t.est",      "cpts",         "keep",
  "smt",        "smt.signal",   "keep",
  "out",        "out.position", "keep",
  "K",          "nb.seg",       "keep",
  "method",     "method",       "keep",
  "test.stat",  "test.stat",    "keep",
  "minseglen",  "minseglen",    "keep"
  ))

toKeep <- which(link.param[, 3] == "keep")

ind=match(link.param[, 1], names(tmp.res))[!is.na(match(link.param[, 1], names(tmp.res)))]
names(tmp.res)[ind] <- link.param[match(names(tmp.res),link.param[, 1]), 2]

if(class==TRUE){
  return(class_input(x, cpttype="nonparametric", method="FPOP", test.stat=test.stat, penalty='Manual',
                     pen.value=pen.value, minseglen=1, out=tmp.res))
}
else{
  return(tmp.res)
}

}, ex=function(){
  x <- c(rnorm(100), rnorm(100)+2)
  std.dev <- mad(diff(x)/sqrt(2))
  x_ <- x/std.dev
  lambda = log(length(x))
  res.l1 <- fpop_intern(x_,  "L1", pen.value=lambda)
  res.Hu <- fpop_intern(x_,  "Huber", pen.value=1.4*lambda, lthreshold=1.345)
  res.Ou <- fpop_intern(x_,  "Outlier", pen.value=2*lambda, lthreshold=3)
  plot(x_, pch=20)
  matlines(data.frame(res.l1$smt.signal, res.Hu$smt.signal, res.Ou$smt.signal), lty=2, lwd=2)
}
)

####################### Wrapper for C

Rob_seg <- function
### Function calling the fpop algorithm with a loss function of the form
### \gamma(X_i, \mu) = (X-i-\mu)^2 if \mu \in [X_i -lthrs, X_i+rthres]
### \gamma(X_i, \mu) = lslope \mu + la0 if \mu \in ]mini, X_i -lthrs]
### \gamma(X_i, \mu) = rslope \mu + ra0 if \mu \in [X_i +lthrs, maxi]
(x,
### A vector of double : the signal to be segmented
lambda,
### Value of the penalty
lthreshold,
### Value of the left threshold
rthreshold=lthreshold,
### Value of the left slope
lslope=0,
### Value of the left threshold
rslope=-lslope,
### Value of the left threshold
mini=min(x),
### Min value for the mean parameter of the segment
maxi=max(x)
### Max value for the mean parameter of the segment
){
  n <- length(x)
  A <- .C("rob_fpop_RtoC", signal=as.double(x), n=as.integer(n),
		lambda=as.double(lambda),
		lthreshold=as.double(lthreshold), rthreshold=as.double(rthreshold),
		lslope=as.double(lslope), rslope=as.double(rslope),
		min=as.double(mini), max=as.double(maxi),
		path=integer(n), cost=double(n) , mean=double(n)
	, PACKAGE="changepoint.np")
    A$t.est <- getPath(A$path, n)
    A$smt <- rep(A$mean[A$t.est], diff(c(0, A$t.est)))
    ## does not always make sense in particular if L1...
    A$out <- ( (x - A$smt) > A$rthreshold) | ( (x - A$smt) < -A$lthreshold )

    A$K <- length(A$t.est)
    #A$J.est <- getCostRob_Seg(A)
    return(A);
### return a list with a vector t.est containing the position of the change-points
}

getPath <- function
### This function is used by the Rob_seg function to recover the best segmentation from 1:n from the C output
(path,
 ### the path vector of the "Rob_seg" function
 i
 ### the last position to consider in the path vector
){
  chaine <- integer(1)
  chaine[1] <- length(path)
  j <- 2
  while(chaine[j-1] > 0){
    chaine[j] <- path[chaine[j-1]]
    j=j+1
  }
  return(rev(chaine)[-1])
  ### return a vector with the best change-points w.r.t. to L2 to go from point 1 to i
}




######################## Standard.

Rob_seg.std <- structure(function
### main function to use fpop for L1, Huber and biweight (outlier) losses
(
x,
### A vector of double : the signal to be segmented
loss="L1",
### loss function (L1, Huber and outlier)
lambda,
### penalty value
lthreshold
### for L1 (none), for Huber typically 1.345 if sd=1, for Outlier typically 3 if sd=1
){

if(loss=="L1"){
  Rob_seg(x, lambda=lambda, lthreshold=0, lslope=-1)-> res
}
if("Huber" == loss){
  Rob_seg(x, lambda=lambda, lthreshold=lthreshold, lslope=-2*lthreshold) -> res
}

if("Outlier"== loss){
  Rob_seg(x, lambda=lambda, lthreshold=lthreshold) -> res
}

return(res)
}, ex=function(){
  x <- c(rnorm(100), rnorm(100)+2)
  std.dev <- mad(diff(x)/sqrt(2))
  x_ <- x/std.dev
  lambda = log(length(x))
  res.l1 <- Rob_seg.std(x_,  "L1", lambda=lambda)
  res.Hu <- Rob_seg.std(x_,  "Huber", lambda=1.4*lambda, lthreshold=1.345)
  res.Ou <- Rob_seg.std(x_,  "Outlier", lambda=2*lambda, lthreshold=3)
  plot(x_, pch=20)
  matlines(data.frame(res.l1$smt, res.Hu$smt, res.Ou$smt), lty=2, lwd=2)
}

)
