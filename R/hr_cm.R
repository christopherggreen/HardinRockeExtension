hr.cm <- 
#
# Simulation function to extend Hardin and Rocke (2005)
# to mcd fractions other than the maximum breakdown point
# case
#
# Author: Christopher G. Green
# Date: 2011 
#
# 2014-02-10 begin finalizing for publication 
# 
# TODO: ensure randomness across blocks in the simulated
#       data; right now each block could potentially use
#       the same samples?
#
function(p,n,N,B=10000,
  mcd.alpha=max.bdp.mcd.alpha,logfile="hr_cm.log") 
{

  n.mcdalpha <- length(mcd.alpha)

  # function to run on each block
  blockfcn <- function(m,b) {
    cat("Running block",m,"\n",file=logfile,append=TRUE)
    # generate observations
    simdata <- array(
      rmvnorm(n*b, mean=rep(0,p), sigma=diag(rep(1,p))),
      dim=c(n,p,b)
    )

    #print(simdata[,,1])

    # compute MCD
    results <- array(NA, dim=c(b,2,n.mcdalpha))
    cat("\t Iteration: ",file=logfile,append=TRUE)
    
    for ( j in 1:b ) {
      if ( (j %% 50)==0 ) cat(j," ")
      for ( a in 1:n.mcdalpha ) {
        mcd <- CovMcd(simdata[,,j],alpha=mcd.alpha[a])
        # get the covariance of the MCD subset
		# 2014-03-01 when alpha == 1, mcd@best is 
		# null b/c CovMcd uses the classical covariance
        mcd.cov <- if ( mcd.alpha[a]==1.0 ) {
			#cov( simdata[,,j] ) 
			# covariance already computed by CovMcd
			# use raw.cov to avoid any potential reweighting
			# see code for covMcd in robustbase
			mcd.cov@raw.cov
		} else {
			cov( simdata[mcd@best,,j] )
		}
        results[j,1,a] <- sum(diag(mcd.cov))
        results[j,2,a] <- sum(diag(mcd.cov)^2)
      }
    }
    cat("\nDone with block",m,"\n\n",file=logfile,append=TRUE)
    results
  }

  B <- min(N,B)
  M <- floor(N/B)

  blockresults <- lapply(1:M, blockfcn,b=B)

  # do last block
  r <- N - M*B
  if (r > 0 ) { blockresults[[M+1]] <- blockfcn(m=M+1,b=r) }

  blockresults <- do.call("abind",c(blockresults,list(along=1)))

  dimnames(blockresults)[[2]] <- c("AVE","SQ")
  dimnames(blockresults)[[3]] <- paste("a",gsub("\\.","",mcd.alpha),sep="")

  # calculate sum of each columns
  # see Joanna Hardin's code cm.R and mcd.est.R
  mcd.ave <- colSums(blockresults[,"AVE",])
  mcd.sq  <- colSums(blockresults[,"SQ" ,])

  # estimate of c is average of all diagonal elements
  mcd.c   <- mcd.ave/(N*p)
  # (sum of squares - n * square of mean)/(n - 1)
  mcd.var <- (mcd.sq - (mcd.c * mcd.c * N * p))/(N * p - 1)
  # estimate of m follows from formula: 2 div by (coef of var)^2
  # coef of var is sqrt of var div by c
  # var of diagonal elements is sum of (s_jj - c)^2, div by
  # Np-1; rearrange numerator
  mcd.m   <- 2*mcd.c*mcd.c/mcd.var

  data.frame(c=mcd.c, m=mcd.m, mcd.alpha=mcd.alpha)

}
