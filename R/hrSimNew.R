hrSimNew <- 
#
# Simulation function to replicate tables on page 942 of 
# Hardin and Rocke (2005)
#
# Author: Christopher G. Green
# Date: 2011 
#
# 2011-07-19 version---changed to use new calibration library
# 2014-01-26 changed default mcd.alpha to maximum breakdown 
#   point case, which in general will be slightly larger than 0.5
#
function(p,n,N,B=10000,alpha=0.05,mcd.alpha=max.bdp.mcd.alpha(n,p)) {

  #require(rrcov, quiet=TRUE)
  #require(mvtnorm, quiet=TRUE)

  blockfcn <- function(m,b) {
    cat("Running block",m,"\n")
    # generate observations
    simdata <- array(
      rmvnorm(n*b, mean=rep(0,p), sigma=diag(rep(1,p))),
      dim=c(n,p,b)
    )

    #print(simdata[,,1])

    # compute MCD
    results <- array(NA, dim=c(b,9))
    cat("\t Iteration: ")
    hr <- hr05CutoffMvnormal(n.obs=n,p.dim=p,
      mcd.alpha=mcd.alpha,signif.alpha=alpha, method="CG")
    
    for ( j in 1:b ) {
      if ( (j%% 50)==0 ) cat(j," ")
      mcd <- covMcd(simdata[,,j],alpha=mcd.alpha)
      # remove scaling: Hardin and Rocke assumed the 
      # covariance matrix was not scaled by the 
      # small-sample correction
      cs           <- mcd$raw.cnp2[1]
      ss           <- mcd$raw.cnp2[2]

      # raw.mah is based on the MCD estimate before reweighting
      # and hence was corrected, so need to remove the small-sample
      # correction: corrected sigma estimate is cs * ss * raw mcd estimate
      # so mahalanobis distances are scaled by 1/(cs*ss)
      # multiply by ss to remove the small-sample correction

      mahdist      <- mcd$raw.mah * ss
      #mahdist2      <- mahalanobis(simdata[,,j], center=mcd$raw.center,
      #  cov=mcd$raw.cov/(ss * cs) )
      #print( sum(abs(mahdist2-mahdist)) )

      # completely uncorrected
      results[j,1] <- mean(mahdist*cs          > qchisq(1-alpha,df=p))
      # consistency corrected
      results[j,2] <- mean(mahdist             > qchisq(1-alpha,df=p))
      # consistency and small-sample
      results[j,3] <- mean(mahdist/ss          > qchisq(1-alpha,df=p))

      results[j,4] <- mean(mahdist*cs          > hr$cutoff.asy)
      results[j,5] <- mean(mahdist             > hr$cutoff.asy)
      results[j,6] <- mean(mahdist/ss          > hr$cutoff.asy)

      results[j,7] <- mean(mahdist*cs          > hr$cutoff.pred)
      results[j,8] <- mean(mahdist             > hr$cutoff.pred)
      results[j,9] <- mean(mahdist/ss          > hr$cutoff.pred)
    }
    cat("\nDone with block",m,"\n\n")
    results
  }

  B <- min(N,B)
  M <- floor(N/B)

  blockresults <- lapply(1:M, blockfcn,b=B)

  # do last block
  r <- N - M*B
  if (r > 0 ) { blockresults[[M+1]] <- blockfcn(i=M+1,b=r) }

  blockresults <- do.call("rbind",blockresults)
  dimnames(blockresults)[[2]] <- c("CHI2.RAW","CHI2.CON","CHI2.SM","CGASY.RAW",
    "CGASY.CON","CGASY.SM","CGPRED.RAW","CGPRED.CON","CGPRED.SM")
  blockresults


}
