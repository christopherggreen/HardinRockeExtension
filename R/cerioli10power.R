Cerioli2010PowerShift <- 
#
# Simulation function to calculate power of Cerioli (2010)
# IRMCD methodology with and without Green and Martin (2017)
# improvements to Hardin and Rocke (2005)
#
# This calculates power to detect a shift in location
#
# Author: Christopher G. Green
# Date: 2017 
#
#
function(#
	cl, # cluster object
	p,  # dimension aka nu in Cerioli (2010)
	n,  # sample size
	N,  # number of simulations
	B=10000, # block size for each simulation
	gamma=0.01, # nominal size/false positive rate for intersection test
	mcd.alpha=max.bdp.mcd.alpha(n,p), # alpha to use for MCD
	tau = 0.05, # contamination rate tau
    #contam.dist=c("rmvnorm","rmvt"), # contaminate by mv normal or mv t
	#contam.df=NULL, # degrees of freedom for mv t contamination
	lambda = seq(0.2,4,0.2), # vector of location translation parameter for contamination
	lgf="" # log file 
) {

  blockfcn <- function(m,b) {
    cat("Running block",m,"\n",file=lgf,append=TRUE)

    # generate b blocks of n observations
	# n (1 - \tau) come from N(0, I)
	# n \tau come from N(\lambda, I)
	#
	# start with a clean data set.

    simdata <- array(
      mvtnorm::rmvnorm(n*b, mean=rep(0,p), sigma=diag(rep(1,p)), 
	  	method="chol"), dim=c(n,p,b)
    )

    print(simdata[,,1])

	# contaminate
    # contaminate rows in each slice with probability tau
    ind <- matrix(rbinom(n=n*b, size=1, prob=tau),nrow=n,ncol=b)
	for ( block in seq_int(1,b) ) {
		rws <- which(ind[,block]==1)
		ncontam <- sum(ind)
		# if ind contains all zeros, do nothing
		if ( ncontam > 0 ) {
			for ( lmb in lambda ) {
				simdata[rws,,block] <- mvtnorm::rmvnorm( ncontam, 
					mean=lambda*rep(1, p),
					sigma=diag(rep(1,p)), method="chol") 
	# calculate FSRMCD, IRMCD with and without HR improvement
	# identify points as outliers, then calculate power as
	# the number of contaminated observations (rws in each slice)
	# that were identified as outliers, divided by the overall
	# number of contaminated observations.
			}
    	}
  	} 

	results <- rep(1,5)
    results
  }

  if ( tau <= 0 || tau >= 0.5 ) 
  	stop("Tau must be non-zero and less than 0.5.")

  B <- min(N,B)
  M <- floor(N/B)

  worklist <- matrix(c(1:M,rep(B,M)),ncol=2) 

  #blockresults <- lapply(1:M, blockfcn,b=B)
  # do last block
  r <- N - M*B
  #if (r > 0) { blockresults[[M+1]] <- blockfcn(i=M+1,b=r) }
  if ( r > 0 ) worklist <- rbind( worklist, c(M+1,r) )
  worklist <- data.frame(t(worklist))
  #print(worklist)

  # parallel version
  # install needed variables on the cluster
  # force evaluation of this now to get the current frame number
  pf <- sys.nframe()
  
  clusterExport(cl = cl, c("M","B","blockfcn","p","n","N","mcd.alpha","gamma","tau","lambda","lgf"),
  envir=sys.frame(pf) )

  blockresults <- parLapply(cl=cl, X=worklist, function(i) blockfcn(m=i[1],b=i[2]))


  blockresults <- do.call("rbind",blockresults)
  #dimnames(blockresults)[[2]] <- c("CHI2.RAW","CHI2.CON","CHI2.SM","CGASY.RAW",
  #  "CGASY.CON","CGASY.SM","CGPRED.RAW","CGPRED.CON","CGPRED.SM")
  blockresults


}

