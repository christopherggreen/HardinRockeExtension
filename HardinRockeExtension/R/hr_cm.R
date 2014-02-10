hr.cm <- function(p,n,N,B=10000,mcd.alpha=0.5,logfile="hr_cm.log") {

	require(rrcov, quiet=TRUE)
	require(mvtnorm, quiet=TRUE)

	n.mcdalpha <- length(mcd.alpha)

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
			if ( (j%% 50)==0 ) cat(j," ")
			for ( a in 1:n.mcdalpha ) {
				mcd <- CovMcd(simdata[,,j],alpha=mcd.alpha[a])
				# get the covariance of the MCD subset
				mcd.cov <- cov( simdata[mcd@best,,j] )
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
	if (r > 0 ) { blockresults[[M+1]] <- blockfcn(i=M+1,b=r) }

	blockresults <- do.call("abind",c(blockresults,list(along=1)))

	dimnames(blockresults)[[2]] <- c("AVE","SQ")
	dimnames(blockresults)[[3]] <- paste("a",gsub("\\.","",mcd.alpha),sep="")

	# calculate sum of each columns
	# see Joanna Hardin's code cm.R and mcd.est.R
	mcd.ave <- colSums(blockresults[,"AVE",])
	mcd.sq  <- colSums(blockresults[,"SQ" ,])

	mcd.c   <- mcd.ave/(N*p)
	mcd.var <- (mcd.sq - (mcd.c * mcd.c * N * p))/(N * p - 1)
	mcd.m   <- 2*mcd.c*mcd.c/mcd.var

	data.frame(c=mcd.c, m=mcd.m, mcd.alpha=mcd.alpha)

}
