table13sim.parallel <- function(cl,p,nn,N,B=250,alpha=c(0.01,0.025,0.05),
    lgf="",mlgf="") {

    # need to push these to cluster init function
    #require(abind,      quietly=TRUE)
    #require(rrcov,      quietly=TRUE)
    #require(mvtnorm,    quietly=TRUE)
    #require(timeSeries, quietly=TRUE)
    #require(cggRcalibration, quietly=TRUE)
    #require(CerioliOutlierDetection, quietly=TRUE)

    maxtries <- 50

    # compute chi-square cutoff value for table1 which will not change in block
    qpalpha.t1 <- qchisq(1-alpha,df=p)
    # compute chi-square cutoff value for table3 which will not change in block
    qpalpha.t3 <- qchisq(1-(alpha/nn),df=p)
	# 2015-01-18 alternate cutoff value given that we are testing a maximum of
	# iid random variables
	qpalpha.t3.alt <- qchisq(1-(alpha^(1./nn)),df=p)

    # compute cutoff values for table1 using new method, these will not change in 
    # the block
    hrmbp.t1 <- sapply(alpha, function(a)
      hr05CutoffMvnormal(nn,p,               signif.alpha=a,method="GM14")$cutoff.pred)
    hrmbp.orig.t1 <- sapply(alpha, function(a)
      hr05CutoffMvnormal(nn,p,               signif.alpha=a,method="HR")$cutoff.pred)
    hr75.t1 <- sapply(alpha, function(a)
      hr05CutoffMvnormal(nn,p,mcd.alpha=0.75,signif.alpha=a,method="GM14")$cutoff.pred)
    hr95.t1 <- sapply(alpha, function(a)
      hr05CutoffMvnormal(nn,p,mcd.alpha=0.95,signif.alpha=a,method="GM14")$cutoff.pred)
    # compute cutoff values for table3 using new method, these will not change in 
    # the block
    hrmbp.t3 <- sapply(alpha, function(a)
      hr05CutoffMvnormal(nn,p,               signif.alpha=a/nn,method="GM14")$cutoff.pred)
    hrmbp.orig.t3 <- sapply(alpha, function(a)
      hr05CutoffMvnormal(nn,p,               signif.alpha=a/nn,method="HR")$cutoff.pred)
    hr75.t3 <- sapply(alpha, function(a)
      hr05CutoffMvnormal(nn,p,mcd.alpha=0.75,signif.alpha=a/nn,method="GM14")$cutoff.pred)
    hr95.t3 <- sapply(alpha, function(a)
      hr05CutoffMvnormal(nn,p,mcd.alpha=0.95,signif.alpha=a/nn,method="GM14")$cutoff.pred)
	# alternate cutoffs for simultaneous case
    hrmbp.t3.alt <- sapply(alpha, function(a)
      hr05CutoffMvnormal(nn,p,               signif.alpha=(a^(1./nn)),
	  	method="GM14")$cutoff.pred)
    hrmbp.orig.t3.alt <- sapply(alpha, function(a)
      hr05CutoffMvnormal(nn,p,               signif.alpha=(a^(1./nn)),
	  	method="HR")$cutoff.pred)
    hr75.t3.alt <- sapply(alpha, function(a)
      hr05CutoffMvnormal(nn,p,mcd.alpha=0.75,signif.alpha=(a^(1./nn)),
	  	method="GM14")$cutoff.pred)
    hr95.t3.alt <- sapply(alpha, function(a)
      hr05CutoffMvnormal(nn,p,mcd.alpha=0.95,signif.alpha=(a^(1./nn)),
	  	method="GM14")$cutoff.pred)

    alphaind <- seq(along=alpha)

    blockfcn <- function(m,b) {
      cat("\nRunning block",m,"\n", file=lgf, append=TRUE)
      # generate observations
      simdata <- tryCatch(array(
        mvtnorm::rmvnorm(nn*b, mean=rep(0,p), sigma=diag(rep(1,p))),
        dim=c(nn,p,b)
      ), error = function(e) e)
      if ( inherits(simdata, "simpleError") ) {
        cat(conditionMessage(simdata),"\n",file=lgf,append=TRUE)
        stop("Failed to generate random data set.")
      }
      #  if ( any(!is.finite(simdata) ) ) {
      #  stop("Missing or infinite values in random number array.")
      #  }

      results <- array(NA, dim=c(b,74,length(alpha)))
      cat("\tIteration: ", file=lgf, append=TRUE)

      cat("\tp = ",p,"nn = ",nn,"b = ",b,"\n", file=lgf, append=TRUE)
      for ( j in 1:b ) {
        if ( (j %% 100)==0 ) cat("\t",j," \n", file=lgf, append=TRUE)

        # calculate average number of rejections at chi-squared(p,1-alpha)
        # level

        # reweighted ogk with beta = 0.9 
        ogk <- tryCatch(CovOgk(simdata[,,j], beta = 0.9), error=function(e) e)
        if ( inherits(ogk,"simpleError") ) {
          cat("\t",conditionMessage(ogk),"\n", file=lgf, append=TRUE)
          stop("ogk was the culprit.")
        }
        for ( k in alphaind ) {
		  # table 1
          # unweighted ogk distances
          results[j, 1,k] <- mean(ogk@raw.mah       > qpalpha.t1[k])
          # reweighted ogk distances
          results[j, 2,k] <- mean(getDistance(ogk)  > qpalpha.t1[k])

          # table 3
          # unweighted ogk distances
          results[j,17,k] <- (max(ogk@raw.mah)       > qpalpha.t3[k])
          # reweighted ogk distances
          results[j,18,k] <- (max(getDistance(ogk))  > qpalpha.t3[k])

		  # new versions with max quantile
          results[j,49,k] <- (max(ogk@raw.mah)       > qpalpha.t3.alt[k])
		  results[j,50,k] <- (max(getDistance(ogk))  > qpalpha.t3.alt[k])
        }
		# reweighted ogk with beta = 1 - (0.1/nn) for table 3
        ogk <- tryCatch(CovOgk(simdata[,,j], beta = (1 - (0.1/nn))), error=function(e) e)
        if ( inherits(ogk,"simpleError") ) {
          cat("\t",conditionMessage(ogk),"\n", file=lgf, append=TRUE)
          stop("ogk was the culprit.")
        }
        for ( k in alphaind ) {
		  # table 3
          # reweighted ogk distances
          results[j,19,k] <- (max(getDistance(ogk))  > qpalpha.t3[k])
		  results[j,51,k] <- (max(getDistance(ogk))  > qpalpha.t3.alt[k])
        }

        # bisquare S-estimator
        sest <- tryCatch(CovSest(simdata[,,j], method="bisquare"), error=function(e) e)
        if ( inherits(sest,"simpleError") ) {
          cat("\t",conditionMessage(sest),"\n", file=lgf, append=TRUE)
          stop("BS S-est was the culprit.")
        }
        for ( k in alphaind ) {
		  # table 1
          results[j, 3,k] <- mean(getDistance(sest) > qpalpha.t1[k])

		  # table 3
          results[j,20,k] <- (max(getDistance(sest)) > qpalpha.t3[k])
          results[j,52,k] <- (max(getDistance(sest)) > qpalpha.t3.alt[k])
        }

        # rocke S-estimator (won't work for n = 50 and p >= 16)
		# using asymptotic rejectin probability of 0.05 per MMY chapter 6

		# 2014-10-19 the occasional error from Rocke seems to be 
		# related to a bad random sample within the routine.
		# can move the random number generator along, and try again
		# this often seems to get us an answer
        sest <- tryCatch(CovSest(simdata[,,j], method="rocke", arp=0.05), 
		  error=function(e) e)
        if ( inherits(sest,"simpleError") ) {
          #cat("\tError was: ",conditionMessage(sest),"\n",file=lgf, append=TRUE)
		  tries <- 1
		  while ( (inherits(sest, "simpleError")) && (tries < maxtries) ) {
            cat("\tError was: ",conditionMessage(sest),"\n",file=lgf, append=TRUE)
			cat("\tRetrying iteration ",j," in this block.\n",file=lgf,append=TRUE)
			runif(100)
            sest <- tryCatch(CovSest(simdata[,,j], method="rocke", arp=0.05), 
		      error=function(e) e)
		    tries <- tries + 1
		  }
		  if ( inherits(sest,"simpleError") ) {
			cat("\tFailed ", tries, "times. Aborting.\n",file=lgf,append=TRUE)
		    stop(paste("Failed to calculate iteration ",j,": Rocke sest was the culprit."))
		  } 
        }
        for ( k in alphaind ) {
		    # table 1
            results[j, 4,k] <- mean(getDistance(sest) > qpalpha.t1[k])

			# table 3
            results[j,21,k] <- (max(getDistance(sest)) > qpalpha.t3[k])
            results[j,53,k] <- (max(getDistance(sest)) > qpalpha.t3.alt[k])
        }

        # mcd(MBP)
        mcdmbp <- tryCatch(CovMcd(simdata[,,j]), error = function(e) e)
        if ( inherits(mcdmbp,"simpleError") ) {
          cat("\t",conditionMessage(mcdmbp),"\n", file=lgf, append=TRUE)
          stop("MCDMBP was the culprit.")
        }
        # hardin and rocke used consistency factor, but only seemed to
        # use the small sample factor for n < 100
        ss  <- mcdmbp@raw.cnp2[2]
        # covMCD adjusts cov by multiplying by consistency and small-sample
        # factors; these therefore appear inverted in the distance metric
        # so to get rid of them we would multiply by ss
        mahdist.hr   <- if (nn < 100) mcdmbp@raw.mah else mcdmbp@raw.mah*ss
        for ( k in alphaind ) {
		  # table 1
          # unweighted mcd
          results[j, 5,k] <- mean(mcdmbp@raw.mah > qpalpha.t1[k])
          # unweighted mcd with revised HR distribution
          results[j, 6,k] <- mean(mcdmbp@raw.mah > hrmbp.t1[k])
          # unweighted mcd with original HR distribution
          results[j,38,k] <- mean(mcdmbp@raw.mah > hrmbp.orig.t1[k])
          # unweighted mcd with no small sample correction
          results[j, 7,k] <- mean(mahdist.hr    > hrmbp.t1[k])
          # unweighted mcd with no small sample correction
          results[j,39,k] <- mean(mahdist.hr    > hrmbp.orig.t1[k])
          # reweighted mcd (qchisq(0.950,p))
          results[j, 8,k] <- mean(mcdmbp@mah     > qpalpha.t1[k])

		  # table 3
          # unweighted mcd
          results[j,22,k] <- (max(mcdmbp@raw.mah) > qpalpha.t3[k])
          # unweighted mcd with revised HR distribution
          results[j,23,k] <- (max(mcdmbp@raw.mah) > hrmbp.t3[k])
          # unweighted mcd with revised HR distribution
          results[j,40,k] <- (max(mcdmbp@raw.mah) > hrmbp.orig.t3[k])
          # unweighted mcd with no small sample correction
          results[j,24,k] <- (max(mahdist.hr)    > hrmbp.t3[k])
          # Hardin-Rocke's original method
          results[j,25,k] <- (max(mahdist.hr)    > hrmbp.orig.t3[k])
          # reweighted mcd (qchisq(0.950,p))
          results[j,26,k] <- (max(mcdmbp@mah)     > qpalpha.t3[k])

          # unweighted mcd
          results[j,54,k] <- (max(mcdmbp@raw.mah) > qpalpha.t3.alt[k])
          # unweighted mcd with revised HR distribution
          results[j,55,k] <- (max(mcdmbp@raw.mah) > hrmbp.t3.alt[k])
          # unweighted mcd with revised HR distribution
          results[j,56,k] <- (max(mcdmbp@raw.mah) > hrmbp.orig.t3.alt[k])
          # unweighted mcd with no small sample correction
          results[j,57,k] <- (max(mahdist.hr)    > hrmbp.t3.alt[k])
          # Hardin-Rocke's original method
          results[j,58,k] <- (max(mahdist.hr)    > hrmbp.orig.t3.alt[k])
          # reweighted mcd (qchisq(0.950,p))
          results[j,59,k] <- (max(mcdmbp@mah)     > qpalpha.t3.alt[k])
        }

        # Reweighted-MCD with Bonferroni adjustment to the weighting
		for ( k in alphaind ) {
          mcdmbp <- CovMcd2(simdata[,,j], 
		    reweighting=list(method="chisquare",quantile=1-(alpha[k]/nn),
			  debug=FALSE))
		  # table 3
		  results[j,27,k] <- (max(mcdmbp@mah)     > qpalpha.t3[k])
		  results[j,60,k] <- (max(mcdmbp@mah)     > qpalpha.t3.alt[k])
		}

        # mcd(0.75)
        mcd75 <- tryCatch(CovMcd(simdata[,,j],alpha=0.75), error=function(e) e)
        if ( inherits(mcd75,"simpleError") ) {
          cat("\t",conditionMessage(mcd75),"\n", file=lgf, append=TRUE)
          stop("MCD75 was the culprit.")
        }
        # hardin and rocke used consistency factor, but only seemed to
        # use the small sample factor for n < 100
        ss  <- mcd75@raw.cnp2[2]
        # covMCD adjusts cov by multiplying by consistency and small-sample
        # factors; these therefore appear inverted in the distance metric
        # so to get rid of them we would multiply by ss
        mahdist.hr   <- if (nn < 100) mcd75@raw.mah else mcd75@raw.mah*ss
        for ( k in alphaind ) {
		  # table 1
          # unweighted mcd
          results[j, 9,k] <- mean(mcd75@raw.mah > qpalpha.t1[k])
          # unweighted mcd with HR distribution
          results[j,10,k] <- mean(mcd75@raw.mah > hr75.t1[k])
          # unweighted mcd with original HR distribution
          results[j,41,k] <- mean(mcd75@raw.mah > hrmbp.orig.t1[k])
          # unweighted mcd with no small sample correction
          results[j,11,k] <- mean(mahdist.hr    > hr75.t1[k])
          # unweighted mcd with no small sample correction
          results[j,42,k] <- mean(mahdist.hr    > hrmbp.orig.t1[k])
          # reweighted mcd (qchisq(0.975,p))
          results[j,12,k] <- mean(mcd75@mah     > qpalpha.t1[k])

		  # table 3
          # unweighted mcd
          results[j,28,k] <- (max(mcd75@raw.mah) > qpalpha.t3[k])
          # unweighted mcd with HR distribution
          results[j,29,k] <- (max(mcd75@raw.mah) > hr75.t3[k])
          # unweighted mcd with revised HR distribution
          results[j,43,k] <- (max(mcd75@raw.mah) > hrmbp.orig.t3[k])
          # unweighted mcd with no small sample correction
          results[j,30,k] <- (max(mahdist.hr)    > hr75.t3[k])
          # unweighted mcd with no small sample correction
          results[j,44,k] <- (max(mahdist.hr)    > hrmbp.orig.t3[k])
          # reweighted mcd (qchisq(0.975,p))
          results[j,31,k] <- (max(mcd75@mah)     > qpalpha.t3[k])

          # unweighted mcd
          results[j,61,k] <- (max(mcd75@raw.mah) > qpalpha.t3.alt[k])
          # unweighted mcd with HR distribution
          results[j,62,k] <- (max(mcd75@raw.mah) > hr75.t3.alt[k])
          # unweighted mcd with revised HR distribution
          results[j,63,k] <- (max(mcd75@raw.mah) > hrmbp.orig.t3.alt[k])
          # unweighted mcd with no small sample correction
          results[j,64,k] <- (max(mahdist.hr)    > hr75.t3.alt[k])
          # unweighted mcd with no small sample correction
          results[j,65,k] <- (max(mahdist.hr)    > hrmbp.orig.t3.alt[k])
          # reweighted mcd (qchisq(0.975,p))
          results[j,66,k] <- (max(mcd75@mah)     > qpalpha.t3.alt[k])
        }
        # Reweighted-MCD with Bonferroni adjustment to the weighting
		for ( k in alphaind ) {
          mcd75 <- CovMcd2(simdata[,,j], alpha=0.75, 
		    reweighting=list(method="chisquare",quantile=1-(alpha[k]/nn),
			  debug=FALSE))
		  # table 3
		  results[j,32,k] <- (max(mcd75@mah)     > qpalpha.t3[k])
		  results[j,67,k] <- (max(mcd75@mah)     > qpalpha.t3.alt[k])
		}

        # mcd(0.95)
        mcd95 <- tryCatch(CovMcd(simdata[,,j],alpha=0.95), error=function(e) e)
        if ( inherits(mcd95,"simpleError") ) {
          cat("\t",conditionMessage(mcd95),"\n", file=lgf, append=TRUE)
          stop("MCD95 was the culprit.")
        }
        ss  <- mcd95@raw.cnp2[2]
        mahdist.hr   <- if (nn < 100) mcd95@raw.mah else mcd95@raw.mah*ss
        for ( k in alphaind ) {
		  # table 1
          # unweighted mcd
          results[j,13,k] <- mean(mcd95@raw.mah > qpalpha.t1[k])
          # unweighted mcd with HR distribution
          results[j,14,k] <- mean(mcd95@raw.mah > hr95.t1[k])
          # unweighted mcd with original HR distribution
          results[j,45,k] <- mean(mcd95@raw.mah > hrmbp.orig.t1[k])
          # unweighted mcd with no small sample correction
          results[j,15,k] <- mean(mahdist.hr    > hr95.t1[k])
          # unweighted mcd with no small sample correction
          results[j,46,k] <- mean(mahdist.hr    > hrmbp.orig.t1[k])
          # reweighted mcd (qchisq(0.975,p))
          results[j,16,k] <- mean(mcd95@mah     > qpalpha.t1[k])

		  # table 3
          # unweighted mcd
          results[j,33,k] <- (max(mcd95@raw.mah) > qpalpha.t3[k])
          # unweighted mcd with HR distribution
          results[j,34,k] <- (max(mcd95@raw.mah) > hr95.t3[k])
          # unweighted mcd with revised HR distribution
          results[j,47,k] <- (max(mcd95@raw.mah) > hrmbp.orig.t3[k])
          # unweighted mcd with no small sample correction
          results[j,35,k] <- (max(mahdist.hr)    > hr95.t3[k])
          # unweighted mcd with no small sample correction
          results[j,48,k] <- (max(mahdist.hr)    > hrmbp.orig.t3[k])
          # reweighted mcd (qchisq(0.975,p))
          results[j,36,k] <- (max(mcd95@mah)     > qpalpha.t3[k])

          # unweighted mcd
          results[j,68,k] <- (max(mcd95@raw.mah) > qpalpha.t3.alt[k])
          # unweighted mcd with HR distribution
          results[j,69,k] <- (max(mcd95@raw.mah) > hr95.t3.alt[k])
          # unweighted mcd with revised HR distribution
          results[j,70,k] <- (max(mcd95@raw.mah) > hrmbp.orig.t3.alt[k])
          # unweighted mcd with no small sample correction
          results[j,71,k] <- (max(mahdist.hr)    > hr95.t3.alt[k])
          # unweighted mcd with no small sample correction
          results[j,72,k] <- (max(mahdist.hr)    > hrmbp.orig.t3.alt[k])
          # reweighted mcd (qchisq(0.975,p))
          results[j,73,k] <- (max(mcd95@mah)     > qpalpha.t3.alt[k])
        }
        # Reweighted-MCD with Bonferroni adjustment to the weighting
		for ( k in alphaind ) {
          mcd95 <- CovMcd2(simdata[,,j], alpha=0.95, 
		    reweighting=list(method="chisquare",quantile=1-(alpha[k]/nn),
			  debug=FALSE))
		  # table 3
		  results[j,37,k] <- (max(mcd95@mah)     > qpalpha.t3[k])
		  results[j,74,k] <- (max(mcd95@mah)     > qpalpha.t3.alt[k])
		}

     } # looping over the block
     cat("\nDone with block",m,"\n\n", file=lgf, append=TRUE)
     results
   } # end blockfcn

   B <- min(N,B)
   M <- floor(N/B)
   
   worklist <- matrix( c(1:M,rep(B,M)), ncol=2 )

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

   clusterExport(cl = cl, c("M","B","blockfcn","p","nn","N","alpha",
     "qpalpha.t1","qpalpha.t3","hrmbp.t1","hrmbp.orig.t1",
	 "hr75.t1","hr95.t1","hrmbp.t3",
	 "hrmbp.orig.t3","hr75.t3","hr95.t3",
	 "qpalpha.t3.alt","hrmbp.t3.alt","hrmbp.orig.t3.alt",
	 "hr75.t3.alt","hr95.t3.alt",
	 "alphaind","maxtries"),
     envir=sys.frame(pf) )

   blockresults <- parLapply(cl=cl, X=worklist, 
     function(i) try(blockfcn(m=i[1],b=i[2]),silent=FALSE) )
   
   # 2011/10/12 occasionally get an error from blockfcn (e.g., a bad data set)
   # rerun any failed blocks
   failed <- which(sapply( blockresults, function(x) inherits(x,"try-error") ))
   if ( length(failed) > 0 ) {

   		.rerun <- function(j,mxt) {
			res <- try(blockfcn(m=j[1],b=j[2]),silent=FALSE) 
			tries <- 1
			while ( (inherits(res, "try-error")) && (tries < maxtries) ) {
				cat("\nRetrying block",j[1],"\n",file=lgf,append=TRUE)
				res <- try(blockfcn(m=j[1],b=j[2]), silent=FALSE)
				tries <- tries + 1
			}
			if ( inherits(res,"try-error") ) {
				cat("Failed ", tries, "times. Aborting.\n",file=lgf,append=TRUE)
				stop(paste("Failed to calculate block",j[1],"\n"))
			} else {
				res
			}
		} # .rerun

   		blockresults[failed] <- parLapply(cl=cl,X=worklist[failed], 
			.rerun, mxt=maxtries )

#   		for ( i in failed ) {
#			tries <- 0
#			while( (inherits(blockresults[[i]], "try-error")) && (tries < maxtries) ) {
#				cat("\nRetrying block",i,"\n",file=mlgf,append=TRUE)
#				blockresults[[i]] <- try(blockfcn(m=worklist[1,i],b=worklist[2,i]), 
#				  silent=FALSE)
#				tries <- tries + 1
#			}
#			if ( inherits(blockresults[[i]],"try-error") ) {
#				cat("Failed ", tries, "times. Aborting.\n",file=mlgf,append=TRUE)
#				stop(paste("Failed to calculate block",i,"\n"))
#			}
#		}
   } # retry block

   # join together individual run results, assign dimnames
   blockresults <- do.call("abind",c(blockresults,list(along=1)))
   dimnames(blockresults)[[2]] <- c(
     # table 1 stuff
     "OGK.T1","ROGK.T1","SEST.BS.T1","SEST.RK.T1",
     "MCDMBP.RAW.T1","MCDMBP.GMRAW.T1","MCDMBP.GMADJ.T1","RMCDMBP.T1",
     "MCD75.RAW.T1","MCD75.GMRAW.T1","MCD75.GMADJ.T1","RMCD75.T1",
     "MCD95.RAW.T1","MCD95.GMRAW.T1","MCD95.GMADJ.T1","RMCD95.T1",
	 # table 3 stuff
     "OGK.T3","ROGK.T3","ROGK.CH.T3",
	 "SEST.BS.T3","SEST.RK.T3",
     "MCDMBP.RAW.T3","MCDMBP.GMRAW.T3","MCDMBP.GMADJ.T3","MCDMBP.HRADJ.T3",
	   "RMCDMBP.T3","RMCDMBP.CH.T3",
     "MCD75.RAW.T3","MCD75.GMRAW.T3","MCD75.GMADJ.T3","RMCD75.T3","RMCD75.CH.T3",
     "MCD95.RAW.T3","MCD95.GMRAW.T3","MCD95.GMADJ.T3","RMCD95.T3","RMCD95.CH.T3",
	 # backfilled stuff
	 "MCDMBP.HRRAW.T1","MCDMBP.HRADJ.T1","MCDMBP.HRRAW.T3",
	 "MCD75.HRRAW.T1","MCD75.HRADJ.T1","MCD75.HRRAW.T3","MCD75.HRADJ.T3",
	 "MCD95.HRRAW.T1","MCD95.HRADJ.T1","MCD95.HRRAW.T3","MCD95.HRADJ.T3",
	 # maximum stuff
     "OGK.T3.ALT","ROGK.T3.ALT","ROGK.CH.T3.ALT",
	 "SEST.BS.T3.ALT","SEST.RK.T3.ALT",
     "MCDMBP.RAW.T3.ALT","MCDMBP.GMRAW.T3.ALT","MCDMBP.HRRAW.T3.ALT",
	 	"MCDMBP.GMADJ.T3.ALT","MCDMBP.HRADJ.T3.ALT",
		"RMCDMBP.T3.ALT","RMCDMBP.CH.T3.ALT",
     "MCD75.RAW.T3.ALT","MCD75.GMRAW.T3.ALT","MCD75.HRRAW.T3.ALT",
	 	"MCD75.GMADJ.T3.ALT","MCD75.HRADJ.T3.ALT",
		"RMCD75.T3.ALT","RMCD75.CH.T3.ALT",
     "MCD95.RAW.T3.ALT","MCD95.GMRAW.T3.ALT","MCD95.HRRAW.T3.ALT",
	 	"MCD95.GMADJ.T3.ALT","MCD95.HRADJ.T3.ALT",
		"RMCD95.T3.ALT","RMCD95.CH.T3.ALT"
   )
 
   dimnames(blockresults)[[3]] <- paste("alpha",alpha,sep="") 
   blockresults
}
