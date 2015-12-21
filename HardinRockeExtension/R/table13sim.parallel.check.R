table13sim.parallel.check <- function(cl,p,nn,N,B=250,alpha=c(0.01,0.025,0.05),
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

    # compute cutoff values for table1 using new method, these will not change in 
    # the block
    gm50.t1 <- sapply(alpha, function(a)
      hr05CutoffMvnormal(nn,p,mcd.alpha=0.50,signif.alpha=a,method="GM14")$cutoff.pred)
    hr50.orig.t1 <- sapply(alpha, function(a)
      hr05CutoffMvnormal(nn,p,               signif.alpha=a,method="HR")$cutoff.pred)
    gmmbp.t1 <- sapply(alpha, function(a)
      hr05CutoffMvnormal(nn,p,               signif.alpha=a,method="GM14")$cutoff.pred)
    hrmbp.orig.t1 <- sapply(alpha, function(a)
      hr05CutoffMvnormal(nn,p,               signif.alpha=a,method="HR")$cutoff.pred)
    # compute cutoff values for table3 using new method, these will not change in 
    # the block
    gm50.t3 <- sapply(alpha, function(a)
      hr05CutoffMvnormal(nn,p,mcd.alpha=0.50,signif.alpha=a/nn,method="GM14")$cutoff.pred)
    hr50.orig.t3 <- sapply(alpha, function(a)
      hr05CutoffMvnormal(nn,p,               signif.alpha=a/nn,method="HR")$cutoff.pred)
    gmmbp.t3 <- sapply(alpha, function(a)
      hr05CutoffMvnormal(nn,p,               signif.alpha=a/nn,method="GM14")$cutoff.pred)
    hrmbp.orig.t3 <- sapply(alpha, function(a)
      hr05CutoffMvnormal(nn,p,               signif.alpha=a/nn,method="HR")$cutoff.pred)

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

      results <- array(NA, dim=c(b,26,length(alpha)))
      cat("\tIteration: ", file=lgf, append=TRUE)

      cat("\tp = ",p,"nn = ",nn,"b = ",b,"\n", file=lgf, append=TRUE)
      for ( j in 1:b ) {
        if ( (j %% 100)==0 ) cat("\t",j," \n", file=lgf, append=TRUE)

        # calculate average number of rejections at chi-squared(p,1-alpha)
        # level

        # mcd(mbp)
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
          results[j, 1,k] <- mean(mcdmbp@raw.mah > qpalpha.t1[k])
          # unweighted mcd with revised HR distribution
          results[j, 2,k] <- mean(mcdmbp@raw.mah > gmmbp.t1[k])
          # unweighted mcd with original HR distribution
          results[j, 3,k] <- mean(mcdmbp@raw.mah > hrmbp.orig.t1[k])
          # unweighted mcd with no small sample correction (revised)
          results[j, 4,k] <- mean(mahdist.hr    > gmmbp.t1[k])
          # unweighted mcd with no small sample correction (original)
          results[j, 5,k] <- mean(mahdist.hr    > hrmbp.orig.t1[k])
          # reweighted mcd (qchisq(0.950,p))
          results[j, 6,k] <- mean(mcdmbp@mah     > qpalpha.t1[k])

		  # table 3
          # unweighted mcd
          results[j, 7,k] <- (max(mcdmbp@raw.mah) > qpalpha.t3[k])
          # unweighted mcd with revised HR distribution
          results[j, 8,k] <- (max(mcdmbp@raw.mah) > gmmbp.t3[k])
          # unweighted mcd with original HR distribution
          results[j, 9,k] <- (max(mcdmbp@raw.mah) > hrmbp.orig.t3[k])
          # unweighted mcd with no small sample correction (revised)
          results[j,10,k] <- (max(mahdist.hr)    > gmmbp.t3[k])
          # unweighted mcd with no small sample correction (original)
          results[j,11,k] <- (max(mahdist.hr)    > hrmbp.orig.t3[k])
          # reweighted mcd (qchisq(0.950,p))
          results[j,12,k] <- (max(mcdmbp@mah)     > qpalpha.t3[k])
        }

        # Reweighted-MCD with Bonferroni adjustment to the weighting
		for ( k in alphaind ) {
          mcdmbp <- CovMcd2(simdata[,,j], 
		    reweighting=list(method="chisquare",quantile=1-(alpha[k]/nn),
			  debug=FALSE))
		  # table 3
		  results[j,13,k] <- (max(mcdmbp@mah)     > qpalpha.t3[k])
		}

        # mcd(50)
        mcd50 <- tryCatch(CovMcd(simdata[,,j],alpha=0.50), error = function(e) e)
        if ( inherits(mcd50,"simpleError") ) {
          cat("\t",conditionMessage(mcd50),"\n", file=lgf, append=TRUE)
          stop("MCDMBP was the culprit.")
        }
        # hardin and rocke used consistency factor, but only seemed to
        # use the small sample factor for n < 100
        ss  <- mcd50@raw.cnp2[2]
        # covMCD adjusts cov by multiplying by consistency and small-sample
        # factors; these therefore appear inverted in the distance metric
        # so to get rid of them we would multiply by ss
        mahdist.hr   <- if (nn < 100) mcd50@raw.mah else mcd50@raw.mah*ss
        for ( k in alphaind ) {
		  # table 1
          # unweighted mcd
          results[j,14,k] <- mean(mcd50@raw.mah > qpalpha.t1[k])
          # unweighted mcd with revised GM distribution
          results[j,15,k] <- mean(mcd50@raw.mah > gm50.t1[k])
          # unweighted mcd with original HR distribution
          results[j,16,k] <- mean(mcd50@raw.mah > hr50.orig.t1[k])
          # unweighted mcd with no small sample correction (revised)
          results[j,17,k] <- mean(mahdist.hr    > gm50.t1[k])
          # unweighted mcd with no small sample correction (original)
          results[j,18,k] <- mean(mahdist.hr    > hr50.orig.t1[k])
          # reweighted mcd (qchisq(0.950,p))
          results[j,19,k] <- mean(mcd50@mah     > qpalpha.t1[k])

		  # table 3
          # unweighted mcd
          results[j,20,k] <- (max(mcd50@raw.mah) > qpalpha.t3[k])
          # unweighted mcd with revised GM distribution
          results[j,21,k] <- (max(mcd50@raw.mah) > gm50.t3[k])
          # unweighted mcd with original HR distribution
          results[j,22,k] <- (max(mcd50@raw.mah) > hr50.orig.t3[k])
          # unweighted mcd with no small sample correction (revised)
          results[j,23,k] <- (max(mahdist.hr)    > gm50.t3[k])
          # unweighted mcd with no small sample correction (original)
          results[j,24,k] <- (max(mahdist.hr)    > hr50.orig.t3[k])
          # reweighted mcd (qchisq(0.950,p))
          results[j,25,k] <- (max(mcd50@mah)     > qpalpha.t3[k])
        }

        # Reweighted-MCD with Bonferroni adjustment to the weighting
		for ( k in alphaind ) {
          mcd50 <- CovMcd2(simdata[,,j], 
		    reweighting=list(method="chisquare",quantile=1-(alpha[k]/nn),
			  debug=FALSE))
		  # table 3
		  results[j,26,k] <- (max(mcd50@mah)     > qpalpha.t3[k])
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
     "qpalpha.t1","qpalpha.t3",
	 "gm50.t1","hr50.orig.t1","gmmbp.t1","hrmbp.orig.t1",
	 "gm50.t3","hr50.orig.t3","gmmbp.t3","hrmbp.orig.t3",
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
     "MCDMBP.RAW.T1","MCDMBP.RAWGM.T1","MCDMBP.RAWHR.T1",
	 "MCDMBP.RAWNOSSGM.T1","MCDMBP.RAWNOSSHR.T1","RMCDMBP.T1",
     "MCDMBP.RAW.T3","MCDMBP.RAWGM.T3","MCDMBP.RAWHR.T3",
	 "MCDMBP.RAWNOSSGM.T3","MCDMBP.RAWNOSSHR.T3","RMCDMBP.T3",
	 "RMCDMBP.CH.T3",
     "MCD50.RAW.T1","MCD50.RAWGM.T1","MCD50.RAWHR.T1",
	 "MCD50.RAWNOSSGM.T1","MCD50.RAWNOSSHR.T1","RMCD50.T1",
     "MCD50.RAW.T3","MCD50.RAWGM.T3","MCD50.RAWHR.T3",
	 "MCD50.RAWNOSSGM.T3","MCD50.RAWNOSSHR.T3","RMCD50.T3",
	 "RMCD50.CH.T3"
   )
 
   dimnames(blockresults)[[3]] <- paste("alpha",alpha,sep="") 
   blockresults
}
