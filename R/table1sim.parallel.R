table1sim.parallel <- function(cl,p,nn,N,B=10000,alpha=c(0.01,0.025,0.05),
    cutoff.method="GM14", lgf="") {

    # need to push these to cluster init function
    #require(abind,      quietly=TRUE)
    #require(rrcov,      quietly=TRUE)
    #require(mvtnorm,    quietly=TRUE)
    #require(timeSeries, quietly=TRUE)
    #require(cggRcalibration, quietly=TRUE)
    #require(CerioliOutlierDetection, quietly=TRUE)

    maxtries <- 50

    blockfcn <- function(m,b) {
      cat("\nRunning block",m,"\n", file=lgf, append=TRUE)
      # generate observations
      simdata <- tryCatch(array(
        mvtnorm::rmvnorm(nn*b, mean=rep(0,p), sigma=diag(rep(1,p))),
        dim=c(nn,p,b)
      ), error = function(e) e)
      if ( inherits(simdata, "simpleError") ) {
        cat(conditionMessage(simdata),"\n",file=logfile,append=TRUE)
        stop("Failed to generate random data set.")
      }
      #  if ( any(!is.finite(simdata) ) ) {
      #  stop("Missing or infinite values in random number array.")
      #  }

      results <- array(NA, dim=c(b,16,length(alpha)))
      cat("\t Iteration: ", file=lgf, append=TRUE)

      # compute chi-square cutoff value which will not change in block
      qpalpha <- qchisq(1-alpha,df=p)
      # compute cutoff values using new method, these will not change in 
      # the block
      hr50 <- sapply(alpha, function(a)
        hr05CutoffMvnormal(nn,p,               signif.alpha=a,method=cutoff.method)$cutoff.pred)
      hr75 <- sapply(alpha, function(a)
        hr05CutoffMvnormal(nn,p,mcd.alpha=0.75,signif.alpha=a,method=cutoff.method)$cutoff.pred)
      hr95 <- sapply(alpha, function(a)
        hr05CutoffMvnormal(nn,p,mcd.alpha=0.95,signif.alpha=a,method=cutoff.method)$cutoff.pred)
      alphaind <- seq(along=alpha)

      cat("p = ",p,"nn = ",nn,"b = ",b,"\n", file=lgf, append=TRUE)
      for ( j in 1:b ) {
        if ( (j%% 50)==0 ) cat(j," ", file=lgf, append=TRUE)

        # calculate average number of rejections at chi-squared(p,1-alpha)
        # level

        # reweighted ogk with beta = 0.9 
        ogk <- tryCatch(CovOgk(simdata[,,j], beta = 0.9), error=function(e) e)
        if ( inherits(ogk,"simpleError") ) {
          cat(conditionMessage(ogk),"\n", file=logfile, append=TRUE)
          stop("ogk was the culprit.")
        }
        for ( k in alphaind ) {
          # unweighted ogk distances
          results[j,1,k] <- mean(ogk@raw.mah       > qpalpha[k])
          # reweighted ogk distances
          results[j,2,k] <- mean(getDistance(ogk)  > qpalpha[k])
        }

        # bisquare S-estimator
        sest <- tryCatch(CovSest(simdata[,,j], method="bisquare"), error=function(e) e)
        if ( inherits(sest,"simpleError") ) {
          cat(conditionMessage(sest),"\n", file=logfile, append=TRUE)
          stop("BS S-est was the culprit.")
        }
        for ( k in alphaind ) {
          results[j,3,k] <- mean(getDistance(sest) > qpalpha[k])
        }

        # rocke S-estimator (won't work for n = 50 and p = 16)
		# using asymptotic rejectin probability of 0.10
        sest <- tryCatch(CovSest(simdata[,,j], method="rocke", arp=0.10), error=function(e) e)
        if ( inherits(sest,"simpleError") ) {
          cat(conditionMessage(sest),"\n",file=lgf, append=TRUE)
          stop("Rocke sest was the culprit.")
        }
        for ( k in alphaind ) {
            results[j,4,k] <- mean(getDistance(sest) > qpalpha[k])
        }

        # mcd(0.50)
        mcd50 <- tryCatch(CovMcd(simdata[,,j]), error = function(e) e)
        if ( inherits(mcd50,"simpleError") ) {
          cat(conditionMessage(mcd50),"\n", file=logfile, append=TRUE)
          stop("MCD50 was the culprit.")
        }
        # hardin and rocke used consistency factor, but only seemed to
        # use the small sample factor for n < 100
        ss  <- mcd50@raw.cnp2[2]
        # covMCD adjusts cov by multiplying by consistency and small-sample
        # factors; these therefore appear inverted in the distance metric
        # so to get rid of them we would multiply by ss
        mahdist.hr   <- if (nn < 100) mcd50@raw.mah else mcd50@raw.mah*ss
        for ( k in alphaind ) {
          # unweighted mcd
          results[j,5,k] <- mean(mcd50@raw.mah > qpalpha[k])
          # unweighted mcd with HR distribution
          results[j,6,k] <- mean(mcd50@raw.mah > hr50[k])
          # unweighted mcd with no small sample correction
          results[j,7,k] <- mean(mahdist.hr    > hr50[k])
          # reweighted mcd (qchisq(0.950,p))
          results[j,8,k] <- mean(mcd50@mah     > qpalpha[k])
        }

        # mcd(0.75)
        mcd75 <- tryCatch(CovMcd(simdata[,,j],alpha=0.75), error=function(e) e)
        if ( inherits(mcd75,"simpleError") ) {
          cat(conditionMessage(mcd75),"\n", file=logfile, append=TRUE)
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
          # unweighted mcd
          results[j, 9,k] <- mean(mcd75@raw.mah > qpalpha[k])
          # unweighted mcd with HR distribution
          results[j,10,k] <- mean(mcd75@raw.mah > hr75[k])
          # unweighted mcd with no small sample correction
          results[j,11,k] <- mean(mahdist.hr    > hr75[k])
          # reweighted mcd (qchisq(0.975,p))
          results[j,12,k] <- mean(mcd75@mah     > qpalpha[k])
        }

        # mcd(0.95)
        mcd95 <- tryCatch(CovMcd(simdata[,,j],alpha=0.95), error=function(e) e)
        if ( inherits(mcd95,"simpleError") ) {
          cat(conditionMessage(mcd95),"\n", file=logfile, append=TRUE)
          stop("MCD95 was the culprit.")
        }
        ss  <- mcd95@raw.cnp2[2]
        mahdist.hr   <- if (nn < 100) mcd95@raw.mah else mcd95@raw.mah*ss
        for ( k in alphaind ) {
          # unweighted mcd
          results[j,13,k] <- mean(mcd95@raw.mah > qpalpha[k])
          # unweighted mcd with HR distribution
          results[j,14,k] <- mean(mcd95@raw.mah > hr95[k])
          # unweighted mcd with no small sample correction
          results[j,15,k] <- mean(mahdist.hr    > hr95[k])
          # reweighted mcd (qchisq(0.975,p))
          results[j,16,k] <- mean(mcd95@mah     > qpalpha[k])
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

   clusterExport(cl = cl, c("M","B","blockfcn","p","nn","N","alpha"),
     envir=sys.frame(pf) )

   blockresults <- parLapply(cl=cl, X=worklist, function(i) blockfcn(m=i[1],b=i[2]))

   # 2011/10/12 occasionally get an error from blockfcn (e.g., a bad data set)
   #blockresults <- lapply(1:M, function(i,bb) try(blockfcn(i,bb), silent=FALSE), b=B)
   # rerun any failed blocks
   #for ( i in 1:M ) {
   #  tries <- 0
   #  while ( (tries < maxtries) && (inherits(blockresults[[i]], "try-error")) ) {
   #    cat("\nRetrying block",i,"\n", file=lgf, append=TRUE)
   #    blockresults[[i]] <- try(blockfcn(i,B), silent=FALSE)
   #    tries <- tries + 1
   #  }
   #  if ( inherits(blockresults[[i]], "try-error") ) 
   #      stop(paste("Failed to calculate block",i,"\n"))  
   #
   #}
 
   #     # do last block
   #     r <- N - M*B
   #     if ( r > 0 ) { 
   #       blockresults[[M+1]] <- try(blockfcn(i=M+1,b=r), silent=FALSE)
   #       tries <- 0
   #    while ( (tries < maxtries) && ( inherits(blockresults[[M+1]], "try-error") ) ) {
   #      cat("\nRetrying block",M+1,"\n", file=lgf, append=TRUE)
   #      blockresults[[M+1]] <- try(blockfcn(i=M+1,b=r), silent=FALSE)
   #      tries <- tries + 1
   #    }
   #    if ( inherits(blockresults[[M+1]], "try-error") ) 
   #      stop(paste("Failed to calculate block",M+1,"\n"))  
   #  }
 
   blockresults <- do.call("abind",c(blockresults,list(along=1)))
   dimnames(blockresults)[[2]] <- c(
     "OGK","ROGK","SEST.BS","SEST.RK",
     "MCD50.RAW","MCD50.HRRAW","MCD50.HRADJ","RMCD50",
     "MCD75.RAW","MCD75.HRRAW","MCD75.HRADJ","RMCD75",
     "MCD95.RAW","MCD95.HRRAW","MCD95.HRADJ","RMCD95"
   )
 
   dimnames(blockresults)[[3]] <- paste("alpha",alpha,sep="") 
   blockresults
}
