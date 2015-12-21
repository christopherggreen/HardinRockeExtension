#### This is originally from the R package
####
####  rrcov : Scalable Robust Estimators with High Breakdown Point
####
#### by Valentin Todorov

##  I would like to thank Peter Rousseeuw and Katrien van Driessen for
##  providing the initial code of this function.

## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, a copy is available at
## http://www.r-project.org/Licenses/

## MM: the way it's set up, *must* be kept in sync with rrcov.control()'s
## defaults --> ./rrcov.control.R  :
## cg 2015-12-10 this is from robustbase version 0.92-5.
covMcd2 <- function(x,
           cor = FALSE,
	       raw.only = FALSE,
           alpha = control$alpha,
           nsamp = control$nsamp,
           nmini = control$nmini, kmini = control$kmini,
           scalefn=control$scalefn, maxcsteps=control$maxcsteps,
           initHsets = NULL, save.hsets = FALSE, names = TRUE,
           seed  = control$seed,
           tolSolve = control$tolSolve, # had 1e-10 hardwired {now 1e-14 default}
           trace = control$trace,
           use.correction = control$use.correction,
           wgtFUN = control$wgtFUN,
           control = rrcov.control(),
		   # cgg 2015-12-10 added this argument
           reweighting = list(method="chisquare",quantile=0.975,debug=FALSE))
{
    logdet.Lrg <- 50 ## <-- FIXME add to  rrcov.control() and then use that
    ##   Analyze and validate the input parameters ...
    if(length(seed) > 0) {
	if(length(seed) < 3 || seed[1L] < 100)
	    stop("invalid 'seed'. Must be compatible with .Random.seed !")
        if(exists(".Random.seed", envir=.GlobalEnv, inherits=FALSE))  {
            seed.keep <- get(".Random.seed", envir=.GlobalEnv, inherits=FALSE)
            on.exit(assign(".Random.seed", seed.keep, envir=.GlobalEnv))
        }
        assign(".Random.seed", seed, envir=.GlobalEnv)
    }

    ## For back compatibility, as some new args did not exist pre 2013-04,
    ## and callers of covMcd() may use a "too small"  'control' list:
    defCtrl <- if(missing(control)) control else rrcov.control()
	# cgg 2015-12-10 add robustbase:::
    if(missing(wgtFUN)) robustbase:::getDefCtrl("wgtFUN", defCtrl)
    if(is.null (nmini)) robustbase:::getDefCtrl("nmini", defCtrl)

    ##   vt::03.02.2006 - added options "best" and "exact" for nsamp
    ##   nsamp will be further analized in the wrapper .fastmcd()
    if(is.numeric(nsamp) && nsamp <= 0)
        stop("Invalid number of trials nsamp = ",nsamp, "!")

    if(is.data.frame(x))
	x <- data.matrix(x, rownames.force=FALSE)
    else if (!is.matrix(x))
        x <- matrix(x, length(x), 1,
                    dimnames = if(names) list(names(x), deparse(substitute(x))))

    if(!names) dimnames(x) <- NULL # (speedup)
    ## drop all rows with missing values (!!) :
    ok <- is.finite(x %*% rep.int(1, ncol(x)))
    x <- x[ok, , drop = FALSE]
    if(!length(dx <- dim(x)))
        stop("All observations have missing values!")
    n <- dx[1]; p <- dx[2]
    if(names) dimn <- dimnames(x)
    ## h(alpha) , the size of the subsamples
    h <- h.alpha.n(alpha, n, p)
    if(n <= p + 1) # ==> floor((n+p+1)/2) > n - 1  -- not Ok
        stop(if (n <= p) # absolute barrier!
             "n <= p -- you can't be serious!"
        else "n == p+1  is too small sample size for MCD")
    ## else
    if(n < 2 * p) { ## p+1 < n < 2p
        warning("n < 2 * p, i.e., possibly too small sample size")
        ## was stop("Need at least 2*(number of variables) observations ")
    }
    ##     jmin <- (n + p + 1) %/% 2
    ##     if(alpha < 1/2) ## FIXME? shouldn't we rather test   'alpha < jmin/n' ?
    ##  stop("The MCD must cover at least", jmin, "observations")
    ## MM: I think this should be sufficient;
    ##     we should even omit the (n < 2p) warning
    if(h > n)
        stop("Sample size n  <  h(alpha; n,p) := size of \"good\" subsample")
    else if(2*h < n)
	warning("subsample size	 h < n/2  may be too small")

    if(is.character(wgtFUN)) {
	if(is.function(mkWfun <- .wgtFUN.covMcd[[wgtFUN]]))
            wgtFUN <- mkWfun(p=p, n=n, control)
    }
    if(!is.function(wgtFUN))
	stop(gettextf("'wgtFUN' must be a function or one of the strings %s.",
		      pasteK(paste0('"',names(.wgtFUN.covMcd),'"'))), domain=NA)

	## cgg 2010-01-15
	## modified this to allow for different quantile values and for an empirical
	## quantile
	## empirical quantile of mahalanobis distances cannot be computed ahead of time
	## so I've replaced the value quantiel with a function that gets the value
	##
	## cgg 2015-12-10
	## reincorporated this into the latest version of covMcd
	## in the future I will probably convert this to a wgtFUN
	## which should remove the need to modify covMcd

	reweighting$method <- match.arg(reweighting$method, c("chisquare","empirical"))
	if ( is.null(reweighting$debug) ) reweighting$debug <- FALSE
	quantiel <- if ( reweighting$method == "chisquare" ) {
	  function() qchisq(reweighting$quantile, p)
	} else {
	  # when called, mah will be bound to the current value of mah
	  function() quantile(mah, reweighting$quantile)  
	}
	## end cgg addition

    ## vt::03.02.2006 - raw.cnp2 and cnp2 are vectors of size 2 and  will
    ##   contain the correction factors (concistency and finite sample)
    ##   for the raw and reweighted estimates respectively. Set them
    ##   initially to 1.  If use.correction is false (not the default),
    ##   the finite sample correction factor will not be used
    ##   (neither for the raw estimates nor for the reweighted ones)
    raw.cnp2 <- cnp2 <- c(1,1)

    ans <- list(call = match.call(), nsamp = nsamp,
                method = sprintf("MCD(alpha=%g ==> h=%d)", alpha, h))

    if(h == n) { ## <==> alpha ~= 1 : Just compute the classical estimates --------
        mcd <- cov(x) #MM: was  cov.wt(x)$cov
        loc <- as.vector(colMeans(x))
        obj <- determinant(mcd, logarithm = TRUE)$modulus[1]
        if ( -obj/p > logdet.Lrg ) {
            ans$cov <- mcd
	    if(names) dimnames(ans$cov) <- list(dimn[[2]], dimn[[2]])
            if (cor)
                ans$cor <- cov2cor(ans$cov)
            ans$center <- loc
            if(names && length(dimn[[2]]))
                names(ans$center) <- dimn[[2]]
            ans$n.obs <- n
            ans$singularity <- list(kind = "classical")
            ## cgg some extra debugging info
			if(trace) cat("classical estimate is singular\n")
			##
            weights <- 1
        }
        else {
            mah <- mahalanobis(x, loc, mcd, tol = tolSolve)
            ## VT:: 01.09.2004 - bug in alpha=1
			# cgg 2015-12-10 changed quantiel from value to function call
			# replacing wgtFUN call for now
            # weights <- wgtFUN(mah) # 0/1
            weights <- as.numeric(mah < quantiel()) # 0/1
            sum.w <- sum(weights)
            ans <- c(ans, cov.wt(x, wt = weights, cor = cor))
            ## cov.wt() -> list("cov", "center", "n.obs", ["wt", "cor"])
            ## Consistency factor for reweighted MCD -- ok for default wgtFUN only: FIXME
            if(sum.w != n) {
			    # cgg add robustbase:: to access internal function
                cnp2[1] <- robustbase::.MCDcons(p, sum.w/n)
                ans$cov <- ans$cov * cnp2[1]
            }
            obj <- determinant(mcd, logarithm = TRUE)$modulus[1]
            if( -obj/p > logdet.Lrg ) {
                ans$singularity <- list(kind = "reweighted.MCD")
				# cgg 2015-12-10 add more debugging info
                if(trace) cat("reweighted MCD is singular\n")
            }
            else {
                mah <- mahalanobis(x, ans$center, ans$cov, tol = tolSolve)
			    # cgg changed quantiel from value to function call
                # weights <- wgtFUN(mah) # 0/1
                weights <- as.numeric(mah < quantiel()) # 0/1
            }
        }

        ans$alpha <- alpha
        ans$quan <- h
        ans$raw.cov <- mcd
        ans$raw.center <- loc
        if(names && !is.null(nms <- dimn[[2]])) {
            names(ans$raw.center) <- nms
            dimnames(ans$raw.cov) <- list(nms,nms)
        }
        ans$crit <- obj # was exp(obj); but log-scale is "robust" against under/overflow
        ans$method <- paste(ans$method,
                            "\nalpha = 1: The minimum covariance determinant estimates based on",
                            n, "observations \nare equal to the classical estimates.")
        ans$mcd.wt <- rep.int(NA, length(ok))
        ans$mcd.wt[ok] <- weights
        if(names && length(dimn[[1]]))
            names(ans$mcd.wt) <- dimn[[1]]
        ans$wt <- NULL
        ans$X <- x
        if(names) {
            if(length(dimn[[1]]))
                dimnames(ans$X)[[1]] <- names(ans$mcd.wt)[ok]
            else
                dimnames(ans$X) <- list(seq(along = ok)[ok], NULL)
        }
        if(trace)
            cat(ans$method, "\n")
        ans$raw.cnp2 <- raw.cnp2
        ans$cnp2 <- cnp2
        class(ans) <- "mcd"
        return(ans)
    } ## end { alpha = 1   <==>   h = n }

    mcd <- if(nsamp == "deterministic") {
	ans$method <- paste("Deterministic", ans$method)
	.detmcd (x, h, hsets.init = initHsets,
		 save.hsets=save.hsets, # full.h=full.h,
		 scalefn=scalefn, maxcsteps=maxcsteps, trace=as.integer(trace),
		 names=names)
    } else {
	ans$method <- paste0("Fast ", ans$method, "; nsamp = ", nsamp,
			     "; (n,k)mini = (", nmini,",",kmini,")")
	# cgg 2015-12-10 add robustbase descriptor
	robustbase:::.fastmcd(x, h, nsamp, nmini, kmini, trace=as.integer(trace))
    }

    ## Compute the consistency correction factor for the raw MCD
    ##  (see calfa in Croux and Haesbroeck)
	# cgg 2015-12-10 add robustbase descriptor
    calpha <- robustbase::.MCDcons(p, h/n)    ## VT::19.3.2007
    correct <- if(use.correction) robustbase::.MCDcnp2(p, n, alpha) else 1.
    raw.cnp2 <- c(calpha, correct)

    if(p == 1) {
        ## ==> Compute univariate location and scale estimates
	ans$method <- paste("Univariate", ans$method)
        scale <- sqrt(calpha * correct) * as.double(mcd$initcovariance)
        center <- as.double(mcd$initmean)
        if(abs(scale - 0) < 1e-07) {
            ans$singularity <- list(kind = "identicalObs", q = h)
            ans$raw.cov <- ans$cov <- matrix(0)
            ans$raw.center <- ans$center <- center
            ans$n.obs <- n
            ans$alpha <- alpha
            ans$quan <- h
            if(names && !is.null(nms <- dimn[[2]][1])) {
                names(ans$raw.center) <- names(ans$center) <- nms
                dimnames(ans$raw.cov) <- dimnames(ans$cov) <- list(nms,nms)
            }
            ans$crit <- -Inf # = log(0)
            weights <- as.numeric(abs(x - center) < 1e-07) # 0 / 1
        } ## end { scale ~= 0 }
        else {
            ## Compute the weights for the raw MCD in case p=1
            # weights <- wgtFUN(((x - center)/scale)^2) # 0/1
			# cgg 2015-12-10 changed quantiel from value to function call
            weights <- as.numeric(((x - center)/scale)^2 < quantiel()) # 0/1
            sum.w <- sum(weights)
            ans <- c(ans, cov.wt(x, wt = weights, cor=cor))

	    if(sum.w != n) {
		# cgg 2015-12-10 add robustbase::
		cdelta.rew <- robustbase::.MCDcons(p, sum.w/n) ## VT::19.3.2007
		correct.rew <- if(use.correction) robustbase::.MCDcnp2.rew(p, n, alpha) else 1.
		cnp2 <- c(cdelta.rew, correct.rew)
		ans$cov <- cdelta.rew * correct.rew * ans$cov
	    }
            ans$alpha <- alpha
            ans$quan <- h
            ans$raw.cov <- as.matrix(scale^2)
            ans$raw.center <- as.vector(center)
            if(names && !is.null(nms <- dimn[[2]][1])) {
                dimnames(ans$raw.cov) <- list(nms,nms)
                names(ans$raw.center) <- nms
            }
	    ans$crit <- ## log(det) =
		log(sum(sort((x - as.double(mcd$initmean))^2, partial = h)[1:h])/max(1,h-1))
            center <- ans$center
            scale <- as.vector(sqrt(ans$cov))
            # weights <- wgtFUN(((x - center)/scale)^2)
			# cgg changed quantiel from value to function call
            weights <- as.numeric(((x - center)/scale)^2 < quantiel())
        } ## end{ scale > 0 }
    } ## end p=1

    else { ## p >= 2 : ---------------------------------------------------------

      ## Apply correction factor to the raw estimates
      ## and use them to compute weights
      mcd$initcovariance <- matrix(calpha * correct * mcd$initcovariance, p,p)
      if(raw.only || mcd$exactfit != 0) {
        ## If not all observations are in general position, i.e. more than
        ## h observations lie on a hyperplane, the program still yields
        ## the MCD location and scatter matrix, the latter being singular
        ## (as it should be), as well as the equation of the hyperplane.

        dim(mcd$coeff) <- c(5, p)
        ans$cov <- ans$raw.cov <- mcd$initcovariance
        ans$center <- ans$raw.center <- as.vector(mcd$initmean)

        if(names && !is.null(nms <- dimn[[2]])) {
            dimnames(ans$cov) <- list(nms, nms)
            names(ans$center) <- nms
        }
        ans$n.obs <- n

	if(raw.only) {
	    ans$raw.only <- TRUE
	} else {
	    ## no longer relevant:
	    ##      if(mcd$exactfit == -1)
	    ##      stop("The program allows for at most ", mcd$kount, " observations.")
	    ##      if(mcd$exactfit == -2)
	    ##      stop("The program allows for at most ", mcd$kount, " variables.")
	    if(!(mcd$exactfit %in% c(1,2,3)))
		stop("Unexpected 'exactfit' code ", mcd$exactfit, ". Please report!")
	    ## new (2007-01) and *instead* of older long 'method' extension;
	    ## the old message is still *printed* via .MCDsingularityMsg()
	    ##
	    ## exactfit is now *passed* to result instead of coded into 'message':
	    ans$singularity <-
		list(kind = "on.hyperplane", exactCode = mcd$exactfit,
		     p = p, h = h, count = mcd$kount, coeff = mcd$coeff[1,])
	}
        ans$alpha <- alpha
        ans$quan <- h
        if(names && !is.null(nms <- dimn[[2]])) {
            names(ans$raw.center) <- nms
            dimnames(ans$raw.cov) <- list(nms,nms)
        }
        ans$crit <- -Inf # = log(0)
        weights <- mcd$weights

      } ## end (raw.only || exact fit)

      else { ## have general position (exactfit == 0) : ------------------------

        ## FIXME? here, we assume that mcd$initcovariance is not singular:
        mah <- mahalanobis(x, mcd$initmean, mcd$initcovariance, tol = tolSolve)
        #weights <- wgtFUN(mah)
		# cgg 2015-12-10 changed quantiel from value to function call
        mcd$weights <- weights <- as.numeric(mah < quantiel())

# cgg
if ( reweighting$debug ) {
	cat("Raw mahalanobis Distances are:")
	print(sort(mah))
	cat("Method is", reweighting$method, "-", reweighting$quantile, "Quantile is", quantiel(),"\n")   
}     

        sum.w <- sum(weights)
        ans <- c(ans, cov.wt(x, wt = weights, cor=cor))
        ## simple check for singularity, much cheaper than determinant() below:
        sing.rewt <- any(apply(ans$cov == 0, 2, all))

        ## Compute and apply the consistency correction factor for
        ## the reweighted cov
        if(!sing.rewt && sum.w != n) {
		# cgg 2015-12-10 add robustbase:: here
	    cdelta.rew <- robustbase::.MCDcons(p, sum.w/n) ## VT::19.3.2007
	    correct.rew <- if(use.correction) robustbase::.MCDcnp2.rew(p, n, alpha) else 1.
	    cnp2 <- c(cdelta.rew, correct.rew)
	    ans$cov <- cdelta.rew * correct.rew * ans$cov
        }

        ##vt:: add also the best found subsample to the result list
        ans$best <- sort(as.vector(mcd$best))

        ans$alpha <- alpha
        ans$quan <- h
        ans$raw.cov <- mcd$initcovariance
        ans$raw.center <- as.vector(mcd$initmean)
        if(names && !is.null(nms <- dimn[[2]])) {
            names(ans$raw.center) <- nms
            dimnames(ans$raw.cov) <- list(nms,nms)
        }
        ans$raw.weights <- weights
        ans$crit <- mcd$mcdestimate # now in log scale!
        ## 'mah' already computed above
        ans$raw.mah <- mah ## mahalanobis(x, ans$raw.center, ans$raw.cov, tol = tolSolve)
        ## Check if the reweighted scatter matrix is singular.
        if(sing.rewt || - determinant(ans$cov, logarithm = TRUE)$modulus[1]/p > logdet.Lrg) {
	    ans$singularity <- list(kind = paste0("reweighted.MCD",
				    if(sing.rewt)"(zero col.)"))
			# cgg 2015-12-10
            if(trace) cat("The reweighted MCD scatter matrix is singular.\n")
	    ans$mah <- mah
        }
        else {
            mah <- mahalanobis(x, ans$center, ans$cov, tol = tolSolve)
            ans$mah <- mah
#if ( reweighting$debug ) {
#	cat("After reweighting mahalanobis Distances are:")
#	print(sort(mah))
#	cat("Method is", reweighting$method, "-", reweighting$quantile, "Quantile is", quantiel(),"\n")   
#}   		
			# cgg changed quantiel from value to function call
            #weights <- wgtFUN(mah)
            weights <- as.numeric(mah < quantiel())
        }
      } ## end{ not exact fit }

    } ## end{ p >= 2 }

    ans$mcd.wt <- rep.int(NA, length(ok))
    ans$mcd.wt[ok] <- weights
    if(names) {
	if(length(dimn[[1]]))
	    names(ans$mcd.wt) <- dimn[[1]]
	if(length(dimn[[1]]))
	    dimnames(x)[[1]] <- names(ans$mcd.wt)[ok]
	else
	    dimnames(x) <- list(seq(along = ok)[ok], NULL)
    }
    ans$X <- x
    ans$wt <- NULL
    if(trace)
        cat(ans$method, "\n")
    ans$raw.cnp2 <- raw.cnp2
    ans$cnp2 <- cnp2
    if(nsamp == "deterministic")
	ans <- c(ans, mcd[c("iBest","n.csteps", if(save.hsets) "initHsets")])
    class(ans) <- "mcd"
    ## warn if we have a singularity:
    if(is.list(ans$singularity))
	warning(paste(strwrap(.MCDsingularityMsg(ans$singularity, ans$n.obs)), collapse="\n"),
		domain=NA)
    ## return
    ans
} ## {covMcd}
