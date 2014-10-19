#### This is originally from the R package
####
####  rrcov : Scalable Robust Estimators with High Breakdown Point
####
#### by Valentin Todorov

covMcd2 <- function(x,
           cor = FALSE,
           alpha = 1/2,
           nsamp = 500,
           nmini = 300, # cgg added 2011-07-03
           seed = NULL,
           trace = FALSE,
           use.correction = TRUE,
           reweighting = list(method="chisquare",quantile=0.975,debug=FALSE),
           control = rrcov.control())
{

    require(robustbase, quietly=TRUE)

    ##   Analyze and validate the input parameters ...

    ## if a control object was supplied, take the option parameters
    ## from it, but if single parameters were passed (not defaults)
    ## they will override the control object.
    ## defCtrl <- rrcov.control()# default control
    # 'control' specified
    if(missing(alpha))  alpha <- control$alpha
    if(missing(nsamp))  nsamp <- control$nsamp
    if(missing(seed))    seed <- control$seed
    if(missing(trace))  trace <- control$trace
    if(missing(use.correction)) use.correction <- control$use.correction

    tolSolve <- control$tolSolve # had 1e-10 hardwired {now defaults to 1e-14}

    if(length(seed) > 0) {
        if(exists(".Random.seed", envir=.GlobalEnv, inherits=FALSE))  {
            seed.keep <- get(".Random.seed", envir=.GlobalEnv, inherits=FALSE)
            on.exit(assign(".Random.seed", seed.keep, envir=.GlobalEnv))
        }
        assign(".Random.seed", seed, envir=.GlobalEnv)
    }

    ##   vt::03.02.2006 - added options "best" and "exact" for nsamp
    ##   nsamp will be further analized in the wrapper .fastmcd()
    if(!missing(nsamp) && is.numeric(nsamp) && nsamp <= 0)
        stop("Invalid number of trials nsamp = ",nsamp, "!")


    if(is.data.frame(x))
        x <- data.matrix(x)
    else if (!is.matrix(x))
        x <- matrix(x, length(x), 1,
                    dimnames = list(names(x), deparse(substitute(x))))

    ## drop all rows with missing values (!!) :
    na.x <- !is.finite(x %*% rep(1, ncol(x)))
    ok <- !na.x
    x <- x[ok, , drop = FALSE]
    dx <- dim(x)
    if(!length(dx))
        stop("All observations have missing values!")
    dimn <- dimnames(x)
    n <- dx[1]
    p <- dx[2]
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
    else if(alpha > 1) stop("alpha must be <= 1")

	 ## cg::01.15.2010
	 ## modified this to allow for different quantile values and for an empirical
	 ## quantile
	 ## empirical quantile of mahalanobis distances cannot be computed ahead of time
	 ## so I've replaced the value quantiel with a function that gets the value
	 ##
	 reweighting$method <- match.arg(reweighting$method, c("chisquare","empirical"))
	 if ( is.null(reweighting$debug) ) reweighting$debug <- FALSE
	 quantiel <- if ( reweighting$method == "chisquare" ) {
		function() qchisq(reweighting$quantile, p)
	 } else {
	 	# when called, mah will be bound to the current value of mah
	 	function() quantile(mah, reweighting$quantile)  
	 }

    ## vt::03.02.2006 - raw.cnp2 and cnp2 are vectors of size 2 and  will
    ##   contain the correction factors (concistency and finite sample)
    ##   for the raw and reweighted estimates respectively. Set them
    ##   initially to 1.  If use.correction is set to FALSE
    ##   (default=TRUE), the finite sample correction factor will not
    ##   be used (neither for the raw estimates nor for the reweighted)
    raw.cnp2 <- cnp2 <- c(1,1)

    ans <- list(method = "Minimum Covariance Determinant Estimator.",
        call = match.call())

    if(alpha == 1) { ## alpha=1: Just compute the classical estimates --------
        mcd <- cov.wt(x)$cov
        loc <- as.vector(colMeans(x))
        obj <- determinant(mcd, log = TRUE)$modulus[1]
        if ( -obj/p > 50 ) {
            ans$cov <- mcd
            dimnames(ans$cov) <- list(dimn[[2]], dimn[[2]])
            if (cor)
                ans$cor <- cov2cor(ans$cov)
            ans$center <- loc
            if(length(dimn[[2]]))
                names(ans$center) <- dimn[[2]]
            ans$n.obs <- n
            ans$singularity <- list(kind = "classical")
            if(trace) cat("classical estimate is singular\n")

            weights <- 1
        }
        else {
            mah <- mahalanobis(x, loc, mcd, tol = tolSolve)
            ## VT:: 01.09.2004 - bug in alpha=1
            weights <- as.numeric(mah < quantiel()) # 0/1
            sum.w <- sum(weights)
            ans <- c(ans, cov.wt(x, wt = weights, cor = cor))
            ## ans$cov <- sum.w/(sum.w - 1) * ans$cov

            ## Consistency factor for reweighted MCD
            if(sum.w != n) {
                cnp2[1] <- robustbase:::MCDcons(p, sum.w/n)
                ans$cov <- ans$cov * cnp2[1]
            }
            if( - (determinant(ans$cov, log = TRUE)$modulus[1] - 0)/p > 50) {
                ans$singularity <- list(kind = "reweighted.MCD")
                if(trace) cat("reweighted MCD is singular\n")
            }
            else {
                mah <- mahalanobis(x, ans$center, ans$cov, tol = tolSolve)
                weights <- as.numeric(mah < quantiel()) # 0/1
            }
        }

        ans$alpha <- alpha
        ans$quan <- h
        ans$raw.cov <- mcd
        ans$raw.center <- loc
        if(!is.null(nms <- dimn[[2]])) {
            names(ans$raw.center) <- nms
            dimnames(ans$raw.cov) <- list(nms,nms)
        }
        ans$crit <- exp(obj)
        ans$method <-
            paste(ans$method,
                  "\nThe minimum covariance determinant estimates based on",
                  n, "observations \nare equal to the classical estimates.")
        ans$mcd.wt <- rep(NA, length(ok))
        ans$mcd.wt[ok] <- weights
        if(length(dimn[[1]]))
            names(ans$mcd.wt) <- dimn[[1]]
        ans$wt <- NULL
        ans$X <- x
        if(length(dimn[[1]]))
            dimnames(ans$X)[[1]] <- names(ans$mcd.wt)[ok]
        else
            dimnames(ans$X) <- list(seq(along = ok)[ok], NULL)
        if(trace)
            cat(ans$method, "\n")
        ans$raw.cnp2 <- raw.cnp2
        ans$cnp2 <- cnp2
        class(ans) <- "mcd"
        return(ans)
    } ## end {alpha=1} --

    mcd <- robustbase:::.fastmcd(x, h, nsamp, nmini, trace=as.integer(trace))

    ## Compute the consistency correction factor for the raw MCD
    ##  (see calfa in Croux and Haesbroeck)
    calpha <- robustbase:::MCDcons(p, h/n)    ## VT::19.3.2007
    correct <- if(use.correction) robustbase:::MCDcnp2(p, n, alpha) else 1.
    raw.cnp2 <- c(calpha, correct)

    if(p == 1) {
        ## ==> Compute univariate location and scale estimates
        ans$method <- "Univariate location and scale estimation."

        scale <- sqrt(calpha * correct) * as.double(mcd$initcovariance)
        center <- as.double(mcd$initmean)
        if(abs(scale - 0) < 1e-07) {
            ans$singularity <- list(kind = "identicalObs", q = h)
            ans$raw.cov <- ans$cov <- matrix(0)
            ans$raw.center <- ans$center <- center
            ans$n.obs <- n
            ans$alpha <- alpha
            ans$quan <- h
            if(!is.null(nms <- dimn[[2]][1])) {
                names(ans$raw.center) <- names(ans$center) <- nms
                dimnames(ans$raw.cov) <- dimnames(ans$cov) <- list(nms,nms)
            }
            ans$crit <- 0
            weights <- as.numeric(abs(x - center) < 1e-07) # 0 / 1
        } ## end { scale ~= 0 }
        else {
            ## Compute the weights for the raw MCD in case p=1
            weights <- as.numeric(((x - center)/scale)^2 < quantiel()) # 0/1
            sum.w <- sum(weights)
            ans <- c(ans, cov.wt(x, wt = weights, cor = cor))
            ## ans$cov <- sum.w/(sum.w - 1) * ans$cov

            ## Apply the correction factor for the reweighted cov
            if(sum.w == n) {
                cdelta.rew <- 1
                correct.rew <- 1
            }
            else {
                cdelta.rew  <- robustbase:::MCDcons(p, sum.w/n) ## VT::19.3.2007
                correct.rew <- if(use.correction) robustbase:::MCDcnp2.rew(p, n, alpha) else 1.
                cnp2 <- c(cdelta.rew, correct.rew)
            }
            ans$cov <- ans$cov * cdelta.rew * correct.rew
            ans$alpha <- alpha
            ans$quan <- h
            ans$raw.cov <- as.matrix(scale^2)
            ans$raw.center <- as.vector(center)
            if(!is.null(nms <- dimn[[2]][1])) {
                dimnames(ans$raw.cov) <- list(nms,nms)
                names(ans$raw.center) <- nms
            }
            ans$crit <- 1/(h - 1) *
                sum(sort((x - as.double(mcd$initmean))^2, partial = h)[1:h])
            center <- ans$center
            scale <- as.vector(sqrt(ans$cov))
            weights <- as.numeric(((x - center)/scale)^2 < quantiel())
        } ## end{ scale > 0 }
    } ## end p=1

    else { ## p >= 2 : ---------------------------------------------------------

      ## Apply correction factor to the raw estimates
      ## and use them to compute weights
      mcd$initcovariance <- calpha * mcd$initcovariance * correct
      dim(mcd$initcovariance) <- c(p, p)

      if(mcd$exactfit != 0) {
        ## If not all observations are in general position, i.e. more than
        ## h observations lie on a hyperplane, the program still yields
        ## the MCD location and scatter matrix, the latter being singular
        ## (as it should be), as well as the equation of the hyperplane.

        dim(mcd$coeff) <- c(5, p)
        ans$cov <- mcd$initcovariance
        ans$center <- as.vector(mcd$initmean)
        if(!is.null(nms <- dimn[[2]])) {
            dimnames(ans$cov) <- list(nms, nms)
            names(ans$center) <- nms
        }
        ans$n.obs <- n

  ## no longer relevant:
  ##      if(mcd$exactfit == -1)
  ##      stop("The program allows for at most ", mcd$kount, " observations.")
  ##      if(mcd$exactfit == -2)
  ##      stop("The program allows for at most ", mcd$kount, " variables.")
        if(!(mcd$exactfit %in% c(1,2)))
            stop("Unexpected 'exactfit' code ", mcd$exactfit, ". Please report!")
        ## new (2007-01) and *instead* of older long 'method' extension;
        ## the old message is stilled *printed* via singularityMessage()
        ##
        ## exactfit is now *passed* to result instead of coded into 'message':
        ans$singularity <-
            list(kind = "on.hyperplane", exactCode = mcd$exactfit,
                 p = p, count = mcd$kount, coeff = mcd$coeff[1,])
        ans$alpha <- alpha
        ans$quan <- h
        ans$raw.cov <- mcd$initcovariance
        ans$raw.center <- as.vector(mcd$initmean)
        if(!is.null(nms <- dimn[[2]])) {
            names(ans$raw.center) <- nms
            dimnames(ans$raw.cov) <- list(nms,nms)
        }
        ans$crit <- 0
        weights <- mcd$weights

      } ## end exact fit <==>  (mcd$exactfit != 0)

      else { ## exactfit == 0 : have general position ------------------------

        ## FIXME: here we assume that mcd$initcovariance is not singular
        ## ----- but it is for data(mortality, package = "riv") !
        mah <- mahalanobis(x, mcd$initmean, mcd$initcovariance, tol = tolSolve)
        mcd$weights <- weights <- as.numeric(mah < quantiel())


if ( reweighting$debug ) {
	cat("Raw mahalanobis Distances are:")
	print(sort(mah))
	cat("Method is", reweighting$method, "-", reweighting$quantile, "Quantile is", quantiel(),"\n")   
}     
        sum.w <- sum(weights)

        ## Compute and apply the consistency correction factor for
        ## the reweighted cov
        if(sum.w == n) {
            cdelta.rew <- 1
            correct.rew <- 1
        }
        else {
            cdelta.rew <- robustbase:::MCDcons(p, sum.w/n) ## VT::19.3.2007
            correct.rew <- if(use.correction) robustbase:::MCDcnp2.rew(p, n, alpha) else 1.
            cnp2 <- c(cdelta.rew, correct.rew)
        }

        ans <- c(ans, cov.wt(x, wt = weights, cor))
        ## ans$cov <- sum.w/(sum.w - 1) * ans$cov
        ans$cov <- ans$cov * cdelta.rew * correct.rew

        ##vt:: add also the best found subsample to the result list
        ans$best <- sort(as.vector(mcd$best))

        ans$alpha <- alpha
        ans$quan <- h
        ans$raw.cov <- mcd$initcovariance
        ans$raw.center <- as.vector(mcd$initmean)
        if(!is.null(nms <- dimn[[2]])) {
            names(ans$raw.center) <- nms
            dimnames(ans$raw.cov) <- list(nms,nms)
        }
        ans$raw.weights <- weights
        ans$crit <- mcd$mcdestimate
        ans$raw.mah <- mahalanobis(x, ans$raw.center, ans$raw.cov, tol = tolSolve)
  
        ## Check if the reweighted scatter matrix is singular.
        if( - (determinant(ans$cov, log = TRUE)$modulus[1] - 0)/p > 50) {
            ans$singularity <- list(kind = "reweighted.MCD")
            if(trace) cat("The reweighted MCD scatter matrix is singular.\n")
            ans$mah <- ans$raw.mah
        }
        else {
            mah <- mahalanobis(x, ans$center, ans$cov, tol = tolSolve)
            ans$mah <- mah

#if ( reweighting$debug ) {
#	cat("After reweighting mahalanobis Distances are:")
#	print(sort(mah))
#	cat("Method is", reweighting$method, "-", reweighting$quantile, "Quantile is", quantiel(),"\n")   
#}   
            weights <- as.numeric(mah < quantiel())
        }
      } ## end{ not exact fit }

    } ## end{ p >= 2 }

    ans$mcd.wt <- rep(NA, length(ok))
    ans$mcd.wt[ok] <- weights
    if(length(dimn[[1]]))
        names(ans$mcd.wt) <- dimn[[1]]
    ans$wt <- NULL
    ans$X <- x
    if(length(dimn[[1]]))
        dimnames(ans$X)[[1]] <- names(ans$mcd.wt)[ok]
    else
        dimnames(ans$X) <- list(seq(along = ok)[ok], NULL)
    if(trace)
        cat(ans$method, "\n")
    ans$raw.cnp2 <- raw.cnp2
    ans$cnp2 <- cnp2
    class(ans) <- "mcd"
    ## even when not 'trace': say if we have singularity!
    if(is.list(ans$singularity))
        cat(strwrap(robustbase:::singularityMsg(ans$singularity, ans$n.obs)), sep ="\n")

    return(ans)
}
