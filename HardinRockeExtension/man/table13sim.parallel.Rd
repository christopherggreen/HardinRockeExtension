\name{table13sim.parallel}
\alias{table13sim.parallel}
\title{Replicate the experiment presented in Cerioli et al. (2009)}
\description{
  Replicates the experiment presented in Cerioli et al. (2009), 
  Tables 1 and 3, for a wider variety of estimators.
}
\usage{
table13sim.parallel(cl, p, nn, N, B = 250, alpha = c(0.01, 0.025, 0.05), lgf = "", mlgf = "", maxtries = 100)
}
\arguments{
  \item{cl}{A cluster object, e.g., returned from \code{makePSOCKcluster}.
    The user must create this object before calling \code{hrSimNewParallel}.}
  \item{p}{The dimension of the data used in each simulated run.}
  \item{nn}{The number of observations used in each simulated run.}
  \item{N}{The number of simulations to run.}
  \item{B}{The batch/block size: the number of simulations to run 
    in each block. This is useful when running very large
    simulation runs (\code{N} very large) where memory is a concern.}
  \item{alpha}{The significance level to use for detecting outliers. Can
    be a vector; the outlier detection tests will be run at each level.}
  %\item{cutoff.method}{String indicating with asymptotic distribution to use
  %  for the MCD-based distances. Valid values are \code{"HR05"} for the
  %  method of Hardin and Rocke (2005), and \code{"GM14"} for the method
  %  of Green and Martin (2014). Default is \code{"GM14"}.}
  \item{lgf}{Path to log file into which logging information should be
    written.}
  \item{mlgf}{ not used at this time }
  \item{maxtries}{The maximum number of times to retry failed blocks. The Rocke
    S-estimator can fail when \eqn{n/p} is small if it gets a bad random sample, 
	so we restart such blocks. Default is 100.}
}
\details{
  This is a work function designed for use in replicating Tables 1 and 3
  of Cerioli et al. (2009), page XXX, but using the asymptotic method
  of Green and Martin (2014) instead of the Hardin-Rocke method. The
  experiment investigates how many false-positives certain Mahalanobis-based
  tests of outlyingness produce, compared to the nominal Type I error
  rate \eqn{\alpha}.

  Internally the simulation function does \code{B} runs at a time. Blocks
  of size \code{B} are distributed across the cluster. Set \code{B} smaller if 
  your machines have less memory or you have lots of cluster nodes.
}
\value{An array of dimension 3:
  \enumerate{
    \item{The results of each of the \code{N} simulation runs appear along the first dimension.}
	\item{The various estimators and tests appear along the second dimension. Results with
	  suffix ``T1'' correspond to Table 1 of Cerioli et al. (2009) (the individual outlier tests)
	  while those with suffix ``T3'' correspond to Table 3 (the simultaneous outlier tests).
	  Currently the 74 columns appear in the following order.
	  \tabular{lll}{
	    Column Name \tab Covariate Estimate \tab Test Statistic \cr
        "OGK.T1" \tab OGK estimate (\eqn{\beta = 0.9}) \tab chi-squared \cr 
        "ROGK.T1" \tab Reweighted OGK estimate \tab chi-squared \cr 
        "SEST.BS.T1" \tab S-estimate using bisquare ? \tab chi-squared \cr 
        "SEST.RK.T1" \tab S-estimate using Rocke ? (\code{arp = 0.05})\tab chi-squared \cr 
        "MCDMBP.RAW.T1" \tab MCD (max. breakdown pt.) \tab chi-squared \cr 
        "MCDMBP.GMRAW.T1" \tab MCD (max. breakdown pt.) \tab Green-Martin \cr 
        "MCDMBP.GMADJ.T1" \tab MCD (max. breakdown pt.) \tab Green-Martin (adj.) \cr 
        "RMCDMBP.T1" \tab Reweighted MCD (max. breakdown pt.) \tab chi-squared \cr 
        "MCD75.RAW.T1" \tab MCD(0.75) \tab chi-squared \cr 
        "MCD75.GMRAW.T1" \tab MCD(0.75) \tab Green-Martin \cr 
        "MCD75.GMADJ.T1" \tab MCD(0.75) \tab Green-Martin (adj.) \cr 
        "RMCD75.T1" \tab Reweighted MCD(0.75) \tab chi-squared \cr 
        "MCD95.RAW.T1" \tab MCD(0.95) \tab chi-squared \cr 
        "MCD95.GMRAW.T1" \tab MCD(0.95) \tab Green-Martin \cr 
        "MCD95.GMADJ.T1" \tab MCD(0.95) \tab Green-Martin (adj.) \cr 
        "RMCD95.T1" \tab Reweighted MCD(0.95) \tab chi-squared \cr 
        "OGK.T3" \tab OGK estimate \tab chi-squared \cr 
        "ROGK.T3" \tab Reweighted OGK estimate \tab chi-squared \cr 
        "ROGK.CH.T3" \tab Reweighted OGK estimate using Bonferroni corrected \eqn{\beta} \tab chi-squared \cr 
        "SEST.BS.T3" \tab S-estimate using bisquare ? \tab chi-squared \cr 
        "SEST.RK.T3" \tab S-estimate using Rocke ? \tab chi-squared \cr 
        "MCDMBP.RAW.T3" \tab MCD (max. breakdown pt.) \tab chi-squared \cr 
        "MCDMBP.GMRAW.T3" \tab MCD (max. breakdown pt.) \tab Green-Martin \cr 
        "MCDMBP.GMADJ.T3" \tab MCD (max. breakdown pt.) \tab Green-Martin (adj.) \cr 
        "MCDMBP.HRADJ.T3" \tab MCD (max. breakdown pt.) \tab Hardin-Rocke (adj.) \cr 
        "RMCDMBP.T3" \tab Reweighted MCD (max. breakdown pt.) \tab chi-squared \cr 
        "RMCDMBP.CH.T3" \tab Reweighted MCD (max. breakdown pt.) with Bonferroni correction \tab chi-squared \cr 
        "MCD75.RAW.T3" \tab MCD(0.75) \tab chi-squared \cr 
        "MCD75.GMRAW.T3" \tab MCD(0.75) \tab Green-Martin \cr 
        "MCD75.GMADJ.T3" \tab MCD(0.75) \tab Green-Martin (adj.) \cr 
        "RMCD75.T3" \tab Reweighted MCD(0.75) \tab chi-squared \cr 
        "RMCD75.CH.T3" \tab Reweighted MCD(0.75) with Bonferroni correction \tab chi-squared \cr 
        "MCD95.RAW.T3" \tab MCD(0.95)  \tab chi-squared \cr 
        "MCD95.GMRAW.T3" \tab MCD(0.95) \tab Green-Martin \cr 
        "MCD95.GMADJ.T3" \tab MCD(0.95) \tab Green-Martin (adj.) \cr 
        "RMCD95.T3" \tab Reweighted MCD(0.95) \tab chi-squared \cr 
        "RMCD95.CH.T3" \tab Reweighted MCD(0.95) with Bonferroni correction \tab chi-squared \cr 
        "MCDMBP.HRRAW.T1" \tab MCD (max. breakdown pt.) \tab Hardin-Rocke \cr 
        "MCDMBP.HRADJ.T1" \tab MCD (max. breakdown pt.) \tab Hardin-Rocke (adj.) \cr 
        "MCDMBP.HRRAW.T3" \tab MCD (max. breakdown pt.) \tab Hardin-Rocke \cr 
        "MCD75.HRRAW.T1" \tab MCD(0.75) \tab Hardin-Rocke \cr 
        "MCD75.HRADJ.T1" \tab MCD(0.75) \tab Hardin-Rocke (adj.) \cr 
        "MCD75.HRRAW.T3" \tab MCD(0.75) \tab Hardin-Rocke \cr 
        "MCD75.HRADJ.T3" \tab MCD(0.75) \tab Hardin-Rocke (adj.) \cr 
        "MCD95.HRRAW.T1" \tab MCD(0.95) \tab Hardin-Rocke \cr 
        "MCD95.HRADJ.T1" \tab MCD(0.95) \tab Hardin-Rocke (adj.) \cr 
        "MCD95.HRRAW.T3" \tab MCD(0.95) \tab Hardin-Rocke \cr 
        "MCD95.HRADJ.T3" \tab MCD(0.95) \tab Hardin-Rocke (adj.) \cr 
        "OGK.T3.ALT" \tab OGK estimate \tab chi-squared \cr 
        "ROGK.T3.ALT" \tab Reweighted OGK estimate \tab chi-squared \cr 
        "ROGK.CH.T3.ALT" \tab Reweighted OGK estimate with Bonferroni correction \tab chi-squared \cr 
        "SEST.BS.T3.ALT" \tab S-estimate with bisquare ? \tab chi-squared \cr 
        "SEST.RK.T3.ALT" \tab S-estimate with Rocke ? \tab chi-squared \cr 
        "MCDMBP.RAW.T3.ALT" \tab MCD (max. breakdown pt.) \tab chi-squared \cr 
        "MCDMBP.GMRAW.T3.ALT" \tab MCD (max. breakdown pt.) \tab Green-Martin \cr 
        "MCDMBP.HRRAW.T3.ALT" \tab MCD (max. breakdown pt.) \tab Hardin-Rocke \cr 
        "MCDMBP.GMADJ.T3.ALT" \tab MCD (max. breakdown pt.) \tab Green-Martin (adj.) \cr 
        "MCDMBP.HRADJ.T3.ALT" \tab MCD (max. breakdown pt.) \tab Hardin-Rocke (adj.) \cr 
        "RMCDMBP.T3.ALT" \tab Reweighted MCD (max. breakdown pt.) \tab chi-squared \cr 
        "RMCDMBP.CH.T3.ALT" \tab Reweighted MCD (max. breakdown pt.) with Bonferroni correction \tab chi-squared \cr 
        "MCD75.RAW.T3.ALT" \tab MCD(0.75) \tab chi-squared \cr 
        "MCD75.GMRAW.T3.ALT" \tab MCD(0.75) \tab Green-Martin \cr 
        "MCD75.HRRAW.T3.ALT" \tab MCD(0.75) \tab Hardin-Rocke \cr 
        "MCD75.GMADJ.T3.ALT" \tab MCD(0.75) \tab Green-Martin (adj.) \cr 
        "MCD75.HRADJ.T3.ALT" \tab MCD(0.75) \tab Hardin-Rocke (adj.) \cr 
        "RMCD75.T3.ALT" \tab Reweighted MCD(0.75) \tab chi-squared \cr 
        "RMCD75.CH.T3.ALT" \tab Reweighted MCD(0.75) with Bonferroni correction \tab chi-squared \cr 
        "MCD95.RAW.T3.ALT" \tab MCD(0.95) \tab chi-squared \cr 
        "MCD95.GMRAW.T3.ALT" \tab MCD(0.95) \tab Green-Martin \cr 
        "MCD95.HRRAW.T3.ALT" \tab MCD(0.95) \tab Hardin-Rocke \cr 
        "MCD95.GMADJ.T3.ALT" \tab MCD(0.95) \tab Green-Martin (adj.) \cr 
        "MCD95.HRADJ.T3.ALT" \tab MCD(0.95) \tab Hardin-Rocke (adj.) \cr 
        "RMCD95.T3.ALT" \tab Reweighted MCD(0.95) \tab chi-squared \cr 
        "RMCD95.CH.T3.ALT" \tab Reweighted MCD(0.95) with Bonferroni correction \tab chi-squared \cr 
      }
	  The adjusted versions of the Hardin-Rocke tests remove the finite
	  sample correction when the sample size is 100 or greater. (WHY DID
	  WE DO THIS)

	  Bonferroni quantile
	  \deqn{ ( 1 - \alpha/nn ) }
	  Signficance level is \eqn{\alpha/nn}.

	  The alt versions use the chi-squared but with the correction
	  \deqn{ ( 1 - \alpha )^(1/nn) }.
	  Significance level is \eqn{1 - ((1-\alpha)^(1/nn))}.
	}
	\item{The specified values of \code{alpha} correspond to the third 
	  dimension; the dimnames will be of the form ``alpha'' + \code{alpha}.}
  }
}
\references{
  Andrea Cerioli, Marco Riani, and Anthony C. Atkinson. Controlling the size of multivariate
  outlier tests with the mcd estimator of scatter. Statistical Computing, 19:341-353, 2009.

  C. G. Green and R. Douglas Martin. An extension of a method of Hardin and Rocke, with
    an application to multivariate outlier detection via the IRMCD method of Cerioli.
    Working Paper, 2014. Available from 
    \url{http://students.washington.edu/cggreen/uwstat/papers/cerioli_extension.pdf}

  J. Hardin and D. M. Rocke. The distribution of robust distances. 
    Journal of Computational and Graphical Statistics, 14:928-946, 2005.
}
\author{
Written and maintained by Christopher G. Green <christopher.g.green@gmail.com>
}
\seealso{
	\code{\link{table1sim.parallel}}, \code{\link{table13sim.parallel.check}},
	\code{\link{CovMcd2}}
}
\examples{
  \dontrun{
    # this runs an experiment
	# assumes a cluster
	# the vignette provides a better recipe for 
	# replicating Cerioli et al. (2009)

    require( parallel                )
    require( CerioliOutlierDetection )
    require( HardinRockeExtension    )
    
	# we use a socket cluster on Windows,
	# change to your preferred method of
	# creating a cluster
    thecluster <- makePSOCKcluster(4)

    N.SIM <- 500
    B.SIM <- 50
    
    # initialize each node
    tmp.rv <- clusterEvalQ( cl = thecluster, {
    
      require(abind,                   quietly=TRUE)
      require(rrcov,                   quietly=TRUE)
      require(mvtnorm,                 quietly=TRUE)
      require(CerioliOutlierDetection, quietly=TRUE)
      require(HardinRockeExtension,    quietly=TRUE)
    
      Sys.sleep(30)
    
      invisible(NULL)
    })
    
    results <- table1sim.parallel(cl=thecluster, p = 4, nn = 300, 
          N=500, B=50, lgf=logfile)
    stopCluster(thecluster)

    # calculate some statistics 
    apply(results,c(2,3),mean),
    apply(results,c(2,3),sd)
  }
} % examples
\keyword{ robust }
\keyword{ multivariate }
