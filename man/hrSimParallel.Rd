\name{hrSimParallel}
\alias{hrSimParallel}
\title{
  Replicate Tables 1--3 of Hardin and Rocke (2005)
}
\description{
  \code{hrSimParallel} is used to replicate Tables 1--3 of Hardin
  and Rocke (2005), page 942.
}
\usage{
hrSimParallel(cl, p, n, N, B = 10000, alpha = 0.05, mcd.alpha = max.bdp.mcd.alpha(n, p), lgf = "")
}
\arguments{
  \item{cl}{A cluster object, e.g., returned from \code{makePSOCKcluster}.
    The user must create this object before calling \code{hrSimNewParallel}.}
  \item{p}{The dimension of the data used in each simulated run.}
  \item{n}{The number of observations used in each simulated run.}
  \item{N}{The number of simulations to run.}
  \item{B}{The batch/block size: the number of simulations to run 
    in each block. This is useful when running very large
    simulation runs (\code{N} very large) where memory is a concern.}
  \item{alpha}{The significance level to use for detecting outliers.}
  \item{mcd.alpha}{The fraction of the data to use in computing the
    MCD. Defaults to the maximum breakdown point fraction.}
  \item{lgf}{Path to log file into which logging information should be
    written.}
}
\details{
  This is a work function designed for use in replicating Tables 1--3 of Hardin
  and Rocke (2005), pages 942-944. Use different values of \code{alpha} to
  replicate each of the tables.

  Internally the simulation function does \code{B} runs at a time. Set 
  \code{B} smaller if your machine has less memory.

  This function performs the same calculation as \code{\link{hrSim}},
  but does so using internal parallelism---multiple blocks of size 
  \code{B} are run in parallel.
}
\value{

  The function returns a matrix with \code{N} rows, one for each simulation,
  and at present, 9 columns: each column reports the fraction of observations
  in a simulation run that exceeded a given threshold (i.e., were flagged as 
  outliers). 
  \enumerate{ 
  \item{The first three test Mahalanobis distances (MD)
  against a chi-squared quantile (prefix is ``CHI2'');}
  \item{the next three test MDs against the asympotic
  cutoff used in Hardin and Rocke (2005) (prefix is ``HRASY''); and}
  \item{the last three test against the cutoff predicted 
  in Hardin and Rocke (2005) (prefix is ``HRPRED'').}
  } 
  
  Within each group of three, the first entry (suffix ``RAW'') uses (raw)
  MDs without the consistency correction or the small sample correction;
  the second entry (suffix ``CON'') uses (raw) MDs without the small
  sample correction; and the third entry (suffix ``SM'') uses the (raw)
  MDs with both correction factors. (It was not clear to the package
  author whether Hardin and Rocke (2005) used these correction factors in 
  their calculations; so all variants were calculated and examined. 
  Empirically, it seems the ``CON'' approach is the best match for their 
  results.)

  Look at the column means of the resulting matrix to see the average 
  fraction of outliers detected (which is an estimate of the Type 1 error rate
  of the procedure, since the simulated data had no outliers).

  The vignette ``HardinRocke'' provides a detailed example of how to
  replicate the data in Hardin and Rocke (2005).
}
\references{
J. Hardin and D. M. Rocke. The distribution of robust distances. 
Journal of Computational and Graphical Statistics, 14:928-946, 2005.
}
\author{
Written and maintained by Christopher G. Green <christopher.g.green@gmail.com>
}
\seealso{
\code{\link{hrSim}}, \code{\link{hrSimNewParallel}}
}
\examples{
  \dontrun{
    # example of how to replicate some of the 
	# calculations in Hardin and Rocke (2005)
    
    # on Windows you must use socket clusters
    # Linux/UNIX supports other types of clusters
    # 
    # Change '4' to reflect the number of 
    # cores/processors you want to use
    require( parallel )
    thecluster <- makePSOCKcluster(4)
  
    # initialize each node
    tmp.rv <- clusterEvalQ( cl = thecluster, {
      require( CerioliOutlierDetection )
      require( HardinRockeExtensionSimulations )

      invisible(NULL)
    })
    
	# compare to Hardin and Rocke, Table 1
    results <- hrSimParallel(cl=thecluster, p = 5, n = 500, 
		N=5000, B=125, lgf="logfile.txt")
	colMeans(results)
    
    stopCluster(thecluster)
	
  }
}
\keyword{ robust }
\keyword{ multivariate }
