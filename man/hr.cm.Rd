\name{hr.cm}
\alias{hr.cm}
\title{ Hardin and Rocke (2005) corrections for general MCD estimates }
\description{
This is a work function to calculate, via simulation, the (reciprocal)
consistency correction ``c'' and the Wishart degrees of freedom ``m'' 
needed in the Hardin and Rocke approximate distribution for the minimum
covariance estimate (MCD) covariance. This function generalizes
their work to MCD using data fractions other than the maximum breakdown 
point case.
}
\usage{
hr.cm(p, n, N, B = 10000, 
  mcd.alpha = max.bdp.mcd.alpha, logfile = "hr_cm.log")
}
\arguments{
  \item{p}{The dimension of the data used in each simulated run.}
  \item{n}{The number of observations used in each simulated run.}
  \item{N}{The number of simulations to run.}
  \item{B}{The batch/block size: the number of simulations to run 
    in each block. This is useful when running very large
    simulation runs (\code{N} very large) where memory is a concern.}
  \item{alpha}{The significance level to use for detecting outliers.}
  \item{mcd.alpha}{The fraction of the data to use in computing
    MCD. Defaults to the maximum breakdown point fraction. Can
  be a vector of values, in which case the calculations will
  be done for each value of \code{mcd.alpha}.}
  \item{logfile}{Name of file into which to write logging information.}
}
\details{
  For each simulated data set, calculate
  the covariance of the MCD subset. Then
  we save the sum of the diagonal entries
  of said covariance matrix and the sum of the 
  squares of the diagonal entries. The estimate
  of ``c'', the reciprocal consistency constant, 
  is then given by the mean of the sums of the 
  diagonal entries. The estimate of ``m'' uses this 
  estimated ``c'' and an estimate of the variance 
  of the diagonal elements.

  Based on Johanna Hardin's \code{cm.R} and 
  \code{mcd.est.R} code.
}
\value{A data frame with the following columns
  \item{c}{The estimate of ``c'', the reciprocal consistency 
    constant}
  \item{m}{The estimated Wishart degrees of freedom ``m''}
  \item{mcd.alpha}{The input data fraction(s) \code{mcd.alpha}}
}
\references{
J. Hardin and D. M. Rocke. The distribution of robust distances. 
Journal of Computational and Graphical Statistics, 14:928-946, 2005.

J. Hardin. R code: to estimate the MCD – \code{mcd.est.r} and 
to estimate c and m – \code{cm.r}. Available
from \url{http://pages.pomona.edu/~jsh04747/Research/cm.r}
and \url{http://pages.pomona.edu/~jsh04747/Research/mcd.est.r}
respectively. 2005. Accessed December 10, 2015.
}
\author{
Written and maintained by Christopher G. Green <christopher.g.green@gmail.com>.
Based on code by Johanna S. Hardin.
}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
  \dontrun{
    # Monte Carlo simulation to compute the
    # Wishart degrees of freedom m associated
    # with the Hardin-Rocke approximation to 
    # the MCD covariance estimate
    #
    # See Hardin and Rocke (2005) for the 
    # original simulation in the maximum
    # breakdown point case. 
    # 
    # This script extends their work to 
    # sampling fractions greater than the
    # maximum breakdown point case 
    # floor( ( n+v+1 )/2 )/n.
    #
	
    require(parallel)
    
	# on Windows we use socket clusters
	# use whatever parallel backend works best for you
	# change '4' to the number of nodes/workers/etc. you
	# want to use

    thecluster <- makePSOCKcluster(4)
    # make this reproducible
    clusterSetRNGStream(cl = thecluster, 2014)
    
    # initialize each node
    tmp.rv <- clusterEvalQ( cl = thecluster, {
      require( CerioliOutlierDetection )
      require( HardinRockeExtensionSimulations )
      N.SIM <- 5000
      B.SIM <- 250
     
      my.pid <- Sys.getpid()
      cat("My pid is ", my.pid, "\n")
      logfile <- paste("Simulation_AllAlphas_Parallel_logfile_",my.pid,".txt",sep="")
      cat("Initialized\n\n", file=logfile)
    
      invisible(NULL)
    })
    
    # build the pairs of sample size n and dimension p
    hr.cm.params <- expand.grid(
          list(
          p=c(3,5,7,10,15,20),
          n=c(50,100,250,500,750,1000)
        )
        )
    # adding more coverage for small sample sizes
    hr.cm.params <- rbind( hr.cm.params, within( 
      expand.grid(list(p=c(3,5,7,10,15,20), ratio=c( 3,5,7,9,11 ) )), 
      {
        n <- p * ratio
        rm(ratio)
      }
    ))
    # remove any duplicates
    hr.cm.params <- unique(hr.cm.params)
    # want to run most expensive cases first
    hr.cm.params <- hr.cm.params[ order( hr.cm.params$n, hr.cm.params$p, decreasing=TRUE ), ]
    
    # add maximum breakdown point case to the params data set
    hr.cm.params[,"mbp"] <- apply( hr.cm.params, 1, function(x) floor( (x[2] + x[1] + 1)/2 )/x[2] )
    
    # want each case to be a column so that we can use parLapply
    hr.cm.params <- data.frame(t(as.matrix(hr.cm.params)))
          
    mcd.alphas <- c(0.55,0.60,0.65,0.70,0.75,0.80,0.85,0.90,0.95,0.99,0.995) 
    clusterExport(cl = thecluster, "hr.cm.params")
    clusterExport(cl = thecluster, "mcd.alphas")
    
    #
    # using parLapply here to prevent simplification of the
    # results (as parApply would attempt to do)
    #
    cat("Starting run at ", format(Sys.time()), "\n")
    
    hr.cm.results.all.pre <- parLapply(cl = thecluster, 
      X = hr.cm.params[1:5], function(pn) {
        cat("Starting case p = ",pn[1]," n = ",pn[2]," at time ", 
          format(Sys.time()), " \n",file=logfile,append=TRUE)
        results <- hr.cm(p = pn[1] , n = pn[2], N=N.SIM, B=B.SIM, 
          mcd.alpha=unique(c(pn[3],mcd.alphas)), logfile=logfile)
        cat("Finished case p = ",pn[1]," n = ",pn[2]," at time ", 
          format(Sys.time()), " \n",file=logfile,append=TRUE)
        data.frame(p=pn[1],n=pn[2],mbp=pn[3],results)
      }
    )
    cat("Run completed at ", format(Sys.time()), "\n")
    
    stopCluster(thecluster)
    
    hr.cm.results.all <- do.call("rbind", hr.cm.results.all.pre )
    
    save("hr.cm.results.all", file="hr.cm.results.all.final.rda")
  }    
}
\keyword{ robust }
\keyword{ multivariate }
