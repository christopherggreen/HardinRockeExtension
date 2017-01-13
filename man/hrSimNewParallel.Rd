\name{hrSimNewParallel}
\alias{hrSimNewParallel}
\title{
  Calculate Tables 1--3 of Hardin and Rocke (2005) using the
  method of Green and Martin (2014)
}
\description{
  \code{hrSimNewParallel} is used to calculate statistics like those presented
  in Tables 1--3 of Hardin and Rocke (2005), page 942, but with the 
  modified asymptotic formulas calculated in Green and Martin (2014).

  \code{hrSimNewParallel} is a parallel version of \code{hrSimNew}.
}
\usage{
hrSimNewParallel(cl, p, n, N, B = 10000, alpha = 0.05, mcd.alpha = max.bdp.mcd.alpha(n, p), lgf = "")
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
  This is a work function designed for use in replicating Tables 1--3 
  of Hardin and Rocke (2005), page 942, but using the asymptotic method of 
  Green and Martin (2014) instead of the Hardin-Rocke method.

  Internally the simulation function does \code{B} runs at a time. Set 
  \code{B} smaller if your machine has less memory.

  This function is nearly identical to \code{\link{hrSim}}, except that
  it uses different cutoff values.

  This function performs the same calculation as \code{\link{hrSimNew}},
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
  cutoff used in Green and Martin (2014) (prefix is ``CGASY''); and}
  \item{the last three test against the cutoff predicted 
  in Green and Martin (2014) (prefix is ``GGPRED'').}
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

}
\references{
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
\code{\link{hrSim}}, \code{\link{hrSimParallel}}, \code{\link{hrSimNew}}
}
\examples{
  \dontrun{
    require( parallel )
    
    # on Windows you must use socket clusters
    # Linux/UNIX supports other types of clusters
    # 
    # Change '4' to reflect the number of 
    # cores/processors you want to use
    thecluster <- makePSOCKcluster(4)
  
    # initialize each node
    tmp.rv <- clusterEvalQ( cl = thecluster, {
      require( CerioliOutlierDetection )
      require( HardinRockeExtensionSimulations )
      N.SIM <- 5000
      B.SIM <- 500

      my.pid <- Sys.getpid()
      cat("My pid is ", my.pid, "\n")
      logfile <- paste("Test_New_Method_Parallel_logfile_",my.pid,".txt",sep="")
      cat("Initialized\n\n", file=logfile)

      invisible(NULL)
    })
    
    hr.cases <- data.frame(t(as.matrix(expand.grid(list(
      p=c(5,10,20),n=c(50,100,500,1000),mcd.alpha=c(0.60,0.75,0.95))))))
    
    hrResults <- lapply(hr.cases, function(pn,clst,ns,bs) {
      cat("Trial p = ",pn[1]," n = ",pn[2],"\n")
      hrSimNewParallel(cl=clst, p = pn[1] , n = pn[2], 
          mcd.alpha=pn[3], N=ns, B=bs, lgf=logfile)
    }, clst=thecluster, ns=N.SIM, bs=B.SIM)

    stopCluster(thecluster)
    
    # calculate column means for each result
    allmeans <- as.data.frame(t(rbind( hr.cases, 
      100*sapply( hrResults, function(x) colMeans(x) )
    )))
    row.names(allmeans) <- NULL
    
    # calculate column stdevs for each result
    allstds <- as.data.frame(t(rbind( hr.cases, 
      100*sapply( hrResults, function(x) apply(x,2,sd) )
    )))
    row.names(allstds) <- NULL
    
    # format for easier comparision to hardin and rocke paper
    reshape(allmeans[,c("p","n","CHI2.CON")  ], 
      direction="wide", idvar="p",timevar="n")
    reshape(allmeans[,c("p","n","CGASY.CON") ], 
      direction="wide", idvar="p",timevar="n")
    reshape(allmeans[,c("p","n","CGPRED.CON")], 
      direction="wide", idvar="p",timevar="n")
    
    reshape(allstds[,c("p","n","CHI2.CON")   ], 
      direction="wide", idvar="p",timevar="n")
    reshape(allstds[,c("p","n","CGASY.CON")  ], 
      direction="wide", idvar="p",timevar="n")
    reshape(allstds[,c("p","n","CGPRED.CON") ], 
      direction="wide", idvar="p",timevar="n")
  }
}
\keyword{ robust }
\keyword{ multivariate }
