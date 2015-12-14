\name{hrSim}
\alias{hrSim}
\title{
  Replicate Tables 1--3 of Hardin and Rocke (2005)
}
\description{
  \code{hrSim} is used to replicate Tables 1--3 of Hardin
  and Rocke (2005), page 942.
}
\usage{
hrSim(p, n, N, B = 10000, alpha = 0.05, mcd.alpha = max.bdp.mcd.alpha(n, p))
}
\arguments{
  \item{p}{The dimension of the data used in each simulated run.}
  \item{n}{The number of observations used in each simulated run.}
  \item{N}{The number of simulations to run.}
  \item{B}{The batch/block size: the number of simulations to run 
    in each block. This is useful when running very large
    simulation runs (\code{N} very large) where memory is a concern.}
  \item{alpha}{The significance level to use for detecting outliers.}
  \item{mcd.alpha}{The fraction of the data to use in computing the
    MCD. Defaults to the maximum breakdown point fraction.}
}
\details{
  This is a work function designed for use in replicating Tables 1--3 
  of Hardin and Rocke (2005), page 942. 

  Internally the simulation function does \code{B} runs at a time. Set 
  \code{B} smaller if your machine has less memory.
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
  author whether Hardin and Rocke (2005) used these correction factors,
  so all variants were calculated and examined. It seems the ``CON''
  approach is the best match for their results. CHECK THIS ASSERTION)

  Look at the column means of the resulting matrix to see the average 
  fraction of outliers detected (which is an estimate of the Type 1 error rate
  of the procedure, since the simulated data had no outliers).
}
\references{
  J. Hardin and D. M. Rocke. The distribution of robust distances. 
  Journal of Computational and Graphical Statistics, 14:928-946, 2005.
}
\author{
Written and maintained by Christopher G. Green <christopher.g.green@gmail.com>
}
\note{
  This function is deprecated; use \code{\link{hrSimParallel}} instead.
}
\seealso{
\code{\link{covMcd}}, \code{\link{hrSimParallel}}, \code{\link{hrSimNew}},
\code{\link{hrSimNewParallel}}
}
\examples{
  require( CerioliOutlierDetection )
  hrResults.p5.n50   <- hrSim(p=5,n=50  ,N=50,B=50)
}
\keyword{ robust }
\keyword{ multivariate }
