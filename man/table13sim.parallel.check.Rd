\name{table13sim.parallel.check}
\alias{table13sim.parallel.check}
\title{Replicate part of the experiment presented in Cerioli et al. (2009)}
\description{
  Replicates a small part the experiment presented in Cerioli et al. (2009), 
  Tables 1 and 3, for the MCD using the maximum breakdown point
  fraction and a fraction of exactly 0.5.
}
\usage{
table13sim.parallel.check(cl, p, nn, N, B = 250, alpha = c(0.01, 0.025, 0.05), lgf = "", mlgf = "")
}
\arguments{
  \item{cl}{A cluster object, e.g., returned from \code{makePSOCKcluster}.
    The user must create this object before calling \code{table13sim.parallel}.}
  \item{p}{The dimension of the data used in each simulated run.}
  \item{nn}{The number of observations used in each simulated run.}
  \item{N}{The number of simulations to run.}
  \item{B}{The batch/block size: the number of simulations to run 
    in each block. This is useful when running very large
    simulation runs (\code{N} very large) where memory is a concern. 
  Blocks are distributed across cluster nodes, so this is also a means
  of controlling the workload on each node.}
  \item{alpha}{The significance level to use for detecting outliers. Can
    be a vector; the outlier detection tests will be run at each level.}
  \item{lgf}{Path to log file into which logging information should be
    written.}
  \item{mlgf}{ not used at this time }
}
\details{
  This is a variant of \code{\link{table13sim.parallel}} designed
  to investigate differences between the outlier detection tests with
  the MCD when the data fraction is (a) the maximum breakdown point and
  (b) exactly 0.5. It also checks whether the small-sample correction
  was used in the results of Hardin and Rocke (2005) and Cerioli et 
  al. (2009).

  This function is not really useful to anyone other than the author,
  and is not supported. \emph{Do not use it.}
}
\value{An array of dimension 3:
  \enumerate{
    \item{The results of each of the \code{N} simulation runs appear along the first dimension.}
  \item{The various estimators and tests appear along the second dimension. Results with
    suffix ``T1'' correspond to Table 1 of Cerioli et al. (2009) (the individual outlier tests)
    while those with suffix ``T3'' correspond to Table 3 (the simultaneous outlier tests).
    Currently the 26 columns appear in the following order.
    \tabular{lll}{
      Column Name \tab Covariate Estimate \tab Test Statistic \cr
        "MCDMBP.RAW.T1" \tab MCD (max. breadown pt.) \tab chi-squared \cr
        "MCDMBP.RAWGM.T1" \tab MCD (max. breadown pt.) \tab Green-Martin \cr
        "MCDMBP.RAWHR.T1" \tab MCD (max. breadown pt.) \tab Hardin-Rocke \cr
        "MCDMBP.RAWNOSSGM.T1" \tab MCD (max. breadown pt.), no small sample correction \tab Green-Martin \cr
        "MCDMBP.RAWNOSSHR.T1" \tab MCD (max. breadown pt.), no small sample correction \tab Hardin-Rocke \cr
        "RMCDMBP.T1" \tab reweighted MCD (max. breadown pt.) \tab chi-squared \cr
        "MCDMBP.RAW.T3" \tab MCD (max. breadown pt.) \tab chi-squared \cr
        "MCDMBP.RAWGM.T3" \tab MCD (max. breadown pt.) \tab Green-Martin \cr
        "MCDMBP.RAWHR.T3" \tab MCD (max. breadown pt.) \tab Hardin-Rocke \cr
        "MCDMBP.RAWNOSSGM.T3" \tab MCD (max. breadown pt.), no small sample correction \tab Green-Martin \cr
        "MCDMBP.RAWNOSSHR.T3" \tab MCD (max. breadown pt.), no small sample correction \tab Hardin-Rocke \cr
        "RMCDMBP.T3" \tab reweighted MCD (max. breadown pt.) \tab chi-squared \cr
        "RMCDMBP.CH.T3" \tab reweighted MCD (max. breadown pt.) with Bonferroni correction \tab chi-squared \cr
        "MCD50.RAW.T1" \tab MCD(0.50) \tab chi-squared \cr
        "MCD50.RAWGM.T1" \tab MCD(0.50) \tab Green-Martin \cr
        "MCD50.RAWHR.T1" \tab MCD(0.50) \tab Hardin-Rocke \cr
        "MCD50.RAWNOSSGM.T1" \tab MCD(0.50), no small sample correction \tab Green-Martin \cr
        "MCD50.RAWNOSSHR.T1" \tab MCD(0.50), no small sample correction \tab Hardin-Rocke \cr
        "RMCD50.T1" \tab reweighted MCD(0.50) \tab chi-squared \cr
        "MCD50.RAW.T3" \tab MCD(0.50) \tab chi-squared \cr
        "MCD50.RAWGM.T3" \tab MCD(0.50) \tab Green-Martin \cr
        "MCD50.RAWHR.T3" \tab MCD(0.50) \tab Hardin-Rocke \cr
        "MCD50.RAWNOSSGM.T3" \tab MCD(0.50), no small sample correction \tab Green-Martin \cr
        "MCD50.RAWNOSSHR.T3" \tab MCD(0.50), no small sample correction \tab Hardin-Rocke \cr
        "RMCD50.T3" \tab reweighted MCD(0.50) \tab chi-squared \cr
        "RMCD50.CH.T3" \tab reweighted MCD(0.50) with Bonferroni correction \tab chi-squared \cr
      }
  } % end item
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
  \code{\link{table1sim.parallel}}
  \code{\link{CovMcd2}}
}
\examples{
  \dontrun{
    #
  }
} % examples
\keyword{ robust }
\keyword{ multivariate }
