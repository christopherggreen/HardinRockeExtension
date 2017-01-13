\name{CovMcd2}
\alias{CovMcd2}
\title{MCD Covariance Estimate}
\description{
  This is a slight modification of the \code{CovMcd} function
  from the \code{rrcov} package for use with Cerioli et al. (2009)
  replication and extension work presented in Green and Martin (2014). 
  The modifications provide more flexibility in computing the 
  reweighted MCD estimate.
}
\usage{
CovMcd(x,
       raw.only=FALSE, alpha=control@alpha, nsamp=control@nsamp,
       scalefn=control@scalefn, maxcsteps=control@maxcsteps,
       initHsets=NULL, save.hsets=FALSE,
       seed=control@seed, trace=control@trace,
       use.correction=control@use.correction,
       control=CovControlMcd(), ...)
}
\arguments{
  \item{x}{See help file for \code{\link[rrcov]{CovMcd}}.}
  \item{raw.only}{See help file for \code{\link[rrcov]{CovMcd}}.}
  \item{alpha}{See help file for \code{\link[rrcov]{CovMcd}}.}
  \item{nsamp}{See help file for \code{\link[rrcov]{CovMcd}}.}
  \item{scalefn}{See help file for \code{\link[rrcov]{CovMcd}}.}
  \item{maxcsteps}{See help file for \code{\link[rrcov]{CovMcd}}.}
  \item{initHsets}{See help file for \code{\link[rrcov]{CovMcd}}.}
  \item{save.hsets}{See help file for \code{\link[rrcov]{CovMcd}}.}
  \item{seed}{See help file for \code{\link[rrcov]{CovMcd}}.}
  \item{trace}{See help file for \code{\link[rrcov]{CovMcd}}.}
  \item{use.correction}{See help file for \code{\link[rrcov]{CovMcd}}.}
  \item{control}{See help file for \code{\link[rrcov]{CovMcd}}.}
  \item{\dots}{See help file for \code{\link[rrcov]{CovMcd}}.}
  \item{reweighting}{A list describing how to perform the reweighting
    when computing the reweighted MCD estimate. The list should have
    the following components.
    \describe{
      \item{method: }{Method to score observations in the reweighting.
        Must be a string from the set (\code{"chisquare","empirical"}). 
        \code{"chisquare"} uses a quantile from a chi-squared distribution,
        while \code{"empirical"} will use a quantile from the empirical
        distribution of the Mahalanobis distances calculated from the 
        unreweighted MCD estimate. Default is \code{"chisquare"}.}
      \item{quantile: }{The quantile to use in scoring. Default is 0.975
        for the chi-squared distribution, which matches how the 
      unmodified \code{covMcd} works.}
      \item{debug: }{TRUE to print additional debugging information; currently 
        prints the Mahalanobis distances (sorted) calculated from the 
        unreweighted MCD estimate. Defaults to FALSE.}
    }
  }
}
\details{
  This function is a modification of \code{CovMcd} from the \code{rrcov}
  package and uses a modified version of \code{covMcd} from the 
  \code{robustbase} package. It works identically to those functions but 
  allows for more flexibility in the reweighting portion of the calculation.
  This was needed to run/extend some of the experiments in Cerioli et al. (2009),
  particularly those presented in Table 3 of that paper, where the reweighting
  threshold is a Bonferroni-corrected quantile.

  The main work is done inside \code{covMcd2}, a private function to the
  \code{HardinRockeExtensionSimulations} package. This function is a modified version of
  \code{covMcd} from the \code{robustbase} package: it replaces the
  quantity \code{quantiel} used in the reweighting with a function of the same
  name that computes the quantile specified in the \code{reweighting}
  argument.
}
\value{
  An S4 object of class \code{\link{CovMcd-class}} which is a subclass of the
  virtual class \code{\link{CovRobust-class}}. See \code{\link[rrcov]{CovMcd}}
  for more details.
}
\references{
  Andrea Cerioli, Marco Riani, and Anthony C. Atkinson. Controlling the size of multivariate
  outlier tests with the mcd estimator of scatter. Statistical Computing, 19:341-353, 2009.

  C. G. Green and R. Douglas Martin. An extension of a method of Hardin and Rocke, with
  an application to multivariate outlier detection via the IRMCD method of Cerioli.
  Working Paper, 2014. Available from 
  \url{http://students.washington.edu/cggreen/uwstat/papers/cerioli_extension.pdf}

  Todorov V & Filzmoser P (2009),
  An Object Oriented Framework for Robust Multivariate Analysis.
  \emph{Journal of Statistical Software}, \bold{32}(3), 1--47.
  URL \url{http://www.jstatsoft.org/v32/i03/}.
}
\author{
  Written and maintained by Christopher G. Green <christopher.g.green@gmail.com>.
}
\seealso{
  \code{\link[rrcov]{CovMcd}}, \code{\link[robustbase]{covMcd}}
}
\note{
  This function is subject to change as the \code{CovMcd} and \code{covMcd}
  functions change. It may also be out of sync with those functions at 
  times.

  The latest version of \code{covMcd} (from \code{robustbase} version 0.92-5)
  has support for a custom weighting function, so at some point in the future
  the need to use a modified version of \code{covMcd} might be eliminated.
}
\examples{
   require(rrcov)
   require(mvtnorm)
   alpha   <- 0.05 # significance level
   nn      <- 100  # number of observations
   p       <- 4    # dimension
   # generate some data
   simdata <- mvtnorm::rmvnorm(nn, mean=rep(0,p), sigma=diag(rep(1,p)))
   CovMcd2(simdata,
      reweighting=list(
      method="chisquare",
      quantile=1-(alpha/nn),
      debug=FALSE
    )
   )
}
\keyword{ robust }
\keyword{ multivariate }
