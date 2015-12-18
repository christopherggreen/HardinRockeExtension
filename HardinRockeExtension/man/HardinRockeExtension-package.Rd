\name{HardinRockeExtension-package}
\alias{HardinRockeExtension-package}
\alias{HardinRockeExtension}
\docType{package}
\title{
  Package for replicating the results of Hardin and Rocke (2005)
  and Green and Martin (2014).
}
\description{
  Package for replicating the results of Hardin and Rocke (2005)
  and Green and Martin (2014).
}
\details{
\tabular{ll}{
Package: \tab HardinRockeExtension\cr
Type: \tab Package\cr
Version: \tab 1.0.13\cr
Date: \tab 2015-12-17\cr
License: \tab GPL-2\cr
}
The \code{hrSim} and \code{hrSimParallel} functions use the
Hardin and Rocke approximation for the distribution of the MCD
to 
}
\author{
Written and maintained by Christopher G. Green <christopher.g.green@gmail.com>
}
\references{
C. G. Green and R. Douglas Martin. An extension of a method of Hardin and Rocke, with
  an application to multivariate outlier detection via the IRMCD method of Cerioli.
  Working Paper, 2014. Available from 
  \url{http://students.washington.edu/cggreen/uwstat/papers/cerioli_extension.pdf}

J. Hardin and D. M. Rocke. The distribution of robust distances. 
Journal of Computational and Graphical Statistics, 14:928-946, 2005.
}
\keyword{ package }
\keyword{ robust }
\keyword{ multivariate }
\seealso{
  \code{\link[CerioliOutlierDetection:CerioliOutlierDetection-package]{CerioliOutlierDetection}} 
}
\examples{
  \dontrun{
      hrSim(p=5,n=50  ,N=5000,B=500)
  }
}
