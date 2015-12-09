\name{hr.cm}
\alias{hr.cm}
\title{ Hardin and Rocke (2005) corrections for general MCD estimates }
\description{
This is a work function to calculate, via simulation, the consistency
correction ``c'' and the Wishart degrees of freedom ``m'' needed in
the Hardin and Rocke approximation for ???. This function generalizes
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
	we calculate the sum of the diagonal entries
	of said covariance matrix and the sum of the 
	squares of the diagonal entries. The estimate
	of ``c'', the consistency constant, is then
	given by the mean of the sums of the diagonal
	entries. The estimate of ``m'' uses this 
	estimated ``c'' and an estimate of the variance 
	of the diagonal elements.

	See Joanna Hardin's code cm.R and mcd.est.R.
}
\value{A data frame with the following columns
  \item{c}{The estimate of ``c'', the consistency constant}
  \item{m}{The estimated Wishart degrees of freedom ``m''}
  \item{mcd.alpha}{The input data fraction(s) \code{mcd.alpha}}
}
\references{
J. Hardin and D. M. Rocke. The distribution of robust distances. 
Journal of Computational and Graphical Statistics, 14:928-946, 2005.
}
\author{
Written and maintained by Christopher G. Green <christopher.g.green@gmail.com>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (p, n, N, B = 10000, mcd.alpha = max.bdp.mcd.alpha, 
    logfile = "hr_cm.log") 
{
    n.mcdalpha <- length(mcd.alpha)
    blockfcn <- function(m, b) {
        cat("Running block", m, "\n", file = logfile, append = TRUE)
        simdata <- array(rmvnorm(n * b, mean = rep(0, p), sigma = diag(rep(1, 
            p))), dim = c(n, p, b))
        results <- array(NA, dim = c(b, 2, n.mcdalpha))
        cat("\t Iteration: ", file = logfile, append = TRUE)
        for (j in 1:b) {
            if ((j\%\%50) == 0) 
                cat(j, " ")
            for (a in 1:n.mcdalpha) {
                mcd <- CovMcd(simdata[, , j], alpha = mcd.alpha[a])
                mcd.cov <- if (mcd.alpha[a] == 1) {
                  mcd@raw.cov
                }
                else {
                  cov(simdata[mcd@best, , j])
                }
                results[j, 1, a] <- sum(diag(mcd.cov))
                results[j, 2, a] <- sum(diag(mcd.cov)^2)
            }
        }
        cat("\nDone with block", m, "\n\n", file = logfile, append = TRUE)
        results
    }
    B <- min(N, B)
    M <- floor(N/B)
    blockresults <- lapply(1:M, blockfcn, b = B)
    r <- N - M * B
    if (r > 0) {
        blockresults[[M + 1]] <- blockfcn(m = M + 1, b = r)
    }
    blockresults <- do.call("abind", c(blockresults, list(along = 1)))
    dimnames(blockresults)[[2]] <- c("AVE", "SQ")
    dimnames(blockresults)[[3]] <- paste("a", gsub("\\.", "", 
        mcd.alpha), sep = "")
    mcd.ave <- colSums(blockresults[, "AVE", ])
    mcd.sq <- colSums(blockresults[, "SQ", ])
    mcd.c <- mcd.ave/(N * p)
    mcd.var <- (mcd.sq - (mcd.c * mcd.c * N * p))/(N * p - 1)
    mcd.m <- 2 * mcd.c * mcd.c/mcd.var
    data.frame(c = mcd.c, m = mcd.m, mcd.alpha = mcd.alpha)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
