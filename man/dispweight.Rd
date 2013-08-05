\encoding{UTF-8}
\name{dispweight}
\alias{dispweight}
\title{Dispersion-based weighting of species counts}

\description{Transform abundance data using dispersion-based weighting of species counts}

\usage{
  dispweight(comm, group, nperm = 1000)
  dispweight2(comm, groups, nperm = 999, nullmodel = "c0_ind2",
      plimit = 0.05, ...)
  dispweight3(formula, data, plimit = 0.05)
}


\arguments{
  \item{comm}{Community data matrix.}
  \item{group, groups}{Factor describing the group structure.}
  \item{formula, data}{Formula where the left-hand side is the
    community data frame and right-hand side gives the explanatory
    variables. The explanatory variables are found in the data frame
    given in \code{data}.}
  \item{nperm}{Number of permutations.}
  \item{nullmodel}{The \code{\link{nullmodel}} used in
    \code{\link{make.commsim}} within \code{groups}. The default
    follows Clarke et al. (2006). }
  \item{plimit}{Downweight species if their \eqn{p}-value is at or
     below this limit. } 
}

\details{
Dispersion weighting can be used when there are \emph{a priori} defined groups 
(eg. location of samples). 

The dispersion index (\eqn{D}) is calculated as ratio between 
variance and mean and averaged across groups. D can then be used as divisor to 
downweight species abundances.

However species should be only downweighted, when there is evidence for 
over-dispersion (\eqn{D > 1}). This is tested using a permutational Chi-squared test. 
If this test shows statistically significant over-dispersion then the species is 
downweighted by \eqn{1/D}, else the divisor is \eqn{1}.

Functions \code{dispweight} and \code{dispweight2} can estimate the
species means and variances within levels of one factor (\code{group},
\code{groups}), and they estimate the significance of overdispersion
using a \code{\link{nullmodel}} simulation within the levels of this
factor. Function \code{dispweight3} finds the species means and their
residual variances with respect to any \code{\link{quasipoisson}}
\code{\link{glm}} model. Function \code{dispweight3} does not use
simulation, but it assesses the significance of overdispersion based
on the approximate Chi-square distribution of Pearson residuals.
}

\value{
A list of class \code{dispweight} with the following items:
    \item{D}{Average Dispersion index.}
    \item{p}{\eqn{p}-value of permutation test that \eqn{D = 1}.}
    \item{weights}{weights applied to community data.}
    \item{transformed}{transformed community data.}
}

\references{
Clarke, K. R., M. G. Chapman, P. J. Somerfield, and H. R. Needham. 2006. Dispersion-based Weighting of Species Counts in Assemblage Analyses. \emph{Marine Ecology Progress Series}, 320, 11–27.
}

\author{ 
  Eduard Szöcs \email{szoe8822@uni-landau.de} and Jari Oksanen
  (\code{dispweight2}).
}

\examples{
data(dune)
data(dune.env)
# calculate weights
dpw <- dispweight(dune, dune.env$Management, nperm = 100)
# transformed community data
dpw$transformed
}


\keyword{multivariate}

