\name{movMF_distribution}
\alias{dmovMF}
\alias{rmovMF}
\alias{vMF}
\title{Mixtures of von Mises-Fisher Distributions}
\description{
  Density and random number generation for finite mixtures of von
  Mises-Fisher distributions.
}
\usage{
dmovMF(x, theta, alpha = 1, log = FALSE)
rmovMF(n, theta, alpha = 1)
}
\arguments{
  \item{x}{a matrix of rows of points on the unit hypersphere.
    Standardized to unit row length if necessary.}
  \item{theta}{a matrix with rows giving the parameters of the mixture
    components.}
  \item{alpha}{a numeric vector with non-negative elements giving the
    mixture probabilities.
    Standardized to sum to one if necessary.}
  \item{log}{a logical; if \code{TRUE} log-densities are computed.}
  \item{n}{an integer giving the number of samples to draw.}  
}
\value{
  For \code{dmovMF}, a numeric vector of (log) density values.

  For \code{rmovMF}, a matrix with \code{n} unit length rows
  representing the samples from the vMF mixture distribution.
}
\details{
  A random \eqn{d}-dimensional unit length vector \eqn{x} has a von
  Mises-Fisher (or Langevin, short: vMF) distribution with parameter
  \eqn{\theta} if its density with respect to the uniform distribution
  on the unit hypersphere is given by
  \deqn{f(x|\theta) = \exp(\theta'x) / {}_0F_1(; d/2; \|\theta\|^2/4),}{
    f(x|theta) = exp(theta'x) / 0F1(; d / 2; ||theta||^2 / 4),
  }
  where \eqn{{}_0F_1}{0F1} is a generalized hypergeometric function
  (e.g.,
  \url{http://en.wikipedia.org/wiki/Generalized_hypergeometric_function})
  and related to the modified Bessel function \eqn{I_\nu} of the first
  kind via
  \deqn{{}_0F_1(; \nu+1; z^2/4) =
    I_\nu(z)\Gamma(\nu+1) / (z/2)^\nu.}{
    0F1(; nu+1; z^2 / 4) = I_nu(z) Gamma(nu + 1) / (z/2)^nu.
  }
  With this parametrization, the von Mises-Fisher family is the natural
  exponential family through the uniform distribution on the unit
  sphere, with cumulant transform
  \deqn{M(\theta) = \log({}_0F_1(; d/2; \|\theta\|^2/4)).}{
    M(theta) = log(0F1(; d / 2; ||theta||^2 / 4)).}

  We note that the vMF distribution is commonly parametrized by the
  \emph{mean direction parameter} \eqn{\mu = \theta / \|\theta\|}{mu =
  theta / ||theta||} (which however is not well-defined if \eqn{\theta =
  0}) and the \emph{concentration parameter} \eqn{\kappa = \|\theta\|}{
  kappa = ||theta||}, e.g., 
  \url{http://en.wikipedia.org/wiki/Von_Mises\%E2\%80\%93Fisher_distribution}
  (which also uses the un-normalized Haar measure on the unit sphere as
  the reference distribution, and hence includes the \dQuote{area} of
  the unit sphere as an additional normalizing constant).

  \code{dmovMF} computes the (log) density of mixtures of vMF
  distributions.

  \code{rmovMF} generates samples from finite mixtures of vMF
  distributions, using Algorithm VM* in Wood (1994) for sampling from
  the vMF distribution.

  Arguments \code{theta} and \code{alpha} are recycled to a common
  number of mixture components.
}
\references{
  A. T. A. Wood (1994).
  Simulation of the von Mises Fisher distribution.
  \emph{Communications in Statistics -- Simulation and Computation},
  \bold{23}(1), 157--164.
}
\examples{
## To simulate from the vMF distribution with mean direction
## proportional to c(1, -1) and concentration parameter 3:
rmovMF(10, 3 * c(1, -1) / sqrt(2))
## To simulate from a mixture of vMF distributions with mean direction
## parameters c(1, 0) and c(0, 1), concentration parameters 3 and 4, and
## mixture probabilities 1/3 and 2/3, respectively:
rmovMF(10, c(3, 4) * rbind(c(1, 0), c(0, 1)), c(1, 2))
}
