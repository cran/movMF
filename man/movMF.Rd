\name{movMF}
\alias{movMF}
\title{Fit Mixtures of von Mises-Fisher Distributions}
\description{
  Fit mixtures of von Mises-Fisher distributions.
}
\usage{
movMF(x, k, control = list(), ...)
}
\arguments{
  \item{x}{a numeric data matrix, with rows corresponding to
    observations.  Standardized to unit row lengths if necessary.
    Can be a dense matrix, a
    \link[slam:matrix]{simple triplet matrix}
    (package \pkg{slam}), or a
    \link[Matrix:dgTMatrix-class]{dgTMatrix}
    (package \pkg{Matrix}).}
  \item{k}{an integer giving the desired number of mixture components
    (classes).}
  \item{control}{a list of control parameters.  See \bold{Details}.}
  \item{\dots}{a list of control parameters (overriding those specified
    in \code{control}).}
}
\value{
  An object of class \code{"movMF"} representing the fitted mixture of
  von Mises-Fisher distributions, which is a list containing at least
  the following components:
  \item{theta}{a matrix with rows giving the fitted parameters of the
    mixture components.}
  \item{alpha}{a numeric vector with the fitted mixture probabilities.}
  See \link{vMF} for the employed parametrization of the von
  Mises-Fisher distribution.
}
\details{
  \code{movMF} returns an object of class \code{"movMF"} representing
  the fitted mixture of von Mises-Fisher distributions model.  Available
  methods for such objects include \code{\link{coef}},
  \code{\link{logLik}}, \code{\link{print}} and \code{\link{predict}}.
  \code{\link{predict}} has an extra \code{type} argument with possible
  values \code{"class_ids"} (default) and \code{"memberships"} for
  indicating hard or soft prediction, respectively.
  
  The mixture of von Mises-Fisher distributions is fitted using EM
  variants as specified by control option \code{E} (specifying the
  E-step employed), with possible values \code{"softmax"} (default),
  \code{"hardmax"} or \code{"stochmax"} where the first two implement
  algorithms soft-moVMF and hard-moVMF of Banerjee et al (2005).  For
  \code{"stochmax"}, class assignments are drawn from the posteriors for
  each observation in the E-step as outlined as SEM in Celeux and
  Govaert (1992). The stopping criterion for this algorithm is by
  default changed to not check for convergence (logical control option
  \code{converge}), but to return the parameters with the maximum
  likelihood encountered. \code{E} may be abbreviated.

  In the M-step, the parameters \eqn{\theta} of the respective component
  distributions are estimated via maximum likelihood, which is
  accomplished by taking \eqn{\theta} proportional to suitable weighted
  sample means \eqn{\bar{x}}{xbar}, with length \eqn{\kappa} solving the
  equation \eqn{A_d(\kappa) = \|\bar{x}\|}{A_d(kappa) = ||xbar||}, where
  \eqn{A_d(\kappa) = I_{d/2}(\kappa) / I_{d/2-1}(\kappa)} with \eqn{I}
  the modified Bessel function of the first kind.  Via control argument
  \code{kappa}, one can specify how to (approximately) solve these
  equations, and whether a common (possibly given) length \eqn{\kappa}
  should be employed.  If \code{kappa} is a number, it gives a common
  length to be employed.  If it is a character string, it specified the
  method to be used for solving the \eqn{\kappa} equation.  Method
  \code{"Banerjee_et_al_2005"} uses the approximation of Banerjee et al
  (2005), \code{"Tanabe_et_al_2007"} the fixed-point iteration of Tanabe
  et al (2007), \code{"Sra_2012"} uses two Newton steps as suggested in
  Sra (2012), \code{"uniroot"} a straightforward call to
  \code{\link{uniroot}}, and finally \code{"Newton"} (default) uses a
  full Newton algorithm.  The lower-cased version of the given
  \code{kappa} specification is matched against the lower-cased names of
  the available methods using \code{\link{pmatch}}.  Finally, to
  indicate using a common (but not given) \eqn{kappa} for all component
  distributions, \code{kappa} should be a list with element \code{common
  = TRUE} (and optionally a character string giving the estimation
  method).

  Additional control parameters are as follows.

  \describe{
    \item{\code{maxiter}}{an integer giving the maximal number of EM
      iterations to be performed.  Default: 100.
    }
    \item{\code{reltol}}{The minimum relative improvement of the
      objective function per iteration. If improvement is less, the EM
      algorithm will stop under the assumption that no further
      significant improvement can be made.  Defaults to
      \code{sqrt(.Machine$double.eps)}.
    }
    \item{\code{ids}}{either a vector of class memberships or
      \code{TRUE} which implies that the class memberships are obtained
      from the attribute named \code{"z"} of the input data; these class
      memberships are used for initializing the EM algorithm and the
      algorithm is stopped after the first iteration.}
    \item{\code{start}}{a specification of the starting values to be
      employed.  Can be a list of matrices giving the memberships of
      objects to components, or of vectors giving component ids
      (numbers from 1 to the given \code{k}).  Can also be a character
      vector with elements \code{"i"} (randomly pick component ids for
      the observations), or one of \code{"p"}, \code{"S"} or \code{"s"}.
      The latter first determine component \dQuote{prototypes}, and then
      determine an optimal \dQuote{fuzzy} membership matrix from the
      implied cosine dissimilarities between observations and
      prototypes.  Prototypes are obtained as follows: for \code{"p"},
      observations are randomly picked.  For \code{"S"}, one takes the
      first prototype to minimize total cosine dissimilarity to the
      observations, and then successively picks observations farthest
      away from the already picked prototypes.  For \code{"s"}, one
      takes a randomly chosen observation as the first prototype, and
      then proceeds as for \code{"S"}. 

      By default, initialization method \code{"p"} is used.

      If several starting values are specified, the EM algorithm is
      performed individually to each starting value, and the best
      solution found is returned.
    }
    \item{\code{nruns}}{an integer giving the number of EM runs to be
      performed.  Default: 1. 
      Only used if \code{start} is not given.
    }
    \item{\code{minalpha}}{a numeric indicating the minimum prior
      probability.  Components falling below this threshold are removed
      during the iteration.  If \eqn{\ge 1}, this is taken as the
      minimal number of observations in a component.
      Default: 0.}

    \item{\code{converge}}{a logical, if \code{TRUE} the EM algorithm is
      stopped if the \code{reltol} criterion is met and the current
      parameter estimate is returned. If \code{FALSE} the EM algorithm
      is run for \code{maxiter} iterations and the parametrizations
      with the maximum likelihood encountered during the EM algorithm is
      returned. Default: \code{TRUE}, changed to \code{FALSE} if
      \code{E="stochmax"}.}

    \item{\code{verbose}}{a logical indicating whether to provide
      some output on algorithmic progress.
      Defaults to \code{getOption("verbose")}.}
  }

  One popular application context of mixtures of von Mises-Fisher
  distributions is text mining, where the data matrices are typically
  very large and sparse.  The provided implementation should be able to
  handle such large corpora with reasonable efficiency by employing
  suitable sparse matrix representations and computations.  In addition,
  straightforward computations of the normalizing constants in the vMF
  densities (see \link[=dmovMF]{movMF_distribution}) by directly
  employing the modified Bessel functions of the first kind are
  computationally infeasible for large \eqn{d} (dimension of the
  observations) and/or values of the concentration parameters.  Instead,
  we use suitably scaled hypergeometric-type power series for computing
  (the logarithms of) the normalizing constants.
}
\references{
  A. Banerjee, I. S. Dhillon, J. Ghosh, and S. Sra (2005).
  Clustering on the unit hypersphere using von Mises-Fisher
  distributions.
  \emph{Journal of Machine Learning Research}, \bold{6}, 1345--1382.
  \url{http://jmlr.csail.mit.edu/papers/v6/banerjee05a.html}.

  G. Celeux, and G. Govaert (1992).
  A Classification EM algorithm for clustering and two stochastic
  versions.
  \emph{Computational Statistics & Data Analysis}, \bold{14}, 315--332.
  doi: 10.1016/0167-9473(92)90042-E.
   
  S. Sra (2012).
  A short note on parameter approximation for von Mises-Fisher
  distributions: and a fast implementation of \eqn{I_s(x)}.
  \emph{Computational Statistics}, \bold{27}, 177--190.
  doi: 10.1007/s00180-011-0232-x.

  A. Tanabe, K. Fukumizu, S. Oba, T. Takenouchi, and S. Ishii.
  Parameter estimation for von Mises-Fisher distributions.
  \emph{Computational Statistics}, \bold{22}, 145--157.
  doi: 10.1007/s00180-007-0030-7.  
}
\examples{
## Generate and fit a "small-mix" data set a la Banerjee et al.
mu <- rbind(c(-0.251, -0.968),
            c(0.399, 0.917))
kappa <- c(4, 4)
theta <- kappa * mu
theta
alpha <- c(0.48, 0.52)
## Generate a sample of size n = 50 from the vMF mixture with the above
## parameters.
set.seed(123)
x <- rmovMF(50, theta, alpha)
## Fit a vMF mixture with the "right" number of components, using 10 EM
## runs.
y2 <- movMF(x, 2, nruns = 10)
## Inspect the fitted parameters:
y2
## Compare the fitted classes to the true ones:
table(True = attr(x, "z"), Fitted = predict(y2))
## To use a common kappa:
y2cv <- movMF(x, 2, nruns = 10, kappa = list(common = TRUE))
## To use a common kappa fixed to the true value of 4:
y2cf <- movMF(x, 2, nruns = 10, kappa = 4)
## Comparing solutions via BIC:
sapply(list(y2, y2cf, y2cv), BIC)
}
