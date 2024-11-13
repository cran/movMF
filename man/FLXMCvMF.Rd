\name{FLXMCvMF}
\alias{FLXMCvMF}
\title{Flexmix Driver for Mixtures of von Mises-Fisher Distributions}
\description{
  This driver for \code{\link[flexmix]{flexmix}} implements
  estimation of mixtures of von Mises-Fisher distributions where the
  data can be stored in a dense or a \link[slam]{simple triplet matrix}
  (package \pkg{slam}) format.
}
\usage{
FLXMCvMF(formula = . ~ ., kappa = NULL)
}
\arguments{
  \item{formula}{a formula which is interpreted relative to the formula
    specified in the call to \code{\link[flexmix]{flexmix}} using
    \code{\link[stats]{update.formula}}.  Only the
    left-hand side (response) of the formula is used. Default is to use
    the original \code{\link[flexmix]{flexmix}} model
    formula.}
  \item{kappa}{see the \code{control} argument of \code{\link{movMF}}}.
}
\value{
  An object of class \code{"FLXMCvMF"}.
}
\author{
  Bettina Grün
}
\examples{
if (requireNamespace("flexmix", quietly = TRUE)) {
## Generate and fit a "small-mix" data set a la Banerjee et al.
mu <- rbind(c(-0.251, -0.968),
            c(0.399, 0.917))
kappa <- c(4, 4)
theta <- kappa * mu
theta
alpha <- c(0.48, 0.52)
## Generate a sample of size n = 50 from the von Mises-Fisher mixture
## with the above parameters.
set.seed(123)
x <- rmovMF(50, theta, alpha)
## Fit a von Mises-Fisher mixture with the "right" number of components,
## using 10 EM runs.
set.seed(123)
y2 <- flexmix::stepFlexmix(x ~ 1, k = 2, model = FLXMCvMF(), verbose = FALSE)
## Inspect the fitted parameters:
y2
## Compare the fitted classes to the true ones:
table(True = attr(x, "z"), Fitted = flexmix::clusters(y2))
## To use a common kappa:
y2cv <- flexmix::stepFlexmix(x ~ 1, k = 2,
  model = FLXMCvMF(kappa = list(common = TRUE)), verbose = FALSE)
## To use a common kappa fixed to the true value of 4:
y2cf <- flexmix::stepFlexmix(x ~ 1, k = 2,
  model = FLXMCvMF(kappa = 4), verbose = FALSE)
## Comparing solutions via BIC:
sapply(list(y2, y2cf, y2cv), BIC)
##  Use a different kappa solver:
set.seed(123)
y2a <- flexmix::stepFlexmix(x ~ 1, k = 2,
  model = FLXMCvMF(kappa = "uniroot"), verbose = FALSE)
y2a
## Using a sparse matrix:
x <- slam::as.simple_triplet_matrix(x)
y2 <- flexmix::stepFlexmix(x ~ 1, k = 2,
  model = FLXMCvMF(), verbose = FALSE)
}
}
\keyword{cluster}
