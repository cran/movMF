\name{NEWS}
\title{News for Package \pkg{movMF}}
\encoding{UTF-8}
\section{Changes in version 0.2-8}{
  \itemize{
    \item \pkg{corpus.useR.2008.abstracts} was removed as suggested
    package and the dataset added to the package to ensure that the
    vignette code can be successfully executed even if the package is
    not available for installation.
  }
}
\section{Changes in version 0.2-7}{
  \itemize{
    \item Performance tweaks for random number generation from the vMF
    distribution.
  }
}
\section{Changes in version 0.2-6}{
  \itemize{
    \item Update URLs.
  }
}
\section{Changes in version 0.2-5}{
  \itemize{
    \item Model driver added to extend package \pkg{flexmix}.
    \item Package \pkg{HSAUR3} now used instead of \pkg{HSAUR2} for
    dataset \code{household}.
  }
}

\section{Changes in version 0.2-4}{
  \itemize{
    \item Minor vignette improvements.
  }
}

\section{Changes in version 0.2-3}{
  \itemize{
    \item A bug in function \code{solve_kappa_uniroot} is now corrected.
  }
}

\section{Changes in version 0.2-2}{
  \itemize{
    \item Minor vignette improvements.
  }
}

\section{Changes in version 0.2-1}{
  \itemize{
    \item The information on publication of the vignette in JSS is now
    also mentioned in the vignette.
  }
}

\section{Changes in version 0.2-0}{
  \itemize{
    \item The package vignette has been published in the Journal of
    Statistical Software, Volume 58, Issue 10
    (\doi{10.18637/jss.v058.i10}), and the paper should be
    used as citation for the package, run \code{citation("movMF")} for
    details.
    \item In the calculation of 0F1 it is now ensured that summation is
    continued as long as the terms are increasing and only then
    numerical convergence is checked as stopping criterion.
    \item For simulating from the vMF distribution the necesssary
    quantities in the rejection sampling are now calculated in a
    numerically more stable way. Thanks to Mike Allerhand for pointing
    out that problems did arise for large values of the precision
    parameter for the previous implementation.
    \item Package \pkg{corpus.useR.2008.abstracts} is now listed in the
    DESCRIPTION file together with the information that is available from
    the additional repository \url{http://datacube.wu.ac.at}.
  }
}

\section{Changes in version 0.1-3}{
  \itemize{
    \item An improved approximation is used to calculate the normalizing
    constants of the log-likelihoods.  See the package vignette for more
    information.
  }
}

\section{Changes in version 0.1-2}{
  \itemize{
    \item For determining the concentration parameters the following
    additonal methods were implemented: \code{"Halley"}, \code{"hybrid"}
    and \code{"Newton_Fourier"}. The default method is now
    \code{"Newton_Fourier"}. 
  }
}

\section{Changes in version 0.1-1}{
  \itemize{
    \item Package \pkg{SnowballC} is now suggested instead of
    \pkg{Snowball}.
  }
}

\section{Changes in version 0.1-0}{
  \itemize{
    \item Package now uses registration for native (C) routines.
    \item For determining the concentration parameters the additional
    method \code{"Song_et_al_2012"} was added which performs 2 Halley
    steps. 
    \item A vignette documenting the use of the package and some
    implementation issues is added.
    \item Function \code{A()} and the functions of its first and second
    derivative now use the Taylor series expansion for small values of
    kappa. 
  }
}

\section{Changes in version 0.0-2}{
  \itemize{
    \item Function \code{A()} modified to have now three methods
    (\code{"PCF"}: Perron continued fraction, \code{"GCF"}: Gauss
    continued fraction, \code{"RH"}: using a ratio of generalized
    hypergeometric series).  By default, method \code{"PCF"} is used.
  }
}

\section{Changes in version 0.0-1}{
  \itemize{
    \item The random initializations for multiple runs are now
    determined within the loop separately for each run. 
  }
}
