\name{NEWS}
\title{News for Package \pkg{movMF}}
\encoding{UTF-8}
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
