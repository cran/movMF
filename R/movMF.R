## Try to estimate movMF via EM.

movMF <-
function(x, k, control = list(), ...)
{
    ## Be nice a la David.
    control <- c(control, list(...))

    if(missing(k)) k <- NULL
    ## Normalize data just in case.
    x <- skmeans:::row_normalize(x)

    n <- nrow(x)
    d <- ncol(x)
    s <- d / 2 - 1

    ## Control parameters.

    maxiter <- control$maxiter
    if(is.null(maxiter))
        maxiter <- 100L

    reltol <- control$reltol
    if(is.null(reltol))
        reltol <- sqrt(.Machine$double.eps)

    verbose <- control$verbose
    if(is.null(verbose))
        verbose <- getOption("verbose")

    E_methods <-
        c("softmax",
          "hardmax",
          "stochmax")
    E <- control$E
    if(is.null(E))
        E <- "softmax"
    else {
        pos <- pmatch(tolower(E), tolower(E_methods))
        if(is.na(pos))
            stop("Invalid E-step method.")
        E <- E_methods[pos]
    }
    do_P <- switch(EXPR = E,
                   "softmax" =
                   function(G, g) exp(G  - g),
                   "hardmax" =
                   function(G, g)
                   .posterior_from_ids(max.col(G), ncol(G)),
                   function(G, g)
                   t(apply(exp(G - g), 1,
                           function(prob) rmultinom(1, 1, prob))))

    kappa_solvers <-
        c("Banerjee_et_al_2005",
          "Tanabe_et_al_2007",
          "Sra_2012",
          "uniroot",
          "Newton")
    kappa <- control$kappa
    if(is.numeric(kappa)) {
        ## CASE A: Use a common given value of kappa.
        kappa_given <- kappa
        do_kappa <- function(norms, alpha)
            rep.int(kappa_given, length(alpha))
        df_kappa <- function(k) 0L
    } else {
        kappa <- as.list(kappa)
        ## Should a common value be used or not?
        pos <- match("common", names(kappa), nomatch = 0L)
        if(pos > 0L) {
            use_common_kappa <- identical(kappa[[pos]], TRUE)
            kappa <- kappa[-pos]
        } else {
            use_common_kappa <- FALSE
        }
        if(length(kappa)) {
            ## Solver specifications.
            kname <- kappa[[1L]]
            kargs <- kappa[-1L]
            pos <- pmatch(tolower(kname), tolower(kappa_solvers))
            if(is.na(pos))
                stop("Invalid kappa solver.")
            kappa <- kappa_solvers[pos]
            kfun <- get(sprintf("solve_kappa_%s", kappa))
            solve_kappa <- function(Rbar)
                do.call(kfun, c(list(Rbar, d), kargs))
        } else {
            ## Default solver.
            solve_kappa <- function(Rbar)
                solve_kappa_Newton(Rbar, d)
        }
        if(use_common_kappa) {
            do_kappa <- function(norms, alpha) 
                rep.int(solve_kappa(sum(norms) / n), length(alpha))
            df_kappa <- function(k) 1L
        } else {
            do_kappa <- function(norms, alpha)
                solve_kappa(norms / (n * alpha))
            df_kappa <- function(k) k
        }
    }

    ## Allow to specify "known" ids for fitting the respective movMF
    ## model by maximum likelihood.  This can be accomplished within our
    ## framework as follows:
    ## * Compute P as the corresponding membership matrix.
    ## * Perform one M step to estimate the parameters.
    ## * Perform one E step to compute the likelihood, but keep P.
    ids <- control$ids
    if(!is.null(ids)) {
        if(identical(ids, TRUE))
            ids <- attr(x, "z")         # Be nice for the case of data
                                        # simulated by rmovMF().
        P0 <- .posterior_from_ids(ids, k)
        do_P <- function(G, g) P0
        start  <- list(P0)
        maxiter <- 1L
        if(!is.null(control$nruns))
            warning("control argument nruns ignored because ids are specified")
        if(!is.null(control$start))
            warning("control argument start ignored because ids are specified")
    } else {
        ## Initialization.
        start <- control$start
        nruns <- control$nruns
        if(is.null(start)) {
            if(is.null(nruns))
                nruns <- 1L
            start <- as.list(rep.int("p", nruns))
        }
        else {
            if(!is.list(start))
                start <- list(start)
            if(!is.null(nruns))
                warning("control argument nruns ignored because start is specified")
        }
    }
    
    nruns <- length(start)

    minalpha <- control$minalpha
    if(is.null(minalpha)) minalpha <- 0
    if(minalpha >= 1)
        minalpha <- minalpha / n

    converge <- control$converge
    if(is.null(converge)) {
        converge <- if(E == "stochmax") FALSE else TRUE
    }
    
    L_opt <- -Inf
    opt_old <- opt <- NULL

    run <- 1L
    logLiks <- vector(length = nruns)
    
    if(verbose && (nruns > 1L))
        message(gettextf("Run: %d", run))

    repeat {
        G <- NULL
        P <- movMF_init(x, k, start[[run]])
        L_old <- -Inf
        iter <- 0L
        logLiks[run] <- tryCatch({
            while(iter < maxiter) {
            ## M step.
                alpha <- colMeans(P)
                while(any(alpha < minalpha)) {
                    if(verbose) 
                        message("*** Removing one component ***")
                    nok <- which.min(alpha)
                    P <- P[, -nok, drop = FALSE]
                    P <- do_P(P, log_row_sums(P))
                    alpha <- colMeans(P)
                    if(!is.null(G))
                        L_old <- sum(log_row_sums(G[, -nok, drop = FALSE]))
                }
                if(any(alpha == 0))
                    stop("Cannot handle empty components")
                M <- skmeans:::g_crossprod(P, x)
                norms <- skmeans:::row_norms(M)
                M <- M / norms
                ## If a cluster contains only identical observations,
                ## Rbar = 1.
                kappa <- do_kappa(norms, alpha)
            
                ## E step.
                G <- cadd(skmeans:::g_tcrossprod(x, kappa * M),
                          log(alpha) -  lH(kappa, s))
                g <- log_row_sums(G)
                L_new <- sum(g)
                if(verbose && (iter %% verbose == 0))
                    message(gettextf("Iteration: %d *** L: %g", iter, L_new))
                if(converge) {
                    if(abs(L_old - L_new) < reltol * (abs(L_old) + reltol)) {
                        L_old <- L_new     
                        break
                    }
                    L_old <- L_new     
                } else if(L_new > L_old) {
                    L_old <- L_new
                    opt_old <- .movMF_object(kappa * M, alpha, L_old, P, iter)
                }
            
                P <- do_P(G, g)
                iter <- iter + 1L
            }
            if(L_old > L_opt) {
                opt <- if(converge)
                    .movMF_object(kappa * M, alpha, L_old, P, iter)
                else opt_old
                L_opt <- L_old
            }
            L_old
        }, error = function(e) {
            if(verbose) {
                ## <NOTE>
                ## Reporting problems is a bit of a mess in cases the
                ## above explicitly called stop(), in which case the
                ## condition call is
                ##   doTryCatch(return(expr), name, parentenv, handler)
                ## Ideally, in these cases we would throw suitably
                ## classed conditions, and provide call/message methods
                ## for them.
                call <- conditionCall(e)
                msg <- conditionMessage(e)
                s <- if(!is.null(call) &&
                        (substring(s <- deparse(call)[1L], 1L, 10L) !=
                         "doTryCatch"))
                    sprintf("Error in %s: %s", s, msg)
                else
                    sprintf("Error: %s", msg)
                message(sprintf("EM algorithm did not converge:\n%s", s))
                ## </NOTE>
            }
            NA
        })
        
        if(run >= nruns) break
        
        run <- run + 1L
        if(verbose)
            message(gettextf("Run: %d", run))
    }

    ## Compute log-likelihood.
    if(is.null(opt))
        stop("the EM algorithm did not converge for any run")
    ll <- L(x, opt$theta, opt$alpha)
    ## Add the "degrees of freedom" (number of (estimated) parameters in
    ## the model): with k the number of classes actually used,
    ##   \mu: k * (d - 1) (as constrained to unit length),
    ##   \kappa: 0, 1 or k depending on whether we use a common given,
    ##        a common (but not given) kappa, or individual kappas.
    ##   \alpha: k - 1    (as constrained to unit sum),
    ## for a total of
    ##   k d - 1 + df_kappa(k)
    k <- length(alpha)
    attr(ll, "df") <- d * k - 1L + df_kappa(k)
    attr(ll, "nobs") <- n
    class(ll) <- "logLik"
    opt$ll <- ll
    opt$details <- list(reltol = reltol,
                        iter = c(iter = opt$iter, maxiter = maxiter),
                        logLiks = logLiks,
                        E = E,
                        kappa = control$kappa,
                        minalpha = minalpha,
                        converge = converge)
    opt
}

## Generator.

.movMF_object <-
function(theta, alpha, L, P, iter)
{
    o <- list(theta = theta, alpha = alpha, L = L, P = P, iter = iter)
    class(o) <- "movMF"
    o
}

## Methods.

print.movMF <-
function(x, ...)
{
    cat("theta:\n")
    print(x$theta)
    cat("alpha:\n")
    print(x$alpha)
    cat("L:\n")
    print(x$L)
    invisible(x)
}

coef.movMF <-
function(object, ...)
    object[c("theta", "alpha")]


logLik.movMF <-
function(object, newdata, ...)
{
  if (missing(newdata))
    return(object$ll)
  else {
    newdata <- skmeans:::row_normalize(newdata)
    return(L(newdata, object$theta, object$alpha))
  }
}

predict.movMF <-
function(object, newdata = NULL,
         type = c("class_ids", "memberships"), ...)
{
    type <- match.arg(type)
    P <- if(is.null(newdata)) object$P else {
        x <- skmeans:::row_normalize(newdata)
        theta <- object$theta
        alpha <- object$alpha
        kappa <- sqrt(rowSums(theta ^ 2))
        ## Could maybe check dimensions.
        ## Same as for E step in movMF().
        d <- nrow(theta)
        s <- d / 2 - 1 
        G <- cadd(skmeans:::g_tcrossprod(x, theta),
                  log(alpha) -  lH(kappa, s))
        g <- log_row_sums(G)
        exp(G - g)
    }
    if(type == "class_ids")
        max.col(P)
    else
        P
}

## Initializer for normalized x.

## Try something similar to what we do for skmeans, but note that
## eventually we need a matrix of posterior probabilities P.

movMF_init <-
function(x, k, start)
{
    if(is.character(start)) {
        if(start %in% c("p", "S", "s")) {
            ## <FIXME>
            ## How should we turn initial centroids into a posterior?
            ## For now, do fuzzy classification with m = 2 and cosine
            ## dissimilarity.
            M <- skmeans:::.skmeans_init_for_normalized_x(x, k, start)
            D <- pmax(1 - skmeans:::g_tcrossprod(x, M), 0)
            clue:::.memberships_from_cross_dissimilarities(D, 2)
            ## </FIXME>
        }
        else if(start == "i") {
            ## Initialize by choosing random class ids, and building a
            ## binary posterior from these.
            ids <- sample.int(k, nrow(x), replace = TRUE)
            .posterior_from_ids(ids, k)
        }
        else 
            stop(gettextf("Invalid control option 'start'"))
    }
    else if(!is.null(dim(start)))
        start
    else {
        ## A vector of class ids, hopefully.
        ids <- match(start, sort(unique(start)))
        .posterior_from_ids(ids, k)
    }
}

## Normalizing constant for the von Mises Fisher distribution.
## (In case the reference density on the unit sphere is 1.)

C <-
function(kappa, d)
{
    s <- d / 2 - 1
    kappa^s / ((2 * pi)^(s + 1) * besselI(kappa, s))
}

## Utility function for scaling the columns of a matrix.

cmult <-
function(A, x)
    A * rep.int(x, rep.int(nrow(A), ncol(A)))


cadd <-
function(A, x)
    A + rep.int(x, rep.int(nrow(A), ncol(A)))

## Utility for computing sums avoiding *flows.

log_row_sums <-
function(x)
{
    M <- x[cbind(seq_len(nrow(x)), max.col(x, "first"))]
    M + log(rowSums(exp(x - M)))
}

## Utility functions for estimating kappa.

## <NOTE>
## Work only for 0 <= Rbar < 1.
## We could try to additionally implement solve_kappa(1, ...) = Inf.
## But the utilities are internal only and the other functions cannot
## gracefully handle kappa = Inf anyways ...
## </NOTE>

solve_kappa_Banerjee_et_al_2005 <-
function(Rbar, d)
{
    if(any(Rbar >= 1))
        stop("Cannot handle infinite concentration parameters")
    Rbar * (d - Rbar ^ 2) / (1 - Rbar ^ 2)
}

solve_kappa_Tanabe_et_al_2007 <-
function(Rbar, d, c = 1, tol = 1e-6)
{
    if(any(Rbar >= 1))
        stop("Cannot handle infinite concentration parameters")
    old <- Rbar * (d - c) / (1 - Rbar ^ 2)
    repeat {
        kappa <- old * Rbar / A(old, d)
        if(max(abs(kappa - old)) < tol) break
        old <- kappa
    }
    kappa
}

solve_kappa_Sra_2012 <-
function(Rbar, d)
{
    ## Initialize using the Banerjee et al approximation.
    kappa <- solve_kappa_Banerjee_et_al_2005(Rbar, d)
    ## Perform two Newton steps.
    kappa <- kappa - (A(kappa, d) - Rbar) / Aprime(kappa, d)
    kappa <- kappa - (A(kappa, d) - Rbar) / Aprime(kappa, d)
    kappa
}

solve_kappa_uniroot <-
function(Rbar, d, tol = 1e-6)
{
    if(any(Rbar >= 1))
        stop("Cannot handle infinite concentration parameters")
    sapply(Rbar,
           function(r)
           uniroot(function(kappa) A(kappa, d) - r,
                   interval = r * (d - c(2, 0)) / (1 - r^2),
                   tol = tol)$root)
}

solve_kappa_Newton <-
function(Rbar, d, c = 1, tol = 1e-6, maxiter = 100)
{
    if(any(Rbar >= 1))
        stop("Cannot handle infinite concentration parameters")
    kappa_0 <- Inf
    kappa <- Rbar * (d - c) / (1 - Rbar ^ 2)
    A <- A(kappa, d) - Rbar
    n <- 0
    while(((abs(A) >= tol) ||
           (abs(kappa_0 - kappa) >= tol * (kappa_0 + tol)))
          && (n < maxiter)) {
        kappa_0 <- kappa
        kappa <- kappa - A / Aprime(kappa, d)
        A <- A(kappa, d) - Rbar
        n <- n + 1
    }
    kappa
}

## Utility function for computing the normalization constant actually
## used (see the implementation notes).

H <-
function(kappa, nu, v0 = 1)
{
    ## Add range checking eventually.
    n <- max(length(kappa), length(nu), length(v0))
    kappa <- rep(kappa, length.out = n)
    nu <- rep(nu, length.out = n)
    v0 <- rep(v0, length.out = n)
    .C("my0F1",
       n,
       as.double(kappa ^ 2 / 4),
       as.double(nu + 1),
       as.double(v0),
       y = double(n))$y
}

lH_asymptotic <-
function(kappa, nu)
{
    ## Add range checking eventually.    
    n <- max(length(kappa), length(nu))
    kappa <- rep(kappa, length.out = n)
    nu <- rep(nu, length.out = n)
    y <- double(n)
    ipk <- (kappa > 0)
    ipn <- (nu > 0)
    ind <- ipk & !ipn
    if(any(ind)) {
        y[ind] <- kappa[ind] - log(2 * pi * kappa[ind]) / 2
    }
    ind <- ipk & ipn
    if(any(ind)) {
        kappa <- kappa[ind]
        nu <- nu[ind]
        u <- sqrt(kappa ^ 2 + nu ^ 2)
        y[ind] <-
           (lgamma(nu + 1) - lgamma(nu + 1/2)
            + (u - nu + nu * log((2 * nu) / (u + nu))) - log(u) / 4)
    }
    y
}
        
lH <-
function(kappa, nu)
{
    ## log(H): see notes.
    n <- max(length(kappa), length(nu))
    kappa <- rep(kappa, length.out = n)
    nu <- rep(nu, length.out = n)
    ## <FIXME>
    ## Maybe use a "quick check" for determining a kappa/nu region
    ## feasible for the series expansion?
    ## </FIXME>
    y <- lH_asymptotic(kappa, nu)
    ind <- y <= 666
    if(any(ind))
        y[ind] <- log(H(kappa[ind], nu[ind]))
    ind <- !ind & (y <= 1333)
    if(any(ind)) {
        v <- y[ind] / 2
        y[ind] <- v + log(H(kappa[ind], nu[ind], exp(-v)))
    }
    y
}

## For testing this: delta should be close to 0.

delta <-
function(kappa, nu)
{
    (besselI(kappa, nu) -
     H(kappa, nu) * (kappa / 2)^nu / gamma(nu + 1))
}

## A and A prime.

A <-
function(kappa, d, method = c("PCF", "GCF", "RH"))
{
    n <- max(length(kappa), length(d))
    kappa <- rep(kappa, length.out = n)
    d <- rep(d, length.out = n)

    method <- match.arg(method)
    if(method == "PCF") {
        .C("mycfP",
           as.integer(n), as.double(kappa), as.double(d / 2),
           y = double(n))$y
    }
    else if(method == "GCF") {
        .C("mycfG",
           as.integer(n), as.double(kappa), as.double(d / 2),
           y = double(n))$y
    }
    else {
        s <- d / 2 - 1
        y <- kappa / d
        a <- lH_asymptotic(kappa, s + 1)
        ind <- (a <= 666)
        if(any(ind))
            y[ind] <- y[ind] *
                (H(kappa[ind], s + 1) /
                 H(kappa[ind], s))
        ind <- !ind & (a <= 1333)
        if(any(ind)) {
            v <- exp(- a[ind] / 2)
            y[ind] <- y[ind] *
                (H(kappa[ind], s + 1, v) /
                 H(kappa[ind], s, v))
        }
        ind <- (a > 1333)
        if(any(ind))
            y[ind] <- y[ind] *
                exp(a[ind] - lH_asymptotic(kappa[ind], s))
        if (any(y >= 1))
          stop("RH evaluation gave infeasible values which are not in the range [0, 1)")
        y
    }
}

Aprime <-
function(kappa, d)
{
    a <- A(kappa, d)
    1 - a ^ 2 - a * (d - 1) / kappa
}

## Log-likelihood

L <-
function(x, theta, alpha)
{
    sum(ldmovMF(x, theta, alpha))
}

## Utility for computing a binary posterior from the class ids.

.posterior_from_ids <-
function(ids, k = NULL)
{
    if(is.null(k))
        k <- max(ids)
    else if(max(ids) < k) {
        k <- max(ids)
        warning("due to initialization number of components reduced to ",
                k)
    }
    else if(max(ids) > k)
        stop("number of components k smaller than those provided for initialization")
    n <- length(ids)
    P <- matrix(0, n, k)
    P[cbind(seq_len(n), ids)] <- 1
    P
}
