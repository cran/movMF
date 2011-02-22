## Try to estimate movMF via EM.

movMF <-
function(x, k, control = list(), ...)
{
    ## Be nice a la David.
    control <- c(control, list(...))

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
                   t(apply(exp(G - g), 1, function(prob) rmultinom(1, 1, prob))))

    kappa_estimators <-
        c("Banerjee_et_al_2005",
          "Tanabe_et_al_2007",
          "Sra_2011",
          "uniroot",
          "Newton")
    kappa <- control$kappa
    if(is.null(kappa))
        do_kappa <- function(Rbar)
            do_kappa_Newton(Rbar, d)
    else {
        kappa <- as.list(kappa)
        kname <- kappa[[1L]]
        kargs <- kappa[-1L]
        pos <- pmatch(tolower(kname), tolower(kappa_estimators))
        if(is.na(pos))
            stop("Invalid kappa estimator.")
        kappa <- kappa_estimators[pos]
        kfun <- get(sprintf("do_kappa_%s", kappa))
        do_kappa <- function(Rbar)
            do.call(kfun, c(list(Rbar, d), kargs))
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
    } else {
        ## Initialization.
        start <- control$start
        if(is.null(start)) {
            nruns <- control$nruns
            if(is.null(nruns))
                nruns <- 1L
            start <- rep.int("p", nruns)
        }
        start <- movMF_init(x, k, start)
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

    if(verbose && (nruns > 1L))
        message(gettextf("Run: %d", run))

    repeat {
        G <- NULL
        P <- start[[run]]
        L_old <- -Inf
        iter <- 1L
        
        while(iter <= maxiter) {
            ## M step.
            alpha <- colMeans(P)
            while(any(alpha < minalpha)) {
                if(verbose) 
                    message("*** Removing one component ***")
                nok <- which.min(alpha)
                P <- P[, -nok, drop = FALSE]
                P <- do_P(P, log_row_sums(P))
                alpha <- colSums(P) / n
                if(!is.null(G))
                    L_old <- sum(log_row_sums(G[, -nok, drop = FALSE]))
            }
            M <- skmeans:::g_crossprod(P, x)
            norms <- skmeans:::row_norms(M)
            Rbar <- norms / (n * alpha)
            M <- M / norms
            kappa <- do_kappa(Rbar)
            
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
        
        if(run >= nruns) break
        
        run <- run + 1L
        if(verbose)
            message(gettextf("Run: %d", run))
    }

    ## Compute log-likelihood.
    ll <- L(x, opt$theta, opt$alpha)
    ## Add the "degrees of freedom" (number of (estimated) parameters in
    ## the model): \mu: k * (d - 1) (as constrained to unit length),
    ## \kappa: k, \alpha: k - 1 (as constrained to unit sum), for a
    ## total of (d + 1) k - 1
    attr(ll, "df") <- (d + 1L) * k - 1L
    class(ll) <- "logLik"
    opt$ll <- ll

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
function(object, ...)
{
    object$ll
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
        if(any(is.na(match(start, c("p", "i", "S", "s")))))
            stop(gettextf("Invalid control option 'start'"))
        out <- vector("list", length(start))
        pos <- which(start %in% c("p", "S", "s"))
        if(length(pos)) {
            ## <FIXME>
            ## How should we turn initial centroids into a posterior?
            ## For now, do fuzzy classification with m = 2 and cosine
            ## dissimilarity.
            out[pos] <-
                lapply(skmeans:::.skmeans_init_for_normalized_x(x, k,
                                                                start[pos]),
                       function(M) {
                           clue:::.memberships_from_cross_dissimilarities(1 -
                                                                          skmeans:::g_tcrossprod(x,
                                                                                                 M),
                                                                          2)
                       })
            ## </FIXME>
        }
        pos <- which(start == "i")
        if(length(pos)) {
            out[pos] <-
                replicate(length(pos), {
                    ## Initialize by choosing random class ids, and
                    ## building a binary posterior from these.
                    ids <- sample.int(k, nrow(x), replace = TRUE)
                    .posterior_from_ids(ids, k)
                })
        }
        out
    } else {
        if(!is.list(start))
            start <- list(start)
        lapply(start,
               function(s) {
                   if(!is.null(dim(s)))
                       s
                   else {
                       ## A vector of class ids, hopefully.
                       ids <- match(s, unique(s))
                       .posterior_from_ids(ids, k)
                   }
               })
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
    M <- x[cbind(seq_len(nrow(x)), max.col(x))]
    M + log(rowSums(exp(x - M)))
}

## Utility functions for estimating kappa.

## <NOTE>
## Work only for 0 <= Rbar < 1.
## We could try to additionally implement do_kappa(1, ...) = Inf.
## But the utilities are internal only and the other functions cannot
## gracefully handle kappa = Inf anyways ...
## </NOTE>

do_kappa_Banerjee_et_al_2005 <-
function(Rbar, d)
{
    Rbar * (d - Rbar ^ 2) / (1 - Rbar ^ 2)
}

do_kappa_Tanabe_et_al_2007 <-
function(Rbar, d, c = 1, tol = 1e-6)
{
    old <- Rbar * (d - c) / (1 - Rbar ^ 2)
    repeat {
        kappa <- old * Rbar / A(old, d)
        if(max(abs(kappa - old)) < tol) break
        old <- kappa
    }
    kappa
}

do_kappa_Sra_2011 <-
function(Rbar, d)
{
    ## Initialize using the Banerjee et al approximation.
    kappa <- do_kappa_Banerjee_et_al_2005(Rbar, d)
    ## Perform two Newton steps.
    kappa <- kappa - (A(kappa, d) - Rbar) / Aprime(kappa, d)
    kappa <- kappa - (A(kappa, d) - Rbar) / Aprime(kappa, d)
    kappa
}

do_kappa_uniroot <-
function(Rbar, d, tol = 1e-6)
{
    sapply(Rbar,
           function(r)
           uniroot(function(kappa) A(kappa, d) - r,
                   interval = r * (d - c(2, 0)) / (1 - r^2),
                   tol = tol)$root)
}

do_kappa_Newton <-
function(Rbar, d, c = 1, tol = 1e-6, maxiter = 100)
{
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
        y[ind] <- kappa[ind] - log(2 * pi * kappa[ind])
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
function(kappa, d, method = c("CF", "RH"))
{
    n <- max(length(kappa), length(d))
    kappa <- rep(kappa, length.out = n)
    d <- rep(d, length.out = n)

    s <- d / 2 - 1
    y <- kappa / d
    
    method <- match.arg(method)
    if(method == "CF") {
        y <- y * .C("mygcf",
                    n,
                    as.double(kappa ^ 2 / 4),
                    as.double(s + 1),
                    y = double(n))$y
    }
    else {
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
    }
    
    y
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
function(ids, k)
{
    n <- length(ids)
    P <- matrix(0, n, k)
    P[cbind(seq_len(n), ids)] <- 1
    P
}
