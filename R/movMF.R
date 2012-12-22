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

    nu <- if (is.null(control$nu)) 0 else control$nu
    
    kappa_solvers <-
        c("Banerjee_et_al_2005",
          "Tanabe_et_al_2007",
          "Sra_2012",
          "Song_et_al_2012",
          "uniroot",
          "Newton")
    kappa <- control$kappa
    if(is.numeric(kappa)) {
        ## CASE A: Use a common given value of kappa.
        if (length(kappa) > 1)
          warning("only the first element of 'kappa' is used for a common given kappa")
        kappa_given <- kappa[1]
        use_common_kappa <- TRUE
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
            if (length(nu) > 1)
              warning("only the first element of 'nu' is used for common kappa")
            nu <- nu[1]
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
                rep.int(solve_kappa(sum(norms) / (n + nu)), length(alpha))
            df_kappa <- function(k) 1L
        } else {
            do_kappa <- function(norms, alpha)
                solve_kappa(norms / (n * alpha + nu))
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
        if (nrow(x) != length(ids))
          stop("length of 'ids' needs to match the number of observations")
        P0 <- .posterior_from_ids(ids, k)
        do_P <- function(G, g) G
        start  <- list(P0)
        maxiter <- 1L
        if(!is.null(control$nruns))
            warning("control argument 'nruns' ignored because ids are specified")
        if(!is.null(control$start))
            warning("control argument 'start' ignored because ids are specified")
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
                warning("control argument 'nruns' ignored because 'start' is specified")
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
        message(gettextf("Run: %d", run),
                domain = NA)

    repeat {
        G <- NULL
        P <- movMF_init(x, k, start[[run]])
        if (!use_common_kappa) {
          if (length(nu) > 1 & length(nu) != ncol(P))
            warning("nu is changed to have length", ncol(P))
          nu <- rep(nu, length.out = ncol(P))
        }
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
                    if (!use_common_kappa)
                      nu <- nu[-nok]
                    P <- do_P(P, log_row_sums(P))
                    alpha <- colMeans(P)
                    if(!is.null(G))
                        L_old <- sum(log_row_sums(G[, -nok, drop = FALSE]))
                }
                if(any(alpha == 0 & nu <= 0))
                    stop("Cannot handle empty components")
                M <- skmeans:::g_crossprod(P, x)
                norms <- skmeans:::row_norms(M)
                M <- M / ifelse(norms > 0, norms, 1)
                ## If a cluster contains only identical observations,
                ## Rbar = 1.
                kappa <- do_kappa(norms, alpha)
            
                ## E step.
                G <- cadd(skmeans:::g_tcrossprod(x, kappa * M),
                          log(alpha) -  lH(kappa, s))
                g <- log_row_sums(G)
                L_new <- sum(g)
                if(verbose && (iter %% verbose == 0))
                    message(gettextf("Iteration: %d *** L: %g",
                                     iter, L_new),
                            domain = NA)
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
                message(gettextf("EM algorithm did not converge:\n%s",
                                 s),
                        domain = NA)
                ## </NOTE>
            }
            NA
        })
        
        if(run >= nruns) break
        
        run <- run + 1L
        if(verbose)
            message(gettextf("Run: %d", run),
                    domain = NA)
    }

    ## Compute log-likelihood.
    if(is.null(opt))
        stop("EM algorithm did not converge for any run")
    k <- length(alpha)
    dimnames(opt$theta)<- list(seq_len(k), colnames(x))
    dimnames(opt$P)<- list(rownames(x), seq_len(k))
    ll <- L(x, opt$theta, opt$alpha)
    ## Add the "degrees of freedom" (number of (estimated) parameters in
    ## the model): with k the number of classes actually used,
    ##   \mu: k * (d - 1) (as constrained to unit length),
    ##   \kappa: 0, 1 or k depending on whether we use a common given,
    ##        a common (but not given) kappa, or individual kappas.
    ##   \alpha: k - 1    (as constrained to unit sum),
    ## for a total of
    ##   k d - 1 + df_kappa(k)
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
            stop("Invalid control option 'start'")
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
      kappa <- ifelse(old == 0, 0, old * (Rbar / A(old, d)))
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
    A <- A(kappa, d)
    kappa <- kappa - (A - Rbar) / Aprime(kappa, d, A)
    A <- A(kappa, d)
    kappa <- kappa - (A - Rbar) / Aprime(kappa, d, A)
    kappa
}

solve_kappa_Song_et_al_2012 <-
function(Rbar, d)
{
    ## Initialize using the Banerjee et al approximation.
    kappa <- solve_kappa_Banerjee_et_al_2005(Rbar, d)
    ## Perform two Halley steps.
    A <- A(kappa, d)
    Adiff <- A - Rbar
    Aprime <- Aprime(kappa, d, A)
    kappa <- kappa - 2 * Adiff * Aprime / (2 * Aprime^2 - Adiff * Adoubleprime(kappa, d, A))
    A <- A(kappa, d)
    Adiff <- A - Rbar
    Aprime <- Aprime(kappa, d, A)
    kappa <- kappa - 2 * Adiff * Aprime / (2 * Aprime^2 - Adiff * Adoubleprime(kappa, d, A))
    kappa
}

solve_kappa_uniroot <-
function(Rbar, d, tol = 1e-6)
{
    if(any(Rbar >= 1))
        stop("Cannot handle infinite concentration parameters")
    
    nu <- d / 2 - 1
    sapply(Rbar,
           function(r) {
             interval <- c(Rinv_lower_Amos_bound(r, nu),
                           Rinv_upper_Amos_bound(r, nu))
             if (abs(diff(interval)) < tol)
               mean(interval)
             else 
               uniroot(function(kappa) A(kappa, d) - r,
                       interval = interval,
                       tol = tol)$root})
}

solve_kappa_Newton <-
function(Rbar, d, tol = 1e-6, maxiter = 100)
{
    if(any(Rbar >= 1))
        stop("Cannot handle infinite concentration parameters")
    
    kappa_0 <- Inf
    kappa <- Rinv_lower_Amos_bound(Rbar, d / 2 - 1)
    A <- A(kappa, d)
    Adiff <- A - Rbar
    n <- 0
    while(((abs(Adiff) >= tol) ||
           (abs(kappa_0 - kappa) >= tol * (kappa_0 + tol)))
          && (n < maxiter)) {
        kappa_0 <- kappa
        kappa <- kappa - Adiff / Aprime(kappa, d, A)
        A <- A(kappa, d)
        Adiff <- A - Rbar
        n <- n + 1
    }
    kappa
}

## Utility functions for computing A, H and its logarithm (see the
## implementation notes).
##
## With R_\nu = I_{\nu+1} / I_\nu the commonly considered modified
## Bessel function ratio and
##
##   log(H_\nu)(\kappa) = \int_0^\kappa R_\nu(t) dt
##
## (so that (log(H_\nu))' = R_\nu and H_\nu(0) = 1), we have
##
##   A_d = I_{d/2} / I_{d/2-1} = R_{d/2-1}
##
## and
##
##   log(f(x|theta)) = theta' x - log(H_{d/2-1})(\|theta\|).
##
## For R_\nu we have Amos-type bounds
##
##   G_{\alpha,\beta}(t) = t / (alpha + sqrt(t^2 + beta^2))
##
## where G_{\nu+1/2,\nu+3/2} is the uniformly best lower Amos-type
## bound, whereas there is no uniformly optimal upper bound, with
## G_{\nu,\nu+2} and G_{\nu+1/2,beta_SS(\nu)} being second-order exact
## and optimal at 0 and infinity, respectively, where
##
##   beta_SS(\nu) = sqrt((\nu+1/2)(\nu+3/2)).
##
## Inverses and anti-derivatives of G are given by
##
##   G_{\alpha,\beta}^{-1}(\rho)
##     = \frac{\rho}{1 - \rho^2}
##       (\alpha + \sqrt{\alpha^2 \rho^2 + \beta^2 (1 - \rho^2)})
##
## and
##
##   S_{\alpha,\beta}(\kappa)
##     = \sqrt{\kappa^2 + \beta^2} - \beta
##       - \alpha \log(\alpha + \sqrt{\kappa^2 + \beta^2})
##       + \alpha \log(\alpha + \beta).
##
## We use
##
##   \max(G_{\nu,\nu+2}^{-1}, G_{\nu+1/2,\beta_{SS}(\nu)}^{-1})
##
## as lower bound and approximation for R_\nu^{-1},
##
##   G_{\nu+1/2,\nu+3/2}^{-1}
##
## as upper bound for R_{\nu}^{-1}, and
##
##   \min(S_{\nu,\nu+2}, S_{\nu+1/2,\beta_{SS}(\nu)})
##
## as approximation for \log(H_\nu).

beta_SS <-
function(nu)
    sqrt((nu + 1/2) * (nu + 3/2))

Ginv <-
function(rho, alpha, beta)
{
    ## Perhaps recycle arguments eventually ...
    
    sigma <- rho^2
    rho * (alpha + sqrt(alpha^2 * sigma + beta^2 * (1 - sigma))) /
        (1 - sigma)
}

Rinv_lower_Amos_bound <-
function(rho, nu)
{
    ## Perhaps recycle arguments eventually ...
    pmax(Ginv(rho, nu, nu + 2),
         Ginv(rho, nu + 1/2, beta_SS(nu)))
}

Rinv_upper_Amos_bound <-
function(rho, nu)
{
    ## Perhaps recycle arguments eventually ...
    Ginv(rho, nu + 1/2, nu + 3/2)
}

S <-
function(kappa, alpha, beta)
{
    u <- sqrt(kappa^2 + beta^2)
    u - beta - alpha * log((alpha + u) / (alpha + beta))
}

## Utility functions for computing H and log(H).

H <-
function(kappa, nu, v0 = 1)
{
    ## Compute v0 H_\nu(\kappa) by direct Taylor series summation.
    
    ## Add range checking eventually.
    n <- max(length(kappa), length(nu), length(v0))
    kappa <- rep(kappa, length.out = n)
    nu <- rep(nu, length.out = n)
    v0 <- rep(v0, length.out = n)
    
    .C(C_my0F1,
       as.integer(n),
       as.double(kappa ^ 2 / 4),
       as.double(nu + 1),
       as.double(v0),
       y = double(n))$y
}

lH_asymptotic <-
function(kappa, nu)
{
    ## Compute a suitable asymptotic approximation to
    ## \log(H_\nu(\kappa)).
    
    ## Add range checking eventually.    
    n <- max(length(kappa), length(nu))
    kappa <- rep(kappa, length.out = n)
    nu <- rep(nu, length.out = n)
    
    y <- double(n)
    ipk <- (kappa > 0)
    ipn <- (nu > 0)
    ind <- ipk & !ipn
    if(any(ind)) {
        ## For \log(H_0) = \log(I_0), use the asymptotic approximation
        ##   I_0(\kappa) \approx e^\kappa / \sqrt{2 \pi \kappa}
        ## (e.g., http://dlmf.nist.gov/10.40).
        y[ind] <- kappa[ind] - log(2 * pi * kappa[ind]) / 2
    }
    ind <- ipk & ipn
    if(any(ind)) {
        ## For \nu > 0, use the Amos-type approximation discussed above.
        kappa <- kappa[ind]
        nu <- nu[ind]
        y[ind] <-
            pmin(S(kappa, nu, nu + 2),
                 S(kappa, nu + 1/2, beta_SS(nu)))
    }
    y
}
        
lH <-
function(kappa, nu)
{
    ## Compute \log(H_\nu(\kappa)) (or an approximation to it).
    ## See the implementation notes for details.
    
    n <- max(length(kappa), length(nu))
    kappa <- rep(kappa, length.out = n)
    nu <- rep(nu, length.out = n)

    y <- lH_asymptotic(kappa, nu)
    ## If the value from the asymptotic approximation is small enough,
    ## we can use direct Taylor series summation.
    ind <- y <= 699.5
    if(any(ind))
        y[ind] <- log(H(kappa[ind], nu[ind]))
    ## For intermediate values of the asymptotic approximation, we use
    ## rescaling and direct Taylor series summation.
    ind <- !ind & (y <= 1399)
    if(any(ind)) {
        v <- y[ind] / 2
        y[ind] <- v + log(H(kappa[ind], nu[ind], exp(-v)))
    }
    ## (Otherwise, we use the asymptotic approximation.)
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
function(kappa, d, method = c("PCF", "GCF", "RH"), tol = 1e-6)
{
    method <- match.arg(method)
    n <- max(length(kappa), length(d))
    kappa <- rep(kappa, length.out = n)
    d <- rep(d, length.out = n)
    A <- vector("numeric", length = n)
    
    index <- kappa >= tol
    if (sum(index)) {
      method <- match.arg(method)
      if(method == "PCF") {
        ## Use the Perron continued fraction for R_{\nu-1}.
        ## Implementation based on Eqn 3.3' in Gautschi and Slavik
        ## (1978).
        ## Note that A_d = R_{d/2-1}.
        A[index] <- .C(C_mycfP,
                       as.integer(sum(index)),
                       as.double(kappa[index]),
                       as.double(d[index] / 2),
                       y = double(sum(index)))$y
      }
      else if(method == "GCF") {
        ## Use the Gauss continued fraction for R_{\nu-1}.
        ## Implementation based on Eqn 3.2' in Gautschi and Slavik
        ## (1978).
        ## Note that A_d = R_{d/2-1}.
        A[index] <- .C(C_mycfG,
                       as.integer(sum(index)),
                       as.double(kappa[index]),
                       as.double(d[index] / 2),
                       y = double(sum(index)))$y
      }
      else {
        ## Compute A_d via a ratio of H functions:
        ##   A_d(\kappa)
        ##     = (kappa/d) * H_{d/2}(\kappa) / H_{d/2-1}(\kappa).
        s <- d[index] / 2 - 1
        kappai <- kappa[index]
        y <- kappai / d[index]
        a <- lH_asymptotic(kappai, s + 1)
        ind <- (a <= 699.5)
        if(any(ind))
            y[ind] <- y[ind] *
                (H(kappai[ind], s + 1) /
                 H(kappai[ind], s))
        ind <- !ind & (a <= 1399)
        if(any(ind)) {
            v <- exp(- a[ind] / 2)
            y[ind] <- y[ind] *
                (H(kappai[ind], s + 1, v) /
                 H(kappai[ind], s, v))
        }
        ind <- (a > 1399)
        if(any(ind))
            y[ind] <- y[ind] *
                exp(a[ind] - lH_asymptotic(kappai[ind], s))
        if(any(y >= 1))
          stop("RH evaluation gave infeasible values which are not in the range [0, 1)")
        A[index] <- y
      }
    }
    if (sum(!index)) {
      di <- d[!index]
      kappai <- kappa[!index]
      A[!index] <- kappai / di - kappai^3 / (di^2 * (di + 2)) + 2 * kappai^5 / (di^3 * (di + 2) * (di + 4))
    }
    A 
}

Aprime <-
function(kappa, d, a = NULL, tol = 1e-6, ...)
{
  n <- max(length(kappa), length(d), length(a))
  kappa <- rep(kappa, length.out = n)
  d <- rep(d, length.out = n)
  a <- rep(a, length.out = n)
  aprime <- vector("numeric", length = n)
  
  index <- kappa >= tol
  if (sum(index)) {
    if (is.null(a)) 
      a <- A(kappa[index], d[index], tol = tol, ...)
    else
      a <- a[index]
    aprime[index] <- 1 - a ^ 2 - a * (d[index] - 1) / kappa[index]
  }
  if (sum(!index)) {
    di <- d[!index]
    aprime[!index] <- 1 / di - 3 / (di^2 * (di + 2)) * kappa[!index]^2 + 10 / (di^3 * (di + 2) * (di + 4)) * kappa[!index]^4
  }
  aprime
}

Adoubleprime <-
function(kappa, d, a = NULL, tol = 1e-6, ...)
{
  n <- max(length(kappa), length(d), length(a))
  kappa <- rep(kappa, length.out = n)
  d <- rep(d, length.out = n)
  a <- rep(a, length.out = n)
  adoubleprime <- vector("numeric", length = n)
  
  index <- kappa >= tol
  if (sum(index)) {
    di <- d[index]
    if (is.null(a)) 
      a <- A(kappa[index], di, tol = tol, ...)
    else
      a <- a[index]
    kappa2 <- kappa[index]^2
    adoubleprime[index] <- 2 * a^3 + 3 * (di - 1) / kappa[index] * a^2 + (di^2 - di - 2 * kappa2) / kappa2 * a - (di - 1) / kappa[index]
  }
  if (sum(!index)) {
    di <- d[!index]
    adoubleprime[!index] <- - 6 / (di^2 * (di + 2)) * kappa[!index] + 40 / (di^3 * (di + 2) * (di + 4)) * kappa[!index]^3
  }
  adoubleprime
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
        warning(gettextf("due to initialization number of components reduced to  %d",
                         k),
                domain = NA)
    }
    else if(max(ids) > k)
        stop("number of components k smaller than those provided for initialization")
    n <- length(ids)
    P <- matrix(0, n, k)
    P[cbind(seq_len(n), ids)] <- 1
    P
}
