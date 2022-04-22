rvMF <-
function(n, theta)
{
    d <- length(theta)
    kappa <- sqrt(sum(theta ^ 2))

    if(kappa == 0) {
        y <- matrix(rnorm(n * d), n, d)
        y <- y / sqrt(rowSums(y ^ 2))
    }
    else if(d == 1) {
        y <- cbind((-1) ^ rbinom(n, 1, 1 / (1 + exp(2 * theta))))
    }
    else {
        ## <FIXME>
        ## Explain in the package vignette eventually ...
        mu <- cbind(theta / kappa)
        w <- rW(n, kappa, d)
        z <- matrix(rnorm(n * d), n, d)
        v <- z - tcrossprod(z %*% mu, mu)
        v <- v / sqrt(rowSums(v ^ 2))
        y <- tcrossprod(w, mu) + sqrt(1 - w^2) * v
        ## </FIXME>
    }

    y
}

rmovMF <-
function(n, theta, alpha = 1)
{
    ## Be nice to users.
    theta <- rbind(theta)
    k <- max(nrow(theta), length(alpha))
    theta <- theta[rep_len(seq_len(nrow(theta)), k), , drop = FALSE]
    alpha <- rep_len(alpha, k)
    alpha <- alpha / sum(alpha)

    y <- matrix(0, n, ncol(theta))
    ind <- as.numeric(cut(runif(n), c(0, cumsum(alpha)),
                          include.lowest = TRUE))
    pos <- split(seq_len(n), ind)
    nms <- names(pos)
    for(i in seq_along(pos)) {
        j <- as.numeric(nms[i])
        p <- pos[[i]]
        y[p, ] <- rvMF(length(p), theta[j, ])
    }

    attr(y, "z") <- ind
    class(y) <- "rmovMF"
    y
}

print.rmovMF <-
function(x, ...)
{
    print.default(matrix(c(x), dim(x)), ...)
    invisible(x)
}
    
rW <-
function(n, kappa, d)
{
    .C(C_rW,
       as.integer(n),
       as.double(kappa),
       as.integer(d),
       y = double(n))$y
}
