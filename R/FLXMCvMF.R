.onLoad <-
function(libname, pkgname) {    
if (requireNamespace("methods", quietly = TRUE) && requireNamespace("flexmix", quietly = TRUE)) {    
    methods::setClass("FLXMCvMF",
                      contains =
                          methods::getClassDef("FLXMCsparse",
                                               where = asNamespace("flexmix")))
    methods::setMethod(flexmix::FLXmstep, methods::signature(model = "FLXMCvMF"),
                       function(model, weights, ...) {
        model@fit(model@x, model@y, weights)
    })
}
}

FLXMCvMF <- function(formula = . ~ ., kappa = NULL) {
  if (is.null(methods::getClassDef("FLXMCvMF", where = asNamespace("movMF")))) {
      stop("The flexmix model driver requires packages methods and flexmix. Please re-install.")
  }
  z <- methods::new("FLXMCvMF", weighted = TRUE, formula = formula,
                    name = "model-based von Mises-Fisher clustering")
  z@preproc.y <- function(x) {
      x <- x / row_norms(x)
      x
  }
  
  z@defineComponent <- function(para) {
      logLik <- function(x, y) {
          dmovMF(y, para$theta, log = TRUE)
      }
      
      predict <- function(x, ...) {
          mu <- para$theta / sqrt(sum(para$theta^2))
          matrix(mu, nrow = nrow(x), ncol = length(mu), byrow = TRUE)
      }
    
      methods::new("FLXcomponent", parameters = list(theta = para$theta), df = para$df,
                   logLik = logLik, predict = predict)
  }

  solve_kappa <- get_solve_kappa(kappa)
  do_kappa <- solve_kappa$do_kappa
  df_kappa <- solve_kappa$df_kappa
  use_common_kappa <- solve_kappa$use_common_kappa
  
  z@fit <- function(x, y, w, ...) {
      M <- skmeans:::g_crossprod(w, y)
      norms <- row_norms(M)
      M <- M / ifelse(norms > 0, norms, 1)
      kappa <- do_kappa(norms, skmeans:::g_col_sums(w), ncol(y))
      theta <- kappa * M
      k <- ncol(w)
      df <- df_kappa(k) / k + ncol(y) - 1
      lapply(seq_len(k), function(K)
          z@defineComponent(list(theta = theta[K, ], df = df)))
  }
  z 
}

