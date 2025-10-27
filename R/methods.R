#' Print method for `endog_copula_fit`
#'
#' @param x An object returned by one of the copula estimators.
#' @param ... Passed to downstream methods.
#' @export
print.endog_copula_fit <- function(x, ...) {
  cat(sprintf("Gaussian copula estimator: %s\n", x$method))
  if (!is.na(x$cdf)) {
    cat(sprintf("CDF transformation: %s\n", x$cdf))
  }
  cat("\nCoefficients:\n")
  print(x$coefficients)
  invisible(x)
}

#' Summary method for `endog_copula_fit`
#'
#' @param object An object returned by one of the copula estimators.
#' @param ... Passed to downstream methods.
#' @export
summary.endog_copula_fit <- function(object, ...) {
  out <- list(
    call = object$call,
    method = object$method,
    cdf = object$cdf,
    coefficients = object$coefficients,
    residuals = object$residuals,
    bootstrap = object$bootstrap
  )
  class(out) <- "summary.endog_copula_fit"
  out
}

#' @export
print.summary.endog_copula_fit <- function(x, ...) {
  cat(sprintf("Gaussian copula estimator: %s\n", x$method))
  if (!is.na(x$cdf)) {
    cat(sprintf("CDF transformation: %s\n", x$cdf))
  }
  cat("\nCall:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  print(x$coefficients)
  invisible(x)
}

#' @export
coef.endog_copula_fit <- function(object, ...) {
  if (is.matrix(object$coefficients)) {
    object$coefficients[, 1]
  } else {
    object$coefficients
  }
}

#' @export
residuals.endog_copula_fit <- function(object, ...) {
  object$residuals
}
