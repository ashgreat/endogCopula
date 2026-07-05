# Build a coefficient table with z statistics and p-values when bootstrap
# standard errors are available; returns NULL otherwise.
coefficient_table <- function(object) {
  coefs <- object$coefficients
  if (!is.matrix(coefs) || !all(c("Estimate", "Std. Error") %in% colnames(coefs))) {
    return(NULL)
  }
  estimate <- coefs[, "Estimate"]
  std_error <- coefs[, "Std. Error"]
  z_value <- estimate / std_error
  p_value <- 2 * stats::pnorm(abs(z_value), lower.tail = FALSE)
  out <- cbind(
    Estimate = estimate,
    `Std. Error` = std_error,
    `z value` = z_value,
    `Pr(>|z|)` = p_value
  )
  rownames(out) <- rownames(coefs)
  out
}

print_coefficients <- function(object, ...) {
  table <- coefficient_table(object)
  if (is.null(table)) {
    print(object$coefficients)
  } else {
    stats::printCoefmat(table, P.values = TRUE, has.Pvalue = TRUE, ...)
  }
  invisible(object)
}

#' Print method for `endog_copula_fit`
#'
#' @param x An object returned by one of the copula estimators.
#' @param ... Passed to [stats::printCoefmat()] when bootstrap standard errors
#'   are available.
#' @return Invisibly returns `x`, after printing the estimator name, the CDF
#'   transformation, and the coefficient table.
#' @examples
#' data(sim_endog)
#' fit <- CopRegPG(y ~ z_endog | x_exog + w_instr, data = sim_endog,
#'                 cdf = "resc.ecdf", nboots = 0)
#' print(fit)
#' @export
print.endog_copula_fit <- function(x, ...) {
  cat(sprintf("Gaussian copula estimator: %s\n", x$method))
  if (!is.na(x$cdf)) {
    cat(sprintf("CDF transformation: %s\n", x$cdf))
  }
  cat("\nCoefficients:\n")
  print_coefficients(x, ...)
  invisible(x)
}

#' Summary method for `endog_copula_fit`
#'
#' @param object An object returned by one of the copula estimators.
#' @param ... Passed to downstream methods.
#' @return An object of class `summary.endog_copula_fit`: a list carrying the
#'   original `call`, `method`, `cdf`, and `coefficients`, plus a
#'   `coefficient_table` with z statistics and p-values (or `NULL` when no
#'   bootstrap standard errors are available), the `residuals`, the matrix of
#'   `bootstrap` replicates, and the number of bootstrap draws `nboots`.
#' @examples
#' data(sim_endog)
#' fit <- CopRegPG(y ~ z_endog | x_exog + w_instr, data = sim_endog,
#'                 cdf = "resc.ecdf", nboots = 0)
#' summary(fit)
#' @export
summary.endog_copula_fit <- function(object, ...) {
  out <- list(
    call = object$call,
    method = object$method,
    cdf = object$cdf,
    coefficients = object$coefficients,
    coefficient_table = coefficient_table(object),
    residuals = object$residuals,
    bootstrap = object$bootstrap,
    nboots = if (is.null(object$bootstrap)) 0L else nrow(object$bootstrap)
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
  if (is.null(x$coefficient_table)) {
    print(x$coefficients)
  } else {
    stats::printCoefmat(x$coefficient_table, P.values = TRUE, has.Pvalue = TRUE, ...)
  }
  if (x$nboots > 0) {
    cat(sprintf("\nStandard errors based on %d bootstrap replicates.\n", x$nboots))
  }
  invisible(x)
}

#' @export
coef.endog_copula_fit <- function(object, ...) {
  if (is.matrix(object$coefficients)) {
    estimates <- object$coefficients[, 1]
    names(estimates) <- rownames(object$coefficients)
    estimates
  } else {
    object$coefficients
  }
}

#' @export
residuals.endog_copula_fit <- function(object, ...) {
  object$residuals
}

#' Confidence intervals for `endog_copula_fit` objects
#'
#' Computes confidence intervals for the estimated coefficients. The default
#' `type = "normal"` uses a normal approximation based on the bootstrap
#' standard errors; `type = "percentile"` uses percentiles of the bootstrap
#' replicates stored in `object$bootstrap`. Both require the model to have
#' been fitted with `nboots > 0`.
#'
#' @param object An object returned by one of the copula estimators.
#' @param parm Coefficients to compute intervals for, either names or indices.
#'   Defaults to all coefficients.
#' @param level Confidence level. Defaults to `0.95`.
#' @param type Either `"normal"` (normal approximation from bootstrap standard
#'   errors) or `"percentile"` (percentile intervals from the bootstrap
#'   replicates).
#' @param ... Currently unused.
#' @return A matrix with one row per coefficient and columns giving the lower
#'   and upper confidence limits, labelled with the corresponding percentiles.
#'   Coefficients without bootstrap information are filled with `NA`.
#' @examples
#' data(sim_endog)
#' fit <- CopRegPG(y ~ z_endog | x_exog + w_instr, data = sim_endog,
#'                 cdf = "resc.ecdf", nboots = 25)
#' confint(fit)
#' confint(fit, parm = "z_endog", type = "percentile")
#' @export
confint.endog_copula_fit <- function(object, parm, level = 0.95,
                                     type = c("normal", "percentile"), ...) {
  type <- match.arg(type)
  if (!is.numeric(level) || length(level) != 1 || level <= 0 || level >= 1) {
    stop("Argument 'level' must be a single number strictly between 0 and 1.", call. = FALSE)
  }
  estimates <- stats::coef(object)
  if (missing(parm)) {
    parm <- names(estimates)
  } else if (is.numeric(parm)) {
    parm <- names(estimates)[parm]
  }
  alpha <- (1 - level) / 2
  probs <- c(alpha, 1 - alpha)
  labels <- paste(format(100 * probs, trim = TRUE, scientific = FALSE, digits = 3), "%")
  out <- matrix(NA_real_, nrow = length(parm), ncol = 2,
                dimnames = list(parm, labels))
  if (type == "percentile") {
    if (is.null(object$bootstrap)) {
      stop(
        "Percentile intervals require bootstrap replicates; refit with 'nboots' > 0.",
        call. = FALSE
      )
    }
    available <- intersect(parm, colnames(object$bootstrap))
    for (name in available) {
      out[name, ] <- stats::quantile(object$bootstrap[, name], probs = probs,
                                     na.rm = TRUE, names = FALSE)
    }
  } else {
    coefs <- object$coefficients
    if (!is.matrix(coefs) || !("Std. Error" %in% colnames(coefs))) {
      stop(
        "Normal-approximation intervals require bootstrap standard errors; refit with 'nboots' > 0.",
        call. = FALSE
      )
    }
    std_errors <- coefs[, "Std. Error"]
    names(std_errors) <- rownames(coefs)
    z <- stats::qnorm(1 - alpha)
    available <- intersect(parm, names(estimates))
    out[available, 1] <- estimates[available] - z * std_errors[available]
    out[available, 2] <- estimates[available] + z * std_errors[available]
  }
  out
}
