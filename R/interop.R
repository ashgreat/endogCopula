#' Panel Gaussian Copula Estimator (delegated)
#'
#' Thin wrapper around `endogCopulaPanel::CopRegML_par()`. The function is only
#' available when the optional `endogCopulaPanel` package is installed.
#'
#' @param formula Two-part formula `y ~ endog | exog` describing the model.
#' @param index Character vector of length two with panel and time identifiers.
#' @param data Data frame containing the referenced variables.
#' @param ecdf Logical; use empirical CDFs (`TRUE`) or kernel estimates (`FALSE`).
#' @param nboots Number of bootstrap replications.
#' @param starting_values Optional numeric vector of starting values.
#' @param method Optimisation method passed to `stats::optim()`.
#' @return Delegates to `endogCopulaPanel::CopRegML_par()`.
#' @export
CopRegML_par <- function(formula, index, data, ecdf = TRUE, nboots = 199,
                         starting_values = NULL, method = "Nelder-Mead") {
  pkg <- "endogCopulaPanel"
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Install 'endogCopulaPanel' to use the panel estimator.", call. = FALSE)
  }
  panel_fun <- getExportedValue(pkg, "CopRegML_par")
  panel_fun(
    formula = formula,
    index = index,
    data = data,
    ecdf = ecdf,
    nboots = nboots,
    starting_values = starting_values,
    method = method
  )
}

#' Bayesian Gaussian Copula Estimator (delegated)
#'
#' Thin wrapper around `endogCopulaBayes::CopRegBayes()`. Requires the optional
#' `endogCopulaBayes` package to be installed.
#'
#' @param data Data frame containing the columns `y`, `z`, and `x`.
#' @param iterations Total number of MCMC iterations.
#' @param burnin Number of initial iterations discarded as burn-in.
#' @param thin Thinning interval applied after burn-in.
#' @param startvalue Optional numeric vector of starting values.
#' @export
CopRegBayes <- function(data, iterations = 10000, burnin = 2000, thin = 10,
                        startvalue = NULL) {
  pkg <- "endogCopulaBayes"
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Install 'endogCopulaBayes' to use the Bayesian estimator.", call. = FALSE)
  }
  bayes_fun <- getExportedValue(pkg, "CopRegBayes")
  bayes_fun(
    data = data,
    iterations = iterations,
    burnin = burnin,
    thin = thin,
    startvalue = startvalue
  )
}
