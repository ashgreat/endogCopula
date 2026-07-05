#' Copula Regression: Park and Gupta (2012)
#'
#' Implements the Gaussian copula approach proposed by Park and Gupta (2012)
#' for correcting endogeneity in linear models. Bootstrap standard errors are
#' obtained via resampling.
#'
#' @param formula A two-part formula of the form
#'   `y ~ endog1 + endog2 | exog1 + exog2`, where the left-hand side lists the
#'   endogenous regressors and the right-hand side after the bar lists
#'   exogenous regressors. Use `- 1` to remove the intercept.
#' @param data A `data.frame` containing the variables referenced in `formula`.
#' @param cdf Character string specifying the cumulative distribution function
#'   estimator. One of `"kde"`, `"ecdf"`, `"resc.ecdf"`, or `"adj.ecdf"`.
#' @param nboots Number of bootstrap replications used to compute standard
#'   errors. Set to zero to skip bootstrapping.
#' @return An object of class `endog_copula_fit`: a list with components
#'   \describe{
#'     \item{`call`}{The matched call.}
#'     \item{`method`}{Character string identifying the estimator.}
#'     \item{`cdf`}{The CDF estimator used for the copula transformation.}
#'     \item{`coefficients`}{Matrix of coefficient estimates with an
#'       `Estimate` column and, when `nboots > 0`, a `Std. Error` column
#'       computed from the bootstrap replicates.}
#'     \item{`residuals`}{Numeric vector of residuals from the structural
#'       part of the model (copula control terms excluded).}
#'     \item{`bootstrap`}{Matrix of bootstrap coefficient replicates, or
#'       `NULL` when `nboots = 0`.}
#'   }
#'   Objects of this class support [print()], [summary()], [coef()],
#'   [residuals()], and [confint()].
#' @details Point estimates are numerically equivalent to the published
#'   reference implementation. Bootstrap standard errors deliberately differ:
#'   the reference code omits the inverse-normal (`qnorm`) transformation of
#'   the copula terms inside bootstrap replicates even though its main fit
#'   applies it. Here every replicate refits the same estimator as the main
#'   fit, so the standard errors describe the estimator actually reported.
#' @references Park, S. and S. Gupta (2012). Handling endogenous regressors by
#'   joint estimation using copulas. *Marketing Science* 31 (4), 567–586.
#' @examples
#' data(sim_endog)
#'
#' fit <- CopRegPG(y ~ z_endog | x_exog + w_instr, data = sim_endog,
#'                 cdf = "resc.ecdf", nboots = 0)
#' summary(fit)
#'
#' if (requireNamespace("ks", quietly = TRUE)) {
#'   fit_kde <- CopRegPG(y ~ z_endog | x_exog + w_instr, data = sim_endog,
#'                       cdf = "kde", nboots = 0)
#'   coef(fit_kde)
#' }
#' @export
CopRegPG <- function(formula, data, cdf = c("kde", "ecdf", "resc.ecdf", "adj.ecdf"),
                     nboots = 199) {
  call <- match.call()
  cdf <- match.arg(cdf)
  if (!is.numeric(nboots) || length(nboots) != 1 || nboots < 0) {
    stop("Argument 'nboots' must be a single non-negative number.", call. = FALSE)
  }
  components <- extract_components(formula, data)
  if (inherits(data, "data.table")) {
    stop("Convert 'data.table' objects to 'data.frame' before estimation.", call. = FALSE)
  }
  vars_needed <- unique(c(components$response, components$endogenous, components$exogenous))
  cleaned <- select_complete_cases(as.data.frame(data), vars_needed)
  validate_numeric_nonconstant(cleaned, components$endogenous)

  fit <- pg_fit(cleaned, components, formula, cdf)

  boot_mat <- NULL
  if (nboots > 0) {
    message("Computing bootstrap standard errors (this may take a while)...")
    expected_names <- rownames(fit$coefficients)
    # Deliberate divergence from the reference: boots1/boots_PG in
    # CopRegPG.R omit the qnorm() step in bootstrap replicates even though
    # the main fit applies it (the other reference estimators apply qnorm()
    # in both). Each replicate here refits the same estimator as the main
    # fit, so standard errors describe the estimator actually reported.
    fit_fun <- function(boot_data) {
      pg_fit(boot_data, components, formula, cdf)$coefficients[, 1]
    }
    boot_mat <- bootstrap_estimates(cleaned, fit_fun, nboots, expected_names)
    ses <- apply(boot_mat, 2, stats::sd)
    fit$coefficients <- cbind(
      Estimate = fit$coefficients[, 1],
      `Std. Error` = ses[expected_names]
    )
  }

  make_result(
    estimates = fit$coefficients,
    residuals = fit$residuals,
    call = call,
    method = "Park and Gupta (2012)",
    cdf = cdf,
    boot_replicates = boot_mat
  )
}

pg_fit <- function(data, components, original_formula, cdf) {
  response <- components$response
  endog <- components$endogenous
  exog <- components$exogenous
  has_intercept <- components$has_intercept

  if (length(exog) == 0) {
    design_info <- pg_design_no_exog(data, response, has_intercept)
    p_matrix <- design_info$endogenous_matrix
  } else {
    design_info <- pg_design_with_exog(data, original_formula, response, endog)
    p_matrix <- design_info$design_matrix[, endog, drop = FALSE]
  }

  transformed <- cdf_transform(p_matrix, cdf)
  transformed <- apply(transformed, 2, stats::qnorm)
  colnames(transformed) <- paste0(colnames(p_matrix), "_cop")

  estimation_matrix <- cbind(design_info$full_matrix, transformed)
  estimation_df <- as.data.frame(estimation_matrix)
  lm_formula <- stats::as.formula(paste(response, "~ . -", response, "- 1"))
  fit <- safe_lm(lm_formula, estimation_df)

  estimates <- stats::coef(fit)
  # Residuals follow the reference implementation: drop the copula control
  # terms and multiply the remaining coefficients into the design matrix by
  # position. lm() backtick-quotes non-syntactic column names (e.g.
  # "(Intercept)" becomes "`(Intercept)`"), so a name-based match would
  # silently drop the intercept from the fitted values.
  beta <- estimates[!grepl("_cop", names(estimates))]
  fitted <- as.numeric(design_info$design_matrix %*% beta)
  residuals <- as.numeric(estimation_matrix[, response] - fitted)
  coeff_mat <- matrix(estimates, ncol = 1)
  rownames(coeff_mat) <- names(estimates)
  colnames(coeff_mat) <- "Estimate"
  list(coefficients = coeff_mat, residuals = residuals)
}

pg_design_no_exog <- function(data, response, has_intercept) {
  if (has_intercept) {
    full_matrix <- stats::model.matrix(~ ., data = data)
  } else {
    full_matrix <- stats::model.matrix(~ . - 1, data = data)
  }
  design_matrix <- full_matrix[, colnames(full_matrix) != response, drop = FALSE]
  check_full_rank(full_matrix, "Perfect fit detected: design matrix has rank deficiency.")
  check_full_rank(design_matrix, "Design matrix is rank deficient.")
  endogenous_matrix <- if (has_intercept) {
    design_matrix[, -1, drop = FALSE]
  } else {
    design_matrix
  }
  list(
    full_matrix = full_matrix,
    design_matrix = design_matrix,
    endogenous_matrix = endogenous_matrix
  )
}

pg_design_with_exog <- function(data, original_formula, response, endogenous) {
  # deparse1() keeps long formulas intact, whereas as.character(formula)[3]
  # truncates the right-hand side at the first deparse line.
  rhs <- gsub("|", "+", deparse1(original_formula[[3]]), fixed = TRUE)
  full_formula <- stats::as.formula(paste(response, "~", rhs))
  full_matrix <- stats::model.matrix(full_formula, data = data)
  full_matrix <- cbind(DependentVar = data[[response]], full_matrix)
  colnames(full_matrix)[1] <- response
  design_matrix <- full_matrix[, colnames(full_matrix) != response, drop = FALSE]
  check_full_rank(full_matrix, "Perfect fit detected: design matrix has rank deficiency.")
  check_full_rank(design_matrix, "Design matrix is rank deficient.")
  list(
    full_matrix = full_matrix,
    design_matrix = design_matrix,
    endogenous = endogenous
  )
}
