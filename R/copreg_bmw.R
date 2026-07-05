#' Copula Regression: Breitung, Mayer, and Wied (2024)
#'
#' Implements the estimator proposed by Breitung, Mayer, and Wied (2024) that
#' combines classical first-stage regressions with copula-based control
#' functions. With exogenous regressors, each endogenous regressor is first
#' regressed on the exogenous ones and the copula control function is built
#' from the CDF-transformed (and normal-quantile mapped) first-stage
#' residuals. Without exogenous regressors the CDF-transformed endogenous
#' regressors enter the control-function regression directly.
#'
#' @inheritParams CopRegPG
#' @inherit CopRegPG return
#' @references Breitung, J., A. Mayer, and D. Wied (2024). Asymptotic properties
#'   of endogeneity corrections using nonlinear transformations. *The
#'   Econometrics Journal* 27 (3), 362–383.
#' @examples
#' data(sim_endog)
#'
#' fit <- CopRegBMW(y ~ z_endog | x_exog + w_instr, data = sim_endog,
#'                  cdf = "resc.ecdf", nboots = 0)
#' summary(fit)
#' @export
CopRegBMW <- function(formula, data, cdf = c("kde", "ecdf", "resc.ecdf", "adj.ecdf"),
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

  fit <- bmw_fit(cleaned, components, formula, cdf)

  boot_mat <- NULL
  if (nboots > 0) {
    message("Computing bootstrap standard errors (this may take a while)...")
    expected_names <- rownames(fit$coefficients)
    fit_fun <- function(boot_data) {
      res <- bmw_fit(boot_data, components, formula, cdf)
      stats::setNames(res$coefficients[, 1], rownames(res$coefficients))
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
    method = "Breitung, Mayer, and Wied (2024)",
    cdf = cdf,
    boot_replicates = boot_mat
  )
}

bmw_fit <- function(data, components, original_formula, cdf) {
  response <- components$response
  endog <- components$endogenous
  exog <- components$exogenous
  has_intercept <- components$has_intercept

  if (length(exog) == 0) {
    # Without exogenous regressors the reference implementation skips the
    # first stage and enters the CDF-transformed endogenous regressors
    # directly (no qnorm).
    design_info <- pg_design_no_exog(data, response, has_intercept)
    transformed <- cdf_transform(design_info$endogenous_matrix, cdf)
    colnames(transformed) <- paste0(colnames(design_info$endogenous_matrix), "_cop")

    estimation_matrix <- cbind(design_info$full_matrix, transformed)
    # Includes the intercept column when present, mirroring the reference
    # residual computation.
    fitted_matrix <- design_info$design_matrix
  } else {
    design_info <- pg_design_with_exog(data, original_formula, response, endog)
    design_matrix <- design_info$design_matrix
    if (has_intercept) {
      design_matrix <- design_matrix[, -1, drop = FALSE]
    }
    regressors <- setdiff(colnames(design_matrix), endog)
    regressors_terms <- format_formula_term(regressors)

    lm_residuals <- matrix(NA_real_, nrow = nrow(design_matrix), ncol = length(endog))
    colnames(lm_residuals) <- endog

    df_full <- as.data.frame(design_info$full_matrix)
    rhs <- if (length(regressors_terms) == 0) "1" else paste(regressors_terms, collapse = " + ")
    for (j in seq_along(endog)) {
      formula_P <- stats::as.formula(paste(endog[j], "~", rhs))
      lm_fit <- safe_lm(formula_P, data = df_full)
      lm_residuals[, j] <- stats::residuals(lm_fit)
    }

    transformed <- cdf_transform(lm_residuals, cdf)
    transformed <- apply(transformed, 2, stats::qnorm)
    colnames(transformed) <- paste0(endog, "_cop")

    estimation_matrix <- cbind(design_info$full_matrix, transformed)
    fitted_matrix <- if (has_intercept) {
      cbind(`(Intercept)` = rep(1, nrow(design_matrix)), design_matrix)
    } else {
      design_matrix
    }
  }

  estimation_df <- as.data.frame(estimation_matrix)
  lm_formula <- stats::as.formula(paste(response, "~ . -", response, "- 1"))
  fit <- safe_lm(lm_formula, estimation_df)

  estimates <- stats::coef(fit)
  beta <- estimates[!grepl("_cop$", names(estimates))]
  # lm() backtick-quotes non-syntactic names such as `(Intercept)`; strip the
  # quotes so the coefficients can be matched to design-matrix columns.
  names(beta) <- sanitise_column_names(names(beta))
  coef_names <- intersect(names(beta), colnames(fitted_matrix))
  fitted <- as.numeric(fitted_matrix[, coef_names, drop = FALSE] %*% beta[coef_names])
  residuals <- estimation_matrix[, response] - fitted
  coeff_mat <- matrix(estimates, ncol = 1)
  rownames(coeff_mat) <- names(estimates)
  colnames(coeff_mat) <- "Estimate"
  list(coefficients = coeff_mat, residuals = residuals)
}
