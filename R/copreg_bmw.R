#' Copula Regression: Breitung, Meyer, and Wied (2024)
#'
#' Implements the estimator proposed by Breitung, Meyer, and Wied (2024) that
#' combines classical first-stage regressions with copula-based control
#' functions.
#'
#' @inheritParams CopRegPG
#' @return An object of class `endog_copula_fit`.
#' @references Breitung, J., A. Meyer, and D. Wied (2024). Asymptotic properties
#'   of endogeneity corrections using nonlinear transformations. *The
#'   Econometrics Journal* 27 (3), 362–383.
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
    k <- nrow(fit$coefficients)
    boot_res <- vapply(
      seq_len(nboots),
      function(i) {
        boot_data <- bootstrap_sample(cleaned)
        res <- bmw_fit(boot_data, components, formula, cdf)
        res$coefficients[, 1]
      },
      numeric(k)
    )
    rownames(boot_res) <- rownames(fit$coefficients)
    boot_mat <- t(boot_res)
    colnames(boot_mat) <- rownames(fit$coefficients)
    ses <- apply(boot_res, 1, stats::sd)
    ses <- ses[rownames(fit$coefficients)]
    fit$coefficients <- cbind(Estimate = fit$coefficients[, 1], `Std. Error` = ses)
  }

  make_result(
    estimates = fit$coefficients,
    residuals = fit$residuals,
    call = call,
    method = "Breitung, Meyer, and Wied (2024)",
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
    design_info <- pg_design_no_exog(data, response, has_intercept)
    transformed <- cdf_transform(design_info$endogenous_matrix, cdf)
    colnames(transformed) <- paste0(colnames(design_info$endogenous_matrix), "_cop")

    estimation_matrix <- cbind(design_info$full_matrix, transformed)
    estimation_df <- as.data.frame(estimation_matrix)
    lm_formula <- stats::as.formula(paste(response, "~ . -", response, "- 1"))
    fit <- safe_lm(lm_formula, estimation_df)

    estimates <- stats::coef(fit)
    beta <- estimates[!grepl("_cop$", names(estimates))]
    coef_names <- intersect(names(beta), colnames(design_info$design_matrix))
    X_beta <- design_info$design_matrix[, coef_names, drop = FALSE]
    fitted <- as.numeric(X_beta %*% beta[coef_names])
    residuals <- estimation_matrix[, response] - fitted
    coeff_mat <- matrix(estimates, ncol = 1)
    rownames(coeff_mat) <- names(estimates)
    colnames(coeff_mat) <- "Estimate"
    list(coefficients = coeff_mat, residuals = residuals)
  } else {
    design_info <- pg_design_with_exog(data, original_formula, response, endog)
    design_matrix <- design_info$design_matrix
    full_matrix <- design_info$full_matrix
    regressors <- setdiff(colnames(design_matrix), endog)
    regressors_terms <- format_formula_term(regressors)

    lm_residuals <- matrix(NA_real_, nrow = nrow(design_matrix), ncol = length(endog))
    colnames(lm_residuals) <- endog

    df_full <- as.data.frame(full_matrix)
    for (j in seq_along(endog)) {
      p_var <- endog[j]
      rhs <- if (length(regressors_terms) == 0) "1" else paste(regressors_terms, collapse = " + ")
      formula_P <- stats::as.formula(paste(p_var, "~", rhs))
      lm_fit <- safe_lm(formula_P, data = df_full)
      lm_residuals[, j] <- stats::residuals(lm_fit)
    }

    transformed <- cdf_transform(lm_residuals, cdf)
    transformed <- apply(transformed, 2, stats::qnorm)
    colnames(transformed) <- paste0(endog, "_cop")

    estimation_matrix <- cbind(full_matrix, transformed)
    estimation_df <- as.data.frame(estimation_matrix)
    lm_formula <- stats::as.formula(paste(response, "~ . -", response, "- 1"))
    fit <- safe_lm(lm_formula, estimation_df)

    estimates <- stats::coef(fit)
    beta <- estimates[!grepl("_cop$", names(estimates))]
    X_beta <- design_info$design_matrix[, names(beta), drop = FALSE]
    fitted <- as.numeric(X_beta %*% beta)
    residuals <- estimation_matrix[, response] - fitted
    coeff_mat <- matrix(estimates, ncol = 1)
    rownames(coeff_mat) <- names(estimates)
    colnames(coeff_mat) <- "Estimate"
    list(coefficients = coeff_mat, residuals = residuals)
  }
}
