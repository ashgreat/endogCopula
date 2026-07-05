# Bootstrap fit functions ------------------------------------------------
#
# Each function computes the CopRegJAMS point estimates on a (resampled)
# data.frame and returns the named coefficient vector. They mirror the
# boots1(), boots_JAMS1() and boots_JAMS2() functions of the reference
# implementation, except that resampling and retrying degenerate draws is
# delegated to bootstrap_estimates() in utils.R, which keeps the bootstrap
# output deterministic in shape.

jams_boot_no_exog <- function(data_cleaned, formula, dependent_var,
                              has_intercept, cdf) {

  design_matrix <- model.matrix(formula, data = data_cleaned)
  YP <- t(as.matrix(design_matrix)) %*% as.matrix(design_matrix)
  if (Matrix::rankMatrix(YP)[1] != ncol(design_matrix)) {
    stop("Bootstrap design matrix is rank deficient", call. = FALSE)
  }

  if (!has_intercept) { full_matrix <- model.matrix(~ . - 1, data = data_cleaned) } else {
    full_matrix <- model.matrix(~ ., data = data_cleaned) }

  # endogenous regressor(s)
  if (has_intercept) { P1 <- design_matrix[, -1] } else { P1 <- design_matrix }
  P1 <- as.matrix(P1)
  P_star1 <- cdf_transform(P1, cdf)

  # rename columns of P_star matrix
  colnames(P_star1) <- paste0(colnames(P1), "_cop")

  # merge columns
  est_matrix <- cbind(full_matrix, P_star1)

  # control function approach
  mod1 <- lm(as.formula(paste(dependent_var, "~ . -", dependent_var, "- 1")),
             data = as.data.frame(est_matrix))

  mod1$coefficients

}

jams_boot_continuous <- function(data_cleaned, full_formula, dependent_var,
                                 independent_P_vars, has_intercept, cdf,
                                 design_matrix1) {

  design_matrix <- model.matrix(full_formula, data = data_cleaned)
  YP <- t(as.matrix(design_matrix)) %*% as.matrix(design_matrix)
  rank_A <- Matrix::rankMatrix(YP)[1]
  if (!(rank_A == ncol(design_matrix) &&
        all(colnames(design_matrix1) %in% colnames(design_matrix)))) {
    stop("Bootstrap design matrix is rank deficient", call. = FALSE)
  }

  full_matrix <- cbind(data_cleaned[[all.vars(full_formula)[1]]], design_matrix)
  colnames(full_matrix)[1] <- all.vars(full_formula)[1]

  if (has_intercept) { P1 <- design_matrix[, -1] } else { P1 <- design_matrix }
  P_star1 <- cdf_transform(P1, cdf)

  # rename columns of P_star matrix
  colnames(P_star1) <- paste0(colnames(P1))

  # generate copula-correction terms
  P_star2 <- P_star1 %*% solve(cor(apply(P_star1, 2, qnorm)))
  P_star2 <- P_star2[, c(1:length(independent_P_vars))]
  P_star2 <- as.matrix(P_star2)

  # rename columns of P_star matrix
  colnames(P_star2) <- paste0(independent_P_vars, "_cop")

  # merge columns
  est_matrix <- cbind(full_matrix, P_star2)

  # control function approach
  mod1 <- lm(as.formula(paste(dependent_var, "~ . -", dependent_var, "- 1")),
             data = as.data.frame(est_matrix))

  mod1$coefficients

}

jams_boot_factor <- function(data_cleaned, full_formula, dependent_var,
                             independent_P_vars, independent_X_vars, cdf,
                             factor_vars, design_matrix1) {

  design_matrix <- model.matrix(full_formula, data = data_cleaned)
  YP <- t(as.matrix(design_matrix)) %*% as.matrix(design_matrix)
  rank_A <- Matrix::rankMatrix(YP)[1]
  if (!(rank_A == ncol(design_matrix) &&
        all(colnames(design_matrix1) %in% colnames(design_matrix)))) {
    stop("Bootstrap design matrix is rank deficient", call. = FALSE)
  }

  # all regressor(s)
  data_cleaned_with_P_star2 <- data_cleaned
  full_matrix_with_P_Pstar2 <- data_cleaned_with_P_star2

  # Loop over each factor variable
  for (var in factor_vars) {

    levels_var <- levels(data_cleaned[[var]])  # Get factor levels

    for (lvl in levels_var) {

      subdat1 <- subset(data_cleaned, data_cleaned[[var]] == lvl)

      if (nrow(subdat1) > 3 &&
          any(sapply(subdat1[independent_P_vars], function(col) length(unique(col)) > 1))) {

        subdat2 <- subdat1[, setdiff(c(independent_P_vars, independent_X_vars), factor_vars)]
        subdat2 <- as.matrix(subdat2)
        subdat2 <- subdat2[, apply(subdat2, 2, function(x) length(unique(x)) > 1), drop = FALSE]

        P_star1 <- cdf_transform(subdat2, cdf)

        # generate copula-correction terms
        is_issue <- tryCatch(
          { solve(cor(apply(P_star1, 2, qnorm))) ; FALSE },
          error = function(e) TRUE,
          warning = function(w) TRUE
        )

        if (is_issue) {

          P_star2 <- matrix(data = 0, nrow = nrow(P_star1), ncol = ncol(P_star1))
          P_star2 <- as.matrix((P_star2))

        } else {

          P_star2 <- P_star1 %*% solve(cor(apply(P_star1, 2, qnorm)))

          if (length(independent_P_vars) > ncol(P_star2)) {
            P_star2 <- as.matrix((P_star2))
          } else {
            P_star2 <- P_star2[, c(1:length(independent_P_vars))]
            P_star2 <- as.matrix((P_star2))
          }

        }

        # rename columns of P_star matrix
        colnames(P_star2) <- paste0(independent_P_vars[c(1:ncol(P_star2))],
                                    "_", var, "_", lvl, "_cop")

        # merge data
        P_star2_full <- matrix(0, nrow = nrow(data_cleaned_with_P_star2), ncol = ncol(P_star2))
        P_star2_full[data_cleaned_with_P_star2[[var]] == lvl, ] <- P_star2
        P_star2_df <- as.data.frame(P_star2_full)
        colnames(P_star2_df) <- colnames(P_star2)

        full_matrix_with_P_Pstar2 <- cbind(full_matrix_with_P_Pstar2, P_star2_df)

      }

    }

  }

  # control function approach
  mod1 <- lm(as.formula(paste(dependent_var, "~ . -",
                              dependent_var)),
             data = full_matrix_with_P_Pstar2)

  Estimates <- mod1$coefficients
  Estimates <- Estimates[!is.na(Estimates)]

  # Align coefficient names with the initial fit (the reference strips
  # "as.factor()" wrappers and backticks before comparing replicate names)
  names(Estimates) <- jams_normalise_names(names(Estimates))

  Estimates

}

jams_normalise_names <- function(x) {
  x <- gsub("as.factor\\((.*?)\\)", "\\1", x)
  gsub("`", "", x)
}

# Liengaard et al. (2024) estimator
CopRegJAMS_impl <- function(formula, data, cdf = c("kde", "ecdf", "resc.ecdf", "adj.ecdf"),
                       nboots = 199) {

  ################################################################################

  # check arguments
  if (rlang::is_formula(formula) == FALSE) { stop("Argument formula is not a formula object", call. = FALSE) }
  if (is.data.frame(data) == FALSE) { stop("Argument data is not a data.frame object", call. = FALSE) }
  if (inherits(data, "data.table")) { stop("Argument data should not be a data.table.", call = FALSE) }
  if (is.numeric(nboots) == FALSE) { stop("Argument nboots is not numeric", call. = FALSE) }

  # seperate endogenous and exogenous regressor(s)
  f1 <- nlme::splitFormula(formula, sep = "|")

  # check if intercept is removed
  has_intercept <- attr(terms(f1[[1]]), "intercept") == 1

  ##############################################################################

  # no exogenous regressors
  if (length(f1) == 1) {

    # no exogenous regressors
    f1P <- f1[[1]]

    # Check if all variables exist in the data
    variables <- all.vars(f1P)
    missing_vars <- setdiff(variables, names(data))
    if(length(missing_vars) > 0) {

      stop(paste("The following variables are missing in the data:",
                 paste(missing_vars, collapse=", ")))

    }

    # Check if all variables are numeric and non-constant
    numeric_vars <- sapply(data[variables], is.numeric)
    constant_vars <- sapply(data[variables], function(x) length(unique(x)) == 1)

    if (!all(numeric_vars)) {
      stop("The following are dummy variables and cannot be endogenous: ", paste(variables[!numeric_vars], collapse = ", "))
    }

    if (any(constant_vars)) {
      stop("The following variables are constant: ", paste(variables[constant_vars], collapse = ", "))
    }

    dependent_var <- all.vars(formula)[1]
    independent_vars <- all.vars(f1P)

    data_cleaned <- data %>%
      dplyr::select(all_of(c(dependent_var, independent_vars))) %>%
      na.omit() %>%
      as.data.frame()

    # Check if design matrix has full column rank
    if (!has_intercept) {

      # Use model.matrix to handle categorical variables
      full_matrix <- model.matrix(~ . - 1, data = data_cleaned)
      design_matrix <- full_matrix[, colnames(full_matrix) != dependent_var, drop = FALSE]

      # Check rank of the full and the design matrix
      rank_A <- rankMatrix(full_matrix)[1]
      rank_B <- rankMatrix(design_matrix)[1]

      if (rank_A != ncol(full_matrix)) {
        stop("Perfect fit")
      }
      if (rank_B != ncol(design_matrix)) {
        stop("Design matrix is rank deficient")
      }

    } else {

      # Add an intercept column to the design matrix
      full_matrix <- model.matrix(~ ., data = data_cleaned)
      design_matrix <- full_matrix[, colnames(full_matrix) != dependent_var, drop = FALSE]

      # Check rank of the full and the design matrix
      rank_A <- rankMatrix(full_matrix)[1]
      rank_B <- rankMatrix(design_matrix)[1]

      if (rank_A != ncol(full_matrix)) {
        stop("Perfect fit")
      }
      if (rank_B != ncol(design_matrix)) {
        stop("Design matrix is rank deficient")
      }

    }

    # endogenous regressor(s)
    if (has_intercept) { P1 <- design_matrix[, -1] } else { P1 <- design_matrix }
    P1 <- as.matrix(P1)
    P_star1 <- cdf_transform(P1, cdf)

    # rename columns of P_star matrix
    colnames(P_star1) <- paste0(colnames(P1), "_cop")

    # merge columns
    est_matrix <- cbind(full_matrix, P_star1)

    # control function approach
    mod1 <- lm(as.formula(paste(dependent_var, "~ . -", dependent_var, "- 1")),
               data = as.data.frame(est_matrix))

    Estimates <- mod1$coefficients

    # obtain residuals
    beta <- Estimates[!grepl("_cop", names(Estimates))]
    fitted_values <- design_matrix %*% beta
    residuals_manual <- est_matrix[, dependent_var] - fitted_values

    # Bootstrapping
    if (nboots > 0) {

      message("Estimation done. Calculating bootstrap standard errors")
      boot_mat <- bootstrap_estimates(
        data = data_cleaned,
        fit_fun = function(d) {
          jams_boot_no_exog(d, formula = formula, dependent_var = dependent_var,
                            has_intercept = has_intercept, cdf = cdf)
        },
        nboots = nboots,
        expected_names = names(Estimates)
      )
      ses <- apply(boot_mat, 2, sd)

      Estimates1 <- cbind(Estimates, ses)
      colnames(Estimates1) <- c("Estimate", "Std. Error")

    } else {

      boot_mat <- NULL
      Estimates1 <- matrix(Estimates, ncol = 1,
                           dimnames = list(names(Estimates), "Estimate"))

    }

  } else {

    ############################################################################

    f1P <- f1[[1]]
    f1X <- f1[[2]]

    variables <- all.vars(f1P)
    missing_vars <- setdiff(variables, names(data))
    if(length(missing_vars) > 0) {

      stop(paste("The following endogenous variables are missing in the data:",
                 paste(missing_vars, collapse=", ")))

    }

    variables <- all.vars(f1X)
    missing_vars <- setdiff(variables, names(data))
    if(length(missing_vars) > 0) {

      stop(paste("The following exogenous variables are missing in the data:",
                 paste(missing_vars, collapse=", ")))

    }

    dependent_var <- all.vars(formula)[1]
    independent_P_vars <- all.vars(f1P)
    independent_X_vars <- all.vars(f1X)

    # Check if all variables are numeric and non-constant
    numeric_vars <- sapply(data[independent_P_vars], is.numeric)
    constant_vars <- sapply(data[variables], function(x) length(unique(x)) == 1)

    if (!all(numeric_vars)) {
      stop("Only continuous variables can be endogenous. The following variables are not numeric: ", paste(independent_P_vars[!numeric_vars], collapse = ", "))
    }

    if (any(constant_vars)) {
      stop("The following variables are constant: ",
           paste(variables[constant_vars], collapse = ", "))
    }

    # all variables
    data_cleaned <- data %>%
      dplyr::select(all_of(c(dependent_var, independent_P_vars,
                             independent_X_vars))) %>% na.omit() %>%
      as.data.frame()

    # Check if design matrix has full column rank
    full_formula <- as.formula(paste(all.vars(formula)[1], "~",
                                     gsub("\\|", "+", as.character(formula)[3])))
    full_matrix <- model.matrix(full_formula, data = data_cleaned)
    full_matrix <- cbind(DependentVar = data_cleaned[[dependent_var]], full_matrix)

    colnames(full_matrix)[1] <- dependent_var
    design_matrix <- full_matrix[, colnames(full_matrix) != dependent_var,
                                 drop = FALSE]

    rank_A <- Matrix::rankMatrix(full_matrix)[1]
    rank_B <- Matrix::rankMatrix(design_matrix)[1]

    if (rank_A != ncol(full_matrix)) {
      stop("Perfect fit")
    }
    if (rank_B != ncol(design_matrix)) {
      stop("Design matrix is rank deficient")
    }

    ############################################################################

    if (sum(sapply(data_cleaned, is.factor)) == 0 &
        !grepl("as.factor\\(", paste(deparse(full_formula), collapse = " "))) {

      # all regressor(s)
      if (has_intercept) { P1 <- design_matrix[, -1] } else { P1 <- design_matrix }
      P_star1 <- cdf_transform(P1, cdf)

      # rename columns of P_star matrix
      colnames(P_star1) <- paste0(colnames(P1))

      # generate copula-correction terms
      P_star2 <- P_star1 %*% solve(cor(apply(P_star1, 2, qnorm)))
      P_star2 <- P_star2[, c(1:length(independent_P_vars))]
      P_star2 <- as.matrix(P_star2)

      # rename columns of P_star matrix
      colnames(P_star2) <- paste0(independent_P_vars, "_cop")

      # merge columns
      est_matrix <- cbind(full_matrix, P_star2)

      # control function approach
      mod1 <- lm(as.formula(paste(dependent_var, "~ . -",
                                  dependent_var, "- 1")),
                 data = as.data.frame(est_matrix))

      Estimates <- mod1$coefficients

      # obtain residuals
      beta <- Estimates[!grepl("_cop", names(Estimates))]
      fitted_values <- design_matrix %*% beta
      residuals_manual <- est_matrix[, dependent_var] - fitted_values

      # Bootstrapping
      if (nboots > 0) {

        message("Estimation done. Calculating bootstrap standard errors")
        boot_mat <- bootstrap_estimates(
          data = data_cleaned,
          fit_fun = function(d) {
            jams_boot_continuous(d, full_formula = full_formula,
                                 dependent_var = dependent_var,
                                 independent_P_vars = independent_P_vars,
                                 has_intercept = has_intercept, cdf = cdf,
                                 design_matrix1 = design_matrix)
          },
          nboots = nboots,
          expected_names = names(Estimates)
        )
        ses <- apply(boot_mat, 2, sd)

        Estimates1 <- cbind(Estimates, ses)
        colnames(Estimates1) <- c("Estimate", "Std. Error")

      } else {

        boot_mat <- NULL
        Estimates1 <- matrix(Estimates, ncol = 1,
                             dimnames = list(names(Estimates), "Estimate"))

      }

    } else {

      ##########################################################################

      # extract factor variables
      formula_str <- paste(deparse(full_formula), collapse = " ")
      factor_vars_formula <- regmatches(formula_str, gregexpr("as.factor\\((.*?)\\)", formula_str))[[1]]
      factor_vars_formula <- gsub("as.factor\\(|\\)", "", factor_vars_formula)

      factor_vars_data <- names(data_cleaned)[sapply(data_cleaned, is.factor)]
      factor_vars <- unique(c(factor_vars_formula, factor_vars_data))

      # all regressor(s)
      data_cleaned_with_P_star2 <- data_cleaned
      full_matrix_with_P_Pstar2 <- full_matrix
      data_cleaned[factor_vars] <- lapply(data_cleaned[factor_vars], as.factor)

      # Loop over each factor variable
      for (var in factor_vars) {

        levels_var <- levels(data_cleaned[[var]])  # Get factor levels

        for (lvl in levels_var) {

          subdat1 <- subset(data_cleaned, data_cleaned[[var]] == lvl)

          if (nrow(subdat1) > 3 &&
              any(sapply(subdat1[independent_P_vars], function(col) length(unique(col)) > 1))) {

            subdat2 <- subdat1[, setdiff(c(independent_P_vars, independent_X_vars), factor_vars)]
            subdat2 <- as.matrix(subdat2)
            subdat2 <- subdat2[, apply(subdat2, 2, function(x) length(unique(x)) > 1), drop = FALSE]

            P_star1 <- cdf_transform(subdat2, cdf)

            # generate copula-correction terms
            is_issue <- tryCatch(
              { solve(cor(apply(P_star1, 2, qnorm))) ; FALSE },
              error = function(e) TRUE,
              warning = function(w) TRUE
            )

            if (is_issue) {

              P_star2 <- matrix(data = 0, nrow = nrow(P_star1), ncol = ncol(P_star1))
              P_star2 <- as.matrix((P_star2))

            } else {

              P_star2 <- P_star1 %*% solve(cor(apply(P_star1, 2, qnorm)))

              if (length(independent_P_vars) > ncol(P_star2)) {
                P_star2 <- as.matrix((P_star2))
              } else {
                P_star2 <- P_star2[, c(1:length(independent_P_vars))]
                P_star2 <- as.matrix((P_star2))
              }

            }

            # rename columns of P_star matrix
            colnames(P_star2) <- paste0(independent_P_vars[c(1:ncol(P_star2))],
                                        "_", var, "_", lvl, "_cop")

            # merge data
            P_star2_full <- matrix(0, nrow = nrow(data_cleaned_with_P_star2), ncol = ncol(P_star2))
            P_star2_full[data_cleaned_with_P_star2[[var]] == lvl, ] <- P_star2
            P_star2_df <- as.data.frame(P_star2_full)
            colnames(P_star2_df) <- colnames(P_star2)

            full_matrix_with_P_Pstar2 <- cbind(full_matrix_with_P_Pstar2, P_star2_df)

          }

        }

      }

      # control function approach
      mod1 <- lm(as.formula(paste(dependent_var, "~ . -",
                                  dependent_var)),
                 data = full_matrix_with_P_Pstar2)

      Estimates <- mod1$coefficients
      Estimates <- Estimates[!is.na(Estimates)]

      # obtain residuals
      beta <- Estimates[!grepl("_cop", names(Estimates))]
      fitted_values <- design_matrix %*% beta
      residuals_manual <- data_cleaned_with_P_star2[, dependent_var] - fitted_values

      # Bootstrapping
      if (nboots > 0) {

        message("Estimation done. Calculating bootstrap standard errors")
        boot_mat <- bootstrap_estimates(
          data = data_cleaned,
          fit_fun = function(d) {
            jams_boot_factor(d, full_formula = full_formula,
                             dependent_var = dependent_var,
                             independent_P_vars = independent_P_vars,
                             independent_X_vars = independent_X_vars,
                             cdf = cdf, factor_vars = factor_vars,
                             design_matrix1 = design_matrix)
          },
          nboots = nboots,
          expected_names = jams_normalise_names(names(Estimates))
        )
        ses <- apply(boot_mat, 2, sd)

        Estimates1 <- cbind(Estimates, ses)
        colnames(Estimates1) <- c("Estimate", "Std. Error")

      } else {

        boot_mat <- NULL
        Estimates1 <- matrix(Estimates, ncol = 1,
                             dimnames = list(names(Estimates), "Estimate"))

      }

    }

  }

  return(list(Estimates1, residuals_manual, boot_mat))

}


# This function implements the copula-based endogeneity correction by
# Liengaard et al. (2024) using the least-squares-based correction function
# approach.
#
# formula = depvar ~ endog_var1 + endog_var2 + ... | exog_var1 + as.factor(exog_var2) + ...
#
# data = as.data.frame(datset)
#
# cdf = c("kde", "ecdf", "resc.ecdf", "adj.ecdf")
# kde is the integral of a density estimator used in Park & Gupta (2012) and Haschka (2022)
# ecdf is the empirical cumulative distribution function (ecdf) with replaced boundary proposed by Becker et al. (2022)
# resc.ecdf is a rescaled ecdf proposed by Qian et al. (2024)
# adj.ecdf is an adjusted ecdf proposed by Liengaard (2024)
#
# CopRegJAMS returns a list of legth 2. First entry are estimates with standard
# errors. Second entry are residuals.


### REFERENCES

# Becker, J.-M., D. Proksch, and C. M. Ringle (2021). Revisiting Gaussian
# copulas to handle endogenous regressors. Journal of the Academy of Marketing
# Science 50, 46–66.
#
# Haschka, R. E (2022). Handling endogenous regressors using copulas: A
# generalisation to linear panel models with fixed effects and correlated
# regressors. Journal of Marketing Research 59(4), 860–881.
#
# Liengaard, B. D., J.-M. Becker, M. Bennedsen, P. Heiler, L. N. Taylor, and
# C. M. Ringle (2024). Dealing with regression models’ endogeneity by means of
# an adjusted estimator for the Gaussian copula approach. Journal of the
# Academy of Marketing Science, 1–21.
#
# Park, S. and S. Gupta (2012). Handling endogenous regressors by joint
# estimation using copulas. Marketing Science 31 (4), 567–586.
#
# Qian, Y., A. Koschmann, and H. Xie (2024). A practical guide to endogeneity
# correction using copulas. NBER Working Paper.
# https://www.nber.org/system/files/workingpapers/w32231/w32231.pdf.


#' Copula Regression: Liengaard et al. (2025)
#'
#' Implements the adjusted estimator for Gaussian copula endogeneity correction
#' proposed by Liengaard et al. (2025). The implementation closely follows the
#' original research code and supports factor-specific transformations.
#'
#' @inheritParams CopRegPG
#' @inherit CopRegPG return
#' @references Liengaard, B. D., J.-M. Becker, M. Bennedsen, P. Heiler, L. N.
#'   Taylor, and C. M. Ringle (2025). Dealing with regression models'
#'   endogeneity by means of an adjusted estimator for the Gaussian copula
#'   approach. *Journal of the Academy of Marketing Science* 53, 279–299.
#' @examples
#' data(sim_endog)
#'
#' fit <- CopRegJAMS(y ~ z_endog | x_exog + w_instr, data = sim_endog,
#'                   cdf = "adj.ecdf", nboots = 0)
#' summary(fit)
#' @importFrom dplyr select all_of
#' @importFrom magrittr %>%
#' @importFrom stats as.formula cor ecdf lm model.matrix na.omit predict qnorm sd terms
#' @importFrom Matrix rankMatrix
#' @rdname CopRegJAMS
#' @export
CopRegJAMS <- function(formula, data, cdf = c("kde", "ecdf", "resc.ecdf", "adj.ecdf"),
                       nboots = 199) {
  call <- match.call()
  cdf <- match.arg(cdf)
  raw <- CopRegJAMS_impl(formula = formula, data = data, cdf = cdf, nboots = nboots)
  if (!is.list(raw) || length(raw) < 2) {
    stop("Internal implementation did not return the expected output.", call. = FALSE)
  }
  estimates <- as.matrix(raw[[1]])
  residuals <- raw[[2]]
  boot_mat <- if (length(raw) >= 3) raw[[3]] else NULL
  make_result(
    estimates = estimates,
    residuals = residuals,
    call = call,
    method = "Liengaard et al. (2025)",
    cdf = cdf,
    boot_replicates = boot_mat
  )
}
