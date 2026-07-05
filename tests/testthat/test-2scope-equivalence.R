# Numerical equivalence of CopReg2sCOPE against the published reference
# implementation in Copula-based-endogeneity-corrections-main/CopReg2sCOPE.R.
#
# Coefficient-name mapping: both implementations estimate the second stage via
# lm(y ~ . - y - 1) on as.data.frame(cbind(y, design_matrix, cop_terms)) with
# identically constructed columns, so coefficient names coincide exactly --
# including backtick-quoted non-syntactic names such as `(Intercept)` and
# factor dummy columns. The single exception is the single-endogenous,
# no-exogenous case: the reference drops the intercept column without
# drop = FALSE (design_matrix[, -1]), losing the column name, so its copula
# term is labelled "_cop" instead of "<var>_cop". That case is compared by
# position (unnamed).
#
# The reference bootstrap (pbsapply over boots1/boots_2sCOPE) only affects the
# "Std. Error" column; the "Estimate" column and the residuals are computed
# before bootstrapping. nboots = 2 is the smallest value the reference code
# tolerates (nboots = 1 makes pbsapply return a vector, breaking
# apply(trapped, 1, sd)), so the reference is run with nboots = 2 and only
# point estimates and residuals are compared.

skip_if_no_reference()

suppressMessages({
  library(dplyr)
  library(copula)
  library(ks)
  library(Matrix)
})

ref_env <- load_reference_functions("CopReg2sCOPE.R")
# The reference calls pbapply::pbsapply unqualified; provide a plain sapply
# stand-in so the tests do not depend on pbapply being installed.
assign("pbsapply", function(X, FUN, ...) sapply(X, FUN, ...), envir = ref_env)

run_reference_2scope <- function(...) {
  res <- NULL
  # The reference prints a progress message; keep test output clean.
  invisible(utils::capture.output(res <- ref_env$CopReg2sCOPE(...)))
  res
}

simulate_scope_data <- function(n = 400, seed = 20240704) {
  set.seed(seed)
  common <- stats::rnorm(n)
  P1 <- exp(0.5 * common + stats::rnorm(n, sd = 0.8))
  P2 <- stats::rchisq(n, df = 3) + 0.4 * common
  X1 <- stats::rnorm(n)
  X2 <- stats::runif(n)
  Fct <- factor(sample(c("a", "b", "c"), n, replace = TRUE))
  eps <- 0.9 * common + stats::rnorm(n, sd = 0.5)
  y <- 1 + 0.8 * P1 - 0.5 * P2 + 0.6 * X1 - 0.4 * X2 +
    0.3 * (Fct == "b") - 0.2 * (Fct == "c") + eps
  data.frame(y = y, P1 = P1, P2 = P2, X1 = X1, X2 = X2, Fct = Fct)
}

expect_scope_equivalent <- function(formula, data, cdf) {
  fit <- CopReg2sCOPE(formula, data = data, cdf = cdf, nboots = 0)
  set.seed(1)
  ref <- run_reference_2scope(formula = formula, data = data, cdf = cdf, nboots = 2)
  expect_equal(
    fit$coefficients[, "Estimate"],
    ref[[1]][, "Estimate"],
    tolerance = 1e-8
  )
  expect_equal(
    as.numeric(fit$residuals),
    as.numeric(ref[[2]]),
    tolerance = 1e-8
  )
}

test_that("CopReg2sCOPE matches reference: 2 endogenous + 2 exogenous", {
  dat <- simulate_scope_data()
  for (cdf in c("resc.ecdf", "adj.ecdf", "ecdf", "kde")) {
    expect_scope_equivalent(y ~ P1 + P2 | X1 + X2, dat, cdf)
  }
})

test_that("CopReg2sCOPE matches reference: endogenous-only", {
  dat <- simulate_scope_data()
  for (cdf in c("resc.ecdf", "adj.ecdf", "ecdf", "kde")) {
    expect_scope_equivalent(y ~ P1 + P2, dat, cdf)
  }
})

test_that("CopReg2sCOPE matches reference: single endogenous, no exogenous", {
  # Reference names the copula term "_cop" here (colname lost by
  # design_matrix[, -1]); the package keeps "P1_cop". Values are identical,
  # so compare by position.
  dat <- simulate_scope_data()
  for (cdf in c("resc.ecdf", "adj.ecdf", "ecdf")) {
    fit <- CopReg2sCOPE(y ~ P1, data = dat, cdf = cdf, nboots = 0)
    set.seed(1)
    ref <- run_reference_2scope(formula = y ~ P1, data = dat, cdf = cdf, nboots = 2)
    expect_equal(
      unname(fit$coefficients[, "Estimate"]),
      unname(ref[[1]][, "Estimate"]),
      tolerance = 1e-8
    )
    expect_equal(
      as.numeric(fit$residuals),
      as.numeric(ref[[2]]),
      tolerance = 1e-8
    )
  }
})

test_that("CopReg2sCOPE matches reference: factor exogenous regressor", {
  dat <- simulate_scope_data()
  for (cdf in c("resc.ecdf", "adj.ecdf", "ecdf")) {
    expect_scope_equivalent(y ~ P1 + P2 | X1 + Fct, dat, cdf)
  }
})
