# Numerical equivalence of CopRegIMA against the published reference
# implementation (CopRegIMA.R in Copula-based-endogeneity-corrections-main).
#
# Notes on the reference algorithm (Haschka 2024) and coefficient mapping:
# * With exogenous regressors, ALL regressors (endogenous + exogenous,
#   including factor dummy columns) are CDF-transformed and mapped through
#   qnorm; the copula control for each endogenous variable is the residual
#   from a no-intercept first-stage regression of its transformed values on
#   the transformed exogenous columns.
# * WITHOUT exogenous regressors the reference adds the raw CDF values as
#   controls (NO qnorm step) -- deliberately different from CopRegPG.
# * The second stage is lm(y ~ . - y - 1) on a data.frame built from the
#   model matrix, so both implementations name the intercept coefficient
#   `(Intercept)` (with literal backticks) and the copula controls
#   <endog>_cop. Names therefore match one-to-one; we still index the
#   package estimates by the reference names to be order-independent.
# * The reference bootstrap runs pbsapply(1:nboots, ...), so nboots = 0
#   would iterate over c(1, 0); the smallest safe value is 2. Point
#   estimates are computed before bootstrapping, so they are unaffected,
#   and we compare only the "Estimate" column (package run uses nboots = 0).

suppressMessages({
  library(dplyr)
  library(pbapply)
  library(copula)
  library(ks)
  library(Matrix)
  library(nlme)
})

make_ima_data <- function(n = 400, seed = 20240704) {
  set.seed(seed)
  err <- rnorm(n)
  dat <- data.frame(
    P1 = exp(0.6 * rnorm(n) + 0.5 * err),        # lognormal, endogenous
    P2 = rchisq(n, df = 3) + 0.4 * err,          # chi-square, endogenous
    X1 = rnorm(n),
    X2 = runif(n, -1, 1),
    fac = factor(sample(c("a", "b", "c"), n, replace = TRUE))
  )
  dat$y <- 1 + 2 * dat$P1 - 1.5 * dat$P2 + 0.8 * dat$X1 - 0.6 * dat$X2 +
    0.5 * (dat$fac == "b") + 1.2 * (dat$fac == "c") + err
  dat
}

# Runs the reference (nboots = 2, smallest its bootstrap tolerates) and the
# package (nboots = 0) on the same data and compares point estimates and
# residuals. Reference progress bars / print() calls are silenced.
expect_ima_matches_reference <- function(ref_env, formula, dat, cdf) {
  old_pbo <- pbapply::pboptions(type = "none")
  on.exit(pbapply::pboptions(old_pbo), add = TRUE)

  invisible(utils::capture.output(
    ref <- ref_env$CopRegIMA(formula = formula, data = dat, cdf = cdf, nboots = 2)
  ))
  pkg <- CopRegIMA(formula, data = dat, cdf = cdf, nboots = 0)

  ref_est <- ref[[1]][, "Estimate"]
  pkg_est <- pkg$coefficients[, "Estimate"]

  expect_setequal(names(pkg_est), names(ref_est))
  expect_equal(pkg_est[names(ref_est)], ref_est, tolerance = 1e-8)
  expect_equal(
    unname(as.numeric(pkg$residuals)),
    unname(as.numeric(ref[[2]])),
    tolerance = 1e-8
  )
  invisible(pkg)
}

test_that("CopRegIMA matches reference: 2 endogenous + 2 exogenous", {
  skip_if_no_reference()
  ref_env <- load_reference_functions("CopRegIMA.R")
  dat <- make_ima_data()
  for (cdf in c("resc.ecdf", "adj.ecdf", "ecdf", "kde")) {
    expect_ima_matches_reference(ref_env, y ~ P1 + P2 | X1 + X2, dat, cdf)
  }
})

test_that("CopRegIMA matches reference: endogenous-only model", {
  skip_if_no_reference()
  ref_env <- load_reference_functions("CopRegIMA.R")
  dat <- make_ima_data(seed = 987)
  # The reference supports one-part formulas: it appends the raw CDF values
  # (no qnorm, no first stage) as control functions.
  for (cdf in c("resc.ecdf", "adj.ecdf", "ecdf")) {
    expect_ima_matches_reference(ref_env, y ~ P1 + P2, dat, cdf)
  }
})

test_that("CopRegIMA matches reference: factor exogenous regressor", {
  skip_if_no_reference()
  ref_env <- load_reference_functions("CopRegIMA.R")
  dat <- make_ima_data(seed = 555)
  # The factor expands to dummy columns in the model matrix; the reference
  # CDF-transforms and qnorm-maps those dummies alongside the continuous
  # exogenous regressors in the first stage.
  for (cdf in c("resc.ecdf", "adj.ecdf")) {
    expect_ima_matches_reference(ref_env, y ~ P1 + P2 | X1 + fac, dat, cdf)
  }
})

test_that("CopRegIMA matches reference: intercept removed", {
  skip_if_no_reference()
  ref_env <- load_reference_functions("CopRegIMA.R")
  dat <- make_ima_data(seed = 321)
  expect_ima_matches_reference(ref_env, y ~ P1 + P2 - 1 | X1 + X2, dat, "resc.ecdf")
})
