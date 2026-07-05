# Numerical equivalence of CopRegBMW against the published reference
# implementation (CopRegBMW.R in Copula-based-endogeneity-corrections-main).
#
# The reference and the package produce identical coefficient names,
# including the backtick-quoted `(Intercept)` term that lm() creates from
# the literal "(Intercept)" column of the estimation matrix, so the name
# mapping between the two implementations is the identity. Copula control
# terms are named "<endog>_cop" in both.
#
# The reference runs its bootstrap unconditionally, so it is called with
# nboots = 1 (the smallest value pbsapply tolerates); only the point
# estimates are compared. The package is called with nboots = 0.

simulate_bmw_data <- function(n = 400) {
  set.seed(20240704)
  z1 <- stats::rnorm(n)
  z2 <- stats::runif(n)
  err <- stats::rnorm(n)
  # Non-normal endogenous regressors correlated with the error
  p1 <- exp(stats::rnorm(n) + 0.5 * err)
  p2 <- stats::rchisq(n, df = 3) + 0.7 * err
  f <- factor(sample(letters[1:3], n, replace = TRUE))
  y <- 1 + 2 * p1 - 1.5 * p2 + 0.8 * z1 + 0.5 * z2 + err
  data.frame(y = y, p1 = p1, p2 = p2, z1 = z1, z2 = z2, f = f)
}

expect_bmw_matches_reference <- function(ref_env, formula, data, cdf) {
  # Silence the reference's print() call and pbapply progress bar
  old_pbo <- pbapply::pboptions(type = "none")
  on.exit(pbapply::pboptions(old_pbo), add = TRUE)
  invisible(capture.output(
    ref <- ref_env$CopRegBMW(formula = formula, data = data, cdf = cdf, nboots = 1)
  ))
  pkg <- CopRegBMW(formula, data = data, cdf = cdf, nboots = 0)

  ref_est <- ref[[1]][, "Estimate"]
  pkg_est <- pkg$coefficients[, "Estimate"]

  expect_setequal(names(pkg_est), names(ref_est))
  expect_equal(pkg_est[names(ref_est)], ref_est, tolerance = 1e-8)
  expect_equal(as.numeric(pkg$residuals), as.numeric(ref[[2]]), tolerance = 1e-8)
}

test_that("CopRegBMW matches reference: 2 endogenous + 2 exogenous", {
  skip_if_no_reference()
  suppressMessages({
    library(dplyr)
    library(copula)
    library(ks)
    library(Matrix)
    library(nlme)
    library(pbapply)
  })
  ref_env <- load_reference_functions("CopRegBMW.R")
  dat <- simulate_bmw_data()

  for (cdf in c("resc.ecdf", "adj.ecdf", "ecdf", "kde")) {
    expect_bmw_matches_reference(ref_env, y ~ p1 + p2 | z1 + z2, dat, cdf)
  }
})

test_that("CopRegBMW matches reference: endogenous-only formula", {
  skip_if_no_reference()
  suppressMessages({
    library(dplyr)
    library(copula)
    library(ks)
    library(Matrix)
    library(nlme)
    library(pbapply)
  })
  ref_env <- load_reference_functions("CopRegBMW.R")
  dat <- simulate_bmw_data()

  # The reference supports zero exogenous regressors: it skips the first
  # stage and enters CDF(P) directly, without the qnorm mapping used in the
  # exogenous branch.
  for (cdf in c("resc.ecdf", "adj.ecdf", "ecdf")) {
    expect_bmw_matches_reference(ref_env, y ~ p1 + p2, dat, cdf)
  }
})

test_that("CopRegBMW matches reference: factor exogenous regressor", {
  skip_if_no_reference()
  suppressMessages({
    library(dplyr)
    library(copula)
    library(ks)
    library(Matrix)
    library(nlme)
    library(pbapply)
  })
  ref_env <- load_reference_functions("CopRegBMW.R")
  dat <- simulate_bmw_data()

  # Factor dummies enter both the first-stage regressions and the outer
  # control-function regression as model-matrix columns (fb, fc).
  for (cdf in c("resc.ecdf", "adj.ecdf")) {
    expect_bmw_matches_reference(ref_env, y ~ p1 + p2 | z1 + f, dat, cdf)
  }
})

test_that("CopRegBMW matches reference: intercept removed", {
  skip_if_no_reference()
  suppressMessages({
    library(dplyr)
    library(copula)
    library(ks)
    library(Matrix)
    library(nlme)
    library(pbapply)
  })
  ref_env <- load_reference_functions("CopRegBMW.R")
  dat <- simulate_bmw_data()

  expect_bmw_matches_reference(ref_env, y ~ p1 + p2 - 1 | z1 + z2, dat, "resc.ecdf")
})
