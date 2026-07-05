# Numerical equivalence of CopRegJAMS against the published reference
# implementation (CopRegJAMS.R in Copula-based-endogeneity-corrections-main).
#
# The reference has no way to switch the bootstrap off (nboots must be >= 1
# because pbsapply(1:nboots, ...) is always evaluated), but the point
# estimates are computed before any resampling takes place. We therefore run
# the reference with nboots = 1 (its smallest tolerable value; the resulting
# standard errors are NA and ignored) and compare ONLY the "Estimate" column
# against the package fit with nboots = 0.
#
# Coefficient-name mapping: both implementations build the estimation data
# frames identically, so coefficient names match exactly:
#   * continuous paths: "X.Intercept." (as.data.frame() mangles the
#     "(Intercept)" column of the model matrix), regressor names, and
#     "<endog>_cop" for the correction terms;
#   * factor path: "(Intercept)", regressor names, backticked factor dummies
#     such as "`as.factor(w)1`", and "<endog>_<var>_<level>_cop" terms.

jams_ref_env <- NULL

load_jams_reference <- function() {
  if (is.null(jams_ref_env)) {
    suppressMessages({
      library(dplyr)
      library(Matrix)
      library(copula)
      library(ks)
      library(nlme)
      library(pbapply)
    })
    jams_ref_env <<- load_reference_functions("CopRegJAMS.R")
  }
  jams_ref_env
}

# Run the reference estimator, swallowing its print() call and the pbapply
# progress bar, and return the coefficient matrix (first list entry).
run_reference_jams <- function(ref, formula, data, cdf) {
  out <- NULL
  capture.output(
    out <- ref$CopRegJAMS(formula = formula, data = data, cdf = cdf, nboots = 1)
  )
  out[[1]]
}

# Seeded data with non-normal endogenous regressors correlated with the error.
make_jams_data <- function(n = 400) {
  set.seed(20260704)
  S <- matrix(c(1.0, 0.4, 0.5,
                0.4, 1.0, 0.5,
                0.5, 0.5, 1.0), nrow = 3)
  Z <- matrix(rnorm(n * 3), nrow = n) %*% chol(S)
  p1 <- exp(Z[, 1])                    # lognormal endogenous regressor
  p2 <- qchisq(pnorm(Z[, 2]), df = 3)  # chi-square endogenous regressor
  eps <- Z[, 3]                        # error correlated with p1 and p2
  x1 <- rnorm(n)
  x2 <- runif(n)
  w <- rep_len(c(0, 1, 2), n)          # three-level discrete exogenous variable
  y <- 1 + 0.8 * p1 - 0.6 * p2 + 0.5 * x1 + 0.3 * x2 +
    0.4 * (w == 1) - 0.2 * (w == 2) + eps
  data.frame(y = y, p1 = p1, p2 = p2, x1 = x1, x2 = x2, w = w)
}

test_that("CopRegJAMS matches reference: 2 endogenous + 2 exogenous continuous", {
  skip_if_no_reference()
  ref <- load_jams_reference()
  dat <- make_jams_data()

  for (cdf in c("resc.ecdf", "adj.ecdf", "ecdf", "kde")) {
    ref_tab <- run_reference_jams(ref, y ~ p1 + p2 | x1 + x2, dat, cdf)
    fit <- CopRegJAMS(y ~ p1 + p2 | x1 + x2, data = dat, cdf = cdf, nboots = 0)

    pkg_est <- fit$coefficients[, "Estimate"]
    ref_est <- ref_tab[, "Estimate"]

    expect_equal(pkg_est, ref_est, tolerance = 1e-8,
                 info = paste("continuous case, cdf =", cdf))
  }
})

test_that("CopRegJAMS matches reference residuals (continuous case)", {
  skip_if_no_reference()
  ref <- load_jams_reference()
  dat <- make_jams_data()

  out <- NULL
  capture.output(
    out <- ref$CopRegJAMS(formula = y ~ p1 + p2 | x1 + x2, data = dat,
                          cdf = "adj.ecdf", nboots = 1)
  )
  fit <- CopRegJAMS(y ~ p1 + p2 | x1 + x2, data = dat, cdf = "adj.ecdf", nboots = 0)
  expect_equal(as.numeric(residuals(fit)), as.numeric(out[[2]]), tolerance = 1e-8)
})

test_that("CopRegJAMS matches reference: endogenous regressors only", {
  skip_if_no_reference()
  ref <- load_jams_reference()
  dat <- make_jams_data()

  for (cdf in c("resc.ecdf", "adj.ecdf")) {
    ref_tab <- run_reference_jams(ref, y ~ p1 + p2, dat, cdf)
    fit <- CopRegJAMS(y ~ p1 + p2, data = dat, cdf = cdf, nboots = 0)

    pkg_est <- fit$coefficients[, "Estimate"]
    ref_est <- ref_tab[, "Estimate"]

    expect_equal(pkg_est, ref_est, tolerance = 1e-8,
                 info = paste("endogenous-only case, cdf =", cdf))
  }
})

test_that("CopRegJAMS matches reference: factor exogenous variable", {
  skip_if_no_reference()
  ref <- load_jams_reference()
  dat <- make_jams_data()

  for (cdf in c("resc.ecdf", "adj.ecdf", "ecdf")) {
    ref_tab <- run_reference_jams(ref, y ~ p1 + p2 | x1 + as.factor(w), dat, cdf)
    fit <- CopRegJAMS(y ~ p1 + p2 | x1 + as.factor(w), data = dat,
                      cdf = cdf, nboots = 0)

    pkg_est <- fit$coefficients[, "Estimate"]
    ref_est <- ref_tab[, "Estimate"]

    expect_equal(pkg_est, ref_est, tolerance = 1e-8,
                 info = paste("factor case, cdf =", cdf))
  }
})
