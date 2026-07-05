# Unit tests for CopRegBMW that do not require the reference implementation.

make_bmw_test_data <- function(n = 80) {
  set.seed(123)
  err <- stats::rnorm(n)
  data.frame(
    y = stats::rnorm(n) + err,
    p1 = exp(stats::rnorm(n) + 0.4 * err),
    p2 = stats::rchisq(n, df = 3) + 0.5 * err,
    z1 = stats::rnorm(n),
    z2 = stats::runif(n),
    f = factor(sample(c("a", "b"), n, replace = TRUE))
  )
}

test_that("CopRegBMW returns an endog_copula_fit with expected structure", {
  dat <- make_bmw_test_data()
  fit <- CopRegBMW(y ~ p1 + p2 | z1 + z2, data = dat, cdf = "resc.ecdf", nboots = 0)

  expect_s3_class(fit, "endog_copula_fit")
  expect_true(is.matrix(fit$coefficients))
  expect_identical(colnames(fit$coefficients), "Estimate")
  expect_true(all(c("p1", "p2", "z1", "z2", "p1_cop", "p2_cop") %in%
                    rownames(fit$coefficients)))
  expect_type(fit$residuals, "double")
  expect_length(fit$residuals, nrow(dat))
  expect_true(all(is.finite(fit$residuals)))
  expect_identical(fit$cdf, "resc.ecdf")
  expect_null(fit$bootstrap)
})

test_that("CopRegBMW supports the endogenous-only branch", {
  dat <- make_bmw_test_data()
  fit <- CopRegBMW(y ~ p1 + p2, data = dat, cdf = "adj.ecdf", nboots = 0)

  expect_s3_class(fit, "endog_copula_fit")
  expect_true(all(c("p1", "p2", "p1_cop", "p2_cop") %in% rownames(fit$coefficients)))
  expect_true(all(is.finite(fit$coefficients[, "Estimate"])))
})

test_that("CopRegBMW handles factor exogenous regressors", {
  dat <- make_bmw_test_data()
  fit <- CopRegBMW(y ~ p1 + p2 | z1 + f, data = dat, cdf = "resc.ecdf", nboots = 0)

  expect_true("fb" %in% rownames(fit$coefficients))
  expect_true(all(is.finite(fit$coefficients[, "Estimate"])))
})

test_that("CopRegBMW gives informative errors", {
  dat <- make_bmw_test_data()

  expect_error(
    CopRegBMW("not a formula", data = dat, cdf = "resc.ecdf", nboots = 0),
    "must be a formula"
  )
  expect_error(
    CopRegBMW(y ~ p1 | z1, data = "not a data frame", cdf = "resc.ecdf", nboots = 0),
    "must be a data.frame"
  )
  expect_error(
    CopRegBMW(y ~ p1 + missing_var | z1, data = dat, cdf = "resc.ecdf", nboots = 0),
    "missing in the data"
  )
  expect_error(
    CopRegBMW(y ~ f | z1, data = dat, cdf = "resc.ecdf", nboots = 0),
    "must be numeric"
  )
  expect_error(
    CopRegBMW(y ~ p1 | z1, data = transform(dat, p1 = 1), cdf = "resc.ecdf", nboots = 0),
    "constant"
  )
  expect_error(
    CopRegBMW(y ~ p1 | z1, data = dat, cdf = "banana", nboots = 0),
    "should be one of"
  )
  expect_error(
    CopRegBMW(y ~ p1 | z1, data = dat, cdf = "resc.ecdf", nboots = -1),
    "non-negative"
  )
})

test_that("CopRegBMW bootstrap adds finite standard errors", {
  dat <- make_bmw_test_data()
  set.seed(99)
  fit <- suppressMessages(
    CopRegBMW(y ~ p1 + p2 | z1 + z2, data = dat, cdf = "resc.ecdf", nboots = 5)
  )

  expect_identical(colnames(fit$coefficients), c("Estimate", "Std. Error"))
  expect_true(all(is.finite(fit$coefficients[, "Std. Error"])))
  expect_true(all(fit$coefficients[, "Std. Error"] > 0))
  expect_true(is.matrix(fit$bootstrap))
  expect_identical(dim(fit$bootstrap), c(5L, nrow(fit$coefficients)))
  expect_identical(colnames(fit$bootstrap), rownames(fit$coefficients))
})
