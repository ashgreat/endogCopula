# Unit tests for CopRegIMA that do not require the reference implementation.

make_ima_unit_data <- function(n = 150, seed = 42) {
  set.seed(seed)
  err <- rnorm(n)
  dat <- data.frame(
    P1 = exp(0.5 * rnorm(n) + 0.5 * err),
    P2 = rchisq(n, df = 4) + 0.4 * err,
    X1 = rnorm(n),
    X2 = runif(n)
  )
  dat$y <- 1 + 2 * dat$P1 - dat$P2 + 0.5 * dat$X1 + 0.3 * dat$X2 + err
  dat
}

test_that("CopRegIMA returns a well-formed endog_copula_fit", {
  dat <- make_ima_unit_data()
  fit <- CopRegIMA(y ~ P1 + P2 | X1 + X2, data = dat, cdf = "resc.ecdf", nboots = 0)

  expect_s3_class(fit, "endog_copula_fit")
  expect_true(is.matrix(fit$coefficients))
  expect_identical(colnames(fit$coefficients), "Estimate")
  expect_true(all(c("P1", "P2", "X1", "X2", "P1_cop", "P2_cop") %in%
                    rownames(fit$coefficients)))
  expect_length(residuals(fit), nrow(dat))
  expect_true(all(is.finite(residuals(fit))))
  expect_identical(fit$method, "Haschka (2024)")
  expect_identical(fit$cdf, "resc.ecdf")
  expect_null(fit$bootstrap)
})

test_that("CopRegIMA supports endogenous-only formulas", {
  dat <- make_ima_unit_data()
  fit <- CopRegIMA(y ~ P1 + P2, data = dat, cdf = "adj.ecdf", nboots = 0)

  expect_s3_class(fit, "endog_copula_fit")
  expect_true(all(c("P1", "P2", "P1_cop", "P2_cop") %in% rownames(fit$coefficients)))
  expect_length(residuals(fit), nrow(dat))
  expect_true(all(is.finite(residuals(fit))))
})

test_that("CopRegIMA rejects invalid inputs with informative errors", {
  dat <- make_ima_unit_data()

  expect_error(
    CopRegIMA(y ~ P1 | X1, data = dat, cdf = "spline", nboots = 0),
    "should be one of"
  )
  expect_error(
    CopRegIMA(y ~ P1 | X1, data = dat, cdf = "resc.ecdf", nboots = -1),
    "non-negative"
  )
  expect_error(
    CopRegIMA(y ~ P1 | X1, data = "not a data frame", cdf = "resc.ecdf", nboots = 0),
    "data.frame"
  )
  expect_error(
    CopRegIMA(y ~ P1 + missing_var | X1, data = dat, cdf = "resc.ecdf", nboots = 0),
    "missing in the data"
  )

  dat$fac <- factor(rep(c("a", "b"), length.out = nrow(dat)))
  expect_error(
    CopRegIMA(y ~ fac | X1, data = dat, cdf = "resc.ecdf", nboots = 0),
    "numeric"
  )

  dat$const <- 1
  expect_error(
    CopRegIMA(y ~ const | X1, data = dat, cdf = "resc.ecdf", nboots = 0),
    "constant"
  )
})

test_that("CopRegIMA bootstrap produces finite standard errors", {
  dat <- make_ima_unit_data(n = 120, seed = 7)
  fit <- suppressMessages(
    CopRegIMA(y ~ P1 | X1 + X2, data = dat, cdf = "resc.ecdf", nboots = 5)
  )

  expect_identical(colnames(fit$coefficients), c("Estimate", "Std. Error"))
  expect_true(all(is.finite(fit$coefficients[, "Std. Error"])))
  expect_true(all(fit$coefficients[, "Std. Error"] > 0))
  expect_true(is.matrix(fit$bootstrap))
  expect_identical(dim(fit$bootstrap), c(5L, nrow(fit$coefficients)))
  expect_identical(colnames(fit$bootstrap), rownames(fit$coefficients))
})
