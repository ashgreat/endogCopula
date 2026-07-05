# Structural and error-handling tests for CopReg2sCOPE. These do not require
# the reference implementation.

simulate_2scope_unit_data <- function(n = 150, seed = 99) {
  set.seed(seed)
  common <- stats::rnorm(n)
  P1 <- exp(0.4 * common + stats::rnorm(n, sd = 0.7))
  P2 <- stats::rchisq(n, df = 4) + 0.3 * common
  X1 <- stats::rnorm(n)
  Fct <- factor(sample(c("lo", "hi"), n, replace = TRUE))
  y <- 1 + 0.7 * P1 - 0.4 * P2 + 0.5 * X1 + 0.8 * common + stats::rnorm(n, sd = 0.5)
  data.frame(y = y, P1 = P1, P2 = P2, X1 = X1, Fct = Fct)
}

test_that("CopReg2sCOPE returns a well-formed endog_copula_fit", {
  dat <- simulate_2scope_unit_data()
  fit <- CopReg2sCOPE(y ~ P1 + P2 | X1, data = dat, cdf = "resc.ecdf", nboots = 0)
  expect_s3_class(fit, "endog_copula_fit")
  expect_true(is.matrix(fit$coefficients))
  expect_identical(colnames(fit$coefficients), "Estimate")
  expect_true(all(c("P1", "P2", "X1", "P1_cop", "P2_cop") %in%
                    rownames(fit$coefficients)))
  expect_true(any(grepl("Intercept", rownames(fit$coefficients))))
  expect_length(fit$residuals, nrow(dat))
  expect_true(all(is.finite(fit$residuals)))
  expect_null(fit$bootstrap)
  expect_identical(fit$method, "Yang et al. (2025)")
  expect_identical(fit$cdf, "resc.ecdf")
})

test_that("CopReg2sCOPE handles the endogenous-only specification", {
  dat <- simulate_2scope_unit_data()
  fit <- CopReg2sCOPE(y ~ P1 + P2, data = dat, cdf = "adj.ecdf", nboots = 0)
  expect_s3_class(fit, "endog_copula_fit")
  expect_true(all(c("P1", "P2", "P1_cop", "P2_cop") %in%
                    rownames(fit$coefficients)))
  expect_length(fit$residuals, nrow(dat))
})

test_that("CopReg2sCOPE rejects invalid inputs informatively", {
  dat <- simulate_2scope_unit_data()
  expect_error(
    CopReg2sCOPE("not a formula", data = dat, cdf = "ecdf", nboots = 0),
    "formula"
  )
  expect_error(
    CopReg2sCOPE(y ~ P1 + P2 | X1, data = "not a data frame", cdf = "ecdf", nboots = 0),
    "data.frame"
  )
  expect_error(
    CopReg2sCOPE(y ~ nonexistent | X1, data = dat, cdf = "ecdf", nboots = 0),
    "missing in the data"
  )
  expect_error(
    CopReg2sCOPE(y ~ Fct | X1, data = dat, cdf = "ecdf", nboots = 0),
    "must be numeric"
  )
  dat$constant <- 1
  expect_error(
    CopReg2sCOPE(y ~ constant | X1, data = dat, cdf = "ecdf", nboots = 0),
    "constant"
  )
  expect_error(
    CopReg2sCOPE(y ~ P1 + P2 | X1, data = dat, cdf = "not-a-cdf", nboots = 0),
    "should be one of"
  )
  expect_error(
    CopReg2sCOPE(y ~ P1 + P2 | X1, data = dat, cdf = "ecdf", nboots = -1),
    "non-negative"
  )
})

test_that("CopReg2sCOPE bootstrap produces finite standard errors", {
  dat <- simulate_2scope_unit_data()
  expect_message(
    fit <- CopReg2sCOPE(y ~ P1 + P2 | X1, data = dat, cdf = "resc.ecdf", nboots = 5),
    "bootstrap"
  )
  expect_identical(colnames(fit$coefficients), c("Estimate", "Std. Error"))
  expect_true(all(is.finite(fit$coefficients[, "Std. Error"])))
  expect_true(all(fit$coefficients[, "Std. Error"] > 0))
  expect_true(is.matrix(fit$bootstrap))
  expect_identical(dim(fit$bootstrap), c(5L, nrow(fit$coefficients)))
  expect_identical(colnames(fit$bootstrap), rownames(fit$coefficients))
  expect_true(all(is.finite(fit$bootstrap)))
})
