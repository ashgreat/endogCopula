# Structural and error-handling tests for CopReg2sCOPEnp. These do not need
# the reference implementation. The endogenous-only branch avoids the 'np'
# dependency entirely, so most tests run everywhere.

sim_2scopenp_unit_data <- function(n, seed) {
  set.seed(seed)
  u <- rnorm(n)
  X1 <- rnorm(n)
  P1 <- exp(0.5 * X1 + 0.5 * u + 0.4 * rnorm(n))
  P2 <- (0.4 * X1 - 0.4 * u + rnorm(n))^2
  y <- 1 + 0.9 * P1 - 0.6 * P2 + 0.5 * X1 + u + 0.3 * rnorm(n)
  data.frame(y = y, P1 = P1, P2 = P2, X1 = X1)
}

test_that("CopReg2sCOPEnp returns a well-formed endog_copula_fit", {
  dat <- sim_2scopenp_unit_data(n = 80, seed = 1)
  fit <- CopReg2sCOPEnp(y ~ P1 + P2, data = dat, cdf = "resc.ecdf", nboots = 0)

  expect_s3_class(fit, "endog_copula_fit")
  expect_true(is.matrix(fit$coefficients))
  expect_identical(colnames(fit$coefficients), "Estimate")
  expect_true(all(c("P1", "P2", "P1_cop", "P2_cop") %in% rownames(fit$coefficients)))
  expect_type(fit$residuals, "double")
  expect_length(fit$residuals, nrow(dat))
  expect_identical(fit$method, "Hu et al. (2025)")
  expect_identical(fit$cdf, "resc.ecdf")
  expect_null(fit$bootstrap)
})

test_that("CopReg2sCOPEnp signals informative errors", {
  dat <- sim_2scopenp_unit_data(n = 60, seed = 2)

  expect_error(
    CopReg2sCOPEnp(y ~ P1 + P2, data = dat, cdf = "unknown", nboots = 0),
    "should be one of"
  )
  expect_error(
    CopReg2sCOPEnp(y ~ P1 + P2, data = dat, nboots = -1),
    "nboots"
  )
  expect_error(
    CopReg2sCOPEnp("not a formula", data = dat, nboots = 0),
    "formula"
  )
  expect_error(
    CopReg2sCOPEnp(y ~ P1 + P2, data = "not a data frame", nboots = 0),
    "data.frame"
  )
  expect_error(
    CopReg2sCOPEnp(y ~ P1 + Pmissing, data = dat, nboots = 0),
    "missing in the data"
  )

  dat$Pfac <- factor(rep(c("a", "b"), length.out = nrow(dat)))
  expect_error(
    CopReg2sCOPEnp(y ~ P1 + Pfac, data = dat, nboots = 0),
    "numeric"
  )

  dat$Pconst <- 1
  expect_error(
    CopReg2sCOPEnp(y ~ P1 + Pconst, data = dat, nboots = 0),
    "constant"
  )
})

test_that("CopReg2sCOPEnp bootstrap adds finite standard errors", {
  dat <- sim_2scopenp_unit_data(n = 60, seed = 3)
  set.seed(42)
  fit <- suppressMessages(
    CopReg2sCOPEnp(y ~ P1 + P2, data = dat, cdf = "resc.ecdf", nboots = 5)
  )

  expect_true(all(c("Estimate", "Std. Error") %in% colnames(fit$coefficients)))
  expect_true(all(is.finite(fit$coefficients[, "Std. Error"])))
  expect_true(is.matrix(fit$bootstrap))
  expect_identical(dim(fit$bootstrap), c(5L, nrow(fit$coefficients)))
  expect_identical(colnames(fit$bootstrap), rownames(fit$coefficients))
})

test_that("CopReg2sCOPEnp np branch returns copula control terms", {
  skip_if_not_installed("np")
  old_opts <- options(np.messages = FALSE)
  on.exit(options(old_opts), add = TRUE)

  dat <- sim_2scopenp_unit_data(n = 60, seed = 4)
  fit <- CopReg2sCOPEnp(y ~ P1 | X1, data = dat, cdf = "resc.ecdf", nboots = 0)

  expect_s3_class(fit, "endog_copula_fit")
  expect_true(all(c("P1", "X1", "P1_cop") %in% rownames(fit$coefficients)))
  expect_true(all(is.finite(fit$coefficients[, "Estimate"])))
})
