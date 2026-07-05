# Structural and error-handling tests for CopRegPG. These do not require the
# reference implementation.

make_pg_unit_data <- function(n = 120, seed = 7) {
  set.seed(seed)
  z <- rnorm(n)
  P1 <- exp(z)
  P2 <- qchisq(pnorm(0.4 * z + sqrt(1 - 0.4^2) * rnorm(n)), df = 4)
  X1 <- rnorm(n)
  eps <- 0.7 * z + rnorm(n)
  y <- 1 + 0.6 * P1 - 0.4 * P2 + 0.5 * X1 + eps
  data.frame(y = y, P1 = P1, P2 = P2, X1 = X1)
}

test_that("CopRegPG returns a well-formed endog_copula_fit", {
  dat <- make_pg_unit_data()
  fit <- CopRegPG(y ~ P1 + P2 | X1, data = dat, cdf = "resc.ecdf", nboots = 0)
  expect_s3_class(fit, "endog_copula_fit")
  expect_identical(fit$method, "Park and Gupta (2012)")
  expect_identical(fit$cdf, "resc.ecdf")
  expect_true(is.matrix(fit$coefficients))
  expect_identical(colnames(fit$coefficients), "Estimate")
  expect_true(all(c("`(Intercept)`", "P1", "P2", "X1", "P1_cop", "P2_cop") %in%
                    rownames(fit$coefficients)))
  expect_length(residuals(fit), nrow(dat))
  expect_true(all(is.finite(residuals(fit))))
  expect_null(fit$bootstrap)
})

test_that("CopRegPG supports endogenous-only formulas", {
  dat <- make_pg_unit_data()
  fit <- CopRegPG(y ~ P1 + P2, data = dat, cdf = "adj.ecdf", nboots = 0)
  expect_s3_class(fit, "endog_copula_fit")
  expect_true(all(c("`(Intercept)`", "P1", "P2", "P1_cop", "P2_cop") %in%
                    rownames(fit$coefficients)))
  expect_length(residuals(fit), nrow(dat))
})

test_that("CopRegPG rejects invalid inputs with informative errors", {
  dat <- make_pg_unit_data()
  expect_error(CopRegPG("y ~ P1", data = dat, cdf = "ecdf"), "formula")
  expect_error(CopRegPG(y ~ P1 | X1, data = as.matrix(dat), cdf = "ecdf"),
               "data.frame")
  expect_error(CopRegPG(y ~ P1 | X1, data = dat, cdf = "banana"))
  expect_error(CopRegPG(y ~ P1 | X1, data = dat, cdf = "ecdf", nboots = -1),
               "non-negative")
  expect_error(CopRegPG(y ~ P1 | X1, data = dat, cdf = "ecdf", nboots = c(1, 2)),
               "non-negative")
  expect_error(CopRegPG(y ~ nope | X1, data = dat, cdf = "ecdf"),
               "missing in the data")
  dt_like <- structure(dat, class = c("data.table", "data.frame"))
  expect_error(CopRegPG(y ~ P1 | X1, data = dt_like, cdf = "ecdf"),
               "data.table")
  dat$const <- 1
  expect_error(CopRegPG(y ~ const | X1, data = dat, cdf = "ecdf", nboots = 0),
               "constant")
  dat$fac <- factor(rep(c("a", "b"), length.out = nrow(dat)))
  expect_error(CopRegPG(y ~ fac | X1, data = dat, cdf = "ecdf", nboots = 0),
               "numeric")
})

test_that("CopRegPG bootstrap produces finite standard errors", {
  dat <- make_pg_unit_data()
  set.seed(11)
  expect_message(
    fit <- CopRegPG(y ~ P1 + P2 | X1, data = dat, cdf = "resc.ecdf", nboots = 5),
    "bootstrap"
  )
  expect_identical(colnames(fit$coefficients), c("Estimate", "Std. Error"))
  expect_true(all(is.finite(fit$coefficients[, "Std. Error"])))
  expect_true(all(fit$coefficients[, "Std. Error"] > 0))
  expect_true(is.matrix(fit$bootstrap))
  expect_identical(dim(fit$bootstrap), c(5L, nrow(fit$coefficients)))
  expect_identical(colnames(fit$bootstrap), rownames(fit$coefficients))
})
