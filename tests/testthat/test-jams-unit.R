# Unit tests for CopRegJAMS that do not require the reference implementation.

make_jams_unit_data <- function(n = 120) {
  set.seed(99)
  common <- rnorm(n)
  p1 <- exp(0.6 * common + rnorm(n))
  p2 <- qchisq(pnorm(0.5 * common + 0.8 * rnorm(n)), df = 4)
  x1 <- rnorm(n)
  x2 <- runif(n)
  w <- rep_len(c(0, 1, 2), n)
  y <- 1 + 0.8 * p1 - 0.6 * p2 + 0.5 * x1 + 0.3 * x2 + common + rnorm(n)
  data.frame(y = y, p1 = p1, p2 = p2, x1 = x1, x2 = x2, w = w)
}

test_that("CopRegJAMS returns an endog_copula_fit with the expected structure", {
  dat <- make_jams_unit_data()
  fit <- CopRegJAMS(y ~ p1 + p2 | x1 + x2, data = dat, cdf = "resc.ecdf",
                    nboots = 0)

  expect_s3_class(fit, "endog_copula_fit")
  expect_true(is.matrix(fit$coefficients))
  expect_identical(colnames(fit$coefficients), "Estimate")
  expect_true(all(c("p1", "p2", "x1", "x2", "p1_cop", "p2_cop") %in%
                    rownames(fit$coefficients)))
  expect_length(as.numeric(residuals(fit)), nrow(dat))
  expect_true(all(is.finite(fit$coefficients[, "Estimate"])))
  expect_identical(fit$method, "Liengaard et al. (2025)")
  expect_identical(fit$cdf, "resc.ecdf")
  expect_null(fit$bootstrap)
})

test_that("CopRegJAMS handles the endogenous-only path", {
  dat <- make_jams_unit_data()
  fit <- CopRegJAMS(y ~ p1 + p2, data = dat, cdf = "adj.ecdf", nboots = 0)

  expect_s3_class(fit, "endog_copula_fit")
  expect_true(all(c("p1", "p2", "p1_cop", "p2_cop") %in%
                    rownames(fit$coefficients)))
})

test_that("CopRegJAMS handles a factor exogenous variable", {
  dat <- make_jams_unit_data()
  fit <- CopRegJAMS(y ~ p1 + p2 | x1 + as.factor(w), data = dat,
                    cdf = "resc.ecdf", nboots = 0)

  expect_s3_class(fit, "endog_copula_fit")
  coef_names <- rownames(fit$coefficients)
  expect_true("(Intercept)" %in% coef_names)
  # stratified correction terms for each factor level
  expect_true(all(c("p1_w_0_cop", "p1_w_1_cop", "p1_w_2_cop") %in% coef_names))
  expect_true(all(is.finite(fit$coefficients[, "Estimate"])))
})

test_that("CopRegJAMS raises informative errors", {
  dat <- make_jams_unit_data()

  expect_error(CopRegJAMS("not a formula", data = dat, cdf = "ecdf", nboots = 0),
               "not a formula object")
  expect_error(CopRegJAMS(y ~ p1 | x1, data = "not a data frame", cdf = "ecdf",
                          nboots = 0),
               "not a data.frame object")
  expect_error(CopRegJAMS(y ~ p1 | x1, data = dat, cdf = "banana", nboots = 0),
               "should be one of")
  expect_error(CopRegJAMS(y ~ p1 | x1, data = dat, cdf = "ecdf", nboots = "a"),
               "nboots is not numeric")
  expect_error(CopRegJAMS(y ~ missing_var | x1, data = dat, cdf = "ecdf",
                          nboots = 0),
               "missing in the data")
  dat$fac <- factor(rep_len(c("a", "b"), nrow(dat)))
  expect_error(CopRegJAMS(y ~ fac | x1, data = dat, cdf = "ecdf", nboots = 0),
               "can be endogenous")
})

test_that("CopRegJAMS bootstrap adds finite standard errors", {
  dat <- make_jams_unit_data()
  set.seed(7)
  fit <- suppressMessages(
    CopRegJAMS(y ~ p1 + p2 | x1 + x2, data = dat, cdf = "resc.ecdf", nboots = 5)
  )

  expect_identical(colnames(fit$coefficients), c("Estimate", "Std. Error"))
  expect_true(all(is.finite(fit$coefficients[, "Std. Error"])))
  expect_true(is.matrix(fit$bootstrap))
  expect_identical(nrow(fit$bootstrap), 5L)
  expect_identical(colnames(fit$bootstrap), rownames(fit$coefficients))
})

test_that("CopRegJAMS bootstrap works on the factor path", {
  dat <- make_jams_unit_data()
  set.seed(11)
  fit <- suppressMessages(
    CopRegJAMS(y ~ p1 + p2 | x1 + as.factor(w), data = dat, cdf = "adj.ecdf",
               nboots = 5)
  )

  expect_identical(colnames(fit$coefficients), c("Estimate", "Std. Error"))
  expect_true(all(is.finite(fit$coefficients[, "Std. Error"])))
  expect_identical(nrow(fit$bootstrap), 5L)
})
