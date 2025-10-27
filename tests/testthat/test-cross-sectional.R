test_that("CopRegPG returns an endog_copula_fit", {
  set.seed(42)
  n <- 40
  dat <- data.frame(
    y = rnorm(n),
    endog = rnorm(n),
    exog1 = rnorm(n),
    exog2 = rnorm(n)
  )
  fit <- CopRegPG(y ~ endog | exog1 + exog2, data = dat, cdf = "resc.ecdf", nboots = 0)
  expect_s3_class(fit, "endog_copula_fit")
  expect_type(residuals(fit), "double")
  expect_true("Estimate" %in% colnames(fit$coefficients))
})
