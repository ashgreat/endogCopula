test_that("simulated dataset works with copula estimators", {
  data(sim_endog, package = "endogCopula")

  ols_fit <- lm(y ~ z_endog + x_exog, data = sim_endog)
  ols_beta <- coef(ols_fit)["z_endog"]

  pg_fit <- CopRegPG(
    formula = y ~ z_endog + x_exog | w_instr,
    data = sim_endog,
    cdf = "resc.ecdf",
    nboots = 0
  )
  scope_fit <- CopReg2sCOPE(
    formula = y ~ z_endog + x_exog | w_instr,
    data = sim_endog,
    cdf = "resc.ecdf",
    nboots = 0
  )

  beta_pg <- pg_fit$coefficients["z_endog", "Estimate"]
  beta_scope <- scope_fit$coefficients["z_endog", "Estimate"]

  expect_true(is.finite(beta_pg))
  expect_true(is.finite(beta_scope))
  expect_lt(abs(beta_pg - beta_scope), 0.5)
  expect_gt(abs(beta_pg - ols_beta), 0.01)
})
