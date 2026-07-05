# Numerical equivalence of CopRegPG against the published reference
# implementation in Copula-based-endogeneity-corrections-main/CopRegPG.R.
#
# Coefficient name mapping: both implementations fit the final
# control-function regression through lm(y ~ . - y - 1, data =
# as.data.frame(est_matrix)), so coefficient names agree: non-syntactic
# columns are backtick-quoted identically by lm() (e.g. "(Intercept)" becomes
# "`(Intercept)`"), and copula terms are "<name>_cop". Coefficients can
# therefore be matched by name. The one exception is the single-endogenous /
# no-exogenous case: the reference extracts the endogenous column via
# design_matrix[, -1] without drop = FALSE, losing the column name, so its
# copula term is named "`_cop`" (with literal backticks) while the package
# keeps "P1_cop". The ref_rename argument maps between the two.
#
# Point-estimate/residual tests ignore the bootstrap: the reference runs it
# unconditionally, so it is called with nboots = 1 (the smallest value its
# sapply/apply pipeline tolerates) and the resulting standard errors are
# discarded. Bootstrap standard errors are NOT compared against the
# reference: boots1/boots_PG in CopRegPG.R omit the qnorm() step that the
# reference main fit applies (the other reference estimators apply qnorm()
# in both places), so the reference bootstrap resamples a different
# estimator than the one it reports. The package deliberately refits the
# same estimator in every replicate; the tests below pin down that contract
# and seeded reproducibility instead.

suppressMessages({
  library(dplyr)
  library(copula)
  library(ks)
  library(Matrix)
  library(nlme)
})

simulate_pg_data <- function(n = 400, seed = 20260704) {
  set.seed(seed)
  z1 <- rnorm(n)
  z2 <- 0.5 * z1 + sqrt(1 - 0.5^2) * rnorm(n)
  eps <- 0.8 * z1 + 0.5 * z2 + rnorm(n)
  P1 <- exp(z1)
  P2 <- qchisq(pnorm(z2), df = 3)
  X1 <- rnorm(n)
  X2 <- runif(n)
  F1 <- factor(sample(c("a", "b", "c"), n, replace = TRUE))
  y <- 1 + 0.8 * P1 - 0.5 * P2 + 0.7 * X1 + 0.3 * X2 +
    0.4 * (F1 == "b") - 0.2 * (F1 == "c") + eps
  data.frame(y = y, P1 = P1, P2 = P2, X1 = X1, X2 = X2, F1 = F1)
}

pg_reference_env <- function() {
  env <- load_reference_functions("CopRegPG.R")
  # The reference bootstrap loop calls pbapply::pbsapply() unqualified; a
  # plain sapply() shim keeps the tests independent of pbapply and silences
  # its progress bar. Bootstrap output is discarded anyway.
  env$pbsapply <- function(X, FUN, ...) sapply(X, FUN, ...)
  env
}

expect_pg_equivalent <- function(ref_env, formula, data, cdf, ref_rename = NULL) {
  set.seed(99)
  pkg <- CopRegPG(formula, data = data, cdf = cdf, nboots = 0)
  set.seed(99)
  ref <- NULL
  # capture.output() swallows the reference's print() progress message.
  invisible(capture.output(
    ref <- ref_env$CopRegPG(formula = formula, data = data, cdf = cdf, nboots = 1)
  ))
  ref_est <- ref[[1]][, "Estimate"]
  if (!is.null(ref_rename)) {
    hits <- names(ref_est) %in% names(ref_rename)
    names(ref_est)[hits] <- ref_rename[names(ref_est)[hits]]
  }
  pkg_est <- pkg$coefficients[, "Estimate"]
  expect_setequal(names(pkg_est), names(ref_est))
  expect_equal(pkg_est[names(ref_est)], ref_est, tolerance = 1e-8)
  expect_equal(as.numeric(residuals(pkg)), as.numeric(ref[[2]]), tolerance = 1e-8)
  invisible(NULL)
}

pg_cdf_methods <- c("resc.ecdf", "adj.ecdf", "ecdf", "kde")

test_that("CopRegPG matches reference: 2 endogenous + 2 exogenous", {
  skip_if_no_reference()
  ref_env <- pg_reference_env()
  dat <- simulate_pg_data()
  for (cdf in pg_cdf_methods) {
    expect_pg_equivalent(ref_env, y ~ P1 + P2 | X1 + X2, dat, cdf)
  }
})

test_that("CopRegPG matches reference: 2 endogenous, no exogenous", {
  skip_if_no_reference()
  ref_env <- pg_reference_env()
  dat <- simulate_pg_data()
  for (cdf in pg_cdf_methods) {
    expect_pg_equivalent(ref_env, y ~ P1 + P2, dat, cdf)
  }
})

test_that("CopRegPG matches reference: 1 endogenous + 2 exogenous", {
  skip_if_no_reference()
  ref_env <- pg_reference_env()
  dat <- simulate_pg_data()
  for (cdf in pg_cdf_methods) {
    expect_pg_equivalent(ref_env, y ~ P1 | X1 + X2, dat, cdf)
  }
})

test_that("CopRegPG matches reference: 1 endogenous, no exogenous", {
  skip_if_no_reference()
  ref_env <- pg_reference_env()
  dat <- simulate_pg_data()
  for (cdf in pg_cdf_methods) {
    # Reference loses the column name here (see header comment): its copula
    # coefficient is "`_cop`" where the package reports "P1_cop".
    expect_pg_equivalent(ref_env, y ~ P1, dat, cdf,
                         ref_rename = c("`_cop`" = "P1_cop"))
  }
})

test_that("CopRegPG matches reference: factor exogenous regressor", {
  skip_if_no_reference()
  ref_env <- pg_reference_env()
  dat <- simulate_pg_data()
  for (cdf in pg_cdf_methods) {
    expect_pg_equivalent(ref_env, y ~ P1 + P2 | X1 + F1, dat, cdf)
  }
})

test_that("CopRegPG bootstrap is seeded-reproducible and refits the main-fit estimator", {
  dat <- simulate_pg_data(n = 150)
  formula <- y ~ P1 + P2 | X1 + X2

  set.seed(2468)
  fit1 <- suppressMessages(
    CopRegPG(formula, data = dat, cdf = "adj.ecdf", nboots = 20)
  )
  set.seed(2468)
  fit2 <- suppressMessages(
    CopRegPG(formula, data = dat, cdf = "adj.ecdf", nboots = 20)
  )
  expect_identical(fit1$coefficients, fit2$coefficients)
  expect_identical(fit1$bootstrap, fit2$bootstrap)

  # Documented divergence from the reference: each replicate must refit the
  # SAME estimator as the main fit (qnorm-transformed copula terms). Replay
  # the seeded resampling sequence manually through pg_fit() and demand
  # replicate-for-replicate agreement with the stored bootstrap matrix.
  components <- extract_components(formula, dat)
  cleaned <- select_complete_cases(dat, unique(c(
    components$response, components$endogenous, components$exogenous
  )))
  set.seed(2468)
  manual <- bootstrap_estimates(
    cleaned,
    function(boot_data) {
      pg_fit(boot_data, components, formula, "adj.ecdf")$coefficients[, 1]
    },
    nboots = 20,
    expected_names = rownames(fit1$coefficients)
  )
  expect_equal(fit1$bootstrap, manual, tolerance = 1e-12)
})
