# Numerical equivalence of CopReg2sCOPEnp against the published reference
# implementation in Copula-based-endogeneity-corrections-main/CopReg2sCOPE-np.R.
#
# Coefficient name mapping: both implementations regress the response on the
# columns of an identical estimation matrix via lm(y ~ . - y - 1), so the
# coefficient names agree exactly. Non-syntactic column names are backticked
# by lm in both implementations, e.g. the intercept coefficient is named
# "`(Intercept)`" and copula terms are named "<var>_cop".
#
# The reference bootstrap runs after point estimation and only feeds the
# "Std. Error" column, so the reference is called with nboots = 1 (the
# smallest value its pbsapply(1:nboots, ...) loop tolerates) and only the
# "Estimate" column is compared.

simulate_2scopenp_data <- function(n, seed, with_exog = TRUE, with_factor = FALSE) {
  set.seed(seed)
  u <- rnorm(n)
  X1 <- rnorm(n)
  X2 <- runif(n, -1, 1)
  F1 <- factor(sample(c("a", "b"), n, replace = TRUE))
  f1_shift <- 0.5 * (F1 == "b")
  # Non-normal endogenous regressors, correlated with the error through u.
  P1 <- exp(0.6 * X1 + 0.3 * f1_shift + 0.5 * u + 0.4 * rnorm(n))
  P2 <- (0.5 * X2 - 0.4 * u + rnorm(n))^2
  y <- 0.5 + 1.2 * P1 - 0.7 * P2 + 0.8 * X1 + 0.6 * X2 + 0.4 * f1_shift +
    u + 0.3 * rnorm(n)
  dat <- data.frame(y = y, P1 = P1, P2 = P2, X1 = X1, X2 = X2)
  if (with_factor) {
    dat$F1 <- F1
  }
  if (!with_exog) {
    dat <- dat[, c("y", "P1", "P2")]
  }
  dat
}

expect_matches_reference <- function(fit_pkg, res_ref, tolerance = 1e-6) {
  pkg_est <- fit_pkg$coefficients[, "Estimate"]
  ref_tab <- res_ref[[1]]
  testthat::expect_setequal(names(pkg_est), rownames(ref_tab))
  ref_est <- ref_tab[names(pkg_est), "Estimate"]
  testthat::expect_equal(unname(pkg_est), unname(ref_est), tolerance = tolerance)
  testthat::expect_equal(
    as.numeric(fit_pkg$residuals),
    as.numeric(res_ref[[2]]),
    tolerance = tolerance
  )
}

test_that("CopReg2sCOPEnp matches reference: 2 endogenous + 2 continuous exogenous", {
  skip_if_no_reference()
  skip_if_not_installed("np")
  suppressWarnings(suppressMessages({
    library(dplyr)
    library(pbapply)
    library(np)
  }))
  old_opts <- options(np.messages = FALSE)
  old_pbo <- pbapply::pboptions(type = "none")
  on.exit({
    options(old_opts)
    pbapply::pboptions(old_pbo)
  }, add = TRUE)

  ref <- load_reference_functions("CopReg2sCOPE-np.R")
  dat <- simulate_2scopenp_data(n = 150, seed = 101)
  f <- y ~ P1 + P2 | X1 + X2

  # np's conditional CDF (not the marginal cdf argument) drives this branch,
  # so a single cdf value suffices; npcdistbw with default settings is
  # deterministic, but both calls are seeded identically for safety.
  set.seed(11)
  fit_pkg <- CopReg2sCOPEnp(f, data = dat, cdf = "resc.ecdf", nboots = 0)
  set.seed(11)
  capture.output(
    res_ref <- ref$CopReg2sCOPEnp(formula = f, data = dat, nboots = 1)
  )

  expect_matches_reference(fit_pkg, res_ref)
})

test_that("CopReg2sCOPEnp matches reference: factor exogenous regressor", {
  skip_if_no_reference()
  skip_if_not_installed("np")
  suppressWarnings(suppressMessages({
    library(dplyr)
    library(pbapply)
    library(np)
  }))
  old_opts <- options(np.messages = FALSE)
  old_pbo <- pbapply::pboptions(type = "none")
  on.exit({
    options(old_opts)
    pbapply::pboptions(old_pbo)
  }, add = TRUE)

  ref <- load_reference_functions("CopReg2sCOPE-np.R")
  dat <- simulate_2scopenp_data(n = 120, seed = 202, with_factor = TRUE)
  dat <- dat[, c("y", "P1", "X1", "F1")]
  f <- y ~ P1 | X1 + F1

  set.seed(22)
  fit_pkg <- CopReg2sCOPEnp(f, data = dat, cdf = "resc.ecdf", nboots = 0)
  set.seed(22)
  capture.output(
    res_ref <- ref$CopReg2sCOPEnp(formula = f, data = dat, nboots = 1)
  )

  # Factor levels expand to model-matrix dummies (e.g. "F1b") in both
  # implementations, so names again agree exactly.
  expect_matches_reference(fit_pkg, res_ref)
})

test_that("CopReg2sCOPEnp matches reference: endogenous-only branch across cdf methods", {
  skip_if_no_reference()
  suppressMessages({
    library(dplyr)
    library(pbapply)
  })
  old_pbo <- pbapply::pboptions(type = "none")
  on.exit(pbapply::pboptions(old_pbo), add = TRUE)

  ref <- load_reference_functions("CopReg2sCOPE-np.R")
  # pobs1 (used for cdf = "adj.ecdf") is defined in the other reference
  # scripts but not in CopReg2sCOPE-np.R; borrow it from CopRegPG.R.
  ref_pg <- load_reference_functions("CopRegPG.R")
  assign("pobs1", get("pobs1", envir = ref_pg), envir = ref)

  dat <- simulate_2scopenp_data(n = 400, seed = 303, with_exog = FALSE)
  f <- y ~ P1 + P2

  for (method in c("resc.ecdf", "adj.ecdf", "ecdf", "kde")) {
    # In the reference the endogenous-only branch reads `cdf` as a free
    # variable (CopReg2sCOPEnp has no cdf argument), so it must be placed in
    # the environment enclosing the reference functions.
    assign("cdf", method, envir = ref)

    set.seed(33)
    fit_pkg <- CopReg2sCOPEnp(f, data = dat, cdf = method, nboots = 0)
    set.seed(33)
    capture.output(
      res_ref <- ref$CopReg2sCOPEnp(formula = f, data = dat, nboots = 1)
    )

    expect_matches_reference(fit_pkg, res_ref)
  }
})
