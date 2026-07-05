test_that("cdf_transform ecdf matches the reference boundary replacement", {
  set.seed(101)
  x <- cbind(a = rnorm(200), b = round(rexp(200), 1))

  result <- cdf_transform(x, "ecdf")

  expect_equal(dim(result), dim(x))
  expect_equal(colnames(result), colnames(x))
  expect_true(all(result > 0 & result < 1))

  for (j in seq_len(ncol(x))) {
    Fhat <- stats::ecdf(x[, j])
    expected <- Fhat(x[, j])
    expected[expected == min(expected)] <- 10e-7
    expected[expected == max(expected)] <- 1 - 10e-7
    expect_equal(unname(result[, j]), expected, tolerance = 0)
  }
})

test_that("cdf_transform resc.ecdf matches copula::pobs", {
  skip_if_not_installed("copula")
  set.seed(102)
  x <- cbind(a = rnorm(150), b = sample(c(1, 2, 2, 5), 150, replace = TRUE))

  result <- cdf_transform(x, "resc.ecdf")

  expect_true(all(result > 0 & result < 1))
  expected <- apply(x, 2, copula::pobs)
  expect_equal(unname(result), unname(expected), tolerance = 0)
})

test_that("cdf_transform adj.ecdf matches the reference pobs1 formula", {
  set.seed(103)
  n <- 150
  x <- cbind(a = rnorm(n), b = round(rnorm(n), 1))

  result <- cdf_transform(x, "adj.ecdf")

  expect_true(all(result > 0 & result < 1))
  for (j in seq_len(ncol(x))) {
    ranks <- rank(x[, j], na.last = "keep", ties.method = "average")
    expected <- ranks * ((n - 1) / (n^2)) + 1 / (2 * n)
    expect_equal(unname(result[, j]), expected, tolerance = 0)
  }
})

test_that("cdf_transform kde matches ks::kcde predictions", {
  skip_if_not_installed("ks")
  set.seed(104)
  x <- cbind(a = rnorm(120))

  result <- cdf_transform(x, "kde")

  expect_true(all(result > 0 & result < 1))
  fit <- ks::kcde(x[, 1])
  expected <- as.numeric(predict(fit, x = x[, 1]))
  expect_equal(unname(result[, 1]), expected, tolerance = 0)
})

test_that("cdf_transform validates its input", {
  expect_error(cdf_transform(1:10, "ecdf"), "must be a matrix")
  x_bad <- matrix(letters[1:10], ncol = 1)
  expect_error(cdf_transform(x_bad, "ecdf"), "numeric")
})

test_that("bootstrap_estimates returns a named matrix of the right shape", {
  set.seed(105)
  dat <- data.frame(x = rnorm(50), y = rnorm(50))
  fit_fun <- function(d) c(mean_x = mean(d$x), mean_y = mean(d$y))

  result <- bootstrap_estimates(dat, fit_fun, nboots = 7,
                                expected_names = c("mean_x", "mean_y"))

  expect_true(is.matrix(result))
  expect_equal(dim(result), c(7L, 2L))
  expect_equal(colnames(result), c("mean_x", "mean_y"))
  expect_false(anyNA(result))
})

test_that("bootstrap_estimates retries when fit_fun errors, returns NA, or wrong names", {
  set.seed(106)
  dat <- data.frame(x = rnorm(30))
  calls <- 0L
  fit_fun <- function(d) {
    calls <<- calls + 1L
    if (calls == 1L) stop("marker failure")
    if (calls == 2L) return(c(estimate = NA_real_))
    if (calls == 3L) return(c(wrong_name = mean(d$x)))
    c(estimate = mean(d$x))
  }

  result <- bootstrap_estimates(dat, fit_fun, nboots = 2,
                                expected_names = "estimate")

  expect_equal(dim(result), c(2L, 1L))
  expect_equal(colnames(result), "estimate")
  expect_false(anyNA(result))
  expect_equal(calls, 5L)
})

test_that("bootstrap_estimates stops after max_tries failures", {
  set.seed(107)
  dat <- data.frame(x = rnorm(30))
  fit_fun <- function(d) stop("always fails")

  expect_error(
    bootstrap_estimates(dat, fit_fun, nboots = 1,
                        expected_names = "estimate", max_tries = 5L),
    "failed after 5 attempts"
  )
})

test_that("reference helper loads only function definitions", {
  skip_if_no_reference()

  ref_env <- load_reference_functions("CopRegPG.R")

  expect_true(is.function(ref_env$pobs1))
  expect_true(is.function(ref_env$CopRegPG))
  expect_true(is.function(ref_env$boots1))
  expect_false(exists("dat1", envir = ref_env, inherits = FALSE))
  expect_false(exists("modPG1", envir = ref_env, inherits = FALSE))

  set.seed(108)
  x <- round(rnorm(75), 1)
  expect_equal(adjusted_pobs(x), ref_env$pobs1(x), tolerance = 0)
})
