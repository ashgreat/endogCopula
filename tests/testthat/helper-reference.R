# Helpers for comparing package output against the published reference
# implementation in Copula-based-endogeneity-corrections-main.
#
# The reference directory is located via the ENDOG_REF_DIR environment
# variable, falling back to a sibling directory of the package. Tests that
# rely on the reference code should call skip_if_no_reference() first.

endog_ref_dir <- function() {
  dir <- Sys.getenv("ENDOG_REF_DIR", unset = "")
  if (nzchar(dir)) {
    return(dir)
  }
  testthat::test_path("..", "..", "..", "Copula-based-endogeneity-corrections-main")
}

skip_if_no_reference <- function() {
  if (!dir.exists(endog_ref_dir())) {
    testthat::skip("Reference implementation directory not available")
  }
  invisible(TRUE)
}

# Parse a reference .R file and evaluate ONLY top-level assignments whose
# right-hand side is a function definition into a fresh environment. This
# keeps the demo code at the bottom of the reference scripts (library(bayesm),
# data(...), hist(...), model runs) from ever executing. The reference
# functions call many symbols unqualified (pobs, kcde, sample_n, %>%,
# rmvnorm, ...), so tests must attach the required packages first via
# suppressMessages(library(...)).
load_reference_functions <- function(filename) {
  path <- file.path(endog_ref_dir(), filename)
  if (!file.exists(path)) {
    testthat::skip(sprintf("Reference file '%s' not available", filename))
  }
  exprs <- parse(file = path, keep.source = FALSE)
  env <- new.env(parent = globalenv())
  for (expr in exprs) {
    if (!is.call(expr) || length(expr) != 3L) {
      next
    }
    op <- expr[[1L]]
    is_assign <- identical(op, as.name("<-")) || identical(op, as.name("="))
    if (!is_assign) {
      next
    }
    rhs <- expr[[3L]]
    if (is.call(rhs) && identical(rhs[[1L]], as.name("function"))) {
      eval(expr, envir = env)
    }
  }
  env
}
