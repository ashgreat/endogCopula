NULL

adjusted_pobs <- function(x, lower.tail = TRUE) {
  if (!is.null(dim(x))) {
    n <- nrow(x)
    ranks <- apply(x, 2, base::rank, na.last = "keep", ties.method = "average")
    U <- ranks * ((n - 1) / (n^2)) + 1 / (2 * n)
  } else {
    n <- length(x)
    ranks <- base::rank(x, na.last = "keep", ties.method = "average")
    U <- ranks * ((n - 1) / (n^2)) + 1 / (2 * n)
  }
  if (lower.tail) {
    U
  } else {
    1 - U
  }
}

pseudo_observations <- function(x) {
  transform_vec <- function(v) {
    n <- length(v)
    ranks <- base::rank(v, na.last = "keep", ties.method = "average")
    (ranks - 0.5) / n
  }
  if (is.null(dim(x))) {
    res <- transform_vec(x)
  } else {
    res <- vapply(seq_len(ncol(x)), function(j) transform_vec(x[, j]),
      numeric(nrow(x)))
    colnames(res) <- colnames(x)
    rownames(res) <- rownames(x)
  }
  pmin(pmax(res, 1e-6), 1 - 1e-6)
}

pbsapply <- function(X, FUN, simplify = TRUE, ...) {
  res <- lapply(X, FUN, ...)
  if (!simplify) {
    return(res)
  }
  simplify2array(res)
}

cdf_transform <- function(x, method = c("kde", "ecdf", "resc.ecdf", "adj.ecdf")) {
  if (!is.matrix(x)) {
    stop("Argument 'x' must be a matrix.", call. = FALSE)
  }
  method <- match.arg(method)
  if (method == "kde" && !requireNamespace("ks", quietly = TRUE)) {
    stop("The 'ks' package is required for method = 'kde'.", call. = FALSE)
  }
  n_cols <- ncol(x)
  output <- matrix(NA_real_, nrow = nrow(x), ncol = n_cols)
  colnames(output) <- colnames(x)
  for (j in seq_len(n_cols)) {
    column <- x[, j]
    if (!is.numeric(column)) {
      stop("All columns supplied to 'cdf_transform' must be numeric.", call. = FALSE)
    }
    values <- switch(
      method,
      kde = {
        fit <- ks::kcde(column)
        as.numeric(stats::predict(fit, x = column))
      },
      `resc.ecdf` = pseudo_observations(column),
      `adj.ecdf` = adjusted_pobs(column),
      {
        Fhat <- stats::ecdf(column)
        Fhat(column)
      }
    )
    values <- pmin(pmax(values, 1e-6), 1 - 1e-6)
    output[, j] <- values
  }
  output
}

bootstrap_sample <- function(data) {
  n <- nrow(data)
  data[sample.int(n = n, size = n, replace = TRUE), , drop = FALSE]
}

sanitise_column_names <- function(x) {
  gsub("`", "", x, fixed = TRUE)
}

format_formula_term <- function(x) {
  vapply(
    x,
    function(term) {
      if (grepl("[^[:alnum:]_.]", term)) {
        paste0("`", term, "`")
      } else {
        term
      }
    },
    character(1)
  )
}

combined_formula <- function(response, endogenous, exogenous, include_intercept = TRUE) {
  rhs_terms <- c(endogenous, exogenous)
  rhs_terms <- rhs_terms[rhs_terms != response]
  rhs <- if (length(rhs_terms) == 0) {
    "1"
  } else {
    paste(rhs_terms, collapse = " + ")
  }
  if (!include_intercept && rhs != "1") {
    rhs <- paste0(rhs, " - 1")
  }
  stats::as.formula(paste(response, "~", rhs))
}

extract_components <- function(formula, data) {
  if (!rlang::is_formula(formula)) {
    stop("Argument 'formula' must be a formula.", call. = FALSE)
  }
  if (!is.data.frame(data)) {
    stop("Argument 'data' must be a data.frame.", call. = FALSE)
  }
  split <- nlme::splitFormula(formula, sep = "|")
  response <- all.vars(formula)[1]
  endogenous <- all.vars(split[[1]])
  exogenous <- if (length(split) > 1) all.vars(split[[2]]) else character(0)
  has_intercept <- attr(stats::terms(split[[1]]), "intercept") == 1
  list(
    response = response,
    endogenous = endogenous,
    exogenous = exogenous,
    has_intercept = has_intercept
  )
}

select_complete_cases <- function(data, vars) {
  missing <- setdiff(vars, colnames(data))
  if (length(missing) > 0) {
    stop(
      sprintf(
        "The following variables are missing in the data: %s",
        paste(missing, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  cleaned <- stats::na.omit(data[, vars, drop = FALSE])
  if (nrow(cleaned) == 0) {
    stop("All observations contain missing values after cleaning.", call. = FALSE)
  }
  rownames(cleaned) <- NULL
  as.data.frame(cleaned)
}

validate_numeric_nonconstant <- function(data, vars, allow_constants = FALSE) {
  if (length(vars) == 0) {
    return(invisible(TRUE))
  }
  numeric_vars <- vapply(vars, function(v) is.numeric(data[[v]]), logical(1))
  if (!all(numeric_vars)) {
    stop(
      sprintf(
        "The following variables must be numeric: %s",
        paste(vars[!numeric_vars], collapse = ", ")
      ),
      call. = FALSE
    )
  }
  if (!allow_constants) {
    constant_vars <- vapply(vars, function(v) length(unique(data[[v]])) == 1L, logical(1))
    if (any(constant_vars)) {
      stop(
        sprintf(
          "The following variables are constant and cannot be used: %s",
          paste(vars[constant_vars], collapse = ", ")
        ),
        call. = FALSE
      )
    }
  }
  invisible(TRUE)
}

check_full_rank <- function(matrix, message) {
  if (Matrix::rankMatrix(matrix)[1] != ncol(matrix)) {
    stop(message, call. = FALSE)
  }
  invisible(TRUE)
}

safe_lm <- function(formula, data) {
  stats::lm(formula, data = data)
}

make_result <- function(estimates, residuals, call, method, cdf = NA_character_,
                        boot_replicates = NULL) {
  structure(
    list(
      call = call,
      method = method,
      cdf = cdf,
      coefficients = estimates,
      residuals = residuals,
      bootstrap = boot_replicates
    ),
    class = "endog_copula_fit"
  )
}
