# Panel Gaussian Copula Estimator (delegated)

Thin wrapper around `endogCopulaPanel::CopRegML_par()`. The function is
only available when the optional `endogCopulaPanel` package is
installed.

## Usage

``` r
CopRegML_par(
  formula,
  index,
  data,
  ecdf = TRUE,
  nboots = 199,
  starting_values = NULL,
  method = "Nelder-Mead"
)
```

## Arguments

- formula:

  Two-part formula `y ~ endog | exog` describing the model.

- index:

  Character vector of length two with panel and time identifiers.

- data:

  Data frame containing the referenced variables.

- ecdf:

  Logical; use empirical CDFs (`TRUE`) or kernel estimates (`FALSE`).

- nboots:

  Number of bootstrap replications.

- starting_values:

  Optional numeric vector of starting values.

- method:

  Optimisation method passed to
  [`stats::optim()`](https://rdrr.io/r/stats/optim.html).

## Value

The fitted panel model as returned by
`endogCopulaPanel::CopRegML_par()`: maximum likelihood coefficient
estimates for the fixed-effects copula panel model together with
bootstrap standard errors when `nboots > 0`. See the documentation in
the `endogCopulaPanel` package for the full structure.

## Examples

``` r
if (FALSE) { # \dontrun{
# Requires the optional 'endogCopulaPanel' package and panel data with
# numeric id and time variables:
fit <- CopRegML_par(y ~ z_endog | x_exog + as.factor(year),
                    index = c("id", "year"), data = panel_data,
                    nboots = 199)
} # }
```
