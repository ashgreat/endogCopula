# Copula Regression: Hu et al. (2025)

Implements the instrument-free two-stage copula control function
estimator proposed by Hu et al. (2025). When exogenous regressors are
present, the first stage estimates conditional distribution functions
nonparametrically, which requires the suggested `np` package.

## Usage

``` r
CopReg2sCOPEnp(
  formula,
  data,
  cdf = c("kde", "ecdf", "resc.ecdf", "adj.ecdf"),
  nboots = 199
)
```

## Arguments

- formula:

  A two-part formula of the form `y ~ endog1 + endog2 | exog1 + exog2`,
  where the left-hand side lists the endogenous regressors and the
  right-hand side after the bar lists exogenous regressors. Use `- 1` to
  remove the intercept.

- data:

  A `data.frame` containing the variables referenced in `formula`.

- cdf:

  Character string specifying the cumulative distribution function
  estimator. One of `"kde"`, `"ecdf"`, `"resc.ecdf"`, or `"adj.ecdf"`.

- nboots:

  Number of bootstrap replications used to compute standard errors. Set
  to zero to skip bootstrapping.

## Value

An object of class `endog_copula_fit`: a list with components

- `call`:

  The matched call.

- `method`:

  Character string identifying the estimator.

- `cdf`:

  The CDF estimator used for the copula transformation.

- `coefficients`:

  Matrix of coefficient estimates with an `Estimate` column and, when
  `nboots > 0`, a `Std. Error` column computed from the bootstrap
  replicates.

- `residuals`:

  Numeric vector of residuals from the structural part of the model
  (copula control terms excluded).

- `bootstrap`:

  Matrix of bootstrap coefficient replicates, or `NULL` when
  `nboots = 0`.

Objects of this class support
[`print()`](https://rdrr.io/r/base/print.html),
[`summary()`](https://rdrr.io/r/base/summary.html),
[`coef()`](https://rdrr.io/r/stats/coef.html),
[`residuals()`](https://rdrr.io/r/stats/residuals.html), and
[`confint()`](https://rdrr.io/r/stats/confint.html).

## References

Hu, X., Y. Qian, and H. Xie (2025). Correcting endogeneity via
instrument-free two-stage nonparametric copula control functions.
National Bureau of Economic Research Working Paper 33607.

## Examples

``` r
data(sim_endog)

# Without exogenous regressors the estimator does not need the 'np' package
fit <- CopReg2sCOPEnp(y ~ z_endog, data = sim_endog,
                      cdf = "resc.ecdf", nboots = 0)
summary(fit)
#> Gaussian copula estimator: Hu et al. (2025)
#> CDF transformation: resc.ecdf
#> 
#> Call:
#> CopReg2sCOPEnp(formula = y ~ z_endog, data = sim_endog, cdf = "resc.ecdf", 
#>     nboots = 0)
#> 
#> Coefficients:
#>                 Estimate
#> `(Intercept)`  1.1946512
#> z_endog        2.1882798
#> z_endog_cop   -0.4709731

# \donttest{
# With exogenous regressors the nonparametric first stage uses 'np'
if (requireNamespace("np", quietly = TRUE)) {
  fit_np <- CopReg2sCOPEnp(y ~ z_endog | x_exog,
                           data = sim_endog[1:150, ],
                           cdf = "resc.ecdf", nboots = 0)
  coef(fit_np)
}
#> `(Intercept)`       z_endog        x_exog   z_endog_cop 
#>     0.9153671    -2.1121331     1.5464786     4.6185621 
# }
```
