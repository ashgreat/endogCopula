# Copula Regression: Haschka (2024)

Implements the estimator proposed by Haschka (2024), which exploits
correlations between regressors to build copula-based control functions.

## Usage

``` r
CopRegIMA(
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

Haschka, R. E. (2024). Robustness of copula-correction models in causal
analysis: Exploiting between-regressor correlation. *IMA Journal of
Management Mathematics* 36 (1), 161–180.

## Examples

``` r
data(sim_endog)

fit <- CopRegIMA(y ~ z_endog | x_exog + w_instr, data = sim_endog,
                 cdf = "resc.ecdf", nboots = 0)
summary(fit)
#> Gaussian copula estimator: Haschka (2024)
#> CDF transformation: resc.ecdf
#> 
#> Call:
#> CopRegIMA(formula = y ~ z_endog | x_exog + w_instr, data = sim_endog, 
#>     cdf = "resc.ecdf", nboots = 0)
#> 
#> Coefficients:
#>                Estimate
#> `(Intercept)` 0.9909636
#> z_endog       1.4347301
#> x_exog        1.5132420
#> w_instr       0.3349117
#> z_endog_cop   0.7203308
```
