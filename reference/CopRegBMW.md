# Copula Regression: Breitung, Mayer, and Wied (2024)

Implements the estimator proposed by Breitung, Mayer, and Wied (2024)
that combines classical first-stage regressions with copula-based
control functions. With exogenous regressors, each endogenous regressor
is first regressed on the exogenous ones and the copula control function
is built from the CDF-transformed (and normal-quantile mapped)
first-stage residuals. Without exogenous regressors the CDF-transformed
endogenous regressors enter the control-function regression directly.

## Usage

``` r
CopRegBMW(
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

Breitung, J., A. Mayer, and D. Wied (2024). Asymptotic properties of
endogeneity corrections using nonlinear transformations. *The
Econometrics Journal* 27 (3), 362–383.

## Examples

``` r
data(sim_endog)

fit <- CopRegBMW(y ~ z_endog | x_exog + w_instr, data = sim_endog,
                 cdf = "resc.ecdf", nboots = 0)
summary(fit)
#> Gaussian copula estimator: Breitung, Mayer, and Wied (2024)
#> CDF transformation: resc.ecdf
#> 
#> Call:
#> CopRegBMW(formula = y ~ z_endog | x_exog + w_instr, data = sim_endog, 
#>     cdf = "resc.ecdf", nboots = 0)
#> 
#> Coefficients:
#>                 Estimate
#> `(Intercept)`  0.9767925
#> z_endog        2.6199074
#> x_exog         1.5141439
#> w_instr       -0.3625772
#> z_endog_cop   -0.4798474
```
