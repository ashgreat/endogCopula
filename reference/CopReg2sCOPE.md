# Copula Regression: Yang et al. (2025)

Implements the two-stage copula generated regressor estimator proposed
by Yang et al. (2025).

## Usage

``` r
CopReg2sCOPE(
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

Yang, F., Y. Qian, and H. Xie (2025). Addressing endogeneity using a
two-stage copula generated regressor approach. *Journal of Marketing
Research* 62(4), 601–623.

## Examples

``` r
data(sim_endog)

fit <- CopReg2sCOPE(y ~ z_endog | x_exog + w_instr, data = sim_endog,
                    cdf = "resc.ecdf", nboots = 0)
summary(fit)
#> Gaussian copula estimator: Yang et al. (2025)
#> CDF transformation: resc.ecdf
#> 
#> Call:
#> CopReg2sCOPE(formula = y ~ z_endog | x_exog + w_instr, data = sim_endog, 
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
