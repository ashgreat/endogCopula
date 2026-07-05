# Copula Regression: Liengaard et al. (2025)

Implements the adjusted estimator for Gaussian copula endogeneity
correction proposed by Liengaard et al. (2025). The implementation
closely follows the original research code and supports factor-specific
transformations.

## Usage

``` r
CopRegJAMS(
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

Liengaard, B. D., J.-M. Becker, M. Bennedsen, P. Heiler, L. N. Taylor,
and C. M. Ringle (2025). Dealing with regression models' endogeneity by
means of an adjusted estimator for the Gaussian copula approach.
*Journal of the Academy of Marketing Science* 53, 279–299.

## Examples

``` r
data(sim_endog)

fit <- CopRegJAMS(y ~ z_endog | x_exog + w_instr, data = sim_endog,
                  cdf = "adj.ecdf", nboots = 0)
summary(fit)
#> Gaussian copula estimator: Liengaard et al. (2025)
#> CDF transformation: adj.ecdf
#> 
#> Call:
#> CopRegJAMS(formula = y ~ z_endog | x_exog + w_instr, data = sim_endog, 
#>     cdf = "adj.ecdf", nboots = 0)
#> 
#> Coefficients:
#>                  Estimate
#> `(Intercept)`  0.97746548
#> z_endog        2.09084002
#> x_exog         1.51392827
#> w_instr       -0.05082422
#> z_endog_cop    0.01732337
```
