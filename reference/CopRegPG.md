# Copula Regression: Park and Gupta (2012)

Implements the Gaussian copula approach proposed by Park and Gupta
(2012) for correcting endogeneity in linear models. Bootstrap standard
errors are obtained via resampling.

## Usage

``` r
CopRegPG(
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

## Details

Point estimates are numerically equivalent to the published reference
implementation. Bootstrap standard errors deliberately differ: the
reference code omits the inverse-normal (`qnorm`) transformation of the
copula terms inside bootstrap replicates even though its main fit
applies it. Here every replicate refits the same estimator as the main
fit, so the standard errors describe the estimator actually reported.

## References

Park, S. and S. Gupta (2012). Handling endogenous regressors by joint
estimation using copulas. *Marketing Science* 31 (4), 567–586.

## Examples

``` r
data(sim_endog)

fit <- CopRegPG(y ~ z_endog | x_exog + w_instr, data = sim_endog,
                cdf = "resc.ecdf", nboots = 0)
summary(fit)
#> Gaussian copula estimator: Park and Gupta (2012)
#> CDF transformation: resc.ecdf
#> 
#> Call:
#> CopRegPG(formula = y ~ z_endog | x_exog + w_instr, data = sim_endog, 
#>     cdf = "resc.ecdf", nboots = 0)
#> 
#> Coefficients:
#>                  Estimate
#> `(Intercept)`  1.00458029
#> z_endog        1.15532197
#> x_exog         1.51341797
#> w_instr       -0.05337007
#> z_endog_cop    1.02440912

if (requireNamespace("ks", quietly = TRUE)) {
  fit_kde <- CopRegPG(y ~ z_endog | x_exog + w_instr, data = sim_endog,
                      cdf = "kde", nboots = 0)
  coef(fit_kde)
}
#> `(Intercept)`       z_endog        x_exog       w_instr   z_endog_cop 
#>    1.01818468    1.01531666    1.51358948   -0.05371194    1.18281589 
```
