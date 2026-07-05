# Print method for `endog_copula_fit`

Print method for `endog_copula_fit`

## Usage

``` r
# S3 method for class 'endog_copula_fit'
print(x, ...)
```

## Arguments

- x:

  An object returned by one of the copula estimators.

- ...:

  Passed to
  [`stats::printCoefmat()`](https://rdrr.io/r/stats/printCoefmat.html)
  when bootstrap standard errors are available.

## Value

Invisibly returns `x`, after printing the estimator name, the CDF
transformation, and the coefficient table.

## Examples

``` r
data(sim_endog)
fit <- CopRegPG(y ~ z_endog | x_exog + w_instr, data = sim_endog,
                cdf = "resc.ecdf", nboots = 0)
print(fit)
#> Gaussian copula estimator: Park and Gupta (2012)
#> CDF transformation: resc.ecdf
#> 
#> Coefficients:
#>                  Estimate
#> `(Intercept)`  1.00458029
#> z_endog        1.15532197
#> x_exog         1.51341797
#> w_instr       -0.05337007
#> z_endog_cop    1.02440912
```
