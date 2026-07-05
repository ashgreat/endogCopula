# Summary method for `endog_copula_fit`

Summary method for `endog_copula_fit`

## Usage

``` r
# S3 method for class 'endog_copula_fit'
summary(object, ...)
```

## Arguments

- object:

  An object returned by one of the copula estimators.

- ...:

  Passed to downstream methods.

## Value

An object of class `summary.endog_copula_fit`: a list carrying the
original `call`, `method`, `cdf`, and `coefficients`, plus a
`coefficient_table` with z statistics and p-values (or `NULL` when no
bootstrap standard errors are available), the `residuals`, the matrix of
`bootstrap` replicates, and the number of bootstrap draws `nboots`.

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
```
