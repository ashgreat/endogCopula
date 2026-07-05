# Confidence intervals for `endog_copula_fit` objects

Computes confidence intervals for the estimated coefficients. The
default `type = "normal"` uses a normal approximation based on the
bootstrap standard errors; `type = "percentile"` uses percentiles of the
bootstrap replicates stored in `object$bootstrap`. Both require the
model to have been fitted with `nboots > 0`.

## Usage

``` r
# S3 method for class 'endog_copula_fit'
confint(object, parm, level = 0.95, type = c("normal", "percentile"), ...)
```

## Arguments

- object:

  An object returned by one of the copula estimators.

- parm:

  Coefficients to compute intervals for, either names or indices.
  Defaults to all coefficients.

- level:

  Confidence level. Defaults to `0.95`.

- type:

  Either `"normal"` (normal approximation from bootstrap standard
  errors) or `"percentile"` (percentile intervals from the bootstrap
  replicates).

- ...:

  Currently unused.

## Value

A matrix with one row per coefficient and columns giving the lower and
upper confidence limits, labelled with the corresponding percentiles.
Coefficients without bootstrap information are filled with `NA`.

## Examples

``` r
data(sim_endog)
fit <- CopRegPG(y ~ z_endog | x_exog + w_instr, data = sim_endog,
                cdf = "resc.ecdf", nboots = 25)
#> Computing bootstrap standard errors (this may take a while)...
confint(fit)
#>                     2.5 %      97.5 %
#> `(Intercept)`  0.94948822 1.059672360
#> z_endog       -0.01926163 2.329905570
#> x_exog         1.46391015 1.562925778
#> w_instr       -0.10831481 0.001574675
#> z_endog_cop   -0.24348216 2.292300395
confint(fit, parm = "z_endog", type = "percentile")
#>             2.5 %  97.5 %
#> z_endog 0.2974968 2.25039
```
