# Bayesian Gaussian Copula Estimator (delegated)

Thin wrapper around `endogCopulaBayes::CopRegBayes()`. Requires the
optional `endogCopulaBayes` package to be installed.

## Usage

``` r
CopRegBayes(
  data,
  iterations = 10000,
  burnin = 2000,
  thin = 10,
  startvalue = NULL
)
```

## Arguments

- data:

  Data frame containing the columns `y`, `z`, and `x`.

- iterations:

  Total number of MCMC iterations.

- burnin:

  Number of initial iterations discarded as burn-in.

- thin:

  Thinning interval applied after burn-in.

- startvalue:

  Optional numeric vector of starting values.

## Value

The MCMC output as returned by `endogCopulaBayes::CopRegBayes()`,
containing posterior draws for the model parameters after burn-in and
thinning. See the documentation in the `endogCopulaBayes` package for
the full structure.

## Examples

``` r
if (FALSE) { # \dontrun{
# Requires the optional 'endogCopulaBayes' package:
data(sim_endog)
bayes_data <- data.frame(y = sim_endog$y, z = sim_endog$z_endog,
                         x = sim_endog$x_exog)
draws <- CopRegBayes(bayes_data, iterations = 10000, burnin = 2000)
} # }
```
