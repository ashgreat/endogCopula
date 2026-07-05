# Simulated Gaussian Copula Endogeneity Data

A synthetic dataset generated to illustrate Gaussian copula endogeneity
corrections. The dependent variable is driven by an endogenous regressor
that is correlated with the structural error via a shared noise
component.

## Usage

``` r
sim_endog
```

## Format

A data frame with 2000 rows and 4 variables:

- y:

  Outcome variable

- z_endog:

  Endogenous regressor that is correlated with the error term

- x_exog:

  Exogenous regressor

- w_instr:

  Instrument used by copula estimators

## Details

The true data-generating process is \$\$y = 1 + 2 z\_{\text{endog}} +
1.5 x\_{\text{exog}} + u,\$\$ where a latent Gaussian factor induces
correlation between \\z\_{\text{endog}}\\ and the structural error
\\u\\. The instrument \\w\_{\text{instr}}\\ shifts the endogenous
regressor but is excluded from the outcome equation.

## Examples

``` r
data(sim_endog)
str(sim_endog)
#> 'data.frame':    2000 obs. of  4 variables:
#>  $ y      : num  0.0945 4.2063 -2.7204 3.1972 2.1352 ...
#>  $ z_endog: num  0.3616 0.8219 -1.5403 1.3684 -0.0897 ...
#>  $ x_exog : num  -0.699 0.996 -0.693 -0.103 0.604 ...
#>  $ w_instr: num  1.482 1.705 -0.934 0.605 -0.524 ...
```
