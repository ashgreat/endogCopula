# endogCopula

`endogCopula` collects recent Gaussian copula based estimators for correcting
endogeneity in linear models. The package unifies code for the Park and Gupta
(2012), Yang et al. (2025), Hu et al. (2025), Breitung et al. (2024), Haschka
(2024), and Liengaard et al. (2025) estimators with consistent interfaces,
documentation, and bootstrap utilities.

## Key features

- A common formula interface, `y ~ endog1 + endog2 | exog1 + exog2`, across
  estimators.
- Shared input validation, copula transformations, and bootstrap helpers.
- S3 methods (`print()`, `summary()`, `coef()`, `residuals()`) for convenient
  post-estimation work.
- Examples and references in the help pages to mirror the original research
  code.

## Getting started

The easiest way to install the development version once the GitHub repository
is available is via

```r
# install.packages("remotes")
remotes::install_github("ashgreat/endogCopula")
```

A minimal example using the Park and Gupta (2012) estimator:

```r
library(endogCopula)

fit <- CopRegPG(
  formula = y ~ price | instrument1 + instrument2,
  data = your_data,
  cdf = "ecdf",
  nboots = 199
)

summary(fit)
```

Please consult the documentation of each estimator (e.g. `?CopRegIMA`) for
assumptions, references, and further examples.
