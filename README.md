# endogCopula

<!-- badges: start -->
[![R-CMD-check](https://github.com/ashgreat/endogCopula/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ashgreat/endogCopula/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Package website: <https://ashgreat.github.io/endogCopula/>

`endogCopula` provides instrument-free, Gaussian copula based corrections for
endogenous regressors in linear models. It unifies six recently published
cross-sectional estimators -- Park and Gupta (2012), Yang et al. (2025), Hu
et al. (2025), Breitung et al. (2024), Haschka (2024), and Liengaard et al.
(2025) -- behind a single two-part formula interface with shared input
validation, a choice of marginal CDF estimators, bootstrap standard errors,
and `print()`, `summary()`, `coef()`, `residuals()`, and `confint()` methods.

## Installation

Install the development version from GitHub with `remotes` or `pak`:

```r
# install.packages("remotes")
remotes::install_github("ashgreat/endogCopula")

# or
# install.packages("pak")
pak::pak("ashgreat/endogCopula")
```

## Quick example

```r
library(endogCopula)
data(sim_endog)

# Naive OLS is biased for the endogenous regressor (true coefficient: 2)
ols <- lm(y ~ z_endog + x_exog, data = sim_endog)
coef(ols)["z_endog"]

# Park and Gupta (2012) copula correction
fit <- CopRegPG(
  y ~ z_endog | x_exog + w_instr,
  data = sim_endog,
  cdf = "resc.ecdf",
  nboots = 199
)
summary(fit)
confint(fit)
```

The formula reads `y ~ endogenous | exogenous`: variables before the bar are
treated as continuous endogenous regressors, variables after the bar as
exogenous. Use `- 1` to drop the intercept and `as.factor()` for dummies.
The `cdf` argument selects the marginal CDF estimator (`"kde"`, `"ecdf"`,
`"resc.ecdf"`, or `"adj.ecdf"`); see `vignette("endogCopula")` for guidance.

## Estimators

| Function            | Estimator  | Reference                                                                 |
|---------------------|------------|---------------------------------------------------------------------------|
| `CopRegPG()`        | PG         | Park and Gupta (2012), *Marketing Science* 31(4), 567–586                 |
| `CopReg2sCOPE()`    | 2sCOPE     | Yang, Qian, and Xie (2025), *Journal of Marketing Research* 62(4), 601–623 |
| `CopReg2sCOPEnp()`  | 2sCOPE-np  | Hu, Qian, and Xie (2025), NBER Working Paper 33607                        |
| `CopRegIMA()`       | IMA        | Haschka (2024), *IMA Journal of Management Mathematics* 36(1), 161–180    |
| `CopRegBMW()`       | BMW        | Breitung, Mayer, and Wied (2024), *The Econometrics Journal* 27(3), 362–383 |
| `CopRegJAMS()`      | JAMS       | Liengaard et al. (2025), *Journal of the Academy of Marketing Science* 53, 279–299 |
| `CopRegML_par()`    | PANEL      | Haschka (2022), *Journal of Marketing Research* 59(4), 860–881            |
| `CopRegBayes()`     | BAYES      | Haschka (2022b), SSRN 4235194                                             |

## Numerical equivalence with the published reference code

The estimator implementations are ports of the replication code published
alongside the papers above ("Copula-based endogeneity corrections"). The
package test suite fits each estimator and the corresponding reference
functions on the same data and asserts numerical equivalence of the
resulting coefficient estimates, so refactoring cannot silently change the
statistical results.

## Companion packages

The panel and Bayesian estimators live in dedicated packages:

- [endogCopulaPanel](https://github.com/ashgreat/endogCopulaPanel) --
  fixed-effects copula panel model of Haschka (2022), exported as
  `CopRegML_par()`.
- [endogCopulaBayes](https://github.com/ashgreat/endogCopulaBayes) --
  Bayesian sampler of Haschka (2022b), exported as `CopRegBayes()`.

When these packages are installed, the wrappers of the same names in
`endogCopula` delegate to them.

## Learning more

- `vignette("endogCopula")` -- getting started, worked example, and guidance
  on choosing the `cdf` argument.
- Help pages of the individual estimators (e.g. `?CopRegJAMS`) for
  assumptions, references, and runnable examples.

## License

MIT © Ashwin Malshe
