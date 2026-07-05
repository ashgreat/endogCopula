# Package index

## Package overview

Overview of the endogCopula package.

- [`endogCopula`](https://ashgreat.github.io/endogCopula/reference/endogCopula-package.md)
  [`endogCopula-package`](https://ashgreat.github.io/endogCopula/reference/endogCopula-package.md)
  : endogCopula: Gaussian Copula Based Endogeneity Corrections

## Cross-sectional copula estimators

Gaussian copula based endogeneity corrections for cross-sectional linear
models.

- [`CopRegPG()`](https://ashgreat.github.io/endogCopula/reference/CopRegPG.md)
  : Copula Regression: Park and Gupta (2012)
- [`CopReg2sCOPE()`](https://ashgreat.github.io/endogCopula/reference/CopReg2sCOPE.md)
  : Copula Regression: Yang et al. (2025)
- [`CopReg2sCOPEnp()`](https://ashgreat.github.io/endogCopula/reference/CopReg2sCOPEnp.md)
  : Copula Regression: Hu et al. (2025)
- [`CopRegBMW()`](https://ashgreat.github.io/endogCopula/reference/CopRegBMW.md)
  : Copula Regression: Breitung, Mayer, and Wied (2024)
- [`CopRegIMA()`](https://ashgreat.github.io/endogCopula/reference/CopRegIMA.md)
  : Copula Regression: Haschka (2024)
- [`CopRegJAMS()`](https://ashgreat.github.io/endogCopula/reference/CopRegJAMS.md)
  : Copula Regression: Liengaard et al. (2025)

## Panel and Bayesian estimators (delegated)

Thin wrappers that dispatch to the companion endogCopulaPanel and
endogCopulaBayes packages.

- [`CopRegML_par()`](https://ashgreat.github.io/endogCopula/reference/CopRegML_par.md)
  : Panel Gaussian Copula Estimator (delegated)
- [`CopRegBayes()`](https://ashgreat.github.io/endogCopula/reference/CopRegBayes.md)
  : Bayesian Gaussian Copula Estimator (delegated)

## Methods

Methods for objects returned by the copula estimators.

- [`print(`*`<endog_copula_fit>`*`)`](https://ashgreat.github.io/endogCopula/reference/print.endog_copula_fit.md)
  :

  Print method for `endog_copula_fit`

- [`summary(`*`<endog_copula_fit>`*`)`](https://ashgreat.github.io/endogCopula/reference/summary.endog_copula_fit.md)
  :

  Summary method for `endog_copula_fit`

- [`confint(`*`<endog_copula_fit>`*`)`](https://ashgreat.github.io/endogCopula/reference/confint.endog_copula_fit.md)
  :

  Confidence intervals for `endog_copula_fit` objects

## Data

Simulated data used in examples and the vignette.

- [`sim_endog`](https://ashgreat.github.io/endogCopula/reference/sim_endog.md)
  : Simulated Gaussian Copula Endogeneity Data
