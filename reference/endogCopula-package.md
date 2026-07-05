# endogCopula: Gaussian Copula Based Endogeneity Corrections

The endogCopula package implements a collection of copula-based
estimators for addressing endogeneity in linear models. It includes
cross-sectional estimators from Park and Gupta (2012)
([`CopRegPG()`](https://ashgreat.github.io/endogCopula/reference/CopRegPG.md)),
Yang et al. (2025)
([`CopReg2sCOPE()`](https://ashgreat.github.io/endogCopula/reference/CopReg2sCOPE.md)),
Hu et al. (2025)
([`CopReg2sCOPEnp()`](https://ashgreat.github.io/endogCopula/reference/CopReg2sCOPEnp.md)),
Breitung et al. (2024)
([`CopRegBMW()`](https://ashgreat.github.io/endogCopula/reference/CopRegBMW.md)),
Haschka (2024)
([`CopRegIMA()`](https://ashgreat.github.io/endogCopula/reference/CopRegIMA.md)),
and Liengaard et al. (2025)
([`CopRegJAMS()`](https://ashgreat.github.io/endogCopula/reference/CopRegJAMS.md)).

## Details

All cross-sectional estimators share the two-part formula interface
`y ~ endog1 + endog2 | exog1 + exog2`, a `cdf` argument selecting the
marginal CDF estimator (`"kde"`, `"ecdf"`, `"resc.ecdf"`, or
`"adj.ecdf"`), and bootstrap standard errors controlled via `nboots`.
Fitted models are returned as `endog_copula_fit` objects with
[`print()`](https://rdrr.io/r/base/print.html),
[`summary()`](https://rdrr.io/r/base/summary.html),
[`coef()`](https://rdrr.io/r/stats/coef.html),
[`residuals()`](https://rdrr.io/r/stats/residuals.html), and
[`confint()`](https://rdrr.io/r/stats/confint.html) methods.

The panel estimator of Haschka (2022) and the Bayesian sampler of
Haschka (2022b) live in the companion packages `endogCopulaPanel` and
`endogCopulaBayes`; the wrappers
[`CopRegML_par()`](https://ashgreat.github.io/endogCopula/reference/CopRegML_par.md)
and
[`CopRegBayes()`](https://ashgreat.github.io/endogCopula/reference/CopRegBayes.md)
delegate to them when those packages are installed.

## See also

Useful links:

- <https://ashgreat.github.io/endogCopula/>

- <https://github.com/ashgreat/endogCopula>

- Report bugs at <https://github.com/ashgreat/endogCopula/issues>

## Author

**Maintainer**: Ashwin Malshe <ashwin@malshe.com>
