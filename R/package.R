#' endogCopula: Gaussian Copula Based Endogeneity Corrections
#'
#' The endogCopula package implements a collection of copula-based estimators
#' for addressing endogeneity in linear models. It includes cross-sectional
#' estimators from Park and Gupta (2012) ([CopRegPG()]), Yang et al. (2025)
#' ([CopReg2sCOPE()]), Hu et al. (2025) ([CopReg2sCOPEnp()]), Breitung et al.
#' (2024) ([CopRegBMW()]), Haschka (2024) ([CopRegIMA()]), and Liengaard et
#' al. (2025) ([CopRegJAMS()]).
#'
#' All cross-sectional estimators share the two-part formula interface
#' `y ~ endog1 + endog2 | exog1 + exog2`, a `cdf` argument selecting the
#' marginal CDF estimator (`"kde"`, `"ecdf"`, `"resc.ecdf"`, or
#' `"adj.ecdf"`), and bootstrap standard errors controlled via `nboots`.
#' Fitted models are returned as `endog_copula_fit` objects with `print()`,
#' `summary()`, `coef()`, `residuals()`, and `confint()` methods.
#'
#' The panel estimator of Haschka (2022) and the Bayesian sampler of Haschka
#' (2022b) live in the companion packages `endogCopulaPanel` and
#' `endogCopulaBayes`; the wrappers [CopRegML_par()] and [CopRegBayes()]
#' delegate to them when those packages are installed.
#'
#' @keywords internal
"_PACKAGE"
