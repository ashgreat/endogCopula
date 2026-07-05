# cran-comments

## Package purpose

endogCopula implements Gaussian copula based estimators for correcting
endogeneity in linear models, including the estimators of Park and Gupta
(2012), Yang, Qian, and Xie (2025), Hu, Qian, and Xie (2025), Breitung,
Mayer, and Wied (2024), Haschka (2024), and Liengaard et al. (2025), along
with wrappers that delegate to companion panel and Bayesian packages.

## Test environments

* local macOS 15 (arm64), R 4.5.0
* GitHub Actions (ubuntu-latest), R release
* GitHub Actions (macos-latest), R release

## R CMD check results

0 errors | 0 warnings | 0 notes

## Submission notes

* This is a new submission.
* All examples that require optional companion packages
  (`endogCopulaPanel`, `endogCopulaBayes`) or the suggested `np` package are
  wrapped in `\donttest{}`/`\dontrun{}` as appropriate; the default
  examples run in well under 5 seconds.
* Point estimates in this package are numerically verified against the
  original published reference implementations for each supported
  estimator; see the package tests.
