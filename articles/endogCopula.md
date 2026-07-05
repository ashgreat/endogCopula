# Getting started with endogCopula

## Copula-based endogeneity correction

Endogeneity arises whenever a regressor is correlated with the error
term of a regression model, for example because of omitted variables,
simultaneity, or measurement error. Ordinary least squares is then
biased and inconsistent. The classical remedy is an instrumental
variable, but valid and strong instruments are notoriously hard to find.
Park and Gupta (2012) proposed an instrument-free alternative: model the
joint distribution of the endogenous regressor and the structural error
with a Gaussian copula. Under the assumption that the endogenous
regressor is non-normally distributed while the error is normal, the
dependence between the two can be recovered from the data itself, and
adding a “copula control function” – a transformation of the endogenous
regressor’s estimated cumulative distribution function (CDF) – to the
regression removes the endogeneity bias.

A lively methodological literature has refined this idea. Becker,
Proksch, and Ringle (2022) revisited the original approach and
highlighted the conditions under which it works well, while Qian,
Koschmann, and Xie (2024) provide a practical guide. Newer estimators
address correlation between endogenous and exogenous regressors, which
the original approach ignores: Yang, Qian, and Xie (2025) introduce a
two-stage copula generated regressor approach (2sCOPE), Hu, Qian, and
Xie (2025) develop a nonparametric two-stage variant (2sCOPE-np),
Breitung, Mayer, and Wied (2024) study corrections based on nonlinear
transformations, Haschka (2024) exploits between-regressor correlation,
and Liengaard et al. (2025) propose an adjusted estimator with improved
finite-sample behaviour.

The `endogCopula` package collects these cross-sectional estimators
behind a single formula interface and provides bootstrap standard errors
plus [`print()`](https://rdrr.io/r/base/print.html),
[`summary()`](https://rdrr.io/r/base/summary.html),
[`coef()`](https://rdrr.io/r/stats/coef.html),
[`residuals()`](https://rdrr.io/r/stats/residuals.html), and
[`confint()`](https://rdrr.io/r/stats/confint.html) methods. The
implementations are tested for numerical equivalence against the
published replication code accompanying the papers.

## Estimators and references

| Function | Estimator | Reference |
|----|----|----|
| [`CopRegPG()`](https://ashgreat.github.io/endogCopula/reference/CopRegPG.md) | PG | Park and Gupta (2012) |
| [`CopReg2sCOPE()`](https://ashgreat.github.io/endogCopula/reference/CopReg2sCOPE.md) | 2sCOPE | Yang, Qian, and Xie (2025) |
| [`CopReg2sCOPEnp()`](https://ashgreat.github.io/endogCopula/reference/CopReg2sCOPEnp.md) | 2sCOPE-np | Hu, Qian, and Xie (2025) |
| [`CopRegIMA()`](https://ashgreat.github.io/endogCopula/reference/CopRegIMA.md) | IMA | Haschka (2024) |
| [`CopRegBMW()`](https://ashgreat.github.io/endogCopula/reference/CopRegBMW.md) | BMW | Breitung, Mayer, and Wied (2024) |
| [`CopRegJAMS()`](https://ashgreat.github.io/endogCopula/reference/CopRegJAMS.md) | JAMS | Liengaard et al. (2025) |
| [`CopRegML_par()`](https://ashgreat.github.io/endogCopula/reference/CopRegML_par.md) | PANEL | Haschka (2022), via `endogCopulaPanel` |
| [`CopRegBayes()`](https://ashgreat.github.io/endogCopula/reference/CopRegBayes.md) | BAYES | Haschka (2022b), via `endogCopulaBayes` |

All cross-sectional estimators use a two-part formula,

``` r

y ~ endog1 + endog2 | exog1 + exog2
```

where the variables before the bar are the continuous endogenous
regressors and the variables after the bar are exogenous. Add `- 1` to
drop the intercept and wrap dummies in
[`as.factor()`](https://rdrr.io/r/base/factor.html).

## A worked example

The package ships a simulated dataset, `sim_endog`, generated from

``` math
y = 1 + 2\, z_{\text{endog}} + 1.5\, x_{\text{exog}} + u,
```

where a latent Gaussian factor makes `z_endog` correlated with the
structural error `u`, and `w_instr` shifts `z_endog` without entering
the outcome equation.

``` r

library(endogCopula)
data(sim_endog)
str(sim_endog)
#> 'data.frame':    2000 obs. of  4 variables:
#>  $ y      : num  0.0945 4.2063 -2.7204 3.1972 2.1352 ...
#>  $ z_endog: num  0.3616 0.8219 -1.5403 1.3684 -0.0897 ...
#>  $ x_exog : num  -0.699 0.996 -0.693 -0.103 0.604 ...
#>  $ w_instr: num  1.482 1.705 -0.934 0.605 -0.524 ...
```

A naive OLS regression that ignores the endogeneity of `z_endog`
produces a biased estimate of its coefficient (true value: 2):

``` r

ols <- lm(y ~ z_endog + x_exog, data = sim_endog)
coef(ols)
#> (Intercept)     z_endog      x_exog 
#>   0.9826393   2.0704182   1.5140852
```

The Park and Gupta (2012) estimator adds a copula control function for
`z_endog`. We use `nboots = 25` here so the vignette builds quickly; in
practice several hundred bootstrap replicates are recommended.

``` r

set.seed(42)
fit_pg <- CopRegPG(
  y ~ z_endog | x_exog + w_instr,
  data = sim_endog,
  cdf = "resc.ecdf",
  nboots = 25
)
summary(fit_pg)
#> Gaussian copula estimator: Park and Gupta (2012)
#> CDF transformation: resc.ecdf
#> 
#> Call:
#> CopRegPG(formula = y ~ z_endog | x_exog + w_instr, data = sim_endog, 
#>     cdf = "resc.ecdf", nboots = 25)
#> 
#> Coefficients:
#>                Estimate Std. Error z value Pr(>|z|)    
#> `(Intercept)`  1.004580   0.035742 28.1065  < 2e-16 ***
#> z_endog        1.155322   0.563230  2.0512  0.04024 *  
#> x_exog         1.513418   0.018161 83.3311  < 2e-16 ***
#> w_instr       -0.053370   0.024242 -2.2015  0.02770 *  
#> z_endog_cop    1.024409   0.607497  1.6863  0.09174 .  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Standard errors based on 25 bootstrap replicates.
```

The 2sCOPE estimator of Yang, Qian, and Xie (2025) additionally accounts
for correlation between the endogenous and the exogenous regressors:

``` r

set.seed(42)
fit_2scope <- CopReg2sCOPE(
  y ~ z_endog | x_exog + w_instr,
  data = sim_endog,
  cdf = "resc.ecdf",
  nboots = 25
)
summary(fit_2scope)
#> Gaussian copula estimator: Yang et al. (2025)
#> CDF transformation: resc.ecdf
#> 
#> Call:
#> CopReg2sCOPE(formula = y ~ z_endog | x_exog + w_instr, data = sim_endog, 
#>     cdf = "resc.ecdf", nboots = 25)
#> 
#> Coefficients:
#>               Estimate Std. Error z value Pr(>|z|)    
#> `(Intercept)` 0.990964   0.029324 33.7933  < 2e-16 ***
#> z_endog       1.434730   0.562956  2.5486  0.01082 *  
#> x_exog        1.513242   0.022765 66.4718  < 2e-16 ***
#> w_instr       0.334912   0.329181  1.0174  0.30896    
#> z_endog_cop   0.720331   0.599489  1.2016  0.22953    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Standard errors based on 25 bootstrap replicates.
```

Comparing the estimates of the endogenous coefficient:

``` r

c(
  OLS      = unname(coef(ols)["z_endog"]),
  PG       = unname(coef(fit_pg)["z_endog"]),
  `2sCOPE` = unname(coef(fit_2scope)["z_endog"])
)
#>      OLS       PG   2sCOPE 
#> 2.070418 1.155322 1.434730
```

Confidence intervals based on the bootstrap replicates are available via
[`confint()`](https://rdrr.io/r/stats/confint.html):

``` r

confint(fit_pg, parm = "z_endog")
#>              2.5 %   97.5 %
#> z_endog 0.05141198 2.259232
confint(fit_pg, parm = "z_endog", type = "percentile")
#>             2.5 %   97.5 %
#> z_endog 0.1874727 2.303848
```

## Choosing the `cdf` argument

The copula control function depends on an estimate of each endogenous
regressor’s marginal CDF, selected via the `cdf` argument:

- `"kde"` – the integral of a kernel density estimator, as used by Park
  and Gupta (2012) and Haschka (2022). Requires the suggested `ks`
  package.
- `"ecdf"` – the empirical CDF with its boundary values replaced,
  proposed by Becker, Proksch, and Ringle (2022).
- `"resc.ecdf"` – a rescaled empirical CDF proposed by Qian, Koschmann,
  and Xie (2024) and used by Haschka (2024) and Yang, Qian, and Xie
  (2025).
- `"adj.ecdf"` – an adjusted empirical CDF proposed by Liengaard et al.
  (2025).

A pragmatic default is to use the CDF estimator from the paper that
introduced the estimator you are applying (e.g. `"adj.ecdf"` with
[`CopRegJAMS()`](https://ashgreat.github.io/endogCopula/reference/CopRegJAMS.md),
`"resc.ecdf"` with
[`CopReg2sCOPE()`](https://ashgreat.github.io/endogCopula/reference/CopReg2sCOPE.md)),
and to check that the conclusions are robust across the alternatives.

## Panel and Bayesian estimators

The fixed-effects panel estimator of Haschka (2022) and the Bayesian
sampler of Haschka (2022b) live in dedicated companion packages:

- [`endogCopulaPanel`](https://github.com/ashgreat/endogCopulaPanel) for
  [`CopRegML_par()`](https://ashgreat.github.io/endogCopula/reference/CopRegML_par.md),
  the maximum likelihood panel estimator.
- [`endogCopulaBayes`](https://github.com/ashgreat/endogCopulaBayes) for
  [`CopRegBayes()`](https://ashgreat.github.io/endogCopula/reference/CopRegBayes.md),
  the MCMC sampler.

When these packages are installed, the wrappers of the same names
exported by `endogCopula` delegate to them.

## References

- Becker, J.-M., D. Proksch, and C. M. Ringle (2022). Revisiting
  Gaussian copulas to handle endogenous regressors. *Journal of the
  Academy of Marketing Science* 50, 46–66.
- Breitung, J., A. Mayer, and D. Wied (2024). Asymptotic properties of
  endogeneity corrections using nonlinear transformations. *The
  Econometrics Journal* 27(3), 362–383.
- Haschka, R. E. (2022). Handling endogenous regressors using copulas: A
  generalisation to linear panel models with fixed effects and
  correlated regressors. *Journal of Marketing Research* 59(4), 860–881.
- Haschka, R. E. (2022b). Bayesian inference for joint estimation models
  using copulas to handle endogenous regressors.
  <https://ssrn.com/abstract=4235194>
- Haschka, R. E. (2024). Robustness of copula-correction models in
  causal analysis: Exploiting between-regressor correlation. *IMA
  Journal of Management Mathematics* 36(1), 161–180.
- Hu, X., Y. Qian, and H. Xie (2025). Correcting endogeneity via
  instrument-free two-stage nonparametric copula control functions. NBER
  Working Paper 33607. <http://www.nber.org/papers/w33607>
- Liengaard, B. D., J.-M. Becker, M. Bennedsen, P. Heiler, L. N. Taylor,
  and C. M. Ringle (2025). Dealing with regression models’ endogeneity
  by means of an adjusted estimator for the Gaussian copula approach.
  *Journal of the Academy of Marketing Science* 53, 279–299.
- Park, S. and S. Gupta (2012). Handling endogenous regressors by joint
  estimation using copulas. *Marketing Science* 31(4), 567–586.
- Qian, Y., A. Koschmann, and H. Xie (2024). A practical guide to
  endogeneity correction using copulas. NBER Working Paper 32231.
  <https://www.nber.org/papers/w32231>
- Yang, F., Y. Qian, and H. Xie (2025). Addressing endogeneity using a
  two-stage copula generated regressor approach. *Journal of Marketing
  Research* 62(4), 601–623.
