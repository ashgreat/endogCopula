#' Simulated Gaussian Copula Endogeneity Data
#'
#' A synthetic dataset generated to illustrate Gaussian copula endogeneity
#' corrections. The dependent variable is driven by an endogenous regressor
#' that is correlated with the structural error via a shared noise component.
#'
#' @format A data frame with 2000 rows and 4 variables:
#' \describe{
#'   \item{y}{Outcome variable}
#'   \item{z_endog}{Endogenous regressor that is correlated with the error term}
#'   \item{x_exog}{Exogenous regressor}
#'   \item{w_instr}{Instrument used by copula estimators}
#' }
#' @details The true data-generating process is
#'   \deqn{y = 1 + 2 z_{\text{endog}} + 1.5 x_{\text{exog}} + u,}
#'   where a latent Gaussian factor induces correlation between
#'   \eqn{z_{\text{endog}}} and the structural error \eqn{u}. The instrument
#'   \eqn{w_{\text{instr}}} shifts the endogenous regressor but is excluded from
#'   the outcome equation.
#' @examples
#' data(sim_endog)
#' str(sim_endog)
"sim_endog"
