#' maczic: Mediation Analysis for Count and Zero-inflated Count Data.
#'
#' Performs causal mediation analysis for count and zero-inflated count data;
#' calculates power of causal mediation effects, direct effects, total effects;
#' implements Instrumental Variable (IV) method to estimate the controlled
#' (natural) direct and mediation effect, and the bootstrap CIs; calculates
#' power of controlled (natural) direct and mediation effect using IV method;
#' performs sensitivity analysis when there is a treatment-induced
#' mediator-outcome confounder by varing coefficient of treatment in the
#' confounder model and plots sensitivity graph.
#'
#' @aliases maczic-package
#'
#' @name maczic
#'
#' @importFrom MASS mvrnorm gamma.shape polr rnegbin glm.nb
#' @importFrom graphics matplot legend
#' @importFrom sandwich vcovHC estfun sandwich
#' @importFrom pscl zeroinfl
#' @importFrom survival survreg.distributions
#' @importFrom stats binomial coef delete.response getCall glm lm median
#'   model.frame model.matrix model.response model.weights plogis pnorm
#'   poisson predict quantile rbinom rgamma rlogis rmultinom rnbinom rnorm
#'   rpois rt runif rweibull terms vcov weighted.mean update
#' @importFrom mediation mediate
#' @importFrom emplik el.test
#' @importFrom BB BBoptim
#' @import mathjaxr
NULL
