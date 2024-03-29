% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mediate_zi.R
\name{mediate_zi}
\alias{mediate_zi}
\title{Mediation Analysis for Count and Zero-Inflated Count Data}
\usage{
mediate_zi(
  model.m,
  model.y,
  sims = 1000,
  boot = FALSE,
  treat = "treat.name",
  mediator = "med.name",
  covariates = NULL,
  outcome = NULL,
  control = NULL,
  conf.level = 0.95,
  control.value = 0,
  treat.value = 1,
  long = TRUE,
  dropobs = FALSE,
  robustSE = FALSE,
  cluster = NULL,
  ...
)
}
\arguments{
\item{model.m}{A fitted model object for mediator.  Can be of class 'lm',
'polr', 'bayespolr', 'glm', 'bayesglm', 'gam', 'rq', or 'survreg'.}

\item{model.y}{A fitted model object for outcome.  Can be of class 'lm',
'polr', 'bayespolr', 'glm', 'bayesglm', 'gam', 'vglm', 'rq', 'survreg',
or 'zeroinfl'}

\item{sims}{Number of Monte Carlo draws for nonparametric bootstrap or
quasi-Bayesian approximation.}

\item{boot}{A logical value. if 'FALSE' a quasi-Bayesian approximation is
used for confidence intervals; if 'TRUE' nonparametric bootstrap will be
used. Default is 'FALSE'.}

\item{treat}{A character string indicating the name of the treatment variable
used in the models.  The treatment can be either binary (integer or a
two-valued factor) or continuous (numeric).}

\item{mediator}{A character string indicating the name of the mediator
variable used in the models.}

\item{covariates}{A list or data frame containing values for a subset of the
pre-treatment covariates in 'model.m' and 'model.y'. If provided, the
function will return the estimates conditional on those covariate values.}

\item{outcome}{A character string indicating the name of the outcome variable
in `model.y'. Only necessary if 'model.y' is of class 'survreg'; otherwise
ignored.}

\item{control}{A character string indicating the name of the control group
indicator. Only relevant if 'model.y' is of class 'gam'. If provided, 'd0',
'z0' and 'n0' are allowed to differ from 'd1', 'z1' and 'n1', respectively.}

\item{conf.level}{Level of the returned two-sided confidence intervals.
Default is to return the 2.5 and 97.5 percentiles of the simulated
quantities.}

\item{control.value}{Value of the treatment variable used as the control
condition. Default is 0.}

\item{treat.value}{Value of the treatment variable used as the treatment
condition. Default is 1.}

\item{long}{A logical value. If 'TRUE', the output will contain the entire
sets of simulation draws of the the average causal mediation effects,
direct effects, proportions mediated, and total effect. Default is 'TRUE'.}

\item{dropobs}{A logical value indicating the behavior when the model frames
of 'model.m' and 'model.y' (and the 'cluster' variable if included) are
composed of different observations. If 'TRUE', models will be re-fitted
using common data rows. If 'FALSE', error is returned. Default is 'FALSE'.}

\item{robustSE}{A logical value. If 'TRUE', heteroskedasticity-consistent
standard errors will be used in quasi-Bayesian simulations. Ignored if
'boot' is 'TRUE' or neither 'model.m' nor 'model.y' has a method for
\code{vcovHC} in the \code{sandwich} package. Default is 'FALSE'.}

\item{cluster}{A variable indicating clusters for standard errors. Note that
this should be a vector of cluster indicators itself, not a character
string for the name of the variable.}

\item{...}{other arguments passed to \code{vcovHC} in the \code{sandwich}
package: typically the 'type' argument, which is ignored if 'robustSE' is
'FALSE'. Arguments to the \code{boot} in the \code{boot} package may also
be passed, e.g. 'parallel' and 'ncpus'.}
}
\value{
\code{mediate} returns an object of class "\code{mediate}", (or
  "\code{mediate.order}" if the outcome model used is 'polr' or 'bayespolr'),
  a list that contains the components listed below.  Some of these elements
  are not available if 'long' is set to 'FALSE' by the user.

  The function \code{summary} (i.e., \code{summary.mediate},
  or \code{summary.mediate.order}) can be used to obtain a table of the
  results.  The function \code{plot} (i.e., \code{plot.mediate}, or
  \code{plot.mediate.order}) can be used to produce a plot of the estimated
  average causal mediation, average direct, and total effects along with
  their confidence intervals.

  \item{d0, d1}{point estimates for average causal mediation effects under
  the control and treatment conditions.}
  \item{d0.ci, d1.ci}{confidence intervals for average causal mediation
  effects. The confidence level is set at the value specified in
  'conf.level'.}
  \item{d0.p, d1.p}{two-sided p-values for average causal mediation effects.}
  \item{d0.sims, d1.sims}{vectors of length 'sims' containing simulation
  draws of average causal mediation effects.}
  \item{z0, z1}{point estimates for average direct effect under the control
  and treatment conditions.}
  \item{z0.ci, z1.ci}{confidence intervals for average direct effects.}
  \item{z0.p, z1.p}{two-sided p-values for average causal direct effects.}
  \item{z0.sims, z1.sims}{vectors of length 'sims' containing simulation
  draws of average direct effects.}
  \item{n0, n1}{the "proportions mediated", or the size of the average causal
  mediation effects relative to the total effect.}
  \item{n0.ci, n1.ci}{confidence intervals for the proportions mediated.}
  \item{n0.p, n1.p}{two-sided p-values for proportions mediated.}
  \item{n0.sims, n1.sims}{vectors of length 'sims' containing simulation
  draws of the proportions mediated.}
  \item{tau.coef}{point estimate for total effect.}
  \item{tau.ci}{confidence interval for total effect.}
  \item{tau.p}{two-sided p-values for total effect.}
  \item{tau.sims}{a vector of length 'sims' containing simulation draws of
  the total effect.}
  \item{d.avg, z.avg, n.avg}{simple averages of d0 and d1, z0 and z1, n0 and
  n1, respectively, which users may want to use as summary values when those
  quantities differ.}
  \item{d.avg.ci, z.avg.ci, n.avg.ci}{confidence intervals for the above.}
  \item{d.avg.p, z.avg.p, n.avg.p}{two-sided p-values for the above.}
  \item{d.avg.sims, z.avg.sims, n.avg.sims}{vectors of length 'sims'
  containing simulation draws of d.avg, z.avg and n.avg, respectively.}
  \item{boot}{logical, the 'boot' argument used. If 'FALSE' a quasi-Bayesian
  approximation was used for confidence intervals; if 'TRUE' nonparametric
  bootstrap was used}
  \item{boot.ci.type}{a character string 'perc' indicating percentile
  bootstrap confidence intervals were estimated if the argument boot = TRUE}
  \item{treat}{a character string indicating the name of the 'treat' variable
  used in the models}
  \item{mediator}{a character string indicating the name of the 'mediator'
  variable used in the models}
  \item{INT}{a logical value indicating whether the model specification
  allows the effects to differ between the treatment and control conditions.}
  \item{conf.level}{the confidence level used. }
  \item{model.y}{the outcome model used.}
  \item{model.m}{the mediator model used.}
  \item{control.value}{value of the treatment variable used as the control
  condition.}
  \item{treat.value}{value of the treatment variable used as the treatment
  condition.}
  \item{nobs}{number of observations in the model frame for 'model.m' and
  'model.y'. May differ from the numbers in the original models input to
  'mediate' if 'dropobs' was 'TRUE'.}
  \item{robustSE}{`TRUE' or `FALSE'.}
  \item{cluster}{the clusters used.}
}
\description{
\loadmathjax{} 'mediate_zi' is modified from the \code{mediate} function in
\code{mediation} package (version 4.0.1) with the new feature
added to handle Zero-inflated count outcomes including Zero-inflated
Poisson and Zero-inflated Negative Binomial outcomes. This function is used
to estimate various quantities for causal mediation analysis, including
average causal mediation effects  (ACME) or indirect effect, average direct
effects (ADE), proportions mediated, and total effect. The Usage, Argument,
Value, and Details of this function are the same as \code{mediate()} in
\code{mediation} package version 4.0.1 (see below) except for the description
 of the 'model.y' argument, which has class 'zeroinfl' added.
}
\details{
This is the workhorse function for estimating causal mediation
  effects for a variety of data types. The average causal mediation effect
  (ACME) represents the expected difference in the potential outcome when the
  mediator took the value that would realize under the treatment condition as
  opposed to the control condition, while the treatment status itself is held
  constant. That is,
  \mjdeqn{\delta(t) \ = \ E[Y(t, M(t_1)) - Y(t, M(t_0))],}{%
        \delta(t) = E[Y(t, M(t1)) - Y(t, M(t0))],}
  where \mjeqn{t, t_1, t_0}{t, t1, t0} are particular values of the treatment
  \mjseqn{T} such that \mjeqn{t_1 \neq t_0}{t1 != t0}, \mjseqn{M(t)} is the potential
  mediator, and \mjseqn{Y(t,m)} is the potential outcome variable. The average
  direct effect (ADE) is defined similarly as,
  \mjdeqn{\zeta(t) \ = \ E[Y(t_1, M(t)) - Y(t_0, M(t))],}{%
        \zeta(t) = E[Y(t1, M(t)) - Y(t0, M(t))],}
  which represents the expected difference in the potential outcome when the
  treatment is changed but the mediator is held constant at the value that
  would realize if the treatment equals \mjseqn{t}. The two quantities on
  average add up to the total effect of the treatment on the outcome,
  \mjseqn{\tau}. See the references for more details.

  When both the mediator model ('model.m') and outcome model ('model.y') are
  normal linear regressions, the results will be identical to the usual LSEM
  method by Baron and Kenny (1986).  The function can, however, accommodate
  other data types including binary, ordered and count outcomes and mediators
  as well as censored outcomes.  Variables can also be modeled
  nonparametrically, semiparametrically, or using quantile regression.

  If it is desired that inference be made conditional on specific values of
  the pre-treatment covariates included in the model, the `covariates'
  argument can be used to set those values as a list or data frame. The list
  may contain either the entire set or any strict subset of the covariates in
  the model; in the latter case, the effects will be averaged over the other
  covariates. The `covariates' argument will be particularly useful when the
  models contain interactions between the covariates and the treatment and/or
  mediator (known as ``moderated mediation'').

  The prior weights in the mediator and outcome models are taken as sampling
  weights and the estimated effects will be weighted averages when non-NULL
  weights are used in fitting 'model.m' and 'model.y'. This will be useful
  when data does not come from a simple random sample, for example.

  As of version 4.0, the mediator model can be of either 'lm', 'glm' (or
  `bayesglm'), 'polr' (or `bayespolr'), 'gam', 'rq', or `survreg' class
  corresponding respectively to the linear regression models,
  generalized linear models, ordered response models, generalized additive
  models, quantile regression models, or parametric duration models. For
  binary response models, the 'mediator' must be a
  numeric variable with values 0 or 1 as opposed to a factor.
  Quasi-likelihood-based inferences are not allowed for the mediator model
  because the functional form must be exactly specified for the estimation
  algorithm to work.  The 'binomial' family can only be used for binary
  response mediators and cannot be used for multiple-trial responses.  This
  is due to conflicts between how the latter type of models are implemented
  in \code{\link{glm}} and how 'mediate' is currently written.

  For the outcome model, the censored regression model fitted via package
  \code{VGAM} (of class 'vglm' with 'family@vfamily' equal to "tobit") can be
  used in addition to the models listed above for the mediator.  The
  'mediate' function is not compatible with censored regression models fitted
  via other packages.  When the quantile regression is used for the outcome
  model ('rq'), the estimated quantities are quantile causal mediation
  effects, quantile direct effects and etc., instead of the average effects.
  If the outcome model is of class 'survreg', the name of the outcome
  variable must be explicitly supplied in the `outcome' argument. This is due
  to the fact that 'survreg' objects do not contain that information in an
  easily extractable form. It should also be noted that for
  \code{\link{survreg}} models, the \code{\link{Surv}} function must be
  directly used in the model formula in the call to the survreg function, and
  that censoring types requiring more than two arguments to Surv (e.g.,
  interval censoring) are not currently supported by 'mediate'.

  The quasi-Bayesian approximation (King et al. 2000) cannot be used if
  'model.m' is of class 'rq' or 'gam', or if 'model.y' is of class 'gam',
  'polr' or 'bayespolr'. In these cases, either an error message is returned
  or use of the nonparametric bootstrap is forced. Users should note that use
  of the nonparametric bootstrap often requires significant computing time,
  especially when 'sims' is set to a large value.

  The 'control' argument must be provided when 'gam' is used for the outcome
  model and user wants to allow ACME and ADE to vary as functions of the
  treatment (i.e., to relax the "no interaction" assumption). Note that the
  outcome model must be fitted via package \code{\link{mgcv}} with
  appropriate formula using \code{\link{s}} constructs (see Imai et al. 2009
  in the references). For other model types, the interaction can be allowed
  by including an interaction term between \eqn{T} and \eqn{M} in the linear
  predictor of the outcome model. As of version 3.0, the 'INT' argument is
  deprecated and the existence of the interaction term is automatically
  detected (except for 'gam' outcome models).

  When the treatment variable is continuous or a factor with multiple levels,
  user must specify the values of \mjeqn{t_1}{t1} and \mjeqn{t_0}{t0} using the
  'treat.value' and 'control.value' arguments, respectively.  The value of
  \mjseqn{t} in the above expressions is set to \mjeqn{t_0}{t0} for 'd0', 'z0',
  etc. and to \mjeqn{t_1}{t1} for 'd1', 'z1', etc.
}
\examples{
# For illustration purposes a small number of simulations are used
# Example: Zero-inflated Count Outcome and Binary Mediator
# Generate example data
n <- 50
# Generate binary treatment variable
z <- rep(0, n)
z[sample(1:n, n/2)] <- 1
# Generate binary covariate variable
x1 <- rbinom(n, 1, p = 0.5)
# Generate continuous covariate variable
x2 <- rnorm(n)
# Create binary mediator
m <- c(1, 1, 0, 1, 1, 1, 1, 0, 1, 1,
       0, 1, 1, 1, 1, 0, 1, 0, 1, 0,
       1, 1, 0, 1, 0, 1, 1, 1, 1, 0,
       1, 1, 0, 1, 0, 1, 1, 0, 1, 1,
       0, 1, 1, 0, 0, 1, 1, 1, 0, 1)
# Create Zero-inflated Count Outcome
y <- c(0, 7, 0, 0, 8, 0, 7, 3, 5, 0,
       5, 0, 0, 1, 10,0, 3, 4, 5, 7,
       9, 7, 0, 6, 2, 4, 0, 0, 0, 6,
       0, 11,0, 5, 6, 0, 0, 12,0, 0,
       13,6, 8, 6, 5, 0, 4, 0, 6, 8)
sdata <- data.frame(z, x1, x2, m, y)
mFit <- glm(m ~ z + x1 + x2, data = sdata, family = binomial)
# Fit with Zero-inflated Poisson model
yzipFit <- zeroinfl(y ~ z + m + x1 + x2, data = sdata)
# Estimation via Quasi-Bayesian approximation
zipMA <- mediate_zi(mFit, yzipFit, sims = 100, treat = "z", mediator = "m")
summary(zipMA)
# Estimation via bootstrap approximation
zipMA_bt <- mediate_zi(mFit, yzipFit, sims = 100, boot = TRUE, treat = "z",
mediator = "m")
summary(zipMA_bt)

}
\references{
Cheng, J., Cheng, N.F., Guo, Z., Gregorich, S., Ismail, A.I.,
  Gansky, S.A (2018) Mediation analysis for count and zero-inflated count
  data. Statistical Methods in Medical Research. 27(9):2756-2774.

  Tingley, D., Yamamoto, T., Hirose, K., Imai, K. and Keele, L.
  (2014). "mediation: R package for Causal Mediation Analysis", Journal of
  Statistical Software, Vol. 59, No. 5, pp. 1-38.

  Imai, K., Keele, L., Tingley, D. and Yamamoto, T. (2011). Unpacking the
  Black Box of Causality: Learning about Causal Mechanisms from Experimental
  and Observational Studies, American Political Science Review, Vol. 105, No.
  4 (November), pp. 765-789.

  Imai, K., Keele, L. and Tingley, D. (2010) A General Approach to Causal
  Mediation Analysis, Psychological Methods, Vol. 15, No. 4 (December), pp.
  309-334.

  Imai, K., Keele, L. and Yamamoto, T. (2010) Identification, Inference, and
  Sensitivity Analysis for Causal Mediation Effects, Statistical Science,
  Vol. 25, No. 1 (February), pp. 51-71.

  Imai, K., Keele, L., Tingley, D. and Yamamoto, T. (2009) "Causal Mediation
  Analysis Using R" in Advances in Social Science Research Using R, ed. H. D.
  Vinod New York: Springer.
}
\seealso{
\code{\link{summary.mediate}}, \code{\link{plot.mediate}}
}
\author{
Nancy Cheng,
  \email{Nancy.Cheng@ucsf.edu}; Jing Cheng,
  \email{Jing.Cheng@ucsf.edu}.
}
