#' Plot Sensitivity Graph
#'
#' 'plot_sensitivity' performs sensitivity analysis when there is a treatment-
#' induced mediator-outcome confounder by varying coefficient of treatment in
#' the confounder model and plots sensitivity graph.
#'
#' @details The function uses the estimated coefficient of treatment in the
#'   fitted confounder model with observed data (see Cheng et al 2018) as the
#'   middle point (mid.beta.u), and calculates 1 unit change in the coefficient
#'   as |2*mid.beta.u/3/(n-1)| where n = n.beta.u if n.beta.u is an odd number
#'   and n = n.beta.u+1 if n.beta.u is an even number, and n.beta.u is the
#'   number of varying coefficient of treatment specified by users.
#'   \code{\link{mediate_zi_vcoef}} is then used to estimate average causal
#'   mediation effects (indirect effect), average direct effects, and total
#'   effect for each value of the treatment coefficient in the confounder
#'   model. The function prints all the effect estimates, and produces the
#'   sensitivity graph.
#' @param model.u A fitted model object for confounder. Can be of class 'lm',
#'   'polr', 'bayespolr', 'glm', 'bayesglm', 'gam', 'rq', or 'survreg'.
#' @param n.beta.u Number of varying coefficient of treatment in the confounder
#'   model. This number is adjusted as n.beta.u + 1 if n.beta.u is an even
#'   number.
#' @param model.m A fitted model object for mediator. Can be of class 'lm',
#'   'polr', 'bayespolr', 'glm', 'bayesglm', 'gam', 'rq', or 'survreg'.
#' @param model.y A fitted model object for outcome. Can be of class 'lm',
#'   'polr', 'bayespolr', 'glm', 'bayesglm', 'gam', 'vglm', 'rq', 'survreg',
#'   or 'zeroinfl'.
#' @param sims Number of Monte Carlo draws for quasi-Bayesian approximation.
#' @param boot A logical value. if 'FALSE' a quasi-Bayesian approximation is
#'   used for confidence intervals; if 'TRUE' nonparametric bootstrap will be
#'   used. Default is 'FALSE'.
#' @param confounder A character string indicating the name of the confounder
#'   variable used in the models.
#' @param treat A character string indicating the name of the treatment
#'   variable used in the models.  The treatment can be either binary (integer
#'   or a two-valued factor) or continuous (numeric).
#' @param mediator A character string indicating the name of the mediator
#'   variable used in the models.
#' @param covariates A list or data frame containing values for a subset of the
#'   pre-treatment covariates in 'model.m' and 'model.y'. If provided, the
#'   function will return the estimates conditional on those covariate values.
#'   Default is NULL.
#' @param outcome A character string indicating the name of the outcome
#'   variable in `model.y'. Only necessary if 'model.y' is of class 'survreg';
#'   otherwise ignored.Default is NULL.
#' @param digits integer indicating the number of decimal places to round
#'   the values to be returned. Default is 3.
#' @param xlab,ylab Labels for x and y axes, as in \code{\link{plot}}.
#'   Default xlab = "Beta.u", ylab = "Effect".
#' @param xlim,ylim Ranges of x and y axes, as in \code{\link{plot}}.
#' The default value, NULL, indicates that the range of the finite values to be
#' plotted should be used.
#' @param main A main title for the plot, as in \code{\link{plot}}.
#' @param type A character string (length 1 vector) or vector of 1-character
#'   strings indicating the type of plot for 3 columns of y, as in
#'   \code{\link{plot}}. Default is 'l' for lines.
#' @param col A vector of strings indicating the colors for 3 lines of y.
#'   Default is ("black","black","black").
#' @param pch A vector of plotting characters or symbols, as in
#'   \code{\link{plot}}. Default is NULL.
#' @param lty A vector of line types, as in \code{\link{plot}}. Default is
#'   c(1,2,3).
#' @param legend.x,legend.y The x and y co-ordinates to be used to position
#'   the legend, see x, y in \code{\link{legend}}. Default legend.y = NULL.
#' @param legend.inset Inset distance(s) from the margins as a fraction of the
#'   plot region when legend is placed by keyword, as \code{inset} in
#'   \code{\link{legend}}. Default is 0.05.
#' @param legend A character or expression vector of 3 to appear in the
#'   legend, as in \code{\link{legend}}. Default is c("Direct Effect",
#'   "Mediation Effect", "Total Effect")
#' @param legend.horiz A logical value. If TRUE, set the legend horizontally
#'   Default is FALSE, which sets the legend vertically.
#' @return \code{plot_sensitivity} produces sensitivity graph and returns an
#'   object of data frame with the following columns:
#'    \item{beta.u}{coefficient of treatment in the confounder model.}
#'   \item{d0, d1}{point estimates for average causal mediation effects under
#'   the control and treatment conditions.}
#'   \item{z0, z1}{point estimates for average direct effect under the control
#'   and treatment conditions.}
#'   \item{d.avg, z.avg}{simple averages of d0 and d1, z0 and z1, respectively.}
#'   \item{tau.coef}{point estimate for total effect.}
#'
#' @author Nancy Cheng,
#'   \email{Nancy.Cheng@@ucsf.edu}; Jing Cheng,
#'   \email{Jing.Cheng@@ucsf.edu}.
#'
#' @seealso \code{\link{mediate_zi_vcoef}}
#'
#' @references
#' Cheng, J., Cheng, N.F., Guo, Z., Gregorich, S., Ismail, A.I.,
#'   Gansky, S.A (2018) Mediation analysis for count and zero-inflated count
#'   data. Statistical Methods in Medical Research. 27(9):2756-2774.
#'
#' Ismail AI, Ondersma S, Willem Jedele JM, et al. (2011) Evaluation of
#'  a brief tailored motivational intervention to prevent early childhood
#'  Community Dentistry and Oral Epidemiology 39: 433–448.
#'
#' @export
#' @examples
#'
#' data("midvd_bt100")
#' uFit <- glm(PDVisit_6 ~ intervention + BrushTimes_W2 + HealthyMeals_W2
#'                         + PDVisit_W2,
#'             family = 'binomial', data = midvd_bt100)
#' mFit <- glm(PBrushBedt_6 ~ intervention + BrushTimes_W2 + HealthyMeals_W2
#'                            + PBrush_W2 + PDVisit_6,
#'             family = 'binomial', data = midvd_bt100)
#' yFit <- zeroinfl(Untreated_W3 ~ intervention + PBrushBedt_6 + BrushTimes_W2
#'                                 + HealthyMeals_W2 + PBrush_W2 + PDVisit_6,
#'                  data = midvd_bt100)
#' # For illustration purposes a small number of simulations are used
#' plot_sensitivity(uFit, n.beta.u = 5, mFit, yFit, sims = 25,
#'                  treat = "intervention", mediator = "PBrushBedt_6",
#'                  confounder = "PDVisit_6",
#'                  main = "Effects on the number of new untreated cavities at 2 years",
#'                  legend.x = "right")
#'
plot_sensitivity <- function(model.u, n.beta.u = 10, model.m, model.y,
                             sims = 1000, boot = FALSE,
                             confounder ="confd.name", treat = "treat.name",
                             mediator = "med.name", covariates = NULL,
                             outcome = NULL, digits=3,
                             xlab = "Beta.u", ylab = "Effect", xlim = NULL,
                             ylim = NULL, main = NULL,  type = "l",
                             col = c("black", "black", "black"), pch = NULL,
                             lty = c(1, 2, 3), legend.x, legend.y = NULL,
                             legend.inset = 0.05,
                             legend = c("Direct Effect", "Mediation Effect",
                                       "Total Effect"),
                             legend.horiz = FALSE) {
  n <- ifelse(n.beta.u %% 2 == 0, n.beta.u + 1, n.beta.u)

  d0 <- vector(length = n)
  d1 <- vector(length = n)
  z0 <- vector(length = n)
  z1 <- vector(length = n)
  d.avg <- vector(length = n)
  z.avg <- vector(length = n)
  tau.coef <- vector(length = n)

  beta.u <- vector(length = n)
  mid.beta.u <- coef(model.u)[treat]
  inc <- abs(2 * mid.beta.u / 3 / (n - 1))
  offset <- -inc * (n + 1) / 2

  for (i in 1:n) {
    delta.beta.u <- offset + inc * i

    ee <- mediate_zi_vcoef(model.u, delta.beta.u, model.m, model.y, sims = sims,
                           boot = boot, treat = treat, mediator = mediator,
                           confounder = confounder)

    beta.u[i] <- mid.beta.u + delta.beta.u
    d0[i] <- ee$d0
    d1[i] <- ee$d1
    z0[i] <- ee$z0
    z1[i] <- ee$z1
    d.avg[i] <- ee$d.avg
    z.avg[i] <- ee$z.avg
    tau.coef[i] <- ee$tau.coef
  }
  results <- data.frame(beta.u, d0, d1, z0, z1, d.avg, z.avg, tau.coef)

  #plot the sensitivity graph
  matplot(results$beta.u, cbind(results$z.avg, results$d.avg, results$tau.coef),
          type = type, col = col, lty = lty, xlab = xlab, ylab = ylab,
          main = main, xlim = xlim, ylim = ylim)
  #add legend
  legend(x = legend.x, y = legend.y, inset = legend.inset, legend = legend,
         col = col, lty = lty, horiz = legend.horiz)

  #round the results
  results<-round(results, digits=digits)

  return(results)
}
