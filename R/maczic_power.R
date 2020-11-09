#' Power Calculation for Mediation Analysis with Count and Zero-Inflated
#' Count Data
#'
#' 'maczic_power' computes powers to detect average causal mediation effects
#' (indirect effect), average direct effects, and total effect. This function
#' uses simulations of 3 optional covariates (binary, normal, and
#' multinomial), mediator (can be binary or continuous), and outcome (can be
#' Normal, Poisson, Negative Binomial, zero-inflated Poisson/Negative Binomial)
#' based on user-specified parameter values.
#'
#' @param nsim Number of simulations.
#' @param nsp Total sample size.
#' @param mtype Type of mediator, either 'binary' or 'continuous'.
#' @param ydist Outcome distribution. Can be 'poisson', 'negbin', 'zip', 'zinb',
#'   or 'normal'.
#' @param size Dispersion parameter for negative binomial outcome distribution.
#'   Default is 1. It's ignored for other outcome distribution.
#' @param ymax Maximum value of outcome allowed.
#' @param boot A logical value. if 'FALSE' a quasi-Bayesian approximation is
#'   used for confidence intervals; if 'TRUE' nonparametric bootstrap will be
#'   used. Default is 'FALSE'.
#' @param sims Number of Monte Carlo draws for nonparametric bootstrap or
#'   quasi-Bayesian approximation. Default is 1000.
#' @param conf.level Level of the returned two-sided confidence intervals.
#'   The default value, 0.95, is to return the 2.5 and 97.5 percentiles of the
#'   simulated quantities.
#' @param lpct0 Low bound for percent of zeros in outcome. Default is 0.
#' @param hpct0 High bound for percent of zeros in outcome. Default is 100.
#' @param px1 Probability of the binary covariate being 1.
#' @param am User-specified value for the intercept in the mediator model.
#' @param bm User-specified value for the treatment coefficient in the mediator
#'   model.
#' @param e1m User-specified value for the coefficient of the binary covariate
#'   variable in the mediator model.
#' @param e2m User-specified value for the coefficient of the continuous
#'   covariate variable in the mediator model.
#' @param e3m User-specified value for the coefficient of the multinomial
#'   covariate variable in the mediator model.
#' @param ag User-specified value for the intercept in the outcome model or in
#'   the count model of zero-inflated Poisson/Negative Binomial outcome.
#' @param bg User-specified value for the treatment coefficient in the outcome
#'   model or in the count model of zero-inflated Poisson/Negative Binomial
#'   outcome.
#' @param gg User-specified value for the mediator coefficient in the outcome
#'   model or in the count model of zero-inflated Poisson/Negative Binomial
#'   outcome.
#' @param e1g User-specified value for the coefficient of the binary covariate
#'   variable in the outcome model or in the count model of zero-inflated
#'   Poisson/Negative Binomial outcome.
#' @param e2g User-specified value for the coefficient of the continuous
#'   covariate variable in the outcome model or in the count model of
#'   zero-inflated Poisson/Negative Binomial outcome.
#' @param e3g User-specified value for the coefficient of the multinomial
#'   covariate variable in the outcome model or in the count model of
#'   zero-inflated Poisson/Negative Binomial outcome.
#' @param delta User-specified value for the treatment-by-mediator interaction
#'   coefficient in the outcome model or in the count model of zero-inflated
#'   Poisson/Negative Binomial outcome model. Default is 0.
#' @param ag2 User-specified value for the intercept in the zero-inflation
#'   model of zero-inflated Poisson/Negative Binomial outcome. Note that this
#'   argument along with the following argument bg2, gg2, e1g2, e2g2, e3g2 and
#'   delta2 only apply to zero-inflated Poisson/Negative Binomial outcome.
#'   Default is 0.
#' @param bg2 User-specified value for the treatment coefficient in the
#'   zero-inflation model of zero-inflated Poisson/Negative Binomial outcome.
#'   Default is 0.
#' @param gg2 User-specified value for the mediator coefficient in the
#'   zero-inflation model of zero-inflated Poisson/Negative Binomial outcome.
#'   Default is 0.
#' @param e1g2 User-specified value for the coefficient of the binary covariate
#'   variable in the zero-inflation model of zero-inflated Poisson/Negative
#'   Binomial outcome. Default is 0.
#' @param e2g2 User-specified value for the coefficient of the continuous
#'   covariate variable in the zero-inflation model of zero-inflated
#'   Poisson/Negative Binomial outcome. Default is 0.
#' @param e3g2 User-specified value for the coefficient of the multinomial
#'   covariate variable in the zero-inflation model of zero-inflated
#'   Poisson/Negative Binomial outcome. Default is 0.
#' @param delta2 User-specified value for the coefficient of
#'   treatment-by-mediator interaction in the zero-inflation model of
#'   zero-inflated Poisson/Negative Binomial outcome. Default is 0.
#'
#' @return \code{maczic_power} returns a data frame with the following
#'   components and prints them out in a matrix format:
#'
#'   \item{te.d0, te.d1}{average true mediation effects under the control and
#'   treatment conditions.}
#'   \item{te.z0, te.z1}{average true direct effects under the control and
#'   treatment conditions.}
#'   \item{te.tau}{average true total effect.}
#'   \item{ee.d0.rej, ee.d1.rej}{power to detect mediation effects under the
#'   control and treatment conditions.}
#'   \item{ee.z0.rej, ee.z1.rej}{power to detect direct effects under the
#'   control and treatment conditions.}
#'   \item{ee.tau.rej}{power to detect total effect.}
#'   \item{mean.y.z0}{mean outcome in control in the simulated data, not
#'   available if outcome is normal.}
#'   \item{mean.y.z1}{mean outcome in treatment in the simulated data, not
#'   available if outcome is normal.}
#'   \item{mean.y.gt0.z0}{mean non-zero outcome in control in the simulated
#'   data, not available if outcome is normal.}
#'   \item{mean.y.gt0.z1}{mean non-zero outcome in treatment in the simulated
#'   data, not available if outcome is normal}
#'   \item{pct0.y.z0}{mean percent zero outcome in control in the simulated
#'   data, not available if outcome is normal.}
#'   \item{pct0.y.z1}{mean percent zero outcome in treatment in the simulated
#'   data, not available if outcome is normal.}
#'
#' @author Nancy Cheng,
#'   \email{Nancy.Cheng@@ucsf.edu}; Jing Cheng,
#'   \email{Jing.Cheng@@ucsf.edu}.
#'
#' @seealso \code{\link{mediate_zi}}
#'
#' @references
#' Cheng, J., Cheng, N.F., Guo, Z., Gregorich, S., Ismail, A.I.,
#'   Gansky, S.A (2018) Mediation analysis for count and zero-inflated count
#'   data. Statistical Methods in Medical Research. 27(9):2756-2774.
#'
#'   Tingley, D., Yamamoto, T., Hirose, K., Imai, K. and Keele, L.
#'   (2014). "mediation: R package for Causal Mediation Analysis", Journal of
#'   Statistical Software, Vol. 59, No. 5, pp. 1-38.
#'
#'   Imai, K., Keele, L., Tingley, D. and Yamamoto, T. (2011). Unpacking the
#'   Black Box of Causality: Learning about Causal Mechanisms from Experimental
#'   and Observational Studies, American Political Science Review, Vol. 105, No.
#'   4 (November), pp. 765-789.
#'
#'   Imai, K., Keele, L. and Tingley, D. (2010) A General Approach to Causal
#'   Mediation Analysis, Psychological Methods, Vol. 15, No. 4 (December), pp.
#'   309-334.
#'
#'   Imai, K., Keele, L. and Yamamoto, T. (2010) Identification, Inference, and
#'   Sensitivity Analysis for Causal Mediation Effects, Statistical Science,
#'   Vol. 25, No. 1 (February), pp. 51-71.
#'
#'   Imai, K., Keele, L., Tingley, D. and Yamamoto, T. (2009) "Causal Mediation
#'   Analysis Using R" in Advances in Social Science Research Using R, ed. H. D.
#'   Vinod New York: Springer.
#'
#' @export
#' @examples
#' # For illustration purposes a small number of simulations are used
#' # Example 1: simulate Poisson outcome with sample size 100, binary mediator
#' # and 2 covariate (binary and normal) variables
#' posOut <- maczic_power(nsim = 8, nsp = 100,  mtype = 'binary',
#'                        sims = 40, ydist = "Poisson", ymax = 60,
#'                        px1 = 0.5, am = 0.2, bm = 0.5,
#'                        e1m = 0.1, e2m = 0.1, e3m = 0,
#'                        ag = 0.1, bg = 0.3, gg = 1,
#'                        e1g = 0.5, e2g = -0.2, e3g = 0, delta = 0.1)
#'
#' # Example 2: simulate zero-inflated Poisson outcome with sample size 80,
#' # continuous mediator and 1 normal covariate variable
#' zipOut <-maczic_power(nsim = 5, nsp = 80, mtype = 'continuous',
#'                        sims=30, ydist = "zip", ymax = 88, hpct0 = 60,
#'                        px1 = 0.5, am = 0.1, bm = 1,
#'                        e1m = 0, e2m = 0.2, e3m = 0,
#'                        ag = 0.6, bg = 0.6, gg = 0.2, e1g = 0, e2g = -0.2,
#'                        e3g = 0, ag2 = -0.7, bg2 = 0.2, gg2 = 0.1,
#'                        e2g2 = 0.1, delta2 = 0.15)
#'
maczic_power <- function(nsim, nsp, mtype, boot = FALSE, sims = 1000,
                         conf.level = 0.95, ydist, size=1, ymax, lpct0 = 0,
                         hpct0 = 100, px1, am, bm, e1m, e2m, e3m,
                         ag, bg, gg, e1g, e2g, e3g, delta = 0,
                         ag2=0, bg2=0, gg2=0, e1g2=0, e2g2=0, e3g2=0,
                         delta2 = 0) {
  sim <- vector(length = nsim)

  te.d0 <- vector(length = nsim)
  te.d1 <- vector(length = nsim)
  te.z0 <- vector(length = nsim)
  te.z1 <- vector(length = nsim)
  te.tau <- vector(length = nsim)

  ee.d0 <- vector(length = nsim)
  ee.d1 <- vector(length = nsim)
  ee.z0 <- vector(length = nsim)
  ee.z1 <- vector(length = nsim)
  ee.tau <- vector(length = nsim)

  ee.d0.lcl <- vector(length = nsim)
  ee.d0.ucl <- vector(length = nsim)
  ee.d0.ci.cover <- vector(length = nsim)

  ee.d1.lcl <- vector(length = nsim)
  ee.d1.ucl <- vector(length = nsim)
  ee.d1.ci.cover <- vector(length = nsim)

  ee.z0.lcl <- vector(length = nsim)
  ee.z0.ucl <- vector(length = nsim)
  ee.z0.ci.cover <- vector(length = nsim)

  ee.z1.lcl <- vector(length = nsim)
  ee.z1.ucl <- vector(length = nsim)
  ee.z1.ci.cover <- vector(length = nsim)

  ee.tau.lcl <- vector(length = nsim)
  ee.tau.ucl <- vector(length = nsim)
  ee.tau.ci.cover <- vector(length = nsim)

  d0.diff <- vector(length = nsim)
  d1.diff <- vector(length = nsim)
  z0.diff <- vector(length = nsim)
  z1.diff <- vector(length = nsim)
  tau.diff <- vector(length = nsim)

  ee.z1.p <- vector(length = nsim)
  ee.z0.p <- vector(length = nsim)
  ee.d1.p <- vector(length = nsim)
  ee.d0.p <- vector(length = nsim)
  ee.tau.p <- vector(length = nsim)

  ee.z1.rej <- vector(length = nsim)
  ee.z0.rej <- vector(length = nsim)
  ee.d1.rej <- vector(length = nsim)
  ee.d0.rej <- vector(length = nsim)
  ee.tau.rej <- vector(length = nsim)

  # without this def, for some unknown reason  R CMD check
  # generates a NOTE - no visible binding for global variable ‘y’
  y <- vector(length = nsim)

  pct0_y <- vector(length = nsim)
  pct0.y.z0 <- vector(length = nsim)
  pct0.y.z1 <- vector(length = nsim)
  mean.y.z0 <- vector(length = nsim)
  mean.y.z1 <- vector(length = nsim)
  mean.y.gt0.z0 <- vector(length = nsim)
  mean.y.gt0.z1 <- vector(length = nsim)

  if (tolower(mtype) != "binary" & tolower(mtype) != "continuous") {
    stop("Unsupported mediator type")
  }

  if (tolower(ydist) != "poisson" & tolower(ydist) != "negbin" &
      tolower(ydist) != "normal" & tolower(ydist) != "zip" &
      tolower(ydist) != "zinb") {
    stop("Unsupported outcome distribution")
  }

  for (i in 1:nsim) {

    writeLines(paste("Simulation:", i))
    sim[i] <- i
    error <- 1
    while (error == 1) {

      simdata1 <- .simulate1(n = nsp, mtype = mtype,
                             ydist = ydist, ymax = ymax, size = size,
                             px1 = px1, am = am, bm = bm, e1m = e1m,
                             e2m = e2m, e3m = e3m,
                             ag = ag, bg = bg, gg = gg,
                             e1g = e1g, e2g = e2g, e3g = e3g,
                             delta = delta,
                             ag2 = ag2, bg2 = bg2, gg2 = gg2,
                             e1g2 = e1g2, e2g2 = e2g2, e3g2 = e3g2,
                             delta2 = delta2)
      if (tolower(mtype) == "binary") {
        # for binary mediator Fit glm with logit link for binary mediator m
        if (e1m == 0 && e2m == 0 && e3m == 0) {
          mFit <- glm(m ~ z, family = "binomial", data = simdata1)
        } else if (e1m != 0 && e2m == 0 && e3m == 0) {
          mFit <- glm(m ~ z + x1, family = "binomial", data = simdata1)
        } else if (e1m == 0 && e2m != 0 && e3m == 0) {
          mFit <- glm(m ~ z + x2, family = "binomial", data = simdata1)
        } else if (e1m == 0 && e2m == 0 && e3m != 0) {
          mFit <- glm(m ~ z + x3, family = "binomial", data = simdata1)
        } else if (e1m != 0 && e2m != 0 && e3m == 0) {
          mFit <- glm(m ~ z + x1 + x2, family = "binomial", data = simdata1)
        } else if (e1m != 0 && e2m == 0 && e3m != 0) {
          mFit <- glm(m ~ z + x1 + x3, family = "binomial", data = simdata1)
        } else if (e1m == 0 && e2m != 0 && e3m != 0) {
          mFit <- glm(m ~ z + x2 + x3, family = "binomial", data = simdata1)
        } else {
          mFit <- glm(m ~ z + x1 + x2 + x3, family = "binomial",
                      data = simdata1)
        }
      }
      else if (tolower(mtype) == "continuous") {
        # for continuous mediator Fit glm with normal continuous mediator m
        if (e1m == 0 && e2m == 0 && e3m == 0) {
          mFit <- glm(m ~ z, family = "gaussian", data = simdata1)
        } else if (e1m != 0 && e2m == 0 && e3m == 0) {
          mFit <- glm(m ~ z + x1, family = "gaussian", data = simdata1)
        } else if (e1m == 0 && e2m != 0 && e3m == 0) {
          mFit <- glm(m ~ z + x2, family = "gaussian", data = simdata1)
        } else if (e1m == 0 && e2m == 0 && e3m != 0) {
          mFit <- glm(m ~ z + x3, family = "gaussian", data = simdata1)
        } else if (e1m != 0 && e2m != 0 && e3m == 0) {
          mFit <- glm(m ~ z + x1 + x2, family = "gaussian", data = simdata1)
        } else if (e1m != 0 && e2m == 0 && e3m != 0) {
          mFit <- glm(m ~ z + x1 + x3, family = "gaussian", data = simdata1)
        } else if (e1m == 0 && e2m != 0 && e3m != 0) {
          mFit <- glm(m ~ z + x2 + x3, family = "gaussian", data = simdata1)
        } else {
          mFit <- glm(m ~ z + x1 + x2 + x3, family = "gaussian",
                      data = simdata1)
        }
      }
      #clear the previous warnings
      assign("last.warning", NULL, envir = baseenv())
      if (tolower(ydist) == "normal") {
        # Fit REG model for outcome y with mediator m
        if (e1g == 0 && e2g == 0 && e3g == 0) {
          if (delta == 0) {
            yFit <- try(glm(y ~ z + m, family = "gaussian",
                            data = simdata1), TRUE)
          } else {
            yFit <- try(glm(y ~ z + m + z * m, family = "gaussian",
                            data = simdata1), TRUE)
          }
        } else if (e1g != 0 && e2g == 0 && e3g == 0) {
          if (delta == 0) {
            yFit <- try(glm(y ~ z + m + x1, family = "gaussian",
                            data = simdata1), TRUE)
          } else {
            yFit <- try(glm(y ~ z + m + z * m + x1, family = "gaussian",
                            data = simdata1), TRUE)
          }
        } else if (e1g == 0 && e2g != 0 && e3g == 0) {
          if (delta == 0) {
            yFit <- try(glm(y ~ z + m + x2, family = "gaussian",
                            data = simdata1), TRUE)
          } else {
            yFit <- try(glm(y ~ z + m + z * m + x2, family = "gaussian",
                            data = simdata1), TRUE)
          }
        } else if (e1g == 0 && e2g == 0 && e3g != 0) {
          if (delta == 0) {
            yFit <- try(glm(y ~ z + m + x3, family = "gaussian",
                            data = simdata1), TRUE)
          } else {
            yFit <- try(glm(y ~ z + m + z * m + x3, family = "gaussian",
                            data = simdata1), TRUE)
          }
        } else if (e1g != 0 && e2g != 0 && e3g == 0) {
          if (delta == 0) {
            yFit <- try(glm(y ~ z + m + x1 + x2, family = "gaussian",
                            data = simdata1), TRUE)
          } else {
            yFit <- try(glm(y ~ z + m + z * m + x1 + x2, family = "gaussian",
                            data = simdata1), TRUE)
          }
        } else if (e1g != 0 && e2g == 0 && e3g != 0) {
          if (delta == 0) {
            yFit <- try(glm(y ~ z + m + x1 + x3, family = "gaussian",
                            data = simdata1), TRUE)
          } else {
            yFit <- try(glm(y ~ z + m + z * m + x1 + x3, family = "gaussian",
                            data = simdata1), TRUE)
          }
        } else if (e1g == 0 && e2g != 0 && e3g != 0) {
          if (delta == 0) {
            yFit <- try(glm(y ~ z + m + x2 + x3, family = "gaussian",
                            data = simdata1), TRUE)
          } else {
            yFit <- try(glm(y ~ z + m + z * m + x2 + x3, family = "gaussian",
                            data = simdata1), TRUE)
          }
        } else {
          if (delta == 0) {
            yFit <- try(glm(y ~ z + m + x1 + x2 + x3, family = "gaussian",
                            data = simdata1), TRUE)
          } else {
            yFit <- try(glm(y ~ z + m + z * m + x1 + x2 + x3,
                            family = "gaussian", data = simdata1), TRUE)
          }
        }
      }
      else if (tolower(ydist) == "poisson") {
        # Fit poisson regression model for outcome y with mediator m
        if (e1g == 0 && e2g == 0 && e3g == 0) {
          if (delta == 0) {
            yFit <- try(glm(y ~ z + m, family = "poisson",
                            data = simdata1), TRUE)
          } else {
            yFit <- try(glm(y ~ z + m + z * m, family = "poisson",
                            data = simdata1), TRUE)
          }
        } else if (e1g != 0 && e2g == 0 && e3g == 0) {
          if (delta == 0) {
            yFit <- try(glm(y ~ z + m + x1, family = "poisson",
                            data = simdata1), TRUE)
          } else {
            yFit <- try(glm(y ~ z + m + z * m + x1, family = "poisson",
                            data = simdata1), TRUE)
          }
        } else if (e1g == 0 && e2g != 0 && e3g == 0) {
          if (delta == 0) {
            yFit <- try(glm(y ~ z + m + x2, family = "poisson",
                            data = simdata1), TRUE)
          } else {
            yFit <- try(glm(y ~ z + m + z * m + x2, family = "poisson",
                            data = simdata1), TRUE)
          }
        } else if (e1g == 0 && e2g == 0 && e3g != 0) {
          if (delta == 0) {
            yFit <- try(glm(y ~ z + m + x3, family = "poisson",
                            data = simdata1), TRUE)
          } else {
            yFit <- try(glm(y ~ z + m + z * m + x3, family = "poisson",
                            data = simdata1), TRUE)
          }
        } else if (e1g != 0 && e2g != 0 && e3g == 0) {
          if (delta == 0) {
            yFit <- try(glm(y ~ z + m + x1 + x2, family = "poisson",
                            data = simdata1), TRUE)
          } else {
            yFit <- try(glm(y ~ z + m + z * m + x1 + x2,
                            family = "poisson", data = simdata1), TRUE)
          }
        } else if (e1g != 0 && e2g == 0 && e3g != 0) {
          if (delta == 0) {
            yFit <- try(glm(y ~ z + m + x1 + x3, family = "poisson",
                            data = simdata1), TRUE)
          } else {
            yFit <- try(glm(y ~ z + m + z * m + x1 + x3, family = "poisson",
                            data = simdata1), TRUE)
          }
        } else if (e1g == 0 && e2g != 0 && e3g != 0) {
          if (delta == 0) {
            yFit <- try(glm(y ~ z + m + x2 + x3, family = "poisson",
                            data = simdata1), TRUE)
          } else {
            yFit <- try(glm(y ~ z + m + z * m + x2 + x3, family = "poisson",
                            data = simdata1), TRUE)
          }
        } else {
          if (delta == 0) {
            yFit <- try(glm(y ~ z + m + x1 + x2 + x3, family = "poisson",
                            data = simdata1), TRUE)
          } else {
            yFit <- try(glm(y ~ z + m + z * m + x1 + x2 + x3,
                            family = "poisson", data = simdata1), TRUE)
          }
        }
      }
      else if (tolower(ydist) == "negbin") {
        # Fit negative binomial regression model for outcome y with mediator m
        if (e1g == 0 && e2g == 0 && e3g == 0) {
          if (delta == 0) {
            yFit <- try(glm.nb(y ~ z + m, data = simdata1), TRUE)
          } else {
            yFit <- try(glm.nb(y ~ z + m + z * m, data = simdata1), TRUE)
          }
        } else if (e1g != 0 && e2g == 0 && e3g == 0) {
          if (delta == 0) {
            yFit <- try(glm.nb(y ~ z + m + x1, data = simdata1), TRUE)
          } else {
            yFit <- try(glm.nb(y ~ z + m + z * m + x1, data = simdata1), TRUE)
          }
        } else if (e1g == 0 && e2g != 0 && e3g == 0) {
          if (delta == 0) {
            yFit <- try(glm.nb(y ~ z + m + x2, data = simdata1), TRUE)
          } else {
            yFit <- try(glm.nb(y ~ z + m + z * m + x2, data = simdata1), TRUE)
          }
        } else if (e1g == 0 && e2g == 0 && e3g != 0) {
          if (delta == 0) {
            yFit <- try(glm.nb(y ~ z + m + x3, data = simdata1), TRUE)
          } else {
            yFit <- try(glm.nb(y ~ z + m + z * m + x3, data = simdata1),
                        TRUE)
          }
        } else if (e1g != 0 && e2g != 0 && e3g == 0) {
          if (delta == 0) {
            yFit <- try(glm.nb(y ~ z + m + x1 + x2, data = simdata1), TRUE)
          } else {
            yFit <- try(glm.nb(y ~ z + m + z * m + x1 + x2,
                               data = simdata1), TRUE)
          }
        } else if (e1g != 0 && e2g == 0 && e3g != 0) {
          if (delta == 0) {
            yFit <- try(glm.nb(y ~ z + m + x1 + x3, data = simdata1), TRUE)
          } else {
            yFit <- try(glm.nb(y ~ z + m + z * m + x1 + x3, data = simdata1),
                        TRUE)
          }
        } else if (e1g == 0 && e2g != 0 && e3g != 0) {
          if (delta == 0) {
            yFit <- try(glm.nb(y ~ z + m + x2 + x3, data = simdata1), TRUE)
          } else {
            yFit <- try(glm.nb(y ~ z + m + z * m + x2 + x3, data = simdata1),
                        TRUE)
          }
        } else {
          if (delta == 0) {
            yFit <- try(glm.nb(y ~ z + m + x1 + x2 + x3, data = simdata1),
                        TRUE)
          } else {
            yFit <- try(glm.nb(y ~ z + m + z * m + x1 + x2 + x3,
                               data = simdata1), TRUE)
          }
        }
      }
      else if (tolower(ydist) == "zip") {
        # Fit ZIP model for outcome y
        if (e1g == 0 && e2g == 0 && e3g == 0 && e1g2 == 0 && e2g2 == 0
            && e3g2 == 0) {
          if (delta == 0 && delta2 == 0) {
            yFit <- try(zeroinfl(y ~ z + m, data = simdata1), TRUE)
          } else {
            yFit <- try(zeroinfl(y ~ z + m + z * m, data = simdata1), TRUE)
          }

        } else if ((e1g != 0 || e1g2 != 0) && (e2g == 0 && e2g2 == 0) &&
                   (e3g == 0 && e3g2 == 0)) {
          if (delta == 0 && delta2 == 0) {
            yFit <- try(zeroinfl(y ~ z + m + x1, data = simdata1), TRUE)
          } else {
            yFit <- try(zeroinfl(y ~ z + m + z * m + x1, data = simdata1),
                        TRUE)
          }
        } else if ((e1g == 0 && e1g2 == 0) && (e2g != 0 || e2g2 != 0) &&
                   (e3g == 0 && e3g2 == 0)) {
          if (delta == 0 && delta2 == 0) {
            yFit <- try(zeroinfl(y ~ z + m + x2, data = simdata1), TRUE)
          } else {
            yFit <- try(zeroinfl(y ~ z + m + z * m + x2, data = simdata1),
                        TRUE)
          }
        } else if ((e1g == 0 && e1g2 == 0) && (e2g == 0 && e2g2 == 0) &&
                   (e3g != 0 || e3g2 != 0)) {
          if (delta == 0 && delta2 == 0) {
            yFit <- try(zeroinfl(y ~ z + m + x3, data = simdata1), TRUE)
          } else {
            yFit <- try(zeroinfl(y ~ z + m + z * m + x3, data = simdata1),
                        TRUE)
          }
        } else if ((e1g != 0 || e1g2 != 0) && (e2g != 0 || e2g2 != 0) &&
                   (e3g == 0 && e3g2 == 0)) {
          if (delta == 0 && delta2 == 0) {
            yFit <- try(zeroinfl(y ~ z + m + x1 + x2, data = simdata1), TRUE)
          } else {
            yFit <- try(zeroinfl(y ~ z + m + z * m + x1 + x2,
                                 data = simdata1), TRUE)
          }
        } else if ((e1g != 0 || e1g2 != 0) && (e2g == 0 && e2g2 == 0) &&
                   (e3g != 0 || e3g2 != 0)) {
          if (delta == 0 && delta2 == 0) {
            yFit <- try(zeroinfl(y ~ z + m + x1 + x3, data = simdata1), TRUE)
          } else {
            yFit <- try(zeroinfl(y ~ z + m + z * m + x1 + x3,
                                 data = simdata1), TRUE)
          }
        } else if ((e1g == 0 && e1g2 == 0) && (e2g != 0 || e2g2 != 0) &&
                   (e3g != 0 || e3g2 != 0)) {
          if (delta == 0 && delta2 == 0) {
            yFit <- try(zeroinfl(y ~ z + m + x2 + x3, data = simdata1), TRUE)
          } else {
            yFit <- try(zeroinfl(y ~ z + m + z * m + x2 + x3,
                                 data = simdata1), TRUE)
          }
        } else {
          if (delta == 0  && delta2 == 0) {
            yFit <- try(zeroinfl(y ~ z + m + x1 + x2 + x3, data = simdata1),
                        TRUE)
          } else {
            yFit <- try(zeroinfl(y ~ z + m + z * m + x1 + x2 + x3,
                                 data = simdata1), TRUE)
          }
        }
      }
      else {
        # Fit ZINB model for outcome y
        if (e1g == 0 && e2g == 0 && e3g == 0 && e1g2 == 0 && e2g2 == 0 &&
            e3g2 == 0) {
          if (delta == 0 && delta2 == 0) {
            yFit <- try(zeroinfl(y ~ z + m, data = simdata1,
                                 dist = "negbin"), TRUE)
          } else {
            yFit <- try(zeroinfl(y ~ z + m + z * m, data = simdata1,
                                 dist = "negbin"), TRUE)
          }

        } else if ((e1g == 0 && e1g2 == 0) && (e2g != 0 || e2g2 != 0)
                   && (e3g == 0 && e3g2 == 0)) {
          if (delta == 0 && delta2 == 0) {
            yFit <- try(zeroinfl(y ~ z + m + x2, data = simdata1,
                                 dist = "negbin"), TRUE)
          } else {
            yFit <- try(zeroinfl(y ~ z + m + z * m + x2, data = simdata1,
                                 dist = "negbin"), TRUE)
          }
        } else if ((e1g == 0 && e1g2 == 0) && (e2g == 0 && e2g2 == 0) &&
                   (e3g != 0 || e3g2 != 0)) {
          if (delta == 0 && delta2 == 0) {
            yFit <- try(zeroinfl(y ~ z + m + x3, data = simdata1,
                                 dist = "negbin"), TRUE)
          } else {
            yFit <- try(zeroinfl(y ~ z + m + z * m + x3, data = simdata1,
                                 dist = "negbin"), TRUE)
          }
        } else if ((e1g != 0 || e1g2 != 0) && (e2g != 0 || e2g2 != 0) &&
                   (e3g == 0 && e3g2 == 0)) {
          if (delta == 0 && delta2 == 0) {
            yFit <- try(zeroinfl(y ~ z + m + x1 + x2, data = simdata1,
                                 dist = "negbin"), TRUE)
          } else {
            yFit <- try(zeroinfl(y ~ z + m + z * m + x1 + x2, data = simdata1,
                                 dist = "negbin"), TRUE)
          }
        } else if ((e1g != 0 || e1g2 != 0) && (e2g == 0 && e2g2 == 0) &&
                   (e3g != 0 || e3g2 != 0)) {
          if (delta == 0 && delta2 == 0) {
            yFit <- try(zeroinfl(y ~ z + m + x1 + x3, data = simdata1,
                                 dist = "negbin"), TRUE)
          } else {
            yFit <- try(zeroinfl(y ~ z + m + z * m + x1 + x3,
                                 data = simdata1, dist = "negbin"), TRUE)
          }
        } else if ((e1g == 0 && e1g2 == 0) && (e2g != 0 || e2g2 != 0) &&
                   (e3g != 0 || e3g2 != 0)) {
          if (delta == 0 && delta2 == 0) {
            yFit <- try(zeroinfl(y ~ z + m + x2 + x3, data = simdata1,
                                 dist = "negbin"), TRUE)
          } else {
            yFit <- try(zeroinfl(y ~ z + m + z * m + x2 + x3,
                                 data = simdata1, dist = "negbin"), TRUE)
          }

        } else {
          if (delta == 0 && delta2 == 0) {
            yFit <- try(zeroinfl(y ~ z + m + x1 + x2 + x3, data = simdata1,
                                 dist = "negbin"), TRUE)
          } else {
            yFit <- try(zeroinfl(y ~ z + m + z * m + x1 + x2 + x3,
                                 data = simdata1, dist = "negbin"), TRUE)
          }
        }
      }


      len <- length(grep("Error", yFit[1]))
      if (len == 0)
        error <- 0
      else error <- 1

      if (inherits(yFit, "try-error"))
        error <- 1
      else {
        if (is.na(det(vcov(yFit, model = "zero"))) ||
            is.na(det(vcov(yFit, model = "count"))))
          error <- 1
        else if (det(vcov(yFit, model = "zero")) <= 0 ||
                 det(vcov(yFit, model = "count")) <= 0)
          error <- 1
        if (length(warnings()) > 0)
          error <- 1

      }

      if (tolower(ydist) == "poisson" || tolower(ydist) == "negbin"
          || tolower(ydist) == "zip" || tolower(ydist) == "zinb") {
        pct0 <- table(simdata1$y)[1] / nsp * 100

        if (pct0 < lpct0 || pct0 > hpct0)
          error <- 1
      }
    }
    if (tolower(ydist) == "poisson" || tolower(ydist) == "negbin"
        || tolower(ydist) == "zip" || tolower(ydist) == "zinb") {
      pct0_y[i] <- pct0
      pct0.y.z0[i] <- sum(simdata1$y == 0 & simdata1$z == 0) /
        sum(simdata1$z == 0) * 100
      pct0.y.z1[i] <- sum(simdata1$y == 0 & simdata1$z == 1) /
        sum(simdata1$z == 1) * 100

      mean_y.by.z <- by(simdata1$y, simdata1$z, mean)
      mean.y.z0[i] <- as.numeric(mean_y.by.z[1])
      mean.y.z1[i] <- as.numeric(mean_y.by.z[2])

      simdata1.y.gt0 <- subset(simdata1, y > 0)
      mean_y.gt0.by.z <- by(simdata1.y.gt0$y, simdata1.y.gt0$z, mean)
      mean.y.gt0.z0[i] <- as.numeric(mean_y.gt0.by.z[1])
      mean.y.gt0.z1[i] <- as.numeric(mean_y.gt0.by.z[2])
    }

    # get true effects
    te <- .compute_te(ds = simdata1)

    te.d0[i] <- te$d0
    te.d1[i] <- te$d1
    te.z0[i] <- te$z0
    te.z1[i] <- te$z1
    te.tau[i] <- te$tau.coef

    # get estimated effects
    if (boot) {
      ee <- mediate_zi(mFit, yFit, sims = sims, boot = TRUE, treat = "z",
                       mediator = "m", conf.level = conf.level)
    } else {
      # bayesian
      ee <- mediate_zi(mFit, yFit, sims = sims, treat = "z",
                       mediator = "m", conf.level = conf.level)
    }
    ee.d0[i] <- ee$d0
    ee.d1[i] <- ee$d1
    ee.z0[i] <- ee$z0
    ee.z1[i] <- ee$z1
    ee.tau[i] <- ee$tau.coef

    ee.d0.lcl[i] <- ee$d0.ci[1]
    ee.d0.ucl[i] <- ee$d0.ci[2]
    ee.d0.ci.cover[i] <- ifelse(ee.d0.lcl[i] <= te.d0[i] &
                                  te.d0[i] <= ee.d0.ucl[i], 1, 0)

    ee.d1.lcl[i] <- ee$d1.ci[1]
    ee.d1.ucl[i] <- ee$d1.ci[2]
    ee.d1.ci.cover[i] <- ifelse(ee.d1.lcl[i] <= te.d1[i] &
                                  te.d1[i] <= ee.d1.ucl[i], 1, 0)

    ee.z0.lcl[i] <- ee$z0.ci[1]
    ee.z0.ucl[i] <- ee$z0.ci[2]
    ee.z0.ci.cover[i] <- ifelse(ee.z0.lcl[i] <= te.z0[i] &
                                  te.z0[i] <= ee.z0.ucl[i], 1, 0)

    ee.z1.lcl[i] <- ee$z1.ci[1]
    ee.z1.ucl[i] <- ee$z1.ci[2]
    ee.z1.ci.cover[i] <- ifelse(ee.z1.lcl[i] <= te.z1[i] &
                                  te.z1[i] <= ee.z1.ucl[i], 1, 0)

    ee.tau.lcl[i] <- ee$tau.ci[1]
    ee.tau.ucl[i] <- ee$tau.ci[2]
    ee.tau.ci.cover[i] <- ifelse(ee.tau.lcl[i] <= te.tau[i] &
                                   te.tau[i] <= ee.tau.ucl[i], 1, 0)

    d0.diff[i] <- ee.d0[i] - te.d0[i]
    d1.diff[i] <- ee.d1[i] - te.d1[i]
    z0.diff[i] <- ee.z0[i] - te.z0[i]
    z1.diff[i] <- ee.z1[i] - te.z1[i]
    tau.diff[i] <- ee.tau[i] - te.tau[i]

    ee.z1.p[i] <- ee$z1.p
    ee.z0.p[i] <- ee$z0.p
    ee.d1.p[i] <- ee$d1.p
    ee.d0.p[i] <- ee$d0.p
    ee.tau.p[i] <- ee$tau.p

    ee.z1.rej[i] <- ifelse(ee.z1.p[i] < 0.05, 1, 0)
    ee.z0.rej[i] <- ifelse(ee.z0.p[i] < 0.05, 1, 0)
    ee.d1.rej[i] <- ifelse(ee.d1.p[i] < 0.05, 1, 0)
    ee.d0.rej[i] <- ifelse(ee.d0.p[i] < 0.05, 1, 0)
    ee.tau.rej[i] <- ifelse(ee.tau.p[i] < 0.05, 1, 0)

  }

  if (tolower(ydist) == "poisson" || tolower(ydist) == "negbin"
      || tolower(ydist) == "zip" || tolower(ydist) == "zinb") {
    simout_long <- data.frame(sim, te.d0, te.d1, te.z0, te.z1, ee.d0, ee.d1,
                              ee.z0, ee.z1, ee.tau, ee.d0.lcl, ee.d0.ucl,
                              ee.d0.ci.cover, ee.d1.lcl, ee.d1.ucl,
                              ee.d1.ci.cover, ee.z0.lcl, ee.z0.ucl,
                              ee.z0.ci.cover, ee.z1.lcl, ee.z1.ucl,
                              ee.z1.ci.cover, ee.tau.lcl, ee.tau.ucl,
                              ee.tau.ci.cover, d0.diff, d1.diff, z0.diff,
                              z1.diff, tau.diff, ee.z1.p, ee.z0.p,
                              ee.d1.p, ee.d0.p, ee.tau.p, ee.z1.rej,
                              ee.z0.rej, ee.d1.rej, ee.d0.rej, ee.tau.rej,
                              mean.y.z0, mean.y.z1, mean.y.gt0.z0,
                              mean.y.gt0.z1, pct0.y.z0, pct0.y.z1)
    simout <- data.frame(te.d0, te.d1, te.z0, te.z1, te.tau, ee.d0.rej,
                         ee.d1.rej, ee.z0.rej, ee.z1.rej, ee.tau.rej,
                         mean.y.z0, mean.y.z1, mean.y.gt0.z0, mean.y.gt0.z1,
                         pct0.y.z0, pct0.y.z1)
    # compute mean TE, rejection rate (power), and mean pct0_y
    mout <- colMeans(simout)
    out <- mout

    names(mout) <- c("True mediation effect in control",
                     "True mediation effect in treatment",
                     "True direct effect in control",
                     "True direct effect in treatment",
                     "True total effect",
                     "Power to detect mediation effect in control",
                     "Power to detect mediation effect in treatment",
                     "Power to detect direct effect in control",
                     "Power to detect direct effect in treatment",
                     "Power to detect total effect",
                     "Mean outcome in control",
                     "Mean outcome in treatment",
                     "Mean non-zero outcome in control",
                     "Mean non-zero outcome in treatment",
                     "Mean percent zero outcome in control",
                     "Mean percent zero outcome in treatment")

    #output as one-column matrix
    mout_matrix <- as.matrix(mout)

  } else { # for normal outcome

    simout_long <- data.frame(sim, te.d0, te.d1, te.z0, te.z1, te.tau, ee.d0,
                              ee.d1, ee.z0, ee.z1, ee.tau, ee.d0.lcl,
                              ee.d0.ucl, ee.d0.ci.cover, ee.d1.lcl, ee.d1.ucl,
                              ee.d1.ci.cover, ee.z0.lcl, ee.z0.ucl,
                              ee.z0.ci.cover, ee.z1.lcl, ee.z1.ucl,
                              ee.z1.ci.cover, ee.tau.lcl, ee.tau.ucl,
                              ee.tau.ci.cover, d0.diff, d1.diff, z0.diff,
                              z1.diff, tau.diff, ee.z1.p, ee.z0.p, ee.d1.p,
                              ee.d0.p, ee.tau.p, ee.z1.rej, ee.z0.rej,
                              ee.d1.rej, ee.d0.rej, ee.tau.rej)
    simout <- data.frame(te.d0, te.d1, te.z0, te.z1, te.tau, ee.d0.rej,
                         ee.d1.rej, ee.z0.rej, ee.z1.rej, ee.tau.rej)

    # compute mean TE, rejection rate (power)
    mout <- colMeans(simout)
    out <- mout
    names(mout) <- c("True mediation effect in control",
                     "True mediation effect in treatment",
                     "True direct effect in control",
                     "True direct effect in treatment",
                     "True total effect",
                     "Power to detect mediation effect in control",
                     "Power to detect mediation effect in treatment",
                     "Power to detect direct effect in control",
                     "Power to detect direct effect in treatment",
                     "Power to detect total effect")

    #output as one-column matrix
    mout_matrix <- as.matrix(mout)
  }
  writeLines("\nResults:")
  prmatrix(mout_matrix, collab = rep("", 1))
  return(out)
}

#' function to get 1 simulation data - which includes z, m0, m1, x1, x2, x3,
#' y0.m0, y0.m1, y1.m0, y1.m1, delta1, delta0, zeta1, zeta0, m, y, this
#' function is for internal use, called by maczic_power(). all the parameters
#' below are passed from the parameters of maczic_power().
#' @param n Total sample size.
#' @param mtype Type of mediator, either 'binary' or 'continuous'.
#' @param ydist Outcome distribution. Can be 'poisson', 'negbin', 'zip', 'zinb',
#'   or 'normal'.
#' @param size Dispersion parameter for negative binomial outcome distribution.
#'   Default is 1. It's ignored for other outcome distribution.
#' @param ymax Maximum value of outcome allowed.
#' @param px1 Probability of the binary covariate being 1.
#' @param am User-specified value for the intercept in the mediator model.
#' @param bm User-specified value for the treatment coefficient in the mediator
#'   model.
#' @param e1m User-specified value for the coefficient of the binary covariate
#'   variable in the mediator model.
#' @param e2m User-specified value for the coefficient of the continuous
#'   covariate variable in the mediator model.
#' @param e3m User-specified value for the coefficient of the multinomial
#'   covariate variable in the mediator model.
#' @param ag User-specified value for the intercept in the outcome model or in
#'   the count model of zero-inflated Poisson/Negative Binomial outcome.
#' @param bg User-specified value for the treatment coefficient in the outcome
#'   model or in the count model of zero-inflated Poisson/Negative Binomial
#'   outcome.
#' @param gg User-specified value for the mediator coefficient in the outcome
#'   model or in the count model of zero-inflated Poisson/Negative Binomial
#'   outcome.
#' @param e1g User-specified value for the coefficient of the binary covariate
#'   variable in the outcome model or in the count model of zero-inflated
#'   Poisson/Negative Binomial outcome.
#' @param e2g User-specified value for the coefficient of the continuous
#'   covariate variable in the outcome model or in the count model of
#'   zero-inflated Poisson/Negative Binomial outcome.
#' @param e3g User-specified value for the coefficient of the multinomial
#'   covariate variable in the outcome model or in the count model of
#'   zero-inflated Poisson/Negative Binomial outcome.
#' @param delta User-specified value for the treatment-by-mediator interaction
#'   coefficient in the outcome model or in the count model of zero-inflated
#'   Poisson/Negative Binomial outcome model. Default is 0.
#' @param ag2 User-specified value for the intercept in the zero-inflation
#'   model of zero-inflated Poisson/Negative Binomial outcome. Note that this
#'   argument along with the following argument bg2, gg2, e1g2, e2g2, e3g2 and
#'   delta2 only apply to zero-inflated Poisson/Negative Binomial outcome.
#'   Default is 0.
#' @param bg2 User-specified value for the treatment coefficient in the
#'   zero-inflation model of zero-inflated Poisson/Negative Binomial outcome.
#'   Default is 0.
#' @param gg2 User-specified value for the mediator coefficient in the
#'   zero-inflation model of zero-inflated Poisson/Negative Binomial outcome.
#'   Default is 0.
#' @param e1g2 User-specified value for the coefficient of the binary covariate
#'   variable in the zero-inflation model of zero-inflated Poisson/Negative
#'   Binomial outcome. Default is 0.
#' @param e2g2 User-specified value for the coefficient of the continuous
#'   covariate variable in the zero-inflation model of zero-inflated
#'   Poisson/Negative Binomial outcome. Default is 0.
#' @param e3g2 User-specified value for the coefficient of the multinomial
#'   covariate variable in the zero-inflation model of zero-inflated
#'   Poisson/Negative Binomial outcome. Default is 0.
#' @param delta2 User-specified value for the coefficient of
#'   treatment-by-mediator interaction in the zero-inflation model of
#'   zero-inflated Poisson/Negative Binomial outcome. Default is 0.
#'
#' @keywords internal
#' @noRd
.simulate1 <- function(n, mtype, ydist, size, ymax,
                       px1, am, bm, e1m, e2m, e3m,
                       ag, bg, gg, e1g, e2g, e3g, delta,
                       ag2, bg2, gg2, e1g2, e2g2, e3g2, delta2) {
  # randomly assign half subjects to z=1 and the other half to z=0
  z <- rep(0, n)
  z[sample(1:n, n / 2)] <- 1
  # Generate binary covariate variable x1
  x1 <- rbinom(n, 1, prob = px1)
  # Generate continous normal covariate variable x2
  x2 <- rnorm(n, 0, 1)
  # Generate multinomial covariate variable x3 from category 1, 2, 3 and 4
  # with equal
  # selection prob
  x3 <- sample(c(1, 2, 3, 4), n, replace = T,
               prob = c(0.25, 0.25, 0.25, 0.25))

  intcpt <- rep(1, n)
  if (tolower(mtype) == "binary") {
    # for binary mediator link function -logit for mediator when z=0
    h_pm0 <- am + e1m * x1 + e2m * x2 + e3m * x3
    # prob of binary mediator
    pm0 <- exp(h_pm0) / (1 + exp(h_pm0))

    ## Generate binary mediator variable m0 with prob=pm0
    m0 <- rbinom(n, 1, prob = pm0)

    # link function -logit for mediator when z=1
    h_pm1 <- am + bm + e1m * x1 + e2m * x2 + e3m * x3
    # prob of binary mediator
    pm1 <- exp(h_pm1) / (1 + exp(h_pm1))

    ## Generate binary mediator variable m1 with prob=pm1
    m1 <- rbinom(n, 1, prob = pm1)
  } else {
    # for continuous normal mediator Generate normal continuous mediator
    # variable m when z=0
    mu <- am + e1m * x1 + e2m * x2 + e3m * x3
    m0 <- rnorm(n, mu, 1)
    # Generate normal continuous mediator variable m when z=1
    m1 <- m0 + bm
  }
  ID <- 1:n
  y0.m0 <- vector(length = n)
  y0.m1 <- vector(length = n)
  y1.m0 <- vector(length = n)
  y1.m1 <- vector(length = n)
  m <- vector(length = n)
  y <- vector(length = n)

  for (k in 1:n) {

    y_gt_ymax <- 1
    while (y_gt_ymax == 1) {

      if (tolower(ydist) == "normal") {
        y0.m0[k] <- .gen_norm(1, x = cbind(intcpt, m0, x1, x2, x3),
                              b = c(ag, gg, e1g, e2g, e3g))
        y0.m1[k] <- .gen_norm(1, x = cbind(intcpt, m1, x1, x2, x3),
                              b = c(ag, gg, e1g, e2g, e3g))
        y1.m0[k] <- .gen_norm(1, x = cbind(intcpt, 1, m0, x1, x2, x3),
                              b = c(ag, bg, gg + delta, e1g, e2g, e3g))
        y1.m1[k] <- .gen_norm(1, x = cbind(intcpt, 1, m1, x1, x2, x3),
                              b = c(ag, bg, gg + delta, e1g, e2g, e3g))

      } else if (tolower(ydist) == "poisson") {
        y0.m0[k] <- .gen_poi(1, x = cbind(intcpt, m0, x1, x2, x3),
                             b = c(ag, gg, e1g, e2g, e3g))
        y0.m1[k] <- .gen_poi(1, x = cbind(intcpt, m1, x1, x2, x3),
                             b = c(ag, gg, e1g, e2g, e3g))
        y1.m0[k] <- .gen_poi(1, x = cbind(intcpt, 1, m0, x1, x2, x3),
                             b = c(ag, bg, gg + delta, e1g, e2g, e3g))
        y1.m1[k] <- .gen_poi(1, x = cbind(intcpt, 1, m1, x1, x2, x3),
                             b = c(ag, bg, gg + delta, e1g, e2g, e3g))
      } else if (tolower(ydist) == "negbin") {
        y0.m0[k] <- .gen_nb(1, x = cbind(intcpt, m0, x1, x2, x3),
                            b = c(ag, gg, e1g, e2g, e3g), size)
        y0.m1[k] <- .gen_nb(1, x = cbind(intcpt, m1, x1, x2, x3),
                            b = c(ag, gg, e1g, e2g, e3g), size)
        y1.m0[k] <- .gen_nb(1, x = cbind(intcpt, 1, m0, x1, x2, x3),
                            b = c(ag, bg, gg + delta, e1g, e2g, e3g), size)
        y1.m1[k] <- .gen_nb(1, x = cbind(intcpt, 1, m1, x1, x2, x3),
                            b = c(ag, bg, gg + delta, e1g, e2g, e3g), size)

      } else if (tolower(ydist) == "zip") {
        y0.m0[k] <- .gen_zip(1, zx = cbind(intcpt, m0, x1, x2, x3),
                             bz = c(ag2, gg2, e1g2, e2g2, e3g2),
                             cx = cbind(intcpt, m0, x1, x2, x3),
                             bc = c(ag, gg, e1g, e2g, e3g))
        y0.m1[k] <- .gen_zip(1, zx = cbind(intcpt, m1, x1, x2, x3),
                             bz = c(ag2, gg2, e1g2, e2g2, e3g2),
                             cx = cbind(intcpt, m1, x1, x2, x3),
                             bc = c(ag, gg, e1g, e2g, e3g))
        y1.m0[k] <- .gen_zip(1, zx = cbind(intcpt, 1, m0, x1, x2, x3),
                             bz = c(ag2, bg2, gg2 + delta2, e1g2, e2g2, e3g2),
                             cx = cbind(intcpt, 1, m0, x1, x2, x3),
                             bc = c(ag, bg, gg + delta, e1g, e2g, e3g))
        y1.m1[k] <- .gen_zip(1, zx = cbind(intcpt, 1, m1, x1, x2, x3),
                             bz = c(ag2, bg2, gg2 + delta2, e1g2, e2g2, e3g2),
                             cx = cbind(intcpt, 1, m1, x1, x2, x3),
                             bc = c(ag, bg, gg + delta, e1g, e2g, e3g))

      } else if (tolower(ydist) == "zinb") {
        y0.m0[k] <- .gen_zinb(1, zx = cbind(intcpt, m0, x1, x2, x3),
                              bz = c(ag2, gg2, e1g2, e2g2, e3g2),
                              cx = cbind(intcpt, m0, x1, x2, x3),
                              bc = c(ag, gg, e1g, e2g, e3g), size)
        y0.m1[k] <- .gen_zinb(1, zx = cbind(intcpt, m1, x1, x2, x3),
                              bz = c(ag2, gg2, e1g2, e2g2, e3g2),
                              cx = cbind(intcpt, m1, x1, x2, x3),
                              bc = c(ag, gg, e1g, e2g, e3g), size)
        y1.m0[k] <- .gen_zinb(1, zx = cbind(intcpt, 1, m0, x1, x2, x3),
                              bz = c(ag2, bg2, gg2 + delta2, e1g2, e2g2, e3g2),
                              cx = cbind(intcpt, 1, m0, x1, x2, x3),
                              bc = c(ag, bg, gg + delta, e1g, e2g, e3g), size)
        y1.m1[k] <- .gen_zinb(1, zx = cbind(intcpt, 1, m1, x1, x2, x3),
                              bz = c(ag2, bg2, gg2 + delta2, e1g2, e2g2, e3g2),
                              cx = cbind(intcpt, 1, m1, x1, x2, x3),
                              bc = c(ag, bg, gg + delta, e1g, e2g, e3g), size)

      } else {
        stop("Unsupported outcome distribution")
      }

      # get observed y
      if (z[k] == 1) {
        m[k] <- m1[k]
        y[k] <- y1.m1[k]
      } else {
        m[k] <- m0[k]
        y[k] <- y0.m0[k]
      }
      if (y[k] <= ymax) {
        y_gt_ymax <- 0
      }
    }
  }
  delta1 <- y1.m1 - y1.m0
  delta0 <- y0.m1 - y0.m0
  zeta1 <- y1.m1 - y0.m1
  zeta0 <- y1.m0 - y0.m0

  simdata <- data.frame(ID, z, m0, m1, x1, x2, x3, y0.m0, y0.m1, y1.m0,
                        y1.m1, delta1, delta0, zeta1, zeta0, m, y)
  simdata
}


#' Generate data y from normal distribution
#'
#' '.gen_norm' is used to generate data y from normal distribution
#'
#' @param n sample size
#' @param x covariates
#' @param b coefficients of corresponding covariates
#' @keywords internal
#' @noRd
.gen_norm <- function(n, x, b) {
  u <- b %*% t(x)
  # Generate outcome y
  y <- vector(length = n)
  y <- rnorm(n, u, 1)
  y
}

#' Generate data y from Poisson distribution
#'
#' '.gen_poi' is used to generate data y from Poisson distribution
#'
#' @param n sample size
#' @param x covariates
#' @param b coefficients of corresponding covariates
#' @keywords internal
#' @noRd
.gen_poi <- function(n, x, b) {
  u <- b %*% t(x)
  # Generate outcome y
  y <- vector(length = n)
  y <- rpois(n, exp(u))
  y
}

#' Generate data y from Negative binomial distribution
#'
#' '.gen_nb' is used to generate data y from Negative binomial distribution
#'
#' @param n sample size
#' @param x covariates
#' @param b coefficients of corresponding covariates
#' @param size dispersion parameter
#' @keywords internal
#'
#' @noRd
.gen_nb <- function(n, x, b, size) {
  u <- b %*% t(x)
  # Generate outcome y
  y <- vector(length = n)
  y <- rnbinom(n = n, mu = exp(u), size = size)
  y
}

#' Generate data y from zero-inflated Poisson (ZIP) distribution
#'
#' '.gen_zip' is used to generate data y from ZIP distribution
#'
#' @param n sample size
#' @param zx covariates for zero inflation model
#' @param bz coefficients of corresponding covariates for zero inflation model
#' @param cx covariates for Poisson count model
#' @param bc coefficients of corresponding covariates for Poisson count model
#'
#' @keywords internal
#' @noRd
.gen_zip <- function(n, zx, bz, cx, bc) {
  # probability of zero inflation is a function of zx
  h_p0 <- bz %*% t(zx)
  p0 <- exp(h_p0) / (1 + exp(h_p0))

  zero_inf <- rbinom(n, 1, prob = p0)

  # count model (Poisson mean) is a function of cx
  h_lambda <- bc %*% t(cx)
  lambda <- exp(h_lambda)

  # Generate outcome y
  y <- vector(length = n)
  for (i in 1:n) {
    if (zero_inf[i])
      y[i] <- 0 else y[i] <- rpois(1, lambda[i])
  }
  y
}
#' Generate data y from zero-inflated Negative Binomial (ZINB) distribution
#'
#' '.gen_zinb' is used to generate data y from ZINB distribution
#'
#' @param n sample size
#' @param zx covariates for zero inflation model
#' @param bz coefficients of corresponding covariates for zero inflation model
#' @param cx covariates for Negative Binomial count model
#' @param bc coefficients of corresponding covariates for Negative Binomial
#'   count model
#' @param size dispersion parameter (the shape parameter of the gamma mixing
#'   distribution)
#'
#' @keywords internal
#' @noRd
.gen_zinb <- function(n, zx, bz, cx, bc, size) {
  # probability of zero inflation is a function of zx
  h_p0 <- bz %*% t(zx)
  p0 <- exp(h_p0) / (1 + exp(h_p0))

  zero_inf <- rbinom(n, 1, prob = p0)

  # Negative binomial mean is a function of cx
  h_mu <- bc %*% t(cx)
  mu <- exp(h_mu)

  # Generate outcome y
  y <- vector(length = n)
  for (i in 1:n) {
    if (zero_inf[i]) {
      y[i] <- 0
    }
    else {
      y[i] <- rnegbin(1, mu = mu[i], theta = size)
    }
  }

  y
}
## function to compute true effects from simulation data
# ds= simulation data with true effects on individual subject level
#' Generate data y from zero-inflated Negative Binomial (ZINB) distribution
#'
#' '.compute_te' is used to compute true effects (true average mediation effect,
#' true average direct effect, and true total effect) from a simulation data,
#' and them as a list.
#'
#' @param ds simulation data generated from .simulate1
#' @keywords internal
#' @noRd
.compute_te <- function(ds) {
  # average mediation effect when fixing treatment at 0
  d0 <- mean(ds$delta0)
  # average mediation effect when fixing treatment at 1
  d1 <- mean(ds$delta1)
  # average direct effect when the mediator is set at its level
  # under treatment = 0
  z0 <- mean(ds$zeta0)
  # average direct effect when the mediator is set at its level
  # under treatment = 1
  z1 <- mean(ds$zeta1)

  # average total effect
  tau.coef <- (d1 + d0 + z1 + z0) / 2
  # proportion via mediation when fixing treatment at 0
  n0 <- d0 / tau.coef
  # proportion via mediation when fixing treatment at 1
  n1 <- d1 / tau.coef
  # average mediation effect
  d.avg <- (d1 + d0) / 2
  # average direct effect
  z.avg <- (z1 + z0) / 2
  # proportion via mediation (Avg.)
  n.avg <- (n0 + n1) / 2

  out <- list(d0 = d0, d1 = d1, z0 = z0, z1 = z1, tau.coef = tau.coef, n0 = n0,
              n1 = n1, d.avg = d.avg, z.avg = z.avg, n.avg = n.avg)
  out
}
