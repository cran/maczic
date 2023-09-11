#' Mediation Analysis for Count and Zero-Inflated Count Data Using Instrumental
#'  Variable Method
#'
#'  \loadmathjax{} 'mediate_iv' is used to estimate controlled or natural
#' direct effect ratio, indirect (mediation) effect ratio for count and
#' zero-inflated count data using Instrumental Variable (IV) method
#' described in Guo, Z., Small, D.S., Gansky, S.A., Cheng, J. (2018).
#'
#' @details 'mediate_iv' is to estimate the controlled and natural direct and
#'   indirect ratios of randomized treatment around and through mediators when
#'   there is unmeasured confounding. This function considers count and
#'   zero-inflated count outcomes, and binary and continuous mediators. The
#'   estimation uses an instrumental variable approach with empirical
#'   likelihood and estimating equations, where no parametric assumption is
#'   needed for the outcome distribution.
#'
#'   \emph{Controlled effect ratio}
#'
#'   Given the generalized linear model for the expected potential outcomes:
#'     \mjsdeqn{f[E(Y(z,m)|x,u)] = \beta_0 + \beta_zz + \beta_mm + \beta_xx + u,}
#'   where \mjseqn{f} is a log-link function, \mjseqn{\beta_0} is the intercept,
#'   \mjseqn{z} is the treatment, \mjseqn{m} is pre-specified mediator value, \mjseqn{x}
#'   are baseline covariates, \mjseqn{\beta_z}, \mjseqn{\beta_m}, \mjseqn{\beta_x} are
#'   their corresponding coefficients, \mjseqn{u} are unmeasured confounders and
#'   \mjseqn{Y(z,m)} is the potential outcome under treatment \mjseqn{z} and mediator
#'   \mjseqn{m}.
#'   Then the controlled direct effect ratio when the treatment changes from
#'   \mjseqn{z^\ast} to \mjseqn{z} while keeping the mediator at \mjseqn{m}:
#'
#'     \mjsdeqn{\frac{E[Y(z,m)|x,u]}{E[Y(z^\ast,m)|x,u]} = \exp[\beta_z(z-z^\ast)]}
#'     and the controlled indirect effect ratio when the mediator changes from
#'     \mjseqn{m^\ast} to \mjseqn{m} while keeping the treatment at \mjseqn{z}:
#'     \mjsdeqn{\frac{E[Y(z,m)|x,u]}{E[Y(z,m^\ast)|x,u]} = \exp[\beta_m(m-m^\ast)]}
#'
#'   \emph{Natural effect ratio}
#'
#'   Given a generalized linear model for the expected potential outcomes and
#'   a model for the mediator:
#'   \mjsdeqn{f[E(Y(z,M^{z^\ast}|x,u))] = \beta_0 + \beta_zz + \beta_mM^{z^\ast} + \beta_xx
#'    + u,}
#'   where \mjseqn{f} is a log-link function, \mjseqn{M^{z^\ast}} is the potential
#'   mediator under treatment \mjseqn{z^\ast}, and \mjseqn{Y(z,M^{z^\ast})} is
#'   the potential outcome under treatment \mjseqn{z} and potential mediator
#'   \mjseqn{M^{z^\ast}}
#'   \mjsdeqn{h[E(M^{z^\ast}|x,u)] = \alpha_0 + \alpha_zz^\ast + \alpha_xx +
#'   \alpha_{IV}z^\ast x^{IV} + u,}
#'   where \mjseqn{h} is an identity function for continuous mediators and logit
#'   function for binary mediators; \mjseqn{\alpha_{IV}} is the coefficient for
#'   instrumental variables \mjseqn{x^{IV}}, \mjseqn{z^\ast x^{IV}}, interaction of randomized
#'   treatment and instrumental variables, and \mjseqn{M^z} is the potential mediator
#'   under treatment \mjseqn{z^\ast}.
#'   Then the natural direct effect ratio is:
#'  \mjsdeqn{\frac{E[Y(z,M^{z^\ast})|x,u]}{E[Y(z^\ast,M^{z^\ast}|x,u]} = \exp[\beta_z(z-z^\ast)].}
#'  And the natural indirect effect ratio with a continuous mediator is:
#'  \mjsdeqn{\frac{E[Y(z,M^z)|x,u]}{E[Y(z,M^{z^\ast})|x,u]} = \exp[\beta_m\alpha_z(z-z^\ast)
#'  + \beta_m\alpha_{IV}x^{IV}(z-z^\ast)].}
#'  The natural indirect effect ratio with a binary mediator is:
#'  \mjsdeqn{\frac{E[Y(z,M^z)|x,u]}{E[Y(z,M^{z^\ast})|x,u]}
#'   = \frac{P(M^z = 1|x,u)\exp(\beta_m) + P(M^z = 0|x, u)}
#'   {P(M^{z^\ast} = 1|x,u) \exp(\beta_m) + P(M^{z^\ast} = 0|x, u)}. }
#'
#'   The parameters are estimated by the maximized empirical likelihood
#'   estimates, which are solutions to the estimating equations that maximizes
#'   the empirical likelihood. See details in Guo, Z., Small, D.S., Gansky,
#'   S.A., Cheng, J. (2018).
#'
#' @param y Outcome variable.
#' @param z Treatment variable.
#' @param m Mediator variable.
#' @param zm.int A logical value. if 'FALSE' the interaction of treatment
#'   by mediator will not be included in the outcome model; if 'TRUE' the interaction
#'   will be included in the outcome model. Default is 'FALSE'.
#' @param ydist A character string indicating outcome distribution which will
#'   be used to fit for getting the initial parameter estimates. Can only take
#'   one of the these three values - 'poisson', 'negbin', or 'neyman type a'
#'   Default is 'poisson'. Use 'neyman type a' for zero-inflated outcomes.
#' @param mtype A character string indicating mediator type, either 'binary'
#'   or 'continuous'. Default is 'binary'.
#' @param x.IV Specifies baseline variable(s) to construct Instrumental Variable(s),
#'   Use \code{cbind()} for 2 or more IVs.
#' @param x.nonIV Specifies baseline variable(s) that are not used to construct
#'   Instrumental Variable(s) as covariates. Use \code{cbind()} for 2 or more non
#'   IVs. Default is NULL.
#' @param tol Numeric tolerance value for parameter convergence. When the square
#'   root of sum of squared difference between two consecutive set of optimized
#'   parameters less than \code{tol}, then the convergence is reached.
#'   Default is 1e-5.
#' @param n.init Number of set of initial parameter values. Default is 3.
#' @param control A list of parameters governing the nonlinear optimization
#'   algorithm behavior (to minimize the profile empirical log-likelihood).
#'   This list is the same as that for \code{BBoptim()}, which is a wrapper
#'   for \code{spg()}. \code{control} default is
#'   \code{list(maxit = 1500, M = c(50, 10), ftol=1e-10, gtol = 1e-05,
#'    maxfeval = 10000, maximize = FALSE, trace = FALSE, triter = 10,
#'    eps = 1e-07, checkGrad=NULL)}.  See \code{spg} for more details.
#' @param conf.level Level of the returned two-sided bootstrap confidence
#'   intervals. The default value, 0.95, is to return the 2.5 and 97.5
#'   percentiles of the simulated quantities.
#' @param sims Number of Monte Carlo draws for nonparametric bootstrap.
#'   Default is 1000.
#'
#' @return \code{mediate_iv} returns a list with the following components and
#'   prints them out in a matrix format:
#'   \item{beta}{parameter estimates from optimized profile empirical
#'   log-likelihood.}
#'   \item{ter, ter.ci}{point estimate for total effect ratio
#'   and its bootstrap confidence interval.}
#'
#'   If the interaction between treatment and mediator is not specified, ie,
#'   zm.int=FALSE, a list of additional components are:
#'   \item{der, der.ci}{point estimate for controlled direct effect ratio and
#'   its bootstrap confidence interval.}
#'   \item{der.nat, der.nat.ci}{point estimate for natural direct effect
#'   ratio and its bootstrap confidence interval.}
#'   \item{ier, ier.ci}{point estimate for controlled indirect effect ratio
#'   and its bootstrap confidence interval.}
#'   \item{ier.nat, ier.nat.ci}{point estimate for natural indirect effect ratio
#'   and its bootstrap confidence interval; They are only applicable to
#'   continuous mediator.}
#'
#'   If the interaction between treatment and mediator is specified, ie,
#'   zm.int=TRUE, which is only supported for binary mediator, a list of
#'   additional components are:
#'   \item{der.m0, der.m0.ci}{point estimate for controlled direct effect ratio
#'   when m = 0 and its bootstrap confidence interval.}
#'   \item{der.m1, der.m1.ci}{point estimate for controlled direct effect ratio
#'   when m = 1 and its bootstrap confidence interval.}
#'   \item{der.nat.m0, der.nat.m0.ci}{point estimate for natural direct effect
#'   ratio when m = 0 and its bootstrap confidence interval.}
#'   \item{der.nat.m1, der.nat.m1.ci}{point estimate for natural direct effect
#'   ratio when m = 1 and its bootstrap confidence interval.}
#'   \item{ier.z0, ier.z0.ci}{point estimate for controlled indirect effect
#'   ratio when z = 0 and its bootstrap confidence interval.}
#'   \item{ier.z1, ier.z1.ci}{point estimate for controlled indirect effect
#'   ratio when z = 1 and its bootstrap confidence interval.}
#'
#' @author Nancy Cheng,
#'   \email{Nancy.Cheng@@ucsf.edu}; Zijian Guo,
#'   \email{zijguo@@stat.rutgers.edu}; Jing Cheng,
#'   \email{Jing.Cheng@@ucsf.edu}.
#'
#' @seealso \code{\link{el.test}}, \code{\link{BBoptim}}, \code{\link{spg}}
#'
#' @references
#' Guo, Z., Small, D.S., Gansky, S.A., Cheng, J. (2018), Mediation analysis
#'   for count and zero-inflated count data without sequential ignorability
#'   and its application in dental studies. Journal of the Royal Statistical
#'   Society, Series C.; 67(2):371-394.
#'
#'  Ismail AI, Ondersma S, Willem Jedele JM, et al. (2011) Evaluation of
#'  a brief tailored motivational intervention to prevent early childhood
#'  caries.
#'  Community Dentistry and Oral Epidemiology 39: 433–448.
#'
#' @export
#' @examples
#' # For illustration purposes a small number of bootstrap iterations are used
#' data("midvd_bt100")
#' # The outcome is Poisson distribution
#' ee <- mediate_iv(y = midvd_bt100$Untreated_W3,
#'                  z = midvd_bt100$intervention,
#'                  m = midvd_bt100$PBrush_6, ydist = "poisson",
#'                  mtype = "binary", x.IV = midvd_bt100$BrushTimes_W2,
#'                  tol = 0.5, n.init = 1,
#'                  control = list(maxit = 15, ftol = 0.5, gtol = 0.5,
#'                                 trace = FALSE),
#'                  sims = 3)
#'
mediate_iv <- function(y, z, m, zm.int = FALSE, ydist = "poisson",
                     mtype = "binary", x.IV, x.nonIV = NULL, tol = 1e-5,
                     n.init = 3, control = list(), conf.level = 0.95,
                     sims = 1000) {
  x.IV <- cbind(x.IV)
  x.nonIV <- cbind(x.nonIV)

  n <- length(y)

  if (sims == 0) {
    stop("Invalid bootstrap iteration number")
  }

  for (b in 1:(sims + 1)) { # Bootstrap loop begins
    if (b != sims + 1) {
      print(paste("Bootstrap iteration:", b), quote = FALSE)
    }

    if (b == sims + 1) {
      #use original data for the last iteration
      index <- 1:n
    }

    repeat {
      #record number of warnings before resampling
      n_warnings_pre <- length(warnings())

      #Resampling Step
      if (b != sims + 1) {
        index <- sample(1:n, n, replace = TRUE)
      }
      if (is.null(x.nonIV)) {
        x.nonIV.b <- NULL
      }
      else {
        x.nonIV.b <- x.nonIV[index, ]
      }

      x.IV.b <- cbind(x.IV[index, ])
      m.b <- m[index]
      z.b <- z[index]
      y.b <- y[index]

      if (tolower(mtype) == "binary") {
        if (is.null(x.nonIV.b)) {
          model.m <- glm(m.b ~ z.b + x.IV.b + I(x.IV.b * z.b), family = binomial)
        }
        else {
          model.m <- glm(m.b ~ z.b + x.IV.b + x.nonIV.b + I(x.IV.b * z.b),
                         family = binomial)
        }

      } else if (tolower(mtype) == "continuous") {
        if (is.null(x.nonIV.b)) {
          model.m <- lm(m.b ~  z.b + x.IV.b + I(x.IV.b * z.b))
        }
        else {
          model.m <- lm(m.b ~  z.b + x.IV.b + x.nonIV.b + I(x.IV.b * z.b))
        }
      }
      else {
        stop("Unsupported mediator type")
      }

      if (!zm.int) { #without z*m interaction
        if (tolower(ydist) == "poisson") {
          if (is.null(x.nonIV.b)) {
            model.y <- glm(y.b ~ z.b + m.b + x.IV.b +
                             resid(model.m, type = "response"),
                           family = poisson)
          } else {
            model.y <- glm(y.b ~ z.b + m.b + x.IV.b + x.nonIV.b +
                             resid(model.m, type = "response"),
                           family = poisson)
          }
        }
        else if (tolower(ydist) == "negbin" ||
                 tolower(ydist) == "neyman type a") {
          if (is.null(x.nonIV.b)) {
            model.y <- glm.nb(y.b ~ z.b + m.b + x.IV.b +
                                resid(model.m, type = "response"))
          }
          else {
            model.y <- glm.nb(y.b ~ z.b + m.b + x.IV.b + x.nonIV.b +
                                resid(model.m, type = "response"))
          }
        }
      }
      else if (zm.int)  { #with z*m interaction
        if (tolower(ydist) == "poisson") {
          if (is.null(x.nonIV.b)) {
            model.y <- glm(y.b ~ z.b + m.b + I(z.b * m.b) + x.IV.b +
                             resid(model.m, type = "response"),
                           family = poisson)
          }
          else {
            model.y <- glm(y.b ~ z.b + m.b + I(z.b * m.b) + x.IV.b + x.nonIV +
                             resid(model.m, type = "response"),
                           family = poisson)
          }
        }
        else if (tolower(ydist) == "negbin" ||
                 tolower(ydist) == "neyman type a") {
          if (is.null(x.nonIV.b)) {
            model.y <- glm.nb(y.b ~ z.b + m.b + I(z.b * m.b) + x.IV.b +
                                resid(model.m, type = "response"))
          }
          else {
            model.y <- glm.nb(y.b ~ z.b + m.b + I(z.b * m.b) + x.IV.b + x.nonIV.b
                              + resid(model.m, type = "response"))
          }
        }
      }

      #record number of warnings again after all the model fittings
      n_warnings_post <- length(warnings())
      #if there are no new warnings then the resampling data for a bootstrap
      #iteration would be kept
      if (n_warnings_post==n_warnings_pre || b == sims + 1) {
        break
      }
    }

    ncoef.model.y <- length(coef(model.y)) - 1

    #run with (n.init) set of initial values
    for (i in 1: n.init) {
      if (i == 1) {
        beta.start.b <- coef(model.y)[1: ncoef.model.y]
      }
      else {
        beta.start.b <- coef(model.y)[1: ncoef.model.y] +
                        rnorm(ncoef.model.y, 0, 0.1)
      }

      ee.i1 <- .mediate_iv1(y = y.b, z = z.b, m = m.b, zm.int = zm.int,
                            ydist = ydist, mtype = mtype, model.m = model.m,
                            x.IV = x.IV.b, x.nonIV = x.nonIV.b,
                            beta.start = beta.start.b, tol = tol,
                            control = control)
      ee.i1.ee <- matrix(unlist(ee.i1[-1]), ncol = length(unlist(ee.i1[-1])))

      ee.i1.beta <- matrix(unlist(ee.i1[1]$beta),
                           ncol = length(unlist(ee.i1[1]$beta)))

      if (i == 1) {
        ee.i <- NULL
        ee.i <- cbind(i, ee.i1.ee)
        beta.i <- NULL
        beta.i <- cbind(i, ee.i1.beta)

      } else {
        ee.i <- rbind(ee.i, cbind(i, ee.i1.ee))
        beta.i <- rbind(beta.i, cbind(i, ee.i1.beta))
      }
      colnames(ee.i) <- c("i", names(unlist(ee.i1[-1])))
      colnames(beta.i) <- c("i", colnames(ee.i1[1]$beta))

    }
    # choose the set of effect ratio estimates that are from the converged
    # betas with smallest neg2llr among those started from the n.init sets
    # of initial values
    ee1 <- subset(ee.i, ee.i[, "neg2llr"] == min(ee.i[, "neg2llr"]))

    beta1 <- subset(beta.i, beta.i[, "i"] == ee1[, "i"])
    ee1 <- subset(ee1, select = -c(i))
    beta1 <- subset(beta1, select = -c(i))

    if (b == 1) {
      ee <- cbind(b, ee1)
      beta.mat <- cbind(b, beta1)
    } else {
      ee <- rbind(ee, cbind(b, ee1))
      beta.mat <- rbind(beta.mat, cbind(b, beta1))
    }

  }  # Bootstrap loop ends

  data.ee <- subset(ee, ee[, "b"] == sims + 1, select = -c(b))
  ter <- as.numeric(data.ee[, "ter"])
  data.beta <- subset(beta.mat, beta.mat[, 'b'] == sims + 1, select = -c(b))

  #get and correct variable names for betas
  beta.names <- gsub("beta.|[()]|.b", "", colnames(ee.i1[1]$beta))
  # correct Z * m interaction name
  beta.names <- gsub("Iz", "z", beta.names)

  colnames(data.beta) <- beta.names

  #compute boostrap confindence interval
  cfl <- 100 * conf.level
  boot.ee <- subset(ee, ee[, "b"] != sims + 1, select = -c(b))

  boot.beta <- subset(beta.mat, beta.mat[, "b"] != sims + 1, select = -c(b))
  colnames(boot.beta) <-  colnames(data.beta)
  low <- (1 - conf.level) / 2
  high <- 1 - low
  ter.ci <- quantile(boot.ee[, "ter"], c(low, high), na.rm = TRUE)
  boot.beta.ci <- apply(boot.beta, 2, quantile, c(low, high))

  #create summary matrix for beta
  beta.smat <- cbind(t(data.beta), t(boot.beta.ci))
  rownames(beta.smat) <- colnames(data.beta)
  colnames(beta.smat) <- c("Estimate", paste(cfl, "% CI Lower", sep = ""),
                           "Upper")

  if (!zm.int) { #no z-m interaction
    if (tolower(mtype) == "binary") {
      der <- as.numeric(data.ee[, "der"])
      der.nat <- as.numeric(data.ee[, "der.nat"])
      ier <- as.numeric(data.ee[, "ier"])
      der.ci <- quantile(boot.ee[, "der"], c(low, high), na.rm = TRUE)
      der.nat.ci <- quantile(boot.ee[, "der.nat"], c(low, high), na.rm = TRUE)
      ier.ci <- quantile(boot.ee[, "ier"], c(low, high), na.rm = TRUE)

      out <- list(beta = data.beta, der = der, der.ci = der.ci,
                  der.nat = der.nat, der.nat.ci = der.nat.ci,
                  ier = ier, ier.ci = ier.ci,
                  ter = ter, ter.ci = ter.ci)
      #create summary matrix for effect ratios
      er.smat <- c(der, der.ci)
      er.smat <- rbind(er.smat, c(der.nat, der.nat.ci))
      er.smat <- rbind(er.smat, c(ier, ier.ci))
      er.smat <- rbind(er.smat, c(ter, ter.ci))

      rownames(er.smat) <- c("Controlled direct effect ratio",
                             "Natural direct effect ratio",
                             "Controlled indirect effect ratio",
                             "Total effect ratio")
      colnames(er.smat) <- c("Estimate", paste(cfl, "% CI Lower", sep = ""),
                             "Upper")

    }
    else if (tolower(mtype) == "continuous") {

      der <- as.numeric(data.ee[, "der"])
      der.nat <- as.numeric(data.ee[, "der.nat"])
      ier <- as.numeric(data.ee[, "ier"])
      ier.nat <- as.numeric(data.ee[, "ier.nat"])

      der.ci <- quantile(boot.ee[, "der"], c(low, high), na.rm = TRUE)
      der.nat.ci <- quantile(boot.ee[, "der.nat"], c(low, high), na.rm = TRUE)
      ier.ci <- quantile(boot.ee[, "ier"], c(low, high), na.rm = TRUE)
      ier.nat.ci <- quantile(boot.ee[, "ier.nat"], c(low, high), na.rm = TRUE)

      out <- list(beta = data.beta, der = der, der.ci = der.ci,
                  der.nat = der.nat, der.nat.ci = der.nat.ci,
                  ier = ier, ier.ci = ier.ci,
                  ier.nat = ier.nat, ier.nat.ci = ier.nat.ci,
                  ter = ter, ter.ci = ter.ci)

      #create summary matrix for effect ratios
      er.smat <- c(der, der.ci)
      er.smat <- rbind(er.smat, c(der.nat, der.nat.ci))
      er.smat <- rbind(er.smat, c(ier, ier.ci))
      er.smat <- rbind(er.smat, c(ier.nat, ier.nat.ci))
      er.smat <- rbind(er.smat, c(ter, ter.ci))

      rownames(er.smat) <- c("Controlled direct effect ratio",
                             "Natural direct effect ratio",
                             "Controlled indirect effect ratio",
                             "Natural indirect effect ratio",
                             "Total effect ratio")

      colnames(er.smat) <- c("Estimate", paste(cfl, "% CI Lower", sep = ""),
                             "Upper")

    } else {
      stop("Unsupported mediator type")
    }
  } else { #with z x m interaction
    if (tolower(mtype) == "binary") {
      der.m0 <- as.numeric(data.ee[, "der.m0"])
      der.m1 <- as.numeric(data.ee[, "der.m1"])
      der.nat.m0 <- as.numeric(data.ee[, "der.nat.m0"])
      der.nat.m1 <- as.numeric(data.ee[, "der.nat.m1"])
      ier.z0 <- as.numeric(data.ee[, "ier.z0"])
      ier.z1 <- as.numeric(data.ee[, "ier.z1"])

      der.m0.ci <- quantile(boot.ee[, "der.m0"], c(low, high), na.rm = TRUE)
      der.m1.ci <- quantile(boot.ee[, "der.m1"], c(low, high), na.rm = TRUE)
      der.nat.m0.ci <- quantile(boot.ee[, "der.nat.m0"], c(low, high),
                                na.rm = TRUE)
      der.nat.m1.ci <- quantile(boot.ee[, "der.nat.m1"], c(low, high),
                                na.rm = TRUE)

      ier.z0.ci <- quantile(boot.ee[, "ier.z0"], c(low, high), na.rm = TRUE)
      ier.z1.ci <- quantile(boot.ee[, "ier.z1"], c(low, high), na.rm = TRUE)

      out <- list(beta = data.beta,
                  der.m0 = der.m0, der.m0.ci = der.m0.ci,
                  der.m1 = der.m1, der.m1.ci = der.m1.ci,
                  der.nat.m0 = der.nat.m0, der.nat.m0.ci = der.nat.m0.ci,
                  der.nat.m1 = der.nat.m1, der.nat.m1.ci = der.nat.m1.ci,
                  ier.z0 = ier.z0, ier.z0.ci = ier.z0.ci,
                  ier.z1 = ier.z1, ier.z1.ci = ier.z1.ci,
                  ter = ter, ter.ci = ter.ci)

      #create summary matrix for effect ratios
      er.smat <- c(der.m0, der.m0.ci)
      er.smat <- rbind(er.smat, c(der.m1, der.m1.ci))
      er.smat <- rbind(er.smat, c(der.nat.m0, der.nat.m0.ci))
      er.smat <- rbind(er.smat, c(der.nat.m1, der.nat.m1.ci))
      er.smat <- rbind(er.smat, c(ier.z0, ier.z0.ci))
      er.smat <- rbind(er.smat, c(ier.z1, ier.z1.ci))
      er.smat <- rbind(er.smat, c(ter, ter.ci))

      rownames(er.smat) <- c("Controlled direct effect ratio when m = 0",
                             "Controlled direct effect ratio when m = 1",
                             "Natural direct effect ratio when m = 0",
                             "Natural direct effect ratio when m = 1",
                             "Controlled indirect effect ratio when z = 0",
                             "Natural indirect effect ratio when z = 1",
                             "Total effect ratio")
      colnames(er.smat) <- c("Estimate", paste(cfl, "% CI Lower", sep = ""),
                             "Upper")

    }
  }

  writeLines("\nMediation Analysis Using IV Method")
  writeLines("\nCoefficient Estimates:")

  beta.smat <- round(beta.smat, 4)
  print(beta.smat)
  writeLines("\nEffect Ratio Estimates:")
  writeLines("\nConfidence Intervals Based on Nonparametric Bootstrap")

  er.smat <- round(er.smat, 4)
  print(er.smat)
  return(out)
}

#' @keywords internal
#' .mediate_iv1 function is called by mediate_iv() for estimating the various
#' effect ratios using Instrumental Variable (IV) method for only 1 sample or
#' 1 bootstrap sample. It has all the same parameters as mediate_iv() except
#' for 'conf.level' and 'sims'. It returns of list of components which includes
#' 'beta' and estimated various effect ratio.
.mediate_iv1 <- function(y, z, m, zm.int, ydist, mtype, model.m = NULL, x.IV,
                       x.nonIV, beta.start, tol, control) {
  x.IV <- cbind(x.IV)
  x <- cbind(x.IV, x.nonIV)
  #get total effect
  if (tolower(ydist) == "poisson") {
    model.y.te <- glm(y ~ z + x, family = poisson)
  }
  else if (tolower(ydist) == "negbin" || tolower(ydist) == "neyman type a") {
    model.y.te <- glm.nb(y ~ z + x)
  }
  else {
    stop("Unsupported outcome distribution")
  }

  #warnings()
  ter <- as.numeric(exp(coef(model.y.te)[2]))

  ## optimize with respect to beta
  optimize_beta <- function(beta) {

    if (!zm.int) { #without z x m interaction
      covas <- cbind(1, z, m, x)
      covas.no.m <- cbind(1, z, x)
      beta.no.m <- beta[-c(3)]
    }
    else { #with z x m interaction
      covas <- cbind(1, z, m, z * m, x)
      covas.no.m <- cbind(1, z, x)
      beta.no.m <- beta[-c(3, 4)]
    }
    ## define the estimating equations
    g11 <- y * exp(-covas %*% beta) - 1
    g21 <- g11 * z
    g31 <- diag(as.vector(g11)) %*% x
    g41 <- diag(as.vector(g11)) %*% x.IV * z

    if (!zm.int) {
      g1 <- y * exp(-beta[3] * m) - exp(covas.no.m %*% beta.no.m)
    } else {
      g1 <- y * exp(-beta[3] * m - beta[4] * z * m) -
            exp(covas.no.m %*% beta.no.m)
    }

    g2 <- g1 * z
    g3 <- diag(as.vector(g1)) %*% x
    g4 <- diag(as.vector(g1)) %*% x.IV * z

    if (tolower(mtype) == "binary") {
      # Chen et al. J of Computational and Graphical Statistics 2008
      tmp <- -1
      g11 <- rbind(g11, colMeans(g11) * (-tmp))
      g21 <- rbind(g21, colMeans(g21) * (-tmp))
      g31 <- rbind(g31, colMeans(g31) * (-tmp))
      g41 <- rbind(g41, colMeans(g41) * (-tmp))

      g1 <- rbind(g1, colMeans(g1) * (-tmp))
      g2 <- rbind(g2, colMeans(g2) * (-tmp))
      g3 <- rbind(g3, colMeans(g3) * (-tmp))
      g4 <- rbind(g4, colMeans(g4) * (-tmp))

    }
    estimating.equations <- cbind(g11, g21, g31, g41, g1, g2, g3, g4)

    n.ees <- ncol(estimating.equations)

    mu <- rep(0, n.ees)
    ## Compute the empirical likelihood ratio
    # with the mean vector fixed at mu
    el.result <- el.test(estimating.equations, mu)
    ## record the -2*loglikelihood of the ratio
    el.neg2llr <- el.result$"-2LLR" + 10^{-6}
    el.neg2llr
  }
  beta.old <- beta.start
  flag <- TRUE
  while (flag) {
    beta.optim <- BBoptim(par = beta.old, fn = optimize_beta,
                          control = control)
    beta.new <- beta.optim$par
    distance <- sqrt(sum((beta.new - beta.old)^2))
    if (distance < tol) {
      flag <- FALSE
    } else {
      beta.old <- beta.new
    }
  }
  #save optimal betas and neg2llr, which is beta.optim$value
  beta.out <- cbind(t(beta.new))
  neg2llr <- beta.optim$value

  beta.z <- as.numeric(beta.out[2])
  beta.m <- as.numeric(beta.out[3])

  #calculate effect ratios
  if (!zm.int) { #without z x m interaction
    der <- exp(beta.z)
    ier <- exp(beta.m)
    if (tolower(mtype) == "binary") {
      der.nat <- der
      out <- list(beta = beta.out, neg2llr = neg2llr, der = der,
                  der.nat = der.nat, ier = ier, ter = ter)

    }
    else if (tolower(mtype) == "continuous") {
      der.nat <- der
      n.x.IV <- ncol(x.IV)
      if (n.x.IV == 1) {
        x.IV.mu <- mean(x.IV)
        ier.nat <- as.numeric(exp(beta.m * (coef(model.m)[2] +
                                  coef(model.m)[4] * x.IV.mu)))
      }
      else {
        x.IV.mu <- colMeans(x.IV)
        ier.nat <- as.numeric(exp(beta.m * (coef(model.m)[2] +
                                  coef(model.m)[(3 + n.x.IV):(2 + 2 * n.x.IV)]
                                  %*% x.IV.mu)))
      }
      out <- list(beta = beta.out, neg2llr = neg2llr, der = der,
                  der.nat = der.nat, ier = ier, ier.nat = ier.nat, ter = ter)
    }
    else {
      stop("Unsupported mediator type")
    }
  }
  else { #with z x m interaction
    beta.zm <- beta.out[4]
    if (tolower(mtype) == "binary") {
      der.m0 <- exp(beta.z)
      der.m1 <- exp(beta.z + beta.zm)
      der.nat.m0 <- der.m0
      der.nat.m1 <- der.m1
      ier.z0 <- exp(beta.m)
      ier.z1 <- exp(beta.m + beta.zm)
      out <- list(beta = beta.out, neg2llr = neg2llr, der.m0 = der.m0,
                  der.m1 =der.m1, der.nat.m0 = der.nat.m0,
                  der.nat.m1 = der.nat.m1, ier.z0 = ier.z0,
                  ier.z1 = ier.z1, ter = ter)

    } else if (tolower(mtype) == "continuous") {
      stop("Treatment-mediator interaction for continuous mediator is not supported")
    }
  }
  return(out)
}
