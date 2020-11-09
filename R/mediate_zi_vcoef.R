#' Mediation Sensitivity Analysis for Count and Zero-Inflated Count Data with
#' a Post-Treatment Confounder
#'
#'  \loadmathjax{} 'mediate_zi_vcoef' is modified from \code{mediate_zi} function with 3
#' confounder-related parameters ('model.u', 'delta.beta.u', and 'confounder')
#' added. It is used to estimate causal mediation effects when there is a
#' treatment-induced mediator-outcome confounder, and the coefficient of treatment
#' in the confounder model is specified by users. Users can perform sensitivity
#' analysis with a range of specified coefficient values when there is a
#' post-treatment confounder.

#' @param model.u A fitted model object for confounder. Can be of class 'lm',
#'   'polr', 'bayespolr', 'glm', 'bayesglm', 'gam', 'rq', or 'survreg'.
#' @param delta.beta.u Sensitivity parameter as difference from the estimated
#'   treatment coefficient in the confounder model (model.u) based on the
#'   observed data.
#' @param confounder A character string indicating the name of the confounder
#'   variable used in the models.
#' @inheritParams mediate_zi
#'
#' @return \code{mediate_zi_vcoef} returns an object of class "\code{mediate}", (or
#'   "\code{mediate.order}" if the outcome model used is 'polr' or 'bayespolr'),
#'   a list that contains the components listed below.  Some of these elements
#'   are not available if 'long' is set to 'FALSE' by the user.
#'
#'   The function \code{summary} (i.e., \code{summary.mediate},
#'   or \code{summary.mediate.order}) can be used to obtain a table of the
#'   results.  The function \code{plot} (i.e., \code{plot.mediate}, or
#'   \code{plot.mediate.order}) can be used to produce a plot of the estimated
#'   average causal mediation, average direct, and total effects along with
#'   their confidence intervals.
#'
#'   \item{d0, d1}{point estimates for average causal mediation effects under
#'   the control and treatment conditions.}
#'   \item{d0.ci, d1.ci}{confidence intervals for average causal mediation
#'   effects. The confidence level is set at the value specified in
#'   'conf.level'.}
#'   \item{d0.p, d1.p}{two-sided p-values for average causal mediation effects.}
#'   \item{d0.sims, d1.sims}{vectors of length 'sims' containing simulation
#'   draws of average causal mediation effects.}
#'   \item{z0, z1}{point estimates for average direct effect under the control
#'   and treatment conditions.}
#'   \item{z0.ci, z1.ci}{confidence intervals for average direct effects.}
#'   \item{z0.p, z1.p}{two-sided p-values for average causal direct effects.}
#'   \item{z0.sims, z1.sims}{vectors of length 'sims' containing simulation
#'   draws of average direct effects.}
#'   \item{n0, n1}{the "proportions mediated", or the size of the average causal
#'   mediation effects relative to the total effect.}
#'   \item{n0.ci, n1.ci}{confidence intervals for the proportions mediated.}
#'   \item{n0.p, n1.p}{two-sided p-values for proportions mediated.}
#'   \item{n0.sims, n1.sims}{vectors of length 'sims' containing simulation
#'   draws of the proportions mediated.}
#'   \item{tau.coef}{point estimate for total effect.}
#'   \item{tau.ci}{confidence interval for total effect.}
#'   \item{tau.p}{two-sided p-values for total effect.}
#'   \item{tau.sims}{a vector of length 'sims' containing simulation draws of
#'   the total effect.}
#'   \item{d.avg, z.avg, n.avg}{simple averages of d0 and d1, z0 and z1, n0 and
#'   n1, respectively, which users may want to use as summary values when those
#'   quantities differ.}
#'   \item{d.avg.ci, z.avg.ci, n.avg.ci}{confidence intervals for the above.}
#'   \item{d.avg.p, z.avg.p, n.avg.p}{two-sided p-values for the above.}
#'   \item{d.avg.sims, z.avg.sims, n.avg.sims}{vectors of length 'sims'
#'   containing simulation draws of d.avg, z.avg and n.avg, respectively.}
#'   \item{boot}{logical, the 'boot' argument used.}
#'   \item{treat}{a character string indicating the name of the 'treat' variable
#'   used.}
#'   \item{mediator}{a character string indicating the name of the 'mediator'
#'   variable used.}
#'   \item{INT}{a logical value indicating whether the model specification
#'   allows the effects to differ between the treatment and control conditions.}
#'   \item{conf.level}{the confidence level used. }
#'   \item{model.y}{the outcome model used.}
#'   \item{model.m}{the mediator model used.}
#'   \item{control.value}{value of the treatment variable used as the control
#'   condition.}
#'   \item{treat.value}{value of the treatment variable used as the treatment
#'   condition.}
#'   \item{nobs}{number of observations in the model frame for 'model.m' and
#'   'model.y'. May differ from the numbers in the original models input to
#'   'mediate' if 'dropobs' was 'TRUE'.}
#'   \item{robustSE}{`TRUE' or `FALSE'.}
#'   \item{cluster}{the clusters used.}
#'
#' @author Nancy Cheng,
#'   \email{Nancy.Cheng@@ucsf.edu}; Jing Cheng,
#'   \email{Jing.Cheng@@ucsf.edu}.
#'
#' @seealso \code{\link{plot_sensitivity}}, \code{\link{mediate_zi}},
#'   \code{\link{summary.mediate}}, \code{\link{plot.mediate}}
#'
#' @references
#' Cheng, J., Cheng, N.F., Guo, Z., Gregorich, S., Ismail, A.I.,
#'   Gansky, S.A (2018) Mediation analysis for count and zero-inflated count
#'   data. Statistical Methods in Medical Research. 27(9):2756-2774.
#'
#'   Ismail AI, Ondersma S, Willem Jedele JM, et al. (2011) Evaluation of
#'   a brief tailored motivational intervention to prevent early childhood caries.
#'   Community Dentistry and Oral Epidemiology 39: 433–448.
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
#' data("midvd_bt100")
#' uFit <- glm (PDVisit_6 ~ intervention + BrushTimes_W2 + HealthyMeals_W2
#'                          + PDVisit_W2,
#'              family = 'binomial', data = midvd_bt100)
#' mFit <- glm (PBrushBedt_6 ~ intervention + BrushTimes_W2 + HealthyMeals_W2
#'                             + PBrush_W2 + PDVisit_6,
#'              family = 'binomial', data = midvd_bt100)
#' yFit <- zeroinfl(Untreated_W3 ~ intervention + PBrushBedt_6 + BrushTimes_W2
#'                                 + HealthyMeals_W2 + PBrush_W2+ PDVisit_6,
#'                  data = midvd_bt100)
#' # For illustration purposes a small number of simulations are used
#' ee <-mediate_zi_vcoef(uFit, delta.beta.u = 0.01, mFit, yFit, sims = 100,
#'                       treat = "intervention", mediator = "PBrushBedt_6",
#'                       confounder ="PDVisit_6")
#' summary(ee)
#'
mediate_zi_vcoef <- function(model.u, delta.beta.u, model.m, model.y,
                             sims = 1000, boot = FALSE,
                             treat = "treat.name", mediator = "med.name",
                             confounder ="confd.name",
                             covariates = NULL, outcome = NULL,
                             control = NULL, conf.level = .95,
                             control.value = 0, treat.value = 1,
                             long = TRUE, dropobs = FALSE,
                             robustSE = FALSE, cluster = NULL, ...) {

  # Warn users who still use INT option
  if (match("INT", names(match.call()), 0L)) {
    warning("'INT' is deprecated - existence of interaction terms is now automatically detected from model formulas")
  }

  # Warning for robustSE and cluster used with boot
  if (robustSE && boot) {
    warning("'robustSE' is ignored for nonparametric bootstrap")
  }

  if (!is.null(cluster) && boot) {
    warning("'cluster' is ignored for nonparametric bootstrap")
  }

  if (robustSE & !is.null(cluster)) {
    stop("choose either `robustSE' or `cluster' option, not both")
  }

  # Drop observations not common to all confounder, mediator and outcome models
  if (dropobs) {
    odata.u <- model.frame(model.u)
    odata.m <- model.frame(model.m)
    odata.y <- model.frame(model.y)
    odata.um <- merge(odata.u, odata.m, sort = FALSE,
                      by = c("row.names", intersect(names(odata.u), names(odata.m))))
    rownames(odata.um) <- odata.um$Row.names
    newdata <- merge(odata.um[, -1L], odata.y, sort = FALSE,
                     by = c("row.names", intersect(names(odata.um), names(odata.y))))

    rownames(newdata) <- newdata$Row.names
    newdata <- newdata[, -1L]
    rm(odata.u, odata.m, odata.um, odata.y)

    call.u <- getCall(model.u)
    call.m <- getCall(model.m)
    call.y <- getCall(model.y)

    call.u$data <- call.m$data <- call.y$data <- newdata
    if (c("(weights)") %in% names(newdata)) {
      call.u$weights <- call.m$weights <- call.y$weights <- model.weights(newdata)
    }
    model.u <- eval.parent(call.u)
    model.m <- eval.parent(call.m)
    model.y <- eval.parent(call.y)
  }

  # Model type indicators
  isGam.y <- inherits(model.y, "gam")
  isGam.u <- inherits(model.u, "gam")
  isGam.m <- inherits(model.m, "gam")
  isGlm.y <- inherits(model.y, "glm")  # Note gam and bayesglm also inherits "glm"
  isGlm.u <- inherits(model.u, "glm")  # Note gam and bayesglm also inherits "glm"
  isGlm.m <- inherits(model.m, "glm")  # Note gam and bayesglm also inherits "glm"
  isLm.y <- inherits(model.y, "lm")    # Note gam, glm and bayesglm also inherit "lm"
  isLm.u <- inherits(model.u, "lm")    # Note gam, glm and bayesglm also inherit "lm"
  isLm.m <- inherits(model.m, "lm")    # Note gam, glm and bayesglm also inherit "lm"
  isVglm.y <- inherits(model.y, "vglm")
  isRq.y <- inherits(model.y, "rq")
  isRq.u <- inherits(model.u, "rq")
  isRq.m <- inherits(model.m, "rq")
  isOrdered.y <- inherits(model.y, "polr")  # Note bayespolr also inherits "polr"
  isOrdered.u <- inherits(model.u, "polr")  # Note bayespolr also inherits "polr"
  isOrdered.m <- inherits(model.m, "polr")  # Note bayespolr also inherits "polr"
  isSurvreg.y <- inherits(model.y, "survreg")
  isSurvreg.u <- inherits(model.u, "survreg")
  isSurvreg.m <- inherits(model.m, "survreg")
  isZeroinfl.y <- inherits(model.y, "zeroinfl")

  # Record family of model.u if glm
  if (isGlm.u) {
    FamilyU <- model.u$family$family
  }
  # Record family of model.m if glm
  if (isGlm.m) {
    FamilyM <- model.m$family$family
  }

  # Record vfamily of model.y if vglm (currently only tobit)
  if (isVglm.y) {
    VfamilyY <- model.y@family@vfamily
  }

  # Warning for unused options
  if (!is.null(control) && !isGam.y) {
    warning("'control' is only used for GAM outcome models - ignored")
    control <- NULL
  }
  if (!is.null(outcome) && !(isSurvreg.y && boot)) {
    warning("'outcome' is only relevant for survival outcome models with bootstrap - ignored")
  }

  # Model frames for U, M and Y models
  u.data <- model.frame(model.u)  # Call.U$data
  m.data <- model.frame(model.m)  # Call.M$data
  y.data <- model.frame(model.y)  # Call.Y$data

  # Numbers of observations and categories
  n.u <- nrow(u.data)
  n.m <- nrow(m.data)
  n.y <- nrow(y.data)

  if (n.u != n.m || n.u != n.y) {
    stop("number of observations do not match between confounder, mediator and outcome models")
  } else {
    n <- n.u
  }
  m <- length(sort(unique(model.frame(model.m)[, 1])))
  u <- length(sort(unique(model.frame(model.u)[, 1])))

  # Extracting weights from models
  weights.u <- model.weights(u.data)
  weights.m <- model.weights(m.data)
  weights.y <- model.weights(y.data)

  if (!is.null(weights.u) && isGlm.u && FamilyU == "binomial") {
    message("weights taken as sampling weights, not total number of trials")
  }
  if (!is.null(weights.m) && isGlm.m && FamilyM == "binomial") {
    message("weights taken as sampling weights, not total number of trials")
  }

  if (is.null(weights.u)) {
    weights.u <- rep(1, nrow(u.data))
  }
  if (is.null(weights.m)) {
    weights.m <- rep(1, nrow(m.data))
  }
  if (is.null(weights.y)) {
    weights.y <- rep(1, nrow(y.data))
  }
  if (!all(weights.u == weights.y)) {
    stop("weights on outcome and confounder models not identical")
  } else if (!all(weights.m == weights.y)) {
    stop("weights on outcome and mediator models not identical")
  }
  else {
    weights <- weights.u
  }

  # Convert character treatment to factor
  if (is.character(u.data[, treat])) {
    u.data[, treat] <- factor(u.data[, treat])
  }
  if (is.character(m.data[, treat])) {
    m.data[, treat] <- factor(m.data[, treat])
  }
  if (is.character(y.data[, treat])) {
    y.data[, treat] <- factor(y.data[, treat])
  }

  # Convert character confounder to factor
  if (is.character(y.data[, confounder])) {
    y.data[, confounder] <- factor(y.data[, confounder])
  }
  # Convert character mediator to factor
  if (is.character(y.data[, mediator])) {
    y.data[, mediator] <- factor(y.data[, mediator])
  }

  # Factor treatment indicator
  isFactorT.u <- is.factor(u.data[, treat])
  isFactorT.m <- is.factor(m.data[, treat])
  isFactorT.y <- is.factor(y.data[, treat])
  if (isFactorT.u != isFactorT.m || isFactorT.u != isFactorT.y) {
    stop("treatment variable types differ in confounder, mediator and outcome models")
  } else {
    isFactorT <- isFactorT.y
  }

  if (isFactorT) {
    t.levels <- levels(y.data[, treat])
    if (treat.value %in% t.levels & control.value %in% t.levels) {
      cat.0 <- control.value
      cat.1 <- treat.value
    } else {
      cat.0 <- t.levels[1]
      cat.1 <- t.levels[2]
      warning("treatment and control values do not match factor levels; using ", cat.0, " and ", cat.1, " as control and treatment, respectively")
    }
  } else {
    cat.0 <- control.value
    cat.1 <- treat.value
  }

  # Factor confounder indicator
  isFactorU <- is.factor(y.data[, confounder])

  if (isFactorU) {
    u.levels <- levels(y.data[, confounder])
  }

  # Factor mediator indicator
  isFactorM <- is.factor(y.data[, mediator])

  if (isFactorM) {
    m.levels <- levels(y.data[, mediator])
  }
  #####################################
  ## Define functions
  #####################################
  indexmax <- function(x) {
    ## Return position of largest element in vector x
    order(x)[length(x)]
  }
  ###########################################################################
  ### CASE I: EVERYTHING EXCEPT ORDERED OUTCOME
  ###########################################################################
  if (!isOrdered.y) {

    #######################################################################
    ## Case I-1: Quasi-Bayesian Monte Carlo
    #######################################################################
    if (!boot) {
      # Error if gam outcome or quantile confounder
      if (isGam.u | isGam.y | isRq.u) {
        stop("'boot' must be 'TRUE' for models used")
      }

      # Get mean and variance parameters for confounder simulations
      if (isSurvreg.u && is.null(survreg.distributions[[model.u$dist]]$scale)) {
        UModel.coef <- c(coef(model.u), log(model.u$scale))
        scalesim.u <- TRUE
      } else {
        beta.u <- coef(model.u)[treat]
        UModel.coef <- coef(model.u)
        UModel.coef[treat] < -beta.u + delta.beta.u
        scalesim.u <- FALSE
      }

      if (isOrdered.u) {
        if (is.null(model.u$Hess)) {
          cat("Confounder model object does not contain 'Hessian';")
        }
        k <- length(UModel.coef)
        UModel.var.cov <- vcov(model.u)[(1:k), (1:k)]
      } else if (isSurvreg.u) {
        UModel.var.cov <- vcov(model.u)
      } else {
        if (robustSE) {
          UModel.var.cov <- vcovHC(model.u, ...)
        } else if (!is.null(cluster)) {
          if (nrow(u.data) != length(cluster)) {
            warning("length of cluster vector differs from # of obs for mediator model")
          }
          dta <- merge(u.data, as.data.frame(cluster), sort = FALSE,
                       by = "row.names")
          fm <- update(model.u, data = dta)
          UModel.var.cov <- sandwich::vcovCL(fm, dta[, ncol(dta)])

        } else {
          UModel.var.cov <- vcov(model.u)
        }
      }

      # Error if gam outcome or quantile mediator
      if (isGam.m | isGam.y | isRq.m) {
        stop("'boot' must be 'TRUE' for models used")
      }
      # Get mean and variance parameters for mediator simulations
      if (isSurvreg.m && is.null(survreg.distributions[[model.m$dist]]$scale)) {
        MModel.coef <- c(coef(model.m), log(model.m$scale))
        scalesim.m <- TRUE
      } else {
        MModel.coef <- coef(model.m)
        scalesim.m <- FALSE
      }

      if (isOrdered.m) {
        if (is.null(model.m$Hess)) {
          cat("Mediator model object does not contain 'Hessian';")
        }
        k <- length(MModel.coef)
        MModel.var.cov <- vcov(model.m)[(1:k), (1:k)]
      } else if (isSurvreg.m) {
        MModel.var.cov <- vcov(model.m)
      } else {
        if (robustSE) {
          MModel.var.cov <- vcovHC(model.m, ...)
        } else if (!is.null(cluster)) {
          if (nrow(m.data) != length(cluster)) {
            warning("length of cluster vector differs from # of obs for mediator model")
          }
          dta <- merge(m.data, as.data.frame(cluster), sort = FALSE,
                       by = "row.names")
          fm <- update(model.m, data = dta)
          MModel.var.cov <- sandwich::vcovCL(fm, dta[, ncol(dta)])
        } else {
          MModel.var.cov <- vcov(model.m)
        }
      }

      # Get mean and variance parameters for outcome simulations
      if (isSurvreg.y && is.null(survreg.distributions[[model.y$dist]]$scale)) {
        YModel.coef <- c(coef(model.y), log(model.y$scale))
        scalesim.y <- TRUE  # indicates if survreg scale parameter is simulated
      } else {
        YModel.coef <- coef(model.y)
        scalesim.y <- FALSE
      }

      if (isRq.y) {
        YModel.var.cov <- summary(model.y, covariance = TRUE)$cov
      } else if (isSurvreg.y) {
        YModel.var.cov <- vcov(model.y)
      } else {
        if (robustSE) {
          YModel.var.cov <- vcovHC(model.y, ...)
        } else if (!is.null(cluster)) {
          if (nrow(y.data) != length(cluster)) {
            warning("length of cluster vector differs from # of obs for outcome model")
          }
          dta <- merge(y.data, as.data.frame(cluster), sort = FALSE,
                       by = "row.names")
          fm <- update(model.y, data = dta)
          YModel.var.cov <- sandwich::vcovCL(fm, dta[, ncol(dta)])

        } else {
          YModel.var.cov <- vcov(model.y)
        }
      }

      # Draw model coefficients from normal
      if (sum(is.na(UModel.coef), is.na(MModel.coef), is.na(YModel.coef)) > 0) {
        stop("NA in model coefficients; rerun models with nonsingular design matrix")
      }

      UModel <- mvrnorm(sims, mu = UModel.coef, Sigma = UModel.var.cov)

      # Draw model coefficients from normal
      if (sum(is.na(MModel.coef), is.na(YModel.coef)) > 0) {
        stop("NA in model coefficients; rerun models with nonsingular design matrix")
      }

      MModel <- mvrnorm(sims, mu = MModel.coef, Sigma = MModel.var.cov)

      if (isZeroinfl.y) {
        model.y.coef.count <- coef(model.y, model = "count")
        model.y.vcov.count <- vcov(model.y, model = "count")
        model.y.coef.zero <- coef(model.y, model = "zero")
        model.y.vcov.zero <- vcov(model.y, model ="zero")

        YModel.coefficients.count <- mvrnorm(sims, mu = model.y.coef.count, Sigma = model.y.vcov.count)
        YModel.coefficients.zero <- mvrnorm(sims, mu = model.y.coef.zero, Sigma = model.y.vcov.zero)

      } else {
        YModel <- mvrnorm(sims, mu = YModel.coef, Sigma = YModel.var.cov)
      }
      if (robustSE && (isSurvreg.u | isSurvreg.m | isSurvreg.y)) {
        warning("`robustSE' ignored for survival models; fit the model with `robust' option instead\n")
      }
      if (!is.null(cluster) && (isSurvreg.u | isSurvreg.m | isSurvreg.y)) {
        warning("`cluster' ignored for survival models; fit the model with 'cluster()' term in the formula\n")
      }

      #####################################
      #  Confounder Predictions
      #####################################
      pred.data.t <- pred.data.c <- u.data

      if (isFactorT) {
        pred.data.t[, treat] <- factor(cat.1, levels = t.levels)
        pred.data.c[, treat] <- factor(cat.0, levels = t.levels)
      } else {
        pred.data.t[, treat] <- cat.1
        pred.data.c[, treat] <- cat.0
      }

      if (!is.null(covariates)) {
        for (p in 1:length(covariates)) {
          vl <- names(covariates[p])
          x <- ifelse(is.factor(pred.data.t[, vl]),
                      factor(covariates[[p]], levels = levels(u.data[, vl])),
                      covariates[[p]])
          pred.data.t[, vl] <- pred.data.c[, vl] <- x
        }
      }

      umat.t <- model.matrix(terms(model.u), data = pred.data.t)
      umat.c <- model.matrix(terms(model.u), data = pred.data.c)

      ### Case I-1-a: GLM Confounder
      if (isGlm.u) {

        muU1 <- model.u$family$linkinv(tcrossprod(UModel, umat.t))
        muU0 <- model.u$family$linkinv(tcrossprod(UModel, umat.c))

        if (FamilyU == "poisson") {
          PredictU1 <- matrix(rpois(sims * n, lambda = muU1), nrow = sims)
          PredictU0 <- matrix(rpois(sims * n, lambda = muU0), nrow = sims)
        } else if (FamilyU == "Gamma") {
          shape <- gamma.shape(model.u)$alpha
          PredictU1 <- matrix(rgamma(n * sims, shape = shape,
                                     scale = muU1/shape), nrow = sims)
          PredictU0 <- matrix(rgamma(n * sims, shape = shape,
                                     scale = muU0/shape), nrow = sims)

        } else if (FamilyU == "binomial") {
          PredictU1 <- matrix(rbinom(n * sims, size = 1,
                                     prob = muU1), nrow = sims)
          PredictU0 <- matrix(rbinom(n * sims, size = 1,
                                     prob = muU0), nrow = sims)
        } else if (FamilyU == "gaussian") {
          sigma <- sqrt(summary(model.u)$dispersion)
          error <- rnorm(sims * n, mean = 0, sd = sigma)
          PredictU1 <- muU1 + matrix(error, nrow = sims)
          PredictU0 <- muU0 + matrix(error, nrow = sims)
          rm(sigma, error)
        } else if (FamilyU == "inverse.gaussian") {
          disp <- summary(model.u)$dispersion
          PredictU1 <- matrix(SuppDists::rinvGauss(n * sims, nu = muU1,
                                                   lambda = 1/disp), nrow = sims)
          PredictU0 <- matrix(SuppDists::rinvGauss(n*sims, nu = muU0,
                                                   lambda = 1/disp), nrow = sims)
        } else {
          stop("unsupported glm family")
        }
        rm(muU1, muU0)

        ### Case I-1-b: Ordered confounder
      } else if (isOrdered.u) {
        if (model.u$method == "logistic") {
          linkfn <- plogis
        } else if (model.u$method == "probit") {
          linkfn <- pnorm
        } else {
          stop("unsupported polr method; use 'logistic' or 'probit'")
        }

        #u.cat <- sort(unique(model.frame(model.u)[, 1]))
        lambda <- model.u$zeta

        umat.t <- umat.t[, -1]
        umat.c <- umat.c[, -1]

        ystar_u1 <- tcrossprod(UModel, umat.t)
        ystar_u0 <- tcrossprod(UModel, umat.c)

        PredictU1 <- matrix(NA, nrow = sims, ncol = n)
        PredictU0 <- matrix(NA, nrow = sims, ncol = n)

        for (i in 1:sims) {

          cprobs_u1 <- matrix(NA, n, u)
          cprobs_u0 <- matrix(NA, n, u)
          probs_u1 <- matrix(NA, n, u)
          probs_u0 <- matrix(NA, n, u)

          for (j in 1:(u - 1)) {  # loop to get category-specific probabilities
            cprobs_u1[, j] <- linkfn(lambda[j] - ystar_u1[i, ])
            cprobs_u0[, j] <- linkfn(lambda[j] - ystar_u0[i, ])
            # cumulative probabilities
            probs_u1[, u] <- 1 - cprobs_u1[, u-1] # top category
            probs_u0[, u] <- 1-cprobs_u0[, u-1] # top category
            probs_u1[, 1] <- cprobs_u1[, 1]     # bottom category
            probs_u0[, 1] <- cprobs_u0[, 1]     # bottom category
          }

          for (j in 2:(u - 1)) {  # middle categories
            probs_u1[, j] <- cprobs_u1[, j] - cprobs_u1[, j - 1]
            probs_u0[, j] <- cprobs_u0[, j] - cprobs_u0[, j - 1]
          }

          draws_u1 <- matrix(NA, n, u)
          draws_u0 <- matrix(NA, n, u)

          for (ii in 1:n) {
            draws_u1[ii, ] <- t(rmultinom(1, 1, prob = probs_u1[ii, ]))
            draws_u0[ii, ] <- t(rmultinom(1, 1, prob = probs_u0[ii, ]))
          }

          PredictU1[i, ] <- apply(draws_u1, 1, indexmax)
          PredictU0[i, ] <- apply(draws_u0, 1, indexmax)
        }
        rm(umat.t, umat.c)

        ### Case I-1-c: Linear
      } else if (isLm.u) {
        sigma <- summary(model.u)$sigma
        error <- rnorm(sims * n, mean = 0, sd = sigma)
        muU1 <- tcrossprod(UModel, umat.t)
        muU0 <- tcrossprod(UModel, umat.c)
        PredictU1 <- muU1 + matrix(error, nrow = sims)
        PredictU0 <- muU0 + matrix(error, nrow = sims)
        rm(error)

        ### Case I-1-d: Survreg
      } else if (isSurvreg.u) {
        dd <- survreg.distributions[[model.u$dist]]
        if (is.null(dd$itrans)) {
          itrans <- function(x) x
        } else {
          itrans <- dd$itrans
        }
        dname <- dd$dist
        if (is.null(dname)) {
          dname <- model.u$dist
        }
        if (scalesim.u) {
          scale <- exp(UModel[, ncol(UModel)])
          lpU1 <- tcrossprod(UModel[, 1:(ncol(UModel) - 1)], umat.t)
          lpU0 <- tcrossprod(UModel[, 1:(ncol(UModel) - 1)], umat.c)
        } else {
          scale <- dd$scale
          lpU1 <- tcrossprod(UModel, umat.t)
          lpU0 <- tcrossprod(UModel, umat.c)
        }
        error <- switch(dname,
                        extreme = log(rweibull(sims * n, shape = 1, scale = 1)),
                        gaussian = rnorm(sims * n),
                        logistic = rlogis(sims * n),
                        t = rt(sims * n, df = dd$parms))
        PredictU1 <- itrans(lpU1 + scale * matrix(error, nrow = sims))
        PredictU0 <- itrans(lpU0 + scale * matrix(error, nrow = sims))
        rm(error)

      } else {
        stop("confounder model is not yet implemented")
      }
      rm(umat.t, umat.c)

      #####################################
      #  Mediator Predictions
      #####################################
      PredictM1 <- matrix(NA, nrow = sims, ncol = n)
      PredictM0 <- matrix(NA, nrow = sims, ncol = n)
      for (jj in 1:sims) {

        pred.data.t <- pred.data.c <- m.data

        if (isFactorT) {
          pred.data.t[, treat] <- factor(cat.1, levels = t.levels)
          pred.data.c[, treat] <- factor(cat.0, levels = t.levels)
        } else {
          pred.data.t[, treat] <- cat.1
          pred.data.c[, treat] <- cat.0
        }

        if (!is.null(covariates)) {
          for (p in 1:length(covariates)) {
            vl <- names(covariates[p])
            x <- ifelse(is.factor(pred.data.t[, vl]),
                        factor(covariates[[p]], levels = levels(m.data[, vl])),
                        covariates[[p]])
            pred.data.t[, vl] <- pred.data.c[, vl] <- x
          }
        }

        if (isFactorU) {
          pred.data.t[, confounder] <- factor(PredictU1[jj, ], levels = 1:u, labels = u.levels)
          pred.data.c[, confounder] <- factor(PredictU0[jj, ], levels = 1:u, labels = u.levels)
        } else {
          pred.data.t[, confounder] <- PredictU1[jj, ]
          pred.data.c[, confounder] <- PredictU0[jj, ]
        }

        mmat.t <- model.matrix(terms(model.m), data = pred.data.t)
        mmat.c <- model.matrix(terms(model.m), data = pred.data.c)

        ### Case I-1-a: GLM Mediator
        if (isGlm.m) {
          muM1 <- model.m$family$linkinv(tcrossprod(MModel[jj, ], mmat.t))
          muM0 <- model.m$family$linkinv(tcrossprod(MModel[jj, ], mmat.c))

          if (FamilyM == "poisson") {
            PredictM1[jj, ] <- rpois(n, lambda = muM1)
            PredictM0[jj, ] <- rpois(n, lambda = muM0)
          } else if (FamilyM == "Gamma") {
            shape <- gamma.shape(model.m)$alpha
            PredictM1[jj, ] <- rgamma(n, shape = shape, scale = muM1/shape)
            PredictM0[jj, ] <- rgamma(n, shape = shape, scale = muM0/shape)
          } else if (FamilyM == "binomial") {
            PredictM1[jj, ] <- rbinom(n, size = 1, prob = muM1)
            PredictM0[jj, ] <- rbinom(n, size = 1, prob = muM0)
          } else if (FamilyM == "gaussian") {
            sigma <- sqrt(summary(model.m)$dispersion)
            error <- rnorm(n, mean = 0, sd = sigma)
            PredictM1[jj, ] <- muM1 + error
            PredictM0[jj, ] <- muM0 + error
            rm(error)
          } else if (FamilyM == "inverse.gaussian") {
            disp <- summary(model.m)$dispersion
            PredictM1[jj, ] <- SuppDists::rinvGauss(n, nu = muM1, lambda = 1 / disp)
            PredictM0[jj, ] <- SuppDists::rinvGauss(n, nu = muM0, lambda = 1 / disp)

          } else {
            stop("unsupported glm family")
          }
          rm(muM1, muM0)
          ### Case I-1-b: Ordered mediator
        } else if (isOrdered.m) {
          if (model.m$method == "logistic") {
            linkfn <- plogis
          } else if (model.m$method == "probit") {
            linkfn <- pnorm
          } else {
            stop("unsupported polr method; use 'logistic' or 'probit'")
          }

          #m.cat <- sort(unique(model.frame(model.m)[, 1]))
          lambda <- model.m$zeta

          mmat.t <- mmat.t[, -1]
          mmat.c <- mmat.c[, -1]

          ystar_m1 <- tcrossprod(MModel[jj, ], mmat.t)
          ystar_m0 <- tcrossprod(MModel[jj, ], mmat.c)

          cprobs_m1 <- matrix(NA, n, m)
          cprobs_m0 <- matrix(NA, n, m)
          probs_m1 <- matrix(NA, n, m)
          probs_m0 <- matrix(NA, n, m)

          for (j in 1:(m - 1)) {  # loop to get category-specific probabilities
            cprobs_m1[, j] <- linkfn(lambda[j] - ystar_m1[jj, ])
            cprobs_m0[, j] <- linkfn(lambda[j] - ystar_m0[jj, ])
            # cumulative probabilities
            probs_m1[, m] <- 1 - cprobs_m1[, m-1] # top category
            probs_m0[, m] <- 1 - cprobs_m0[, m-1] # top category
            probs_m1[, 1] <- cprobs_m1[, 1]     # bottom category
            probs_m0[, 1] <- cprobs_m0[, 1]     # bottom category
          }

          for (j in 2:(m - 1)){  # middle categories
            probs_m1[, j] <- cprobs_m1[, j] - cprobs_m1[, j-1]
            probs_m0[, j] <- cprobs_m0[, j] - cprobs_m0[, j-1]
          }

          draws_m1 <- matrix(NA, n, m)
          draws_m0 <- matrix(NA, n, m)

          for (ii in 1:n) {
            draws_m1[ii, ] <- t(rmultinom(1, 1, prob = probs_m1[ii, ]))
            draws_m0[ii, ] <- t(rmultinom(1, 1, prob = probs_m0[ii, ]))
          }

          PredictM1[jj, ] <- apply(draws_m1, 1, indexmax)
          PredictM0[jj, ] <- apply(draws_m0, 1, indexmax)
          rm(mmat.t, mmat.c)

          ### Case I-1-c: Linear
        } else if (isLm.m) {
          sigma <- summary(model.m)$sigma
          error <- rnorm(n, mean = 0, sd = sigma)
          muM1 <- tcrossprod(MModel[jj, ], mmat.t)
          muM0 <- tcrossprod(MModel[jj, ], mmat.c)
          PredictM1[jj, ] <- muM1 + error
          PredictM0[jj, ] <- muM0 + error
          rm(error)

          ### Case I-1-d: Survreg
        } else if (isSurvreg.m) {
          dd <- survreg.distributions[[model.m$dist]]
          if (is.null(dd$itrans)) {
            itrans <- function(x) x
          } else {
            itrans <- dd$itrans
          }
          dname <- dd$dist
          if (is.null(dname)) {
            dname <- model.m$dist
          }
          if (scalesim.m) {
            scale <- exp(MModel[jj, ncol(MModel)])
            lpM1 <- tcrossprod(MModel[jj, 1:(ncol(MModel) - 1)], mmat.t)
            lpM0 <- tcrossprod(MModel[jj, 1:(ncol(MModel) - 1)], mmat.c)
          } else {
            scale <- dd$scale
            lpM1 <- tcrossprod(MModel[jj, ], mmat.t)
            lpM0 <- tcrossprod(MModel[jj, ], mmat.c)
          }
          error <- switch(dname,
                          extreme = log(rweibull(n, shape = 1, scale = 1)),
                          gaussian = rnorm(n),
                          logistic = rlogis(n),
                          t = rt(n, df = dd$parms))
          PredictM1[jj, ] <- itrans(lpM1 + scale * error)
          PredictM0[jj, ] <- itrans(lpM0 + scale * error)
          rm(error)

        } else {
          stop("mediator model is not yet implemented")
        }
        rm(mmat.t, mmat.c)
      }#end of loop - for(jj in 1:sims)

      #####################################
      ##  Outcome Predictions
      #####################################
      effects.tmp <- array(NA, dim = c(n, sims, 4))

      for (e in 1:4) {
        tt <- switch(e, c(1, 1, 1, 0), c(0, 0, 1, 0), c(1, 0, 1, 1), c(1, 0, 0, 0))
        Pr1 <- matrix(NA, nrow = n, ncol = sims)
        Pr0 <- matrix(NA, nrow = n, ncol = sims)

        for (j in 1:sims) {
          pred.data.t <- pred.data.c <- y.data

          if (!is.null(covariates)) {
            for (p in 1:length(covariates)) {
              vl <- names(covariates[p])
              x <- ifelse(is.factor(pred.data.t[, vl]),
                          factor(covariates[[p]], levels = levels(y.data[, vl])),
                          covariates[[p]])
              pred.data.t[, vl] <- pred.data.c[, vl] <- x
            }
          }

          # Set treatment values
          cat.t <- ifelse(tt[1], cat.1, cat.0)
          cat.c <- ifelse(tt[2], cat.1, cat.0)
          cat.t.ctrl <- ifelse(tt[1], cat.0, cat.1)
          cat.c.ctrl <- ifelse(tt[2], cat.0, cat.1)
          if (isFactorT) {
            pred.data.t[, treat] <- factor(cat.t, levels = t.levels)
            pred.data.c[, treat] <- factor(cat.c, levels = t.levels)
            if (!is.null(control)) {
              pred.data.t[, control] <- factor(cat.t.ctrl, levels = t.levels)
              pred.data.c[, control] <- factor(cat.c.ctrl, levels = t.levels)
            }
          } else {
            pred.data.t[, treat] <- cat.t
            pred.data.c[, treat] <- cat.c
            if (!is.null(control)) {
              pred.data.t[, control] <- cat.t.ctrl
              pred.data.c[, control] <- cat.c.ctrl
            }
          }

          # Set confounder values
          PredictUt <- PredictU1[j, ] * tt[1] + PredictU0[j, ] * (1 - tt[1])
          PredictUc <- PredictU1[j, ] * tt[2] + PredictU0[j, ] * (1 - tt[2])

          if (isFactorU) {
            pred.data.t[, confounder] <- factor(PredictUt, levels = 1:u, labels = u.levels)
            pred.data.c[, confounder] <- factor(PredictUc, levels = 1:u, labels = u.levels)
          } else {
            pred.data.t[, confounder] <- PredictUt
            pred.data.c[, confounder] <- PredictUc
          }

          # Set mediator values
          PredictMt <- PredictM1[j, ] * tt[3] + PredictM0[j, ] * (1 - tt[3])
          PredictMc <- PredictM1[j, ] * tt[4] + PredictM0[j, ] * (1 - tt[4])
          if (isFactorM) {
            pred.data.t[, mediator] <- factor(PredictMt, levels = 1:m, labels = m.levels)
            pred.data.c[, mediator] <- factor(PredictMc, levels = 1:m, labels = m.levels)
          } else {
            pred.data.t[, mediator] <- PredictMt
            pred.data.c[, mediator] <- PredictMc
          }

          ymat.t <- model.matrix(terms(model.y), data = pred.data.t)
          ymat.c <- model.matrix(terms(model.y), data = pred.data.c)

          if (!isZeroinfl.y) {
            #get covariate data including intercept that match design matrix
            # ie: rearrange data as the order of intercept covariate1 covariate3 ...
            ymat.t <- model.matrix(terms(model.y), data = pred.data.t)
            ymat.c <- model.matrix(terms(model.y), data = pred.data.c)
          }

          if (isVglm.y) {
            if (VfamilyY == "tobit") {
              Pr1.tmp <- ymat.t %*% YModel[j, -2]
              Pr0.tmp <- ymat.c %*% YModel[j, -2]
              Pr1[, j] <- pmin(pmax(Pr1.tmp, model.y@misc$Lower), model.y@misc$Upper)
              Pr0[, j] <- pmin(pmax(Pr0.tmp, model.y@misc$Lower), model.y@misc$Upper)
            } else {
              stop("outcome model is in unsupported vglm family")
            }
          } else if (scalesim.y) {
            Pr1[, j] <- t(as.matrix(YModel[j, 1:(ncol(YModel) - 1)])) %*% t(ymat.t)
            Pr0[, j] <- t(as.matrix(YModel[j, 1:(ncol(YModel) - 1)])) %*% t(ymat.c)
          } else if (isZeroinfl.y) {
            mf.t <- model.frame(delete.response(model.y$terms$full), pred.data.t, xlev = model.y$levels)
            mf.c <- model.frame(delete.response(model.y$terms$full), pred.data.c, xlev = model.y$levels)
            X.t <- model.matrix(delete.response(model.y$terms$count), mf.t, contrasts = model.y$contrasts$count)
            X.c <- model.matrix(delete.response(model.y$terms$count), mf.c, contrasts = model.y$contrasts$count)
            Z.t <- model.matrix(delete.response(model.y$terms$zero),  mf.t, contrasts = model.y$contrasts$zero)
            Z.c <- model.matrix(delete.response(model.y$terms$zero),  mf.c, contrasts = model.y$contrasts$zero)
            mu.t <- exp(X.t %*% YModel.coefficients.count[j, ])[, 1]
            mu.c <- exp(X.c %*% YModel.coefficients.count[j, ])[, 1]

            p.t <- model.y$linkinv(Z.t %*% YModel.coefficients.zero[j, ])[, 1]
            p.c <- model.y$linkinv(Z.c %*% YModel.coefficients.zero[j, ])[, 1]

            Pr1[, j] <- (1 - p.t) * mu.t
            Pr0[, j] <- (1 - p.c) * mu.c

            rm(X.t, X.c, Z.t, Z.c, p.t, p.c, mu.t, mu.c)

          } else {
            Pr1[, j] <- t(as.matrix(YModel[j, ])) %*% t(ymat.t)
            Pr0[, j] <- t(as.matrix(YModel[j, ])) %*% t(ymat.c)
          }
          if (isZeroinfl.y)
            rm(pred.data.t, pred.data.c)
          else rm(ymat.t, ymat.c, pred.data.t, pred.data.c)
        }

        if (isGlm.y) {
          Pr1 <- apply(Pr1, 2, model.y$family$linkinv)
          Pr0 <- apply(Pr0, 2, model.y$family$linkinv)
        } else if (isSurvreg.y) {
          dd <- survreg.distributions[[model.y$dist]]
          if (is.null(dd$itrans)) {
            itrans <- function(x) x
          } else {
            itrans <- dd$itrans
          }
          Pr1 <- apply(Pr1, 2, itrans)
          Pr0 <- apply(Pr0, 2, itrans)
        }

        effects.tmp[, , e] <- Pr1 - Pr0
        rm(Pr1, Pr0)
      }
      if (isZeroinfl.y) {
        rm(PredictU1, PredictU0, PredictM1, PredictM0, MModel, YModel.coefficients.count, YModel.coefficients.zero, UModel)
      }
      else  rm(PredictU1, PredictU0, PredictM1, PredictM0, YModel, MModel, UModel)

      delta.1 <- t(as.matrix(apply(effects.tmp[, , 1], 2, weighted.mean, w = weights)))
      delta.0 <- t(as.matrix(apply(effects.tmp[, , 2], 2, weighted.mean, w = weights)))
      zeta.1 <- t(as.matrix(apply(effects.tmp[, , 3], 2, weighted.mean, w = weights)))
      zeta.0 <- t(as.matrix(apply(effects.tmp[, , 4], 2, weighted.mean, w = weights)))
      rm(effects.tmp)
      tau <- (zeta.1 + delta.0 + zeta.0 + delta.1) / 2
      nu.0 <- delta.0 / tau
      nu.1 <- delta.1 / tau
      delta.avg <- (delta.1 + delta.0) / 2
      zeta.avg <- (zeta.1 + zeta.0) / 2
      nu.avg <- (nu.1 + nu.0) / 2

      d0 <- mean(delta.0)
      d1 <- mean(delta.1)
      z1 <- mean(zeta.1)
      z0 <- mean(zeta.0)
      tau.coef <- mean(tau)
      n0 <- median(nu.0)
      n1 <- median(nu.1)
      d.avg <- (d0 + d1) / 2
      z.avg <- (z0 + z1) / 2
      n.avg <- (n0 + n1) / 2

      ########################################################################
      ## Case I-2: Nonparametric Bootstrap
      ########################################################################
    } else {

      Call.U <- getCall(model.u)
      Call.M <- getCall(model.m)
      Call.Y <- getCall(model.y)

      # Storage
      delta.1 <- matrix(NA, sims, 1)
      delta.0 <- matrix(NA, sims, 1)
      zeta.1 <- matrix(NA, sims, 1)
      zeta.0 <- matrix(NA, sims, 1)
      tau <- matrix(NA, sims, 1)

      # Bootstrap loop begins
      for (b in 1:(sims + 1)) {
        #check if there is an error when fitting the model
        bError <- 1
        while (bError == 1) {

          index <- sample(1:n, n, replace = TRUE)

          if (b == sims + 1) {  # in the final run, use the original data
            index <- 1:n
          }

          if (isSurvreg.u) {
            if (ncol(model.u$y) > 2) {
              stop("unsupported censoring type")
            }
            uname <- names(u.data)[1]
            if (substr(uname, 1, 4) != "Surv") {
              stop("refit the survival model with `Surv' used directly in model formula")
            }
            nc <- nchar(confounder)
            eventname <- substr(uname, 5 + nc + 3, nchar(uname) - 1)
            if (nchar(eventname) == 0) {
              u.data.tmp <- data.frame(u.data,
                                       as.numeric(u.data[, 1L][, 1L]))
              names(u.data.tmp)[c(1L, ncol(u.data) + 1)] <- c(uname, confounder)
            } else {
              u.data.tmp <- data.frame(u.data,
                                       as.numeric(u.data[, 1L][, 1L]),
                                       as.numeric(model.u$y[, 2]))
              names(u.data.tmp)[c(1L, ncol(u.data) + (1:2))] <- c(uname, confounder, eventname)
            }
            Call.U$data <- u.data.tmp[index, ]
          } else {
            Call.U$data <- u.data[index, ]
          }

          if (isSurvreg.m) {
            if (ncol(model.m$y) > 2) {
              stop("unsupported censoring type")
            }
            mname <- names(m.data)[1]
            if (substr(mname, 1, 4) != "Surv") {
              stop("refit the survival model with `Surv' used directly in model formula")
            }
            nc <- nchar(mediator)
            eventname <- substr(mname, 5 + nc + 3, nchar(mname) - 1)
            if (nchar(eventname) == 0) {
              m.data.tmp <- data.frame(m.data,
                                       as.numeric(m.data[, 1L][, 1L]))
              names(m.data.tmp)[c(1L, ncol(m.data) + 1)] <- c(mname, mediator)
            } else {
              m.data.tmp <- data.frame(m.data,
                                       as.numeric(m.data[, 1L][, 1L]),
                                       as.numeric(model.m$y[, 2]))
              names(m.data.tmp)[c(1L, ncol(m.data) + (1:2))] <- c(mname, mediator, eventname)
            }
            Call.M$data <- m.data.tmp[index, ]
          } else {
            Call.M$data <- m.data[index, ]
          }

          if (isSurvreg.y) {
            if (ncol(model.y$y) > 2) {
              stop("unsupported censoring type")
            }
            yname <- names(y.data)[1]
            if (substr(yname, 1, 4) != "Surv") {
              stop("refit the survival model with `Surv' used directly in model formula")
            }
            if (is.null(outcome)) {
              stop("`outcome' must be supplied for survreg outcome with boot")
            }
            nc <- nchar(outcome)
            eventname <- substr(yname, 5 + nc + 3, nchar(yname) - 1)
            if (nchar(eventname) == 0) {
              y.data.tmp <- data.frame(y.data,
                                       as.numeric(y.data[, 1L][, 1L]))
              names(y.data.tmp)[c(1L, ncol(y.data) + 1)] <- c(yname, outcome)
            } else {
              y.data.tmp <- data.frame(y.data,
                                       as.numeric(y.data[, 1L][, 1L]),
                                       as.numeric(model.y$y[, 2]))
              names(y.data.tmp)[c(1L, ncol(y.data) + (1:2))] <- c(yname, outcome, eventname)
            }
            Call.Y$data <- y.data.tmp[index, ]
          } else {
            Call.Y$data <- y.data[index, ]
          }

          Call.U$weights <- u.data[index, "(weights)"]
          Call.M$weights <- m.data[index, "(weights)"]
          Call.Y$weights <- y.data[index, "(weights)"]

          if (isOrdered.u && length(unique(y.data[index, confounder])) != u) {
            stop("insufficient variation on confounder")
          }

          # Refit Models with Resampled Data
          new.fit.U <- eval.parent(Call.U)
          new.fit.M <- eval.parent(Call.M)
          new.fit.Y <- try(eval.parent(Call.Y), TRUE)

          bLen <- length(grep("Error", new.fit.Y[1]))
          if (bLen == 0) bError <- 0 #no error
          else bError <- 1 #error

        }

        #####################################
        #  Confounder Predictions
        #####################################
        pred.data.t <- pred.data.c <- u.data

        if (isFactorT) {
          pred.data.t[, treat] <- factor(cat.1, levels = t.levels)
          pred.data.c[, treat] <- factor(cat.0, levels = t.levels)
        } else {
          pred.data.t[, treat] <- cat.1
          pred.data.c[, treat] <- cat.0
        }

        if (!is.null(covariates)) {
          for (p in 1:length(covariates)) {
            vl <- names(covariates[p])
            x <- ifelse(is.factor(pred.data.t[, vl]),
                        factor(covariates[[p]], levels = levels(u.data[, vl])),
                        covariates[[p]])
            pred.data.t[, vl] <- pred.data.c[, vl] <- x
          }
        }

        ### Case I-2-a: GLM Confounder (including GAMs)
        if (isGlm.u) {

          muU1 <- predict(new.fit.U, type = "response", newdata = pred.data.t)
          muU0 <- predict(new.fit.U, type = "response", newdata = pred.data.c)

          if (FamilyU == "poisson") {
            PredictU1 <- rpois(n, lambda = muU1)
            PredictU0 <- rpois(n, lambda = muU0)
          } else if (FamilyU == "Gamma") {
            shape <- gamma.shape(new.fit.U)$alpha
            PredictU1 <- rgamma(n, shape = shape, scale = muU1/shape)
            PredictU0 <- rgamma(n, shape = shape, scale = muU0/shape)

          } else if (FamilyU == "binomial") {
            PredictU1 <- rbinom(n, size = 1, prob = muU1)
            PredictU0 <- rbinom(n, size = 1, prob = muU0)
          } else if (FamilyU == "gaussian") {
            sigma <- sqrt(summary(new.fit.U)$dispersion)
            error <- rnorm(n, mean = 0, sd = sigma)
            PredictU1 <- muU1 + error
            PredictU0 <- muU0 + error
            rm(error)
          } else if (FamilyU == "inverse.gaussian") {
            disp <- summary(new.fit.U)$dispersion
            PredictU1 <- SuppDists::rinvGauss(n, nu = muU1, lambda = 1/disp)
            PredictU0 <- SuppDists::rinvGauss(n, nu = muU0, lambda = 1/disp)

          } else {
            stop("unsupported glm family")
          }
          rm(muU1, muU0)

          ### Case I-2-b: Ordered Confounder
        } else if (isOrdered.u) {
          probs_u1 <- predict(new.fit.U, newdata = pred.data.t, type = "probs")
          probs_u0 <- predict(new.fit.U, newdata = pred.data.c, type = "probs")
          draws_u1 <- matrix(NA, n, u)
          draws_u0 <- matrix(NA, n, u)
          for (ii in 1:n) {
            draws_u1[ii, ] <- t(rmultinom(1, 1, prob = probs_u1[ii, ]))
            draws_u0[ii, ] <- t(rmultinom(1, 1, prob = probs_u0[ii, ]))
          }
          PredictU1 <- apply(draws_u1, 1, indexmax)
          PredictU0 <- apply(draws_u0, 1, indexmax)

          ### Case I-2-c: Quantile Regression for Confounder
        } else if (isRq.u) {
          # Use inverse transform sampling to predict U
          call.new <- new.fit.U$call
          call.new$tau <- runif(n)
          newfits <- eval.parent(call.new)
          tt <- delete.response(terms(new.fit.U))
          u.t <- model.frame(tt, pred.data.t, xlev = new.fit.U$xlevels)
          u.c <- model.frame(tt, pred.data.c, xlev = new.fit.U$xlevels)
          X.t <- model.matrix(tt, u.t, contrasts = new.fit.U$contrasts)
          X.c <- model.matrix(tt, u.c, contrasts = new.fit.U$contrasts)
          rm(call.new, tt, u.t, u.c)
          PredictU1 <- rowSums(X.t * t(newfits$coefficients))
          PredictU0 <- rowSums(X.c * t(newfits$coefficients))
          rm(newfits, X.t, X.c)

          ### Case I-2-d: Linear
        } else if (isLm.u) {
          sigma <- summary(new.fit.U)$sigma
          error <- rnorm(n, mean = 0, sd = sigma)
          PredictU1 <- predict(new.fit.U, type = "response",
                               newdata = pred.data.t) + error
          PredictU0 <- predict(new.fit.U, type = "response",
                               newdata = pred.data.c) + error
          rm(error)

          ### Case I-2-e: Survreg
        } else if (isSurvreg.u) {
          dd <- survreg.distributions[[new.fit.U$dist]]
          if (is.null(dd$itrans)) {
            itrans <- function(x) x
          } else {
            itrans <- dd$itrans
          }
          dname <- dd$dist
          if (is.null(dname)) {
            dname <- new.fit.U$dist
          }
          scale <- new.fit.U$scale
          lpU1 <- predict(new.fit.U, newdata = pred.data.t, type = "linear")
          lpU0 <- predict(new.fit.U, newdata = pred.data.c, type = "linear")
          error <- switch(dname,
                          extreme = log(rweibull(n, shape = 1, scale = 1)),
                          gaussian = rnorm(n),
                          logistic = rlogis(n),
                          t = rt(n, df = dd$parms))
          PredictU1 <- as.numeric(itrans(lpU1 + scale * error))
          PredictU0 <- as.numeric(itrans(lpU0 + scale * error))
          rm(error)

        } else {
          stop("confounder model is not yet implemented")
        }

        #####################################
        # Mediator Predictions
        #####################################
        pred.data.t <- pred.data.c <- m.data

        if (isFactorT) {
          pred.data.t[, treat] <- factor(cat.1, levels = t.levels)
          pred.data.c[, treat] <- factor(cat.0, levels = t.levels)
        } else {
          pred.data.t[, treat] <- cat.1
          pred.data.c[, treat] <- cat.0
        }

        if (!is.null(covariates)) {
          for (p in 1:length(covariates)) {
            vl <- names(covariates[p])
            x <- ifelse(is.factor(pred.data.t[, vl]),
                        factor(covariates[[p]], levels = levels(m.data[, vl])),
                        covariates[[p]])
            pred.data.t[, vl] <- pred.data.c[, vl] <- x
          }
        }
        if (isFactorU) {
          pred.data.t[, confounder] <- factor(PredictU1, levels = 1:u, labels = u.levels)
          pred.data.c[, confounder] <- factor(PredictU0, levels = 1:u, labels = u.levels)
        } else {
          pred.data.t[, confounder] <- PredictU1
          pred.data.c[, confounder] <- PredictU0
        }

        ### Case II-a: GLM Mediator (including GAMs)
        if (isGlm.m) {

          muM1 <- predict(new.fit.M, type = "response", newdata = pred.data.t)
          muM0 <- predict(new.fit.M, type = "response", newdata = pred.data.c)

          if (FamilyM == "poisson") {
            PredictM1 <- rpois(n, lambda = muM1)
            PredictM0 <- rpois(n, lambda = muM0)
          } else if (FamilyM == "Gamma") {
            shape <- gamma.shape(model.m)$alpha
            PredictM1 <- rgamma(n, shape = shape, scale = muM1/shape)
            PredictM0 <- rgamma(n, shape = shape, scale = muM0/shape)

          } else if (FamilyM == "binomial") {
            PredictM1 <- rbinom(n, size = 1, prob = muM1)
            PredictM0 <- rbinom(n, size = 1, prob = muM0)
          } else if (FamilyM == "gaussian") {
            sigma <- sqrt(summary(model.m)$dispersion)
            error <- rnorm(n, mean = 0, sd = sigma)
            PredictM1 <- muM1 + error
            PredictM0 <- muM0 + error
            rm(error)
          } else if (FamilyM == "inverse.gaussian") {
            disp <- summary(model.m)$dispersion
            PredictM1 <- SuppDists::rinvGauss(n, nu = muM1, lambda = 1/disp)
            PredictM0 <- SuppDists::rinvGauss(n, nu = muM0, lambda = 1/disp)

          } else {
            stop("unsupported glm family")
          }

          ### Case II-b: Ordered Mediator
        } else if (isOrdered.m) {
          probs_m1 <- predict(new.fit.M, type = "probs", newdata = pred.data.t)
          probs_m0 <- predict(new.fit.M, type = "probs", newdata = pred.data.c)
          draws_m1 <- matrix(NA, n, m)
          draws_m0 <- matrix(NA, n, m)

          for (ii in 1:n) {
            draws_m1[ii, ] <- t(rmultinom(1, 1, prob = probs_m1[ii, ]))
            draws_m0[ii, ] <- t(rmultinom(1, 1, prob = probs_m0[ii, ]))
          }
          PredictM1 <- apply(draws_m1, 1, indexmax)
          PredictM0 <- apply(draws_m0, 1, indexmax)

          ### Case II-c: Quantile Regression for Mediator
        } else if (isRq.m) {
          # Use inverse transform sampling to predict M
          call.new <- new.fit.M$call
          call.new$tau <- runif(n)
          newfits <- eval.parent(call.new)
          tt <- delete.response(terms(new.fit.M))
          m.t <- model.frame(tt, pred.data.t, xlev = new.fit.M$xlevels)
          m.c <- model.frame(tt, pred.data.c, xlev = new.fit.M$xlevels)
          X.t <- model.matrix(tt, m.t, contrasts = new.fit.M$contrasts)
          X.c <- model.matrix(tt, m.c, contrasts = new.fit.M$contrasts)
          rm(call.new, tt, m.t, m.c)
          PredictM1 <- rowSums(X.t * t(newfits$coefficients))
          PredictM0 <- rowSums(X.c * t(newfits$coefficients))
          rm(newfits, X.t, X.c)

          ### Case II-d: Linear
        } else if (isLm.m) {
          sigma <- summary(new.fit.M)$sigma
          error <- rnorm(n, mean = 0, sd = sigma)
          PredictM1 <- predict(new.fit.M, type = "response",
                               newdata = pred.data.t) + error
          PredictM0 <- predict(new.fit.M, type = "response",
                               newdata = pred.data.c) + error
          rm(error)

          ### Case I-2-e: Survreg
        } else if (isSurvreg.m) {
          dd <- survreg.distributions[[new.fit.M$dist]]
          if (is.null(dd$itrans)) {
            itrans <- function(x) x
          } else {
            itrans <- dd$itrans
          }
          dname <- dd$dist
          if (is.null(dname)) {
            dname <- new.fit.M$dist
          }
          scale <- new.fit.M$scale
          lpM1 <- predict(new.fit.M, newdata = pred.data.t, type = "linear")
          lpM0 <- predict(new.fit.M, newdata = pred.data.c, type = "linear")
          error <- switch(dname,
                          extreme = log(rweibull(n, shape = 1, scale = 1)),
                          gaussian = rnorm(n),
                          logistic = rlogis(n),
                          t = rt(n, df = dd$parms))
          PredictM1 <- as.numeric(itrans(lpM1 + scale * error))
          PredictM0 <- as.numeric(itrans(lpM0 + scale * error))
          rm(error)

        } else {
          stop("mediator model is not yet implemented")
        }

        #####################################
        #  Outcome Predictions
        #####################################
        effects.tmp <- matrix(NA, nrow = n, ncol = 4)

        for (e in 1:4) {
          tt <- switch(e, c(1, 1, 1, 0), c(0, 0, 1, 0), c(1, 0, 1, 1), c(1, 0, 0, 0))
          pred.data.t <- pred.data.c <- y.data

          if (!is.null(covariates)) {
            for (p in 1:length(covariates)) {
              vl <- names(covariates[p])
              x <- ifelse(is.factor(pred.data.t[, vl]),
                          factor(covariates[[p]], levels = levels(y.data[, vl])),
                          covariates[[p]])
              pred.data.t[, vl] <- pred.data.c[, vl] <- x
            }
          }

          # Set treatment values
          cat.t <- ifelse(tt[1], cat.1, cat.0)
          cat.c <- ifelse(tt[2], cat.1, cat.0)
          cat.t.ctrl <- ifelse(tt[1], cat.0, cat.1)
          cat.c.ctrl <- ifelse(tt[2], cat.0, cat.1)
          if (isFactorT) {
            pred.data.t[, treat] <- factor(cat.t, levels = t.levels)
            pred.data.c[, treat] <- factor(cat.c, levels = t.levels)
            if (!is.null(control)) {
              pred.data.t[, control] <- factor(cat.t.ctrl, levels = t.levels)
              pred.data.c[, control] <- factor(cat.c.ctrl, levels = t.levels)
            }
          } else {
            pred.data.t[, treat] <- cat.t
            pred.data.c[, treat] <- cat.c
            if (!is.null(control)) {
              pred.data.t[, control] <- cat.t.ctrl
              pred.data.c[, control] <- cat.c.ctrl
            }
          }

          # Set confounder values
          PredictUt <- PredictU1 * tt[1] + PredictU0 * (1 - tt[1])
          PredictUc <- PredictU1 * tt[2] + PredictU0 * (1 - tt[2])

          if (isFactorU) {
            pred.data.t[, confounder] <- factor(PredictUt, levels = 1:u, labels=u.levels)
            pred.data.c[, confounder] <- factor(PredictUc, levels = 1:u, labels=u.levels)
          } else {
            pred.data.t[, confounder] <- PredictUt
            pred.data.c[, confounder] <- PredictUc
          }
          # set mediator values
          PredictMt <- PredictM1 * tt[3] + PredictM0 * (1 - tt[3])
          PredictMc <- PredictM1 * tt[4] + PredictM0 * (1 - tt[4])
          if (isFactorM) {
            pred.data.t[, mediator] <- factor(PredictMt, levels = 1:m, labels = m.levels)
            pred.data.c[, mediator] <- factor(PredictMc, levels = 1:m, labels = m.levels)
          } else {
            pred.data.t[, mediator] <- PredictMt
            pred.data.c[, mediator] <- PredictMc
          }

          if (isRq.y) {
            pr.1 <- predict(new.fit.Y, type = "response",
                            newdata = pred.data.t, interval = "none")
            pr.0 <- predict(new.fit.Y, type = "response",
                            newdata = pred.data.c, interval = "none")
          } else {
            pr.1 <- predict(new.fit.Y, type = "response",
                            newdata = pred.data.t)
            pr.0 <- predict(new.fit.Y, type = "response",
                            newdata = pred.data.c)
          }

          if (isZeroinfl.y) {
            mf.t <- model.frame(delete.response(new.fit.Y$terms$full), pred.data.t, xlev = new.fit.Y$levels)
            mf.c <- model.frame(delete.response(new.fit.Y$terms$full), pred.data.c, xlev = new.fit.Y$levels)
            X.t <- model.matrix(delete.response(new.fit.Y$terms$count), mf.t, contrasts = new.fit.Y$contrasts$count)
            X.c <- model.matrix(delete.response(new.fit.Y$terms$count), mf.c, contrasts = new.fit.Y$contrasts$count)
            Z.t <- model.matrix(delete.response(new.fit.Y$terms$zero),  mf.t, contrasts = new.fit.Y$contrasts$zero)
            Z.c <- model.matrix(delete.response(new.fit.Y$terms$zero),  mf.c, contrasts = new.fit.Y$contrasts$zero)
            mu.t <- exp(X.t %*% new.fit.Y$coefficients$count)[, 1]
            mu.c <- exp(X.c %*% new.fit.Y$coefficients$count)[, 1]

            p.t <- new.fit.Y$linkinv(Z.t %*% new.fit.Y$coefficients$zero)[, 1]
            p.c <- new.fit.Y$linkinv(Z.c %*% new.fit.Y$coefficients$zero)[, 1]

            pr.1 <- (1 - p.t) * mu.t
            pr.0 <- (1 - p.c) * mu.c

            rm(mf.t, mf.c, X.t, X.c, Z.t, Z.c, mu.t, mu.c, p.t, p.c)
          }

          pr.mat <- as.matrix(cbind(pr.1, pr.0))
          effects.tmp[,e] <- pr.mat[, 1] - pr.mat[, 2]

          rm(pred.data.t, pred.data.c, pr.1, pr.0, pr.mat)
        }

        # Compute all QoIs
        if (b == sims + 1) {
          d1 <- weighted.mean(effects.tmp[, 1], weights)
          d0 <- weighted.mean(effects.tmp[, 2], weights)
          z1 <- weighted.mean(effects.tmp[, 3], weights)
          z0 <- weighted.mean(effects.tmp[, 4], weights)
        } else {
          delta.1[b] <- weighted.mean(effects.tmp[, 1], weights)
          delta.0[b] <- weighted.mean(effects.tmp[, 2], weights)
          zeta.1[b] <- weighted.mean(effects.tmp[, 3], weights)
          zeta.0[b] <- weighted.mean(effects.tmp[, 4], weights)
        }
      }  # bootstrap loop ends
      rm(effects.tmp)
      tau.coef <- (d1 + d0 + z1 + z0)/2
      n0 <- d0/tau.coef
      n1 <- d1/tau.coef
      d.avg <- (d1 + d0)/2
      z.avg <- (z1 + z0)/2
      n.avg <- (n0 + n1)/2

      tau <- (delta.1 + delta.0 + zeta.1 + zeta.0)/2
      nu.0 <- delta.0/tau
      nu.1 <- delta.1/tau
      delta.avg <- (delta.0 + delta.1)/2
      zeta.avg <- (zeta.0 + zeta.1)/2
      nu.avg <- (nu.0 + nu.1)/2

    }  # nonpara boot branch ends

    ########################################################################
    ## Compute Outputs and Put Them Together
    ########################################################################

    low <- (1 - conf.level)/2
    high <- 1 - low
    d0.ci <- quantile(delta.0, c(low, high), na.rm = TRUE)
    d1.ci <- quantile(delta.1, c(low, high), na.rm = TRUE)
    tau.ci <- quantile(tau, c(low, high), na.rm = TRUE)
    z1.ci <- quantile(zeta.1, c(low, high), na.rm = TRUE)
    z0.ci <- quantile(zeta.0, c(low, high), na.rm = TRUE)
    n0.ci <- quantile(nu.0, c(low, high), na.rm = TRUE)
    n1.ci <- quantile(nu.1, c(low, high), na.rm = TRUE)
    d.avg.ci <- quantile(delta.avg, c(low, high), na.rm = TRUE)
    z.avg.ci <- quantile(zeta.avg, c(low, high), na.rm = TRUE)
    n.avg.ci <- quantile(nu.avg, c(low, high), na.rm = TRUE)

    # p-values
    d0.p <- 2 * sum(sign(delta.0) != sign(median(delta.0)))/sims
    d1.p <- 2 * sum(sign(delta.1) != sign(median(delta.1)))/sims
    d.avg.p <- 2 * sum(sign(delta.avg) != sign(median(delta.avg)))/sims
    z0.p <- 2 * sum(sign(zeta.0) != sign(median(zeta.0)))/sims
    z1.p <- 2 * sum(sign(zeta.1) != sign(median(zeta.1)))/sims
    z.avg.p <- 2 * sum(sign(zeta.avg) != sign(median(zeta.avg)))/sims
    n0.p <- 2 * sum(sign(nu.0) != sign(median(nu.0)))/sims
    n1.p <- 2 * sum(sign(nu.1) != sign(median(nu.1)))/sims
    n.avg.p <- 2 * sum(sign(nu.avg) != sign(median(nu.avg)))/sims
    tau.p <- 2 * sum(sign(tau) != sign(median(tau)))/sims

    # Detect whether models include T-M interaction
    INT <- paste(treat, confounder, sep = ":") %in%
                 attr(terms(model.y), "term.labels") |
           paste(confounder, treat, sep = ":") %in%
                 attr(terms(model.y), "term.labels")
    if (!INT & isGam.y) {
      INT <- !isTRUE(all.equal(d0, d1))  # if gam, determine empirically
    }

    if (long) {
      out <- list(d0 = d0, d1 = d1, d0.ci = d0.ci, d1.ci = d1.ci,
                  d0.p = d0.p, d1.p = d1.p,
                  d0.sims = delta.0, d1.sims = delta.1,
                  z0 = z0, z1 = z1, z0.ci = z0.ci, z1.ci = z1.ci,
                  z0.p = z0.p, z1.p = z1.p,
                  z0.sims = zeta.0, z1.sims = zeta.1,
                  n0 = n0, n1 = n1, n0.ci = n0.ci, n1.ci = n1.ci,
                  n0.p = n0.p, n1.p = n1.p,
                  n0.sims = nu.0, n1.sims = nu.1,
                  tau.coef = tau.coef, tau.ci = tau.ci, tau.p = tau.p,
                  tau.sims = tau,
                  d.avg = d.avg, d.avg.p = d.avg.p, d.avg.ci = d.avg.ci,
                  d.avg.sims = delta.avg,
                  z.avg = z.avg, z.avg.p = z.avg.p, z.avg.ci = z.avg.ci,
                  z.avg.sims = zeta.avg,
                  n.avg = n.avg, n.avg.p = n.avg.p, n.avg.ci = n.avg.ci,
                  n.avg.sims = nu.avg,
                  boot = boot, treat = treat, confounder = confounder,
                  covariates = covariates,
                  INT = INT, conf.level = conf.level,
                  model.y = model.y, model.u = model.u,
                  control.value = control.value, treat.value = treat.value,
                  nobs = n, sims = sims)
    } else {
      out <- list(d0 = d0, d1 = d1, d0.ci = d0.ci, d1.ci = d1.ci,
                  d0.p = d0.p, d1.p = d1.p,
                  z0 = z0, z1 = z1, z0.ci = z0.ci, z1.ci = z1.ci,
                  z0.p = z0.p, z1.p = z1.p,
                  n0 = n0, n1 = n1, n0.ci = n0.ci, n1.ci = n1.ci,
                  n0.p = n0.p, n1.p = n1.p,
                  tau.coef = tau.coef, tau.ci = tau.ci, tau.p = tau.p,
                  d.avg = d.avg, d.avg.p = d.avg.p, d.avg.ci = d.avg.ci,
                  z.avg = z.avg, z.avg.p = z.avg.p, z.avg.ci = z.avg.ci,
                  n.avg = n.avg, n.avg.p = n.avg.p, n.avg.ci = n.avg.ci,
                  boot = boot, treat = treat, confounder = confounder,
                  covariates = covariates,
                  INT = INT, conf.level = conf.level,
                  model.y = model.y, model.u = model.u,
                  control.value = control.value, treat.value = treat.value,
                  nobs = n, sims = sims)
    }
    class(out) <- "mediate"
    out

    ############################################################################
    ### CASE II: ORDERED OUTCOME
    ############################################################################
    } else {
    if (boot != TRUE) {
      warning("ordered outcome model can only be used with nonparametric bootstrap - option forced")
      boot <- TRUE
    }

    n.ycat <- length(unique(model.response(y.data)))

    # Storage
    delta.1 <- matrix(NA, sims, n.ycat)
    delta.0 <- matrix(NA, sims, n.ycat)
    zeta.1 <- matrix(NA, sims, n.ycat)
    zeta.0 <- matrix(NA, sims, n.ycat)
    tau <- matrix(NA, sims, n.ycat)

    # Bootstrap loop begins
    for (b in 1:(sims + 1)) {

      # Resampling Step
      index <- sample(1:n, n, replace = TRUE)
      if (b == sims + 1) {  # use original data for the last iteration
        index <- 1:n
      }

      Call.U <- model.u$call
      Call.M <- model.m$call
      Call.Y <- model.y$call

      if (isSurvreg.u) {
        if (ncol(model.u$y) > 2) {
          stop("unsupported censoring type")
        }
        uname <- names(u.data)[1]
        if (substr(uname, 1, 4) != "Surv") {
          stop("refit the survival model with `Surv' used directly in model formula")
        }
        nc <- nchar(confounder)
        eventname <- substr(uname, 5 + nc + 3, nchar(uname) - 1)
        if (nchar(eventname) == 0) {
          u.data.tmp <- data.frame(u.data,
                                   as.numeric(u.data[, 1L][, 1L]))
          names(u.data.tmp)[c(1L, ncol(u.data) + 1)] <- c(uname, confounder)
        } else {
          u.data.tmp <- data.frame(u.data,
                                   as.numeric(u.data[, 1L][, 1L]),
                                   as.numeric(model.u$y[, 2]))
          names(u.data.tmp)[c(1L, ncol(u.data) + (1:2))] <- c(uname, confounder, eventname)
        }
        Call.U$data <- u.data.tmp[index, ]
      } else {
        Call.U$data <- u.data[index, ]
      }
      if (isSurvreg.m) {
        if (ncol(model.m$y) > 2) {
          stop("unsupported censoring type")
        }
        mname <- names(m.data)[1]
        if (substr(mname, 1, 4) != "Surv") {
          stop("refit the survival model with `Surv' used directly in model formula")
        }
        nc <- nchar(mediator)
        eventname <- substr(mname, 5 + nc + 3, nchar(mname) - 1)
        if (nchar(eventname) == 0) {
          m.data.tmp <- data.frame(m.data,
                                   as.numeric(m.data[, 1L][, 1L]))
          names(m.data.tmp)[c(1L, ncol(m.data) + 1)] <- c(mname, mediator)
        } else {
          m.data.tmp <- data.frame(m.data,
                                   as.numeric(m.data[, 1L][, 1L]),
                                   as.numeric(model.m$y[, 2]))
          names(m.data.tmp)[c(1L, ncol(m.data) + (1:2))] <- c(mname, mediator,
                                                              eventname)
        }
        Call.M$data <- m.data.tmp[index, ]
      } else {
        Call.M$data <- m.data[index, ]
      }

      Call.Y$data <- y.data[index, ]
      Call.U$weights <- u.data[index, "(weights)"]
      Call.M$weights <- m.data[index, "(weights)"]
      Call.Y$weights <- y.data[index, "(weights)"]

      new.fit.U <- eval.parent(Call.U)
      new.fit.M <- eval.parent(Call.M)
      new.fit.Y <- eval.parent(Call.Y)

      if (isOrdered.u && length(unique(y.data[index, confounder])) != u) {
        # Modify the coefficients when confounder has empty cells
        coefnames.y <- names(model.y$coefficients)
        coefnames.new.y <- names(new.fit.Y$coefficients)
        new.fit.Y.coef <- rep(0, length(coefnames.y))
        names(new.fit.Y.coef) <- coefnames.y
        new.fit.Y.coef[coefnames.new.y] <- new.fit.Y$coefficients
        new.fit.Y$coefficients <- new.fit.Y.coef
      }

      #####################################
      # Confounder Predictions
      #####################################
      pred.data.t <- pred.data.c <- u.data

      if (isFactorT) {
        pred.data.t[, treat] <- factor(cat.1, levels = t.levels)
        pred.data.c[, treat] <- factor(cat.0, levels = t.levels)
      } else {
        pred.data.t[, treat] <- cat.1
        pred.data.c[, treat] <- cat.0
      }

      if (!is.null(covariates)) {
        for (p in 1:length(covariates)) {
          vl <- names(covariates[p])
          x <- ifelse(is.factor(pred.data.t[, vl]),
                      factor(covariates[[p]], levels = levels(u.data[, vl])),
                      covariates[[p]])
          pred.data.t[, vl] <- pred.data.c[, vl] <- x
        }
      }

      ### Case II-a: GLM Confounder (including GAMs)
      if (isGlm.u) {

        muU1 <- predict(new.fit.U, type = "response", newdata = pred.data.t)
        muU0 <- predict(new.fit.U, type = "response", newdata = pred.data.c)

        if (FamilyU == "poisson") {
          PredictU1 <- rpois(n, lambda = muU1)
          PredictU0 <- rpois(n, lambda = muU0)
        } else if (FamilyU == "Gamma") {
          shape <- gamma.shape(model.u)$alpha
          PredictU1 <- rgamma(n, shape = shape, scale = muU1/shape)
          PredictU0 <- rgamma(n, shape = shape, scale = muU0/shape)
          rm(shape)
        } else if (FamilyU == "binomial") {
          PredictU1 <- rbinom(n, size = 1, prob = muU1)
          PredictU0 <- rbinom(n, size = 1, prob = muU0)
        } else if (FamilyU == "gaussian") {
          sigma <- sqrt(summary(model.u)$dispersion)
          error <- rnorm(n, mean = 0, sd = sigma)
          PredictU1 <- muU1 + error
          PredictU0 <- muU0 + error
          rm(sigma, error)
        } else if (FamilyU == "inverse.gaussian") {
          disp <- summary(model.u)$dispersion
          PredictU1 <- SuppDists::rinvGauss(n, nu = muU1, lambda = 1/disp)
          PredictU0 <- SuppDists::rinvGauss(n, nu = muU0, lambda = 1/disp)
          rm(disp)
        } else {
          stop("unsupported glm family")
        }
        rm(muU1, muU0)
        ### Case II-b: Ordered Confounder
      } else if (isOrdered.u) {
        probs_u1 <- predict(new.fit.U, type = "probs", newdata = pred.data.t)
        probs_u0 <- predict(new.fit.U, type = "probs", newdata = pred.data.c)
        draws_u1 <- matrix(NA, n, u)
        draws_u0 <- matrix(NA, n, u)

        for (ii in 1:n) {
          draws_u1[ii, ] <- t(rmultinom(1, 1, prob = probs_u1[ii, ]))
          draws_u0[ii, ] <- t(rmultinom(1, 1, prob = probs_u0[ii, ]))
        }
        PredictU1 <- apply(draws_u1, 1, indexmax)
        PredictU0 <- apply(draws_u0, 1, indexmax)
        rm(probs_u1, probs_u0, draws_u1, draws_u0)
        ### Case II-c: Quantile Regression for Confounder
      } else if (isRq.u) {
        # Use inverse transform sampling to predict U
        call.new <- new.fit.U$call
        call.new$tau <- runif(n)
        newfits <- eval.parent(call.new)
        tt <- delete.response(terms(new.fit.U))
        u.t <- model.frame(tt, pred.data.t, xlev = new.fit.U$xlevels)
        u.c <- model.frame(tt, pred.data.c, xlev = new.fit.U$xlevels)
        X.t <- model.matrix(tt, u.t, contrasts = new.fit.U$contrasts)
        X.c <- model.matrix(tt, u.c, contrasts = new.fit.U$contrasts)
        rm(tt, u.t, u.c)
        PredictU1 <- rowSums(X.t * t(newfits$coefficients))
        PredictU0 <- rowSums(X.c * t(newfits$coefficients))
        rm(newfits, X.t, X.c)

        ### Case II-d: Linear
      } else if (isLm.u) {
        sigma <- summary(new.fit.U)$sigma
        error <- rnorm(n, mean = 0, sd = sigma)
        PredictU1 <- predict(new.fit.U, type = "response",
                             newdata = pred.data.t) + error
        PredictU0 <- predict(new.fit.U, type = "response",
                             newdata = pred.data.c) + error
        rm(sigma, error)

        ### Case I-2-e: Survreg
      } else if (isSurvreg.u) {
        dd <- survreg.distributions[[new.fit.U$dist]]
        if (is.null(dd$itrans)) {
          itrans <- function(x) x
        } else {
          itrans <- dd$itrans
        }
        dname <- dd$dist
        if (is.null(dname)) {
          dname <- new.fit.U$dist
        }
        scale <- new.fit.U$scale
        lpU1 <- predict(new.fit.U, newdata = pred.data.t, type = "linear")
        lpU0 <- predict(new.fit.U, newdata = pred.data.c, type = "linear")
        error <- switch(dname,
                        extreme = log(rweibull(n, shape = 1, scale = 1)),
                        gaussian = rnorm(n),
                        logistic = rlogis(n),
                        t = rt(n, df = dd$parms))
        PredictU1 <- as.numeric(itrans(lpU1 + scale * error))
        PredictU0 <- as.numeric(itrans(lpU0 + scale * error))
        rm(dd, dname, scale, lpU1, lpU0, error)

      } else {
        stop("confounder model is not yet implemented")
      }
      #####################################
      # Mediator Predictions
      #####################################
      pred.data.t <- pred.data.c <- m.data

      if (isFactorT) {
        pred.data.t[, treat] <- factor(cat.1, levels = t.levels)
        pred.data.c[, treat] <- factor(cat.0, levels = t.levels)
      } else {
        pred.data.t[, treat] <- cat.1
        pred.data.c[, treat] <- cat.0
      }

      if (!is.null(covariates)) {
        for (p in 1:length(covariates)) {
          vl <- names(covariates[p])
          x <- ifelse(is.factor(pred.data.t[, vl]),
                      factor(covariates[[p]], levels = levels(m.data[, vl])),
                      covariates[[p]])
          pred.data.t[, vl] <- pred.data.c[, vl] <- x
        }
      }
      if (isFactorU) {
        pred.data.t[, confounder] <- factor(PredictU1, levels = 1:u, labels = u.levels)
        pred.data.c[, confounder] <- factor(PredictU0, levels = 1:u, labels = u.levels)
      } else {
        pred.data.t[, confounder] <- PredictU1
        pred.data.c[, confounder] <- PredictU0
      }

      ### Case II-a: GLM Mediator (including GAMs)
      if (isGlm.m) {

        muM1 <- predict(new.fit.M, type = "response", newdata = pred.data.t)
        muM0 <- predict(new.fit.M, type = "response", newdata = pred.data.c)

        if (FamilyM == "poisson") {
          PredictM1 <- rpois(n, lambda = muM1)
          PredictM0 <- rpois(n, lambda = muM0)
        } else if (FamilyM == "Gamma") {
          shape <- gamma.shape(model.m)$alpha
          PredictM1 <- rgamma(n, shape = shape, scale = muM1/shape)
          PredictM0 <- rgamma(n, shape = shape, scale = muM0/shape)
          rm(shape)
        } else if (FamilyM == "binomial") {
          PredictM1 <- rbinom(n, size = 1, prob = muM1)
          PredictM0 <- rbinom(n, size = 1, prob = muM0)
        } else if (FamilyM == "gaussian") {
          sigma <- sqrt(summary(model.m)$dispersion)
          error <- rnorm(n, mean = 0, sd = sigma)
          PredictM1 <- muM1 + error
          PredictM0 <- muM0 + error
          rm(sigma, error)
        } else if (FamilyM == "inverse.gaussian") {
          disp <- summary(model.m)$dispersion
          PredictM1 <- SuppDists::rinvGauss(n, nu = muM1, lambda = 1/disp)
          PredictM0 <- SuppDists::rinvGauss(n, nu = muM0, lambda = 1/disp)
          rm(disp)
        } else {
          stop("unsupported glm family")
        }
        rm(muM1, muM0)
        ### Case II-b: Ordered Mediator
      } else if (isOrdered.m) {
        probs_m1 <- predict(new.fit.M, type = "probs", newdata = pred.data.t)
        probs_m0 <- predict(new.fit.M, type = "probs", newdata = pred.data.c)
        draws_m1 <- matrix(NA, n, m)
        draws_m0 <- matrix(NA, n, m)

        for (ii in 1:n) {
          draws_m1[ii, ] <- t(rmultinom(1, 1, prob = probs_m1[ii, ]))
          draws_m0[ii, ] <- t(rmultinom(1, 1, prob = probs_m0[ii, ]))
        }
        PredictM1 <- apply(draws_m1, 1, indexmax)
        PredictM0 <- apply(draws_m0, 1, indexmax)
        rm(probs_m1, probs_m0, draws_m1, draws_m0)
        ### Case II-c: Quantile Regression for Mediator
      } else if (isRq.m) {
        # Use inverse transform sampling to predict M
        call.new <- new.fit.M$call
        call.new$tau <- runif(n)
        newfits <- eval.parent(call.new)
        tt <- delete.response(terms(new.fit.M))
        m.t <- model.frame(tt, pred.data.t, xlev = new.fit.M$xlevels)
        m.c <- model.frame(tt, pred.data.c, xlev = new.fit.M$xlevels)
        X.t <- model.matrix(tt, m.t, contrasts = new.fit.M$contrasts)
        X.c <- model.matrix(tt, m.c, contrasts = new.fit.M$contrasts)
        rm(call.new, tt, m.t, m.c)
        PredictM1 <- rowSums(X.t * t(newfits$coefficients))
        PredictM0 <- rowSums(X.c * t(newfits$coefficients))
        rm(newfits, X.t, X.c)

        ### Case II-d: Linear
      } else if (isLm.m) {
        sigma <- summary(new.fit.M)$sigma
        error <- rnorm(n, mean = 0, sd = sigma)
        PredictM1 <- predict(new.fit.M, type = "response",
                             newdata = pred.data.t) + error
        PredictM0 <- predict(new.fit.M, type = "response",
                             newdata = pred.data.c) + error
        rm(sigma, error)

        ### Case I-2-e: Survreg
      } else if (isSurvreg.m) {
        dd <- survreg.distributions[[new.fit.M$dist]]
        if (is.null(dd$itrans)) {
          itrans <- function(x) x
        } else {
          itrans <- dd$itrans
        }
        dname <- dd$dist
        if (is.null(dname)) {
          dname <- new.fit.M$dist
        }
        scale <- new.fit.M$scale
        lpM1 <- predict(new.fit.M, newdata = pred.data.t, type = "linear")
        lpM0 <- predict(new.fit.M, newdata = pred.data.c, type = "linear")
        error <- switch(dname,
                        extreme = log(rweibull(n, shape = 1, scale = 1)),
                        gaussian = rnorm(n),
                        logistic = rlogis(n),
                        t = rt(n, df = dd$parms))
        PredictM1 <- as.numeric(itrans(lpM1 + scale * error))
        PredictM0 <- as.numeric(itrans(lpM0 + scale * error))
        rm(dd, dname, scale, lpM1, lpM0, error)

      } else {
        stop("mediator model is not yet implemented")
      }

      #####################################
      #  Outcome Predictions
      #####################################
      effects.tmp <- array(NA, dim = c(n, n.ycat, 4))
      for (e in 1:4) {
        tt <- switch(e, c(1, 1, 1, 0), c(0, 0, 1, 0), c(1, 0, 1, 1), c(1, 0, 0, 0))
        pred.data.t <- pred.data.c <- y.data

        if (!is.null(covariates)) {
          for (p in 1:length(covariates)) {
            vl <- names(covariates[p])
            x <- ifelse(is.factor(pred.data.t[, vl]),
                        factor(covariates[[p]], levels = levels(y.data[, vl])),
                        covariates[[p]])
            pred.data.t[, vl] <- pred.data.c[, vl] <- x
          }
        }

        # Set treatment values
        cat.t <- ifelse(tt[1], cat.1, cat.0)
        cat.c <- ifelse(tt[2], cat.1, cat.0)
        cat.t.ctrl <- ifelse(tt[1], cat.0, cat.1)
        cat.c.ctrl <- ifelse(tt[2], cat.0, cat.1)
        if (isFactorT) {
          pred.data.t[, treat] <- factor(cat.t, levels = t.levels)
          pred.data.c[, treat] <- factor(cat.c, levels = t.levels)
          if (!is.null(control)) {
            pred.data.t[, control] <- factor(cat.t.ctrl, levels = t.levels)
            pred.data.c[, control] <- factor(cat.c.ctrl, levels = t.levels)
          }
        } else {
          pred.data.t[, treat] <- cat.t
          pred.data.c[, treat] <- cat.c
          if (!is.null(control)) {
            pred.data.t[, control] <- cat.t.ctrl
            pred.data.c[, control] <- cat.c.ctrl
          }
        }

        # Set confounder values
        PredictUt <- PredictU1 * tt[1] + PredictU0 * (1 - tt[1])
        PredictUc <- PredictU1 * tt[2] + PredictU0 * (1 - tt[2])

        if (isFactorU) {
          pred.data.t[, confounder] <- factor(PredictUt, levels = 1:u, labels = u.levels)
          pred.data.c[, confounder] <- factor(PredictUc, levels = 1:u, labels = u.levels)
        } else {
          pred.data.t[, confounder] <- PredictUt
          pred.data.c[, confounder] <- PredictUc
        }
        # set mediator values
        PredictMt <- PredictM1 * tt[3] + PredictM0 * (1 - tt[3])
        PredictMc <- PredictM1 * tt[4] + PredictM0 * (1 - tt[4])
        if (isFactorM) {
          pred.data.t[, mediator] <- factor(PredictMt, levels = 1:m, labels = m.levels)
          pred.data.c[, mediator] <- factor(PredictMc, levels = 1:m, labels = m.levels)
        } else {
          pred.data.t[, mediator] <- PredictMt
          pred.data.c[, mediator] <- PredictMc
        }

        probs_p1 <- predict(new.fit.Y, newdata = pred.data.t, type = "probs")
        probs_p0 <- predict(new.fit.Y, newdata = pred.data.c, type = "probs")
        effects.tmp[, , e] <- probs_p1 - probs_p0
        rm(pred.data.t, pred.data.c, probs_p1, probs_p0)
      }

      # Compute all QoIs
      if (b == sims + 1) {
        d1 <- apply(effects.tmp[, , 1], 2, weighted.mean, w = weights)
        d0 <- apply(effects.tmp[, , 2], 2, weighted.mean, w = weights)
        z1 <- apply(effects.tmp[, , 3], 2, weighted.mean, w = weights)
        z0 <- apply(effects.tmp[, , 4], 2, weighted.mean, w = weights)
      } else {
        delta.1[b, ] <- apply(effects.tmp[, , 1], 2, weighted.mean, w = weights)
        delta.0[b, ] <- apply(effects.tmp[, , 2], 2, weighted.mean, w = weights)
        zeta.1[b, ] <- apply(effects.tmp[, , 3], 2, weighted.mean, w = weights)
        zeta.0[b, ] <- apply(effects.tmp[, , 4], 2, weighted.mean, w = weights)
      }

    }  # Bootstrap loop ends

    tau.coef <- (d1 + d0 + z1 + z0)/2
    tau <- (zeta.1 + zeta.0 + delta.0 + delta.1)/2

    ########################################################################
    ## Compute Outputs and Put Them Together
    ########################################################################
    low <- (1 - conf.level)/2
    high <- 1 - low
    d0.ci <- apply(delta.0, 2, quantile, c(low, high))
    d1.ci <- apply(delta.1, 2, quantile, c(low, high))
    tau.ci <- apply(tau, 2, quantile, c(low, high))
    z1.ci <- apply(zeta.1, 2, quantile, c(low, high))
    z0.ci <- apply(zeta.0, 2, quantile, c(low, high))

    # p-values
    d0.p <- 2 * apply(delta.0, 2,
                      function(x) sum(sign(x) != sign(median(x)))/sims)
    d1.p <- 2 * apply(delta.1, 2,
                      function(x) sum(sign(x) != sign(median(x)))/sims)
    z0.p <- 2 * apply(zeta.0, 2,
                      function(x) sum(sign(x) != sign(median(x)))/sims)
    z1.p <- 2 * apply(zeta.1, 2,
                      function(x) sum(sign(x) != sign(median(x)))/sims)
    tau.p <- 2 * apply(tau, 2,
                       function(x) sum(sign(x) != sign(median(x)))/sims)

    # Detect whether models include T-M interaction
    INT <- paste(treat, confounder, sep = ":") %in% attr(model.y$terms,
                                                         "term.labels") |
      paste(confounder, treat, sep = ":") %in% attr(model.y$terms,
                                                    "term.labels")

    if (long) {
      out <- list(d0 = d0, d1 = d1, d0.ci = d0.ci, d1.ci = d1.ci,
                  d0.p = d0.p, d1.p = d1.p,
                  d0.sims = delta.0, d1.sims = delta.1,
                  tau.coef = tau.coef, tau.ci = tau.ci, tau.p = tau.p,
                  z0 = z0, z1 = z1, z0.ci = z0.ci, z1.ci = z1.ci,
                  z0.p = z0.p, z1.p = z1.p,
                  z1.sims = zeta.1, z0.sims = zeta.0, tau.sims = tau,
                  boot = boot, treat = treat, confounder = confounder,
                  covariates = covariates,
                  INT = INT, conf.level = conf.level,
                  model.y = model.y, model.u = model.u,
                  control.value = control.value, treat.value = treat.value, nobs = n, sims = sims)
    } else {
      out <- list(d0 = d0, d1 = d1, d0.ci = d0.ci, d1.ci = d1.ci,
                  d0.p = d0.p, d1.p = d1.p,
                  tau.coef = tau.coef, tau.ci = tau.ci, tau.p = tau.p,
                  z0 = z0, z1 = z1, z0.ci = z0.ci, z1.ci = z1.ci,
                  z0.p = z0.p, z1.p = z1.p,
                  boot = boot, treat = treat, confounder = confounder,
                  covariates = covariates,
                  INT = INT, conf.level = conf.level,
                  model.y = model.y, model.u = model.u,
                  control.value = control.value, treat.value = treat.value, nobs = n, sims = sims)
    }
    class(out) <- "mediate.order"
    out
  }
}
