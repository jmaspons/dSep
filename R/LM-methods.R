## p-values ----

#' Exctract p-value of the predictors
#'
#' @name pValues
#' @param x the result of a model from lm, glm, gls, pgls, MCMCglmm or brm functions
#' @return A named vector with p-values for each predictor
#' @importFrom stats coef
#' @export
pValues<- function(x) UseMethod("pValues")

#' @rdname pValues
#' @export
pValues.default<- function(x) coef(summary(x))[,4]
#' @rdname pValues
#' @export
pValues.lm<- function(x) coef(summary(x))[,"Pr(>|t|)"]
#' @rdname pValues
#' @export
pValues.pgls<- function(x) coef(caper::summary.pgls(x))[,"Pr(>|t|)"]
#' @rdname pValues
#' @export
pValues.phylolm<- function(x) coef(phylolm::summary.phylolm(x))[,"p.value"]
#' @rdname pValues
#' @export
pValues.phyloglm<- function(x) coef(phylolm::summary.phyloglm(x))[,"p.value"]
#' @rdname pValues
#' @export
pValues.MCMCglmm<- function(x) MCMCglmm::summary.MCMCglmm(x)$solutions[,"pMCMC"]
#' @rdname pValues
#' @export
pValues.gls<- function(x) coef(nlme:::summary.gls(x))[,"p-value"]
#' @rdname pValues
#' @export
pValues.brmsfit<- function(x){
  # Does the model contain posterior samples?
  if (!inherits(x$fit, "stanfit") || !length(x$fit@sim)) return(NA_real_)

  coefs<- brms::as.mcmc(x, pars="^b_")

  # lapply(coefs, colMeans) # post.mean
  # lapply(coefs, coda::HPDinterval) # "l-95% CI", "u-95% CI"
  # lapply(coefs, coda::effectiveSize) # eff.samp

  res<- lapply(coefs, function(y){ # pMCMC
    ## code from MCMCglmm::summary.MCMCglmm
    2 * pmax(0.5 / dim(y)[1],
             pmin( colSums(y > 0) / dim(y)[1], 1 - colSums(y > 0 ) / dim(y)[1]) )
  })
  res<- do.call(rbind, res)
  colnames(res)<- colnames(coefs[[1]])

  return(colMeans(res))
}
#' @rdname pValues
#' @export
pValues.lme <- function(x) summary(x)$tTable[,"p-value"]


## scoef ----

#' @importFrom stats coef
#' @export
scoef<- function(object, ...) UseMethod("scoef")

#' @export
scoef.default<- function(object, ...) coef(object, ...) # phylolm, glm, lm, gls
#' @export
scoef.pgls<- function(object, ...) caper::coef.pgls(object, ...)
#' @export
scoef.MCMCglmm<- function(object, ...) MCMCglmm::posterior.mode(object$Sol)
#' @export
scoef.brmsfit<- function(object, ...){
  # Does the model contain posterior samples?
  if (!inherits(object$fit, "stanfit") || !length(object$fit@sim)) return(NA_real_)
  drop(t(brms::fixef(object, ...)))
}
#' @export
scoef.lme <- function(object, ...) object$coefficients$fixed


## nobs ----

#' @importFrom stats nobs
#' @export
nobs.pgls<- function(object, ...) caper::nobs.pgls(object, ...)
#' @export
nobs.phylolm<- function(object, ...) phylolm::nobs.phylolm(object, ...)
#' @export
nobs.phyloglm<- function(object, ...) phylolm::nobs.phylolm(object, ...)
#' @export
nobs.MCMCglmm<- function(object, ...) length(object$error.term)
#' @export
nobs.lme<- function(object, ...) nlme:::nobs.lme(object, ...)
