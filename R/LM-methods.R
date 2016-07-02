## p-values ----

#' Exctract p-value of the predictors
#'
#' @name pValues
#' @param x the result of a model from lm, glm, gls, pgls, MCMCglmm or brm functions
#' @return A named vector with p-values for each predictor
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
pValues.pgls<- function(x) stats:::coef.default(caper::summary.pgls(x))[,"Pr(>|t|)"]
#' @rdname pValues
#' @export
pValues.MCMCglmm<- function(x) MCMCglmm::summary.MCMCglmm(x)$solutions[,"pMCMC"]
#' @rdname pValues
#' @export
pValues.gls<- function(x) coef(nlme:::summary.gls(x))[,"p-value"]
#' @rdname pValues
#' @export
pValues.brmsfit<- function(x){
  coefs<- brms::as.mcmc(x)
  coefs<- brms::as.mcmc(lapply(coefs, function(y) y[,grep("^b_", colnames(y))]))

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

## scoef ----

#' @importFrom caper coef.pgls
#' @export
scoef<- function(object, ...) UseMethod("scoef")
#' @export
scoef.default<- function(object, ...) coef(object, ...)
#' @export
scoef.MCMCglmm<- function(object, ...) MCMCglmm::posterior.mode(object$Sol)
#' @export
scoef.brmsfit<- function(object, ...) drop(t(brms::fixef(object, ...)))


#' @importFrom caper nobs.pgls
#' @export
nobs.MCMCglmm<- function(object, ...) length(object$error.term)

