#' Path coefficients
#' It is advised to standardize the input data (see \code{\link{scale}}). By doing this, the standardized coefficients represents the relative
#' strength of each causal relationship in the model.
#'
#' @name pathCoef
#'
#' @param x a \linkS4class{graph} or a list of graph objects.
#' @param FUN a function or a the name of the function to test the conditional independences.
#' Currently tested with \code{\link[stats]{lm}}, \code{\link[stats]{glm}}, \code{\link[nlme]{gls}}, \code{\link[caper]{pgls}}, \code{\link[MCMCglmm]{MCMCglmm}} and \code{\link[brms]{brm}}.
#' @param formulaArg argument name from FUN that accepts the formula parameter of the lineal model.
#' @param nobs sample size of the dataset. If omited, the function try to extract from the model.
#' @param cl the number of CPU cores or a cluster object to run the models in parallel. Cluster object can be defined with \code{\link[parallel]{makeCluster}}
#' in package \code{parallel} or \code{\link[snow]{makeCluster}} from \code{snow} package.
#' @param ... parameters passed to FUN. Parameters must be named following the FUN arguments (e.g. data=data.frame()).
#' @return a pathCoef object.
#' @author Joan Maspons <\email{j.maspons@@creaf.uab.cat}>
#' @examples
#' ## Dummy data
#' g1<- gRbase::dag(~a:c:d + b:d)
#' g2<- gRbase::dag(~a:c:d + b:d:a)
#' g<- list(m1=g1, m2=g2)
#' d<- data.frame(a=rnorm(100), b=rnorm(100), c=rnorm(100), d=rnorm(100))
#'
#' p.a<- pathCoef(g1, FUN="lm", nobs=nrow(d), data=d)
#' plot(p.b<- pathCoef(list(m1=g1, m2=g2), FUN="lm", nobs=nrow(d), data=d))
#' @import graph
#' @export
pathCoef<- function(x, FUN="lm", formulaArg="formula", cl, alpha=0.05, ...) UseMethod("pathCoef")

#' @rdname pathCoef
#' @export
pathCoef.dSep<- function(x, ...){
  pathCoef(x$graph, FUN=x$FUN, formulaArg=x$formulaArg, ...)
}

#' @rdname pathCoef
#' @export
pathCoef.graph<- function(x, ...){
  pathCoef(list(x), ...)
}

#' @rdname pathCoef
#' @export
pathCoef.list<- function(x, FUN="lm", formulaArg="formula", cl, alpha=0.05, ...){
  formulasG<- lapply(x, function(y){
    edgName<- graph::edgeNames(y)
    # nod<- graph::nodes(y)
    # inEdg<- graph::inEdges(y)
    # inEdg<- inEdg[sapply(inEdg, length) > 0]
    # resVar<- names(inEdg)

    formulas<- strsplit(edgName, "~")

    formulas<- lapply(formulas, function(z){
      as.formula(paste(z[2], "~", z[1]), env=.GlobalEnv)
    })

    return(formulas)
  })
  names(formulasG)<- names(x)

  formulas<- unique(unlist(formulasG))
  varsM<- data.frame(do.call(rbind, strsplit(as.character(formulas), "\\s*~\\s*")), stringsAsFactors=FALSE)
  colnames(varsM)<- c("response", "predictor")
  varsM$edgeName<- paste0(varsM$predictor, "~", varsM$response)


  args0<- list(...)
  args<- lapply(formulas, function(y, args){
    res<- c(y, args)
    names(res)[1]<- formulaArg
    return(res)
  }, args=args0)


  ## Run models if cl parameter is present
  message(length(args), " models to run")
  if (missing(cl)){
    m<- lapply(args, function(y){
      try(do.call(FUN, y))
    })
  }else{
    numericCl<- FALSE
    if (is.numeric(cl)){
      cl<- parallel::makeCluster(cl)
      numericCl<- TRUE
    }
    message("Running on ", length(cl), " cores")
    parallel::clusterExport(cl=cl, "FUN", envir=environment())
    m<- parallel::parLapply(cl=cl, X=args, fun=function(y){
      try(do.call(FUN, y))
    })
    if (numericCl) parallel::stopCluster(cl)
  }
  # message("All models done!")

  ## Extract coefficients and p-values
  coefM<- lapply(m, function(y){
    if (inherits(y, "try-error")) return(NA)
    scoef(y)[2] # Intercept: [1]
  })
  pValM<- lapply(m, function(y){
    if (inherits(y, "try-error")) return(NA)
    pValues(y)[2] # Intercept: [1]
  })
  names(coefM)<- names(pValM)<- varsM$response

  varsM$coefficients<- unlist(coefM)
  varsM$p.value<- unlist(pValM)

  ## Match models with every graph edge and set edgeData
  for (i in seq_along(x)){ # Graph loop
    edgName<- graph::edgeNames(x[[i]]) # tail~head
    tmpVars<- as.data.frame(do.call(rbind, strsplit(edgName, "~")), stringsAsFactors=FALSE)
    colnames(tmpVars)<- c("predictor", "response")
    tmpVars<- merge(varsM, tmpVars)

    graph::edgeDataDefaults(x[[i]], "coefficients")<- NA_real_
    graph::edgeDataDefaults(x[[i]], "p.value")<- NA_real_

    for (j in 1:nrow(tmpVars)){ # Edge loop
      to<- tmpVars$response[j]
      from<- tmpVars$predictor[j]
      graph::edgeData(x[[i]], from=from, to=to, attr="coefficients")<- tmpVars$coefficients[j]
      graph::edgeData(x[[i]], from=from, to=to, attr="p.value")<- tmpVars$p.value[j]
    }
    # edgeData(x[[i]]); edgeData(x[[i]], attr="p.value")
  }

  out<- list(graph=x, res=varsM, models=m, FUN=FUN, args=args0, formulaArg=formulaArg, alpha=alpha)
  class(out)<- "pathCoef"

  return(out)
}

## Generic methods ----

#' @rdname pathCoef
#' @export
summary.pathCoef<- function(x, ....){
  out<- list(res=x$res, FUN=x$FUN, alpha=x$alpha)
  class(out)<- "summary.pathCoef"
  out
}

#' @export
print.summary.pathCoef<- function(x, ...){
  print(x$res, digits=3, ...)
}

#' @export
print.pathCoef<- function(x, ...){
  print(x$res, digits=3, ...)
}

