# Functions for path analysis using D-separation method
# Gonzalez-Voyer, Alejandro, and Achaz Von Hardenberg. 2014. “An Introduction to Phylogenetic Path Analysis.” In Modern Phylogenetic Comparative Methods and Their Application in Evolutionary Biology, edited by László Zsolt Garamszegi, 29. Berlin, Heidelberg: Springer Berlin Heidelberg. doi:10.1007/978-3-662-43550-2.
# http://www.mpcm-evolution.org/practice/online-practical-material-chapter-8/chapter-8-2-step-step-guide-phylogenetic-path-analysis-using-d-sep-method-rhinograds-example

## p-values and coef ----
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


## d-separation ----

#' Calculate the number of conditional independences
#'
#' @param V number of vertices
#' @param A number of arrows
#' @export
condNum<- function(V, A) {
  (factorial(V)/(2 * factorial(V - 2))) - A
}

# C statistic
#
# @param x vector with the relevant p-values
C<- function(x){
  return (-2 * sum(log(x)))
}

# Goodness of the graph model
#
# @param C value of the C-statistics
# @param q number of parameters
# @param n sample size
CICc <- function(C, q, n){
  C + 2 * q * (n / (n - 1 - q))
}


#' Conditional independences from a directed acyclid graphs (DAG)
#'
#' @examples
#' condIndep(g<- dag(~BV:season:interCV + FS:season, forceCheck=TRUE))
#' @import graph
#' @export
condIndep<- function(x, orderResponse){
  if (!inherits(x, "graph")) stop("x must be an object inheriting from graph class.")
  if (require(gRbase) & inherits(x, "graphNEL")){
    if (!gRbase::is.DAG(x)) stop("input must be a directed acyclical graph. See ?dag")
  }else{
    message("gRbase package not installed. No checks ensure that the graph is a directed acyclical graph (DAC).\n\tProceed under your own risk.")
  }
  nod<- graph::nodes(x)
  inEdg<- graph::inEdges(x)
  outEdg<- graph::edges(x)
  # edgeL(x)
  # adj(x, 1)

  # Find all pairs of non adjacent nodes. Conditionally independent (d-separated) pairs of nodes
  dSepPair<- lapply(nod, function(x){
    tmp<- nod[!nod %in% c(x, inEdg[[x]], outEdg[[x]])]
    if (length(tmp) < 1) return(data.frame(x=tmp, y=tmp))
    return (data.frame(x, y=tmp))
  })
  # names(dSepPair)<- nod
  dSepPair<- do.call(rbind, dSepPair)
  # Remove duplicates
  dSepPair<- t(apply(dSepPair, 1, sort))
  dSepPair<- as.data.frame(dSepPair, stringsAsFactors=FALSE)
  dSepPair<- dSepPair[!duplicated(dSepPair),]

  # Find all the nodes with an arrow pointing to any conditionally independent nodes. Causal nodes of two d-separated nodes
  dSepNod<- unique(unlist(dSepPair))
  causalNod<- lapply(dSepNod, function(x){
    inEdg[[x]]
  })
  names(causalNod)<- dSepNod
  # Match d-separated nodes with all causal nodes of the pair
  ## CHECK: ensure that dSepPair is a data.frame. Apply simplifies to vectors when there is only one case.
  dSepCausal<- apply(dSepPair, 1, function(x){
    res<- lapply(x, function(y) causalNod[[y]])
    return (sort(unique(unlist(res))))
  })

  conditionalIndependences<- lapply(1:nrow(dSepPair), function(i){
    if (length(dSepCausal) > i) causal<- dSepCausal[[i]] else causal<- character()
    return (list(dSep=dSepPair[i,], causal=causal))
  })

  model<- lapply(conditionalIndependences, function(cI, orderResponse){
    selResponse<- 1
    if (!missing(orderResponse)){
      selResponse<- which.min(match(cI$dSep, orderResponse))
      if (length(selResponse) == 0) selResponse<- 1
    }

    pred<- paste(sort(as.character(c(cI$dSep[-selResponse], cI$causal))), collapse=" + ")
    mod<- as.formula(paste(cI$dSep[selResponse], "~", pred), env=.GlobalEnv)
    indTest<- cI$dSep[,-selResponse]
    return (list(formula=mod, independentVar=indTest))
  })

  res<- list(conditionalIndependences=conditionalIndependences, model=model)
  class(res)<- "conditionalIndependences"
  return(res)
}


## Test models ----

#' Path analysis usding d-separation
#'
#' @name dSep
#'
#' @param x a \linkS4class{graph} or a list of graph objects representing hypothesis about the correlation structure of the variables. Hypothesis must be represended by directed acyclical graphs (DAG).
#' @param FUN a function or a the name of the function to test the conditional independences.
#' Currently tested with \code{\link[stats]{lm}}, \code{\link[stats]{glm}}, \code{\link[nlme]{gls}}, \code{\link[caper]{pgls}}, \code{\link[MCMCglmm]{MCMCglmm}} and \code{\link[brms]{brm}}.
#' @param formulaArg argument name from FUN that accepts the formula parameter.
#' @param nobs sample size of the dataset.
#' @param cl the number of CPU cores or a cluster object to run the models in parallel. Cluster object can be defined with \code{\link[parallel]{makeCluster}}
#' in package \code{parallel} or \code{\link[snow]{makeCluster}} from \code{snow} package.
#' @param pathCoef if \code{TRUE} calculates the path coefficient with \code{\link{pathCoef}}. It is advised to standardize the
#' input data (see \code{\link{scale}}). By doing this, the standardized coefficients represents the relative
#' strength of each causal relationship in the model.
#' @param ... parameters passed to FUN. Parameters must be named following the FUN arguments (e.g. data=data.frame()).
#' @return a dSep object.
#' @author Joan Maspons <\email{j.maspons@@creaf.uab.cat}>
#' @references Gonzalez-Voyer, Alejandro, and Achaz Von Hardenberg. 2014. “An Introduction to Phylogenetic Path Analysis.” In Modern Phylogenetic Comparative Methods and Their Application in Evolutionary Biology, edited by László Zsolt Garamszegi, 29. Berlin, Heidelberg: Springer Berlin Heidelberg.
#' \url{http://www.mpcm-evolution.org/practice/online-practical-material-chapter-8/chapter-8-2-step-step-guide-phylogenetic-path-analysis-using-d-sep-method-rhinograds-example}
#' @examples
#' ## Dummy data
#' g1<- gRbase::dag(~a:c:d + b:d)
#' g2<- gRbase::dag(~a:c:d + b:d:a)
#' g<- list(m1=g1, m2=g2)
#' d<- data.frame(a=rnorm(100), b=rnorm(100), c=rnorm(100), d=rnorm(100))
#'
#' ## Use of different functions
#' m.lm<- dSep(g1, FUN="lm", nobs=nrow(d), data=d)
#' m.glm<- dSep(g, FUN="glm", nobs=nrow(d), data=d)
#'
#' if (require("nlme"))
#'   m.gls<- dSep(g1, FUN=gls, formulaArg="model", nobs=nrow(d), data=d)
#'
#' if (require(caper)){
#'   data(shorebird, package="caper")
#'   shorebird.data[,2:5]<- scale(log(shorebird.data[,2:5]))
#'   shorebird<- caper::comparative.data(shorebird.tree, shorebird.data, 'Species')
#'   g3<- gRbase::dag(~Egg.Mass:M.Mass:F.Mass)
#'   g4<- gRbase::dag(~Egg.Mass:M.Mass:Cl.size + Cl.size:F.Mass:M.Mass)
#'   m.pgls<- dSep(list(mPhy1=g3, mPhy2=g4), FUN=caper::pgls, nobs=nrow(d), cl=3, data=shorebird, lambda='ML')
#' }
#'
#' if (require("brms"))
#'   m.brmfit<- dSep(g1, FUN=brms::brm, nobs=nrow(d), cl=2, pathCoef=FALSE, data=d)
#'
#' if (require("MCMCglmm"))
#'   m.MCMCglmm<- dSep(g1, FUN=MCMCglmm::MCMCglmm, formulaArg="fixed", nobs=nrow(d), cl=2, data=d, verbose=FALSE)
#'
#' ##
#' plot(m.glm)
#' @export
dSep<- function(x, FUN="lm", formulaArg="formula", nobs, cl, pathCoef=TRUE, ...) UseMethod("dSep")

#' @rdname dSep
#' @export
dSep.graph<- function(x, FUN="lm", formulaArg="formula", nobs, cl, pathCoef=TRUE, ...){
  dSep(list(x), FUN=FUN, formulaArg=formulaArg, nobs=nobs, cl=cl, pathCoef=pathCoef, ...)
}

#' @rdname dSep
#' @export
dSep.pathCoef<- function(x, FUN="lm", formulaArg="formula", nobs, cl, pathCoef=TRUE, ...){
  out<- dSep(list(x$graph), FUN=FUN, formulaArg=formulaArg, nobs=nobs, cl=cl, pathCoef=FALSE, ...)
  out$pathCoefficients<- x
  class(out)<- "dSep"

  return(out)
}

# FUN="lm"; formulaArg="formula"; nobs=nrow(d); cl=1; args0<- list(data=d)
#' @rdname dSep
#' @import graph
#' @importFrom caper pgls
#' @export
dSep.list<- function(x, FUN="lm", formulaArg="formula", nobs, cl, pathCoef=TRUE, ...){
  if (is.null(names(x))) names(x)<- seq_along(x)
  q<- sapply(x, function(y) graph::numNodes(y) + graph::numEdges(y))
  FUN<- match.fun(FUN)

  condInd<- lapply(x, condIndep, orderResponse=match.arg(orderResponse))
  formulas<- unique(unlist(lapply(condInd, function(y) lapply(y$model, function(z) z$formula))))

  args0<- list(...)
  args<- lapply(formulas, function(y, args){
    res<- c(y, args)
    names(res)[1]<- formulaArg
    return(res)
  }, args=args0)


  # Run models in parallel if cl parameter is present
  message(length(args), " models to run")
  if (missing(cl)){
    m<- lapply(args, function(y){
      try(do.call(FUN, y))
    })
  }else{
    numericCl<- FALSE
    if (is.numeric(cl)){
      cl0<- cl
      cl<- parallel::makeCluster(cl)
      numericCl<- TRUE
    }
    message("Running on ", length(cl), " cores")
    parallel::clusterExport(cl=cl, "FUN", envir=environment())
    m<- parallel::parLapply(cl=cl, X=args, fun=function(y){
      try(do.call(FUN, y))
    })
    if (numericCl){
      parallel::stopCluster(cl)
      cl<- cl0
    }
  }
  # message("All models done!")

  ## Match models with conditional inferences and store p-values
  pVal<- lapply(m, function(y){
    if (inherits(y, "try-error")) return(NA)
    pValues(y)
  })

  pValCondInd<- lapply(condInd, function(y){ # Graph loop
    sapply(y$model, function(z){ # Conditional independences loop
      sel<- which(as.character(formulas) %in% as.character(list(z$formula)))
      if (all(is.na(pVal))) return(NA)
      pVal[[sel]][z$independentVar]
    })
  })

  for (i in seq_along(condInd)){ # Graph loop
    for (j in seq_along(condInd[[i]]$model)){ # Conditional independences loop
      condInd[[i]]$model[[j]]$p.value<- pValCondInd[[i]][j]
    }
  }


  ## Stats
  Cx<- lapply(pValCondInd, C)
  Cx.pval<- lapply(Cx, function(p){
    if (all(is.na(p))) return(NA)
    1 - pchisq(p, 2 * length(p))
  })

  ## Try to extract sample size from the models if missing nobs parameter
  if (missing(nobs)) nobs<- nobs(m[[which(!sapply(m, inherits, "try-error"))[1]]])

  CICc.x<- lapply(seq_along(Cx), function(i, Cx, q, nobs) CICc(C=Cx[[i]], q=q[[i]], n=nobs), C=Cx, q=q, nobs=nobs)

  res<- data.frame(q=q, C=as.numeric(Cx), p.value=as.numeric(Cx.pval), CICc=as.numeric(CICc.x), row.names=names(x))
  res<- res[order(res$CICc),]
  res$delta<- res$CICc - min(res$CICc)
  res$logLik<- - res$delta / 2
  likelihood<- exp(res$logLik)
  res$weight<- likelihood / sum(likelihood)

  out<- list(graph=x, res=res, models=m, FUN=FUN, args=args0, formulaArg=formulaArg, conditionalIndependences=condInd, nobs=nobs)
  class(out)<- "dSep"

  if (pathCoef){
    message("Estimating path coefficients...")

## TODO: try
# match.call()
# match(pathCoef, expand.dots=TRUE)
## http://www.r-bloggers.com/function-argument-lists-and-missing/
argList<-  as.list(match.call(expand.dots = TRUE)[-1])
# Enforce inclusion of non-optional arguments
argList$x<- NULL
argList$pathCoef<- NULL
# do.call(pathCoef, argList)

    out$pathCoefficients<- pathCoef.dSep(out, cl=cl, ...)
    class(out)<- "dSep"
  }

  return(out)
}


#' Path coefficients
#'
#' @name pathCoef
#'
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
  varsM$label<- paste0(round(varsM$coefficients, 2), ifelse(varsM$p.value < alpha, " *", ""))
  varsM$label<- gsub("NANA", "NA", varsM$label)

  ## Match models with every graph edge and set edgeData
  for (i in seq_along(x)){ # Graph loop
    edgName<- graph::edgeNames(x[[i]]) # tail~head
    tmpVars<- as.data.frame(do.call(rbind, strsplit(edgName, "~")), stringsAsFactors=FALSE)
    colnames(tmpVars)<- c("predictor", "response")
    tmpVars<- merge(varsM, tmpVars)

    graph::edgeDataDefaults(x[[i]], "coefficients")<- NA_real_
    graph::edgeDataDefaults(x[[i]], "p.value")<- NA_real_
    graph::edgeDataDefaults(x[[i]], "label")<- ""

    for (j in 1:nrow(tmpVars)){ # Edge loop
      to<- tmpVars$response[j]
      from<- tmpVars$predictor[j]
      graph::edgeData(x[[i]], from=from, to=to, attr="coefficients")<- tmpVars$coefficients[j]
      graph::edgeData(x[[i]], from=from, to=to, attr="p.value")<- tmpVars$p.value[j]
      graph::edgeData(x[[i]], from=from, to=to, attr="label")<- tmpVars$label[j]
    }
    # edgeData(x[[i]]); edgeData(x[[i]], attr="label")
  }

  out<- list(graph=x, res=varsM, models=m, FUN=FUN, args=args0, formulaArg=formulaArg, alpha=alpha)
  class(out)<- "pathCoef"

  return(out)
}


## Generics ----

as.data.frame.conditionalIndependences<- function(x, ...){
  condIndp<- sapply(x$conditionalIndependences, function(y){
    paste("(", paste(y$dSep, collapse=", "), "){", paste(y$causal, collapse=", "), "}", sep="")
  })
  model<- sapply(x$model, function(y){
    y$formula
  })
  indepVar<- sapply(x$model, function(y){
    y$independentVar
  })

  pVal<- lapply(x$model, function(y){
    y$p.value
  })
  pVal<- unlist(pVal)
  if (is.null(pVal)) pVal<- NA

  res<- data.frame(condIndp=condIndp, formula=as.character(model), indepVar=as.character(indepVar), p.value=pVal)

  return(res)
}

#' @export
print.conditionalIndependences<- function(x, ...){
  out<- dSep:::as.data.frame.conditionalIndependences(x)
  names(out)<- gsub("condIndp", "condIndp*", names(out))
  if (all(is.na(out$p.value))){
    out<- out[,!names(out) %in% "p.value"]
  }

  print(out, ...)
  cat("\t* (d-separated variables){causal variables}\n")

  invisible(x)
}

summary.list.conditionalIndependences<- function(x, ....){
  condInd<- lapply(x, function(y){
    as.data.frame.conditionalIndependences(y)
  })
  condInd<- unique(do.call(rbind, condInd))
  rownames(condInd)<- NULL
  condInd
}

#' @rdname dSep
#' @export
summary.dSep<- function(x, ....){
  conditionalIndependences<- summary.list.conditionalIndependences(x$conditionalIndependences)
  pathCoefficients<- summary.pathCoef(x$pathCoefficients)

  out<- list(res=x$res, FUN=x$FUN, conditionalIndependences=conditionalIndependences, nobs=nobs)
    if (!is.null(x$pathCoefficients)){
    out<- c(out, list(pathCoefficients=pathCoefficients))
  }

  class(out)<- c("summary.dSep")
  out
}

#' @export
print.summary.dSep<- function(x, ...){
  print(x$res, digits=3, ...)

  cat("\nConditional independence statements:\n")
  names(x$conditionalIndependences)<- gsub("condIndp", "condIndp*", names(x$conditionalIndependences))
  print(x$conditionalIndependences, digits=3, ...)
  cat("\t* (d-separated variables){causal variables}\n")

  if (!is.null(x$pathCoefficients)){
    cat("\nPath coefficients:\n")
    print(x$pathCoefficients, ...)
  }
  invisible(x)
}

#' @export
print.dSep<- function(x, ...){
  print(x$res, digits=3, ...)

  cat("\nConditional independence statements:\n")
  condInd<- lapply(x$conditionalIndependences, function(y){
    tmp<- as.data.frame.conditionalIndependences(y)
    names(tmp)<- gsub("condIndp", "condIndp*", names(tmp))
    tmp
  })
  print(condInd, ...)
  cat("\t* (d-separated variables){causal variables}\n")

  if (!is.null(x$pathCoefficients)){
    cat("\nPath coefficients:\n")
    print(x$pathCoefficients, ...)
  }
}

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

## Plots ----

#' @rdname dSep
#' @import graph Rgraphviz
#' @export
plot.dSep<- function(x, y, ...){
  if (!require(Rgraphviz)){
    stop("You need to install Rgraphviz to plot the results:\n",
         "\tsource('https://bioconductor.org/biocLite.R'\n",
         "\tbiocLite('Rgraphviz')")
  }
  if (!is.null(x$pathCoefficients)){
    invisible(plot(x$pathCoefficients, dSep=x, ...))
  }else{
    g<- x$graph
    parOri<- par(no.readonly=TRUE)
    modelName<- names(g)
    nG<- length(g)
    nRows<- floor(sqrt(nG))
    nCols<- ceiling(sqrt(nG))

    par(mfrow=c(nRows, nCols))

    tmp<- lapply(1:nG, function(i){
      plot(x=g[[i]], main=paste0(modelName[i], "\np-value=", round(x$res$p.value[i], 3), "\tCICc=", round(x$res$CICc[i], 3)), ...)
    })

    par(parOri)

    invisible(tmp)
  }
}

#' @rdname pathCoef
#' @import graph Rgraphviz
#' @export
plot.pathCoef<- function(x, y, alpha=0.05, lty=c(signif=1, nonSignif=2), legend=TRUE, dSep, ...){
  if (!require(Rgraphviz)){
    stop("You need to install Rgraphviz:\n",
         "\tsource('https://bioconductor.org/biocLite.R'\n",
         "\tbiocLite('Rgraphviz')")
  }
  g<- x$graph
  parOri<- par(no.readonly=TRUE)
  modelName<- names(g)
  nG<- length(g)
  nRows<- floor(sqrt(nG))
  nCols<- ceiling(sqrt(nG))

  par(mfrow=c(nRows, nCols))
  eAttrs<- lapply(1:nG, function(i, dSep){
    ## Edge type
    ePvalue<- as.numeric(unlist(graph::edgeWeights(g[[i]], attr="p.value", type.checker=is.numeric)))
    ePvalue<- ePvalue[setdiff(seq(along=ePvalue), Rgraphviz::removedEdges(g[[i]]))]
    eSignif<- ePvalue < alpha
    eLty<- eSignif
    eLty[eSignif]<- lty[1]
    eLty[!eSignif | is.na(eSignif)]<- lty[2]

    ## Labels
    eCoef<- as.numeric(unlist(graph::edgeWeights(g[[i]], attr="coefficients", type.checker=is.numeric)))
    eCoef<- eCoef[setdiff(seq(along=eCoef), Rgraphviz::removedEdges(g[[i]]))]
    eLab<- paste0(round(eCoef, 2), ifelse(eSignif, " *", ""))
    eLab<- gsub("NANA", "NA", eLab)

    ## Edge width
    # ?? eAttrs$penwidth
    # ((to_max - to_min) * (self - from_min)) / (from_max - from_min) + to_min
    ## TODDO: scale in the range [1-10]??
    # eLwd<- unlist(graph::edgeWeights(g[[i]], attr="coefficients", type.checker=is.numeric))
    # eLwd<- eLwd[setdiff(seq(along=eLwd), Rgraphviz::removedEdges(g[[i]]))]

    names(eLab)<- names(eLty)<- graph::edgeNames(g[[i]])
    # names(eLwd)<- graph::edgeNames(g[[i]])

    eAttrs<- list(label=eLab, lty=eLty)
    # eAttrs<- list(label=eLab, lty=eLty, lwd=eLwd)


    if (missing(dSep)){
      plot(x=g[[i]], edgeAttrs=eAttrs, main=modelName[i], ...)
    }else{
      mName<- modelName[i]
      plot(x=g[[i]], edgeAttrs=eAttrs, main=paste0(modelName[i], "\np-value=", round(dSep$res[mName, "p.value"], 3), "\tCICc=", round(dSep$res[mName,"CICc"], 3)), ...)
    }
    eAttrs
  }, dSep=dSep)
  names(eAttrs)<- names(g)

  if (legend){
    legend("bottomright", c("Significant", "Non significant"), lty=lty, title=as.expression(bquote( alpha == .(alpha))))
  }

  par(parOri)

  invisible(eAttrs)
}

