# Functions for path analysis using D-separation method
# Gonzalez-Voyer, Alejandro, and Achaz Von Hardenberg. 2014. “An Introduction to Phylogenetic Path Analysis.” In Modern Phylogenetic Comparative Methods and Their Application in Evolutionary Biology, edited by László Zsolt Garamszegi, 29. Berlin, Heidelberg: Springer Berlin Heidelberg. doi:10.1007/978-3-662-43550-2.
# http://www.mpcm-evolution.org/practice/online-practical-material-chapter-8/chapter-8-2-step-step-guide-phylogenetic-path-analysis-using-d-sep-method-rhinograds-example


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

## conditionalIndependences class ----

#' Conditional independences from a directed acyclid graphs (DAG)
#'
#' @param x a \linkS4class{graph}.
#' @param orderResponse character vector with the order of precedence of the variables when choosing a response variable for the d-separation test.
#' @return a \code{conditionalIndependences} object.
#' @examples
#' g<- gRbase::dag(~BV:season:interCV + FS:season, forceCheck=TRUE)
#' condIndep(g)
#' @export
condIndep<- function(x, orderResponse){
  if (!inherits(x, "graph")) stop("x must be an object inheriting from graph class.")
  if (requireNamespace("gRbase") & inherits(x, "graphNEL")){
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
    mod<- stats::as.formula(paste(cI$dSep[selResponse], "~", pred), env=.GlobalEnv)
    indTest<- cI$dSep[,-selResponse]
    return (list(formula=mod, independentVar=indTest))
  })

  res<- list(conditionalIndependences=conditionalIndependences, model=model)
  class(res)<- "conditionalIndependences"
  return(res)
}


graph.conditionalIndependences<- function(x, ...){
  dSepPair<- lapply(x$conditionalIndependences, function(y) y$dSep)
  dSepPair<- do.call(rbind, dSepPair)
  dSepPair$indepVar<- sapply(x$model, function(y) y$independentVar)
  nod<- sort(unique(unlist(dSepPair)))

  edg<- as.list(nod)
  names(edg)<- nod

  edg<- lapply(edg, function(n){
    tail<- dSepPair$indepVar %in% n
    if (any(tail)){
      tmp<- dSepPair[tail,]
      return(setdiff(tmp[[1]], n))
    }else{
      return(character())
    }
  })
  g<- graph::graphNEL(nodes=nod, edgeL=edg, edgemode="directed")

  pVal<- sapply(x$model, function(y) y$p.value)

  ## Add p-value to edges of the d-separated variables
  if (any(!is.null(pVal))){
    dSepPair$p.value<- pVal

    edgName<- graph::edgeNames(g) # tail~head
    tmpVars<- as.data.frame(do.call(rbind, strsplit(edgName, "~")), stringsAsFactors=FALSE)
    colnames(tmpVars)<- c("predictor", "response")
    dSepPair<- merge(dSepPair, tmpVars)

    # graph::edgeDataDefaults(g, "coefficients")<- NA_real_
    graph::edgeDataDefaults(g, "p.value")<- NA_real_

    for (j in 1:nrow(dSepPair)){ # Edge loop
      to<- dSepPair$response[j]
      from<- dSepPair$predictor[j]
      # graph::edgeData(g, from=from, to=to, attr="coefficients")<- as.numeric(dSepPair$coefficients[j])
      graph::edgeData(g, from=from, to=to, attr="p.value")<- as.numeric(dSepPair$p.value[j])
    }
    # graph::edgeData(g); graph::edgeData(g, attr="p.value")
  }


  return (g)
}

## dSep class ----

#' Path analysis usding d-separation
#'
#' @name dSep
#'
#' @param x a \linkS4class{graph} or a list of graph objects representing hypothesis about the correlation structure of the variables. Hypothesis must be represended by directed acyclical graphs (DAG).
#' @param FUN a function or a the name of the function to test the conditional independences.
#' Currently tested with \code{\link[stats]{lm}}, \code{\link[stats]{glm}}, \code{\link[nlme]{gls}}, \code{\link[caper]{pgls}}, \code{\link[MCMCglmm]{MCMCglmm}} and \code{\link[brms]{brm}}.
#' @param formulaArg argument name from FUN that accepts the formula parameter.
#' @param orderResponse parameter passed to \code{\link{condIndep}}.
#' @param n sample size of the dataset. If omited the function tries to guess from the models using \code{\link{nobs}} and custom functions.
#' @param cl the number of CPU cores or a cluster object to run the models in parallel. Cluster object can be defined with \code{\link[parallel]{makeCluster}}
#' in package \code{parallel} or \code{\link[snow]{makeCluster}} from \code{snow} package.
#' @param pathCoef if \code{TRUE} calculates the path coefficient with \code{\link{pathCoef}}. It is advised to standardize the
#' input data (see \code{\link{scale}}). By doing this, the standardized coefficients represents the relative
#' strength of each causal relationship in the model.
#' @param alpha significance level.
#' @param ... parameters passed to FUN. Parameters must be named following the FUN arguments (e.g. data=data.frame()).
#' @return a dSep object.
#' @author Joan Maspons <\email{j.maspons@@creaf.uab.cat}>
#' @references Gonzalez-Voyer, Alejandro, and Achaz Von Hardenberg. 2014. "An Introduction to Phylogenetic Path Analysis." In Modern Phylogenetic Comparative Methods and Their Application in Evolutionary Biology, edited by Laszlo Zsolt Garamszegi, 29. Berlin, Heidelberg: Springer Berlin Heidelberg.
#' \url{http://www.mpcm-evolution.org/practice/online-practical-material-chapter-8/chapter-8-2-step-step-guide-phylogenetic-path-analysis-using-d-sep-method-rhinograds-example}
#' @examples
#' ## Dummy data
#' g1<- gRbase::dag(~a:c:d + b:d)
#' g2<- gRbase::dag(~a:c:d + b:d:a)
#' g<- list(m1=g1, m2=g2)
#' d<- data.frame(a=rnorm(100), b=rnorm(100), c=rnorm(100), d=rnorm(100))
#'
#' ## Use of different functions
#' m.lm<- dSep(g1, FUN="lm", n=nrow(d), data=d)
#' m.glm<- dSep(g, FUN="glm", n=nrow(d), data=d)
#'
#' if (require("nlme"))
#'   m.gls<- dSep(g1, FUN=gls, formulaArg="model", n=nrow(d), data=d)
#'
#' if (require(caper)){
#'   data(shorebird, package="caper")
#'   shorebird.data[,2:5]<- scale(log(shorebird.data[,2:5]))
#'   shorebird<- comparative.data(shorebird.tree, shorebird.data, 'Species')
#'   g3<- gRbase::dag(~Egg.Mass:M.Mass:F.Mass)
#'   g4<- gRbase::dag(~Egg.Mass:M.Mass:Cl.size + Cl.size:F.Mass:M.Mass)
#'   m1.pgls<- dSep(list(mPhy1=g3, mPhy2=g4), FUN=caper::pgls, n=nrow(d), cl=2, data=shorebird, lambda='ML')
#'
#'   rhino.dat <- read.csv("http://mpcm-evolution.org/OPM/Chapter8_OPM/download/rhino.csv")
#'   rhino.tree <- read.tree("http://mpcm-evolution.org/OPM/Chapter8_OPM/download/rhino.tree")
#'   com.dat<- caper::comparative.data(rhino.tree, rhino.dat, SP, vcv=TRUE, vcv.dim=3, warn.dropped=TRUE)
#'   m<- list()
#'   m$h1<- gRbase::dag(~LS:BM + NL:BM + DD:NL + RS:DD)
#'   m$h2<- gRbase::dag(~LS:BM + NL:BM + DD:NL + RS:DD:LS)
#'   m$h3<- gRbase::dag(~LS:BM + NL:BM + DD:NL + RS:NL)
#'   m$h4<- gRbase::dag(~LS:BM + NL:BM + DD:NL + RS:BM:NL)
#'   m$h5<- gRbase::dag(~LS:BM + NL:BM + DD:NL + RS:BM:NL:DD)
#'   m$h6<- gRbase::dag(~LS:BM + NL:BM:RS + DD:NL + RS:BM)
#'   m$h7<- gRbase::dag(~LS:BM + NL:BM:RS + DD:NL + RS:BM:LS)
#'   m$h8<- gRbase::dag(~LS:BM + NL:BM:RS + DD:NL)
#'   m$h9<- gRbase::dag(~LS:BM + NL:BM:RS + DD:NL + RS:LS)
#'
#'   m2.pgls<- dSep(m, FUN=pgls, cl=2, orderResponse=c("BM", "RS", "LS", "NL", "DD"), data=com.dat, lambda="ML")
#' }
#'
#' if (require("brms"))
#'   m.brmfit<- dSep(g1, FUN=brm, n=nrow(d), cl=2, pathCoef=FALSE, data=d)
#'
#' if (require("MCMCglmm"))
#'   m.MCMCglmm<- dSep(g1, FUN=MCMCglmm, formulaArg="fixed", n=nrow(d), cl=2, data=d, verbose=FALSE)
#'
#' ##
#' plot(m.glm)
#' @export
dSep<- function(x, FUN="lm", orderResponse, formulaArg="formula", n, cl, pathCoef=TRUE, ...) UseMethod("dSep")

#' @rdname dSep
#' @export
dSep.graph<- function(x, FUN="lm", orderResponse, formulaArg="formula", n, cl, pathCoef=TRUE, ...){
  dSep(list(x), FUN=FUN, orderResponse=orderResponse, formulaArg=formulaArg, n=n, cl=cl, pathCoef=pathCoef, ...)
}

#' @rdname dSep
#' @export
dSep.pathCoef<- function(x, FUN="lm", orderResponse, formulaArg="formula", n, cl, pathCoef=TRUE, ...){
  out<- dSep(list(x$graph), FUN=FUN, orderResponse=orderResponse, formulaArg=formulaArg, n=n, cl=cl, pathCoef=FALSE, ...)
  out$pathCoefficients<- x
  class(out)<- "dSep"

  return(out)
}

# FUN="lm"; formulaArg="formula"; n=nrow(d); cl=1; args0<- list(data=d)
#' @rdname dSep
#' @export
dSep.list<- function(x, FUN="lm", orderResponse, formulaArg="formula", n, cl, pathCoef=TRUE, ...){
  if (is.null(names(x))) names(x)<- seq_along(x)
  q<- sapply(x, function(y) graph::numNodes(y) + graph::numEdges(y))

  FUN<- match.fun(FUN)

  condInd<- lapply(x, condIndep, orderResponse=orderResponse)
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
    m<- parallel::parLapply(cl=cl, X=args, function(y){
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
      ## TODO: use grep to catch categorical predictors. Use min p-val of all levels?
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
    1 - stats::pchisq(p, 2 * length(p))
  })

  ## Try to extract sample size from the models if n parameter is missing
  if (missing(n)){
    selEg<- which(!sapply(m, inherits, "try-error"))[1]
    n<- ifelse(is.na(selEg), NA, nobs(m[[selEg]]))
  }

  CICc.x<- lapply(seq_along(Cx), function(i, Cx, q, n) CICc(C=Cx[[i]], q=q[[i]], n=n), C=Cx, q=q, n=n)

  res<- data.frame(q=q, C=as.numeric(Cx), p.value=as.numeric(Cx.pval), CICc=as.numeric(CICc.x), row.names=names(x))
  res<- res[order(res$CICc),]
  res$delta<- res$CICc - min(res$CICc)
  res$logLik<- - res$delta / 2
  likelihood<- exp(res$logLik)
  res$weight<- likelihood / sum(likelihood)

  out<- list(graph=x, res=res, models=m, FUN=FUN, args=args0, formulaArg=formulaArg, conditionalIndependences=condInd, n=n)
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


## Generic methods ----

## conditionalIndependences

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
  out<- as.data.frame.conditionalIndependences(x)
  names(out)<- gsub("condIndp", "condIndp*", names(out))
  if (all(is.na(out$p.value))){
    out<- out[,!names(out) %in% "p.value"]
  }

  print(out, ...)
  cat("\t* (d-separated variables){causal variables}\n")

  invisible(x)
}

summary.list.conditionalIndependences<- function(object, ...){
  condInd<- lapply(object, function(x){
    as.data.frame.conditionalIndependences(x)
  })
  condInd<- unique(do.call(rbind, condInd))
  rownames(condInd)<- NULL
  condInd
}


## dSep

#' @rdname dSep
#' @export
summary.dSep<- function(object, ...){
  conditionalIndependences<- summary.list.conditionalIndependences(object$conditionalIndependences)
  pathCoefficients<- summary.pathCoef(object$pathCoefficients)

  out<- list(res=object$res, FUN=object$FUN, conditionalIndependences=conditionalIndependences, n=object$n)
  if (!is.null(object$pathCoefficients)){
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

