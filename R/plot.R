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


