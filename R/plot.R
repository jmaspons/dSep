## Plots ----

#' @export
plot<- graph::plot

#' @rdname dSep
#'
#' @param x a dSep object.
#' @param plotCoef if \code{TRUE} plot the path coefficients.
#' @param plotdSep if \code{TRUE} plot the d-separation paths. By default in red.
#' @param legend if \code{TRUE} add a legend.
#' @param lty a vector of 2 elements with the line type of significant and non significant paths.
#' @param color a vector of 2 elements with the color of the causal paths and the color of the d-separation paths.
#' @param alpha significance level.
#' @export
#' @S3method plot dSep
plot.dSep<- function(x, y, plotCoef=TRUE, plotdSep=FALSE, legend=TRUE, lty=c(signif=1, nonSignif=2), color=c(causal="black", dsep="red"), alpha=0.05, ...){
  if (!requireNamespace("Rgraphviz")){
    stop("You need to install Rgraphviz to plot the results:\n",
         "\tsource('https://bioconductor.org/biocLite.R'\n",
         "\tbiocLite('Rgraphviz')")
  }

  if (!plotCoef){
    g<- x$graph
  }else if (!is.null(x$pathCoefficients) & plotCoef){
    g<- x$pathCoefficients$graph
  }else if (is.null(x$pathCoefficients) & plotCoef){
    g<- x$graph
    warning("dSep object doesn't have path coefficients. Run the models again with the parameter pathCoef=TRUE.")
  }

  if (plotdSep){
    g<- lapply(seq_along(g), function(i){
      merge(g[[i]], graph.conditionalIndependences(x$conditionalIndependences[[i]]))
    })
    names(g)<- names(x$graph)
  }

  parOri<- graphics::par(no.readonly=TRUE)
  modelName<- names(g)
  nG<- length(g)
  nRows<- floor(sqrt(nG))
  nCols<- ceiling(sqrt(nG))

  graphics::par(mfrow=c(nRows, nCols))

  out<- lapply(1:nG, function(i){
    edgeAttributes<- edgeAttrs(g[[i]], alpha=alpha, lty=lty, color=color)
    title<- paste0(modelName[i], "\np-value=", round(x$res$p.value[i], 3), "\tCICc=", round(x$res$CICc[i], 3))
    plot(x=g[[i]], main=title, edgeAttrs=edgeAttributes, ...)
  })
  names(out)<- names(g)

  if (legend){
    legendText<- c("Significant", "Non significant")
    ncols<- 1 ## remove parameter not used
    if (plotdSep){
      # ncols<- 2 ## remove parameter not used
      legendText<- c(legendText, "d-separation")
      lty<- c(lty, lty[1])
      col<- c(rep(color[1], 2), color[2])
    }else{
      col<- rep(color[1], 2)
    }
    legend("bottomright", legendText, lty=lty, col=col, ncol=ncols, title=as.expression(bquote( alpha == .(alpha) )))
  }
  graphics::par(parOri)

  invisible(out)
}


#' @rdname pathCoef
#'
#' @param x a \code{pathCoef} object.
#' @param plotCoef if \code{TRUE} plot the path coefficients.
#' @param legend if \code{TRUE} add a legend.
#' @param lty a vector of 2 elements with the line type of significant and non significant paths.
#' @param alpha significance level.
#' @export
#' @S3method plot pathCoef
plot.pathCoef<- function(x, y, plotCoef=TRUE, legend=TRUE, lty=c(signif=1, nonSignif=2), alpha=0.05, ...){
  if (!requireNamespace("Rgraphviz")){
    stop("You need to install Rgraphviz to plot the results:\n",
         "\tsource('https://bioconductor.org/biocLite.R'\n",
         "\tbiocLite('Rgraphviz')")
  }
  g<- x$graph
  parOri<- graphics::par(no.readonly=TRUE)
  modelName<- names(g)
  nG<- length(g)
  nRows<- floor(sqrt(nG))
  nCols<- ceiling(sqrt(nG))

  graphics::par(mfrow=c(nRows, nCols))
  out<- lapply(1:nG, function(i){
    edgeAttributes<- edgeAttrs(g[[i]], alpha=alpha, lty=lty)
    if (!plotCoef) edgeAttributes$label<- NULL
    plot(x=g[[i]], main=modelName[i], edgeAttrs=edgeAttributes, ...) # TODO: Rgraphviz::plot?
  })
  names(out)<- names(g)

  if (legend){
    legend("bottomright", c("Significant", "Non significant"), lty=lty, title=as.expression(bquote( alpha == .(alpha))))
  }

  graphics::par(parOri)

  invisible(out)
}


## Helpers ----
## TODO: think API
# getDefaultAttrs()
#' Edge attributes for ploting dSep and pathCoef classes
#'
#' @param g a \linkS4class{graph}.
#' @param alpha significance level.
#' @param lty a vector of 2 elements with the line type of significant and non significant paths.
#' @param color a vector of 2 elements with the color of the causal edges and the color of the d-separation edges.
#' @param ... other valid edge attributes.
#' @importFrom graph edgeNames edgeWeights
edgeAttrs<- function(g, alpha=0.05, lty=1:2, color=1:2, ...){
  eAttrs<- list()
  attrName<- names(edgeDataDefaults(g))

  ## Edge type: significant p-value
  if ("p.value" %in% attrName){
    ePvalue<- as.numeric(unlist(graph::edgeWeights(g, attr="p.value", type.checker=is.numeric)))
    ePvalue<- ePvalue[setdiff(seq(along=ePvalue), Rgraphviz::removedEdges(g))]
    eSignif<- ePvalue < alpha
    eLty<- eSignif
    eLty[eSignif]<- lty[1]
    eLty[!eSignif | is.na(eSignif)]<- lty[2]
    names(eLty)<- graph::edgeNames(g)
    eAttrs<- c(eAttrs, list(lty=eLty))
  }

  ## Labels
  if ("coefficients" %in% attrName){
    eCoef<- as.numeric(unlist(graph::edgeWeights(g, attr="coefficients", type.checker=is.numeric)))
    eCoef<- eCoef[setdiff(seq(along=eCoef), Rgraphviz::removedEdges(g))]
    eLab<- paste0(round(eCoef, 2), ifelse(eSignif, " *", ""))
    eLab<- gsub("NANA", "NA", eLab)
    names(eLab)<- graph::edgeNames(g)
    eAttrs<- c(eAttrs, list(label=eLab))
  }

  ## Edge width
  # ?? eAttrs$penwidth
  # ((to_max - to_min) * (self - from_min)) / (from_max - from_min) + to_min
  # TODO: scale in the range [1-10]??
  # eLwd<- unlist(graph::edgeWeights(g, attr="coefficients", type.checker=is.numeric))
  # eLwd<- eLwd[setdiff(seq(along=eLwd), Rgraphviz::removedEdges(g))]

  ## Edge color
  if ("graph" %in% attrName){
    if (is.numeric(color)) color<- grDevices::palette()[color]
    eGraph<- as.character(unlist(graph::edgeWeights(g, attr="graph", type.checker=is.character)))
    eGraph<- eGraph[setdiff(seq(along=eGraph), Rgraphviz::removedEdges(g))]
    eColor<- eGraph
    eColor[eColor %in% "x"]<- color[1]
    eColor[eColor %in% "y"]<- color[2]
    names(eColor)<- graph::edgeNames(g)

    ## Remove label for d-separated edges (has graph attribute == g2)
    if ("label" %in% names(eAttrs)){
      eAttrs$label[eGraph %in% "y"]<- ""
    }

    eAttrs<- c(eAttrs, list(color=eColor))
  }

  c(eAttrs, list(...))
}


#' Merge graphs
#'
#' @param x \linkS4class{graph}
#' @param y \linkS4class{graph}
#' @details If the same edge attribute is present in both, the value from the second one is omited.
#' TODO: merge nodeData and graphData
#' @importFrom graph edgeData edgeDataDefaults edgeNames join
#' @export
merge.graph<- function(x, y, ...){
  ## Topological union
  xy<- graph::join(x, y)

  ## Set edgeDataDefaults (union of x and y)
  eAttrG1<- graph::edgeDataDefaults(x)
  eAttrG2<- graph::edgeDataDefaults(y)
  eAttrL<- c(eAttrG1, eAttrG2) ## Use a list to keep NA mode (data.frame unify them to character if any)
  eAttr<- unique(data.frame(name=names(eAttrL), value=unlist(eAttrL), stringsAsFactors=FALSE))
  eAttr$x<- sapply(eAttr$name, function(n) n %in% names(eAttrG1))
  eAttr$y<- sapply(eAttr$name, function(n) n %in% names(eAttrG2))

  ## TODO: handle cases where x and y have the same attribute with different default value
  ## By default the value from x is used
  for (i in 1:nrow(eAttr)){
    graph::edgeDataDefaults(xy, eAttr$name[i])<- eAttrL[[eAttr$name[i]]]
  }

  # Attribute to save the original graph of the edges
  graph::edgeDataDefaults(xy, "graph")<- ""


  ## Set edgeData
  edg<- graph::edgeNames(xy) # tail~head
  edg<- data.frame(do.call(rbind, strsplit(edg, "~")), edg, stringsAsFactors=FALSE)
  colnames(edg)<- c("tail", "head", "name")
  edgName1<- graph::edgeNames(x) # tail~head
  edgName2<- graph::edgeNames(y) # tail~head

  if (length(dupEdges<- intersect(edgName1, edgName2)) > 0)
    warning("Edge ", dupEdges, " present in both graphs. If the same edge attribute is present in both, the value from the second one is omited.")

  for (i in 1:nrow(edg)){ # Edge loop
    to<- edg$head[i]
    from<- edg$tail[i]
    oriGraph<- character()
    if (edg$name[i] %in% edgName1) oriGraph<- c(oriGraph, "x")
    if (edg$name[i] %in% edgName2) oriGraph<- c(oriGraph, "y")
    graph::edgeData(xy, from=from, to=to, attr="graph")<- paste(oriGraph, collapse="_")

    for (j in 1:nrow(eAttr)){ # Attribute loop
      if (edg$name[i] %in% edgName1 & eAttr$x[j]){
        val<- graph::edgeData(x, from=from, to=to, attr=eAttr$name[j])
        graph::edgeData(xy, from=from, to=to, attr=eAttr$name[j])<- val
      }else if (edg$name[i] %in% edgName2 & eAttr$y[j]){      ## If edgeName exist in x omit y
        val<- graph::edgeData(y, from=from, to=to, attr=eAttr$name[j])
        graph::edgeData(xy, from=from, to=to, attr=eAttr$name[j])<- val
      }
    }
  }
  # edgeData(xy); edgeData(xy, attr="label")

  return(xy)
}
