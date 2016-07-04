## Run dSep() ----
# with FUN parameter as lm, glm, gls, pgls, MCMCglmm and brm

g1<- gRbase::dag(~a:c:d + b:d)
g2<- gRbase::dag(~a:c:d + b:d:a)
g<- list(m1=g1, m2=g2)
d<- data.frame(a=rnorm(100), b=rnorm(100), c=rnorm(100), d=rnorm(100))

context("dSep-FUN=lm")

test_that("lm works", {
  d.lm<- dSep(g1, FUN="lm", nobs=nrow(d), data=d)
  expect_is(d.lm, "dSep")
  d.lm<- dSep(g, FUN="lm", nobs=nrow(d), data=d)
  expect_is(d.lm, "dSep")

  d.lm<- dSep(g1, FUN="lm", nobs=nrow(d), cl=2, data=d)
  expect_is(d.lm, "dSep")
  d.lm<- dSep(g, FUN="lm", nobs=nrow(d), cl=2, data=d)
  expect_is(d.lm, "dSep")

  cl<- parallel::makeCluster(2, type="FORK")
  d.lm<- dSep(g1, FUN="lm", nobs=nrow(d), cl=cl, data=d)
  expect_is(d.lm, "dSep")
  d.lm<- dSep(g, FUN="lm", nobs=nrow(d), cl=cl, data=d)
  expect_is(d.lm, "dSep")
  parallel::stopCluster(cl)
})

context("dSep-FUN=glm")

test_that("glm works", {
  d.glm<- dSep(g1, FUN="glm", nobs=nrow(d), data=d)
  expect_is(d.glm, "dSep")
  d.glm<- dSep(g, FUN="glm", nobs=nrow(d), data=d)
  expect_is(d.glm, "dSep")

  d.glm<- dSep(g1, FUN="glm", nobs=nrow(d), cl=2, data=d)
  expect_is(d.glm, "dSep")
  d.glm<- dSep(g, FUN="glm", nobs=nrow(d), cl=2, data=d)
  expect_is(d.glm, "dSep")

  cl<- parallel::makeCluster(2, type="FORK")
  d.glm<- dSep(g1, FUN="glm", nobs=nrow(d), cl=cl, data=d)
  expect_is(d.glm, "dSep")
  d.glm<- dSep(g, FUN="glm", nobs=nrow(d), cl=cl, data=d)
  expect_is(d.glm, "dSep")
  parallel::stopCluster(cl)
})

context("dSep-FUN=gls")

test_that("gls works", {
  skip_if_not_installed("nlme")
  # d.gls<- dSep(g1, FUN=nlme::gls, formulaArg="model", nobs=nrow(d), data=d)
  # expect_is(d.gls, "dSep")
  d.gls<- dSep(g, FUN=nlme::gls, formulaArg="model", nobs=nrow(d), data=d)
  expect_is(d.gls, "dSep")

  # d.gls<- dSep(g1, FUN=nlme::gls, formulaArg="model", nobs=nrow(d), cl=2, data=d)
  # expect_is(d.gls, "dSep")
  # d.gls<- dSep(g, FUN=nlme::gls, formulaArg="model", nobs=nrow(d), cl=2, data=d)
  # expect_is(d.gls, "dSep")

  # cl<- parallel::makeCluster(2, type="FORK")
  # d.gls<- dSep(g1, FUN=nlme::gls, formulaArg="model", nobs=nrow(d), cl=cl, data=d)
  # expect_is(d.gls, "dSep")
  # d.gls<- dSep(g, FUN=nlme::gls, formulaArg="model", nobs=nrow(d), cl=cl, data=d)
  # expect_is(d.gls, "dSep")
  # parallel::stopCluster(cl)
})

context("dSep-FUN=pgls")

test_that("pgls works", {
  skip_if_not_installed("caper")
  data(shorebird, package="caper")
  shorebird.data[,2:5]<- scale(log(shorebird.data[,2:5]))
  shorebird<- caper::comparative.data(phy=shorebird.tree, data=shorebird.data, names.col=Species, )
  g3<- gRbase::dag(~Egg.Mass:M.Mass:F.Mass)
  g4<- gRbase::dag(~Egg.Mass:M.Mass:Cl.size + Cl.size:F.Mass:M.Mass)
  gPhy<- list(mPhy1=g3, mPhy2=g4)

  # d.pgls<- dSep(g3, FUN=caper::pgls, nobs=nrow(shorebird$data), data=shorebird)
  # expect_is(d.pgls, "dSep")
  # d.pgls<- dSep(gPhy, FUN=caper::pgls, nobs=nrow(shorebird$data), data=shorebird)
  # expect_is(d.pgls, "dSep")

  # d.pgls<- dSep(g3, FUN=caper::pgls, nobs=nrow(shorebird$data), cl=2, data=shorebird)
  # expect_is(d.pgls, "dSep")
  # d.pgls<- dSep(gPhy, FUN=caper::pgls, nobs=nrow(shorebird$data), cl=2, data=shorebird)
  # expect_is(d.pgls, "dSep")

  cl<- parallel::makeCluster(2, type="FORK")
  # d.pgls<- dSep(g3, FUN=caper::pgls, nobs=nrow(shorebird$data), cl=cl, data=shorebird)
  # expect_is(d.pgls, "dSep")
  d.pgls<- dSep(gPhy, FUN="pgls", nobs=nrow(shorebird$data), cl=cl, data=shorebird)
  expect_is(d.pgls, "dSep")
  parallel::stopCluster(cl)
})

context("dSep-FUN=brm")

test_that("brm works", {
  skip_if_not_installed("brms")
  library(brms)
  # d.brmfit<- dSep(g1, FUN=brms::brm, nobs=nrow(d), pathCoef=FALSE, data=d)
  # expect_is(d.brmfit, "dSep")
  # d.brmfit<- dSep(g, FUN=brms::brm, nobs=nrow(d), pathCoef=FALSE,  data=d)
  # expect_is(d.brmfit, "dSep")
  #
  # d.brmfit<- dSep(g1, FUN=brms::brm, nobs=nrow(d), cl=2, pathCoef=FALSE, data=d)
  # expect_is(d.brmfit, "dSep")
  # d.brmfit<- dSep(g, FUN=brms::brm, nobs=nrow(d), cl=2, pathCoef=FALSE, data=d)
  # expect_is(d.brmfit, "dSep")

  # cl<- parallel::makeCluster(2, type="FORK")
  # d.brmfit<- dSep(g1, FUN=brms::brm, nobs=nrow(d), cl=cl, pathCoef=FALSE, data=d)
  # expect_is(d.brmfit, "dSep")
  # d.brmfit<- dSep(g, FUN=brms::brm, nobs=nrow(d), cl=cl, pathCoef=FALSE, data=d)
  # expect_is(d.brmfit, "dSep")
  # parallel::stopCluster(cl)
})

context("dSep-FUN=MCMCglmm")

test_that("MCMCglmm works", {
  skip_if_not_installed("MCMCglmm")
  # d.MCMCglmm<- dSep(g1, FUN=MCMCglmm::MCMCglmm, formulaArg="fixed", nobs=nrow(d), data=d, verbose=FALSE)
  # expect_is(d.MCMCglmm, "dSep")
  # d.MCMCglmm<- dSep(g, FUN=MCMCglmm::MCMCglmm, formulaArg="fixed", nobs=nrow(d), data=d, verbose=FALSE)
  # expect_is(d.MCMCglmm, "dSep")
  #
  # d.MCMCglmm<- dSep(g1, FUN=MCMCglmm::MCMCglmm, formulaArg="fixed", nobs=nrow(d), cl=2, data=d, verbose=FALSE)
  # expect_is(d.MCMCglmm, "dSep")
  # d.MCMCglmm<- dSep(g, FUN=MCMCglmm::MCMCglmm, formulaArg="fixed", nobs=nrow(d), cl=2, data=d, verbose=FALSE)
  # expect_is(d.MCMCglmm, "dSep")

  cl<- parallel::makeCluster(2, type="FORK")
  # d.MCMCglmm<- dSep(g1, FUN=MCMCglmm::MCMCglmm, formulaArg="fixed", nobs=nrow(d), cl=cl, data=d, verbose=FALSE)
  # expect_is(d.MCMCglmm, "dSep")
  d.MCMCglmm<- dSep(g, FUN=MCMCglmm::MCMCglmm, formulaArg="fixed", nobs=nrow(d), cl=cl, data=d, verbose=FALSE)
  expect_is(d.MCMCglmm, "dSep")
  parallel::stopCluster(cl)
})

## dSep methods ----
context("Methods for dSep objects")

load(system.file(c("extdata/sampleModels.RData"), package="dSep"), verbose=FALSE)

test_that("print works", {
  lapply(D, function(x) return(x))
  lapply(D, function(x) return(x$pathCoefficients))
})

test_that("summary works", {
  # lapply(D, function(x) summary(x))
  # lapply(D, function(x) summary(x$pathCoefficients))
  # lapply(D, function(x) dSep:::summary.list.conditionalIndependences(x$conditionalIndependences))
  lapply(D, function(x) expect_is(summary(x), "summary.dSep"))
  lapply(D[sapply(D, function(x) !is.null(x$pathCoefficients))], function(x) expect_is(summary(x$pathCoefficients), "summary.pathCoef"))
  lapply(D, function(x) expect_is(dSep:::summary.list.conditionalIndependences(x$conditionalIndependences), "data.frame"))
})

test_that("plot works", {
  skip_if_not_installed("Rgraphviz")

  lapply(D, function(x) ifelse(inherits(x, "dSep"), try(plot(x)), "ERROR"))
  lapply(D, function(x) ifelse(inherits(x$pathCoefficients, "pathCoef"), try(plot(x)), "ERROR"))

  lapply(D, function(x) ifelse(inherits(x, "dSep"), try(plot(x, plotCoef=FALSE, plotdSep=FALSE)), "ERROR"))
  lapply(D, function(x) ifelse(inherits(x, "dSep"), try(plot(x, plotCoef=TRUE, plotdSep=FALSE)), "ERROR"))
  lapply(D, function(x) ifelse(inherits(x, "dSep"), try(plot(x, plotCoef=FALSE, plotdSep=TRUE)), "ERROR"))
  lapply(D, function(x) ifelse(inherits(x, "dSep"), try(plot(x, plotCoef=TRUE, plotdSep=TRUE)), "ERROR"))
})

