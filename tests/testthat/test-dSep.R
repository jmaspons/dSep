## Run dSep() ----
# with FUN parameter as lm, glm, gls, pgls, MCMCglmm and brm

g1<- gRbase::dag(~a:c:d + b:d)
g2<- gRbase::dag(~a:c:d + b:d:a)
g<- list(m1=g1, m2=g2)
d<- data.frame(a=rnorm(100), b=rnorm(100), c=rnorm(100), d=rnorm(100))

if (require(caper) | require(phylolm)){
  rhino.dat <- read.csv("http://mpcm-evolution.org/OPM/Chapter8_OPM/download/rhino.csv")
  rhino.tree <- ape::read.tree("http://mpcm-evolution.org/OPM/Chapter8_OPM/download/rhino.tree")
  com.dat<- caper::comparative.data(rhino.tree, rhino.dat, SP, vcv=TRUE, vcv.dim=3, warn.dropped=TRUE)
  m<- list()
  m$h1<- gRbase::dag(~LS:BM + NL:BM + DD:NL + RS:DD)
  m$h2<- gRbase::dag(~LS:BM + NL:BM + DD:NL + RS:DD:LS)
  # m$h3<- gRbase::dag(~LS:BM + NL:BM + DD:NL + RS:NL)
  # m$h4<- gRbase::dag(~LS:BM + NL:BM + DD:NL + RS:BM:NL)
  # m$h5<- gRbase::dag(~LS:BM + NL:BM + DD:NL + RS:BM:NL:DD)
  # m$h6<- gRbase::dag(~LS:BM + NL:BM:RS + DD:NL + RS:BM)
  # m$h7<- gRbase::dag(~LS:BM + NL:BM:RS + DD:NL + RS:BM:LS)
  # m$h8<- gRbase::dag(~LS:BM + NL:BM:RS + DD:NL)
  # m$h9<- gRbase::dag(~LS:BM + NL:BM:RS + DD:NL + RS:LS)
}

## lm ----
context("dSep-FUN=lm")

test_that("lm works", {
  d.lm<- dSep(g1, FUN="lm", n=nrow(d), data=d)
  expect_is(d.lm, "dSep")
  d.lm<- dSep(g, FUN="lm", n=nrow(d), data=d)
  expect_is(d.lm, "dSep")

  d.lm<- dSep(g1, FUN="lm", n=nrow(d), cl=2, data=d)
  expect_is(d.lm, "dSep")
  d.lm<- dSep(g, FUN="lm", n=nrow(d), cl=2, data=d)
  expect_is(d.lm, "dSep")

  cl<- parallel::makeCluster(2, type="FORK")
  d.lm<- dSep(g1, FUN="lm", n=nrow(d), cl=cl, data=d)
  expect_is(d.lm, "dSep")
  d.lm<- dSep(g, FUN="lm", n=nrow(d), cl=cl, data=d)
  expect_is(d.lm, "dSep")
  parallel::stopCluster(cl)
})

## glm ----
context("dSep-FUN=glm")

test_that("glm works", {
  d.glm<- dSep(g1, FUN="glm", n=nrow(d), data=d)
  expect_is(d.glm, "dSep")
  d.glm<- dSep(g, FUN="glm", n=nrow(d), data=d)
  expect_is(d.glm, "dSep")

  d.glm<- dSep(g1, FUN="glm", n=nrow(d), cl=2, data=d)
  expect_is(d.glm, "dSep")
  d.glm<- dSep(g, FUN="glm", n=nrow(d), cl=2, data=d)
  expect_is(d.glm, "dSep")

  cl<- parallel::makeCluster(2, type="FORK")
  d.glm<- dSep(g1, FUN="glm", n=nrow(d), cl=cl, data=d)
  expect_is(d.glm, "dSep")
  d.glm<- dSep(g, FUN="glm", n=nrow(d), cl=cl, data=d)
  expect_is(d.glm, "dSep")
  parallel::stopCluster(cl)
})

## gls ----
context("dSep-FUN=gls")

test_that("gls works", {
  skip_if_not_installed("nlme")
  # d.gls<- dSep(g1, FUN=nlme::gls, formulaArg="model", n=nrow(d), data=d)
  # expect_is(d.gls, "dSep")
  d.gls<- dSep(g, FUN=nlme::gls, formulaArg="model", n=nrow(d), data=d)
  expect_is(d.gls, "dSep")

  # d.gls<- dSep(g1, FUN=nlme::gls, formulaArg="model", n=nrow(d), cl=2, data=d)
  # expect_is(d.gls, "dSep")
  # d.gls<- dSep(g, FUN=nlme::gls, formulaArg="model", n=nrow(d), cl=2, data=d)
  # expect_is(d.gls, "dSep")

  # cl<- parallel::makeCluster(2, type="FORK")
  # d.gls<- dSep(g1, FUN=nlme::gls, formulaArg="model", n=nrow(d), cl=cl, data=d)
  # expect_is(d.gls, "dSep")
  # d.gls<- dSep(g, FUN=nlme::gls, formulaArg="model", n=nrow(d), cl=cl, data=d)
  # expect_is(d.gls, "dSep")
  # parallel::stopCluster(cl)
})

## pgls ----

context("dSep-FUN=pgls")

test_that("pgls works", {
  skip_if_not_installed("caper")
  # data(shorebird, package="caper")
  # shorebird.data[,2:5]<- scale(log(shorebird.data[,2:5]))
  # shorebird<- caper::comparative.data(phy=shorebird.tree, data=shorebird.data, names.col=Species)
  # g3<- gRbase::dag(~Egg.Mass:M.Mass:F.Mass)
  # g4<- gRbase::dag(~Egg.Mass:M.Mass:Cl.size + Cl.size:F.Mass:M.Mass)
  # gPhy<- list(mPhy1=g3, mPhy2=g4)

  # d.pgls<- dSep(g3, FUN=caper::pgls, n=nrow(shorebird$data), data=shorebird)
  # expect_is(d.pgls, "dSep")
  # d.pgls<- dSep(gPhy, FUN=caper::pgls, n=nrow(shorebird$data), data=shorebird)
  # expect_is(d.pgls, "dSep")

  # d.pgls<- dSep(g3, FUN=caper::pgls, n=nrow(shorebird$data), cl=2, data=shorebird)
  # expect_is(d.pgls, "dSep")
  # d.pgls<- dSep(gPhy, FUN=caper::pgls, n=nrow(shorebird$data), cl=2, data=shorebird)
  # expect_is(d.pgls, "dSep")

  # cl<- parallel::makeCluster(2, type="FORK")
  # d.pgls<- dSep(g3, FUN=caper::pgls, n=nrow(shorebird$data), cl=cl, data=shorebird)
  # expect_is(d.pgls, "dSep")
  # d.pgls<- dSep(gPhy, FUN="pgls", n=nrow(shorebird$data), cl=cl, data=shorebird)
  # expect_is(d.pgls, "dSep")
  # parallel::stopCluster(cl)


  ## Example 2
  # d.pgls<- dSep(m, FUN=caper::pgls, cl=4, orderResponse=c("BM", "RS", "LS", "NL", "DD"), data=com.dat, lambda="ML")
  # expect_is(d.pgls, "dSep")
})

## phylolm ----

context("dSep-FUN=phylolm")

test_that("phylolm works", {
  skip_if_not_installed("phylolm")
  # data(shorebird, package="caper")
  # shorebird.data[,2:5]<- scale(log(shorebird.data[,2:5]))
  # shorebird<- caper::comparative.data(phy=shorebird.tree, data=shorebird.data, names.col=Species)
  # g3<- gRbase::dag(~Egg.Mass:M.Mass:F.Mass)
  # g4<- gRbase::dag(~Egg.Mass:M.Mass:Cl.size + Cl.size:F.Mass:M.Mass)
  # gPhy<- list(mPhy1=g3, mPhy2=g4)

  # d.phylolm<- dSep(g3, FUN=phylolm:phylolm, n=nrow(shorebird$data), data=shorebird)
  # expect_is(d.phylolm, "dSep")
  # d.phylolm<- dSep(gPhy, FUN=phylolm:phylolm, n=nrow(shorebird$data), data=shorebird)
  # expect_is(d.phylolm, "dSep")

  # d.phylolm<- dSep(g3, FUN=phylolm:phylolm, n=nrow(shorebird$data), cl=2, data=shorebird)
  # expect_is(d.phylolm, "dSep")
  # d.phylolm<- dSep(gPhy, FUN=phylolm:phylolm, n=nrow(shorebird$data), cl=2, data=shorebird)
  # expect_is(d.phylolm, "dSep")

  # cl<- parallel::makeCluster(2, type="FORK")
  # d.phylolm<- dSep(g3, FUN=phylolm:phylolm, n=nrow(shorebird$data), cl=cl, data=shorebird)
  # expect_is(d.phylolm, "dSep")
  # d.phylolm<- dSep(gPhy, FUN="pgls", n=nrow(shorebird$data), cl=cl, data=shorebird)
  # expect_is(d.phylolm, "dSep")
  # parallel::stopCluster(cl)


  ## Example 2
  d.phylolm<- dSep(m, FUN=phylolm::phylolm, cl=4, orderResponse=c("BM", "RS", "LS", "NL", "DD"), data=com.dat$data, phy=com.dat$phy)
  # expect_is(d.phylolm, "dSep")
})
## brm ----

context("dSep-FUN=brm")

test_that("brm works", {
  skip_if_not_installed("brms")
  # library(brms)
  # d.brmfit<- dSep(g1, FUN=brms::brm, n=nrow(d), pathCoef=FALSE, data=d)
  # expect_is(d.brmfit, "dSep")
  # d.brmfit<- dSep(g, FUN=brms::brm, n=nrow(d), pathCoef=FALSE,  data=d)
  # expect_is(d.brmfit, "dSep")
  #
  # d.brmfit<- dSep(g1, FUN=brms::brm, n=nrow(d), cl=2, pathCoef=FALSE, data=d)
  # expect_is(d.brmfit, "dSep")
  # d.brmfit<- dSep(g, FUN=brms::brm, n=nrow(d), cl=2, pathCoef=FALSE, data=d)
  # expect_is(d.brmfit, "dSep")

  # cl<- parallel::makeCluster(2, type="FORK")
  # d.brmfit<- dSep(g1, FUN=brms::brm, n=nrow(d), cl=cl, pathCoef=FALSE, data=d)
  # expect_is(d.brmfit, "dSep")
  # d.brmfit<- dSep(g, FUN=brms::brm, n=nrow(d), cl=cl, pathCoef=FALSE, data=d)
  # expect_is(d.brmfit, "dSep")
  # parallel::stopCluster(cl)
})

## MCMCglmm ----

context("dSep-FUN=MCMCglmm")

test_that("MCMCglmm works", {
  skip_if_not_installed("MCMCglmm")
  # d.MCMCglmm<- dSep(g1, FUN=MCMCglmm::MCMCglmm, formulaArg="fixed", n=nrow(d), data=d, verbose=FALSE)
  # expect_is(d.MCMCglmm, "dSep")
  # d.MCMCglmm<- dSep(g, FUN=MCMCglmm::MCMCglmm, formulaArg="fixed", n=nrow(d), data=d, verbose=FALSE)
  # expect_is(d.MCMCglmm, "dSep")
  #
  # d.MCMCglmm<- dSep(g1, FUN=MCMCglmm::MCMCglmm, formulaArg="fixed", n=nrow(d), cl=2, data=d, verbose=FALSE)
  # expect_is(d.MCMCglmm, "dSep")
  # d.MCMCglmm<- dSep(g, FUN=MCMCglmm::MCMCglmm, formulaArg="fixed", n=nrow(d), cl=2, data=d, verbose=FALSE)
  # expect_is(d.MCMCglmm, "dSep")

  # cl<- parallel::makeCluster(2, type="FORK")
  # d.MCMCglmm<- dSep(g1, FUN=MCMCglmm::MCMCglmm, formulaArg="fixed", n=nrow(d), cl=cl, data=d, verbose=FALSE)
  # expect_is(d.MCMCglmm, "dSep")
  # d.MCMCglmm<- dSep(g, FUN=MCMCglmm::MCMCglmm, formulaArg="fixed", n=nrow(d), cl=cl, data=d, verbose=FALSE)
  # expect_is(d.MCMCglmm, "dSep")
  # parallel::stopCluster(cl)
})


## dSep methods ----
context("Methods for dSep objects")

# load(system.file(c("extdata/sampleModels.RData"), package="dSep"), verbose=FALSE)

## Calls with missing "n"
D<- list()
suppressMessages(D$lm<- dSep(g, FUN="lm", data=d))
suppressMessages(D$glm<- dSep(g, FUN="glm", pathCoef=FALSE, data=d))
suppressMessages(D$gls<- dSep(g1, FUN=nlme::gls, formulaArg="model", data=d))
suppressMessages(D$pgls<- dSep(m, FUN=caper::pgls, cl=parallel::detectCores(), orderResponse=c("BM", "RS", "LS", "NL", "DD"), data=com.dat, lambda="ML"))
suppressMessages(D$phylolm<- dSep(m, FUN=phylolm::phylolm, cl=parallel::detectCores(), orderResponse=c("BM", "RS", "LS", "NL", "DD"), data=com.dat$data, phy=com.dat$phy))
## TODO ERROR: The model does not contain posterior samples. brm
# suppressMessages(D$brm0<- dSep(g1, FUN=brms::brm, cl=parallel::detectCores(), pathCoef=FALSE, data=d))
# suppressMessages(D$brm<- dSep(m, FUN=brms::brm, cl=parallel::detectCores(), pathCoef=FALSE, data=rhino.dat))
suppressMessages(D$MCMCglmm<- dSep(g1, FUN=MCMCglmm::MCMCglmm, formulaArg="fixed", n=nrow(d), cl=parallel::detectCores(), data=d, verbose=FALSE))

sapply(D, function(x) class(x))

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

# Save models for test-LM_API.R
# M<- lapply(D, function(x) x$model[[1]])
# M$brm<- brms::brm(count ~ log_Age_c + log_Base4_c * Trt_c + (1|patient) + (1|visit) + (1|obs), data=brms::epilepsy, family = poisson(),
#                   prior = c(set_prior("student_t(5,0,10)", class = "b"), set_prior("cauchy(0,2)", class = "sd")), cluster=parallel::detectCores())
#
# save(M, file="inst/extdata/sampleModels.RData")
