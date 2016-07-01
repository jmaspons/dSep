## Run dSep() ----
# with FUN parameter as lm, glm, gls, pgls, MCMCglmm and brm
context("dSep-FUN")

## TODO: Rename context
## TODO: Add more tests

g1<- gRbase::dag(~a:c:d + b:d)
g2<- gRbase::dag(~a:c:d + b:d:a)
g<- list(m1=g1, m2=g2)
d<- data.frame(a=rnorm(100), b=rnorm(100), c=rnorm(100), d=rnorm(100))

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

test_that("gls works", {
  skip_if_not_installed("nlme")
  d.gls<- dSep(g1, FUN=nlme::gls, formulaArg="model", nobs=nrow(d), data=d)
  expect_is(d.gls, "dSep")
  d.gls<- dSep(g, FUN=nlme::gls, formulaArg="model", nobs=nrow(d), data=d)
  expect_is(d.gls, "dSep")

  d.gls<- dSep(g1, FUN=nlme::gls, formulaArg="model", nobs=nrow(d), cl=2, data=d)
  expect_is(d.gls, "dSep")
  d.gls<- dSep(g, FUN=nlme::gls, formulaArg="model", nobs=nrow(d), cl=2, data=d)
  expect_is(d.gls, "dSep")

  cl<- parallel::makeCluster(2, type="FORK")
  d.gls<- dSep(g1, FUN=nlme::gls, formulaArg="model", nobs=nrow(d), cl=cl, data=d)
  expect_is(d.gls, "dSep")
  d.gls<- dSep(g, FUN=nlme::gls, formulaArg="model", nobs=nrow(d), cl=cl, data=d)
  expect_is(d.gls, "dSep")
  parallel::stopCluster(cl)
})

test_that("pgls works", {
  skip_if_not_installed("caper")
  data(shorebird, package="caper")
  shorebird.data[,2:5]<- scale(log(shorebird.data[,2:5]))
  shorebird<- caper::comparative.data(shorebird.tree, shorebird.data, 'Species')
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
  d.pgls<- dSep(gPhy, FUN=caper::pgls, nobs=nrow(shorebird$data), cl=cl, data=shorebird)
  expect_is(d.pgls, "dSep")
  parallel::stopCluster(cl)
})

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

  cl<- parallel::makeCluster(2, type="FORK")
  d.brmfit<- dSep(g1, FUN=brms::brm, nobs=nrow(d), cl=cl, pathCoef=FALSE, data=d)
  expect_is(d.brmfit, "dSep")
  # d.brmfit<- dSep(g, FUN=brms::brm, nobs=nrow(d), cl=cl, pathCoef=FALSE, data=d)
  expect_is(d.brmfit, "dSep")
  parallel::stopCluster(cl)
})

test_that("MCMCglmm works", {
  skip_if_not_installed("MCMCglmm")
  d.MCMCglmm<- dSep(g1, FUN=MCMCglmm::MCMCglmm, formulaArg="fixed", nobs=nrow(d), data=d, verbose=FALSE)
  expect_is(d.MCMCglmm, "dSep")
  d.MCMCglmm<- dSep(g, FUN=MCMCglmm::MCMCglmm, formulaArg="fixed", nobs=nrow(d), data=d, verbose=FALSE)
  expect_is(d.MCMCglmm, "dSep")

  d.MCMCglmm<- dSep(g1, FUN=MCMCglmm::MCMCglmm, formulaArg="fixed", nobs=nrow(d), cl=2, data=d, verbose=FALSE)
  expect_is(d.MCMCglmm, "dSep")
  d.MCMCglmm<- dSep(g, FUN=MCMCglmm::MCMCglmm, formulaArg="fixed", nobs=nrow(d), cl=2, data=d, verbose=FALSE)
  expect_is(d.MCMCglmm, "dSep")

  cl<- parallel::makeCluster(2, type="FORK")
  d.MCMCglmm<- dSep(g1, FUN=MCMCglmm::MCMCglmm, formulaArg="fixed", nobs=nrow(d), cl=cl, data=d, verbose=FALSE)
  expect_is(d.MCMCglmm, "dSep")
  d.MCMCglmm<- dSep(g, FUN=MCMCglmm::MCMCglmm, formulaArg="fixed", nobs=nrow(d), cl=cl, data=d, verbose=FALSE)
  expect_is(d.MCMCglmm, "dSep")
  parallel::stopCluster(cl)
})

## Common API for lm, glm, gls, pgls, MCMCglmm and brm ----
# library(caper)
# library(nlme)
# library(MCMCglmm)
# library(brms)

D<- list(lm=d.lm, glm=d.glm, gls=d.gls, pgls=d.pgls, MCMCglmm=d.MCMCglmm, brmsfit=d.brmfit)
M<- lapply(D, function(x) x$model[[1]])

test_that("models have the right classes",{
  sapply(D, function(x) class(x))
  sapply(D, function(x) { expect_is(x, "dSep")})

  sapply(M, function(x) class(x))
  sapply(M, function(x) { expect_is(x, c("lm", "glm", "gls", "pgls", "MCMCglmm", "brmsfit")); summary(x); })
})

context("Genergic methods for lm, glm, gls, pgls, MCMCglmm and brm")

test_that("pValues() works",{
  sapply(M, function(x) try(pValues(x)))
  lapply(M, function(x) expect_type(try(pValues(x)), "double"))
})

test_that("scoef() works",{
  sapply(M, function(x) try(scoef(x)))
  lapply(M, function(x) expect_type(try(scoef(x)), "double"))
})

test_that("nobs() works",{
  sapply(M, function(x) try(nobs(x))) ## Fail for MCMCglmm
  lapply(M, function(x) expect_type(try(nobs(x)), "double"))
})

## dSep methods ----
context("dSep-methods")

test_that("print works", {
  lapply(D, function(x) return(x))
  lapply(D, function(x) return(x$pathCoefficients))
})

test_that("plot works", {
  skip_if_not_installed("Rgraphviz")

  lapply(D, function(x) try(plot(x)))
  lapply(D, function(x) try(plot(x$pathCoefficients)))
})

