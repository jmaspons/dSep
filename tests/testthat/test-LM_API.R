## Run and save dummy models ----
# lm, glm, gls, pgls, MCMCglmm and brm
## WARNING: test_that runs in an isolated environment thus variables won't be available outside the functions

# g1<- gRbase::dag(~a:c:d + b:d)
# g2<- gRbase::dag(~a:c:d + b:d:a)
# g<- list(m1=g1, m2=g2)
# d<- data.frame(a=rnorm(100), b=rnorm(100), c=rnorm(100), d=rnorm(100))
#
# cl<- parallel::makeCluster(2, type="FORK")
#
# d.lm<- dSep(g1, FUN="lm", n=nrow(d), data=d)
# d.glm<- dSep(g1, FUN="glm", n=nrow(d), data=d)
# d.gls<- dSep(g, FUN=nlme::gls, formulaArg="model", n=nrow(d), data=d)
#
# data(shorebird, package="caper")
# shorebird.data[,2:5]<- scale(log(shorebird.data[,2:5]))
# shorebird<- caper::comparative.data(phy=shorebird.tree, data=shorebird.data, names.col=Species, )
# g3<- gRbase::dag(~Egg.Mass:M.Mass:F.Mass)
# g4<- gRbase::dag(~Egg.Mass:M.Mass:Cl.size + Cl.size:F.Mass:M.Mass)
# gPhy<- list(mPhy1=g3, mPhy2=g4)
# d.pgls<- dSep(gPhy, FUN="pgls", n=nrow(shorebird$data), cl=cl, data=shorebird)
#
# library(brms)
# d.brmfit<- dSep(g1, FUN=brms::brm, n=nrow(d), cl=cl, pathCoef=FALSE, data=d)
# d.MCMCglmm<- dSep(g, FUN=MCMCglmm::MCMCglmm, formulaArg="fixed", n=nrow(d), cl=cl, data=d, verbose=FALSE)
#
# parallel::stopCluster(cl)
# D<- list(lm=d.lm, glm=d.glm, gls=d.gls, pgls=d.pgls, MCMCglmm=d.MCMCglmm, brmsfit=d.brmfit)
# M<- lapply(D, function(x) x$model[[1]])
# sapply(D, function(x) class(x))
# sapply(M, function(x) class(x))
# save(D, M, file="inst/extdata/sampleModels.RData", compress=TRUE)

## Common API for lm, glm, gls, pgls, MCMCglmm and brm ----
# library(caper)
# library(nlme)
# library(MCMCglmm)
# library(brms)


context("lm, glm, gls, pgls, MCMCglmm and brm models")

load(system.file(c("extdata/sampleModels.RData"), package="dSep"), verbose=FALSE)

test_that("models have the right classes",{
  sapply(D, function(x) { expect_is(x, "dSep")})
  sapply(M, function(x) { expect_is(x, c("lm", "glm", "gls", "pgls", "MCMCglmm", "brmsfit")); summary(x); })
})

## Generic methods for lineal models ----
context("Genergic methods for lineal models")

test_that("pValues() works",{
  # sapply(M, function(x) try(pValues(x)))
  lapply(M, function(x) expect_type(try(pValues(x)), "double"))
})

test_that("scoef() works",{
  # sapply(M, function(x) try(scoef(x)))
  lapply(M, function(x) expect_type(try(scoef(x)), "double"))
})

# library(caper)
# library(nlme)
# library(MCMCglmm)
# library(brms)

test_that("nobs() works",{
  # sapply(M, function(x) try(nobs(x))) ## Fail for MCMCglmm
  lapply(M, function(x) expect_equal(try(mode(nobs(x))), "numeric"))
})

