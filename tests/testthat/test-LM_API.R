## Common API for lm, glm, gls, pgls, MCMCglmm and brm ----
# library(caper)
# library(nlme)
# library(MCMCglmm)
# library(brms)


context("lm, glm, gls, pgls, MCMCglmm and brm models")


# D comes from test-dSep.R. Tests ordered alphabetically.
M<- lapply(D, function(x) x$model[[1]])
sapply(D, function(x) class(x))
sapply(M, function(x) class(x))

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

