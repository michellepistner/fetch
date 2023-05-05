library(ALDEx2)
library(plyr)
library(dplyr)
library(fetch)

test_that("globalFetchconds", {
  data(selex)
  #subset for efficiency
  selex <- selex[1201:1600,]
  conds <- c(rep("NS", 7), rep("S", 7))

  set.seed(1)
  expect_warning(fetch(selex, conds, denom = "all", test = "global", mc.samples = 3))
  expect_warning(fetch(selex, conds, denom = "all", test = "global", mc.samples = 3, gl.test = "anova"))
})

test_that("globalFetchmatrix", {
  data(selex)
  #subset for efficiency
  selex <- selex[1201:1600,]
  covariates <- data.frame("A" = sample(0:1, 14, replace = TRUE),
                           "B" = c(rep(0, 7), rep(1, 7)))
  mm <- model.matrix(~ A + B, covariates)
  
  set.seed(1)
  expect_warning(fetch(selex, mm, denom = "all", test = "global", mc.samples = 3))
  expect_warning(fetch(selex, mm, denom = "all", test = "global", mc.samples = 3, gl.test = "anova"))
  
})

test_that("localFetch", {
  set.seed(1)
  Y = rmultinom(10, 10000, prob = seq(1:100))
  X = matrix(sample(1:2,10,replace=TRUE),ncol = 10)
  out =fetch(Y, c(X), denom = "iqlr", test = "t")

  expect_equal(out[1,1], -5.000539551)
})

test_that("localFetchvsAldex", {
  set.seed(1)
  Y = rmultinom(10, 10000, prob = seq(1:15))
  X = matrix(sample(1:2,10,replace=TRUE),ncol = 10)
  out =fetch(Y, c(X), denom = "all", test = "t")

  out2 = ALDEx2::aldex(Y,c(X),bayesEst = FALSE)
  sig2 = ifelse(out2[,11] > 0.05,0,1)
  sig = ifelse(out[,11] > 0.05,0,1)

  expect_equal(sig, sig2)
})

