context("sbhm")

# These tests implicitly test that Stan and JAGS give similar results

# Need to allow room for error. The book states that 200 samples were taken from
# the posterior distributions, so tail quantiles will be variable. We could find
# an alternative version online, but that might be subject to change with
# different builds of Stan. Presumably the authors checked the Stan output against
# the printed version that comes from other software (presumably WinBUGS).

# Compare with a published version: BDA 3rd Ed., page 123
bda <- matrix(ncol=5, byrow=TRUE,
              c(-2, 7, 10, 16, 31, 
                -5, 3, 8, 12, 23,
                -11, 2, 7, 11, 19,
                -7, 4, 8, 11, 21,
                -9, 1, 5, 10, 18,
                -7, 2, 6, 10, 28,
                -1, 7, 10, 15, 26,
                -6, 3, 8, 13, 33))

test_that("8 schools example is reproduced accurately with single chain", {
  Sys.setenv("R_TESTS" = "") # Prevent R CMD check failure with parallel

  set.seed(17010702)

  has.stan <- require(rstan)

  if (has.stan){
    smod <- sbhm(schools$effect, schools$se, nchains=1)
    ss <- get.sbhm.summary(smod)[1:8, ]

    # Quartiles should be less variable than the 95% intervals
    expect_equal(unname(ss[, 2:4] / bda[, 2:4]), matrix(1, ncol=3, nrow=8), tolerance=.15,
                 label="Single chain Stan output matches BDA")
    # Allow more divergence for the 95% limits
    expect_equal(unname(ss[, c(1, 5)] / bda[, c(1, 5)]), matrix(1, ncol=2, nrow=8), tolerance=.2,
                 label="Single chain Stan output matches BDA")
  }

  jmod <- sbhm(schools$effect, schools$se, nchains=1, engine="jags")
  sj <- get.sbhm.summary(jmod)[1:8, ]

  # Quartiles should be less variable than the 95% intervals
  expect_equal(unname(sj[, 2:4] / bda[, 2:4]), matrix(1, ncol=3, nrow=8), tolerance=.15,
               label="Single chain Stan output matches BDA")

  # Allow more divergence for the 95% limits
  expect_equal(unname(sj[, c(1, 5)] / bda[, c(1, 5)]), matrix(1, ncol=2, nrow=8), tolerance=.2,
               label="Single chain Stan output matches BDA")
})

test_that("8 schools example is reproduced accurately with multiple chains", {
  set.seed(17610407)

  has.stan <- require(rstan)

  if  (has.stan){
    smod <- sbhm(schools$effect, schools$se, nchains=4)
    ss <- get.sbhm.summary(smod)[1:8, ]
    
    # Quartiles should be less variable than the 95% intervals
    expect_equal(unname(ss[, 2:4] / bda[, 2:4]), matrix(1, ncol=3, nrow=8), tolerance=.15,
                 label="Multiple chain Stan output matches BDA")
    # Allow more divergence for the 95% limits
    expect_equal(unname(ss[, c(1, 5)] / bda[, c(1, 5)]), matrix(1, ncol=2, nrow=8), tolerance=.2,
                 label="Multiple chain Stan output matches BDA")
  }

  jmod <- sbhm(schools$effect, schools$se, nchains=4, engine="jags")
  sj <- get.sbhm.summary(jmod)[1:8, ]

  # Quartiles should be less variable than the 95% intervals
  expect_equal(unname(sj[, 2:4] / bda[, 2:4]), matrix(1, ncol=3, nrow=8), tolerance=.15,
               label="Multiple chain Stan output matches BDA")

  # Allow more divergence for the 95% limits
  expect_equal(unname(sj[, c(1, 5)] / bda[, c(1, 5)]), matrix(1, ncol=2, nrow=8), tolerance=.2,
               label="Multiple chain Stan output matches BDA")
})
