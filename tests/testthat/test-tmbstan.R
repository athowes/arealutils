test_that("the constant model fits using tmbstan", {
  capture.output(fit <- constant_tmbstan(mw, nsim_warm = 200, nsim_iter = 400, cores = 2))
  expect_s4_class(fit, "stanfit")
})

test_that("the Besag model fits using tmbstan", {
  capture.output(fit <- besag_tmbstan(mw, nsim_warm = 200, nsim_iter = 400, cores = 2))
  expect_s4_class(fit, "stanfit")
})

test_that("the BYM2 model fits using tmbstan", {
  capture.output(fit <- bym2_tmbstan(mw, nsim_warm = 200, nsim_iter = 400, cores = 2))
  expect_s4_class(fit, "stanfit")
})

test_that("the FCK model fits using tmbstan", {
  capture.output(fit <- fck_tmbstan(mw, nsim_warm = 200, nsim_iter = 400, cores = 2))
  expect_s4_class(fit, "stanfit")
})

test_that("the FIK model fits using tmbstan", {
  capture.output(fit <- fik_tmbstan(mw, nsim_warm = 200, nsim_iter = 400, cores = 2))
  expect_s4_class(fit, "stanfit")
})

test_that("the CK model fits using tmbstan", {
  capture.output(fit <- ck_tmbstan(mw, nsim_warm = 200, nsim_iter = 400, cores = 2))
  expect_s4_class(fit, "stanfit")
})

test_that("the IK model fits using tmbstan", {
  capture.output(fit <- ik_tmbstan(mw, nsim_warm = 200, nsim_iter = 400, cores = 2))
  expect_s4_class(fit, "stanfit")
})