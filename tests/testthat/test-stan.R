test_that("the constant model fits using Stan", {
  capture.output(fit <- constant_stan(mw, nsim_warm = 200, nsim_iter = 400, cores = 2))
  expect_s4_class(fit, "stanfit")
})

test_that("the Besag model fits using Stan", {
  capture.output(fit <- besag_stan(mw, nsim_warm = 200, nsim_iter = 400, cores = 2))
  expect_s4_class(fit, "stanfit")
})

test_that("the BYM2 model fits using Stan", {
  capture.output(fit <- bym2_stan(mw, nsim_warm = 200, nsim_iter = 400, cores = 2))
  expect_s4_class(fit, "stanfit")
})

test_that("the FCK model fits using Stan", {
  capture.output(fit <- fck_stan(mw, nsim_warm = 200, nsim_iter = 400, cores = 2))
  expect_s4_class(fit, "stanfit")
})

test_that("the FIK model fits using Stan", {
  capture.output(fit <- fik_stan(mw, nsim_warm = 200, nsim_iter = 400, cores = 2))
  expect_s4_class(fit, "stanfit")
})

test_that("the CK model fits using Stan", {
  capture.output(fit <- ck_stan(mw, nsim_warm = 200, nsim_iter = 400, cores = 2))
  expect_s4_class(fit, "stanfit")
})

test_that("the IK model fits using Stan", {
  capture.output(fit <- ik_stan(mw, nsim_warm = 200, nsim_iter = 400, cores = 2))
  expect_s4_class(fit, "stanfit")
})