test_that("the constant model fits using TMB", {
  capture.output(fit <- constant_tmb(mw, its = 100))
  expect_s3_class(fit, "sdreport")
})

test_that("the Besag model fits using TMB", {
  capture.output(fit <- besag_tmb(mw, its = 100))
  expect_s3_class(fit, "sdreport")
})

test_that("the BYM2 model fits using TMB", {
  capture.output(fit <- bym2_tmb(mw, its = 100))
  expect_s3_class(fit, "sdreport")
})

test_that("the FCK model fits using TMB", {
  capture.output(fit <- fck_tmb(mw, its = 100))
  expect_s3_class(fit, "sdreport")
})

test_that("the FIK model fits using TMB", {
  capture.output(fit <- fik_tmb(mw, its = 100))
  expect_s3_class(fit, "sdreport")
})

test_that("the CK model fits using TMB", {
  capture.output(fit <- ck_tmb(mw, its = 100))
  expect_s3_class(fit, "sdreport")
})

test_that("the IK model fits using TMB", {
  capture.output(fit <- ik_tmb(mw, its = 100))
  expect_s3_class(fit, "sdreport")
})