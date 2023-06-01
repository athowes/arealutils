test_that("the constant model fits using aghq", {
  capture.output(fit <- constant_aghq(mw, its = 100))
  expect_s3_class(fit, "aghq")
})

test_that("the Besag model fits using aghq", {
  capture.output(fit <- besag_aghq(mw, its = 100))
  expect_s3_class(fit, "aghq")
})

test_that("the BYM2 model fits using aghq", {
  capture.output(fit <- bym2_aghq(mw, its = 100))
  expect_s3_class(fit, "aghq")
})

test_that("the FCK model fits using aghq", {
  capture.output(fit <- fck_aghq(mw, its = 100))
  expect_s3_class(fit, "aghq")
})

test_that("the FIK model fits using aghq", {
  capture.output(fit <- fik_aghq(mw, its = 100))
  expect_s3_class(fit, "aghq")
})

test_that("the CK model fits using aghq", {
  capture.output(fit <- ck_aghq(mw, its = 100))
  expect_s3_class(fit, "aghq")
})

test_that("the IK model fits using aghq", {
  capture.output(fit <- ik_aghq(mw, its = 100))
  expect_s3_class(fit, "aghq")
})