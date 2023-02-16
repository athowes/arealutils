test_that("the constant model fits using R-INLA", {
  fit <- constant_inla(mw, cores = 2)
  expect_s3_class(fit, "inla")
})

test_that("the Besag model fits using R-INLA", {
  fit <- besag_inla(mw, cores = 2)
  expect_s3_class(fit, "inla")
})

test_that("the BYM2 model fits using R-INLA", {
  fit <- bym2_inla(mw, cores = 2)
  expect_s3_class(fit, "inla")
})

test_that("the FCK model fits using R-INLA", {
  fit <- fck_inla(mw, cores = 2)
  expect_s3_class(fit, "inla")
})

test_that("the FIK model fits using R-INLA", {
  fit <- fik_inla(mw, cores = 2)
  expect_s3_class(fit, "inla")
})
