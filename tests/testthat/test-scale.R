test_that("get_scale, riebler_gv, scale_gmrf_precision agree", {
  Q_singular <- matrix(c(1, -1, -1, 1), nrow = 2, ncol = 2) # Besag precision for connected pair
  Q_singular_inv <- matrix(c(0.25, -0.25, -0.25, 0.25), nrow = 2, ncol = 2)
  expect_equal(get_scale(Q_singular), riebler_gv(Q_singular_inv), tolerance = 1e-3)
  expect_equal(riebler_gv(Q_singular_inv), scale_gmrf_precision(Q_singular)$scales, tolerance = 1e-3)
})

test_that("scale_gmrf_precision works as expected for two connected components", {
  # Besag precision for a pair of connected pairs
  Q_singular <- matrix(c(1, -1, 0, 0, -1, 1, 0, 0, 0, 0, 1, -1, 0, 0, -1, 1), nrow = 4, ncol = 4)
  Q_scaled <- 0.25 * Q_singular
  expect_equal(scale_gmrf_precision(Q_singular)$scales, c(0.25, 0.25), tolerance = 1e-3)
  expect_equal(scale_gmrf_precision(Q_singular)$Q, Q_scaled, tolerance = 1e-3)
})
