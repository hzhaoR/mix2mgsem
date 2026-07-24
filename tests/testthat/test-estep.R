test_that("EStep returns posterior cluster membership probabilities", {
  loglik <- matrix(
    c(-100, -106,
      -108, -101,
      -103, -103,
      -100, -101),
    nrow = 4,
    byrow = TRUE
  )

  out <- mix2mgsem:::EStep(
    pi_ks = c(0.5, 0.5),
    ngroup = 4,
    nclus = 2,
    loglik = loglik
  )

  expect_true(is.matrix(out))
  expect_equal(dim(out), c(4, 2))
  expect_equal(rowSums(out), rep(1, 4), tolerance = 1e-12)

  expect_gt(out[1, 1], out[1, 2])
  expect_gt(out[2, 2], out[2, 1])
  expect_equal(out[3, 1], out[3, 2], tolerance = 1e-12)
  expect_gt(out[4, 1], out[4, 2])
})

test_that("EStep keeps prior probabilities when log-likelihoods are equal", {
  out <- mix2mgsem:::EStep(
    pi_ks = c(0.25, 0.75),
    ngroup = 3,
    nclus = 2,
    loglik = matrix(
      c(-10, -10,
        -20, -20,
        -30, -30),
      nrow = 3,
      byrow = TRUE
    )
  )

  expected <- matrix(
    c(0.25, 0.75,
      0.25, 0.75,
      0.25, 0.75),
    nrow = 3,
    byrow = TRUE
  )

  expect_equal(out, expected, tolerance = 1e-12)
})
