test_that("MixMix_Step1 requires data to be a data frame", {
  expect_error(
    MixMix_Step1(
      data = matrix(1:6, nrow = 3),
      step1model = "f =~ y1 + y2",
      group = "group",
      MM.nclus = 1,
      MM.design = NULL,
      invar_loadings = NULL,
      markers = NULL
    ),
    "`data` must be a data frame"
  )
})

test_that("MixMix_Step1 gives informative error for missing group column", {
  dat <- data.frame(
    y1 = 1:4,
    y2 = 2:5
  )

  expect_error(
    MixMix_Step1(
      data = dat,
      step1model = "f =~ y1 + y2",
      group = "country",
      MM.nclus = 1,
      MM.design = NULL,
      invar_loadings = NULL,
      markers = NULL
    ),
    "`group` must be a single column name"
  )
})

test_that("MixMix_Step1 requires one selected measurement model cluster number", {
  dat <- data.frame(
    group = rep(1:2, each = 3),
    y1 = 1:6,
    y2 = 2:7
  )

  expect_error(
    MixMix_Step1(
      data = dat,
      step1model = "f =~ y1 + y2",
      group = "group",
      MM.nclus = 1:2,
      MM.design = NULL,
      invar_loadings = NULL,
      markers = NULL
    ),
    "`MM.nclus`"
  )
})

test_that("MixMix_Step2 validates required Step 1 output structure", {
  expect_error(
    MixMix_Step2(
      s1out = list(),
      step2model = "F2 ~ F1",
      nclus = 1,
      nfree = 0
    ),
    "`s1out` is missing required components"
  )
})

test_that("MixMix_Step2 validates nclus and nfree after Step 1 structure is present", {
  fake_s1out <- list(
    ngroups = 2,
    N_gs = c(20, 20),
    S_unbiased = list(diag(2), diag(2)),
    vars = c("y1", "y2"),
    lat_var = c("F1", "F2"),
    mmgfa_output = list(
      lambda_gs = list(diag(2)),
      cluster_memb = matrix(c(1, 1), nrow = 2),
      theta_gs = list(diag(2), diag(2))
    )
  )

  expect_error(
    MixMix_Step2(
      s1out = fake_s1out,
      step2model = "F2 ~ F1",
      nclus = 0,
      nfree = 0
    ),
    "`nclus` must be a single positive integer"
  )

  expect_error(
    MixMix_Step2(
      s1out = fake_s1out,
      step2model = "F2 ~ F1",
      nclus = 1,
      nfree = -1
    ),
    "`nfree` must be a single non-negative integer"
  )
})
