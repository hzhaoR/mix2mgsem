test_that("MixMix_Step1 requires data to be a data frame", {
  expect_error(
    MixMix_Step1(
      data = matrix(1:6, nrow = 3),
      group = "group",
      step1model = "f =~ y1 + y2",
      MM.nclus = 1,
      MM.design = NULL,
      markers = NULL
    ),
    "`data` must be a data frame"
  )
})

test_that("MixMix_Step1 requires one measurement model cluster number", {
  dat <- data.frame(
    group = rep(1:2, each = 3),
    y1 = 1:6,
    y2 = 2:7
  )

  expect_error(
    MixMix_Step1(
      data = dat,
      group = "group",
      step1model = "f =~ y1 + y2",
      MM.nclus = 1:2,
      MM.design = NULL,
      markers = NULL
    ),
    "`MM.nclus`"
  )
})

test_that("MixMix_Step2 requires one structural relations cluster number", {
  expect_error(
    MixMix_Step2(
      s1out = list(),
      step2model = "F2 ~ F1",
      nclus = 1:2,
      nfree = 0
    ),
    "`nclus` must be a single positive integer"
  )
})

test_that("MixMix_Step2 requires one non-missing nfree value", {
  expect_error(
    MixMix_Step2(
      s1out = list(),
      step2model = "F2 ~ F1",
      nclus = 1,
      nfree = NA
    ),
    "`nfree` must be a single non-negative integer"
  )
})
