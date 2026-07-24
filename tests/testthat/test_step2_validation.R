test_that("Step2 recovers SR-clusters and betas", {
  s1out <- readRDS(testthat::test_path(
    "fixtures",
    "s1output19346.rds"
  ))

  truth <- readRDS(testthat::test_path(
    "fixtures",
    "truth_step2_19346.rds"
  ))

  s2out <- MixMix_Step2(
    s1out = s1out,
    step2model = truth$step2model,
    nclus = truth$nclus,
    nfree = truth$nfree,
    seed = 100,
    printing = FALSE
  )

  step2_eval <- evaluate_step2_validation(
    posteriors = s2out$posteriors,
    beta_ks = s2out$param$beta_ks,
    true_cluster = truth$true_cluster,
    true_beta = truth$true_beta,
    step2model = truth$step2model
  )

  expect_equal(step2_eval$clustering$misclassification_error, 0)
  expect_lt(step2_eval$clustering$mean_uncertainty, 0.01)

  expect_lt(step2_eval$beta$rmse, 0.05)
  expect_lt(step2_eval$beta$max_abs_error, 0.10)
})
