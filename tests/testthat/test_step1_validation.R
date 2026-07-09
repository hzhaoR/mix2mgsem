test_that("Step 1 validation recovers MM-clusters and factor loadings", {
  dat <- readRDS(testthat::test_path(
    "fixtures",
    "data19346.rds"
  ))

  truth <- readRDS(testthat::test_path(
    "fixtures",
    "truth_step1_19346.rds"
  ))
  data <- as.data.frame(dat$SimData)
  # data <- data[,c(21,1:20)]

  s1out <- MixMix_Step1(
    data = data,
    step1model = truth$step1model,
    MM.cluster.spec = c("loadings"),
    MM.nclus = length(unique(truth$true_cluster)),
    MM.design = truth$true_design,
    invar_loadings = rep(c(1,0,0,1,1), 4),
    markers = truth$markers
  )

  # s1out <- readRDS(testthat::test_path(
  #   "fixtures",
  #   "s1output19346.rds"
  # ))

  step1_eval <- evaluate_step1_validation(
    posteriors = s1out$mmgfa_output$cluster_memb,
    lambda_ks = s1out$mmgfa_output$lambda_gs,
    cov_eta = s1out$cov_eta,

    true_cluster = truth$true_cluster,
    true_lambda = truth$true_lambda,
    true_cov_eta = truth$true_cov_eta,

    step1model = truth$step1model,
    observed_vars = s1out$vars,
    latent_vars = s1out$lat_var,
    markers = truth$markers,
    include_marker_loadings = FALSE
  )

  expect_equal(step1_eval$clustering$misclassification_error, 0)
  expect_true(step1_eval$clustering$correct_clustering)
  expect_lt(step1_eval$clustering$mean_uncertainty, 0.01)

  expect_lt(step1_eval$loading$rmse, 0.05)
  expect_lt(step1_eval$loading$max_abs_error, 0.10)

  # expect_lt(step1_eval$factor_covariance$rmse, 0.10)
  # expect_lt(step1_eval$factor_covariance$max_abs_error, 0.20)
})
