#' Call lavaan's internal sample-statistics log-likelihood helper
#'
#' Internal wrapper around lavaan's sample-statistics log-likelihood function.
#' This avoids direct use of the `:::` operator in the main estimation code.
#'
#' @noRd
lav_samplestats_loglik <- function(...) {
  utils::getFromNamespace("lav_mvnorm_loglik_samplestats", "lavaan")(...)
}

