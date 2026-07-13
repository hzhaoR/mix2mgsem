#' Select the number of SR-clusters for Step 2
#'
#' Fit [MixMix_Step2()] for several candidate numbers of
#' SR-clusters and collect the resulting information criteria in one overview.
#'
#' @param s1out Output from [MixMix_Step1()].
#' @param step2model A lavaan model syntax string specifying the structural relations among latent variables.
#' @param nclus Integer vector containing at least two candidate numbers of SR-clusters.
#' @param nfree Integer scalar indicating the number of free non-marker loadings in the Step 1 measurement model.
#'   Used for observed-data information-criterion calculations.
#' @param seed Optional integer seed for random starts.
#' @param single Logical. If `FALSE`, Step 2 uses `s1out$cov_eta`. If `TRUE`,
#'   Step 2 uses `s1out$cov_eta_fs`. The single-indicator option is experimental.
#' @param max_it Maximum number of EM iterations per random start.
#' @param nstarts Number of random starts.
#' @param printing Logical; whether to print iteration progress.
#' @param partition Character string specifying the initial partition type.
#'   Currently `"hard"` and `"soft"` are supported.
#' @param Endo2Cov Logical; whether to allow covariance between endogenous 2 variables.
#' @param allG Logical; whether the endogenous covariances are group-specific (`TRUE`) or cluster-specific (`FALSE`).
#' @param fit Character string indicating the likelihood target for model estimation. Currently `"factors"` is the supported option.
#'
#' @return An object of class `mix2mgsem_step2_selection`, containing:
#' \describe{
#'   \item{fits}{Named list of successful [MixMix_Step2()] fits. Failed fits are stored as `NULL`.}
#'   \item{overview}{Data frame containing the fit and information criteria for every candidate number of SR-clusters.}
#'   \item{selected}{Numbers of SR-clusters selected by the information criteria.}
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' s2_selection <- MixMix_Step2_select(
#'   s1out = s1out,
#'   step2model = step2_model,
#'   nclus = c(1:5),
#'   nfree = 4,
#'   seed = 100
#' )
#'
#' s2_selection$overview
#'
#' plot(
#'   s2_selection,
#'   criteria = c("BIC_G", "BIC_N", "AIC", "AIC3")
#' )
#' }
MixMix_Step2_select <- function(s1out, step2model, nclus, nfree,
                                seed = NULL, single = FALSE,
                                max_it = 10000L, nstarts = 50L, printing = FALSE,
                                partition = "hard", Endo2Cov = TRUE, allG = TRUE, fit = "factors"
) {
  valid_nclus <- is.numeric(nclus) &&
    all(is.finite(nclus)) &&
    all(nclus >= 1L & nclus %% 1L == 0L) &&
    length(unique(nclus)) >= 2L

  if (!valid_nclus) {
    stop(
      "`nclus` must contain at least two distinct positive integers.",
      call. = FALSE
    )
  }

  nclus <- sort(unique(as.integer(nclus)))

  if (any(nclus > s1out$ngroups)) {
    stop(
      "Candidate values in `nclus` cannot exceed the number of groups (",
      s1out$ngroups,
      ").",
      call. = FALSE
    )
  }

  fits <- stats::setNames(
    vector("list", length(nclus)),
    as.character(nclus)
  )

  overview_rows <- list()
  model_errors <- character()

  for (k in nclus) {
    if (printing) {
      message("Fitting Step 2 model with ", k, " cluster(s).")
    }

    current_fit <- tryCatch(
      MixMix_Step2(
        s1out = s1out,
        step2model = step2model,
        nclus = k,
        nfree = nfree,
        seed = seed,
        userStart = NULL,
        single = single,
        max_it = max_it,
        nstarts = nstarts,
        printing = printing,
        partition = partition,
        Endo2Cov = Endo2Cov,
        allG = allG,
        fit = fit
      ),
      error = identity
    )

    if (inherits(current_fit, "error")) {
      model_errors[as.character(k)] <- conditionMessage(current_fit)

      if (printing) {
        message(
          "Model with ", k, " cluster(s) failed: ",
          model_errors[as.character(k)]
        )
      }
      next
    }

    fits[[as.character(k)]] <- current_fit

    overview_rows[[as.character(k)]] <- with(
      current_fit,
      data.frame(
        nclus = as.integer(k),
        logLik_observed = as.numeric(logLik$obs_loglik),
        nrpar_observed = as.integer(NrPar$Obs.nrpar),
        BIC_G_observed = as.numeric(model_sel$BIC$observed$BIC_G),
        BIC_N_observed = as.numeric(model_sel$BIC$observed$BIC_N),
        AIC_observed = as.numeric(model_sel$AIC$observed),
        AIC3_observed = as.numeric(model_sel$AIC3$observed),
        logLik_factors = as.numeric(logLik$loglik),
        nrpar_factors = as.integer(NrPar$Fac.nrpar),
        BIC_G_factors = as.numeric(model_sel$BIC$Factors$BIC_G),
        BIC_N_factors = as.numeric(model_sel$BIC$Factors$BIC_N),
        AIC_factors = as.numeric(model_sel$AIC$Factors),
        AIC3_factors = as.numeric(model_sel$AIC3$Factors),
        time_minutes = as.numeric(step2_time, units = "mins")
      )
    )
  }

  if (length(overview_rows) == 0L) {
    error_details <- paste(
      paste0(
        "nclus = ",
        names(model_errors),
        ": ",
        unname(model_errors)
      ),
      collapse = "; "
    )

    stop(
      "All candidate Step 2 models failed. ",
      error_details,
      call. = FALSE
    )
  }

  overview <- do.call(
    rbind,
    unname(overview_rows)
  )

  rownames(overview) <- NULL

  if (nrow(overview) == 1L) {
    warning(
      "Only one candidate model was fitted successfully; the information ",
      "criteria cannot be meaningfully compared.",
      call. = FALSE
    )
  }

  select_minimum <- function(column) {
    values <- overview[[column]]
    valid <- is.finite(values)

    if (!any(valid)) {
      return(integer(0))
    }

    minimum <- min(values[valid])

    overview$nclus[valid & values == minimum]
  }

  criteria <- list(
    observed = c(
      BIC_G = "BIC_G_observed",
      BIC_N = "BIC_N_observed",
      AIC = "AIC_observed",
      AIC3 = "AIC3_observed"
    ),
    factors = c(
      BIC_G = "BIC_G_factors",
      BIC_N = "BIC_N_factors",
      AIC = "AIC_factors",
      AIC3 = "AIC3_factors"
    )
  )

  selected <- lapply(criteria, function(columns) {
    stats::setNames(
      lapply(unname(columns), select_minimum),
      names(columns)
    )
  })

  structure(
    list(
      call = match.call(),
      fits = fits,
      overview = overview,
      selected = selected,
      errors = model_errors
    ),
    class = "mix2mgsem_step2_selection"
  )
}

#' Print a Step 2 model selection object
#'
#' @param x An object returned by [MixMix_Step2_select()].
#' @param ... Additional arguments passed to `print.data.frame()`.
#'
#' @return `x`, invisibly.
#' @export
print.mix2mgsem_step2_selection <- function(x, ...) {
  cat("Step 2 model selection\n\n")

  display_columns <- c(
    "nclus",
    "BIC_G_observed",
    "BIC_N_observed",
    "AIC_observed",
    "AIC3_observed"
  )

  print(
    x$overview[, display_columns, drop = FALSE],
    row.names = FALSE,
    ...
  )

  cat("\nSelected number(s) of SR-clusters\n\n")

  selected_table <- data.frame(
    Criterion = names(x$selected$observed),
    Observed = vapply(
      x$selected$observed,
      paste,
      collapse = ", ",
      FUN.VALUE = character(1)
    ),
    Factors = vapply(
      x$selected$factors,
      paste,
      collapse = ", ",
      FUN.VALUE = character(1)
    ),
    check.names = FALSE
  )

  print(
    selected_table,
    row.names = FALSE
  )

  if (length(x$errors) > 0L) {
    cat("\nFailed candidate models\n")
    print(x$errors)
  }

  invisible(x)
}

#' Plot Step 2 model-selection criteria
#'
#' Plot information criterion values against the candidate number of SR-clusters.
#'
#' @param x An object returned by [MixMix_Step2_select()].
#' @param level Character string indicating whether to plot criteria based on
#'   the observed-data likelihood (`"observed"`) or factor likelihood (`"factors"`).
#' @param criteria Character vector specifying the information criteria to
#'   plot. Available options are `"BIC_G"`, `"BIC_N"`, `"AIC"`, and `"AIC3"`.
#' @param ... Additional graphical arguments passed to [graphics::plot()].
#'
#' @return The plotted model selection values, invisibly.
#' @export
plot.mix2mgsem_step2_selection <- function(
    x,
    level = c("observed", "factors"),
    criteria = c("BIC_G", "AIC"),
    ...
) {
  level <- match.arg(level)

  available_criteria <- c("BIC_G", "BIC_N", "AIC", "AIC3")

  if (
    length(criteria) < 1L ||
    any(!criteria %in% available_criteria)
  ) {
    stop(
      "`criteria` must contain one or more of: ",
      paste(available_criteria, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  criteria <- unique(criteria)
  criterion_columns <- paste0(criteria, "_", level)

  old_par <- graphics::par(
    mfrow = c(1, length(criteria))
  )

  on.exit(
    graphics::par(old_par),
    add = TRUE
  )

  for (j in seq_along(criteria)) {
    graphics::plot(
      x = x$overview$nclus,
      y = x$overview[[criterion_columns[j]]],
      type = "b",
      xaxt = "n",
      xlab = "Number of SR-clusters",
      ylab = criteria[j],
      main = paste(criteria[j], "per model"),
      ...
    )

    graphics::axis(
      side = 1,
      at = x$overview$nclus
    )
  }

  invisible(
    x$overview[
      ,
      c("nclus", criterion_columns),
      drop = FALSE
    ]
  )
}
