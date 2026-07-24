#' Select the number of SR-clusters for Step 2
#'
#' Fit [MixMix_Step2()] for several numbers of SR-clusters and compare the model fit.
#'
#' @inheritParams MixMix_Step2
#' @param nclus Integer vector containing at least two candidate numbers of SR-clusters.
#' @param ... Additional arguments passed to [MixMix_Step2()].
#'
#' @return An object of class `mix2mgsem_step2_selection` containing:
#' \describe{
#'   \item{fits}{A list of [MixMix_Step2()] fits. Failed fits are stored as `NULL`.}
#'   \item{overview}{Data frame containing the fit and information criteria for the fitted models.}
#'   \item{selected}{Numbers of SR-clusters selected by each information criterion.}
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
MixMix_Step2_select <- function(s1out, step2model, nclus, ...) {

  nclus <- sort(unique(as.integer(nclus)))

  if (any(nclus > s1out$ngroups)) {
    stop(
      "Candidate values in `nclus` cannot exceed the number of groups (", s1out$ngroups, ").",
      call. = FALSE
    )
  }

  fits <- vector("list", length(nclus))
  names(fits) <- as.character(nclus)

  overview_rows <- vector("list", length(nclus))
  errors <- character(0)
  n_success <- 0L

  for (i in seq_along(nclus)) {
    k <- nclus[i]

    if (printing) {
      message("Fitting Step 2 model with nclus = ", k, ".")
    }

    s2out <- tryCatch(
      MixMix_Step2(
        s1out = s1out,
        step2model = step2model,
        nclus = k,
        nfree = nfree,
        seed = seed,
        # userStart = NULL,
        max_it = max_it,
        nstarts = nstarts,
        printing = printing,
        partition = partition,
        Endo2Cov = Endo2Cov,
        allG = allG
      ),
      error = function(e) e
    )

    # Continue with the other number of clusters when one model fails
    if (inherits(current_fit, "error")) {
      errors[as.character(k)] <- conditionMessage(s2out)

      if (printing) {
        message(
          "Step 2 model with nclus = ", k, " failed: ", errors[as.character(k)])
      }
      next
    }

    fits[[i]] <- s2out

    overview_rows[[n_success]] <- data.frame(
      nclus = k,
      logLik_observed = s2out$logLik$obs_loglik,
      nrpar_observed = s2out$NrPar$Obs.nrpar,
      BIC_G_observed = s2out$model_sel$BIC$observed$BIC_G,
      BIC_N_observed = s2out$model_sel$BIC$observed$BIC_N,
      AIC_observed = s2out$model_sel$AIC$observed,
      AIC3_observed = s2out$model_sel$AIC3$observed,
      logLik_factors = s2out$logLik$loglik,
      nrpar_factors = s2out$NrPar$Fac.nrpar,
      BIC_G_factors = s2out$model_sel$BIC$Factors$BIC_G,
      BIC_N_factors = s2out$model_sel$BIC$Factors$BIC_N,
      AIC_factors = s2out$model_sel$AIC$Factors,
      AIC3_factors = s2out$model_sel$AIC3$Factors,
        time_minutes = as.numeric(step2_time, units = "mins")
    )
  }

  overview <- do.call(rbind, overview_rows[seq_len(n_success)])
  rownames(overview) <- NULL

  if (n_success == 1L) {
    warning(
      "Only one model was fitted successfully; the information criteria cannot be meaningfully compared.",
      call. = FALSE
    )
  }

  # Return the number of clusters (min BIC/AIC)
  select_minimum <- function(column) {
    values <- overview[[column]]
    valid <- is.finite(values)

    if (!any(valid)) {return(integer(0))}

    minimum <- min(values[valid])

    overview$nclus[valid & values == minimum]
  }

  selected <- list(
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

  selected <- list(
    observed = list(
      BIC_G = select_minimum("BIC_G_observed"),
      BIC_N = select_minimum("BIC_N_observed"),
      AIC = select_minimum("AIC_observed"),
      AIC3 = select_minimum("AIC3_observed")
    ),

    factors = list(
      BIC_G = select_minimum("BIC_G_factors"),
      BIC_N = select_minimum("BIC_N_factors"),
      AIC = select_minimum("AIC_factors"),
      AIC3 = select_minimum("AIC3_factors")
    )
  )

  result <- list(
    fits = fits,
    overview = overview,
    selected = selected,
    errors = errors
  )

  class(result) <- "mix2mgsem_step2_selection"

  result
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
  cat("Observed data information criteria\n\n")

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

#' Plot Step 2 model selection criteria
#'
#' Plot information criterion values against the number of SR-clusters.
#'
#' @param x An object returned by [MixMix_Step2_select()].
#' @param level Character string indicating whether to plot criteria based on the observed data likelihood (`"observed"`) or factor likelihood (`"factors"`).
#' @param criteria Character vector specifying the information criteria to plot. Available options are `"BIC_G"`, `"BIC_N"`, `"AIC"`, and `"AIC3"`.
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
