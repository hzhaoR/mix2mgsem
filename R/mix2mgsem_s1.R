#' Run Step 1: measurement model clustering

#' Fit the first step of the double-mixture multigroup SEM method.
#' In the first step, groups are clustered based on (a subset of) measurement model
#' parameters using mixture multigroup factor analysis via `mixmgfa`.
#' The relevant parameters are extracted to obtain group-specific latent factor covariance matrices
#' that can be passed to [MixMix_Step2()] for the second step structural relations clustering.
#'
#' @param data A data frame containing the observed item variables and a group identifier.
#' @param step1model A lavaan model syntax string specifying the measurement model.
#' @param group Character string giving the name of the grouping variable in `data`. Defaults to `"group"`.
#' @param MM.cluster.spec Character vector passed to `mixmgfa::mixmgfa()`; specifies which measurement parameters are cluster-specific.
#' @param MM.nclus Integer scalar indicating the number of measurement clusters.
#' @param MM.maxiter Maximum number of iterations passed to `mixmgfa::mixmgfa()`.
#' @param MM.nruns Number of random starts passed to `mixmgfa::mixmgfa()`.
#' @param MM.design Design object passed to `mixmgfa::mixmgfa()`. Loading matrix (with ncol = nfactors) indicating position of zero loadings with '0' and non-zero loadings with '1'.
#' @param invar_loadings Invariant loadings specification passed to `mixmgfa::mixmgfa()`.
#' @param markers Marker variables specification passed to `mixmgfa::ScaleRotateMixmgfa_pinvar()`.
#' @param seed Integer seed used before fitting the measurement model.
#'
#' @return A named list containing:
#' \describe{
#'   \item{cov_eta}{An array of group-specific factor covariance matrices.}
#'   \item{ngroups}{Number of groups.}
#'   \item{N_gs}{Group sample sizes.}
#'   \item{S_unbiased}{Observed sample covariance matrices per group.}
#'   \item{vars}{Observed variable names.}
#'   \item{lat_var}{Latent variable names.}
#'   \item{mmgfa_output}{Selected outputs from the Step 1 `mixmgfa` fit.}
#'   \item{step1_time}{Elapsed Step 1 computation time in minutes.}
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' step1_model <- '
#'   F1 =~ x1 + x2 + x3 + x4 + x5
#'   F2 =~ z1 + z2 + z3 + z4 + z5
#'   F3 =~ m1 + m2 + m3 + m4 + m5
#'   F4 =~ y1 + y2 + y3 + y4 + y5
#' '
#'
#' step1_out <- MixMix_Step1(
#'   data = dat,
#'   step1model = step1_model,
#'   group = "country",
#'   MM.cluster.spec=c("loadings"),
#'   MM.nclus = 2,
#'   MM.design = design,
#'   invar_loadings = invar_loadings,
#'   markers = markers
#' )
#' }


MixMix_Step1 <- function(data, step1model, group = "group", MM.cluster.spec = c("loadings"),
                         MM.nclus, MM.maxiter = 10000, MM.nruns = 50,
                         MM.design, invar_loadings = NULL, markers, seed = 100){

  start_time_step1 <- Sys.time()

  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.", call. = FALSE)
  }

  if (!is.character(group) || length(group) != 1L || !group %in% names(data)) {
    stop("`group` must be a single column name in `data`.", call. = FALSE)
  }

  if (length(MM.nclus) != 1L) {
    stop(
      "This development version expects `MM.nclus` to be a single selected number of measurement model clusters.",
      call. = FALSE
    )
  }
  # centered data
  g_name <- as.character(unique(data[, group]))
  vars <- lavaan::lavNames(lavaan::lavaanify(step1model, auto = TRUE)) #observed var
  lat_var <- lavaan::lavNames(lavaan::lavaanify(step1model, auto = TRUE), "lv") #latent var

  centered <- data
  group.idx <- match(data[,group], g_name)
  group.sizes <- tabulate(group.idx)
  group.means <- rowsum.default(as.matrix(data[,vars]),
                                group = group.idx, reorder = FALSE,
                                na.rm = FALSE)/group.sizes
  centered[,vars] <- data[,vars] - group.means[group.idx, ,drop = FALSE]

  N_gs <- group.sizes
  nfactors <- length(lat_var)
  ngroups <- length(N_gs)

  # sample covariance matrix per group
  S_unbiased <- lapply(X = unique(centered[, group]), FUN = function(x) {stats::cov(centered[centered[, group] == x, vars])})

  set.seed(seed)
  # MixMG-CFA: loadings cluster-specific, 1-6 clusters;
  mixmgfa_args <- list(
    data = centered,
    N_gs = N_gs,
    nfactors = nfactors,
    cluster.spec = MM.cluster.spec,
    nsclust = MM.nclus,
    maxiter = MM.maxiter,
    nruns = MM.nruns,
    design = MM.design
  )

  if (!is.null(invar_loadings)) {
    if (!"invar_loadings" %in% names(formals(mixmgfa::mixmgfa))) {
      stop(
        "`invar_loadings` was supplied, but the installed version of ",
        "`mixmgfa::mixmgfa()` does not support this argument. ",
        "Use `invar_loadings = NULL` with the public version of `mixmgfa`, ",
        "or install a development version of `mixmgfa` that supports ",
        "`invar_loadings`.",
        call. = FALSE
      )
    }

    mixmgfa_args$invar_loadings <- invar_loadings
  }

  output1 <- do.call(mixmgfa::mixmgfa, mixmgfa_args)

  # MixMG-CFA: rescaling factors using marker variables
  output2 <- .scale_rotate_mixmgfa(output1, N_gs = N_gs, cluster.spec = MM.cluster.spec,
                                nsclust = MM.nclus, design = MM.design, rescale=1, markers = markers,
                                rotation=0,targetT=0,targetW=0)

  # the list of cluster names from MixMG-CFA solutions
  cluster_names <- grep("^\\d+\\.clusters$", names(output2[["MMGFAsolutions"]]), value = TRUE)
  uniquevar_key <- grep("uniquevariances$", names(output2[["MMGFAsolutions"]][[cluster_names]]), value = TRUE)
  # relevant parameters
  factor_cov_list<- output2[["MMGFAsolutions"]][[cluster_names]][["group.and.clusterspecific.factorcovariances"]]
  cluster_memb_list <- output2[["MMGFAsolutions"]][[cluster_names]][["clustermemberships"]]
  lambda_list <- output2[["MMGFAsolutions"]][[cluster_names]][["clusterspecific.loadings"]]
  theta_list <- output2[["MMGFAsolutions"]][[cluster_names]][[uniquevar_key]] #group or cluster

  # weighted sum approach to get group- and cluster-specific fcov
  weighted_sum_mg <- array(data = NA, dim = c(nfactors, nfactors, ngroups))

  for (j in 1:ngroups) {
      # multiply each cluster-specific matrix by its corresponding cluster membership and sum the results
      weighted_sum_mg[,,j] <- Reduce(`+`, Map(`*`, factor_cov_list[j, ], cluster_memb_list[j, ]))
    }
  cov_eta <- weighted_sum_mg
  dimnames(cov_eta)[[1]] <- dimnames(cov_eta)[[2]] <- lat_var

  end_time_step1 <- Sys.time()  # End time for Step 1
  step1_time <- difftime(end_time_step1, start_time_step1, units = "mins")

  return(list(cov_eta = cov_eta,
              ngroups = ngroups,
              N_gs = N_gs,
              S_unbiased   = S_unbiased,
              vars         = vars, #var names
              lat_var      = lat_var, #latent variable names
              mmgfa_output = list(output1 = output1,
                                  output2 = output2,
                                  cluster_memb = cluster_memb_list,
                                  factor_cov = factor_cov_list, #gk
                                  lambda_gs = lambda_list, #gk
                                  theta_gs = theta_list), #g
              step1_time = step1_time))

}
