# Step 2: structural relations clustering
# adapted from Andres Felipe Perez Alonso's Mixture multigroup SEM code (mmgsem)

#' E-step for posterior structural relations cluster memberships
#'
#' Internal helper used by [MixMix_Step2()] to update posterior classification probabilities
#' from prior cluster probabilities and group-cluster log-likelihoods.
#'
#' @param pi_ks Numeric vector of prior structural relations cluster probabilities.
#' @param ngroup Integer number of groups.
#' @param nclus Integer number of structural relations clusters.
#' @param loglik Numeric matrix of group-cluster log-likelihood values.
#'
#' @return A numeric matrix of posterior cluster membership probabilities.
#' Rows correspond to groups and columns correspond to structural relations clusters.
#'
#' @noRd

EStep <- function(pi_ks, ngroup, nclus, loglik){
  max_g <-rep(0,ngroup)
  z_gks <- matrix(NA,nrow = ngroup,ncol = nclus)

  for(g in 1:ngroup){
    for(k in 1:nclus){
      z_gks[g,k] <- log(pi_ks[k])+loglik[g,k]
    }
    max_g[g] <- max(z_gks[g,])
    z_gks[g,] <- exp(z_gks[g,]-rep(max_g[g],nclus))
  }
  z_gks <- diag(1/apply(z_gks,1,sum))%*%z_gks

  return(z_gks)
}

#' Run Step 2: structural relations clustering
#'
#' Fit the second step of the double-mixture multigroup SEM method.
#' Given the group-specific factor covariance matrices obtained from
#' [MixMix_Step1()], groups are clustered according to their structural
#' relations among latent variables.
#'
#' @param s1out Output from [MixMix_Step1()].
#' @param step2model A lavaan model syntax string specifying the structural relations among latent variables.
#' @param nclus Integer scalar indicating the number of structural relations clusters.
#'@param nfree Integer scalar indicating the number of free non-marker loadings in the Step 1 measurement model.
#'   Used for observed-data information-criterion calculations.
#' @param seed Optional integer seed for random starts.
#' @param userStart Optional user-defined cluster membership matrix (dimensions: number of groups × number of clusters) with binary values
#'   (1 or 0; 1 indicates cluster membership). This can be used to initialize the EM algorithm when the user has prior information about cluster memberships.
#'   There must be only one 1 for each row. Skips random starts.
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
#' @return A named list containing:
#' \describe{
#'   \item{posteriors}{Posterior cluster membership probabilities.}
#'   \item{final_fit}{Final lavaan fit object (containing all group-cluster combinations).}
#'   \item{param}{Estimated group- and cluster-specific residual covariance matrices and cluster-specific regression matrices.}
#'   \item{logLik}{Log-likelihood values from the final solution and random starts.}
#'   \item{model_sel}{Model selection indices including BIC, AIC, and AIC3.}
#'   \item{NrPar}{Number of estimated parameters used for model selection calculations.}
#'   \item{step2_time}{Elapsed Step 2 computation time in minutes.}
#' }
#'
#' @details
#' This function implements an EM algorithm for structural relations clustering
#' using Step 1 factor covariance matrices as input. The implementation
#' is adapted from Andres Felipe Perez Alonso's mixture multigroup SEM code
#' (`mmgsem`) and extended for the two-step double-mixture multigroup SEM workflow.
#'
#' @export
#' @importFrom lavaan sem parTable lavaan lavInspect
#'
#' @examples
#' \dontrun{
#' step2_model <- '
#'   F4 ~ F1 + F3
#'   F3 ~ F1 + F2'
#'
#' step2_out <- MixMix_Step2(
#'   s1out = step1_out,
#'   step2model = step2_model,
#'   nclus = 2,
#'   nfree = 4,
#'   seed = 123
#' )
#' }

MixMix_Step2 <- function(s1out,  step2model, nclus, nfree,
                         seed = NULL, userStart = NULL, single = FALSE,
                         max_it = 10000L, nstarts = 50L, printing = FALSE,
                         partition = "hard", Endo2Cov = TRUE, allG = TRUE, fit = "factors"
){
  start_time_step2 <- Sys.time()

  required_s1 <- c("ngroups", "N_gs", "S_unbiased", "vars", "lat_var", "mmgfa_output")
  missing_s1 <- setdiff(required_s1, names(s1out))
  if (length(missing_s1) > 0L) {
    stop(
      "`s1out` is missing required components: ",
      paste(missing_s1, collapse = ", "),
      call. = FALSE
    )
  }

  required_mmgfa <- c("lambda_gs", "cluster_memb", "theta_gs")
  missing_mmgfa <- setdiff(required_mmgfa, names(s1out[["mmgfa_output"]]))
  if (length(missing_mmgfa) > 0L) {
    stop(
      "`s1out$mmgfa_output` is missing required components: ",
      paste(missing_mmgfa, collapse = ", "),
      call. = FALSE
    )
  }

  if (!is.numeric(nclus) || length(nclus) != 1L || is.na(nclus) || nclus < 1L) {
    stop("`nclus` must be a single positive integer.", call. = FALSE)
  }

  if (!is.numeric(nfree) || length(nfree) != 1L || is.na(nfree) || nfree < 0L) {
    stop("`nfree` must be a single non-negative integer.", call. = FALSE)
  }

  # from step 1
  ngroups   <- s1out$ngroups
  N_gs      <- s1out$N_gs
  lambda_ks <- s1out[["mmgfa_output"]][["lambda_gs"]]  # cluster-specific loadings
  s1_memb   <- s1out[["mmgfa_output"]][["cluster_memb"]]  # MM-cluster memberships
  theta_gs  <- s1out[["mmgfa_output"]][["theta_gs"]]
  S_unbiased <- s1out$S_unbiased
  vars      <- s1out$vars
  lat_var   <- s1out$lat_var
  nfactors  <- length(lat_var)

  # the other approach: single-indicator approach (not used for now)
  if(single==FALSE){
    cov_eta   <- s1out$cov_eta} else {cov_eta <- s1out$cov_eta_fs}
  cov_eta <- lapply(1:ngroups, function(i) cov_eta[,,i])

  # to avoid local maxima: multi-start
  # initialize objects
  results_nstarts <- vector(mode = "list", length = nstarts)
  z_gks_nstarts <- vector(mode = "list", length = nstarts)
  loglik_nstarts <- numeric(nstarts)

  if (!is.null(seed)) {
    set.seed(seed)
  }

  # function to reorder matrices (used later)
  reorder <- function(x) {
    x <- x[c(exog, endog), c(exog, endog)]
    return(x)
  }

  # Do a fake sem() to obtain the correct settings to use in Step 2
  # just a single sample cov!
  fake <- sem(
    model = step2model, sample.cov = rep(cov_eta[1], nclus),
    sample.nobs = rep(sum(N_gs), nclus), do.fit = FALSE,
    baseline = FALSE,
    h1 = FALSE, check.post = FALSE,
    loglik = FALSE,
    sample.cov.rescale = FALSE,
    fixed.x = TRUE
  )
  FakeprTbl <- parTable(fake)
  fake@Options$do.fit <- TRUE
  fake@Options$se     <- "none"
  fake@ParTable$start <- NULL
  fake@ParTable$est   <- NULL
  fake@ParTable$se    <- NULL
  fake@Options$start  <- "default"

  # Get the labels of the endogenous 1 and 2 factors
  endog1 <- lat_var[(lat_var %in% FakeprTbl$rhs[which(FakeprTbl$op == "~")]) &
                      (lat_var %in% FakeprTbl$lhs[which(FakeprTbl$op == "~")])]
  endog2 <- lat_var[!c(lat_var %in% FakeprTbl$rhs[which(FakeprTbl$op == "~")]) &
                      (lat_var %in% FakeprTbl$lhs[which(FakeprTbl$op == "~")])]
  endog <- c(endog1, endog2)
  exog <- lat_var[!c(lat_var %in% endog)]

  # endog2 <- lavNames(fake, "ov.y")
  # endog1 <- lavNames(fake, "eqs.y")
  # both.idx <- which(endog1 %in% endog2)
  # if (length(both.idx) > 0L) {
  #   endog1 <- endog1[-both.idx]
  # }
  # endog <- c(endog1, endog2)
  # exog <- lavNames(fake, "lv")[!c(lavNames(fake, "lv") %in% endog)]

  # Do a fake model per endo LV (to avoid bias due to reconstruction of group-specific endo variances)
  # Please note that "lv" is used to identify each model per endo LV
  # Please also note that this is used AFTER the first iteration. In the first iteration we start with only one model
  fake_lv <- vector(mode = "list", length = length(endog))

  for(lv in 1:length(endog)){
    # Create a parameter table per endogenous latent variables
    # Select the current latent variable
    this_lv <- endog[endog %in% endog[lv]]

    # Keep the (co)variances of the other latent variables
    var_not_this_lv <- which(FakeprTbl$lhs != this_lv & FakeprTbl$op == "~~")

    # Get the new parameter table per endogenous latent variable
    prTbl_lv <- FakeprTbl[c(which(FakeprTbl$lhs == this_lv), var_not_this_lv), ]

    # Run the model per endo latent variable
    fake_lv[[lv]] <- sem(
      model = prTbl_lv, sample.cov = rep(cov_eta[1], nclus),
      sample.nobs = rep(sum(N_gs), nclus), do.fit = FALSE,
      baseline = FALSE,
      h1 = FALSE, check.post = FALSE,
      loglik = FALSE,
      sample.cov.rescale = FALSE,
      fixed.x = TRUE
    )

    fake_lv[[lv]]@Options$do.fit <- TRUE
    fake_lv[[lv]]@ParTable$start <- NULL
    fake_lv[[lv]]@ParTable$est <- NULL
    fake_lv[[lv]]@ParTable$se <- NULL
    fake_lv[[lv]]@Options$start  <- "default"
  }

  # Re-order (order of columns and rows) cov_eta to make sure later computations are comparing correct matrices
  cov_eta <- lapply(1:ngroups, function(x) {reorder(cov_eta[[x]])})

  # Multi-start
  for (s in 1:nstarts) {
    tryCatch({
      if (printing == T){print(paste("Start", s, "-----------------"))}

      # Random Start
      if (!is.null(userStart)) {
        # In case the user inputs a pre-defined start, use it for z_gks
        z_gks <- userStart
      } else if (partition == "hard") {
        # Create initial random partition. Hard partition. (z_gks)
        cl <- 0
        while(cl < 1){ # "while loop" to make sure all clusters get at least one group
          z_gks <- t(replicate(ngroups, sample(x = c(rep(0, (nclus - 1)), 1))))
          cl <- min(colSums(z_gks))
        }
      } else if (partition == "soft") {
        z_gks <- matrix(data = runif(n = c(nclus*ngroups)), ncol = nclus, nrow = ngroups)
        z_gks <- z_gks/rowSums(z_gks)
      }

      # Initialize psi_gks
      psi_gks <- matrix(data = list(NA), nrow = ngroups, ncol = nclus)

      # Initialize weighted z_gks (weighted for each endo LV)
      # Done to avoid bias - NOT NECESSARY in the first iteration
      z_gks_lv <- vector(mode = "list", length = length(endog))
      for (lv in 1:length(endog)){
        z_gks_lv[[lv]] <- matrix(data = NA, ncol = nclus, nrow = ngroups)
      }

      # Prepare objects for the while loop
      i <- 0 # iteration initialization
      prev_LL <- 0 # previous loglikelihood initialization
      diff_LL <- 1 # Set a diff of 1 just to start the while loop
      log_test <- T # TEMPORARY - To check if there is decreasing loglikelihood

      # Run full-convergence multi-start
      while (diff_LL > 1e-6 & i < max_it & isTRUE(log_test)) {

        i <- i + 1
        pi_ks <- colMeans(z_gks) # Prior probabilities
        N_gks <- z_gks * N_gs # Sample size per group-cluster combination
        # N_gks <- as.vector(N_gks)

        # Weight posteriors for each endogenous factor - Not necessary in the first iteration
        # To avoid bias when allG == T. That is, when the endogenous variances are group-specific
        if (isTRUE(allG) & i > 1){
          for (lv in 1:length(endog)){
            for(k in 1:nclus){
              for(g in 1:ngroups){
                # Correct bias by dividing by the correct endo LV
                z_gks_lv[[lv]][g, k] <- z_gks[g, k]/psi_gks[[g, k]][endog[lv], endog[lv]]
              }
            }
          }
        }

        # M-Step --------

        # Trick to avoid slow multi-group estimation
        # Get a weighted averaged covariance matrix for each cluster
        if (isFALSE(allG) | i == 1){
          # Do this when allG is False OR when it is True and we are in the first iteration
          # For the first iteration there is no weighted z_gks
          COV <- vector("list", length = nclus)

          for (k in 1:nclus) {
            # create 'averaged' sample cov for this cluster
            this_nobs <- z_gks[,k] * N_gs
            this_w <- this_nobs/sum(this_nobs)
            tmp <- lapply(seq_along(cov_eta), function(g) {
              cov_eta[[g]] * this_w[g]
            })
            COV[[k]] <- Reduce("+", tmp)
            COV[[k]] <- 0.5 * (COV[[k]] + t(COV[[k]])) # Force the matrix to be symmetric
          }
        } else if (i > 1){
          # After the first iteration
          # Get one weighted cluster-specific COV per endo LV
          COV_lv <- vector("list", length = length(endog))
          for(lv in 1:length(endog)){
            COV_lv[[lv]] <- vector("list", length = nclus)
          }

          for(lv in 1:length(endog)){
            for (k in 1:nclus) {
              # create 'averaged' sample cov for this cluster
              this_nobs <- z_gks_lv[[lv]][, k] * N_gs
              this_w <- this_nobs/sum(this_nobs)
              tmp <- lapply(seq_along(cov_eta), function(g) {
                cov_eta[[g]] * this_w[g]
              })
              COV_lv[[lv]][[k]] <- Reduce("+", tmp)
              COV_lv[[lv]][[k]] <- 0.5 * (COV_lv[[lv]][[k]] + t(COV_lv[[lv]][[k]]))
            }
          }
        }

        # PARAMETER ESTIMATION
        # Call lavaan to estimate the structural parameters
        # the 'groups' are the clusters
        # Note: this makes all resulting parameters to be cluster-specific (it is reconstructed later)

        if (isFALSE(allG)| i == 1){
          # Do this when allG is False OR when it is True and we are in the first iteration
          # For the first iteration, perform the full structural model estimation
          s2out <- lavaan(slotOptions     = fake@Options,
                          slotParTable    = fake@ParTable,
                          sample.cov      = COV,
                          sample.nobs     = rep(sum(N_gs), nclus)
                          #slotModel       = slotModel,
                          #slotData        = fake@Data,
                          #slotSampleStats = fake@SampleStats
          )
        } else if (i > 1){
          # After the first iteration
          # Run structural estimation once per endo LV
          s2out <- vector(mode = "list", length = length(endog))
          for(lv in 1:length(endog)){
            s2out[[lv]] <- lavaan(slotOptions     = fake_lv[[lv]]@Options,
                                  slotParTable    = fake_lv[[lv]]@ParTable,
                                  sample.cov      = COV_lv[[lv]],
                                  sample.nobs     = rep(sum(N_gs), nclus)
                                  #slotModel       = slotModel,
                                  #slotData        = fake@Data,
                                  #slotSampleStats = fake@SampleStats
            )
          }
        }


        # start <- partable(s2out)$est

        # compute loglikelihood for all group/cluster combinations
        # Initialize matrices to store loglikelihoods
        loglik_gks  <- matrix(data = 0, nrow = ngroups, ncol = nclus)
        loglik_gksw <- matrix(data = 0, nrow = ngroups, ncol = nclus)
        gk <- 0

        # Prepare Sigma
        # Initialize the object for estimating Sigma
        Sigma <- matrix(data = list(NA), nrow = ngroups, ncol = nclus)
        I <- diag(length(lat_var)) # Identity matrix based on number of latent variables. Used later

        # Extract cluster-specific parameters from step 2
        if (isFALSE(allG) | i == 1){
          # Do this when allG is False OR when it is True and we are in the first iteration
          # In iteration one, there is only model (s2out) from which we extract the parameters.
          if (nclus == 1){
            EST_s2 <- lavInspect(s2out, "est", add.class = TRUE, add.labels = TRUE)
            beta_ks <- EST_s2[["beta"]]
            psi_ks <- EST_s2[["psi"]]
          } else if (nclus != 1){
            EST_s2 <- lavInspect(s2out, "est", add.class = TRUE, add.labels = TRUE)
            beta_ks <- lapply(EST_s2, "[[", "beta") # Does not work with only one cluster
            psi_ks <- lapply(EST_s2, "[[", "psi")
          }

          # Re order for correct comparisons
          if (nclus == 1){
            beta_ks <- reorder(beta_ks)
            psi_ks <- reorder(psi_ks)
          } else if (nclus != 1){
            beta_ks <- lapply(1:nclus, function(x) {reorder(beta_ks[[x]])}) # Does not work with only one cluster
            psi_ks <- lapply(1:nclus, function(x) {reorder(psi_ks[[x]])})
          }

        } else if (i > 1){
          # After iteration 1 we have several models (s2out[[lv]]) from which we extract the parameters
          # Extract the beta matrices per model (one per endo LV)
          # Initialize lists to store the parameters
          EST_s2_lv <- vector(mode = "list", length = length(endog))
          beta_ks_lv <- vector(mode = "list", length = length(endog))
          psi_ks_lv <- vector(mode = "list", length = length(endog))
          for (lv in 1:length(endog)){
            if (nclus == 1){
              EST_s2_lv[[lv]] <- lavInspect(s2out[[lv]], "est", add.class = TRUE, add.labels = TRUE)
              beta_ks_lv[[lv]] <- EST_s2_lv[[lv]][["beta"]]
              psi_ks_lv[[lv]] <- EST_s2_lv[[lv]][["psi"]]
            } else if (nclus != 1){
              EST_s2_lv[[lv]] <- lavInspect(s2out[[lv]], "est", add.class = TRUE, add.labels = TRUE)
              beta_ks_lv[[lv]] <- lapply(EST_s2_lv[[lv]], "[[", "beta") # Does not work with only one cluster
              psi_ks_lv[[lv]] <- lapply(EST_s2_lv[[lv]], "[[", "psi")
            }
          }

          # Combine all the beta matrices into just one per cluster

          # Start with an empty beta
          # beta <- matrix(data = 0, nrow = length(lat_var), ncol = length(lat_var))
          # colnames(beta) <- rownames(beta) <- lat_var
          # beta_ks <- lapply(X = seq_along(beta_ks), FUN = function(k){beta_ks[[k]] * 0})

          # beta_ks will contain the regression parameters per cluster
          # beta_ks_lv contains regressions per cluster AND per model of each endo latent variables

          for(k in 1:nclus){
            for(lv in 1:length(endog)){
              # Select current endogenous latent variable
              this_lv <- endog[lv]
              col.idx <- colnames(beta_ks_lv[[lv]][[k]])

              # Extract the regression coefficients of each endogenous latent variables
              if (nclus == 1){
                beta_ks[this_lv, col.idx] <- beta_ks_lv[[lv]][this_lv, col.idx]
              } else if (nclus != 1){
                beta_ks[[k]][this_lv, col.idx] <- beta_ks_lv[[lv]][[k]][this_lv, col.idx] # Does not work with only 1 cluster
              }
            }
          }

          # Re-order the beta matrices to make sure we are comparing the correct matrices
          if (nclus == 1){
            beta_ks <- reorder(beta_ks)
          } else if (nclus != 1){
            beta_ks <- lapply(1:nclus, function(x) {reorder(beta_ks[[x]])}) # Does not work with only one cluster
          }
        }

        for (k in 1:nclus) {
          # Previous matrices were only cluster-specific. We have to reconstruct the group-specific matrices (psi and sigma)

          ## Save the cluster-specific psi and beta
          ## ifelse() in case of only 1 cluster
          ifelse(test = (nclus == 1), yes = (psi <- psi_ks), no = (psi <- psi_ks[[k]]))
          ifelse(test = (nclus == 1), yes = (beta <- beta_ks), no = (beta <- beta_ks[[k]]))
          for (g in 1:ngroups) {
            # Reconstruct psi and sigma so they are group- and cluster-specific again.
            # Replace the group-specific part of psi
            # Exogenous (co)variance is always group-specific
            psi[exog, exog] <- cov_eta[[g]][exog, exog] # Replace the group-specific part

            # If the user required group-specific endogenous covariances (allG = T), do:
            if (allG == T){
              # Take into account the effect of the cluster-specific regressions
              # cov_eta[[g]] = solve(I - beta) %*% psi %*% t(solve(I - beta))
              # If we solve for psi, then:
              solved_psi <- ((I - beta) %*% cov_eta[[g]] %*% t((I - beta))) # Extract group-specific endog cov

              # Replace endog 1
              g_endog1_cov <- solved_psi[endog1, endog1]
              if(length(endog1) > 1){ # Remove cov between endog 1 variables
                g_endog1_cov[row(g_endog1_cov) != col(g_endog1_cov)] <- 0
              }
              psi[endog1, endog1] <- g_endog1_cov

              # Replace endog 2
              g_endog2_cov <- solved_psi[endog2, endog2] # Extract group-specific endog cov
              psi[endog2, endog2] <- g_endog2_cov
            }

            # If required by the user, set to 0 the covariance between endog 2 factors
            if (isFALSE(Endo2Cov)){
              psi[endog2, endog2][row(psi[endog2, endog2]) != col(psi[endog2, endog2])] <- 0
            }

            # Store for future check
            psi_gks[[g, k]] <- psi

            # Get log-likelihood by comparing factor covariance matrix of step 1 (cov_eta) and step 2 (Sigma)

            if (fit == "factors"){
              # Estimate Sigma (factor covariance matrix of step 2)
              Sigma[[g, k]] <- solve(I - beta) %*% psi %*% t(solve(I - beta))
              Sigma[[g, k]] <- 0.5 * (Sigma[[g, k]] + t(Sigma[[g, k]])) # Force to be symmetric

              # Estimate the loglikelihood
              loglik_gk <- lavaan:::lav_mvnorm_loglik_samplestats(
                sample.mean = rep(0, nrow(cov_eta[[1]])),
                sample.nobs = N_gs[g], # Use original sample size to get the correct loglikelihood
                # sample.nobs = N_gks[g, k],
                sample.cov  = cov_eta[[g]], # Factor covariance matrix from step 1
                Mu          = rep(0, nrow(cov_eta[[1]])),
                Sigma       = Sigma[[g, k]] # Factor covariance matrix from step 2
              )
              # } else if (fit == "observed"){
              #   # browser()
              #   S_biased <- S_unbiased[[g]] * (N_gs[[g]] - 1)/N_gs[[g]]
              #   Sigma[[g, k]] <- lambda_gs[[g]] %*% solve(I - beta) %*% psi %*% t(solve(I - beta)) %*% t(lambda_gs[[g]]) + theta_gs[[g]]
              #   Sigma[[g, k]] <- 0.5 * (Sigma[[g, k]] + t(Sigma[[g, k]]))
              #   # Sigma[[g, k]][lower.tri(Sigma[[g, k]])] <- t(Sigma[[g, k]])[lower.tri(Sigma[[g, k]])]
              #   loglik_gk <- lavaan:::lav_mvnorm_loglik_samplestats(
              #     sample.mean = rep(0, length(vars)),
              #     sample.nobs = N_gs[g], # Use original sample size to get the correct loglikelihood
              #     # sample.nobs = N_gks[g, k],
              #     sample.cov  = S_biased, # Item (observed) covariance matrix from step 1
              #     Mu          = rep(0, length(vars)),
              #     Sigma       = Sigma[[g, k]] # Item (observed) covariance matrix from step 2
              #   )

            }

            loglik_gks[g, k] <- loglik_gk
            loglik_gksw[g, k] <- log(pi_ks[k]) + loglik_gk # weighted loglik
          } # ngroups
        } # cluster

        # Get total loglikelihood
        # First, deal with arithmetic underflow by subtracting the maximum value per group
        max_gs <- apply(loglik_gksw, 1, max) # Get max value per row
        minus_max <- sweep(x = loglik_gksw, MARGIN = 1, STATS = max_gs, FUN = "-") # Subtract the max per row
        exp_loglik <- exp(minus_max) # Exp before summing for total loglikelihood
        loglik_gsw <- log(apply(exp_loglik, 1, sum)) # Sum exp_loglik per row and then take the log again
        LL <- sum((loglik_gsw + max_gs)) # Add the maximum again and then sum them all for total loglikelihood

        # Now, do E-step
        E_out <- EStep(
          pi_ks = pi_ks, ngroup = ngroups,
          nclus = nclus, loglik = loglik_gks
        )

        z_gks <- E_out
        diff_LL <- abs(LL - prev_LL)
        log_test <- prev_LL < LL | isTRUE(all.equal(prev_LL, LL))
        if (i == 1){log_test <- T}
        if (log_test == F){print(paste("Start", s, "; Iteration", i, "-------")); print(paste("Difference", LL - prev_LL))}#; browser()}
        log_test <- T
        prev_LL <- LL
        if (printing == T){print(i); print(LL)}
      }

      results_nstarts[[s]] <- s2out
      z_gks_nstarts[[s]] <- z_gks
      loglik_nstarts[s] <- LL

    }, error = function(e) {
      # Code to handle the error, such as printing a message or logging the error
      cat("Error in nstarts", s, ":", e$message, "\n")
    })
  }

  # Get best fit and z_gks based on the loglikelihood
  best_idx <- which.max(loglik_nstarts)
  s2out <- results_nstarts[[best_idx]]
  LL <- loglik_nstarts[best_idx]
  z_gks <- z_gks_nstarts[[best_idx]]
  colnames(z_gks) <- paste("Cluster", seq_len(nclus))

  # Extract matrices from final step 2 output
  if (isFALSE(allG)){
    if (nclus == 1){
      EST_s2 <- lavInspect(s2out, "est", add.class = TRUE, add.labels = TRUE) # Estimated matrices step 2
      beta_ks <- EST_s2[["beta"]]
      psi_ks <- EST_s2[["psi"]]
    } else if (nclus != 1){
      EST_s2 <- lavInspect(s2out, "est", add.class = TRUE, add.labels = TRUE) # Estimated matrices step 2
      beta_ks <- lapply(EST_s2, "[[", "beta")
      psi_ks <- lapply(EST_s2, "[[", "psi")
    }
  } else if (isTRUE(allG)){
    EST_s2_lv <- vector(mode = "list", length = length(endog))
    beta_ks_lv <- vector(mode = "list", length = length(endog))
    psi_ks_lv <- vector(mode = "list", length = length(endog))
    for (lv in 1:length(endog)){
      if (nclus == 1){
        EST_s2_lv[[lv]] <- lavInspect(s2out[[lv]], "est", add.class = TRUE, add.labels = TRUE)
        beta_ks_lv[[lv]] <- EST_s2_lv[[lv]][["beta"]]
        psi_ks_lv[[lv]] <- EST_s2_lv[[lv]][["psi"]]
      } else if (nclus != 1){
        EST_s2_lv[[lv]] <- lavInspect(s2out[[lv]], "est", add.class = TRUE, add.labels = TRUE)
        beta_ks_lv[[lv]] <- lapply(EST_s2_lv[[lv]], "[[", "beta") # Does not work with only one cluster
        psi_ks_lv[[lv]] <- lapply(EST_s2_lv[[lv]], "[[", "psi")
      }
    }

    # Re-construct beta_ks
    for(k in 1:nclus){
      for(lv in 1:length(endog)){
        this_lv <- endog[lv]
        col.idx <- colnames(beta_ks_lv[[lv]][[k]])
        if (nclus == 1){
          beta_ks[this_lv, col.idx] <- beta_ks_lv[[lv]][this_lv, col.idx]
        } else if (nclus != 1){
          beta_ks[[k]][this_lv, col.idx] <- beta_ks_lv[[lv]][[k]][this_lv, col.idx]
        }
      }
    }
  }

  # Re-order betas
  if (nclus == 1){
    beta_ks <- reorder(beta_ks)
  } else if (nclus != 1){
    beta_ks <- lapply(1:nclus, function(x) {reorder(beta_ks[[x]])}) # Does not work with only one cluster
  }

  # Get the group- and cluster-specific psi_gks matrices
  psi_gks <- matrix(data = list(NA), nrow = ngroups, ncol = nclus)

  for (k in 1:nclus) {
    ifelse(test = (nclus == 1), yes = (psi_k <- psi_ks), no = (psi_k <- psi_ks[[k]]))
    ifelse(test = (nclus == 1), yes = (beta <- beta_ks), no = (beta <- beta_ks[[k]]))
    for (g in 1:ngroups) {
      psi_gks[[g, k]] <- psi_k
      psi_gks[[g, k]][exog, exog] <- cov_eta[[g]][exog, exog]

      # If the user required group-specific endogenous covariances (allG = T), do:
      if (allG == T){
        # Take into account the effect of the cluster-specific regressions
        # cov_eta[[g]] = solve(I - beta) %*% psi %*% t(solve(I - beta))
        # If we solve for psi, then:
        g_endog1_cov <- ((I - beta) %*% cov_eta[[g]] %*% t((I - beta)))[endog1, endog1] # Extract group-specific endog cov
        if(length(endog1) > 1){
          g_endog1_cov[row(g_endog1_cov) != col(g_endog1_cov)] <- 0
        }
        #browser()
        psi_gks[[g, k]][endog1, endog1] <- g_endog1_cov

        g_endog2_cov <- ((I - beta) %*% cov_eta[[g]] %*% t((I - beta)))[endog2, endog2] # Extract group-specific endog cov
        psi_gks[[g, k]][endog2, endog2] <- g_endog2_cov
      }

      # If required by the user, endog2 covariances are set to 0
      if (isFALSE(Endo2Cov)){
        offdiag <- row(psi_gks[[g, k]][endog2, endog2]) != col(psi_gks[[g, k]][endog2, endog2])
        psi_gks[[g, k]][endog2, endog2][offdiag] <- 0
      }
    } # groups
  } # cluster

  # MODEL SELECTION ------------------------------------------------
  # Get observed data log-likelihood using Kim's code (for model selection purposes)
  H <- length(lambda_ks)
  K <- nclus
  Sigma_ghk <- array(vector("list", ngroups * H * K),
                     dim = c(ngroups, H, K))  # covariance per group × MM × SR
  Obs.loglik_ghk <- array(0, dim = c(ngroups, H, K))     # loglik per group × MM × SR
  Obs.loglik_ghkw <- array(0, dim = c(ngroups, H, K))    # weighted loglik per group × MM × SR
  pi_ks <- colMeans(z_gks)
  pi_hs <- colMeans(s1_memb)
  for (g in 1:ngroups) {
    S_biased <- S_unbiased[[g]] * (N_gs[[g]] - 1) / N_gs[[g]]

    for (h in 1:H) {   # loop over MM-clusters
      lambda_g_h <- lambda_ks[[h]]  # factor loadings for MM-cluster h
      psi_g_h <- theta_gs[[g]]        # residuals for MM-cluster h

      for (k in 1:K) { # loop over SR-clusters
        if (K == 1) {
          beta <- beta_ks
        } else {
          beta <- beta_ks[[k]]
        }
        # latent variance for SR cluster
        var_eta <- solve(I - beta) %*% psi_gks[[g, k]] %*% t(solve(I - beta))

        # observed covariance for this MM × SR combination
        Sigma_ghk[[g, h, k]] <- lambda_g_h %*% var_eta %*% t(lambda_g_h) + psi_g_h

        # compute loglik for this combination
        Obs.loglik_ghk[g, h, k] <- lavaan:::lav_mvnorm_loglik_samplestats(
          sample.mean = rep(0, length(vars)),
          sample.nobs  = N_gs[g],
          sample.cov  = S_biased,
          Mu          = rep(0, length(vars)),
          Sigma       = Sigma_ghk[[g, h, k]]
        )

        # add mixture weights
        Obs.loglik_ghkw[g, h, k] <- log(pi_hs[h]) + log(pi_ks[k]) + Obs.loglik_ghk[g, h, k]
      }
    }
  }


  # Get total observed loglikelihood
  Obs.LL_g <- numeric(ngroups)
  for (g in 1:ngroups) {
    # flatten MM × SR combinations
    temp_mat <- as.vector(Obs.loglik_ghkw[g,,])
    # subtract max per group for numerical stability
    max_val <- max(temp_mat)
    Obs.LL_g[g] <- log(sum(exp(temp_mat - max_val))) + max_val
  }
  # total log-likelihood
  Obs.LL <- sum(Obs.LL_g)

  Q <- length(lat_var)
  J <- length(vars)

  # Structural parameters
  ifelse(test = (nclus == 1), yes = (n_reg <- sum(beta_ks != 0)), no = (n_reg <- sum(beta_ks[[1]] != 0)))
  Q_exo <- length(exog)
  n_cov_exo <- ((Q_exo * (Q_exo + 1)) / 2)
  Q_endo1 <- length(endog1)
  Q_endo2 <- length(endog2)
  n_cov_endo2 <- ((Q_endo2 * (Q_endo2 + 1)) / 2)

  # Measurement parameters
  n_res <- sum(theta_gs[[1]] != 0) #residuals as group-specific
  n_load <- sum(lambda_ks[[1]] != 0)
  n_free <- as.integer(nfree)

  # Get the correct number of free parameters depending on the possible combinations
  # browser()
  if (allG == F){ # Is endogenous covariance group-specific?
    nr_par_factors <- (nclus - 1) + (n_reg * nclus) + (n_cov_exo * ngroups) + (Q_endo1 * nclus) + (n_cov_endo2 * nclus)
    nr_pars <- (nclus - 1) + (n_reg * nclus) + (n_cov_exo * ngroups) + (Q_endo1 * nclus) + (n_cov_endo2 * nclus) + (n_res * ngroups) + (n_load - Q - n_free) + (n_free * length(lambda_ks)) #group-specific uniqvar + invariant loadings + cluster-specific loadings
    # nr_pars_cl <- (nclus - 1) + (n_reg * nclus) + (n_cov_exo * ngroups) + (Q_endo1 * nclus) + (n_cov_endo2 * nclus) + (n_res * ngroups) + (n_load - Q - n_free) + (n_free*2* length(lambda_ks)) #group-specific uniqvar + invariant loadings + cluster-specific loadings + crossloadings

  } else if (allG == T){
    nr_par_factors <- (nclus - 1) + (n_reg * nclus) + (n_cov_exo * ngroups) + (Q_endo1 * ngroups) + (n_cov_endo2 * ngroups)
    nr_pars <- (nclus - 1) + (n_reg * nclus) + (n_cov_exo * ngroups) + (Q_endo1 * ngroups) + (n_cov_endo2 * ngroups) + (n_res * ngroups) + (n_load - Q - n_free) + (n_free * length(lambda_ks))
    # nr_pars_cl <- (nclus - 1) + (n_reg * nclus) + (n_cov_exo * ngroups) + (Q_endo1 * nclus) + (n_cov_endo2 * nclus) + (n_res * ngroups) + (n_load - Q - n_free) + (n_free*2* length(lambda_ks)) #group-specific uniqvar + invariant loadings + cluster-specific loadings + crossloadings
  }

  # Calculate BIC
  # Factors
  BIC_N <- -2 * LL + nr_par_factors * log(sum(N_gs))
  BIC_G <- -2 * LL + nr_par_factors * log(ngroups)

  # Observed LL
  Obs.BIC_N <- -2 * Obs.LL + nr_pars * log(sum(N_gs)) #parameters only step 2 (factors)
  Obs.BIC_G <- -2 * Obs.LL + nr_pars * log(ngroups)


  # Calculate AIC (and AIC3).
  # Factors
  AIC  <- (-2 * LL) + (nr_par_factors * 2)
  AIC3 <- (-2 * LL) + (nr_par_factors * 3)

  # Observed
  Obs.AIC  <- (-2 * Obs.LL) + (nr_pars * 2)
  Obs.AIC3 <- (-2 * Obs.LL) + (nr_pars * 3)

  # Re order matrices so that we get them in the following order:
  # (1) Exogenous latent variables
  # (2) Endogenous latent variables: independent and dependent variables at the same time
  # (3) Endogenous latent variables: only dependent variables

  # Reoder psi_ks and beta_ks by using the reorder function in the lapply function
  gro_clu <- ngroups * nclus
  psi_gks <- array(lapply(1:gro_clu, function(x) {reorder(psi_gks[[x]])}), dim = c(ngroups, nclus))
  if (nclus == 1){
    beta_ks <- reorder(beta_ks)
  } else if (nclus != 1){
    beta_ks <- lapply(1:nclus, function(x) {reorder(beta_ks[[x]])}) # Does not work with only one cluster
  }

  names(beta_ks) <- paste("Cluster", seq_len(nclus))

  end_time_step2 <- Sys.time()
  step2_time <- difftime(end_time_step2, start_time_step2, units = "mins")

  return(list(
    posteriors    = z_gks,
    final_fit     = s2out, # Final fit of step 2 (contains all group-cluster combinations)
    param         = list(psi_gks = psi_gks, beta_ks = beta_ks),
    logLik        = list(loglik        = LL, # Final logLik of the model (its meaning depends on arguments "fit" and "est_method")
                         loglik_gksw   = loglik_gksw, # Weighted logLik per group-cluster combinations
                         runs_loglik   = loglik_nstarts, # loglik for each start
                         obs_loglik    = Obs.LL),
    model_sel     = list(BIC        = list(observed = list(BIC_N = Obs.BIC_N, BIC_G = Obs.BIC_G),
                                           Factors = list(BIC_N = BIC_N, BIC_G = BIC_G)),
                         AIC        = list(observed = Obs.AIC, Factors = AIC),
                         AIC3       = list(observed = Obs.AIC3, Factors = AIC3)),
    NrPar         = list(Obs.nrpar = nr_pars, Fac.nrpar = nr_par_factors),
    step2_time    = step2_time
  ))
}

