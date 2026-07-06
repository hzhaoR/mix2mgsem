# Step 1: measurement model clustering
MixMix_Step1 <- function(data, step1model, group = "group", MM.cluster.spec=c("loadings"),
                         MM.nclus = c(1:6), MM.maxiter = 10000, MM.nruns = 50, 
                         MM.design, invar_loadings, markers, seed = 100){
  
  start_time_step1 <- Sys.time()  
  
  # centered data
  g_name <- as.character(unique(data[, group]))
  vars <- lavaan::lavNames(lavaanify(step1model, auto = TRUE)) #observed var
  lat_var <- lavaan::lavNames(lavaanify(step1model, auto = TRUE), "lv") #latent var
  
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
  S_unbiased <- lapply(X = unique(centered[, group]), FUN = function(x) {cov(centered[centered[, group] == x, vars])})
  
  set.seed(seed)
  # MixMG-CFA: loadings cluster-specific, 1-6 clusters;
  output1 <- mixmgfa::mixmgfa(data = centered, N_gs = N_gs, nfactors = nfactors, 
                    cluster.spec = MM.cluster.spec, nsclust = MM.nclus, 
                    maxiter = MM.maxiter, nruns = MM.nruns, design = MM.design, invar_loadings = invar_loadings)
  
  # MixMG-CFA: rescaling factors using marker variables
  output2 <- mixmgfa::ScaleRotateMixmgfa_pinvar(output1, N_gs = N_gs,
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
