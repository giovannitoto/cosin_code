# ---------------------------------------------------------------------------- #

get_metrics <- function(A, B, metrics = "all") {
  # -------------------------------------------------------------------------- #
  #   A : reference matrix
  #   B : data matrix
  # out : vector containing the metrics to maximize
  # -------------------------------------------------------------------------- #
  confusionMatrix_metrics <- c("Accuracy", "Sensitivity", "Specificity", "Precision", "F1")
  valid_metrics <- c("MSE", "RMSE", "MAE", confusionMatrix_metrics)
  if("all" %in% metrics) {
    metrics <- valid_metrics
  } else {
    metrics <- intersect(metrics, valid_metrics)
  }
  # output
  out <- rep(NA, length(metrics))
  names(out) <- metrics
  # MSE
  if("MSE" %in% metrics) out["MSE"] <- -1 * mean((A - B)^2)
  # RMSE
  if("RMSE" %in% metrics) out["RMSE"] <- -1 * sqrt(mean((A - B)^2))
  # MAE
  if("MAE" %in% metrics) out["MAE"] <- -1 * mean(abs(A - B))
  # note: there is -1 because we want to maximize the output values
  # caret::confusionMatrix
  if(any(metrics %in% confusionMatrix_metrics)) {
    # confusion matrix
    A <- factor(A == 0, levels = c("TRUE", "FALSE"))
    B <- factor(B == 0, levels = c("TRUE", "FALSE"))
    cm <- caret::confusionMatrix(data = B, reference = A)
    # accuracy
    if("Accuracy" %in% metrics) out["Accuracy"] <- cm$overall["Accuracy"]
    # sensibility
    if("Sensitivity" %in% metrics) out["Sensitivity"] <- cm$byClass["Sensitivity"]
    # specificity
    if("Specificity" %in% metrics) out["Specificity"] <- cm$byClass["Specificity"]
    # precision
    if("Precision" %in% metrics) out["Precision"] <- cm$byClass["Precision"]
    # F1
    if("F1" %in% metrics) out["F1"] <- cm$byClass["F1"]
  }
  # output
  out[is.na(out)] <- 0
  return(out)
}

# ---------------------------------------------------------------------------- #

metrics_MCMC <- function(out_MCMC, sim, METRICS, verbose = Inf) {
  metrics_tot <- matrix(-Inf, nrow = length(out_MCMC$numFactors),
                        ncol = length(METRICS) * (sim$k + 1))
  rownames(metrics_tot) <- 1:nrow(metrics_tot)
  colnames(metrics_tot) <- 1:ncol(metrics_tot)
  for (m in 1:length(METRICS)) {
    colnames(metrics_tot)[(m - 1) * (sim$k + 1) + 1:(sim$k + 1)] <- paste(METRICS[m], c(1:sim$k, "tot"), sep="_")
  }
  colnames(metrics_tot)
  n <- nrow(out_MCMC$eta[[1]])
  p <- nrow(out_MCMC$lambda[[1]])
  for (it in 1:length(out_MCMC$numFactors)) {
    # number of factors
    k <- ncol(out_MCMC$lambda[[it]])
    # get C_1, ..., C_k
    C_list <- list()
    for (h in 1:k) {
      if(sum(out_MCMC$lambda[[it]][, h]) != 0) {
        C_list <- append(C_list, list(tcrossprod(out_MCMC$eta[[it]][, h], out_MCMC$lambda[[it]][, h])))
      }
    }
    # number of active factors
    kstar <- length(C_list)
    # add matrices of zeros if kstar is less than the true number of factor (sim$k)
    if (kstar < sim$k) {
      C_list[[kstar + 1]] <- matrix(0, nrow = n, ncol = p)
    }
    # identify all possible permutations
    comb_list <- gtools::permutations(n = max(kstar, sim$k), r = sim$k)
    comb_list[comb_list > kstar] <- kstar + 1
    comb_list <- comb_list[!duplicated(comb_list), ]
    if(it %% verbose == 0) cat(it, ":", nrow(comb_list), "\n")
    # compute metrics for each possible permutation
    for (i in 1:nrow(comb_list)) {
      metrics <- sapply(1:sim$k, function(h) get_metrics(sim$C[[h]], C_list[[comb_list[i, h]]], metrics = METRICS))
      if(is.null(dim(metrics))) metrics <- matrix(metrics, nrow = 1, ncol = sim$k)
      metrics_sum <- apply(metrics, 1, sum)
      for (m in 1:length(METRICS)) {
        # save m-th metric if it is better than the previous ones
        if(metrics_sum[m] > metrics_tot[it, (m - 1) * (sim$k + 1) + sim$k + 1]) {
          metrics_tot[it, (m - 1) * (sim$k + 1) + sim$k + 1] <- metrics_sum[m]
          metrics_tot[it, (m - 1) * (sim$k + 1) + 1:sim$k] <- metrics[m, ]
        }
      }
    }
  }
  metrics_tot <- abs(metrics_tot)
  return(metrics_tot)
}

# ---------------------------------------------------------------------------- #

get_RMSE_glmpca <- function(mod, sim) {
  # -------------------------------------------------------------------------- #
  # mod : glmpca::glmpca output
  # sim : simulation
  # out : (sim$k+1)-dimensional vector
  # -------------------------------------------------------------------------- #
  out <- rep(Inf, sim$k + 1)
  # add columns to factors and loadings if they are less than sim$k
  if(ncol(mod$factors) < sim$k) {
    fact_diff <- sim$k - ncol(mod$factors)
    factors <- cbind(mod$factors, matrix(0, nrow = nrow(mod$factors), ncol = fact_diff))
    loadings <- cbind(mod$loadings, matrix(0, nrow = nrow(mod$loadings), ncol = fact_diff))
  } else {
    factors <- mod$factors
    loadings <- mod$loadings
  }
  # identify all possible permutations
  comb_list <- gtools::permutations(n = ncol(factors), r = sim$k)
  # compute metrics for each possible permutation
  for (i in 1:nrow(comb_list)) {
    metrics <- sapply(1:sim$k, function(h) sqrt(mean((tcrossprod(factors[,comb_list[i, h]], loadings[,comb_list[i, h]]) - sim$C[[h]])^2)))
    metrics <- c(metrics, sum(metrics))
    if(metrics[3] < out[3]) out <- metrics
  }
  return(out)
}

# ---------------------------------------------------------------------------- #
