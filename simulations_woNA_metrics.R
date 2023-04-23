# ---------------------------------------------------------------------------- #

# R --slave --no-restore --file=simulations_metrics.R

# ---------------------------------------------------------------------------- #

rm(list = ls())
source("src/metrics.R")

# ---------------------------------------------------------------------------- #

# save the files in a folder called "FOLDER"
FOLDER <- "simulations"

# max number of cores
maxCores <- 8

# ---------------------------------------------------------------------------- #

METRICS <- list(c("RMSE", "F1"), "RMSE", "RMSE")

np_list <-  list(c(50, 100), c(200, 100), c(200, 1000))
sigma_list <- c(1/10, 1)
s_list <- c(3, 2, 1)
R_list <- 1:50

# ---------------------------------------------------------------------------- #

library(foreach)
numCores <- min(parallel::detectCores() - 1, maxCores)
myCluster <- parallel::makeCluster(numCores,type = "PSOCK")  # "FORK"
doParallel::registerDoParallel(cl = myCluster)

# ---------------------------------------------------------------------------- #
# MCMC WITHOUT META-COVARIATES
# ---------------------------------------------------------------------------- #

for(s in s_list) {
  # grid of all possible combinations of n, p, sigma, scenario
  sim_list <- expand.grid(r = R_list, sigma = sigma_list, np = np_list)
  # filename
  table_filename <- paste0("table_metrics_nometa_woNA_s", s)
  # compute metrics
  table_metrics <- foreach(i = 1:nrow(sim_list), .combine = 'rbind', .inorder = FALSE) %dopar% {
    # ------------------------------------------------------------------------ #
    # get n, p, sigma, scenario (s), replication (r), seed
    n <- sim_list$np[i][[1]][1]
    p <- sim_list$np[i][[1]][2]
    sigma <- sim_list$sigma[i]
    r <- sim_list$r[i]
    # ------------------------------------------------------------------------ #
    sim_file <- paste("sim_s", s, "_n", n, "_p", p, "_sigma", sigma, "_r", r, "_woNA.RDS", sep = "")
    sim_load <- paste("sim_s", s, "_n", n, "_p", p, "_sigma", sigma, "_r", r, "_na0.25.RDS", sep = "")
    sim <- readRDS(file.path(FOLDER, sim_load))
    etaLambda_sim <- tcrossprod(sim$eta, sim$lambda)
    # ------------------------------------------------------------------------ #
    # MCMC with meta-covariates (scenario 1 and 3)
    # MCMC
    out_MCMC <- readRDS(file.path(FOLDER, "results", paste0("s", s), paste0("res_nometa_", sim_file)))
    sp <- length(out_MCMC$numFactors)
    # Metrics
    metrics_meta <- metrics_MCMC(out_MCMC, sim, METRICS[[s]], verbose = Inf)
    metrics_meta <- apply(metrics_meta, 2, mean)
    # kstar
    kstar <- mean(out_MCMC$numFactors)
    # Mean Absolute Deviation (MAD), MAE & RMSE
    rmse_etaLambda <- 0
    mae_Y <- 0
    mad_Y <- 0
    for (it in 1:sp) {
      etaLambda_it <- tcrossprod(out_MCMC$eta[[it]], out_MCMC$lambda[[it]])
      Y_it <- pmin(floor(exp(etaLambda_it)), sim$hyperparameters$y_max)
      rmse_etaLambda <- rmse_etaLambda + sqrt(mean((etaLambda_it - etaLambda_sim)^2)) / sp
      mae_Y <- mae_Y + mean(abs(Y_it - sim$Y)) / sp
      mad_Y <- mad_Y + mean(abs(Y_it - mean(Y_it))) / sp
    }
    # Benchmark
    sd_cont <- sapply(sim$C, sd)
    sd_etaLambda <- sd(tcrossprod(sim$eta, sim$lambda))
    # create row
    met <- c(n, p, sigma, s, r, metrics_meta, mad_Y, mae_Y, rmse_etaLambda, kstar, sd_cont, sd_etaLambda, 0)
    # update txt
    cat(met, sep=",", file=paste0(table_filename,".txt"), append = TRUE)
    cat("\n", file=paste0(table_filename,".txt"), append = TRUE)
    # update table
    met
    # ------------------------------------------------------------------------ #
  }
  # save metrics
  saveRDS(table_metrics, file.path(FOLDER, "metrics", paste0(table_filename,".RDS")))
  # -------------------------------------------------------------------------- #
}

# ---------------------------------------------------------------------------- #
# MCMC WITH META-COVARIATES
# ---------------------------------------------------------------------------- #

s_list <- setdiff(s_list, 2)  # there are no meta-covariates in scenario 2

for(s in s_list) {
  # grid of all possible combinations of n, p, sigma, scenario
  sim_list <- expand.grid(r = R_list, sigma = sigma_list, np = np_list)
  # filename
  table_filename <- paste0("table_metrics_meta_woNA_p1000_s", s)
  # compute metrics
  table_metrics <- foreach(i = 1:nrow(sim_list), .combine = 'rbind', .inorder = FALSE) %dopar% {
    # ------------------------------------------------------------------------ #
    # get n, p, sigma, scenario (s), replication (r), seed
    n <- sim_list$np[i][[1]][1]
    p <- sim_list$np[i][[1]][2]
    sigma <- sim_list$sigma[i]
    r <- sim_list$r[i]
    # ------------------------------------------------------------------------ #
    sim_file <- paste("sim_s", s, "_n", n, "_p", p, "_sigma", sigma, "_r", r, "_woNA.RDS", sep = "")
    sim_load <- paste("sim_s", s, "_n", n, "_p", p, "_sigma", sigma, "_r", r, "_na0.25.RDS", sep = "")
    sim <- readRDS(file.path(FOLDER, sim_load))
    etaLambda_sim <- tcrossprod(sim$eta, sim$lambda)
    # ------------------------------------------------------------------------ #
    # MCMC with meta-covariates (scenario 1 and 3)
    # MCMC
    out_MCMC <- readRDS(file.path(FOLDER, "results", paste0("s", s), paste0("res_meta_", sim_file)))
    sp <- length(out_MCMC$numFactors)
    # Metrics
    metrics_meta <- metrics_MCMC(out_MCMC, sim, METRICS[[s]], verbose = Inf)
    metrics_meta <- apply(metrics_meta, 2, mean)
    # kstar
    kstar <- mean(out_MCMC$numFactors)
    # Mean Absolute Deviation (MAD), MAE & RMSE
    rmse_etaLambda <- 0
    mae_Y <- 0
    mad_Y <- 0
    for (it in 1:sp) {
      etaLambda_it <- tcrossprod(out_MCMC$eta[[it]], out_MCMC$lambda[[it]])
      Y_it <- pmin(floor(exp(etaLambda_it)), sim$hyperparameters$y_max)
      rmse_etaLambda <- rmse_etaLambda + sqrt(mean((etaLambda_it - etaLambda_sim)^2)) / sp
      mae_Y <- mae_Y + mean(abs(Y_it - sim$Y)) / sp
      mad_Y <- mad_Y + mean(abs(Y_it - mean(Y_it))) / sp
    }
    # Benchmark
    sd_cont <- sapply(sim$C, sd)
    sd_etaLambda <- sd(tcrossprod(sim$eta, sim$lambda))
    # create row
    met <- c(n, p, sigma, s, r, metrics_meta, mad_Y, mae_Y, rmse_etaLambda, kstar, sd_cont, sd_etaLambda, 1)
    # update txt
    cat(met, sep=",", file=paste0(table_filename,".txt"), append = TRUE)
    cat("\n", file=paste0(table_filename,".txt"), append = TRUE)
    # update table
    met
    # ------------------------------------------------------------------------ #
  }
  # save metrics
  saveRDS(table_metrics, file.path(FOLDER, "metrics", paste0(table_filename,".RDS")))
  # -------------------------------------------------------------------------- #
}

parallel::stopCluster(cl = myCluster)

# ---------------------------------------------------------------------------- #
