# ---------------------------------------------------------------------------- #
# aGLM-PCA and GLM-PCA
# ---------------------------------------------------------------------------- #
# replications: R = 50
#   n     p  |           scenario
#  50   100  |  full + col1
# 200   100  |  full + row1
# 200  1000  |  full + row1 + col1
# ---------------------------------------------------------------------------- #

# R --slave --no-restore --file=simulations_glmpca.R

# ---------------------------------------------------------------------------- #

rm(list = ls())

# import functions from other R files
source("src/aglmpca.R")
source("src/metrics.R")

# save the files in a folder called "FOLDER"
FOLDER <- "simulations"

# max number of cores
maxCores <- 8

# ---------------------------------------------------------------------------- #

alpha_list <- 2:8

np_list <-  list(c(50, 100), c(200, 100), c(200, 1000))
sigma_list <- c(1/10, 1)
s_list <- c(1, 2, 3)
na_prop_list <- c(0.25)
R_list <- 1:50

# metrics_table: (n, p, sigma, r, best_mae, best_alpha, rmse_cont, mae_Y)

# ---------------------------------------------------------------------------- #

library(foreach)
numCores <- min(parallel::detectCores() - 1, maxCores)
myCluster <- parallel::makeCluster(numCores,type = "PSOCK")  # "FORK"
doParallel::registerDoParallel(cl = myCluster)

# ---------------------------------------------------------------------------- #
# aGLM-PCA and GLM-PCA WITHOUT METACOVARIATES
# ---------------------------------------------------------------------------- #

for(s in s_list) {
 # grid of all possible combinations of n, p, sigma, scenario
 sim_list <- expand.grid(r = R_list, na_prop = na_prop_list, sigma = sigma_list, np = np_list)
 # filename
 table_filename <- paste0("table_metrics_nometa_glmpca_s", s)
 # estimtae models and compute metrics
 table_metrics <- foreach(i = 1:nrow(sim_list), .combine = 'rbind', .inorder = FALSE) %dopar% {
   # -------------------------------------------------------------------------- #
   # get n, p, sigma, proportion of NA (na_prop) and replication (r)
   n <- sim_list$np[i][[1]][1]
   p <- sim_list$np[i][[1]][2]
   sigma <- sim_list$sigma[i]
   na_prop <- sim_list$na_prop[i]
   r <- sim_list$r[i]
   # -------------------------------------------------------------------------- #
   # file name
   sim_file <- paste("sim_s", s, "_n", n, "_p", p, "_sigma", sigma, "_r", r, "_na", na_prop, ".RDS", sep = "")
   # read data
   sim <- readRDS(file.path(FOLDER, sim_file))
   # -------------------------------------------------------------------------- #
   # best alpha and MAE
   best_alpha <- NA
   best_mae <- Inf
   # approximation of glm-pca (poi) without meta-covariates
   for(alpha in alpha_list) {
     Y_aglmpca <- aglmpca(t(sim$Y_na), L = alpha, pca = "ppca", completeY = TRUE)
     Y_aglmpca <- t(Y_aglmpca)
     mae <- mean(abs(sim$Y[is.na(sim$Y_na)] - Y_aglmpca[is.na(sim$Y_na)]))
     if(mae < best_mae) {
       best_alpha <- alpha
       best_mae <- mae
     }
   }
   # glm-pca (poi) without meta-covariates
   set.seed(28)
   mod <- glmpca::glmpca(t(sim$Y), L = best_alpha, fam = "poi", ctl=list("maxIter" = 3000))
   # metrics
   rmse_cont <- get_RMSE_glmpca(mod, sim)
   rmse_etaLambda <- sqrt(mean(tcrossprod(as.matrix(mod$factors), as.matrix(mod$loadings)) - tcrossprod(sim$eta, sim$lambda))^2)
   mae_Y <- mean(abs(t(predict(mod)) - sim$Y))
   # -------------------------------------------------------------------------- #
   # create row: (n, p, sigma, r, best_alpha, best_mae, mae_Y, rmse_cont, rmse_etaLambda, meta)
   met <- c(n, p, sigma, r, best_alpha, best_mae, mae_Y, rmse_cont, rmse_etaLambda, 0)
   # update txt
   cat(met, sep=",", file = file.path(FOLDER, "glmpca", paste0(table_filename,".txt")), append = TRUE)
   cat("\n",         file = file.path(FOLDER, "glmpca", paste0(table_filename,".txt")), append = TRUE)
   # update table
   met
   # -------------------------------------------------------------------------- #
 }
 # save metrics
 saveRDS(table_metrics, file.path(FOLDER, "glmpca", paste0(table_filename,".RDS")))
}

# ---------------------------------------------------------------------------- #
# GLM-PCA WITH METACOVARIATES
# ---------------------------------------------------------------------------- #

s_list <- setdiff(s_list, 2)  # there are no meta-covariates in scenario 2

for(s in s_list) {
  # grid of all possible combinations of n, p, sigma, scenario
  sim_list <- expand.grid(r = R_list, na_prop = na_prop_list, sigma = sigma_list, np = np_list)
  # filename
  table_filename <- paste0("table_metrics_meta_glmpca_s", s)
  # get best_alpha
  tm_wo <- paste0("table_metrics_nometa_glmpca_s", s)
  tm_wo <- as.data.frame(readRDS(file.path(FOLDER, "glmpca", paste0(tm_wo,".RDS"))))
  colnames(tm_wo)[1:6] <- c("n", "p", "sigma", "r", "alpha_a", "MAE_a")
  # estimtae models and compute metrics
  table_metrics <- foreach(i = 1:nrow(sim_list), .combine = 'rbind', .inorder = FALSE) %dopar% {
    # -------------------------------------------------------------------------- #
    # get n, p, sigma, proportion of NA (na_prop) and replication (r)
    n <- sim_list$np[i][[1]][1]
    p <- sim_list$np[i][[1]][2]
    sigma <- sim_list$sigma[i]
    na_prop <- sim_list$na_prop[i]
    r <- sim_list$r[i]
    # -------------------------------------------------------------------------- #
    # file name
    sim_file <- paste("sim_s", s, "_n", n, "_p", p, "_sigma", sigma, "_r", r, "_na", na_prop, ".RDS", sep = "")
    # read data
    sim <- readRDS(file.path(FOLDER, sim_file))
    # read MCMC and save meta-covariates
    wB <- readRDS(file.path(FOLDER, "results", paste0("s",s), paste0("res_meta_", sim_file)))$wB
    # -------------------------------------------------------------------------- #
    # get best alpha identified by aGLM-PCA without metacovariates
    best_alpha <- tm_wo[tm_wo$n==n & tm_wo$p==p & tm_wo$sigma==sigma & tm_wo$r==r, "alpha_a"]
    best_mae <- tm_wo[tm_wo$n==n & tm_wo$p==p & tm_wo$sigma==sigma & tm_wo$r==r, "MAE_a"]
    # glm-pca (poi) with meta-covariates
    set.seed(28)
    mod <- glmpca::glmpca(t(sim$Y), L = best_alpha, fam = "poi", Z = wB, ctl=list("maxIter" = 3000))
    # metrics
    rmse_cont <- get_RMSE_glmpca(mod, sim)
    rmse_etaLambda <- sqrt(mean(tcrossprod(as.matrix(mod$factors), as.matrix(mod$loadings)) - tcrossprod(sim$eta, sim$lambda))^2)
    mae_Y <- mean(abs(t(predict(mod)) - sim$Y))
    # -------------------------------------------------------------------------- #
    # create row: (n, p, sigma, r, best_alpha, best_mae, mae_Y, rmse_cont, rmse_etaLambda, meta)
    met <- c(n, p, sigma, r, best_alpha, best_mae, mae_Y, rmse_cont, rmse_etaLambda, 1)
    # update txt
    cat(met, sep=",", file = file.path(FOLDER, "glmpca", paste0(table_filename,".txt")), append = TRUE)
    cat("\n",         file = file.path(FOLDER, "glmpca", paste0(table_filename,".txt")), append = TRUE)
    # update table
    met
    # -------------------------------------------------------------------------- #
  }
  # save metrics
  saveRDS(table_metrics, file.path(FOLDER, "glmpca", paste0(table_filename,".RDS")))
}

# ---------------------------------------------------------------------------- #
