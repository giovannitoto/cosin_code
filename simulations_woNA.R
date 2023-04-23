# ---------------------------------------------------------------------------- #
# replications: R = 50
#   n     p  |           scenario
#  50   100  |  full + col1
# 200   100  |  full + row1
# 200  1000  |  full + row1 + col1
# ---------------------------------------------------------------------------- #

# R --slave --no-restore --file=simulations_woNA.R

# ---------------------------------------------------------------------------- #

rm(list = ls())

# ---------------------------------------------------------------------------- #

get_output <- function(lambda, eta, wB, hyper, na_prop) {
  n <- nrow(eta)
  p <- nrow(lambda)
  k <- ncol(lambda)
  C_list <- list()
  for (h in 1:k) {
    C_list[[h]] <- eta[, h] %*% t(lambda[, h])
  }
  Z <- eta %*% t(lambda) + matrix(rnorm(n * p, 0, hyper$sigma), nrow = n, ncol = p)
  y <- pmin(floor(exp(Z)), hyper$y_max)
  output <- list("k" = k, "lambda" = lambda, "eta" = eta, "C" = C_list,
                 "Y" = y, "wB" = wB, "hyperparameters" = hyper)
  return(output)
}

# ---------------------------------------------------------------------------- #

# save the files in a folder called "FOLDER"
FOLDER <- "simulations"

# max number of cores
maxCores <- 8

# ---------------------------------------------------------------------------- #

y_max <- Inf
sigma_lambda <- 1
sigma_eta <- 1
sigma0 <- 0.05
alpha <- 5

hyper <- list("y_max" = y_max, "alpha" = alpha, "sigma_lambda" = sigma_lambda,
              "sigma_eta" = sigma_eta, "sigma0" = sigma0)

np_list <-  list(c(50, 100), c(200, 100), c(200, 1000))
sigma_list <- c(1/10, 1)
s_list <- c(1, 2, 3)
na_prop_list <- c(0.25)
R_list <- 1:50

# grid of all possible combinations of n, p, sigma, scenario
sim_list <- expand.grid(r = R_list, na_prop = na_prop_list, sigma = sigma_list, np = np_list, s = s_list)
# seeds of the simulations
sim_list$seed <- 1:nrow(sim_list) + 5000

# ---------------------------------------------------------------------------- #

library(foreach)
numCores <- min(parallel::detectCores() - 1, maxCores)
myCluster <- parallel::makeCluster(numCores,type = "PSOCK")  # "PSOCK", "FORK"
doParallel::registerDoParallel(cl = myCluster)

foreach(i = 1:nrow(sim_list), .combine = 'rbind', .inorder = FALSE) %dopar% {
  # -------------------------------------------------------------------------- #
  # get n, p, sigma, scenario (s), replication (r), seed
  n <- sim_list$np[i][[1]][1]
  p <- sim_list$np[i][[1]][2]
  sigma <- sim_list$sigma[i]
  s <- sim_list$s[i]
  na_prop <- sim_list$na_prop[i]
  r <- sim_list$r[i]
  seed <- sim_list$seed[i]
  # -------------------------------------------------------------------------- #
  # file name
  sim_file <- paste("sim_s", s, "_n", n, "_p", p, "_sigma", sigma, "_r", r, ".RDS", sep = "")
  # hyperparameters
  hyper$sigma <- sigma
  # set seed
  set.seed(seed)
  # ---------------------------------------------------------------------- #
  if(s == 1) {
    # lambda
    l1 <- rnorm(p, 0, sigma_lambda)
    l2 <- c(rep(1, p/2), rep(0, p/2))
    lambda <- cbind(l1, l2)
    # eta
    e1 <- rnorm(n, 0, sigma_eta)
    e2 <- rnorm(n, 0, sigma_eta)
    eta <- cbind(e1, e2)
    # meta-covariates (see lambda)
    wB1 <- c(rep(1, p/2), rep(0, p/2))
    wB <- data.frame("wB1" = as.factor(wB1))
  } else if(s == 2) {
    # lambda
    l1 <- rnorm(p, 0, sigma_lambda)
    l2 <- rnorm(p, 0, sigma_lambda)
    lambda <- cbind(l1, l2)
    # eta
    e1 <- rnorm(n, 0, sigma_eta)
    e2 <- c(rep(1, n/2), rnorm(n/2, 0, sigma0))
    eta <- cbind(e1, e2)
    # meta-covariates (see lambda)
    wB <- NULL
  } else if(s == 3) {
    # lambda
    l1 <- rnorm(p, 0, sigma_lambda)
    l2 <- rnorm(p, 0, sigma_lambda)
    l3 <- c(rep(1, p/2), rnorm(p/2, 0, sigma0))
    lambda <- cbind(l1, l2, l3)
    # eta
    e1 <- rnorm(n, 0, sigma_eta)
    e2 <- c(rep(1, n/2), rnorm(n/2, 0, sigma0))
    e3 <- rnorm(n, 0, sigma_eta)
    eta <- cbind(e1, e2, e3)
    # meta-covariates (see lambda)
    wB1 <- c(rep(1, p/2), rep(0, p/2))
    wB <- data.frame("wB1" = as.factor(wB1))
  }
  # get output
  sim <- get_output(lambda, eta, wB, hyper, na_prop)
  sim$hyperparameters$seed <- seed
  # save data
  saveRDS(sim, file.path(FOLDER, sim_file))
  # -------------------------------------------------------------------------- #
  # Adaptive Gibbs Sampler without meta-covariates
  out_MCMC <- GGIF::AGS_SIS(Y = sim$Y, y_max = sim$hyperparameters$y_max,
                            X_mean = NULL, X_cov = NULL, W = NULL,
                            seed = 28, stdx = TRUE, stdw = TRUE,
                            kinit = NULL, kmax = NULL, kval = 4,
                            nrun = 20000, thin = 2, start_adapt = 100,
                            b0 = 1, b1 = 5*10^(-4),
                            sd_b = 1, sd_mu = 1, sd_beta = 1,
                            a_theta = 1, b_theta = 1,
                            as = 1, bs = 1,
                            p_constant = 0.5,
                            alpha = alpha, output = c("eta", "lambda"),
                            verbose = FALSE)
  saveRDS(out_MCMC, file.path(FOLDER, "results", paste0("s",s), paste0("res_nometa_", sim_file)))
  # -------------------------------------------------------------------------- #
  if(!is.null(sim$wB)) {
    # Adaptive Gibbs Sampler with meta-covariates
    out_MCMC <- GGIF::AGS_SIS(Y = sim$Y, y_max = sim$hyperparameters$y_max,
                              X_mean = NULL, X_cov = sim$X_cov, W = NULL,
                              seed = 28, stdx = TRUE, stdw = TRUE,
                              kinit = NULL, kmax = NULL, kval = 4,
                              nrun = 20000, thin = 2, start_adapt = 100,
                              b0 = 1, b1 = 5*10^(-4),
                              sd_b = 1, sd_mu = 1, sd_beta = 1,
                              a_theta = 1, b_theta = 1,
                              as = 1, bs = 1,
                              p_constant = 0.5,
                              alpha = alpha, output = c("eta", "lambda"),
                              verbose = FALSE)
    saveRDS(out_MCMC, file.path(FOLDER, "results", paste0("s",s), paste0("res_meta_", sim_file)))
  }
  # -------------------------------------------------------------------------- #
}

parallel::stopCluster(cl = myCluster)

# ---------------------------------------------------------------------------- #
