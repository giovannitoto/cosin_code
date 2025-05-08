# ---------------------------------------------------------------------------- #
# replications: R = 50
#   n     p  |            scenario
#  50   100  |  full + col1
# 200   100  |  full + row1
# 200  1000  |  full + row1 + col1
# ---------------------------------------------------------------------------- #

# R --slave --no-restore --file=simulations_clustering.R &

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
                 "y" = y, "wB" = wB, "hyperparameters" = hyper)
  return(output)
}

get_ARI <- function(groups, scores, Y, J, method = "kmeans") {
  # method = c("kmeans","seurat")
  if (method == "kmeans") {
    ARI <- pdfCluster::adj.rand.index(groups, kmeans(scale(scores), centers = J, nstart = 10)$cluster)
  } else {
    Sco <- Seurat::CreateSeuratObject(counts = Y)
    rownames(scores) <- paste0("Cell",1:nrow(scores))
    Sco[["model"]] <- Seurat::CreateDimReducObject(embeddings = scores, key = "model_")
    Sco <- Seurat::FindNeighbors(Sco, reduction = "model", dims = 1:ncol(scores))
    Sco <- Seurat::FindClusters(Sco)
    ARI <- pdfCluster::adj.rand.index(groups, Sco$seurat_clusters)
  }
  return(ARI)
}

# ---------------------------------------------------------------------------- #

# max number of cores
maxCores <- 25

# ---------------------------------------------------------------------------- #

y_max <- Inf
sigma_lambda <- 1
sigma_eta <- 1
sigma0 <- 0.05

clustMethod <- "kmeans"  # c("kmeans","seurat")

table_filename <- paste0("table_metric_clustering_", clustMethod)

# ---------------------------------------------------------------------------- #

hyper <- list("y_max" = y_max, "sigma_lambda" = sigma_lambda,
              "sigma_eta" = sigma_eta, "sigma0" = sigma0)

np_list <-  list(c(50, 100), c(200, 100), c(200, 1000))
sigma_list <- c(1/10, 1)
s_list <- c(1, 2, 3)
J_list <- c(5, 10)
R_list <- 1:50

# grid of all possible combinations of n, p, sigma, scenario
sim_list <- expand.grid(r = R_list, sigma = sigma_list, J = J_list, np = np_list, s = s_list)
# seeds of the simulations
sim_list$seed <- 1:nrow(sim_list)

# filter combinations
sim_list <- sim_list[sim_list$r <= 30,]
sim_list <- sim_list[sim_list$s == 3,]

# ---------------------------------------------------------------------------- #

library(foreach)
numCores <- min(parallel::detectCores() - 1, maxCores)
myCluster <- parallel::makeCluster(numCores,type = "PSOCK")  # "PSOCK", "FORK"
doParallel::registerDoParallel(cl = myCluster)

table_metrics <- foreach(i = 1:nrow(sim_list), .combine = 'rbind', .inorder = FALSE) %dopar% {
  # -------------------------------------------------------------------------- #
  # get n, p, sigma, scenario (s), replication (r), seed
  n <- sim_list$np[i][[1]][1]
  p <- sim_list$np[i][[1]][2]
  sigma <- sim_list$sigma[i]
  J <- sim_list$J[i]
  s <- sim_list$s[i]
  r <- sim_list$r[i]
  seed <- sim_list$seed[i]
  # -------------------------------------------------------------------------- #
  # hyperparameters
  hyper$sigma <- sigma
  # set seed
  set.seed(seed)
  # -------------------------------------------------------------------------- #
  # true cell groups
  groups <- rep(1:J, each = n/J)
  # eta and lambda
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
    eta_means <- rnorm(J, 0, sigma_eta*3)
    e1 <- rnorm(n, rep(eta_means, each = n/J), sigma_eta)
    e2 <- c(rep(1, n/J), rnorm((J-1)*n/J, 0, sigma0))
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
    eta_means <- rnorm(J*2, 0, sigma_eta*3)
    e1 <- rnorm(n, rep(eta_means[1:J], each = n/J), sigma_eta)
    e2 <- c(rep(1, n/J), rnorm((J-1)*n/J, 0, sigma0))
    e3 <- rnorm(n, rep(eta_means[-c(1:J)], each = n/J), sigma_eta)
    eta <- cbind(e1, e2, e3)
    # meta-covariates (see lambda)
    wB1 <- c(rep(1, p/2), rep(0, p/2))
    wB <- data.frame("wB1" = as.factor(wB1))
  }
  # get output
  sim <- get_output(lambda, eta, wB, hyper, na_prop)
  sim$hyperparameters$seed <- seed
  # transpose
  sim$y <- t(sim$y)
  # remove unexpressed genes
  sim$y <- sim$y[rowSums(sim$y) > 0, ]
  # assign rownames and colnames
  rownames(sim$y) <- paste0("Gene", 1:nrow(sim$y))
  colnames(sim$y) <- paste0("Cell", 1:ncol(sim$y))
  # -------------------------------------------------------------------------- #
  # init ARI
  ARI <- list()
  # -------------------------------------------------------------------------- #
  # PCA on log-counts
  # -------------------------------------------------------------------------- #
  Y_log <- log(t(sim$y) / colSums(sim$y) * 10000 + 1)  # n x p
  # PCA
  res <- prcomp(Y_log, scale. = TRUE)
  # proportion of variance explained by each principal component
  var_explained <- res$sdev^2 / sum(res$sdev^2)
  # number of factors explaining at least 90% variance
  K <- which(cumsum(var_explained) >= 0.90)[1]
  # ARI
  ARI[["pca_K"]] <- get_ARI(groups = groups, scores = res$x[,1:K], Y = sim$y, J = J, method = clustMethod)
  ARI[["pca"]] <- get_ARI(groups = groups, scores = res$x[,1:sim$k], Y = sim$y, J = J, method = clustMethod)
  # -------------------------------------------------------------------------- #
  # COSIN without meta-covariates
  # -------------------------------------------------------------------------- #
  res <- cosin::cosin(y = t(sim$y), wT = NULL, wB = NULL, x = NULL,
                      y_max = sim$hyperparameters$y_max,
                      stdwT = TRUE, stdwB = TRUE, stdx = TRUE,
                      sd_gammaT = 1, sd_gammaB = 1, sd_beta = 1,
                      a_theta = 1, b_theta = 1, a_sigma = 1, b_sigma = 1,
                      alpha = K, p_constant = 0.5,
                      kinit = NULL, kmax = NULL, kval = 4,
                      nrun = 20000, thin = 2,
                      start_adapt = 100, b0 = 1, b1 = 5*10^(-4),
                      seed = 28, output = "all",
                      verbose = TRUE)
  # scores
  scores <- cosin::lposterior(res, parameters = "eta")$eta_max
  # ARI
  ARI[["cosin"]] <- get_ARI(groups = groups, scores = scores, Y = sim$y, J = J, method = clustMethod)
  # -------------------------------------------------------------------------- #
  # GLM-PCA
  # -------------------------------------------------------------------------- #
  # GLM-PCA with K factors
  tryCatch({
    res <- glmpca::glmpca(sim$y, L = K, fam = "poi", ctl=list("maxIter" = 3000))
    # compute ARI
    ARI[["glmpca_K"]] <- get_ARI(groups = groups, scores = as.matrix(res$factors), Y = sim$y, J = J, method = clustMethod) 
  }, error = function(e) { ARI[["glmpca_K"]] <- NA } )
  if(!("glmpca_K" %in% names(ARI))) { ARI[["glmpca_K"]] <- NA }
  # GLM-PCA with sim$k factors
  tryCatch({
    res <- glmpca::glmpca(sim$y, L = sim$k, fam = "poi", ctl=list("maxIter" = 3000))
    # compute ARI
    ARI[["glmpca"]] <- get_ARI(groups = groups, scores = as.matrix(res$factors), Y = sim$y, J = J, method = clustMethod)
  }, error = function(e) { ARI[["glmpca"]] <- NA })
  if(!("glmpca" %in% names(ARI))) { ARI[["glmpca"]] <- NA }
  # -------------------------------------------------------------------------- #
  # fastglmpca
  # -------------------------------------------------------------------------- #
  # fastglmpca with K factors
  tryCatch({
    res <- fastglmpca::init_glmpca_pois(sim$y, K = K)
    res <- fastglmpca::fit_glmpca_pois(sim$y, fit0 = res, control=list("maxiter" = 3000))
    # compute ARI
    ARI[["fastglmpca_K"]] <- get_ARI(groups = groups, scores = res$V, Y = sim$y, J = J, method = clustMethod)
  }, error = function(e) { ARI[["fastglmpca_K"]] <- NA } )
  if(!("fastglmpca_K" %in% names(ARI))) { ARI[["fastglmpca_K"]] <- NA }
  # fastglmpca with sim$k factors
  tryCatch({
    res <- fastglmpca::init_glmpca_pois(sim$y, K = sim$k)
    res <- fastglmpca::fit_glmpca_pois(sim$y, fit0 = res, control=list("maxiter" = 3000))
    # compute ARI
    ARI[["fastglmpca"]] <- get_ARI(groups = groups, scores = res$V, Y = sim$y, J = J, method = clustMethod)
  }, error = function(e) { ARI[["fastglmpca"]] <- NA } )
  if(!("fastglmpca" %in% names(ARI))) { ARI[["fastglmpca"]] <- NA }
  # -------------------------------------------------------------------------- #
  # scGBM
  # -------------------------------------------------------------------------- #
  # scGBM with K factors
  tryCatch({
    res <- scGBM::gbm.sc(sim$y, M = K, infer.beta = TRUE)
    # compute ARI
    ARI[["scGBM_K"]] <- get_ARI(groups = groups, scores = res$scores, Y = sim$y, J = J, method = clustMethod)
  }, error = function(e) { ARI[["scGBM_K"]] <- NA } )
  if(!("scGBM_K" %in% names(ARI))) { ARI[["scGBM_K"]] <- NA }
  # scGBM with sim$k factors
  tryCatch({
    res <- scGBM::gbm.sc(sim$y, M = sim$k, infer.beta = TRUE)
    # compute ARI
    ARI[["scGBM"]] <- get_ARI(groups = groups, scores = res$scores, Y = sim$y, J = J, method = clustMethod)
  }, error = function(e) { ARI[["scGBM"]] <- NA } )
  if(!("scGBM" %in% names(ARI))) { ARI[["scGBM"]] <- NA }
  # -------------------------------------------------------------------------- #
  # NewWave
  # -------------------------------------------------------------------------- #
  se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = sim$y))
  # NewWave with K factors
  tryCatch({
    res <- NewWave::newWave(se, K = K, verbose = FALSE)
    # compute ARI
    ARI[["newWave_K"]] <- get_ARI(groups = groups, scores = SingleCellExperiment::reducedDims(res)$newWave, Y = sim$y, J = J, method = clustMethod)
  }, error = function(e) { ARI[["newWave_K"]] <- NA } )
  if(!("newWave_K" %in% names(ARI))) { ARI[["newWave_K"]] <- NA }
  # NewWave with sim$k factors
  tryCatch({
    res <- NewWave::newWave(se, K = sim$k, verbose = FALSE)
    # compute ARI
    ARI[["newWave"]] <- get_ARI(groups = groups, scores = SingleCellExperiment::reducedDims(res)$newWave, Y = sim$y, J = J, method = clustMethod)
  }, error = function(e) { ARI[["newWave"]] <- NA } )
  if(!("newWave" %in% names(ARI))) { ARI[["newWave"]] <- NA }
  # -------------------------------------------------------------------------- #
  # create row
  met <- c(n, p, sigma, J, s, r, sim$k, K, ARI[["cosin"]],
           ARI[["pca_K"]], ARI[["glmpca_K"]], ARI[["fastglmpca_K"]], ARI[["scGBM_K"]], ARI[["newWave_K"]],
           ARI[["pca"]], ARI[["glmpca"]], ARI[["fastglmpca"]], ARI[["scGBM"]], ARI[["newWave"]])
  # update txt
  cat(met, "\n", sep=",", file=paste0(table_filename,".txt"), append = TRUE)
  # update table
  met
  # -------------------------------------------------------------------------- #
}

# ---------------------------------------------------------------------------- #

# save metrics
saveRDS(table_metrics, file.path("results", "clustering", paste0(table_filename,".RDS")))

# ---------------------------------------------------------------------------- #

parallel::stopCluster(cl = myCluster)

# ---------------------------------------------------------------------------- #
