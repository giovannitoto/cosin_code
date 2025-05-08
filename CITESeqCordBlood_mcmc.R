### ------------------------------------------------------------------------ ###
### CITEseq Cord Blood - MCMC
### ------------------------------------------------------------------------ ###

rm(list=ls())

# ---------------------------------------------------------------------------- #
# DATA
# ---------------------------------------------------------------------------- #

sce <- readRDS("data/CITESeqCordBlood/CITEseqCordBlood_an.RDS")

y <- t(as.matrix(counts(sce)))
wT <- as.data.frame(rowData(sce)[,1:2])
wB <- as.data.frame(rowData(sce)[,-c(1:2)])
x <- as.data.frame(colData(sce)[,3:4])

wTformula <- as.formula("~ length + gc_content")
wBformula <- as.formula(paste("~", paste(colnames(wB),collapse = " + ")))
xFormula <- as.formula(paste("~", paste(colnames(x),collapse = " + ")))

# ---------------------------------------------------------------------------- #
# MCMC
# ---------------------------------------------------------------------------- #

file_name_burn <- "cosin_CITE_burnin"
file_name_mcmc <- "cosin_CITE"

# run MCMC in block of 500 iterations
m_nrun <- 500

# burning MCMC: first 500 iterations
mcmc <- cosin::cosin(y = y, wT = wT, wB = wB, x = x, y_max = Inf,
                     wTformula = wTformula, wBformula = wBformula, xFormula = xFormula,
                     stdwT = TRUE, stdwB = TRUE, stdx = TRUE,
                     sd_gammaT = 1, sd_gammaB = 1, sd_beta = 1/3,
                     a_theta = 1, b_theta = 1, a_sigma = 1, b_sigma = 1,
                     alpha = 10, p_constant = NULL,
                     kinit = NULL, kmax = NULL, kval = 4,
                     nrun = m_nrun, burn = 500, thin = 2,
                     start_adapt = 0, b0 = 1, b1 = 5*10^(-4),
                     seed = 28, output = "all",
                     verbose = TRUE)
# save last iteration
start <- mcmc$last
# delete MCMC
rm(mcmc)

# burnin MCMC: additional 4500 iterations for a total of 5000 burnin iterations
m_list_burn <- 1:9
for(m in m_list_burn) {
  # MCMC
  mcmc <- cosin::cosin(y = y, wT = wT, wB = wB, x = x, y_max = Inf,
                       wTformula = wTformula, wBformula = wBformula, xFormula = xFormula,
                       stdwT = TRUE, stdwB = TRUE, stdx = TRUE,
                       sd_gammaT = 1, sd_gammaB = 1, sd_beta = 1/3,
                       a_theta = 1, b_theta = 1, a_sigma = 1, b_sigma = 1,
                       alpha = 10, p_constant = NULL,
                       kinit = NULL, kmax = NULL, kval = 4,
                       nrun = m_nrun, burn = m_nrun, thin = 2,
                       start_adapt = 0, b0 = 1, b1 = 5*10^(-4),
                       seed = 28, start = start, last = TRUE,
                       output = "all", verbose = TRUE)
  # save MCMC
  saveRDS(mcmc, file.path("results", "CITESeqCordBlood", paste0(file_name_burn,"_", m, ".RDS")))
  # save last iteration
  start <- mcmc$last
  # delete MCMC
  rm(mcmc)
}

# MCMC: 10000 iterations (5000 effective samples considering thin=2)
m_list <- 1:20
for(m in m_list) {
  # MCMC
  mcmc <- cosin::cosin(y = y, wT = wT, wB = wB, x = x, y_max = Inf,
                       wTformula = wTformula, wBformula = wBformula, xFormula = xFormula,
                       stdwT = TRUE, stdwB = TRUE, stdx = TRUE,
                       sd_gammaT = 1, sd_gammaB = 1, sd_beta = 1/3,
                       a_theta = 1, b_theta = 1, a_sigma = 1, b_sigma = 1,
                       alpha = 10, p_constant = NULL,
                       kinit = NULL, kmax = NULL, kval = 4,
                       nrun = m_nrun, burn = 0, thin = 2,
                       start_adapt = 0, b0 = 1, b1 = 5*10^(-4),
                       seed = 28, start = start, last = TRUE,
                       output = "all", verbose = TRUE)
  # save MCMC
  saveRDS(mcmc, file.path("results", "CITESeqCordBlood", paste0(file_name_mcmc,"_", m, ".RDS")))
  # save last iteration
  start <- mcmc$last
  # delete MCMC
  rm(mcmc)
}

# ---------------------------------------------------------------------------- #
# remark : originally each of the following blocks was a different R script
# ---------------------------------------------------------------------------- #
# POSTERIOR ESTIMATES (I)
# ---------------------------------------------------------------------------- #

# useful quantities
n <- nrow(y)
p <- ncol(y)

# parameters to estimate
par2estimate <- c("sigmacol", "omega")

k_max <- 0
lpost_loglik <- c()
lpost_mode <- list("lposterior_max" = -Inf)
lpost_mean <- setNames(as.list(rep(0, length(par2estimate))), par2estimate)
C1_mean <- C2_mean <- list()

draws <- 0
for(m in m_list) {
  # load mcmc
  mcmc <- readRDS(file.path("results", "CITESeqCordBlood", paste0(file_name_mcmc,"_", m, ".RDS")))
  # update number of draws
  draws <- draws + length(mcmc$numFactors)
  # compute loglik and posterior modes
  lpost_mode_m <- cosin::lposterior(out_MCMC = mcmc)
  # save loglik
  lpost_loglik <- c(lpost_loglik, lpost_mode_m$lposterior)
  # save posterior modes if it has the highest loglik
  if(lpost_mode_m$lposterior_max > lpost_mode$lposterior_max) {
    lpost_mode <- lpost_mode_m
  }
  lpost_mean_m <- cosin::posterior_mean(out_MCMC = mcmc, parameters = par2estimate)
  lpost_mean <- as.list(mapply(`+`, lpost_mean_m, lpost_mean))
  # update the maximum number of factors kmax; contributions are ordered by frobenius norm
  contributions_m <- cosin::contributions(out_MCMC = mcmc, reference = NULL)
  k_m <- length(contributions_m$C1)
  if(k_m > kmax) {
    C1_mean <- c(C1_mean, rep(list(matrix(0, n,p)), km-kmax) )
    C2_mean <- c(C2_mean, rep(list(matrix(0, n,p)), km-kmax) )
    kmax <- k_m
  } else if(k_m < kmax) {
    C1_mean <- c(contributions_m$C1, rep(list(matrix(0, n,p)), kmax-k_m))
    C2_mean <- c(contributions_m$C2, rep(list(matrix(0, n,p)), kmax-k_m))
  }
  C1_sum <- mapply("+", C1_mean, contributions_m$C1, SIMPLIFY = FALSE)
  C2_sum <- mapply("+", C2_mean, contributions_m$C2, SIMPLIFY = FALSE)
}
lpost_mode$lposterior <- lpost_loglik
lpost_mean <- lapply(lpost_mean, function(par) par/draws)
C1_mean <- lapply(C1_sum, function(par) par/draws)
C2_mean <- lapply(C2_sum, function(par) par/draws)

# save posterior estimates
saveRDS(lpost_mode, file.path("results", "CITESeqCordBlood", "CITEseqCordBlood_lposterior.RDS"))
saveRDS(lpost_mean, file.path("results", "CITESeqCordBlood", "CITEseqCordBlood_posterior_mean.RDS"))
saveRDS(C1_mean, file.path("results", "CITESeqCordBlood", "CITEseqCordBlood_contributions_mean.RDS"))
saveRDS(C2_mean, file.path("results", "CITESeqCordBlood", "CITEseqCordBlood_square_contributions_mean.RDS"))

# ---------------------------------------------------------------------------- #
# POSTERIOR ESTIMATES (II)
# ---------------------------------------------------------------------------- #

# parameter to save the entire chain. They cannot be too much big.
par2save <- c("beta", "gammaT")

lpost_mean <- setNames(as.list(rep(0, length(par2save))), par2save)
parameter_chains <- setNames(lapply(par2save, function(x) list()), par2save)

for(m in m_list) {
  # load mcmc
  mcmc <- readRDS(file.path("results", "CITESeqCordBlood", paste0(file_name_mcmc,"_", m, ".RDS")))
  # MCMC chains of the parameters of interest
  parameter_chains_m = mcmc[par2save]
  # compute posterior mean
  lpost_mean_m <- cosin::posterior_mean(out_MCMC = mcmc, parameters = par2save)
  lpost_mean <- as.list(mapply(`+`, lpost_mean_m  , lpost_mean))
  for(par in par2save){
    parameter_chains[[par]] <- c(parameter_chains[[par]], parameter_chains_m[[par]])
  }
}
lpost_mean <- lapply(lpost_mean, function(par) par/length(m_list))

# save posterior estimates
saveRDS(parameter_chains, file.path("results", "CITESeqCordBlood", "CITEseqCordBlood_beta_chain.RDS"))
saveRDS(lpost_mean, file.path("results", "CITESeqCordBlood", "CITEseqCordBlood_beta_posterior_mean.RDS"))

# ---------------------------------------------------------------------------- #
# POSTERIOR ESTIMATES (III)
# ---------------------------------------------------------------------------- #

par2estimate <- c("sigmacol", "omega")

library(parallel)
res_mcmc <- mclapply(m_list, FUN = function(m){
  # load mcmc
  mcmc <- readRDS(file.path("results", "CITESeqCordBlood", paste0(file_name_mcmc,"_", m, ".RDS")))
  # compute loglik and posterior modes
  lpost_mode_m <- cosin::lposterior(out_MCMC = mcmc)
  # compute posterior mean
  lpost_mean_m <- cosin::posterior_mean(out_MCMC = mcmc, parameters = par2estimate)
  return(list(lpost_mode = lpost_mode_m, lpost_mean = lpost_mean_m))
}, mc.cores = length(m_list))

chain_max <- which.max(lapply(res_mcmc, function(x) x$lpost_mode$lposterior_max))
res_mcmc[[chain_max]]$lpost_mode$iteration_max  <- res_mcmc[[chain_max]]$lpost_mode$iteration_max + 250*(chain_max-1)
lpost_mode <- res_mcmc[[chain_max]]$lpost_mode

# save lposterior mode
saveRDS(lpost_mode, file.path("results", "CITESeqCordBlood", "CITEseqCordBlood_lposterior.RDS"))

# compute lposterior mean
lpost_mean <- setNames(as.list(rep(0, length(par2estimate))), par2estimate)
for(m in m_list) {
  lpost_mean <- as.list(mapply(`+`, res_mcmc[[m]]$lpost_mean , lpost_mean))
}
lpost_mean <- lapply(lpost_mean, function(par) par/length(m_list))
# save lposterior mean
saveRDS(lpost_mean, file.path("results", "CITESeqCordBlood", "CITEseqCordBlood_posterior_mean.RDS"))

# ---------------------------------------------------------------------------- #
# POSTERIOR ESTIMATES (IV)
# ---------------------------------------------------------------------------- #

par2estimate <- c("omega_pcorr")

library(parallel)
res_mcmc <- mclapply(m_list, FUN = function(m){
  # load mcmc
  mcmc <- readRDS(file.path("results", "CITESeqCordBlood", paste0(file_name_mcmc,"_", m, ".RDS")))
  # compute posterior mean
  lpost_mean_m <- cosin::posterior_mean(out_MCMC = mcmc, parameters = par2estimate)

  return(list(lpost_mean = lpost_mean_m))
}, mc.cores = length(m_list))

# compute lposterior mean
lpost_mean <- setNames(as.list(rep(0, length(par2estimate))), par2estimate)
for(m in m_list){
  lpost_mean$omega_pcorr <- res_mcmc[[m]]$lpost_mean$omega_pcorr + lpost_mean$omega_pcorr
}
lpost_mean <- lapply(lpost_mean, function(par) par/length(m_list))

# save lposterior mean
saveRDS(lpost_mean, file.path("results", "CITESeqCordBlood", "CITEseqCordBlood_pcorr_mean.RDS"))

# ---------------------------------------------------------------------------- #
