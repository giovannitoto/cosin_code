# ---------------------------------------------------------------------------- #

aglmpca <- function(Y, L, meta = NULL, pca = "ppca", completeY = TRUE) {
  # -------------------------------------------------------------------------- #
  # APPROXIMATION OF GLM-PCA / POISSON FAMILY + PEARSON RESIDUALS + PPCA
  # -------------------------------------------------------------------------- #
  # We follow Townes et al. (2019) applying principal component analysis on
  # Pearson residuals of the null Poisson model as fast approximation of GLM-PCA,
  # naming such approach aGLM-PCA.
  # -------------------------------------------------------------------------- #
  #         Y = pxn matrix with genes in the rows and samples in the columns
  #         L = number of PCs (integer)
  #      meta = pxq matrix of metacovariates
  #       pca = c("nipals", "ppca", "svdImpute")
  # completeY = c(TRUE, FALSE)
  # -------------------------------------------------------------------------- #
  #    nipals : tolerant to amounts of missing values not more than 5%
  #      ppca : tolerant to amounts of missing values between 10% to 15%
  # svdImpute : tolerant to amounts of missing values greater than 10%
  # -------------------------------------------------------------------------- #
  # adaptation of the code available here:
  # https://code.bioconductor.org/browse/scry/blob/RELEASE_3_16/R/nullResiduals.R
  # -------------------------------------------------------------------------- #
  if(is.null(meta)) {
    sz <- colSums(Y, na.rm = TRUE)
    lsz <- log(sz)
    sz <- exp(lsz-mean(lsz))
    lambda <- rowSums(Y, na.rm = TRUE) / sum(sz)
    mhat <- outer(lambda, sz)
    res <- (Y - mhat) / sqrt(mhat)
    res[is.na(res)] <- 0  # case of 0/0
    res[is.na(Y)] <- NA   # keep NAs of original data
  } else {
    stop("[meta] Not implemented yet!")
  }
  # res must be a pxn matrix (same dim as Y)
  pc <- pcaMethods::pca(res, nPcs = L, method = pca)
  if(completeY) {
    mres <- pc@completeObs
  } else {
    mres <- pc@scores %*% t(pc@loadings)
  }
  Y_woNA <- mres * sqrt(mhat) + mhat
  Y_woNA[Y_woNA < 0] <- 0
  return(Y_woNA)
}

# ---------------------------------------------------------------------------- #
