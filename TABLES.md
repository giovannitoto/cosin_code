# Tables

Tables generated by `simulations_NA_metrics.R` and `simulations_woNA_metrics.R` have the following columns:

| column                  | description                            
|-------------------------|----------------------------------------
| `n`                     | Number of cells
| `p`                     | Number of genes
| `sigma`                 | Standard deviation of errors
| `s`                     | Scenario (1,2,3)
| `r`                     | Replication (1,...,50)
| `RMSE_1`, `RMSE_2`, ... | MCMC average of Root Mean Square Error computed on each contribution
| `RMSE_tot`              | Sum of `RMSE_1`, `RMSE_2`, ...
| `F1_1`, `F1_2`, ...     | MCMC average of F1 computed on each contribution
| `F1_tot`                | Sum of `F1_1`, `F1_2`, ...
| `MAD_Y`                 | MCMC average of Mean Absolute Deviation computed on the gene expression matrix
| `MAE_Y`                 | MCMC average of Mean Absolute Error computed on the gene expression matrix
| `RMSE_el`               | MCMC average of Root Mean Square Error computed on the posterior mean of `eta%*%t(lambda)`
| `kstar`                 | Average number of active latent factors
| `sd1`, `sd2`, ...       | Standard deviation of each real contribution
| `sd_tot`                | Sum of `sd1`, `sd2`, ...
| `meta`                  | Whether the model was estimated using meta-covariates (0,1)

**remark:** `MAD_Y` and `MAE_Y` are computed out-of-sample in the tables generated by `simulations_NA_metrics.R`.


Tables generated by `simulations_glmpca.R` have the following columns:

| column                  | description                            |
|-------------------------|----------------------------------------|
| `n`                     | Number of cells
| `p`                     | Number of genes
| `sigma`                 | Standard deviation of errors
| `s`                     | Scenario (1,2,3)
| `r`                     | Replication (1,...,50)
| `alpha_a`               | Number of latent factors for which GLM-PCA approximation reaches the lowest out-of-sample MAE (this value is `MAE_a`)
| `MAE_a`                 | Out-of-sample Mean Absolute Error computed on the gene expression matrix with GLM-PCA approximation
| `MAE`                   | Mean Absolute Error computed on the gene expression matrix with GLM-PCA
| `RMSE_1`, `RMSE_2`, ... | MCMC average of the Root Mean Square Error computed on each contribution with GLM-PCA
| `RMSE_tot`              | Sum of `RMSE_1`, `RMSE_2`, ...
| `rmse_etaLambda`        | Root Mean Square Error computed on the posterior mean of `eta%*%t(lambda)` with GLM-PCA
| `meta`                  | Whether the model was estimated using meta-covariates (0,1)

**remark:** `alpha_a` and `MAE_a` are obtained from GLM-PCA approximation applied to synthetic data sets with NAs; the remaining metrics from GLM-PCA applied to synthetic data sets without NAs.