# ---------------------------------------------------------------------------- #
#  COSIN : c(n, p, sigma, s, r, metrics_cont, mad_Y, mae_Y, rmse_etaLambda, kstar, sd_cont, sd_etaLambda, meta/nometa)
# ---------------------------------------------------------------------------- #
# GLMPCA : (n, p, sigma, r, best_alpha, best_mae, mae_Y, rmse_cont, rmse_etaLambda, meta)
# ---------------------------------------------------------------------------- #

rm(list=ls())
setwd("simulations/metrics")
library(dplyr)

# ---------------------------------------------------------------------------- #
# SCENARIO 1
# ---------------------------------------------------------------------------- #

metrics_s1 <- rbind(readRDS("table_metrics_meta_NA_p100_s1.RDS"),
                    readRDS("table_metrics_nometa_NA_p100_s1.RDS"),
                    readRDS("table_metrics_meta_NA_p1000_s1.RDS"),
                    readRDS("table_metrics_nometa_NA_p1000_s1.RDS"))

colnames(metrics_s1) <- c("n", "p", "sigma", "s", "r", "RMSE_1", "RMSE_2", "RMSE_tot",
                          "F1_1", "F1_2", "F1_tot", "MAD_Y", "MAE_Y", "RMSE_el",
                          "kstar", "sd_1", "sd_2", "sd_tot", "meta")

metrics_s1 <- as.data.frame(metrics_s1)
metrics_s1$s <- NULL
metrics_s1$r <- NULL

metrics_s1 <- metrics_s1 %>% group_by(n, p, meta, sigma) %>%
  summarise_all(.funs = c("median" = median, "IQR" = IQR))
# order columns
metrics_s1 <- metrics_s1[, c(1:4, rbind(5:17, 5:17+13))]
# change colnames
colnames(metrics_s1) <- colnames(metrics_s1) %>%
  stringr::str_replace_all("_median", " median") %>%
  stringr::str_replace("_IQR", " IQR")
# table latex
# print(xtable::xtable(metrics_s1, digits = 2), include.rownames = FALSE)

View(round(metrics_s1[,c(1:4,19:20)], 4))  # MAE

# without NA
metrics_s1_woNA <- rbind(readRDS("table_metrics_meta_woNA_p100_s1.RDS"),
                         readRDS("table_metrics_nometa_woNA_p100_s1.RDS"),
                         readRDS("table_metrics_meta_woNA_p1000_s1.RDS"),
                         readRDS("table_metrics_nometa_woNA_p1000_s1.RDS"))

colnames(metrics_s1_woNA) <- c("n", "p", "sigma", "s", "r", "RMSE_1", "RMSE_2", "RMSE_tot",
                               "F1_1", "F1_2", "F1_tot", "MAD_Y", "MAE_Y", "RMSE_el",
                               "kstar", "sd_1", "sd_2", "sd_tot", "meta")

metrics_s1_woNA <- as.data.frame(metrics_s1_woNA)
metrics_s1_woNA$s <- NULL
metrics_s1_woNA$r <- NULL

# filter columns
metrics_s1_woNA <- metrics_s1_woNA[, c(1:3,4:5,17)]

metrics_s1_woNA <- metrics_s1_woNA %>% group_by(n, p, meta, sigma) %>%
  summarise_all(.funs = c("median" = median, "IQR" = IQR))

View(round(metrics_s1_woNA, 2))  # RMSE
round(as.matrix(metrics_s1_woNA), 2)

# ---------------------------------------------------------------------------- #
# SCENARIO 2
# ---------------------------------------------------------------------------- #

metrics_s2 <- rbind(readRDS("table_metrics_nometa_NA_p100_s2.RDS"),
                    readRDS("table_metrics_nometa_NA_p1000_s2.RDS"))

colnames(metrics_s2) <- c("n", "p", "sigma", "s", "r", "RMSE_1", "RMSE_2", "RMSE_tot",
                          "MAD_Y", "MAE_Y", "RMSE_el", "kstar", "sd_1", "sd_2", "sd_tot", "meta")

metrics_s2 <- as.data.frame(metrics_s2)
metrics_s2$s <- NULL
metrics_s2$r <- NULL

metrics_s2 <- metrics_s2 %>% group_by(n, p, meta, sigma) %>%
  summarise_all(.funs = c("median" = median, "IQR" = IQR))
# order columns
metrics_s2 <- metrics_s2[, c(1:4, rbind(5:14, 5:14+10))]
# change colnames
colnames(metrics_s2) <- colnames(metrics_s2) %>%
  stringr::str_replace_all("_median", " median") %>%
  stringr::str_replace("_IQR", " IQR")
# table latex
# print(xtable::xtable(metrics_s2, digits = 2), include.rownames = FALSE)

View(round(metrics_s2[,c(1:4,13:14)], 4))  # MAE

# without NA
metrics_s2_woNA <- rbind(readRDS("table_metrics_nometa_woNA_p100_s2.RDS"),
                         readRDS("table_metrics_nometa_woNA_p1000_s2.RDS"))

colnames(metrics_s2_woNA) <- c("n", "p", "sigma", "s", "r", "RMSE_1", "RMSE_2", "RMSE_tot",
                          "MAD_Y", "MAE_Y", "RMSE_el", "kstar", "sd_1", "sd_2", "sd_tot", "meta")

metrics_s2_woNA <- as.data.frame(metrics_s2_woNA)
metrics_s2_woNA$s <- NULL
metrics_s2_woNA$r <- NULL

# filter columns
metrics_s2_woNA <- metrics_s2_woNA[, c(1:3,4:5,14)]

metrics_s2_woNA <- metrics_s2_woNA %>% group_by(n, p, meta, sigma) %>%
  summarise_all(.funs = c("median" = median, "IQR" = IQR))

View(round(metrics_s2_woNA, 2))  # RMSE
round(as.matrix(metrics_s2_woNA), 2)

# ---------------------------------------------------------------------------- #
# SCENARIO 3
# ---------------------------------------------------------------------------- #

metrics_s3 <- rbind(readRDS("table_metrics_meta_NA_p100_s3.RDS"),
                    readRDS("table_metrics_nometa_NA_p100_s3.RDS"),
                    readRDS("table_metrics_meta_NA_p1000_s3.RDS"),
                    readRDS("table_metrics_nometa_NA_p1000_s3.RDS"))

colnames(metrics_s3) <- c("n", "p", "sigma", "s", "r", "RMSE_1", "RMSE_2",
                          "RMSE_3", "RMSE_tot", "MAD_Y", "MAE_Y", "RMSE_el",
                          "kstar", "sd_1", "sd_2", "sd_3", "sd_tot", "meta")

metrics_s3 <- as.data.frame(metrics_s3)
metrics_s3$s <- NULL
metrics_s3$r <- NULL

metrics_s3 <- metrics_s3 %>% group_by(n, p, meta, sigma) %>%
  summarise_all(.funs = c("median" = median, "IQR" = IQR))
# order columns
metrics_s3 <- metrics_s3[, c(1:4, rbind(5:16, 5:16+12))]
# change colnames
colnames(metrics_s3) <- colnames(metrics_s3) %>%
  stringr::str_replace_all("_median", " median") %>%
  stringr::str_replace("_IQR", " IQR")
# table latex
# print(xtable::xtable(metrics_s3, digits = 2), include.rownames = FALSE)

View(round(metrics_s3[,c(1:4,15:16)], 4))  # MAE

# without NA
metrics_s3_woNA <- rbind(readRDS("table_metrics_meta_woNA_p100_s3.RDS"),
                         readRDS("table_metrics_nometa_woNA_p100_s3.RDS"),
                         readRDS("table_metrics_meta_woNA_p1000_s3.RDS"),
                         readRDS("table_metrics_nometa_woNA_p1000_s3.RDS"))

colnames(metrics_s3_woNA) <- c("n", "p", "sigma", "s", "r", "RMSE_1", "RMSE_2",
                          "RMSE_3", "RMSE_tot", "MAD_Y", "MAE_Y", "RMSE_el",
                          "kstar", "sd_1", "sd_2", "sd_3", "sd_tot", "meta")

metrics_s3_woNA <- as.data.frame(metrics_s3_woNA)
metrics_s3_woNA$s <- NULL
metrics_s3_woNA$r <- NULL

# filter columns
metrics_s3_woNA <- metrics_s3_woNA[, c(1:3,4:6,16)]

metrics_s3_woNA <- metrics_s3_woNA %>% group_by(n, p, meta, sigma) %>%
  summarise_all(.funs = c("median" = median, "IQR" = IQR))

View(round(metrics_s3_woNA, 2))  # RMSE
round(as.matrix(metrics_s3_woNA), 2)


# ---------------------------------------------------------------------------- #
# SCENARIO 1
# ---------------------------------------------------------------------------- #

# without metacovariates
metrics_glmpca_s1 <- readRDS("table_metrics_nometa_glmpca_s1.RDS")
colnames(metrics_glmpca_s1) <- c("n", "p", "sigma", "r", "alpha_a", "MAE_a",
                                 "MAE", "RMSE_1", "RMSE_2", "RMSE_tot",
                                 "rmse_etaLambda", "meta")
metrics_glmpca_s1 <- as.data.frame(metrics_glmpca_s1)
metrics_glmpca_s1$r <- NULL
metrics_glmpca_s1 <- metrics_glmpca_s1 %>% group_by(n, p, meta, sigma) %>%
  summarise_all(.funs = c("median" = median, "IQR" = IQR))
View(round(metrics_glmpca_s1, 2))  # all
round(as.matrix(metrics_glmpca_s1[,c(1:4,6,13)]), 4)       # MAE
round(as.matrix(metrics_glmpca_s1[,c(1:4,8:9,15:16)]), 2)  # RMSE

# with metacovariates
metrics_glmpca_s1 <- readRDS("table_metrics_meta_glmpca_s1.RDS")
colnames(metrics_glmpca_s1) <- c("n", "p", "sigma", "r", "alpha_a", "MAE_a",
                                 "MAE", "RMSE_1", "RMSE_2", "RMSE_tot",
                                 "rmse_etaLambda", "meta")
metrics_glmpca_s1 <- as.data.frame(metrics_glmpca_s1)
metrics_glmpca_s1$r <- NULL
metrics_glmpca_s1 <- metrics_glmpca_s1 %>% group_by(n, p, meta, sigma) %>%
  summarise_all(.funs = c("median" = median, "IQR" = IQR))
View(round(metrics_glmpca_s1, 2))  # all
round(as.matrix(metrics_glmpca_s1[,c(1:4,8:9,15:16)]), 2)  # RMSE

# ---------------------------------------------------------------------------- #
# SCENARIO 2
# ---------------------------------------------------------------------------- #

# without metacovariates
metrics_glmpca_s2 <- readRDS("table_metrics_nometa_glmpca_s2.RDS")
colnames(metrics_glmpca_s2) <- c("n", "p", "sigma", "r", "alpha_a", "MAE_a",
                                 "MAE", "RMSE_1", "RMSE_2", "RMSE_tot",
                                 "rmse_etaLambda", "meta")
metrics_glmpca_s2 <- as.data.frame(metrics_glmpca_s2)
metrics_glmpca_s2$r <- NULL
metrics_glmpca_s2 <- metrics_glmpca_s2 %>% group_by(n, p, meta, sigma) %>%
  summarise_all(.funs = c("median" = median, "IQR" = IQR))
View(round(metrics_glmpca_s2, 2))  # all
round(as.matrix(metrics_glmpca_s2[,c(1:4,6,13)]), 4)       # MAE
round(as.matrix(metrics_glmpca_s2[,c(1:4,8:9,15:16)]), 2)  # RMSE

# ---------------------------------------------------------------------------- #
# SCENARIO 3
# ---------------------------------------------------------------------------- #

# without metacovariates
metrics_glmpca_s3 <- readRDS("table_metrics_nometa_glmpca_s3.RDS")
colnames(metrics_glmpca_s3) <- c("n", "p", "sigma", "r", "alpha_a", "MAE_a",
                                 "MAE", "RMSE_1", "RMSE_2", "RMSE_3", "RMSE_tot",
                                 "rmse_etaLambda", "meta")
metrics_glmpca_s3 <- as.data.frame(metrics_glmpca_s3)
metrics_glmpca_s3$r <- NULL
metrics_glmpca_s3 <- metrics_glmpca_s3 %>% group_by(n, p, meta, sigma) %>%
  summarise_all(.funs = c("median" = median, "IQR" = IQR))
View(round(metrics_glmpca_s3, 2))  # all
round(as.matrix(metrics_glmpca_s3[,c(1:4,6,14)]), 4)        # MAE
round(as.matrix(metrics_glmpca_s3[,c(1:4,8:10,16:18)]), 2)  # RMSE

# with metacovariates
metrics_glmpca_s3 <- readRDS("table_metrics_meta_glmpca_s3.RDS")
colnames(metrics_glmpca_s3) <- c("n", "p", "sigma", "r", "alpha_a", "MAE_a",
                                 "MAE", "RMSE_1", "RMSE_2", "RMSE_3", "RMSE_tot",
                                 "rmse_etaLambda", "meta")
metrics_glmpca_s3 <- as.data.frame(metrics_glmpca_s3)
metrics_glmpca_s3$r <- NULL
metrics_glmpca_s3 <- metrics_glmpca_s3 %>% group_by(n, p, meta, sigma) %>%
  summarise_all(.funs = c("median" = median, "IQR" = IQR))
View(round(metrics_glmpca_s3, 2))  # all
round(as.matrix(metrics_glmpca_s3[,c(1:4,8:10,16:18)]), 2)  # RMSE

# ---------------------------------------------------------------------------- #
