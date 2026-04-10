library(dplyr)
library(readr)
library(tidyr)
library(maps)
library(ggplot2)
library(gridExtra)
library(writexl)

data_dir <- "data/Hengl-openlandmap-trainingdata"
sol_chem_pnts_horizons <- readRDS(file.path(data_dir, "sol_chem.pnts_horizons.rds")) %>% as_tibble()

# write locations into a file for AEE fetch
cols_select <- c("site_key", "longitude_decimal_degrees", "latitude_decimal_degrees",
"hzn_top", "hzn_bot", "hzn_depth",
"clay_total", "silt_total", "sand_total", "organic_carbon", "oc_d",    
"total_carbon_ncs", "total_nitrogen_ncs", "ph_kcl", "ph_h2o",
"ph_cacl2", "sum_of_cations_cec_pH_8_2", "cec_nh4_ph_7",
"ecec_base_plus_aluminum", "total_frag_wt_pct_gt_2_mm_ws", "bulk_density_oven_dry",
"ca_ext", "mg_ext", "na_ext", "k_ext",
"ec_water_extractable", "ec_predict_one_to_two")

hzn_depth_threshold = 20
sol_chem_pnts_horizons_select = sol_chem_pnts_horizons %>%
    dplyr::select(cols_select) %>%
    dplyr::filter(hzn_depth <= hzn_depth_threshold)


sol_chem_pnts_horizons_uniqlocs <- sol_chem_pnts_horizons %>%
    dplyr::mutate(id = row_number()) %>%
    dplyr::select(id, cols_select) %>%
    dplyr::mutate(id = as.character(id),
                  longitude = as.character(longitude_decimal_degrees),
                  latitude = as.character(latitude_decimal_degrees)) %>%
    dplyr::select(id, latitude, longitude) %>%
    dplyr::distinct(latitude, longitude, .keep_all = TRUE)

write.csv(sol_chem_pnts_horizons_uniqlocs, "data/sol_chem_pnts_horizons_uniqlocs.csv", quote = TRUE, row.names = FALSE)

n <- nrow(sol_chem_pnts_horizons_uniqlocs)
split_idx <- cut(seq_len(n), breaks = 3, labels = FALSE)
for (i in 1:3) {
    write.csv(sol_chem_pnts_horizons_uniqlocs[split_idx == i, ],
              paste0("data/sol_chem_pnts_horizons_uniqlocs_part", i, ".csv"),
              quote = TRUE, row.names = FALSE)
}

#writexl::write_xlsx(head(sol_chem_pnts_horizons, 5000), "data/sol_chem_pnts_horizons_5000.xlsx")

#p_hzn_depth <- ggplot(sol_chem_pnts_horizons, aes(x = hzn_depth)) +
#    geom_histogram(bins = 50, fill = "steelblue", color = "white") +
#    theme_minimal() +
#    labs(title = "Distribution of hzn_depth", x = "hzn_depth", y = "Count")
#ggsave("data/hzn_depth_distribution.pdf", plot = p_hzn_depth, width = 8, height = 5)

library(glmnet)
library(caret)

cols2predict = colnames(sample2alldata1)[c(66:253, 255:302)]
X = as.matrix(sample2alldata1[, AE_cols])
rownames(X) = sample2alldata1 %>% pull("SampleSiteCode")
#X <- scale(X)

lasso.results = matrix(nrow = length(cols2predict), ncol = 1)
rownames(lasso.results) = cols2predict
colnames(lasso.results) = c("r2")
for(c in 1:length(cols2predict)) {
  y = sample2alldata1 %>% pull(cols2predict[c])
  names(y) = sample2alldata1 %>% pull(SampleSiteCode)
  
  if(sum(y, na.rm = T) == 0) {
    next
  }
  sites2include = names(which(is.na(y) == F))
  X.sites2include = X[sites2include, ]
  y.sites2include = y[sites2include]
  # 3. Use cross-validation to find the optimal lambda (best penalty parameter)
  # cv.glmnet automatically fits the model across a range of lambda values.
  # alpha = 1 specifies Lasso regression (alpha = 0 is Ridge, between 0 and 1 is Elastic Net).
  set.seed(123) # for reproducibility
  cv_model <- cv.glmnet(X.sites2include, y.sites2include, alpha = 1, nfolds = 10) # 10-fold cross-validation is common
  best_lambda <- cv_model$lambda.min
  best_model <- glmnet(X.sites2include, y.sites2include, alpha = 1, lambda = best_lambda)
  # Plot the cross-validation results to visualize MSE vs. log(lambda)
  predicted = predict(best_model, s = best_lambda, newx = X.sites2include)[,1]
  measured = y.sites2include

  toplot = data.frame(site = sites2include, measured = y.sites2include, predicted = predicted) %>% 
    tibble::as_tibble() %>%
    dplyr::left_join(vegetation, by = c("site" = "SampleSiteCode")) %>%
    dplyr::left_join(site2cluster, by = c("site" = "site"))

  min = min(c(y.sites2include, predicted))
  max = max(c(y.sites2include, predicted))


  ggplot(toplot, aes(x=y.sites2include, y=predicted)) + 
    #geom_point(aes(color = factor(`plant-group2`))) +
    geom_point(aes(color = factor(`cluster`))) +
    #geom_point(color = "blue") +
    geom_abline(intercept = 0, slope = 1, color = "gray", linetype = "dashed", size = 0.3) +
    xlim(min, max) + 
    ylim(min, max) + 
    theme_bw() +
    xlab("measured") +
    ylab("predicted")

  conifer_sites = toplot %>% filter(`plant-group2` %in% c("spruce", "fir", "pine")) %>% pull(site)
  aspen_sites = toplot %>% filter(`plant-group2` %in% c("aspen")) %>% pull(site)
  other_sites = setdiff(toplot %>% pull(site), c(conifer_sites, aspen_sites))
  meadow_sites = toplot %>% filter(`Vegetation Type` == "Meadow") %>% pull(site)
  shrub_sites = toplot %>% filter(`Vegetation Type` == "Shrubs") %>% pull(site)
  
  sst_conifer = sum((measured[conifer_sites] - mean(measured[conifer_sites]))^2)
  sse_conifer = sum((predicted[conifer_sites] - measured[conifer_sites])^2)
  rsq_conifer = 1-sse_conifer/sst_conifer
  
  sst_aspen = sum((measured[aspen_sites] - mean(measured[aspen_sites]))^2)
  sse_aspen = sum((predicted[aspen_sites] - measured[aspen_sites])^2)
  rsq_aspen = 1-sse_aspen/sst_aspen
   
  sst_meadow = sum((measured[meadow_sites] - mean(measured[meadow_sites]))^2)
  sse_meadow = sum((predicted[meadow_sites] - measured[meadow_sites])^2)
  rsq_meadow = 1-sse_meadow/sst_meadow
 
  sst_shrub = sum((measured[shrub_sites] - mean(measured[shrub_sites]))^2)
  sse_shrub = sum((predicted[shrub_sites] - measured[shrub_sites])^2)
  rsq_shrub = 1-sse_shrub/sst_shrub
 
  sst <- sum((measured - mean(measured))^2)
  sse <- sum((predicted - measured)^2)
  rsq <- 1 - sse/sst
  lasso.results[cols2predict[c], "r2"] = rsq
  cat(cols2predict[c], "\t", rsq, "\n")
}