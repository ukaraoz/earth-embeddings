# variance partitioning analysis of NEON sites with 
# different sets of variables
# 1. AE embedding similarity
# 2. Soil biogeo
# 3. Community aggregated traits
# 4. Compare to all guilds (first clustered based on microTrait similarity, 
#    EcoSIM+ guilds

library(easypackages)
packages <- c("dplyr", "tibble", "magrittr", "tidyr", "labdsv", "ggplot2", "phyloseq", "ape",
              "ecodist", "doParallel", "vegan", "DESeq2", "pheatmap", "Hmisc", "GGally", "ggfortify", "ggforce",
              "usedist", "lsa", "dendextend", "ComplexHeatmap", "stringr")
libraries(packages)

base = "/Users/ukaraoz/Work/EastRiver/synthesis"
source(file.path(base, "scripts/misc.R"))
outdir = file.path(base, "results/")

#####
# 1. AE embeddings
#####
alphaearth = read_alphaearth(AE_file = file.path(base, "data/alphaearth/ersynthesis_metadata_alphaEarth.csv"),
                             AE_similarity_file = file.path(base, "data/alphaearth/ersynthesis_metadata_AEF_cosine_similarity_matrix.csv"))
AE_cols = paste0("A", formatC(seq(0,63,1), width = 2, format = "d", flag = "0"))

alphaearth_embeddings = alphaearth$embeddings %>%
  filter(str_detect(sample, "-ER18-SC-SDNA")) %>%
  mutate(sample = str_replace_all(sample, "-SC-SDNA", "")) %>%
  select(sample, AE_cols)

alphaearth_dist = alphaearth$embedding_cosine_dist
alphaearth_cosine = alphaearth$embedding_cosine

alphaearth_neonsamples = grep("ER18", attr(alphaearth_dist, "Labels"), value=T)
alphaearth_cosine_dist_neon = usedist::dist_subset(alphaearth_dist, alphaearth_neonsamples)

#####
# 2. Community aggregated trait matrix
#####
sample2weighted_traits_file = file.path(base, "results/neon-bins_sample2weighted_traits.txt")
# "307-ER18" not in the AE list, remove here
sample2weighted_traits = read.table(sample2weighted_traits_file, check.names = F, header = T, sep = "\t") %>% tibble::as_tibble() %>% 
  filter(SampleSiteCode != "307-ER18") %>%
  tibble::column_to_rownames(var = "SampleSiteCode")

select_cols = names(which(apply(as.matrix(sample2weighted_traits),2,sum) != 0))
sample2weighted_traits = sample2weighted_traits[, select_cols]

traitdata_dist = dist(scale(sample2weighted_traits))

traitdata_dist_rho = dist(1 - abs(cor(t(sample2weighted_traits), method = "spearman")))

#####
# 3. biogeochemical data
#####
biogeoremote_temp = read_biogeoremote(biogeodata_file = file.path(base, "data/biogeo_remote/NEON_soil_bgc_master_11.28.2023.csv"),
                                      remotesensingdata_file = file.path(base, "data/biogeo_remote/NEON_RS_master_updated_11.28.2023.csv"),
                                      leaftraitsdata_file = file.path(base, "data/biogeo_remote/NEON_LeafTraits_master_11.28.2023.csv"),
                                      variables_file = file.path(base, "data/biogeo_remote/biogeo-remote-variables.txt"),
                                      plants_file = file.path(base, "data/biogeo_remote/plant-metadata.txt"))
biogeoremote = biogeoremote_temp$data
biogeoremote_vars = biogeoremote_temp$variables
biogeochemical_vars = biogeoremote_vars %>% 
  dplyr::filter(`variable-source` == "biogeochemical") %>%
  dplyr::pull(`variable-new`)

vegetation_vars = biogeoremote_vars %>%
  dplyr::filter(`variable-source` == "vegetation" & `variable-keep` == 1) %>%
  dplyr::pull(`variable-new`)


hyperspectral_vars = biogeoremote_vars %>% 
  dplyr::filter(`variable-source` == "hyperspectral") %>%
  dplyr::pull(`variable-new`)

lidar_vars = biogeoremote_vars %>% 
  dplyr::filter(`variable-source` == "LiDAR") %>%
  dplyr::filter(!`variable-new` %in% "Drought Sensitivity") %>%
  dplyr::pull(`variable-new`)

# select variables, drop NAs
biogeochemical = biogeoremote %>% 
  dplyr::select(SampleSiteCode, biogeochemical_vars) %>% # 438
  drop_na %>% #385
  column_to_rownames(var = "SampleSiteCode")
biogeochemical_wSamplecolumn = biogeochemical %>% rownames_to_column("SampleSiteCode")

vegetation = biogeoremote %>%
  dplyr::select(SampleSiteCode, vegetation_vars)

hyperspectral = biogeoremote %>% 
  dplyr::select(SampleSiteCode, hyperspectral_vars) %>% # 438
  drop_na %>% #403
  column_to_rownames(var = "SampleSiteCode")
hyperspectral_wSamplecolumn = hyperspectral %>% rownames_to_column("SampleSiteCode")

lidar = biogeoremote %>% 
  dplyr::select(SampleSiteCode, lidar_vars) %>% # 438
  drop_na %>% #426
  column_to_rownames(var = "SampleSiteCode")
lidar_wSamplecolumn = lidar %>% rownames_to_column("SampleSiteCode")

biogeochemical_std = round(decostand(biogeochemical, method = "standardize"), 4)
# 385 x 385
biogeochemical_dist <- round(vegdist(biogeochemical_std, method = "euclidian", upper = TRUE), 4)
biogeochemical_dist_rho = dist(1 - abs(cor(t(biogeochemical_std), method = "spearman")))


hyperspectral_std = round(decostand(hyperspectral, method = "standardize"), 4)
# 403 x 403
hyperspectral_dist <- round(vegdist(hyperspectral_std, method = "euclidian", upper = TRUE), 4)
hyperspectral_dist_rho = dist(1 - abs(cor(t(hyperspectral_std), method = "spearman")))

lidar_std = round(decostand(lidar, method = "standardize"), 4)
# 403 x 403
lidar_dist <- round(vegdist(lidar_std, method = "euclidian", upper = TRUE), 4)
lidar_dist_rho = dist(1 - abs(cor(t(lidar_std), method = "spearman")))


# intersect
# 249
alphaearth_samples = grep("ER18", attr(alphaearth_dist, "Labels"), value=T)
traitdata_samples = rownames(sample2weighted_traits)
biogeochemical_samples = attr(biogeochemical_dist, "Labels")
hyperspectral_samples = attr(hyperspectral_dist, "Labels")
lidar_samples = attr(lidar_dist, "Labels")

variable_sets = c("alphaearth", "traitdata", "biogeochemical", "hyperspectral", "lidar")
mantel_results = data.frame(cbind(t(combn(variable_sets,2)), mantel_r = "", mantel_p = ""))
for(i in 1:nrow(mantel_results)) {
  cat(i, "\n")
  dist1 = get(paste0(mantel_results[i,1], "_dist"))
  dist2 = get(paste0(mantel_results[i,2], "_dist"))
  samples1 = get(paste0(mantel_results[i,1], "_samples"))
  samples2 = get(paste0(mantel_results[i,2], "_samples"))
  mantel_temp = mantel(usedist::dist_subset(dist1, intersect(samples1, samples2)),
                       usedist::dist_subset(dist2, intersect(samples1, samples2)),
                       permutations = 10000)
  mantel_results[i, "mantel_r"] = mantel_temp$statistic
  mantel_results[i, "mantel_p"] = mantel_temp$signif
}
mantel_results %>% tibble::as_tibble() %>% 
  select(-mantel_p) %>% 
  pivot_wider(names_from = V1, values_from = mantel_r)

alphaearth_dist_rho = alphaearth_dist
mantel_results_rho = data.frame(cbind(t(combn(variable_sets,2)), mantel_r = "", mantel_p = ""))
for(i in 1:nrow(mantel_results)) {
  cat(i, "\n")
  dist1 = get(paste0(mantel_results_rho[i,1], "_dist_rho"))
  dist2 = get(paste0(mantel_results_rho[i,2], "_dist_rho"))
  samples1 = get(paste0(mantel_results_rho[i,1], "_samples"))
  samples2 = get(paste0(mantel_results_rho[i,2], "_samples"))
  mantel_temp = mantel(usedist::dist_subset(dist1, intersect(samples1, samples2)),
                       usedist::dist_subset(dist2, intersect(samples1, samples2)),
                       permutations = 10000)
  mantel_results_rho[i, "mantel_r"] = mantel_temp$statistic
  mantel_results_rho[i, "mantel_p"] = mantel_temp$signif
}
mantel_results_rho %>% tibble::as_tibble() %>% 
  select(-mantel_p) %>% 
  pivot_wider(names_from = V1, values_from = mantel_r)

##################
# 2. cluster sites based on AE cosine similarity
##################
alphaearth_dist_neon = usedist::dist_subset(alphaearth_dist, alphaearth_neonsamples)
hclust_ae_sites = hclust(alphaearth_dist_neon, method = "ward.D")

#### 
# 2.1. measure embedding distance variance as a function of clusters
#### 
# vary number of clusters and measure explained variance
nclusters = c(2, 5,10,15, 18, 20, 40, 60, 80, 100, 150, 200, 249)

adonis_results_temp =
  parallel::mclapply(1:length(nclusters),
    function(i) {
       v = cutree(hclust_ae_sites, nclusters[i])
       site2cluster = data.frame(site = names(v), cluster = paste0("cluster_", v)) %>% as_tibble()
       rownames(site2cluster) = names(v)
       adonis_results = vegan::adonis2(alphaearth_dist_neon ~ cluster, data = site2cluster, perm = 100)
       result = data.frame(`number of clusters` = nclusters[i],
                           `percent variance` = adonis_results["Model", "R2"]*100,
                           check.names = F)
    },
  mc.cores = ncores)
adonis_results = do.call("rbind", adonis_results_temp)
p = ggplot(adonis_results, aes(x=`number of clusters`, y=`percent variance`)) +
  geom_line() + geom_point(size=0.8) +
  #xlim(1, 500) +
  geom_hline(yintercept = 60, color = "red", linewidth = 0.3) +
  geom_hline(yintercept = 70, color = "red", linewidth = 0.3) +
  geom_hline(yintercept = 90, color = "red", linewidth = 0.3) +
  #scale_x_continuous(breaks = c(seq(0, attr(alphaearth_dist_neon, "Size")-1, by = 5), attr(alphaearth_dist_neon, "Size"))) +
  scale_x_continuous(breaks = c(0, nclusters)) +
  scale_y_continuous(breaks = seq(0, 100, by = 10)) +
  #ggtitle("variance across clusters") +
  xlab("number of clusters") +
  ylab("% explained variance") +
  theme(text = element_text(size = 18),
        axis.text = element_text(size = 10, angle = 90))
pdf_outfile = file.path(outdir, "embedding-dist-variance-acrossclusters.pdf")
ggsave(p, height = 8, width = 8, filename = pdf_outfile)

#### 
# 2.2. cut for 10 clusters at 80% of variance
# summarize trait means across clusters
####
v = cutree(hclust_ae_sites, 10)
site2cluster = data.frame(site = names(v), cluster = paste0("cluster_", v)) %>% as_tibble()
rownames(site2cluster) = names(v)

# write site2cluster results
site2cluster2write = site2cluster %>%
  dplyr::mutate(cluster = paste0("AEE-cosdist-80_", cluster))
write.table(site2cluster2write, 
            file = file.path(base, "results/site2AEE-cosdist-80_cluster.txt"),
            row.names = F, col.names = T, quote = F, sep = "\t")


# read aggregate trait data again so it is a tibble
sample2alldata = read.table(sample2weighted_traits_file, check.names = F, header = T, sep = "\t") %>% tibble::as_tibble() %>% 
  filter(SampleSiteCode != "307-ER18") %>%
  dplyr::left_join(site2cluster, by = c("SampleSiteCode" = "site")) %>%
  # merge with biogeochemical
  dplyr::left_join(biogeochemical_wSamplecolumn, by = c("SampleSiteCode" = "SampleSiteCode")) %>%
  dplyr::left_join(hyperspectral_wSamplecolumn, by = c("SampleSiteCode" = "SampleSiteCode")) %>%
  dplyr::left_join(lidar_wSamplecolumn, by = c("SampleSiteCode" = "SampleSiteCode"))

  # remove undetected traits
sample2weighted_traits = sample2alldata %>% 
  dplyr::select(c("SampleSiteCode", "cluster", select_cols))

cluster2traitprofile = sample2weighted_traits %>%
    dplyr::select(-c("SampleSiteCode")) %>%
    dplyr::group_by(cluster) %>%
    dplyr::summarise(across(everything(), mean))

cluster2meandata = sample2alldata %>%
    dplyr::select(-c("SampleSiteCode")) %>%
    dplyr::group_by(cluster) %>%
    dplyr::summarise(across(everything(), .f = mean, na.rm = T))

# for each trait, test differences across clusters
trait_matrix_p = matrix(nrow = 0, ncol = 2)
all_traits = colnames(sample2weighted_traits)[3:ncol(sample2weighted_traits)]
for(c in 1:length(all_traits)) {
  temp = sample2weighted_traits %>%
    dplyr::select(c("cluster", all_traits[c])) %>%
    dplyr::mutate(cluster = factor(cluster)) %>%
    as.data.frame
  pvalue = kruskal.test(get(colnames(temp)[2]) ~ cluster, data = temp)$p.value
  # turned off temporarily for kbase
  #effsize = temp %>% rstatix::kruskal_effsize(get(colnames(temp)[2]) ~ guild) %>% dplyr::select(effsize)
  trait_matrix_p = rbind(trait_matrix_p,
                         c(all_traits[c], pvalue))
  #cat(c, "\t", pvalue, "\n")
}
pthreshold = 0.05
colnames(trait_matrix_p) = c("trait", "pvalue")
trait_matrix_p = data.frame(trait_matrix_p,
                            pvalue.star = gtools::stars.pval(as.numeric(trait_matrix_p[,"pvalue"]))) %>%
  dplyr::filter(pvalue <= pthreshold)

toplot = cluster2traitprofile %>% 
  dplyr::select(cluster, trait_matrix_p %>% pull(trait)) %>%
  tibble::column_to_rownames(var = "cluster") %>% 
  as.data.frame()
# sort columns alphabetically
toplot = toplot[, sort(colnames(toplot))]
column_split = sub(":.*", "", colnames(toplot))
colnames(toplot) = sub("Resource Acquisition:|Resource Use:|Stress Tolerance:", "", colnames(toplot))
p_main = ComplexHeatmap::Heatmap(scale(toplot),
                                 col = circlize::colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
                                 #col = colorRampPalette(RColorBrewer::brewer.pal(n = 7, name = "Reds"))(20),
                                 column_split = column_split, column_gap = unit(c(0.3, 0.3), "inches"),
                                 cluster_rows = FALSE, cluster_columns = FALSE,
                                 show_row_names = TRUE, show_column_names = TRUE,
                                 column_names_side = "top",
                                 show_row_dend = TRUE, show_column_dend = FALSE,
                                 show_heatmap_legend = TRUE,
                                 row_title = "satellite embedding clusters",
                                 row_title_gp = grid::gpar(fontsize = 15),
                                 column_title = "community aggregated microbial traits",
                                 column_title_gp = grid::gpar(fontsize = 15),
                                 column_title_side = "bottom", # bottom does overlap with the labels
                                 #row_labels = row_labels,
                                 row_names_side = "left",
                                 row_names_gp = grid::gpar(fontsize = 12),
                                 column_names_gp = grid::gpar(fontsize = 10),
                                 heatmap_height = unit(6, "inches"),
                                 heatmap_width = unit(11, "inches"))
heatmap_outfile = file.path(outdir, "cluster2traitprofile_plot.png")
png(filename = heatmap_outfile, width = 14, height = 22, units = "in", res = 800)
ComplexHeatmap::draw(p_main, padding = unit(c(0, 0, 0, 0), "mm"),
                     column_title_gp=grid::gpar(fontsize=38))
dev.off()

toplot_metadata= cluster2meandata %>% 
  dplyr::select(cluster, biogeochemical_vars, hyperspectral_vars, lidar_vars) %>%
  tibble::column_to_rownames(var = "cluster") %>% 
  as.data.frame()
column_split = c(rep("Soil Biogeochemistry", length(biogeochemical_vars)), 
                 rep("Hyperspectral", length(hyperspectral_vars)),
                 rep("Lidar", length(lidar_vars)))
p_metadata = ComplexHeatmap::Heatmap(scale(toplot_metadata),
                                 col = circlize::colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
                                 #col = colorRampPalette(RColorBrewer::brewer.pal(n = 7, name = "Reds"))(20),
                                 column_split = column_split, column_gap = unit(c(0.3, 0.3), "inches"),
                                 cluster_rows = FALSE, cluster_columns = FALSE,
                                 show_row_names = TRUE, show_column_names = TRUE,
                                 column_names_side = "top",
                                 show_row_dend = TRUE, show_column_dend = FALSE,
                                 show_heatmap_legend = TRUE,
                                 #row_title = "",
                                 row_title_gp = grid::gpar(fontsize = 15),
                                 #column_title = "",
                                 column_title_gp = grid::gpar(fontsize = 15),
                                 column_title_side = "bottom", # bottom does overlap with the labels
                                 #row_labels = row_labels,
                                 row_names_side = "left",
                                 row_names_gp = grid::gpar(fontsize = 12),
                                 column_names_gp = grid::gpar(fontsize = 10),
                                 heatmap_height = unit(6, "inches"),
                                 heatmap_width = unit(11, "inches"))
heatmap_outfile = file.path(outdir, "cluster2metadata_plot.png")
png(filename = heatmap_outfile, width = 14, height = 22, units = "in", res = 800)
ComplexHeatmap::draw(p_metadata, padding = unit(c(0, 0, 0, 0), "mm"),
                     column_title_gp=grid::gpar(fontsize=38))
dev.off()


##################
# 3. regress AEE to predict trait and other data
##################
sample2alldata1 = sample2alldata %>%
  dplyr::left_join(alphaearth_embeddings, by = c("SampleSiteCode" = "sample")) %>%
  dplyr::select(SampleSiteCode, AE_cols, everything())
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
#1:188 microbial
lasso.results.microbial.plot = lasso.results %>% tibble::as_tibble(rownames = "variable") %>%
  dplyr::slice(1:188) %>%
  filter(!is.na(r2)) %>%
  arrange(desc(r2)) %>%
  filter(r2>0.1) %>%
  ggplot(aes(x=reorder(variable, r2), y=r2)) + 
  geom_bar(stat = "identity", fill = "#CC5500") + 
  coord_flip() +
  theme(
    axis.title.y = element_blank(),  # Removes the y-axis title
    axis.text.y = element_text(size = 16),
    axis.text.x = element_text(size = 16)
  )

png(filename = file.path(outdir, "lasso_microbial.png"), width = 22, height = 12, units = "in", res = 200)
print(lasso.results.microbial.plot)
dev.off()

lasso.results.others.plot = lasso.results %>% tibble::as_tibble(rownames = "variable") %>%
  dplyr::slice(189:208) %>%
  filter(!is.na(r2)) %>%
  arrange(desc(r2)) %>%
  filter(r2>0.1) %>%
  ggplot(aes(x=reorder(variable, r2), y=r2)) + 
  geom_bar(stat = "identity", fill = "#CC5500") + 
  coord_flip() +
  theme(
    axis.title.y = element_blank(),  # Removes the y-axis title
    axis.text.y = element_text(size = 16),
    axis.text.x = element_text(size = 16)
  )

png(filename = file.path(outdir, "lasso_others.png"), width = 22, height = 12, units = "in", res = 200)
print(lasso.results.others.plot)
dev.off()

# write to file
towrite = lasso.results %>% tibble::as_tibble(rownames = "variable") %>%
  dplyr::slice(1:188) %>%
  filter(!is.na(r2)) %>%
  arrange(desc(r2))
write.table(towrite, file = file.path(outdir, "AEE-lassoreg_models-R2.microbial.xls"), row.names = F, col.names = T, sep = "\t", quote = F)

towrite = lasso.results %>% tibble::as_tibble(rownames = "variable") %>%
  dplyr::slice(189:236) %>%
  filter(!is.na(r2)) %>%
  arrange(desc(r2))
write.table(towrite, file = file.path(outdir, "AEE-lassoreg_models-R2.others.xls"), row.names = F, col.names = T, sep = "\t", quote = F)


##################
# 4. iTag data - incomplete
##################
library(phyloseq)
load("/Users/ukaraoz/Work/EastRiver/neon_itag/16S/results/batch1and2_2018and2022/neon_itag.16S.phyloseq.RData")

sample_data(phyloseq)[, "sample.id"]

###END###

# Create a list to hold dendrograms
dend_list <- dendlist(hclust_AE_dend, hclust_genomes_dend)
dendlist(hclust_AE_dend, hclust_genomes_dend) %>%
  untangle(method = "step1side") %>% # Find the best alignment layout
  tanglegram()

heatmap(as.matrix(sample2weighted_traits[, colSums]), scale = "col")

alphaearth_cosine_dist_intersect = usedist::dist_subset(alphaearth_cosine_dist, alphaearth_neonsamples)
heatmap(alphaearth_cosine_dist_intersect)



##### TEST AREA
data(varespec)
data(varechem)

# CCA
vare.cca <- cca(varespec, varechem)

veg.dist <- vegdist(varespec) # Bray-Curtis
env.dist <- vegdist(scale(varechem), "euclid")
mantel(veg.dist, env.dist)
mantel(veg.dist, env.dist, method="spear")

## Common but bad way: use all variables you happen to have in your
## environmental data matrix
vare.cca
plot(vare.cca)
## Formula interface and a better model
vare.cca <- cca(varespec ~ Al + P*(K + Baresoil), data=varechem)
vare.cca
plot(vare.cca)
## Partialling out and negative components of variance
cca(varespec ~ Ca, varechem)
cca(varespec ~ Ca + Condition(pH), varechem)
## RDA
data(dune)
data(dune.env)
dune.Manure <- rda(dune ~ Manure, dune.env)
plot(dune.Manure) 

