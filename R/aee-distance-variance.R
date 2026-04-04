library(dplyr)
#library(lsa) # lsa::cosine
library(coop)
library(GUniFrac)
library(ggplot2)
library(scales)
# --- Config ---
INPUT_CSV  <- "data/locations_with_vectors.csv"
BAND_PREFIX = "A"

# --- Load data ---
# 38472 Gbs
embeddings = read.table(INPUT_CSV, sep = ",", header = T, check.names = F, stringsAsFactors = F) %>% as_tibble()

# Deredundant: keep one row per unique (lat, lon)
# 12445 deredundant locations
embeddings = embeddings %>% distinct(lat, lon, .keep_all = TRUE)

cols = paste0(BAND_PREFIX, sprintf("%02d", 0:63))
# embeddings_small = embeddings %>% slice_head(n = 1000) %>%
#     dplyr::select(all_of(cols)) %>%
#     as.matrix %>% t()
colnames = embeddings %>% pull(biosample)
embeddings_matrix =  embeddings %>%
     dplyr::select(all_of(cols)) %>%
     as.matrix %>% t()
colnames(embeddings) = colnames
alphaearth_cosine_sim <- coop::cosine(embeddings_matrix, inverse = F) %>% as.matrix()
alphaearth_cosine_dist <- 1 - alphaearth_cosine_sim
rownames(alphaearth_cosine_dist) = colnames
colnames(alphaearth_cosine_dist) = colnames

hclust_ae_sites = hclust(as.dist(alphaearth_cosine_dist), method = "ward.D")

nclusters = c(5, 10,15, 18, 20, 40, 60, 80, 100, 150, 200, 250, 500)
ncores = 2
adonis_results_temp =
  parallel::mclapply(1:length(nclusters),
    function(i) {
        cat(i, "\n")
        v = cutree(hclust_ae_sites, nclusters[i])
        site2cluster = data.frame(site = names(v), cluster = paste0("cluster_", v)) %>% as_tibble()
        rownames(site2cluster) = names(v)
        adonis_results = GUniFrac::dmanova(alphaearth_cosine_dist ~ cluster, data = site2cluster)
        #vegan::adonis2(alphaearth_cosine_dist ~ cluster, data = site2cluster, perm = 1)
        result = data.frame(`number of clusters` = nclusters[i],
                            `percent variance` = adonis_results$aov.tab["Model","R2"]*100,
                            #`percent variance` = adonis_results["Model", "R2"]*100,
                            check.names = F)
    },
  mc.cores = 2)

OUTPUT_PNG        <- "data/aee-dist-variance.png"
OUTPUT_HIST_PNG   <- "data/aee-cosine-sim-histogram.png"

adonis_results = do.call("rbind", adonis_results_temp)
colnames(adonis_results) = c("number of clusters", "percent variance")
p = ggplot(adonis_results, aes(x=`number of clusters`, y=`percent variance`)) +
  geom_line() + geom_point(size=5) +
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
  theme(text = element_text(size = 30),
        axis.text = element_text(size = 20, angle = 90))

ggsave(OUTPUT_PNG, plot = p, width = 14, height = 7, dpi = 400)

# --- Histogram of pairwise cosine similarities (upper triangle) ---
sim_vals <- alphaearth_cosine_sim[upper.tri(alphaearth_cosine_sim, diag = FALSE)]

hist_df <- data.frame(cosine_sim = sim_vals)
n_pairs  <- length(sim_vals)

p_hist <- ggplot(hist_df, aes(x = cosine_sim)) +
  geom_density(fill = "#4a90d9", color = "#2c5f8a", alpha = 0.7, linewidth = 0.6) +
  scale_x_continuous(breaks = seq(-1, 1, by = 0.1)) +
  labs(
    title    = sprintf("Pairwise alphaEarth cosine similarities\n(n = %s pairs)", scales::comma(n_pairs)),
    x        = "Cosine similarity",
    y        = "Density"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title   = element_text(hjust = 0.5, face = "bold", size = 18),
    axis.title   = element_text(size = 16),
    axis.text.x  = element_text(size = 13, angle = 45, hjust = 1),
    axis.text.y  = element_text(size = 13),
    panel.grid.minor = element_blank()
  )

ggsave(OUTPUT_HIST_PNG, plot = p_hist, width = 12, height = 7, dpi = 300)
cat(sprintf("Histogram saved to %s\n", OUTPUT_HIST_PNG))
