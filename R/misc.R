read_alphaearth <- function(AE_file, AE_similarity_file) {
  #base = "/Users/ukaraoz/Work/EastRiver/synthesis"
  #AE_file = file.path(base, "data/alphaearth/ersynthesis_metadata_alphaEarth.csv")
  #AE_similarity_file = file.path(base, "data/alphaearth/ersynthesis_metadata_AEF_cosine_similarity_matrix.csv")

  # Fabio
  # 1. Downloaded the 64-element AlphaEarth vectors for each point along with their latitude and longitude.
  # 2. Computed the pairwise cosine distances to generate a cosine similarity matrix.
  # 3. Constructed a network using the cosine similarities as edge weights, with nodes representing the points.
  # 4. Pruned the network to retain only statistically significant edges (in network science this is known as the disparity filter).
  # 5. Identified clusters of points using the Louvain community detection method (widely used in network science).
  # cluster-->community column
  require(dplyr)
  
  embeddings = read.table(AE_file, sep = ",", header = T, check.names = F, stringsAsFactors = F) %>% as_tibble()
  embedding_cosine = read.table(AE_similarity_file, sep = ",", header = T, comment.char = "", quote = "", check.names = F, stringsAsFactors = F)
  embedding_cosine = round(embedding_cosine, 6)
  # get rid of "-SC-SDNA" for NEON samples
  rownames(embedding_cosine) = sub("-SC-SDNA", "", rownames(embedding_cosine))
  colnames(embedding_cosine) = sub("-SC-SDNA", "", colnames(embedding_cosine))
  embedding_cosine_dist = 1-embedding_cosine

  #alphaearth_embeddings = alphaearth$embeddings %>%
  #  filter(sample %in% paste0(alphaearth_neonsamples, "-SC-SDNA")) %>% 
  #  select(sample, AE_cols) %>% 
  #  mutate(sample = str_replace(sample, "-SC-SDNA", "")) %>%
  #  tibble::column_to_rownames(var = "sample")
  #my_alphaearth_cosine_dist_neon <- lsa::cosine(t(as.matrix(alphaearth_embeddings)))


  return(list(embeddings = embeddings, 
              embedding_cosine = embedding_cosine %>% as.dist(),
              embedding_cosine_dist = embedding_cosine_dist %>% as.dist()))
}


