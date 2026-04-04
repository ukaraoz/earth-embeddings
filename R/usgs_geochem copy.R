library(dplyr)
library(arrow)

base <- "/Users/ukaraoz/Work/earth_embeddings/data"

tbl_parameters = read.csv(file.path(base, "Geochem_DataRelease_2025_11_18.csv/tbl_parameter.csv"), sep = ",") %>% as_tibble() 
grep("aluminum", parameter_desc %>% pull())

# Soluble soil Al content, "Al_pct"     "Al_AM"      "Al_pct_ALL"
Soil Al(OH)3 content
Soil AlPO4 content
Soluble soil Ca content, "Ca_pct"     "Ca_AM"      "Ca_pct_ALL"
Soil CaCO3 content, "CCO3_pct"     "CCO3_AM"      "CCO3_pct_ALL"
Soil CaHPO4 content, -
Soil apatite content, -
Soil CaSO4 content,  "SO4_pct"     "SO4_AM"      "SO4_pct_ALL"
Soluble soil Cl content, "Cl_pct"     "Cl_AM"      "Cl_pct_ALL"

Cation exchange capacity
Soluble soil Fe content
Soil Fe(OH)3 content
Soil FePO4 content
Atmospheric CH4
Soluble soil K content
Soluble soil MG content
Soluble soil Na content
Total soil NH4 concentration
Total soil NO3 concentration
Atmospheric CO2
Total soil organic carbon
Total soil organic nitrogen
Total soil organic phosphorus
POC (part of SOC)
Total soil H2PO4 concentration
Sand content
Silt content
Soluble soil SO4 content

ecosim_parameters = 

file = "USGS_soil_geochemfeatures.csv"

data = read.csv(file.path(base, file), sep = ",") %>% as_tibble()
var_cols = colnames(data)[50:334]


##############
parquet_dir <- file.path(base, "Geochem_DataRelease_2025_11_18.parquet")
csv_dir     <- file.path(base, "Geochem_DataRelease_2025_11_18.csv")

dir.create(csv_dir, showWarnings = FALSE)

parquet_files <- list.files(parquet_dir, pattern = "\\.parquet$", full.names = TRUE)
parquet_files = parquet_files[4:length(parquet_files)]
for (pf in parquet_files) {
  stem    <- tools::file_path_sans_ext(basename(pf))
  out_csv <- file.path(csv_dir, paste0(stem, ".csv"))
  cat(sprintf("Converting %s -> %s\n", basename(pf), out_csv))
  read_parquet(pf) %>% write.csv(out_csv, row.names = FALSE)
}

cat(sprintf("Done. %d files written to %s\n", length(parquet_files), csv_dir))
