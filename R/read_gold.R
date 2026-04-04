library(dplyr)
library(readxl)
library(sf)
library(rnaturalearth)


base = "/Users/ukaraoz/Work/earth_embeddings/data"
goldDataFile = file.path(base, "goldData_032026.xlsx")

Study = read_excel(goldDataFile, sheet = "Study", col_types = "text")
Biosample = read_excel(goldDataFile, sheet = "Biosample", col_types = "text")
Organism = read_excel(goldDataFile, sheet = "Organism",col_types = "text")
SequencingProject = read_excel(goldDataFile, sheet = "Sequencing Project", col_types = "text")
AnalysisProject = read_excel(goldDataFile, sheet = "Analysis Project", col_types = "text")

towrite = SequencingProject %>% 
  dplyr::filter(`SEQUENCING STRATEGY` == "Metagenome") %>% 
  dplyr::left_join(Biosample, by = c("BIOSAMPLE GOLD ID" = "BIOSAMPLE GOLD ID")) %>%
  dplyr::filter(`BIOSAMPLE ECOSYSTEM` == "Environmental" & 
                `BIOSAMPLE ECOSYSTEM CATEGORY` == "Terrestrial" & 
                `BIOSAMPLE ECOSYSTEM TYPE` == "Soil") %>%
  dplyr::select(`BIOSAMPLE GOLD ID`, `BIOSAMPLE NAME`, `BIOSAMPLE LATITUDE`, `BIOSAMPLE LONGITUDE`) %>%
  dplyr::rename(latitude = `BIOSAMPLE LATITUDE`, longitude = `BIOSAMPLE LONGITUDE`) %>%
  dplyr::mutate(latitude = as.numeric(latitude), longitude = as.numeric(longitude))
towrite_kml = towrite %>%
  dplyr::filter(is.na(latitude) == F & is.na(longitude) == F) 
  st_as_sf(coords = c("longitude", "latitude"), # Specify X and Y column names
           crs = 4326)
st_write(towrite_kml, file.path(base, "goldData_032026.kml"), driver = "KML", delete_dsn = TRUE)


write.table(towrite %>% mutate(across(everything(), as.character)), file.path(base, "locations.csv"), sep = ",", quote = T, row.names = F)