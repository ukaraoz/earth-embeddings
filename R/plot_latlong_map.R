library(ggplot2)
library(maps)
library(dplyr)

# --- Config ---
INPUT_CSV  <- "data/locations.csv"
OUTPUT_PNG <- "data/locations_map.png"
OUTPUT_KML <- "data/locations.kml"
LAT_COL    <- "latitude"
LON_COL    <- "longitude"
POINT_SIZE <- 1.2
POINT_COLOR <- "#e63946"
POINT_ALPHA <- 0.6

# --- Load data ---
locs <- read.csv(INPUT_CSV, stringsAsFactors = FALSE)


if (!all(c(LAT_COL, LON_COL) %in% colnames(locs))) {
  stop(sprintf("Expected columns '%s' and '%s' not found in %s.\nFound: %s",
               LAT_COL, LON_COL, INPUT_CSV, paste(colnames(locs), collapse = ", ")))
}

locs <- locs %>%
  rename(lat = all_of(LAT_COL), lon = all_of(LON_COL)) %>%
  filter(!is.na(lat), !is.na(lon),
         lat >= -90, lat <= 90,
         lon >= -180, lon <= 180)

# --- Aggregate by lat/lon ---
other_cols <- setdiff(colnames(locs), c("lat", "lon"))

locs <- locs %>%
  group_by(lat, lon) %>%
  summarise(
    n = n(),
    across(all_of(other_cols), ~ paste(unique(.x), collapse = "; ")),
    .groups = "drop"
  )

cat(sprintf("Plotting %d unique locations (aggregated from %s)\n", nrow(locs), INPUT_CSV))

# --- World base map ---
world <- map_data("world")

p <- ggplot() +
  geom_polygon(data = world,
               aes(x = long, y = lat, group = group),
               fill = "#d9d9d9", color = "#aaaaaa", linewidth = 0.2) +
  geom_point(data = locs,
             aes(x = lon, y = lat, size = n),
             color = POINT_COLOR, alpha = POINT_ALPHA) +
  scale_size_continuous(name = "count", range = c(POINT_SIZE, POINT_SIZE * 5)) +
  coord_fixed(ratio = 1.3, xlim = c(-180, 180), ylim = c(-90, 90), expand = FALSE) +
  labs(
    title = sprintf("Sample locations (%d unique, %d total)", nrow(locs), sum(locs$n)),
    x = "Longitude", y = "Latitude"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_line(color = "#eeeeee"),
    panel.background = element_rect(fill = "#f0f4f8", color = NA),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave(OUTPUT_PNG, plot = p, width = 28, height = 14, dpi = 400)
cat(sprintf("Map saved to %s\n", OUTPUT_PNG))

# --- Write KML ---
desc_cols <- setdiff(colnames(locs), c("lat", "lon"))

placemarks <- vapply(seq_len(nrow(locs)), function(i) {
  row  <- locs[i, ]
  name <- if ("BIOSAMPLE.GOLD.ID" %in% colnames(locs)) row[["BIOSAMPLE.GOLD.ID"]] else sprintf("loc_%d", i)
  desc <- paste(
    sapply(desc_cols, function(col) sprintf("<b>%s:</b> %s", col, row[[col]])),
    collapse = "<br/>"
  )
  sprintf(
    '    <Placemark>\n      <name>%s</name>\n      <description><![CDATA[%s]]></description>\n      <Point><coordinates>%.6f,%.6f,0</coordinates></Point>\n    </Placemark>',
    name, desc, row[["lon"]], row[["lat"]]
  )
}, character(1))

kml <- c(
  '<?xml version="1.0" encoding="UTF-8"?>',
  '<kml xmlns="http://www.opengis.net/kml/2.2">',
  '  <Document>',
  sprintf('    <name>%s</name>', basename(INPUT_CSV)),
  placemarks,
  '  </Document>',
  '</kml>'
)

writeLines(kml, OUTPUT_KML)
cat(sprintf("KML saved to %s (%d placemarks)\n", OUTPUT_KML, nrow(locs)))
