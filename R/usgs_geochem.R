library(dplyr)
library(readr)
library(tidyr)
library(maps)
library(ggplot2)
library(gridExtra)

usgs_geochem <- read_csv("data/USGS_soil_geochemfeatures.csv")
variable_cols = c("Ag_ppm", "Al_pct", "As_ppm", "Au_ppm", "B_ppm", "Ba_ppm", "Be_ppm", "Bi_ppm", "Br_ppm", "C_pct", "Ca_pct", "CCO3_pct", "Cd_ppm", "Ce_ppm",
"Cl_pct", "CO2_pct", "Co_ppm", "Corg_pct", "Cr_ppm", "Cs_ppm", "Cu_ppm", "Dy_ppm", "Er_ppm", "Eu_ppm", "F_pct", "Fe2O3_pct", "Fe_pct", "FeO_pct", 
"Ga_ppm","Gd_ppm", "Ge_ppm", "Hf_ppm", "Hg_ppm", "Ho_ppm", "I_ppm", "In_ppm", "Ir_ppm", "K_pct", "La_ppm", "Li_ppm", "LOI_pct", "Lu_ppm", "Mg_pct",
"Mn_pct", "Mo_ppm", "Na_pct", "Nb_ppm", "Nd_ppm", "Ni_ppm", "Os_ppm", "P_pct", "Pb_ppm", "Pd_ppm", "Pr_ppm", "Pt_ppm", "Rb_ppm", "Re_ppm", "Rh_ppm",
"Ru_ppm", "S_pct", "Sb_ppm", "Sc_ppm", "Se_ppm", "Si_pct", "Sm_ppm", "Sn_ppm", "SO4_pct", "Sr_ppm", "Sulfide_pct", "Ta_ppm", "Tb_ppm", "Te_ppm",
"Th_ppm", "Ti_pct", "Tl_ppm", "Tm_ppm", "U_ppm", "V_ppm", "W_ppm", "Y_ppm", "Yb_ppm", "Zn_ppm", "Zr_ppm")

us_map <- map_data("state")
results = parallel::mclapply(1:length(variable_cols),
    function(r) {
        cat(r, "\n")
        temp = usgs_geochem %>%
            select(latitude, longitude, variable_cols[r]) %>%
            filter(if_all(variable_cols[r], ~ !is.na(.))) %>%
            count(latitude, longitude)
        title = paste0(variable_cols[r], ":", nrow(temp), " unique locations.")
        map = ggplot() +
            geom_polygon(data = us_map, aes(x = long, y = lat, group = group), fill = "gray90", colour = "gray80") +
            geom_point(data = temp, aes(x = longitude, y = latitude), color = "red", size = 0.5, alpha = 0.5) +
            coord_fixed(1.3, xlim = c(-125, -66), ylim = c(24, 50)) +
            theme_minimal() +
            ggtitle(title)
        returnList = list(counts = c(variable_cols[r], nrow(temp)),
                          plot = map)
    },
    mc.cores = 4)

geochem_nlocs = do.call(rbind, lapply(results, `[[`, 1))
colnames(geochem_nlocs) = c("parameter", "n_uniqlatlong")
write.table(geochem_nlocs, "data/geochem_nlocs.tsv", sep = "\t", quote = FALSE, row.names = FALSE)


all_maps = lapply(results, `[[`, 2)
all_maps_multipage = gridExtra::marrangeGrob(all_maps, nrow = 4, ncol = 2)
ggsave("data/geochem_locs_maps.pdf", dpi = 150, limitsize = T, all_maps_multipage)

