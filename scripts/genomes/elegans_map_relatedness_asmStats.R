library(cowplot)
library(readr)
library(dplyr)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)

####################################################################################################
####################################################################################################

# WORLD MAP OF COLLECTION SITES

####################################################################################################
####################################################################################################

# Strain list
strains <- readr::read_tsv("../../processed_data/genome_resources/wild_strain_genome_stats.tsv") %>%
  dplyr::select(Strain) %>%
  dplyr::rename(strain = Strain) %>% 
  dplyr::pull()

# Assessing if IPR term is enriched in specific geo locations
geo_initial <- readr::read_tsv("../../processed_data/genome_resources/elegans_isotypes_sampling_geo.tsv")

# Isolation site of each wild strain
geo <- geo_initial %>%
  dplyr::select(isotype, lat, long, geo) %>%
  dplyr::mutate(hifi_sequence = ifelse(isotype %in% strains, "yes", "no")) %>%
  dplyr::mutate(lat = as.numeric(lat),
                long = as.numeric(long))
# SIX STRAINS ARE MISSING LAT AND LONG INFORMATION #################################################################################################################
# TWO OF THEM ARE SELECTED FOR HIFI SEQUENCING

# Blueprint of a world map
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

# Plotting where isotype strains have been collected and which ones were selected for sequencing
MAP <- ggplot() +
  geom_sf(data = world, fill = "white", color = "black", linewidth = 0.2) +
  geom_point(data = geo, aes(x = long, y = lat, color = hifi_sequence), size = 2, alpha = 0.85) +
  scale_color_manual(values = c("yes" = "red", "no" = "black")) +
  coord_sf(expand = FALSE) +
  theme_minimal() +
  labs(x = NULL, y = NULL, color = "Selected for HiFi sequencing") +
  theme(
    legend.title = element_text(size = 11, color = 'black'),
    legend.text = element_text(size = 11, color = 'black'),
    legend.position = "inside",
    legend.position.inside = c(0.1,0.3),
    axis.text = element_blank(),
    panel.grid = element_blank()) 
MAP 


####################################################################################################
####################################################################################################

# RELATEDNESS MATRIX

####################################################################################################
####################################################################################################

library(revgeo)
library(countrycode)

div1 <- c("ECA1693","ECA3493","ECA3025","ECA724","ECA3678","ECA3656","ECA3654",
          "ECA3157","ECA3671","ECA1195","ECA1193","ECA3673","ECA1191","ECA3012",
          "ECA3011","ECA3601","ECA3524","ECA1187","ECA1252","ECA3672","ECA2897",
          "ECA1253","ECA1238","ECA2882","ECA1232","ECA1196","ECA1242","ECA1240",
          "ECA1251","ECA1185","ECA3222","ECA1194","ECA3256","ECA1192","ECA347",
          "ECA2367","ECA3497","ECA1843","ECA363","ECA2360","ECA1821","ECA2043",
          "ECA741","ECA740","ECA3496","QX1211","ECA36","JU3226","ECA742","ECA229",
          "ECA1825","ECA2950","ECA2942","XZ1514","ECA3005","ECA1761","ECA2993",
          "ECA2951","ECA1717","ECA2291")

div2 <- c("ECA1287","ECA730","ECA1237","ECA723","ECA1284","ECA3668","ECA1278",
          "ECA1276","ECA191","ECA722","ECA2475","ECA396","ECA746","ECA1261",
          "ECA1891","ECA1286","ECA1289","ECA1229","ECA1202","ECA3027","ECA1225",
          "ECA1223","ECA1281","ECA745","ECA2251","ECA2443","ECA1212","ECA744",
          "ECA1969","ECA2482","ECA2417","ECA2081","ECA1216","ECA1228","ECA3605",
          "ECA1293","ECA2073","ECA1247","ECA1919","ECA1208","ECA1283","ECA2489",
          "ECA2473","ECA1997","ECA1206")

div3 <- c("ECA2949","ECA2947","ECA2151","ECA1751","ECA1769","ECA2944","ECA1409",
          "ECA2948")

div4 <- c("ECA2183","ECA2187","ECA1793","ECA2960","ECA2970","ECA2929","ECA2336",
          "ECA2941","ECA1885","ECA3006","ECA1851","ECA2982","ECA2988","ECA2334",
          "ECA2956","ECA2963","ECA1925","ECA2281","ECA2952","ECA2109","ECA2946",
          "ECA2945","ECA2974","ECA1725","ECA2159","ECA2155")

div5 <- c("ECA2091","ECA701")

div6 <- c("ECA2199")

div7 <- c("XZ1516","ECA2191","ECA2195","ECA3088","ECA1713","ECA1493")

div8 <- c("ECA369")

dfdiv <- list(
  "RG2" = div1,
  "RG3" = div2,
  "RG4" = div3,
  "RG5" = div4,
  "RG6" = div5,
  "RG7" = div6,
  "RG8" = div7,
  "RG9" = div8) |>
  stack() |>                      # converts list to data frame (values + group)
  dplyr::rename(STRAIN = values,
                RG = ind) %>%
  dplyr::mutate(RG=as.character(RG))

#ECA1693, ECA2968, ECA1186, ECA3038, ECA2949, ECA3025 could be problematic/misclasified/mixtures

# div4 <- c("ECA1825","ECA2950","ECA2942","XZ1514","ECA3005","ECA1761","ECA2993",
#           "ECA2951","ECA1717","ECA2183","ECA2187","ECA1793","ECA2960","ECA2970",
#           "ECA2929","ECA2336","ECA2941","ECA1885")

tree <- ape::read.tree("/vast/eande106/data/c_elegans/WI/tree/20250625/WI.20250625.hard-filter.isotype.min4.tree")

isos <- readr::read_tsv(file="/vast/eande106/data/c_elegans/WI/concordance/20250625/isotype_groups.tsv") %>%
  dplyr::group_by(isotype_ref_strain) %>%
  dplyr::summarise(count=n()) %>%
  dplyr::rename(isotype=isotype_ref_strain) %>%
  dplyr::left_join(dfdiv,by=c("isotype"="STRAIN")) %>%
  dplyr::mutate(RG=ifelse(is.na(RG),"RG1",RG)) 


rn_isos <- c("ECA243","ECA246","ECA248","ECA250","ECA251","ECA259")

#read pairwise similarity estimates
conc <- readr::read_tsv("/vast/eande106/data/c_elegans/WI/concordance/20250625/gtcheck.txt") %>%
  dplyr::filter(i %in% isos$isotype  & j %in% isos$isotype)

#construct symmetric pairwise similarity matrix
concordance_matrix <- conc %>%
  dplyr::mutate(concordance = (sites - discordance) / sites) %>%
  dplyr::select(i, j, concordance) %>%
  dplyr::bind_rows(
    conc %>%
      dplyr::mutate(concordance = (sites - discordance) / sites) %>%
      dplyr::select(i = j, j = i, concordance)
  ) %>%
  tidyr::pivot_wider(names_from = i, values_from = concordance) %>%
  tibble::column_to_rownames("j") %>%
  as.matrix()

all_iso <- colnames(concordance_matrix)[colnames(concordance_matrix) %in% sort(unique(isos$isotype)) | colnames(concordance_matrix) %in% rn_isos]
concordance_matrix <- concordance_matrix[all_iso, all_iso]
diag(concordance_matrix) <- 1

geo_dat <- readr::read_csv("/vast/eande106/projects/Nicolas/c.elegans/RGs/20250625_c_elegans_strain_data.csv") %>%
  dplyr::select(strain,latitude,longitude) %>%
  dplyr::filter(strain %in% isos$isotype) 


geo_dat$state <- revgeo(long   = geo_dat$longitude,
                        lat    = geo_dat$latitude,
                        output = "addr2")

geo_dat$country <- revgeo(long   = geo_dat$longitude,
                          lat    = geo_dat$latitude,
                          output = "addr1")
geo_dat2 <- geo_dat 

geo_dat3 <- geo_dat2 %>% dplyr::select(-state) |> unnest_wider(country) %>% dplyr::select(-housenumber,-street,-city,-zip) %>%
  dplyr::mutate(state=ifelse(is.na(latitude)|is.na(longitude),NA,state)) %>%
  dplyr::mutate(country=ifelse(is.na(latitude)|is.na(longitude),NA,country))

geo_dat3$continent <- countrycode(geo_dat3$country,
                                  origin      = "country.name",
                                  destination = "continent")

geo_dat3$region <- countrycode(geo_dat3$country, "country.name", "region")
geo_dat3$subregion <- countrycode(geo_dat3$country, "country.name", "region23")

geo_colors <- c("Hawaii"="#66C2A5", 
                "Australia"="#FC8D62", 
                "Central America"="#8DA0CB",
                "South America"="#E78AC3", 
                "Africa"="#A6D854", 
                "Caribbean"="#FFD92F",
                "Taiwan" = "#E5C494", 
                "North America" = "#A65628", 
                "Europe" = "#E41A1C",
                "Asia" = "#399EB8", 
                "Pacific" = "#611EA1",
                "New Zealand" = "green",
                "Atlantic" = "purple", 
                "Oceania" ="#DB7779", 
                "Micronesia" = "#E7298A",
                "Indonesia" = "#7570B3", 
                "Malay Archipelago" = "#4110B3", 
                "Central America"="#8DA0CB",
                "South America"="#E78AC3", 
                "North America" = "#A65628", 
                "Cosmopolitan" = "gray30",
                "unknown" = 'grey')


rg_colors <- c("RG1"="#ff0000",
               "RG2"="#0000ff",
               "RG3"="#ff8200",
               "RG4"="#00ff00",
               "RG5"="#ff037e",
               "RG6"="#16537e",
               "RG7"="#6a329f",
               "RG8"="#f3c588",
               "RG9"="#ffdd02")

#                TD3="#00f4c2", KD = "#8b3700", Temperate="#0000ff", TS1="#9fc5e8", AD="#ff8200", TH="#aeb400",FM="#000000")

geo_dat4 <- geo_dat3 %>%
  dplyr::mutate(geo=ifelse(state=="Hawaii",state,
                           ifelse(country %in% names(geo_colors),country,
                                  ifelse(region=="North America",region,
                                         ifelse(subregion %in% names(geo_colors),subregion,continent))))) %>%
  dplyr::mutate(geo=ifelse(is.na(geo),"unknown",geo)) %>%
  dplyr::select(strain,latitude,longitude,geo) 

all_dat <- isos %>% dplyr::left_join(geo_dat4,by=c("isotype"="strain")) %>%
  dplyr::mutate(abslat=abs(latitude), geo_color = geo_colors[geo], rg_color = rg_colors[RG])

annotation_df <- as.data.frame(all_dat)
annotation_df$abslat <- as.numeric(annotation_df$abslat)
rownames(annotation_df) <- all_dat$isotype

# Assume annotation_df$abslat exists and is numeric
abslat_vec <- annotation_df$abslat
names(abslat_vec) <- rownames(annotation_df)
abslat_vec <- abslat_vec[rownames(concordance_matrix)]  # align order

geo_vec <- annotation_df$geo
names(geo_vec) <- rownames(annotation_df)
geo_vec <- geo_vec[rownames(concordance_matrix)] 

rg_vec <-annotation_df$RG
names(rg_vec) <- rownames(annotation_df)
rg_vec <- rg_vec[rownames(concordance_matrix)]

rg_vec <-annotation_df$RG
names(rg_vec) <- rownames(annotation_df)
rg_vec <- rg_vec[rownames(concordance_matrix)]

rg_cols <- annotation_df %>% dplyr::select(RG,rg_color) %>% dplyr::distinct(RG,.keep_all = T)
rg_col_vec <- rg_cols$rg_color
names(rg_col_vec) <- rg_cols$RG

lat_col_fun <- circlize::colorRamp2(
  seq(0, 60, length.out = 5),
  c("#D73027", "#EDC948" , "#FFFFE0", "#4575B4", "#000080"),
)

geo_vec_filtered <- factor(geo_vec, levels = c(setdiff(geo_vec, c("Cosmopolitan", "unknown")),c("Cosmopolitan", "unknown")))

row_annot <- columnAnnotation(
  #abslat = abslat_vec,
  `Geo.` = geo_vec_filtered,
  Group = rg_vec,
  #subgroup = sublin_vec,
  col = list(#abslat = lat_col_fun,
    `Geo.` = geo_colors,
    Group= rg_col_vec),
  #subgroup= sublin_col_vec),
  annotation_legend_param = list(
    #abslat = list(title = "Absolute\nlatitude", ncol = 1),
    Group = list(title= "Relatedness\ngroup", ncol = 2, fontsize = 14, title_gp = grid::gpar(fontsize = 9), labels_gp = grid::gpar(fontsize = 8)),
    #subgroup = list(title= "Tropical\nsubgroup", ncol = 1, fontsize = 14, title_gp = grid::gpar(fontsize = 9), labels_gp = grid::gpar(fontsize = 8)),
    `Geo.` = list(title = "Geographic\nregion", ncol = 2,   title_gp = grid::gpar(fontsize = 9), labels_gp = grid::gpar(fontsize = 8))
  ),
  annotation_name_side = "left",
  na_col = "white"
  #width = grid::unit(15, "cm")
)

bottom_annot <- rowAnnotation(
  `Abs.Lat.` = abslat_vec,
  col = list(`Abs.Lat.` = lat_col_fun),
  annotation_legend_param = list(
    `Abs.Lat.` = list(title = "Absolute\nlatitude", ncol = 1,title_gp = grid::gpar(fontsize = 9), labels_gp = grid::gpar(fontsize = 8))
  ),
  annotation_name_side = "top"
  #height = unit(15, "cm")
)

conc_test <-  conc %>% dplyr::mutate(conc=(sites-discordance)/sites)

phylo_vals <- as.vector(conc_test$conc)
phylo_vals <- phylo_vals[!is.na(phylo_vals)]  # remove NAs

# Set the breakpoints manually
breaks <- c(
  min(phylo_vals),
  quantile(phylo_vals, 0.25),
  median(phylo_vals),
  quantile(phylo_vals, 0.75),
  1
)

phylo_col_fun <- circlize::colorRamp2(
  breaks,
  #seq(min(phylo_vals, na.rm = TRUE), max(phylo_vals, na.rm = TRUE), length.out = 5),
  c("#4575B4", "#87C6C2", "#FFFFE0","#F4D166","#D73027")
)

heatmap_grob <- grid::grid.grabExpr({
  ComplexHeatmap::draw(
    ComplexHeatmap::Heatmap(
      concordance_matrix,
      name = "Genetic\nsimilarity",
      col = phylo_col_fun,
      cluster_rows = T,
      cluster_columns = T,
      clustering_distance_rows = function(m) as.dist(1 - m),
      clustering_distance_columns = function(m) as.dist(1 - m),
      clustering_method_rows = "complete",
      clustering_method_columns = "complete",
      show_row_names = T,
      show_column_names = F,
      row_names_gp = grid::gpar(fontsize = 1),
      #column_names_gp = grid::gpar(fontsize = 1),
      right_annotation = bottom_annot,
      bottom_annotation = row_annot,
      heatmap_legend_param = list(ncol = 2,
                                  title_gp = grid::gpar(fontsize = 9),
                                  labels_gp = grid::gpar(fontsize = 8)),
    ),
    merge_legend = F,
    heatmap_legend_side = "right",
    annotation_legend_side = "bottom"
  )
})

ggmap <- ggplotify::as.ggplot(heatmap_grob)

# ggsave(plot = ggmap, filename = "/vast/eande106/projects/Nicolas/c.elegans/RGs/CE_heatmap_cc_byGeoLat_20251014.png", width =12, height =12,bg = "white",device = "png",units = "in",dpi = 900)
# ggsave(plot = ggmap, filename = "/vast/eande106/projects/Nicolas/c.elegans/RGs/CE_heatmap_cc_byGeoLat_20251014.pdf", width =12, height =12,bg = "white",device = "pdf",units = "in",dpi = 900)


m <- concordance_matrix

# ensure same row/col set and same order
all_ids <- sort(union(rownames(m), colnames(m)))
m <- m[all_ids, all_ids]

# force symmetry where one side is missing
m2 <- m
m2[is.na(m2)] <- t(m2)[is.na(m2)]

# set diagonal to 1 (typical for similarity); adjust if you want NA
diag(m2) <- 1

hc <- hclust(as.dist(1 - m2), method = "complete")
ord <- hc$order

m_ord <- m2[ord, ord]

m_tri <- m_ord
m_tri[upper.tri(m_tri, diag = FALSE)] <- NA 

ids_ord <- rownames(m_ord)   # same as rownames(m_tri)

annotation_df <- as.data.frame(all_dat)
annotation_df$abslat <- as.numeric(annotation_df$abslat)
rownames(annotation_df) <- all_dat$isotype

# reorder annotation_df to match the heatmap order
annotation_df <- annotation_df[ids_ord, , drop = FALSE]

# aligned vectors
abslat_vec <- annotation_df$abslat; names(abslat_vec) <- ids_ord
geo_vec    <- annotation_df$geo;    names(geo_vec)    <- ids_ord
rg_vec     <- annotation_df$RG;     names(rg_vec)     <- ids_ord

# color vectors (no ordering needed, just named mapping)
rg_cols <- annotation_df %>% dplyr::select(RG, rg_color) %>% dplyr::distinct(RG, .keep_all = TRUE)
rg_col_vec <- rg_cols$rg_color
names(rg_col_vec) <- rg_cols$RG

lat_col_fun <- circlize::colorRamp2(
  seq(0, 60, length.out = 5),
  c("#D73027", "#EDC948", "#FFFFE0", "#4575B4", "#000080")
)

geo_vec_filtered <- factor(
  geo_vec,
  levels = c(setdiff(unique(geo_vec), c("Cosmopolitan", "unknown")), "Cosmopolitan", "unknown")
)

# rebuild annotations using the reordered vectors
row_annot <- ComplexHeatmap::columnAnnotation(
  `Geo.` = geo_vec_filtered,
  `Abs.Lat.` = abslat_vec,
  Group  = rg_vec,
  
  col = list(
    `Geo.` = geo_colors,
    Group  = rg_col_vec,
    `Abs.Lat.` = lat_col_fun
  ),
  annotation_legend_param = list(
    Group = list(title = "Relatedness\ngroup", ncol = 2,
                 title_gp = grid::gpar(fontsize = 10,fontface = "bold"), labels_gp = grid::gpar(fontsize = 9)),
    `Geo.` = list(title = "Geographic\nregion", ncol = 2,
                  title_gp = grid::gpar(fontsize = 10,fontface = "bold"), labels_gp = grid::gpar(fontsize = 9)),
    `Abs.Lat.` = list(title = "Absolute\nlatitude", ncol = 1,
                      title_gp = grid::gpar(fontsize = 10,fontface = "bold"), labels_gp = grid::gpar(fontsize = 9))
  ),
  annotation_name_side = "left",
  na_col = "white"
)

bottom_annot <- ComplexHeatmap::rowAnnotation(
  `Abs.Lat.` = abslat_vec,
  col = list(`Abs.Lat.` = lat_col_fun),
  annotation_legend_param = list(
    `Abs.Lat.` = list(title = "Absolute\nlatitude", ncol = 1,
                      title_gp = grid::gpar(fontsize = 9), labels_gp = grid::gpar(fontsize = 8))
  ),
  annotation_name_side = "top"
)


keep_samples <- c(
  "AB1", "CB4852", "CB4856", "CX11254", "CX11264",
  "CX11314", "DL238", "ECA1185", "ECA1186", "ECA1187",
  "ECA1195", "ECA1202", "ECA1208", "ECA1228", "ECA1237",
  "ECA1243", "ECA1255", "ECA1260", "ECA1286", "ECA1287",
  "ECA1409", "ECA1413", "ECA1493", "ECA1693", "ECA1725",
  "ECA1751", "ECA1757", "ECA1761", "ECA1769", "ECA1825",
  "ECA1843", "ECA1851", "ECA1889", "ECA191", "ECA1943",
  "ECA1997", "ECA2081", "ECA2109", "ECA2111", "ECA2151",
  "ECA2187", "ECA2199", "ECA2281", "ECA2336", "ECA2367",
  "ECA2377", "ECA2417", "ECA2452", "ECA2473", "ECA248",
  "ECA250", "ECA251", "ECA2529", "ECA2555", "ECA2561",
  "ECA2565", "ECA2581", "ECA259", "ECA2607", "ECA2676",
  "ECA2948", "ECA2952", "ECA2968", "ECA3005", "ECA3012",
  "ECA3025", "ECA3035", "ECA3088", "ECA347", "ECA36",
  "ECA369", "ECA594", "ECA701", "ECA703", "ECA706",
  "ECA722", "ECA723", "ECA730", "ECA738", "ECA741",
  "ECA742", "ECA768", "ECA923", "ED3017", "ED3049",
  "ED3073", "EG4725", "JT11398", "JU1400", "JU1409",
  "JU1581", "JU2001", "JU2526", "JU258", "JU2600",
  "JU2829", "JU2841", "JU310", "JU311", "JU3128",
  "JU3166", "JU3226", "JU3282", "JU346", "JU394",
  "JU440", "JU775", "JU782", "LKC34", "MY1",
  "MY10", "MY16", "MY2147", "MY2212", "MY23",
  "MY2693", "NIC1", "NIC1119", "NIC1775", "NIC1790",
  "NIC1792", "NIC1810", "NIC195", "NIC199", "NIC2",
  "NIC259", "NIC526", "PX179", "QG4018", "QG4021",
  "QX1791", "QX1794", "RC301", "TWN2530", "WN2001",
  "WN2082", "WN2117", "XZ1513", "XZ1515", "XZ1516"
)



rn <- rownames(m_tri)

row_labels <- rep("", length(rn))
row_labels[rn %in% keep_samples] <- "→"

heatmap_grob2 <- grid::grid.grabExpr({
  ComplexHeatmap::draw(
    ComplexHeatmap::Heatmap(
      m_tri,
      name = "Genetic\nsimilarity",
      col = phylo_col_fun,
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      show_row_names = TRUE,          # must be logical
      row_labels = row_labels,        # <-- put the vector here
      show_column_names = FALSE,
      row_names_side = "left",
      row_names_gp = grid::gpar(fontsize = 12),
      na_col = "white",
      #left_annotation = bottom_annot
      bottom_annotation = row_annot
    ),
    merge_legends=T
  )
})

ggmap2 <- ggplotify::as.ggplot(heatmap_grob2)













####################################################################################################
####################################################################################################

# GENOME ASSEMBLY STATS

####################################################################################################
####################################################################################################
# Reading in genome stats
stats <- readr::read_tsv("../../processed_data/genome_resources/wild_strain_genome_stats.tsv") %>%
  dplyr::select(Strain,`Number of contigs`,`Genome size`, `Contig N50`, `Contig L90`, `Contig N90`, `Coverage`, `BUSCO completeness (genome)`) %>%
  dplyr::rename(strain = Strain, n_contigs = `Number of contigs`, contig_bp = `Genome size`, ctg_N50 = `Contig N50`, ctg_L90 = `Contig L90`, ctg_N90 = `Contig N90`)

plot_df <- stats %>%
  transmute(
    strain,
    contig_bp,
    ctg_N50_mb = ctg_N50 / 1e6,
    ctg_N90_mb = ctg_N90 / 1e6
  ) %>%
  pivot_longer(
    cols = c(contig_bp, ctg_N50_mb, ctg_N90_mb),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = factor(
      metric,
      levels = c("contig_bp", "ctg_N50_mb", "ctg_N90_mb"),
      labels = c("Assembly size (bp)", "Contig N50 (Mb)", "Contig N90 (Mb)")
    )
  )

plot_df <- stats %>% 
  dplyr::transmute(strain,
    Assembly = contig_bp / 1e6,
    N50 = ctg_N50 / 1e6,
    N90 = ctg_N90 / 1e6) %>%
  tidyr::pivot_longer(cols = c(Assembly, N50, N90),names_to = "metric",values_to = "value")

av_gasm <- plot_df %>% dplyr::filter(metric == "Assembly") %>% dplyr::summarize(mean_asm = mean(value))

top <- ggplot(plot_df, aes(x = metric, y = value)) +
  geom_point(position = position_jitter(width = 0.25), size = 0.3, color = 'black') +
  geom_boxplot(aes(fill = metric),  outlier.shape = NA, alpha = 0.5, fatten = 2) +
  scale_fill_manual(values = c("Assembly" = "blue", "N50" = "green", "N90" = "orange")) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"),
    legend.position = 'none',
    axis.title = element_blank(),
    axis.text.y  = element_text(size = 10, color = 'black'),
    axis.text.x = element_blank(),
    plot.margin = margin(l = 20, t = 2, r =2, b = 7.5),
    axis.ticks.x  = element_blank()
  ) +
  labs(y = "Megabases") +
  coord_cartesian(ylim = c(100,115))
top

av_N50 <- plot_df %>% dplyr::filter(metric == "N50") %>% dplyr::summarize(mean_n50 = mean(value))
av_N90 <- plot_df %>% dplyr::filter(metric == "N90") %>% dplyr::summarize(mean_n90 = mean(value))

bottom <- ggplot(plot_df, aes(x = metric, y = value)) +
  geom_point(position = position_jitter(width = 0.25), size = 0.3, color = 'black') +
  geom_boxplot(aes(fill = metric), outlier.shape = NA, alpha = 0.5, fatten = 2) +
  scale_fill_manual(values = c("Asembly" = "blue", "N50" = "seagreen", "N90" = "orange")) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
    axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"),
    legend.position = 'none',
    plot.margin = margin(l = 20, b = 2, r = 2, t = 7.5),
    axis.text.x  = element_text(size = 10, color = 'black'),
    axis.title.y = element_blank(),
    axis.text.y  = element_text(size = 10, color = 'black')
  ) +
  labs(y = "Megabases") +
  coord_cartesian(ylim = c(0,20))
bottom

aligned <- cowplot::align_plots(top, bottom, align = "v", axis = "lr")

base <- cowplot::plot_grid(cowplot::plot_grid(
  aligned[[1]],aligned[[2]],
  nrow = 2) + draw_label("Megabases", x=0.007, y=0.5, vjust= 1.5, angle=90, size = 12, color = 'black'))

STATS <- cowplot::ggdraw(base) +
  draw_line(x = c(0.145, 0.16), y = c(0.535, 0.540), size = 0.6) +
  draw_line(x = c(0.145, 0.16), y = c(0.47, 0.475), size = 0.6)
STATS






# Concatenate to make final plot
bottom_row <- cowplot::plot_grid(MATRIX, STATS, ncol = 2, rel_widths = c(1,1), labels = c("b","c"))

final_plot <- cowplot::plot_grid(
  MAP, bottom_row,
  ncol = 1,
  rel_heights = c(1,1),
  labels = c("a")
)



# Save the plot
# ggsave("../../figures/strain_selection_genome_stats.png", final_plot, width = 7.5, height = 9, dpi = 600)






