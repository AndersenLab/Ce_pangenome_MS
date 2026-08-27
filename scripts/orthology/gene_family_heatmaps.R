library(readr)
library(dplyr)
library(ggplot2)
library(ggtree)
library(ComplexHeatmap)
library(grid)

# ======================================================================================================================================================================================== #
# Load in orthogroup matrix and classify gene sets
# ======================================================================================================================================================================================== #
og_matrix <- readr::read_tsv("../../processed_data/orthology/OG_relations_matrix_count.tsv")

# ======================================================================================================================================================================================== #
# Looking at how the accesory genome clusters with strain relatedness 
# ======================================================================================================================================================================================== #
strainCol_c2_u <- setdiff(colnames(og_matrix), "Orthogroup")

# 1 / 142 for private gene set frequency
private_freq <- (1 / length(strainCol_c2_u))

# Gene set classification
pav_mat <- og_matrix %>% dplyr::mutate(across(2:(ncol(.)), ~ ifelse(. >= 1, 1, .))) %>% 
  dplyr::mutate(across(2:ncol(.), ~ifelse(is.na(.), 0, .))) %>%
  dplyr::mutate(sum = rowSums(across(-1, ~ ., .names = NULL), na.rm = TRUE)) %>%
  dplyr::mutate(freq = (sum / length(strainCol_c2_u))) %>%
  dplyr::mutate(
    class = case_when(
      freq == 1 ~ "core",
      freq > private_freq & freq < 1 ~ "accessory",
      freq == private_freq ~ "private",
      TRUE ~ "undefined")) %>%
  dplyr::select(-freq) %>%
  dplyr::mutate(class = factor(class, levels = c("core", "accessory", "private"))) %>%
  dplyr::arrange(desc(sum)) %>%
  dplyr::select(-sum) %>%
  dplyr::mutate(Orthogroup = factor(Orthogroup, levels = unique(Orthogroup)))

# ======================================================================================================================================================================================== #
# Prepare matrix for gheatmap - needs to be a matrix with rownames = tip labels
# ======================================================================================================================================================================================== #
# Convert to wide matrix format with strains as rows
pav_matrix_for_gheatmap <- pav_mat %>%
  dplyr::select(-class) %>%
  tidyr::pivot_longer(cols = all_of(strainCol_c2_u), names_to = "strain", values_to = "presence") %>%
  tidyr::pivot_wider(names_from = Orthogroup, values_from = presence) %>%
  tibble::column_to_rownames("strain") %>%
  as.matrix()

# Reorder columns by gene class (core -> accessory -> private)
og_order <- pav_mat %>% dplyr::pull(Orthogroup) %>% as.character()
pav_matrix_for_gheatmap <- pav_matrix_for_gheatmap[, og_order]

# ======================================================================================================================================================================================== #
# Load tree
# ======================================================================================================================================================================================== #
tree_file <- "../../processed_data/genome_resources/trees/BUSCO_supermatrix.fa.contree"
busco_tree <- read.tree(tree_file)
busco_tree_scaled <- ape::compute.brlen(busco_tree, method = "Grafen")

# ======================================================================================================================================================================================== #
# Order strains according to BUSCO-inferred phylogeny
# ======================================================================================================================================================================================== #
# Get strain order from tree (visual order, top to bottom)
strain_order_by_y <- ggtree(busco_tree_scaled)$data %>%
  dplyr::filter(isTip) %>%
  dplyr::arrange(desc(y)) %>%
  dplyr::pull(label)

# Reorder matrix rows to match tree order
pav_matrix_ordered <- pav_matrix_for_gheatmap[strain_order_by_y, ]

# Verify the strain ordering is correct
identical(rownames(pav_matrix_ordered), strain_order_by_y)  # Should be TRUE

# ======================================================================================================================================================================================== #
# Create annotation for columns (orthogroups)
# ======================================================================================================================================================================================== #
col_annotation <- pav_mat %>%
  dplyr::select(Orthogroup, class) %>%
  dplyr::distinct() %>%
  tibble::column_to_rownames("Orthogroup")

# OG name matching
col_annotation <- col_annotation[colnames(pav_matrix_ordered), , drop = FALSE]
ann_colors <- list(class = c("core" = "green4", "accessory" = "#DB6333", "private" = "magenta3"))

# ======================================================================================================================================================================================== #
# Create annotation for row (collection site)
# ======================================================================================================================================================================================== #
geo_initial <- readr::read_tsv("../../processed_data/genome_resources/isotypes/elegans_isotypes_sampling_geo.tsv")
hawaii_islands <- readr::read_tsv("../../processed_data/genome_resources/isotypes/elegans_isotypes_sampling_geo_hawaii_islands.tsv") %>% dplyr::select(isotype,collection_island_Hawaii)
WSs <- readr::read_tsv("../../tables/wild_strain_genome_stats.tsv") %>% dplyr::select(Strain) %>% dplyr::rename(strain = Strain) %>% dplyr::pull()

# Adding Hawaiian island resolution
geo <- geo_initial %>%
  dplyr::left_join(hawaii_islands, by = "isotype") %>%
  dplyr::mutate(geo = ifelse(geo == "Hawaii",collection_island_Hawaii,geo)) %>%
  dplyr::select(isotype, lat, long, geo) %>%
  dplyr::filter(isotype %in% WSs) %>%
  dplyr::select(strain = isotype, Geo = geo) %>%
  dplyr::add_row(strain = "N2", Geo = "Europe") %>%
  dplyr::add_row(strain = "CGC1", Geo = "Europe")

row_annotation <- geo %>%
  tibble::column_to_rownames("strain")

# Make sure row annotation matches matrix row order
row_annotation <- row_annotation[rownames(pav_matrix_ordered), , drop = FALSE]

# Check for any missing strains
setdiff(rownames(pav_matrix_ordered), geo$strain)  # Should be character(0)

# ======================================================================================================================================================================================== #
# Combine ALL annotation colors into ONE list
# ======================================================================================================================================================================================== #
ann_colors <- list(
  # Column annotation colors
  class = c(
    "core" = "green4", 
    "accessory" = "#DB6333", 
    "private" = "magenta3"
  ),
  # Row annotation colors
  Geo = c(
    "Asia" = "lightblue",
    "South America" = "forestgreen",
    "Indian Ocean" = "chocolate",
    "Big Island" = "black", 
    "Molokai" = "#66C2A5", 
    "Maui" = "yellow", 
    "Oahu" = "brown", 
    "Kauai" = "purple", 
    "Africa" = "green", 
    "North America" = "pink", 
    "Europe" = "#E41A1C", 
    "Atlantic" = "blue", 
    "Oceania" = "cyan", 
    "unknown" = "gray"))

# ======================================================================================================================================================================================== #
# Looking at heatmap of all orthogroups
# ======================================================================================================================================================================================== #
# Column annotation (top)
col_ha <- HeatmapAnnotation(
  class = col_annotation$class,
  col = list(class = ann_colors$class),
  show_legend = FALSE,
  show_annotation_name = FALSE)

# Row annotation (right side)
row_ha <- ComplexHeatmap::rowAnnotation(
  Geo = row_annotation$Geo,
  col = list(Geo = ann_colors$Geo),
  show_legend = TRUE,
  show_annotation_name = FALSE,
  annotation_legend_param = list(
    Geo = list(
      title_gp = grid::gpar(fontface = "plain"))))

pav_heatmap <- ComplexHeatmap::Heatmap(
  pav_matrix_ordered,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = FALSE,
  col = c("0" = "white", "1" = "black"),
  top_annotation = col_ha,
  right_annotation = row_ha,
  row_names_gp = grid::gpar(fontsize = 5),
  column_title = "Orthogroups",
  column_title_side = "bottom",
  
  border = TRUE,
  show_heatmap_legend = FALSE)

ComplexHeatmap::draw(pav_heatmap, padding = grid::unit(c(5, 5, 5, 5), "mm"))















# ======================================================================================================================================================================================== #
# Looking at GPCR clustering
# ======================================================================================================================================================================================== #
# Load in all genes and the orthogroups they contribute to
all_genes_ogs <- readr::read_tsv("../../processed_data/orthology/all_genes_class_OGs.tsv")

# Load in GPCR genes among all wild strains
gpcrs <- readr::read_tsv("../../processed_data/genome_resources/annotation/140_wild_strains_IPR_gpcrs.tsv") %>%
  dplyr::left_join(all_genes_ogs, by = c("strain","gene")) %>% dplyr::select(-class) %>% dplyr::distinct(Orthogroup) %>% dplyr::pull(Orthogroup)


# Filtering for GPCR-containing orthogroups
pav_mat_GPCRS <- pav_mat %>% dplyr::filter(Orthogroup %in% gpcrs) 

pav_matrix_for_gheatmap_GPCR <- pav_mat_GPCRS %>%
  dplyr::select(-class) %>%
  tidyr::pivot_longer(cols = all_of(strainCol_c2_u), names_to = "strain", values_to = "presence") %>%
  tidyr::pivot_wider(names_from = Orthogroup, values_from = presence) %>%
  tibble::column_to_rownames("strain") %>%
  as.matrix()

# Reorder columns by gene class (core -> accessory -> private)
og_order <- pav_mat_GPCRS %>% dplyr::pull(Orthogroup) %>% as.character()
pav_matrix_for_gheatmap_GPCR <- pav_matrix_for_gheatmap_GPCR[, og_order]

# ======================================================================================================================================================================================== #
# Order strains according to BUSCO-inferred phylogency
# ======================================================================================================================================================================================== #
# Reorder matrix rows to match tree order
pav_matrix_ordered_GPCRs <- pav_matrix_for_gheatmap_GPCR[strain_order_by_y, ]

# Verify the strain ordering is correct
identical(rownames(pav_matrix_ordered_GPCRs), strain_order_by_y)  # Should be TRUE

# ======================================================================================================================================================================================== #
# Create annotation for columns (orthogroups)
# ======================================================================================================================================================================================== #
col_annotation_GPCRs<- pav_mat_GPCRS %>%
  dplyr::select(Orthogroup, class) %>%
  dplyr::distinct() %>%
  tibble::column_to_rownames("Orthogroup")

# OG name matching
col_annotation_GPCRs <- col_annotation_GPCRs[colnames(pav_matrix_ordered_GPCRs), , drop = FALSE]

# ======================================================================================================================================================================================== #
# Looking at heatmap
# ======================================================================================================================================================================================== #
# Column annotation (top)
col_ha_GPCRs <- HeatmapAnnotation(
  class = col_annotation_GPCRs$class,
  col = list(class = ann_colors$class),
  show_legend = FALSE,
  show_annotation_name = FALSE)

pav_heatmap_GPCRs <- ComplexHeatmap::Heatmap(
  pav_matrix_ordered_GPCRs,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = FALSE,
  col = c("0" = "white", "1" = "black"),
  top_annotation = col_ha_GPCRs,
  right_annotation = row_ha,
  row_names_gp = grid::gpar(fontsize = 5),
  column_title = "GPCR orthogroups",
  column_title_side = "bottom",
  
  border = TRUE,
  show_heatmap_legend = FALSE)

ComplexHeatmap::draw(pav_heatmap_GPCRs, padding = grid::unit(c(5, 5, 5, 5), "mm"))






# ======================================================================================================================================================================================== #
# Looking at F-box clustering
# ======================================================================================================================================================================================== #
# Load in GPCR genes among all wild strains
fbox <- readr::read_tsv("../../processed_data/genome_resources/annotation/140_wild_strains_IPR_fBox.tsv") %>%
  dplyr::left_join(all_genes_ogs, by = c("strain","gene")) %>% dplyr::select(-class) %>% dplyr::distinct(Orthogroup) %>% dplyr::pull(Orthogroup)

# Filtering for GPCR-containing orthogroups
pav_mat_FBOX <- pav_mat %>% dplyr::filter(Orthogroup %in% fbox) 

pav_matrix_for_gheatmap_FBOX <- pav_mat_FBOX %>%
  dplyr::select(-class) %>%
  tidyr::pivot_longer(cols = all_of(strainCol_c2_u), names_to = "strain", values_to = "presence") %>%
  tidyr::pivot_wider(names_from = Orthogroup, values_from = presence) %>%
  tibble::column_to_rownames("strain") %>%
  as.matrix()

# Reorder columns by gene class (core -> accessory -> private)
og_order <- pav_mat_FBOX %>% dplyr::pull(Orthogroup) %>% as.character()
pav_matrix_for_gheatmap_FBOX <- pav_matrix_for_gheatmap_FBOX[, og_order]

# ======================================================================================================================================================================================== #
# Order strains according to BUSCO-inferred phylogency
# ======================================================================================================================================================================================== #
# Reorder matrix rows to match tree order
pav_matrix_ordered_FBOX <- pav_matrix_for_gheatmap_FBOX[strain_order_by_y, ]

# Verify the strain ordering is correct
identical(rownames(pav_matrix_ordered_FBOX), strain_order_by_y)  # Should be TRUE

# ======================================================================================================================================================================================== #
# Create annotation for columns (orthogroups)
# ======================================================================================================================================================================================== #
col_annotation_FBOX<- pav_mat_FBOX %>%
  dplyr::select(Orthogroup, class) %>%
  dplyr::distinct() %>%
  tibble::column_to_rownames("Orthogroup")

# OG name matching
col_annotation_FBOX <- col_annotation_FBOX[colnames(pav_matrix_ordered_FBOX), , drop = FALSE]

# ======================================================================================================================================================================================== #
# Looking at heatmap 
# ======================================================================================================================================================================================== #
# Column annotation (top)
col_ha_FBOX <- HeatmapAnnotation(
  class = col_annotation_FBOX$class,
  col = list(class = ann_colors$class),
  show_legend = FALSE,
  show_annotation_name = FALSE)

pav_heatmap_FBOX <- ComplexHeatmap::Heatmap(
  pav_matrix_ordered_FBOX,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = FALSE,
  col = c("0" = "white", "1" = "black"),
  top_annotation = col_ha_FBOX,
  right_annotation = row_ha,
  row_names_gp = grid::gpar(fontsize = 5),
  column_title = "F-box orthogroups",
  column_title_side = "bottom",
  
  border = TRUE,
  show_heatmap_legend = FALSE)

ComplexHeatmap::draw(pav_heatmap_FBOX, padding = grid::unit(c(5, 5, 5, 5), "mm"))






# ======================================================================================================================================================================================== #
# Looking at C-type Lectins clustering
# ======================================================================================================================================================================================== #
# Load in GPCR genes among all wild strains
lectins <- readr::read_tsv("../../processed_data/genome_resources/annotation/140_wild_strains_IPR_CtypeLectins.tsv") %>%
  dplyr::left_join(all_genes_ogs, by = c("strain","gene")) %>% dplyr::select(-class) %>% dplyr::distinct(Orthogroup) %>% dplyr::pull(Orthogroup)

# Filtering for GPCR-containing orthogroups
pav_mat_CLECTIN <- pav_mat %>% dplyr::filter(Orthogroup %in% lectins) 

pav_matrix_for_gheatmap_CLECTIN <- pav_mat_CLECTIN %>%
  dplyr::select(-class) %>%
  tidyr::pivot_longer(cols = all_of(strainCol_c2_u), names_to = "strain", values_to = "presence") %>%
  tidyr::pivot_wider(names_from = Orthogroup, values_from = presence) %>%
  tibble::column_to_rownames("strain") %>%
  as.matrix()

# Reorder columns by gene class (core -> accessory -> private)
og_order <- pav_mat_CLECTIN %>% dplyr::pull(Orthogroup) %>% as.character()
pav_matrix_for_gheatmap_CLECTIN <- pav_matrix_for_gheatmap_CLECTIN[, og_order]

# ======================================================================================================================================================================================== #
# Order strains according to BUSCO-inferred phylogency
# ======================================================================================================================================================================================== #
# Reorder matrix rows to match tree order
pav_matrix_ordered_CLECTIN <- pav_matrix_for_gheatmap_CLECTIN[strain_order_by_y, ]

# Verify the strain ordering is correct
identical(rownames(pav_matrix_ordered_CLECTIN), strain_order_by_y)  # Should be TRUE

# ======================================================================================================================================================================================== #
# Create annotation for columns (orthogroups)
# ======================================================================================================================================================================================== #
col_annotation_CLECTIN<- pav_mat_CLECTIN %>%
  dplyr::select(Orthogroup, class) %>%
  dplyr::distinct() %>%
  tibble::column_to_rownames("Orthogroup")

# OG name matching
col_annotation_CLECTIN <- col_annotation_CLECTIN[colnames(pav_matrix_ordered_CLECTIN), , drop = FALSE]

# ======================================================================================================================================================================================== #
# Looking at heatmap 
# ======================================================================================================================================================================================== #
# Column annotation (top)
col_ha_CLECTIN <- HeatmapAnnotation(
  class = col_annotation_CLECTIN$class,
  col = list(class = ann_colors$class),
  show_legend = FALSE,
  show_annotation_name = FALSE)

pav_heatmap_CLECTIN <- ComplexHeatmap::Heatmap(
  pav_matrix_ordered_CLECTIN,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = FALSE,
  col = c("0" = "white", "1" = "black"),
  top_annotation = col_ha_CLECTIN,
  right_annotation = row_ha,
  row_names_gp = grid::gpar(fontsize = 5),
  column_title = "C-type lectin orthogroups",
  column_title_side = "bottom",
  
  border = TRUE,
  show_heatmap_legend = FALSE)

ComplexHeatmap::draw(pav_heatmap_CLECTIN, padding = grid::unit(c(5, 5, 5, 5), "mm"))



# ======================================================================================================================================================================================== #
# Looking at NHR clustering
# ======================================================================================================================================================================================== #
# Load in GPCR genes among all wild strains
nhrs <- readr::read_tsv("../../processed_data/genome_resources/annotation/140_wild_strains_IPR_nhr.tsv") %>%
  dplyr::left_join(all_genes_ogs, by = c("strain","gene")) %>% dplyr::select(-class) %>% dplyr::distinct(Orthogroup) %>% dplyr::pull(Orthogroup)

# Filtering for GPCR-containing orthogroups
pav_mat_NHR <- pav_mat %>% dplyr::filter(Orthogroup %in% nhrs) 

pav_matrix_for_gheatmap_NHR <- pav_mat_NHR %>%
  dplyr::select(-class) %>%
  tidyr::pivot_longer(cols = all_of(strainCol_c2_u), names_to = "strain", values_to = "presence") %>%
  tidyr::pivot_wider(names_from = Orthogroup, values_from = presence) %>%
  tibble::column_to_rownames("strain") %>%
  as.matrix()

# Reorder columns by gene class (core -> accessory -> private)
og_order <- pav_mat_NHR %>% dplyr::pull(Orthogroup) %>% as.character()
pav_matrix_for_gheatmap_NHR <- pav_matrix_for_gheatmap_NHR[, og_order]

# ======================================================================================================================================================================================== #
# Order strains according to BUSCO-inferred phylogency
# ======================================================================================================================================================================================== #
# Reorder matrix rows to match tree order
pav_matrix_ordered_NHR <- pav_matrix_for_gheatmap_NHR[strain_order_by_y, ]

# Verify the strain ordering is correct
identical(rownames(pav_matrix_ordered_NHR), strain_order_by_y)  # Should be TRUE

# ======================================================================================================================================================================================== #
# Create annotation for columns (orthogroups)
# ======================================================================================================================================================================================== #
col_annotation_NHR<- pav_mat_NHR %>%
  dplyr::select(Orthogroup, class) %>%
  dplyr::distinct() %>%
  tibble::column_to_rownames("Orthogroup")

# OG name matching
col_annotation_NHR <- col_annotation_NHR[colnames(pav_matrix_ordered_NHR), , drop = FALSE]

# ======================================================================================================================================================================================== #
# Looking at heatmap
# ======================================================================================================================================================================================== #
# Column annotation (top)
col_ha_NHR <- HeatmapAnnotation(
  class = col_annotation_NHR$class,
  col = list(class = ann_colors$class),
  show_legend = FALSE,
  show_annotation_name = FALSE)

pav_heatmap_NHR <- ComplexHeatmap::Heatmap(
  pav_matrix_ordered_NHR,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = FALSE,
  col = c("0" = "white", "1" = "black"),
  top_annotation = col_ha_NHR,
  right_annotation = row_ha,
  row_names_gp = grid::gpar(fontsize = 5),
  column_title = "NHR orthogroups",
  column_title_side = "bottom",
  
  border = TRUE,
  show_heatmap_legend = FALSE)

ComplexHeatmap::draw(pav_heatmap_NHR, padding = grid::unit(c(5, 5, 5, 5), "mm"))



# ======================================================================================================================================================================================== #
# Looking at cytochrome P450s clustering
# ======================================================================================================================================================================================== #
# Load in GPCR genes among all wild strains
cytos <- readr::read_tsv("../../processed_data/genome_resources/annotation/140_wild_strains_IPR_cytochromeP450.tsv") %>%
  dplyr::left_join(all_genes_ogs, by = c("strain","gene")) %>% dplyr::select(-class) %>% dplyr::distinct(Orthogroup) %>% dplyr::pull(Orthogroup)

# Filtering for GPCR-containing orthogroups
pav_mat_CCYTO <- pav_mat %>% dplyr::filter(Orthogroup %in% cytos) 

pav_matrix_for_gheatmap_CCYTO <- pav_mat_CCYTO %>%
  dplyr::select(-class) %>%
  tidyr::pivot_longer(cols = all_of(strainCol_c2_u), names_to = "strain", values_to = "presence") %>%
  tidyr::pivot_wider(names_from = Orthogroup, values_from = presence) %>%
  tibble::column_to_rownames("strain") %>%
  as.matrix()

# Reorder columns by gene class (core -> accessory -> private)
og_order <- pav_mat_CCYTO %>% dplyr::pull(Orthogroup) %>% as.character()
pav_matrix_for_gheatmap_CCYTO <- pav_matrix_for_gheatmap_CCYTO[, og_order]

# ======================================================================================================================================================================================== #
# Order strains according to BUSCO-inferred phylogency
# ======================================================================================================================================================================================== #
# Reorder matrix rows to match tree order
pav_matrix_ordered_CCYTO <- pav_matrix_for_gheatmap_CCYTO[strain_order_by_y, ]

# Verify the strain ordering is correct
identical(rownames(pav_matrix_ordered_CCYTO), strain_order_by_y)  # Should be TRUE

# ======================================================================================================================================================================================== #
# Create annotation for columns (orthogroups)
# ======================================================================================================================================================================================== #
col_annotation_CCYTO<- pav_mat_CCYTO %>%
  dplyr::select(Orthogroup, class) %>%
  dplyr::distinct() %>%
  tibble::column_to_rownames("Orthogroup")

# OG name matching
col_annotation_CCYTO <- col_annotation_CCYTO[colnames(pav_matrix_ordered_CCYTO), , drop = FALSE]

# ======================================================================================================================================================================================== #
# Looking at heatmap 
# ======================================================================================================================================================================================== #
# Column annotation (top)
col_ha_CCYTO <- HeatmapAnnotation(
  class = col_annotation_CCYTO$class,
  col = list(class = ann_colors$class),
  show_legend = FALSE,
  show_annotation_name = FALSE)

pav_heatmap_CCYTO <- ComplexHeatmap::Heatmap(
  pav_matrix_ordered_CCYTO,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = FALSE,
  col = c("0" = "white", "1" = "black"),
  top_annotation = col_ha_CCYTO,
  right_annotation = row_ha,
  row_names_gp = grid::gpar(fontsize = 5),
  column_title = "Cytochrome P450 orthogroups",
  column_title_side = "bottom",
  
  border = TRUE,
  show_heatmap_legend = FALSE)

ComplexHeatmap::draw(pav_heatmap_CCYTO, padding = grid::unit(c(5, 5, 5, 5), "mm"))






























# ======================================================================================================================================================================================== #
# Looking at clustering of only OGs that are found in at least one HDR
# ======================================================================================================================================================================================== #
# Load in all genes and the orthogroups they contribute to
all_genes_ogs <- readr::read_tsv("../../processed_data/orthology/all_genes_class_OGs.tsv")

# Load in HDR OGs among all wild strains
hdr_ogs <- readr::read_tsv("../../processed_data/hdr_liftover/hdr_genes_OG_class.tsv") %>%
  dplyr::distinct(Orthogroup) %>% dplyr::pull()


# Filtering for GPCR-containing orthogroups
pav_mat_HDRs <- pav_mat %>% dplyr::filter(Orthogroup %in% hdr_ogs) 

pav_matrix_for_gheatmap_HDRs <- pav_mat_HDRs %>%
  dplyr::select(-class) %>%
  tidyr::pivot_longer(cols = all_of(strainCol_c2_u), names_to = "strain", values_to = "presence") %>%
  tidyr::pivot_wider(names_from = Orthogroup, values_from = presence) %>%
  tibble::column_to_rownames("strain") %>%
  as.matrix()

# Reorder columns by gene class (core -> accessory -> private)
og_order <- pav_mat_HDRs %>% dplyr::pull(Orthogroup) %>% as.character()
pav_matrix_for_gheatmap_HDRs <- pav_matrix_for_gheatmap_HDRs[, og_order]

# ======================================================================================================================================================================================== #
# Order strains according to BUSCO-inferred phylogency
# ======================================================================================================================================================================================== #
# Reorder matrix rows to match tree order
pav_matrix_ordered_HDRs <- pav_matrix_for_gheatmap_HDRs[strain_order_by_y, ]

# Verify the strain ordering is correct
identical(rownames(pav_matrix_ordered_HDRs), strain_order_by_y)  # Should be TRUE

# ======================================================================================================================================================================================== #
# Create annotation for columns (orthogroups)
# ======================================================================================================================================================================================== #
col_annotation_HDRss<- pav_mat_HDRs %>%
  dplyr::select(Orthogroup, class) %>%
  dplyr::distinct() %>%
  tibble::column_to_rownames("Orthogroup")

# OG name matching
col_annotation_HDRs <- col_annotation_HDRs[colnames(pav_matrix_ordered_HDRs), , drop = FALSE]

# ======================================================================================================================================================================================== #
# Looking at heatmap of all orthogroups
# ======================================================================================================================================================================================== #
# Column annotation (top)
col_ha_HDRss <- HeatmapAnnotation(
  class = col_annotation_HDRs$class,
  col = list(class = ann_colors$class),
  show_legend = FALSE,
  show_annotation_name = FALSE)

pav_heatmap_HDRs <- ComplexHeatmap::Heatmap(
  pav_matrix_ordered_HDRs,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = FALSE,
  col = c("0" = "white", "1" = "black"),
  top_annotation = col_ha_HDRss,
  right_annotation = row_ha,
  row_names_gp = grid::gpar(fontsize = 5),
  column_title = "HDRs orthogroups",
  column_title_side = "bottom",
  
  border = TRUE,
  show_heatmap_legend = FALSE)

ComplexHeatmap::draw(pav_heatmap_HDRs, padding = grid::unit(c(5, 5, 5, 5), "mm"))



