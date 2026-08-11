library(readr)
library(dplyr)
library(ggplot2)
library(stringr)
library(cowplot)
library(ggtree)

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

pav_long <- pav_mat %>%
  tidyr::pivot_longer(
    cols = all_of(strainCol_c2_u),
    names_to = "strain",
    values_to = "presence") %>%
  dplyr::mutate(
    class_presence = case_when(
      presence == 1 ~ paste0(class, "_present"),
      presence == 0 ~ paste0(class, "_absent"))) %>%
  dplyr::group_by(strain,class) %>%
  dplyr::mutate(class_count=sum(presence)) %>%
  dplyr::ungroup() 

strain_order_acc <- pav_long %>%
  dplyr::filter(class=="accessory") %>%
  dplyr::arrange(class_count) %>%
  dplyr::distinct(strain,.keep_all = T) %>%
  dplyr::pull(strain) 

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

# Verify strain matching
setdiff(busco_tree$tip.label, rownames(pav_matrix_for_gheatmap)) 
setdiff(rownames(pav_matrix_for_gheatmap), busco_tree$tip.label)

# ======================================================================================================================================================================================== #
# Create plot with gheatmap
# ======================================================================================================================================================================================== #

# Base tree
p_tree <- ggtree(busco_tree_scaled) + 
  geom_tiplab(size = 1, align = TRUE, linesize = 0.2, color = 'black')

# Fix 1: Correct color mapping (low = absent = white, high = present = black)
heatmap_all <- gheatmap(
  p_tree, 
  pav_matrix_for_gheatmap,
  offset = 0.15,              
  width = 3,                 
  colnames = FALSE,
  color = NA,
  low = "white",             # 0 = absent = white
  high = "black"             # 1 = present = black
) +
  scale_fill_gradient(
    low = "white", 
    high = "black", 
    guide = "none"
  ) +
  theme(legend.position = "none",
        panel.border = element_rect(color = 'black', fill = NA))
# heatmap_all


# Colored bars above OGs indicating gene set classification
og_strip <- pav_mat %>%
  dplyr::select(Orthogroup, class) %>%
  dplyr::distinct() %>%
  dplyr::mutate(Orthogroup = factor(Orthogroup, levels = og_order))

p_strip <- ggplot(og_strip, aes(x = Orthogroup, y = 1, fill = class), alpha = 0.5) +
  geom_tile() +
  scale_fill_manual(values = c(core = "green4", accessory = "#DB6333", private = "magenta3"), guide = "none") +
  theme_void() +
  theme(
    plot.margin = margin(t = 5, r = 5, b = 0, l = 5))



# Create empty spacer for tree portion
p_spacer <- ggplot() + theme_void()

# Combine spacer + strip
p_strip_row <- cowplot::plot_grid(
  p_spacer,
  p_strip + theme(legend.position = "none"),
  ncol = 2,
  rel_widths = c(0.30, 0.72)  # Adjust these to match your tree/heatmap ratio
)

# Stack strip on top of gheatmap
final_heatmap_all <- cowplot::plot_grid(
  p_strip_row,
  heatmap_all,
  ncol = 1,
  rel_heights = c(0.08, 1),
  align = "h") + theme(panel.background = element_rect(fill = "white", color = NA))
# final_heatmap_all

ggsave("../../figures/supplementary/PAV_OG_matrix.png", final_heatmap_all, width = 7.5, height = 6, dpi = 600)


# 
# 
# # ---------------------------------------------------------
# # 2. Read phylogeny
# # ---------------------------------------------------------
# 
# tree_file <- "../../processed_data/genome_resources/trees/BUSCO_supermatrix.fa.contree"
# 
# busco_tree <- ape::read.tree(tree_file)
# 
# # Optional: make branch lengths more visually uniform
# # busco_tree_scaled <- ape::compute.brlen(
# #   busco_tree,
# #   method = "Grafen"
# # )
# 
# 
# # ---------------------------------------------------------
# # 3. Create strain × orthogroup matrix
# # ---------------------------------------------------------
# 
# pav_matrix <- pav_mat %>%
#   select(-class) %>%
#   tibble::column_to_rownames("Orthogroup") %>%
#   t() %>%
#   as.matrix()
# 
# 
# # Make sure tree and matrix contain exactly the same strains
# setdiff(busco_tree$tip.label, rownames(pav_matrix))
# setdiff(rownames(pav_matrix), busco_tree$tip.label)
# 
# # Reorder rows to exactly match the phylogeny
# pav_matrix <- pav_matrix[
#   busco_tree$tip.label,
#   ,
#   drop = FALSE
# ]
# 
# 
# # ---------------------------------------------------------
# # 4. Orthogroup classification annotation
# # ---------------------------------------------------------
# # og_annotation <- pav_mat %>%
# #   dplyr::select(Orthogroup, class) %>%
# #   dplyr::tibble::column_to_rownames("Orthogroup")
# # 
# # og_annotation <- og_annotation[
# #   colnames(pav_matrix),
# #   ,
# #   drop = FALSE
# # ]
# # 
# # colnames(og_annotation) <- "Class"
# 
# 
# # ---------------------------------------------------------
# # 5. Tree + PAV heatmap
# # ---------------------------------------------------------
# 
# p_PAV <- ggtree(busco_tree)
# 
# p_PAV$data <- ggtree::fortify(
#   busco_tree
# )
# 
# heatmap_all <- gheatmap(
#   p_PAV,
#   pav_matrix)
#   # offset = 0.02,
#   # width = 0.8,
#   # colnames = FALSE,
#   # color = "black")
#   # scale_fill_manual(values = c("0" = "white", "1" = "black"), guide = "none")
# heatmap_all
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# # ======================================================================================================================================================================================== #
# # Looking at GPCR clustering
# # ======================================================================================================================================================================================== #
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 