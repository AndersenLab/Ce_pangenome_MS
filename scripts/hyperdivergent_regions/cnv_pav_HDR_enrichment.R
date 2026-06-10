library(ggplot2)
library(dplyr)
library(data.table)
library(readr)
library(rstatix)
library(ggpubr)
  
# Calculating stats for every wild strain and appending to a master dataframe
OG_enrichment <- function(ws_hdr_ogs, strains, og_matrix_relationships, single_copy_ogs) {
  
  cnv_pav_results_df = as.data.frame(matrix(ncol = 11, nrow = 140))
  names(cnv_pav_results_df) = c("strain", "CNV_inHDR", "PAV_inHDR", "one_to_one_inHDR", "single_copy_OGs_inHDRs", "CNV_nonHDR", "PAV_nonHDR", "one_to_one_nonHDR", "single_copy_OGs_nonHDRs", "HDR_OG_count", "nonHDR_OG_count")
  
  # To prevent error-ing out because of no OG type being present
  extract_count <- function(df, relation_type) {
    val <- df %>% dplyr::filter(relation == relation_type) %>% dplyr::pull(count)
    ifelse(length(val) == 0, 0, val)
  }
  
  
  for (i in 1:length(strains)) {
    soi <- strains[i]
    print(paste0("On strain: ", soi, ". ", i, "/140."))
    
    cnv_pav_results_df[i,1] = soi
    
    # HDRs stats
    hdr_ogs <- ws_hdr_ogs %>% dplyr::filter(strain == soi) %>% dplyr::distinct(Orthogroup) %>% dplyr::pull()
    cnv_pav_results_df[i,10] = length(hdr_ogs)
  
    hdr_og_relations <- og_matrix_relationships %>% dplyr::select(Orthogroup, N2, dplyr::all_of(soi)) %>% dplyr::filter(Orthogroup %in% hdr_ogs) %>% filter(!(is.na(N2) & is.na(.data[[soi]]))) %>%
      dplyr::mutate(
        relation = dplyr::case_when(
          is.na(N2) | is.na(.data[[soi]]) ~ "PAV",
          N2 == .data[[soi]]              ~ "one_to_one",
          TRUE                            ~ "CNV" )) %>%
      dplyr::count(relation, name = "count")
    
    sg_ogs <- og_matrix_relationships %>% dplyr::filter(Orthogroup %in% hdr_ogs) %>% dplyr::filter(Orthogroup %in% sc_ogs) %>% dplyr::distinct(Orthogroup) %>% dplyr::pull()
  
    cnv_hdr <- hdr_og_relations %>% dplyr::filter(relation == "CNV") %>% dplyr::pull(count)
    pav_hdr <- hdr_og_relations %>% dplyr::filter(relation == "PAV") %>% dplyr::pull(count)
    one_to_one <- hdr_og_relations %>% dplyr::filter(relation == "one_to_one") %>% dplyr::pull(count)
    
    cnv_pav_results_df[i,2] = extract_count(hdr_og_relations, "CNV")
    cnv_pav_results_df[i,3] = extract_count(hdr_og_relations, "PAV")
    cnv_pav_results_df[i,4] = extract_count(hdr_og_relations, "one_to_one")
    cnv_pav_results_df[i,5] = length(sg_ogs)
    
    
    # outside of HDRs stats
    nonhdr_og_relations <- og_matrix_relationships %>% dplyr::select(Orthogroup, N2, dplyr::all_of(soi)) %>% dplyr::filter(!Orthogroup %in% hdr_ogs) %>% filter(!(is.na(N2) & is.na(.data[[soi]]))) %>%
      dplyr::mutate(
        relation = dplyr::case_when(
          is.na(N2) | is.na(.data[[soi]]) ~ "PAV",
          N2 == .data[[soi]]              ~ "one_to_one",
          TRUE                            ~ "CNV" )) %>%
      dplyr::count(relation, name = "count")
    
    nonhdr_ogs <- og_matrix_relationships %>% dplyr::select(Orthogroup, N2, dplyr::all_of(soi)) %>% dplyr::filter(!Orthogroup %in% hdr_ogs) %>% filter(!(is.na(N2) & is.na(.data[[soi]])))
    sg_ogs_nonHDR <- og_matrix_relationships %>% dplyr::filter(!Orthogroup %in% hdr_ogs) %>% dplyr::filter(Orthogroup %in% single_copy_ogs) %>% dplyr::distinct(Orthogroup) %>% dplyr::pull()
    
    
    cnv_pav_results_df[i,11] = nrow(nonhdr_ogs)
    
    cnv_hdr <- nonhdr_og_relations %>% dplyr::filter(relation == "CNV") %>% dplyr::pull(count)
    pav_hdr <- nonhdr_og_relations %>% dplyr::filter(relation == "PAV") %>% dplyr::pull(count)
    one_to_one <- nonhdr_og_relations %>% dplyr::filter(relation == "one_to_one") %>% dplyr::pull(count)
    
    cnv_pav_results_df[i,6] = extract_count(nonhdr_og_relations, "CNV")
    cnv_pav_results_df[i,7] = extract_count(nonhdr_og_relations, "PAV")
    cnv_pav_results_df[i,8] = extract_count(nonhdr_og_relations, "one_to_one")
    cnv_pav_results_df[i,9] = length(sg_ogs_nonHDR)
  }
  return(cnv_pav_results_df)
}

# Preparing data:
ws_hdr_ogs <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/ws_HDR_liftover/hdr_genes_OG_class.tsv")
strains <- ws_hdr_ogs %>% dplyr::distinct(strain) %>% dplyr::pull()
all_relations <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/ws_HDR_liftover/OG_relations_matrix_count.tsv") %>%
  dplyr::rename_with(~ sub("_count$", "", .x))
sc_ogs <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/orthology/elegans/orthofinder/64_core/OrthoFinder/Results_Dec07/Orthogroups/Orthogroups_SingleCopyOrthologues.txt", col_names = "OGs") %>%
  dplyr::pull(OGs)

# Running orthogroup enrichment
og_enrich_results <- OG_enrichment(ws_hdr_ogs, strains, all_relations, sc_ogs)

### Diagnostic plots
# # non-normalized results plotted
# plot_df <- og_enrich_results %>%
#   tidyr::pivot_longer(
#     cols = -strain,
#     names_to = "metric",
#     values_to = "value") %>%
#   dplyr::mutate(region = case_when(
#       stringr::str_detect(metric, "nonHDR") ~ "non-HDR",
#       stringr::str_detect(metric, "HDR") ~ "HDR",
#       TRUE ~ NA)) %>%
#   dplyr::filter(!is.na(region)) %>%
#   dplyr::mutate(region = factor(region, levels = c("HDR", "non-HDR")),
#                 metric = ifelse(metric == "HDR_OG_count", "OG_count", 
#                                 ifelse(metric == "nonHDR_OG_count", "OG_count", metric)),
#                 stat = metric %>%
#                   stringr::str_remove("_inHDRs?$") %>%
#                   stringr::str_remove("_nonHDRs?$"),
#                 stat = ifelse(stat == "one_to_one","one-to-one", 
#                               ifelse(stat == "single_copy_OGs", "single-copy OGs",
#                                      ifelse(stat == "OG_count", "OGs", stat))),
#                 stat = factor(stat, levels = c("CNV","PAV","one-to-one","single-copy OGs","OGs")))
# 
# ggplot(plot_df, aes(x = stat, y = value, fill = region)) +
#   geom_boxplot(outlier.size = 0.6, width = 0.7, position = position_dodge(width = 0.75), outlier.shape = NA, alpha = 0.5) +
#   geom_point(position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.75), size = 1.25) +
#   scale_fill_manual(values = c("HDR" = "red", "non-HDR" = "blue")) +
#   theme_bw() +
#   theme(
#     axis.title.x = element_blank(),
#     axis.text.x = element_text(size = 22, color = 'black', face = 'bold'),
#     panel.grid.major.x = element_blank(),
#     legend.title = element_blank(),
#     legend.text = element_text(size = 20, color = 'black'),
#     axis.text.y = element_text(size = 18, color = 'black'),
#     axis.title.y = element_text(size = 22, color = 'black', face = 'bold')
#   ) +
#   scale_y_log10() +
#   labs(y = "Orthogroup count", fill = "Region")
# 
# 
# # normalizing for # of orthogroups in HDRs versus not in HDRs
# plot_df_norm <- og_enrich_results %>%
#   dplyr::rename(inHDR_OG_count = HDR_OG_count) %>%
#   dplyr::mutate(
#     dplyr::across(dplyr::matches("inHDR"), ~ .x / inHDR_OG_count),
#     dplyr::across(dplyr::matches("nonHDR"), ~ .x / nonHDR_OG_count)) %>%
#   tidyr::pivot_longer(
#     cols = -strain,
#     names_to = "metric",
#     values_to = "value") %>%
#   dplyr::mutate(region = case_when(
#     stringr::str_detect(metric, "nonHDR") ~ "non-HDR",
#     stringr::str_detect(metric, "HDR") ~ "HDR",
#     TRUE ~ NA)) %>%
#   dplyr::filter(!is.na(region)) %>%
#   dplyr::mutate(region = factor(region, levels = c("HDR", "non-HDR")),
#                 metric = ifelse(metric == "inHDR_OG_count", "OG_count", 
#                                 ifelse(metric == "nonHDR_OG_count", "OG_count", metric)),
#                 stat = metric %>%
#                   stringr::str_remove("_inHDRs?$") %>%
#                   stringr::str_remove("_nonHDRs?$"),
#                 stat = ifelse(stat == "one_to_one","n-to-n", 
#                               ifelse(stat == "single_copy_OGs", "single-copy OGs",
#                                      ifelse(stat == "OG_count", "OGs", stat))),
#                 stat = factor(stat, levels = c("CNV","PAV","n-to-n","single-copy OGs","OGs"))) %>%
#   dplyr::filter(stat != "OGs")
# 
# 
# ggplot(plot_df_norm, aes(x = stat, y = value, fill = region)) +
#   geom_boxplot(outlier.size = 0.6, width = 0.7, position = position_dodge(width = 0.75), outlier.shape = NA, alpha = 0.5) +
#   geom_point(position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.75), size = 1.25) +
#   scale_fill_manual(values = c("HDR" = "red", "non-HDR" = "blue")) +
#   theme_bw() +
#   theme(
#     axis.title.x = element_blank(),
#     axis.text.x = element_text(size = 26, color = 'black', face = 'bold'),
#     panel.grid.major.x = element_blank(),
#     legend.box.background = element_rect(color = "black", size = 1),
#     legend.position  = 'inside',
#     legend.position.inside = c(0.945,0.945),
#     legend.title = element_blank(),
#     legend.key.size = unit(1.5, "cm"),
#     legend.text = element_text(size = 24, color = 'black'),
#     axis.text.y = element_text(size = 18, color = 'black'),
#     axis.title.y = element_text(size = 26, color = 'black', face = 'bold')
#   ) +
#   scale_y_log10() +
#   labs(y = "Orthogroup count", fill = "Region")


#### Calculating stats and adding to plot: 
# Paired Wilcoxon signed-rank test
plot_df_norm <- og_enrich_results %>%
  dplyr::rename(inHDR_OG_count = HDR_OG_count) %>%
  dplyr::mutate(
    dplyr::across(dplyr::matches("inHDR"), ~ .x / inHDR_OG_count),
    dplyr::across(dplyr::matches("nonHDR"), ~ .x / nonHDR_OG_count)) %>%
  tidyr::pivot_longer(
    cols = -strain,
    names_to = "metric",
    values_to = "value") %>%
  dplyr::mutate(region = case_when(
    stringr::str_detect(metric, "nonHDR") ~ "non-HDR",
    stringr::str_detect(metric, "HDR") ~ "HDR",
    TRUE ~ NA)) %>%
  dplyr::filter(!is.na(region)) %>%
  dplyr::mutate(region = factor(region, levels = c("HDR", "non-HDR")),
                metric = ifelse(metric == "inHDR_OG_count", "OG_count",
                                ifelse(metric == "nonHDR_OG_count", "OG_count", metric)),
                stat = metric %>%
                  stringr::str_remove("_inHDRs?$") %>%
                  stringr::str_remove("_nonHDRs?$"),
                stat = ifelse(stat == "one_to_one","n-to-n",
                              ifelse(stat == "single_copy_OGs", "single-copy OGs",
                                     ifelse(stat == "OG_count", "OGs", stat))),
                stat = factor(stat, levels = c("CNV","PAV","n-to-n","single-copy OGs","OGs"))) %>%
  dplyr::filter(stat != "OGs") %>%
  dplyr::mutate(value = ifelse(value == 'NaN', 0, value)) # to account for potential 0 dividided by 0 instances

wilcox_results <- plot_df_norm  %>%
  dplyr::group_by(stat) %>%
  wilcox_test(value ~ region, paired = TRUE) %>%
  adjust_pvalue(method = "BH") %>% # Benjamini-Hochberg correction
  add_significance()

y_pos <- plot_df_norm %>%
  dplyr::group_by(stat) %>%
  dplyr::summarise(y.position = max(value, na.rm = TRUE) * 1.03, .groups = "drop")

wilcox_results <- wilcox_results %>%
  dplyr::left_join(y_pos, by = "stat")

all_plt <- ggplot(plot_df_norm, aes(x = stat, y = value, fill = region)) +
  geom_boxplot(outlier.size = 0.6, width = 0.7, position = position_dodge(width = 0.75), outlier.shape = NA, alpha = 0.5) +
  geom_point(position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.75), size = 1.25) +
  scale_fill_manual(values = c("HDR" = "red", "non-HDR" = "blue")) +
  stat_pvalue_manual(wilcox_results, label = "p.adj.signif", x = "stat", y.position = "y.position", inherit.aes = FALSE, color = 'black', size = 8) +
  scale_y_log10() +
  labs(y = "Proportion of orthogroups", fill = "Region", title = "ALL ORTHOGROUPS") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 26, color = 'black'),
    panel.grid.major.x = element_blank(),
    legend.box.background = element_rect(color = "black", linewidth = 1),
    legend.position  = 'inside',
    panel.border = element_rect(color = 'black', fill =NA),
    legend.position.inside = c(0.95, 0.945),
    legend.title = element_blank(),
    legend.key.size = unit(1.5, "cm"),
    plot.title = element_text(size = 20, color = 'black', face = 'bold', hjust = 0.5),
    legend.text = element_text(size = 24, color = 'black'),
    axis.text.y = element_text(size = 18, color = 'black'),
    axis.title.y = element_text(size = 26, color = 'black', face = 'bold')
  )  +
  # scale_x_discrete(labels = c(
  #   "CNV" = "CNV",
  #   "PAV" = "PAV",
  #   "n-to-n" = expression(italic(n) * "-to-" * italic(n)),
  #   "single-copy OGs" = "single-copy OGs"
  # ))
  scale_x_discrete(labels = c(
  "CNV" = expression(bold("CNV")),
  "PAV" = expression(bold("PAV")),
  "n-to-n" = expression(bold(bolditalic(n) * "-to-" * bolditalic(n))),
  "single-copy OGs" = expression(bold("single-copy OGs"))
  ))
all_plt



all_plt_noscOG <- ggplot(plot_df_norm %>% dplyr::filter(stat != "single-copy OGs"), aes(x = stat, y = value, fill = region)) +
  geom_boxplot(outlier.size = 0.6, width = 0.7, position = position_dodge(width = 0.75), outlier.shape = NA, alpha = 0.5) +
  geom_point(position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.75), size = 1.25) +
  scale_fill_manual(values = c("HDR" = "red", "non-HDR" = "blue")) +
  stat_pvalue_manual(wilcox_results %>% dplyr::filter(stat != "single-copy OGs"), label = "p.adj.signif", x = "stat", y.position = "y.position", inherit.aes = FALSE, color = 'black', size = 8) +
  # scale_y_log10() +
  labs(y = "Proportion of orthogroups", fill = "Region") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 28, color = 'black'),
    # panel.grid.major.x = element_blank(),
    panel.grid = element_blank(),
    legend.box.background = element_rect(color = "black", linewidth = 1),
    legend.position  = 'inside',
    panel.border = element_rect(color = 'black', fill =NA),
    legend.position.inside = c(0.9, 0.1),
    legend.title = element_blank(),
    legend.key.size = unit(1.5, "cm"),
    # plot.title = element_text(size = 20, color = 'black', face = 'bold', hjust = 0.5),
    legend.text = element_text(size = 24, color = 'black'),
    axis.text.y = element_text(size = 18, color = 'black'),
    axis.title.y = element_text(size = 32, color = 'black', face = 'bold')
  )  +
  # scale_x_discrete(labels = c(
  #   "CNV" = "CNV",
  #   "PAV" = "PAV",
  #   "n-to-n" = expression(italic(n) * "-to-" * italic(n)),
  #   "single-copy OGs" = "single-copy OGs"
  # ))
  scale_x_discrete(labels = c(
    "CNV" = expression(bold("copy-number variant (CNV)")),
    "PAV" = expression(bold("presence-absence variant (PAV)")),
    "n-to-n" = expression(bold(bolditalic(n) * "-to-" * bolditalic(n))))) +
  coord_cartesian(ylim = c(-0.0001, 1.0001)) 
all_plt_noscOG

cnv_inset <- ggplot(plot_df_norm %>% dplyr::filter(stat != "single-copy OGs"), aes(x = stat, y = value, fill = region)) +
  geom_boxplot(outlier.size = 0.6, width = 0.7, position = position_dodge(width = 0.75), outlier.shape = NA, alpha = 0.5) +
  geom_point(position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.75), size = 1.25) +
  scale_fill_manual(values = c("HDR" = "red", "non-HDR" = "blue")) +
  stat_pvalue_manual(wilcox_results %>% dplyr::filter(stat != "single-copy OGs"), label = "p.adj.signif", x = "stat", y.position = "y.position", inherit.aes = FALSE, color = 'black', size = 8) +
  # scale_y_log10() +
  labs(y = "Proportion of orthogroups", fill = "Region") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 28, color = 'black'),
    # panel.grid.major.x = element_blank(),
    panel.grid = element_blank(),
    legend.box.background = element_rect(color = "black", linewidth = 1),
    legend.position  = 'inside',
    panel.border = element_rect(color = 'black', fill =NA),
    legend.position.inside = c(0.9, 0.1),
    legend.title = element_blank(),
    legend.key.size = unit(1.5, "cm"),
    # plot.title = element_text(size = 20, color = 'black', face = 'bold', hjust = 0.5),
    legend.text = element_text(size = 24, color = 'black'),
    axis.text.y = element_text(size = 18, color = 'black'),
    axis.title.y = element_text(size = 32, color = 'black', face = 'bold')
  )  +
  # scale_x_discrete(labels = c(
  #   "CNV" = "CNV",
  #   "PAV" = "PAV",
  #   "n-to-n" = expression(italic(n) * "-to-" * italic(n)),
  #   "single-copy OGs" = "single-copy OGs"
  # ))
  scale_x_discrete(labels = c(
    "CNV" = expression(bold("copy-number variant (CNV)")),
    "PAV" = expression(bold("presence-absence variant (PAV)")),
    "n-to-n" = expression(bold(bolditalic(n) * "-to-" * bolditalic(n))))) +
  coord_cartesian(ylim = c(-0.0001, 0.1)) 
cnv_inset



#####################################################################################################
# Adding on gene-class specific variation
#####################################################################################################
### GPCRs #########################################
ws_gpcrs <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/ws_HDR_liftover/140_IPR_gpcrs.tsv") %>% dplyr::mutate(gpcr = TRUE) # see line 1664 in orthogroup_vis.R

gpcr_ws_hdr_ogs <- ws_hdr_ogs %>% dplyr::left_join(ws_gpcrs, by = c("strain","gene")) %>% dplyr::filter(gpcr == TRUE)

all_ws_genes_class_og <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/ws_HDR_liftover/all_genes_class_OGs.tsv")

all_ws_gpcr_OGs <- all_ws_genes_class_og %>% dplyr::left_join(ws_gpcrs, by = c("strain","gene")) %>% dplyr::filter(gpcr == TRUE) %>% dplyr::distinct(Orthogroup) %>% dplyr::pull()

gpcr_all_relations <- all_relations %>% dplyr::filter(Orthogroup %in% all_ws_gpcr_OGs)


# Running orthogroup enrichment
og_enrich_results_GPCRs <- OG_enrichment(gpcr_ws_hdr_ogs, strains, gpcr_all_relations, sc_ogs)


#### Calculating stats and adding to plot: 
# Paired Wilcoxon signed-rank test
plot_df_norm_gpcrs <- og_enrich_results_GPCRs %>%
  dplyr::rename(inHDR_OG_count = HDR_OG_count) %>%
  dplyr::mutate(
    dplyr::across(dplyr::matches("inHDR"), ~ .x / inHDR_OG_count),
    dplyr::across(dplyr::matches("nonHDR"), ~ .x / nonHDR_OG_count)) %>%
  tidyr::pivot_longer(
    cols = -strain,
    names_to = "metric",
    values_to = "value") %>%
  dplyr::mutate(region = case_when(
    stringr::str_detect(metric, "nonHDR") ~ "non-HDR",
    stringr::str_detect(metric, "HDR") ~ "HDR",
    TRUE ~ NA)) %>%
  dplyr::filter(!is.na(region)) %>%
  dplyr::mutate(region = factor(region, levels = c("HDR", "non-HDR")),
                metric = ifelse(metric == "inHDR_OG_count", "OG_count",
                                ifelse(metric == "nonHDR_OG_count", "OG_count", metric)),
                stat = metric %>%
                  stringr::str_remove("_inHDRs?$") %>%
                  stringr::str_remove("_nonHDRs?$"),
                stat = ifelse(stat == "one_to_one","n-to-n",
                              ifelse(stat == "single_copy_OGs", "single-copy OGs",
                                     ifelse(stat == "OG_count", "OGs", stat))),
                stat = factor(stat, levels = c("CNV","PAV","n-to-n","single-copy OGs","OGs"))) %>%
  dplyr::filter(stat != "OGs") %>%
  dplyr::mutate(value = ifelse(value == 'NaN', 0, value))

wilcox_results_gpcrs <- plot_df_norm_gpcrs  %>%
  dplyr::group_by(stat) %>%
  wilcox_test(value ~ region, paired = TRUE) %>%
  adjust_pvalue(method = "BH") %>% # Benjamini-Hochberg correction
  add_significance()

y_pos_gpcrs <- plot_df_norm_gpcrs %>%
  dplyr::group_by(stat) %>%
  dplyr::summarise(y.position = max(value, na.rm = TRUE) * 1.03, .groups = "drop")

wilcox_results_gp <- wilcox_results_gpcrs %>%
  dplyr::left_join(y_pos_gpcrs, by = "stat")

gpcrs_plt <- ggplot(plot_df_norm_gpcrs, aes(x = stat, y = value, fill = region)) +
  geom_boxplot(outlier.size = 0.6, width = 0.7, position = position_dodge(width = 0.75), outlier.shape = NA, alpha = 0.5) +
  geom_point(position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.75), size = 1.25) +
  scale_fill_manual(values = c("HDR" = "red", "non-HDR" = "blue")) +
  stat_pvalue_manual(wilcox_results_gp, label = "p.adj.signif", x = "stat", y.position = "y.position", inherit.aes = FALSE, color = 'black', size = 8) +
  scale_y_log10() +
  labs(y = "Proportion of orthogroups", fill = "Region", title = "GPCRs") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 26, color = 'black'),
    panel.grid.major.x = element_blank(),
    legend.box.background = element_rect(color = "black", size = 1),
    legend.position  = 'inside',
    panel.border = element_rect(color = 'black', fill = NA),
    legend.position.inside = c(0.945, 0.945),
    legend.title = element_blank(),
    plot.title = element_text(size = 20, color = 'black', face = 'bold', hjust = 0.5),
    legend.key.size = unit(1.5, "cm"),
    legend.text = element_text(size = 24, color = 'black'),
    axis.text.y = element_text(size = 18, color = 'black'),
    axis.title.y = element_text(size = 26, color = 'black', face = 'bold')
  )  +
  scale_x_discrete(labels = c(
    "CNV" = expression(bold("CNV")),
    "PAV" = expression(bold("PAV")),
    "n-to-n" = expression(bold(bolditalic(n) * "-to-" * bolditalic(n))),
    "single-copy OGs" = expression(bold("single-copy OGs"))
  )) 
gpcrs_plt









### F-box genes ############################################
ws_fbox <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/ws_HDR_liftover/140_IPR_fBox.tsv") %>% dplyr::mutate(fbox = TRUE) # see line 1664 in orthogroup_vis.R

fbox_ws_hdr_ogs <- ws_hdr_ogs %>% dplyr::left_join(ws_fbox, by = c("strain","gene")) %>% dplyr::filter(fbox == TRUE)

all_ws_fbox_OGs <- all_ws_genes_class_og %>% dplyr::left_join(ws_fbox, by = c("strain","gene")) %>% dplyr::filter(fbox == TRUE) %>% dplyr::distinct(Orthogroup) %>% dplyr::pull()

fbox_all_relations <- all_relations %>% dplyr::filter(Orthogroup %in% all_ws_fbox_OGs)


# Running orthogroup enrichment
og_enrich_results_FBOX <- OG_enrichment(fbox_ws_hdr_ogs, strains, fbox_all_relations, sc_ogs)


#### Calculating stats and adding to plot: 
# Paired Wilcoxon signed-rank test
plot_df_norm_fbox <- og_enrich_results_FBOX %>%
  dplyr::rename(inHDR_OG_count = HDR_OG_count) %>%
  dplyr::mutate(
    dplyr::across(dplyr::matches("inHDR"), ~ .x / inHDR_OG_count),
    dplyr::across(dplyr::matches("nonHDR"), ~ .x / nonHDR_OG_count)) %>%
  tidyr::pivot_longer(
    cols = -strain,
    names_to = "metric",
    values_to = "value") %>%
  dplyr::mutate(region = case_when(
    stringr::str_detect(metric, "nonHDR") ~ "non-HDR",
    stringr::str_detect(metric, "HDR") ~ "HDR",
    TRUE ~ NA)) %>%
  dplyr::filter(!is.na(region)) %>%
  dplyr::mutate(region = factor(region, levels = c("HDR", "non-HDR")),
                metric = ifelse(metric == "inHDR_OG_count", "OG_count",
                                ifelse(metric == "nonHDR_OG_count", "OG_count", metric)),
                stat = metric %>%
                  stringr::str_remove("_inHDRs?$") %>%
                  stringr::str_remove("_nonHDRs?$"),
                stat = ifelse(stat == "one_to_one","n-to-n",
                              ifelse(stat == "single_copy_OGs", "single-copy OGs",
                                     ifelse(stat == "OG_count", "OGs", stat))),
                stat = factor(stat, levels = c("CNV","PAV","n-to-n","single-copy OGs","OGs"))) %>%
  dplyr::filter(stat != "OGs") %>%
  dplyr::mutate(value = ifelse(value == 'NaN', 0, value))

wilcox_results_fbox <- plot_df_norm_fbox  %>%
  dplyr::group_by(stat) %>%
  wilcox_test(value ~ region, paired = TRUE) %>%
  adjust_pvalue(method = "BH") %>% # Benjamini-Hochberg correction
  add_significance()

y_pos_fbox <- plot_df_norm_fbox %>%
  dplyr::group_by(stat) %>%
  dplyr::summarise(y.position = max(value, na.rm = TRUE) * 1.03, .groups = "drop")

wilcox_results_final <- wilcox_results_fbox %>%
  dplyr::left_join(y_pos_fbox, by = "stat")

fbox_plt <- ggplot(plot_df_norm_fbox, aes(x = stat, y = value, fill = region)) +
  geom_boxplot(outlier.size = 0.6, width = 0.7, position = position_dodge(width = 0.75), outlier.shape = NA, alpha = 0.5) +
  geom_point(position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.75), size = 1.25) +
  scale_fill_manual(values = c("HDR" = "red", "non-HDR" = "blue")) +
  stat_pvalue_manual(wilcox_results_final, label = "p.adj.signif", x = "stat", y.position = "y.position", inherit.aes = FALSE, color = 'black', size = 8) +
  scale_y_log10() +
  labs(y = "Proportional orthogroup count", fill = "Region", title = "F-box") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 26, color = 'black'),
    panel.grid.major.x = element_blank(),
    legend.box.background = element_rect(color = "black", size = 1),
    legend.position  = 'inside',
    panel.border = element_rect(color = 'black', fill = NA),
    legend.position.inside = c(0.945, 0.945),
    legend.title = element_blank(),
    plot.title = element_text(size = 20, color = 'black', face = 'bold', hjust = 0.5),
    legend.key.size = unit(1.5, "cm"),
    legend.text = element_text(size = 24, color = 'black'),
    axis.text.y = element_text(size = 18, color = 'black'),
    axis.title.y = element_text(size = 26, color = 'black', face = 'bold')
  )  +
  scale_x_discrete(labels = c(
    "CNV" = expression(bold("CNV")),
    "PAV" = expression(bold("PAV")),
    "n-to-n" = expression(bold(bolditalic(n) * "-to-" * bolditalic(n))),
    "single-copy OGs" = expression(bold("single-copy OGs"))
  )) 
fbox_plt













### C-type lectins ############################################
ws_lectin <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/ws_HDR_liftover/140_IPR_CtypeLectins.tsv") %>% dplyr::mutate(lectin = TRUE) # see line 1664 in orthogroup_vis.R

lectin_ws_hdr_ogs <- ws_hdr_ogs %>% dplyr::left_join(ws_lectin, by = c("strain","gene")) %>% dplyr::filter(lectin == TRUE)

all_ws_lectin_OGs <- all_ws_genes_class_og %>% dplyr::left_join(ws_lectin, by = c("strain","gene")) %>% dplyr::filter(lectin == TRUE) %>% dplyr::distinct(Orthogroup) %>% dplyr::pull()

lectin_all_relations <- all_relations %>% dplyr::filter(Orthogroup %in% all_ws_lectin_OGs)


# Running orthogroup enrichment
og_enrich_results_LECTIN <- OG_enrichment(lectin_ws_hdr_ogs, strains, lectin_all_relations, sc_ogs)


#### Calculating stats and adding to plot: 
# Paired Wilcoxon signed-rank test
plot_df_norm_lectin <- og_enrich_results_LECTIN %>%
  dplyr::rename(inHDR_OG_count = HDR_OG_count) %>%
  dplyr::mutate(
    dplyr::across(dplyr::matches("inHDR"), ~ .x / inHDR_OG_count),
    dplyr::across(dplyr::matches("nonHDR"), ~ .x / nonHDR_OG_count)) %>%
  tidyr::pivot_longer(
    cols = -strain,
    names_to = "metric",
    values_to = "value") %>%
  dplyr::mutate(region = case_when(
    stringr::str_detect(metric, "nonHDR") ~ "non-HDR",
    stringr::str_detect(metric, "HDR") ~ "HDR",
    TRUE ~ NA)) %>%
  dplyr::filter(!is.na(region)) %>%
  dplyr::mutate(region = factor(region, levels = c("HDR", "non-HDR")),
                metric = ifelse(metric == "inHDR_OG_count", "OG_count",
                                ifelse(metric == "nonHDR_OG_count", "OG_count", metric)),
                stat = metric %>%
                  stringr::str_remove("_inHDRs?$") %>%
                  stringr::str_remove("_nonHDRs?$"),
                stat = ifelse(stat == "one_to_one","n-to-n",
                              ifelse(stat == "single_copy_OGs", "single-copy OGs",
                                     ifelse(stat == "OG_count", "OGs", stat))),
                stat = factor(stat, levels = c("CNV","PAV","n-to-n","single-copy OGs","OGs"))) %>%
  dplyr::filter(stat != "OGs") %>%
  dplyr::mutate(value = ifelse(value == 'NaN', 0, value))

wilcox_results_lectin <- plot_df_norm_lectin  %>%
  dplyr::group_by(stat) %>%
  wilcox_test(value ~ region, paired = TRUE) %>%
  adjust_pvalue(method = "BH") %>% # Benjamini-Hochberg correction
  add_significance()

y_pos_lectin <- plot_df_norm_lectin %>%
  dplyr::group_by(stat) %>%
  dplyr::summarise(y.position = max(value, na.rm = TRUE) * 1.03, .groups = "drop")

wilcox_results_ct <- wilcox_results_lectin %>%
  dplyr::left_join(y_pos_lectin, by = "stat")

lectins_plt <- ggplot(plot_df_norm_lectin, aes(x = stat, y = value, fill = region)) +
  geom_boxplot(outlier.size = 0.6, width = 0.7, position = position_dodge(width = 0.75), outlier.shape = NA, alpha = 0.5) +
  geom_point(position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.75), size = 1.25) +
  scale_fill_manual(values = c("HDR" = "red", "non-HDR" = "blue")) +
  stat_pvalue_manual(wilcox_results_ct, label = "p.adj.signif", x = "stat", y.position = "y.position", inherit.aes = FALSE, color = 'black', size = 8) +
  scale_y_log10() +
  labs(y = "Proportional orthogroup count", fill = "Region", title = "C-type lectins") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 26, color = 'black'),
    panel.grid.major.x = element_blank(),
    legend.box.background = element_rect(color = "black", size = 1),
    legend.position  = 'inside',
    panel.border = element_rect(color = 'black', fill = NA),
    legend.position.inside = c(0.945, 0.945),
    legend.title = element_blank(),
    plot.title = element_text(size = 20, color = 'black', face = 'bold', hjust = 0.5),
    legend.key.size = unit(1.5, "cm"),
    legend.text = element_text(size = 24, color = 'black'),
    axis.text.y = element_text(size = 18, color = 'black'),
    axis.title.y = element_text(size = 26, color = 'black', face = 'bold')
  )  +
  scale_x_discrete(labels = c(
    "CNV" = expression(bold("CNV")),
    "PAV" = expression(bold("PAV")),
    "n-to-n" = expression(bold(bolditalic(n) * "-to-" * bolditalic(n))),
    "single-copy OGs" = expression(bold("single-copy OGs"))
  )) 
lectins_plt






### Cytochrome P450s ############################################
ws_cyto <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/ws_HDR_liftover/140_IPR_cytochromeP450.tsv") %>% dplyr::mutate(cyto = TRUE) # see line 1664 in orthogroup_vis.R

cyto_ws_hdr_ogs <- ws_hdr_ogs %>% dplyr::left_join(ws_cyto, by = c("strain","gene")) %>% dplyr::filter(cyto == TRUE)

all_ws_cyto_OGs <- all_ws_genes_class_og %>% dplyr::left_join(ws_cyto, by = c("strain","gene")) %>% dplyr::filter(cyto == TRUE) %>% dplyr::distinct(Orthogroup) %>% dplyr::pull()

cyto_all_relations <- all_relations %>% dplyr::filter(Orthogroup %in% all_ws_cyto_OGs)


# Running orthogroup enrichment
og_enrich_results_cyto <- OG_enrichment(cyto_ws_hdr_ogs, strains, cyto_all_relations, sc_ogs)


#### Calculating stats and adding to plot: 
# Paired Wilcoxon signed-rank test
plot_df_norm_cyto <- og_enrich_results_cyto %>%
  dplyr::rename(inHDR_OG_count = HDR_OG_count) %>%
  dplyr::mutate(
    dplyr::across(dplyr::matches("inHDR"), ~ .x / inHDR_OG_count),
    dplyr::across(dplyr::matches("nonHDR"), ~ .x / nonHDR_OG_count)) %>%
  tidyr::pivot_longer(
    cols = -strain,
    names_to = "metric",
    values_to = "value") %>%
  dplyr::mutate(region = case_when(
    stringr::str_detect(metric, "nonHDR") ~ "non-HDR",
    stringr::str_detect(metric, "HDR") ~ "HDR",
    TRUE ~ NA)) %>%
  dplyr::filter(!is.na(region)) %>%
  dplyr::mutate(region = factor(region, levels = c("HDR", "non-HDR")),
                metric = ifelse(metric == "inHDR_OG_count", "OG_count",
                                ifelse(metric == "nonHDR_OG_count", "OG_count", metric)),
                stat = metric %>%
                  stringr::str_remove("_inHDRs?$") %>%
                  stringr::str_remove("_nonHDRs?$"),
                stat = ifelse(stat == "one_to_one","n-to-n",
                              ifelse(stat == "single_copy_OGs", "single-copy OGs",
                                     ifelse(stat == "OG_count", "OGs", stat))),
                stat = factor(stat, levels = c("CNV","PAV","n-to-n","single-copy OGs","OGs"))) %>%
  dplyr::filter(stat != "OGs") %>%
  dplyr::mutate(value = ifelse(value == 'NaN', 0, value))

wilcox_results_cyto <- plot_df_norm_cyto  %>%
  dplyr::group_by(stat) %>%
  wilcox_test(value ~ region, paired = TRUE) %>%
  adjust_pvalue(method = "BH") %>% # Benjamini-Hochberg correction
  add_significance()

y_pos_cyto <- plot_df_norm_cyto %>%
  dplyr::group_by(stat) %>%
  dplyr::summarise(y.position = max(value, na.rm = TRUE) * 1.03, .groups = "drop")

wilcox_results_cyto <- wilcox_results_cyto %>%
  dplyr::left_join(y_pos_cyto, by = "stat")

cytos_plt <- ggplot(plot_df_norm_cyto, aes(x = stat, y = value, fill = region)) +
  geom_boxplot(outlier.size = 0.6, width = 0.7, position = position_dodge(width = 0.75), outlier.shape = NA, alpha = 0.5) +
  geom_point(position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.75), size = 1.25) +
  scale_fill_manual(values = c("HDR" = "red", "non-HDR" = "blue")) +
  stat_pvalue_manual(wilcox_results_ct, label = "p.adj.signif", x = "stat", y.position = "y.position", inherit.aes = FALSE, color = 'black', size = 8) +
  scale_y_log10() +
  labs(y = "Proportional orthogroup count", fill = "Region", title = "Cytochrome P450s") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 26, color = 'black'),
    # panel.grid.major.x = element_blank(),
    legend.box.background = element_rect(color = "black", size = 1),
    legend.position  = 'inside',
    panel.border = element_rect(color = 'black', fill = NA),
    legend.position.inside = c(0.945, 0.945),
    legend.title = element_blank(),
    plot.title = element_text(size = 20, color = 'black', face = 'bold', hjust = 0.5),
    legend.key.size = unit(1.5, "cm"),
    legend.text = element_text(size = 24, color = 'black'),
    axis.text.y = element_text(size = 18, color = 'black'),
    axis.title.y = element_text(size = 26, color = 'black', face = 'bold')
  )  +
  scale_x_discrete(labels = c(
    "CNV" = expression(bold("CNV")),
    "PAV" = expression(bold("PAV")),
    "n-to-n" = expression(bold(bolditalic(n) * "-to-" * bolditalic(n))),
    "single-copy OGs" = expression(bold("single-copy OGs"))
  )) 
cytos_plt










############################################################################################################################################ 
################# Plotting all gene classes and then faceting by class ###################
############################################################################################################################################ 
final_cnv_pav_geneclass <- plot_df_norm %>% dplyr::mutate(gene_class = "All genes") %>% 
  dplyr::bind_rows(plot_df_norm_cyto %>% dplyr::mutate(gene_class = "Cytochrome P450s"),
                   plot_df_norm_fbox %>% dplyr::mutate(gene_class = "F-box genes"),
                   plot_df_norm_gpcrs %>% dplyr::mutate(gene_class = "GPCRs"),
                   plot_df_norm_lectin %>% dplyr::mutate(gene_class = "C-type lectins")) %>%
  dplyr::filter(stat == "CNV" | stat == "PAV" | stat == "n-to-n") %>%
  dplyr::mutate(gene_class = factor(gene_class, levels = c("All genes", "F-box genes", "GPCRs", "C-type lectins", "Cytochrome P450s")))

concatenated_stats <- wilcox_results %>% dplyr::mutate(gene_class = "All genes") %>% 
  dplyr::bind_rows(wilcox_results_cyto %>% dplyr::mutate(gene_class = "Cytochrome P450s"),
                   wilcox_results_final %>% dplyr::mutate(gene_class = "F-box genes"),
                   wilcox_results_gp %>% dplyr::mutate(gene_class = "GPCRs"),
                   wilcox_results_ct %>% dplyr::mutate(gene_class = "C-type lectins")) %>%
  dplyr::filter(stat == "CNV" | stat == "PAV" | stat == "n-to-n") %>%
  dplyr::mutate(gene_class = factor(gene_class, levels = c("All genes", "F-box genes", "GPCRs", "C-type lectins", "Cytochrome P450s")))


all_plt_class <- ggplot(final_cnv_pav_geneclass, aes(x = stat, y = value, fill = region)) +
  geom_boxplot(outlier.size = 0.6, width = 0.7, position = position_dodge(width = 0.75), outlier.shape = NA, alpha = 0.5) +
  geom_point(position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.75), size = 1.25) +
  scale_fill_manual(values = c("HDR" = "red", "non-HDR" = "blue")) +
  stat_pvalue_manual(concatenated_stats, label = "p.adj.signif", x = "stat", y.position = "y.position", inherit.aes = FALSE, color = 'black', size = 8) +
  # scale_y_log10() +
  facet_wrap(~gene_class, nrow = 1) +
  labs(y = "Proportion of orthogroups", fill = "Region") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 26, color = 'black'),
    panel.grid = element_blank(),
    legend.box.background = element_rect(color = "black", size = 1),
    # legend.position  = 'inside',
    panel.border = element_rect(color = 'black', fill = NA),
    # legend.position.inside = c(0.945, 0.945),
    legend.title = element_blank(),
    strip.text = element_text(size = 26, color = 'black', face = 'bold'),
    plot.title = element_text(size = 20, color = 'black', face = 'bold', hjust = 0.5),
    legend.key.size = unit(1.5, "cm"),
    legend.text = element_text(size = 24, color = 'black'),
    axis.text.y = element_text(size = 18, color = 'black'),
    axis.title.y = element_text(size = 26, color = 'black', face = 'bold')
  )  +
  scale_x_discrete(labels = c(
      "CNV" = expression(bold("CNV")),
      "PAV" = expression(bold("PAV")),
      "n-to-n" = expression(bold(bolditalic(n) * "-to-" * bolditalic(n))))) 
  # coord_cartesian(ylim = c(-0.0001, 1.0001))
all_plt_class
                     


fig_plot <- ggplot(final_cnv_pav_geneclass, aes(x = stat, y = value, fill = region)) +
  geom_boxplot(outlier.size = 0.6, width = 0.7, position = position_dodge(width = 0.75), outlier.shape = NA, alpha = 0.5) +
  geom_point(position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.75), size = 0.01) +
  scale_fill_manual(values = c("HDR" = "red", "non-HDR" = "blue")) +
  stat_pvalue_manual(concatenated_stats, label = "p.adj.signif", x = "stat", y.position = "y.position", inherit.aes = FALSE, color = 'black', size = 4) +
  # scale_y_log10() +
  facet_wrap(~gene_class, nrow = 1) +
  labs(y = "Proportional orthogroup count", fill = "Region") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 10, color = 'black'),
    panel.grid = element_blank(),
    legend.box.background = element_rect(color = "black", size = 1),
    legend.position  = 'none',
    panel.border = element_rect(color = 'black', fill = NA),
    strip.text = element_text(size = 9.5, color = 'black'),
    # plot.title = element_text(size = 20, color = 'black', hjust = 0.5),
    legend.text = element_text(size = 10, color = 'black'),
    axis.text.y = element_text(size = 10, color = 'black'),
    axis.title.y = element_text(size = 10, color = 'black')
  )  +
  scale_x_discrete(labels = c(
    # "CNV" = expression(bold("CNV")),
    # "PAV" = expression(bold("PAV")),
    "n-to-n" = expression(italic("n") * "-to-" * italic("n"))))
# coord_cartesian(ylim = c(-0.0001, 1.0001))
fig_plot


# ggsave("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/plots/cnv_pav_enrichmentPlot.png", fig_plot, dpi = 600, width = 7.5, height = 4)










############################################################################################################################################ 
############################ Looking at averages among all wild strains ############################ 
############################################################################################################################################ 
all_stats <- og_enrich_results_GPCRs %>% dplyr::mutate(type = "GPCRs") %>%
  dplyr::bind_rows((og_enrich_results %>% dplyr::mutate(type = "ALL"))) %>%
  dplyr::bind_rows((og_enrich_results_FBOX %>% dplyr::mutate(type = "FBOX"))) %>%
  dplyr::bind_rows((og_enrich_results_LECTIN %>% dplyr::mutate(type = "C_type_lectins"))) %>% 
  dplyr::bind_rows((og_enrich_results_cyto %>% dplyr::mutate(type = "Cytochrome_P450s"))) %>% 
  dplyr::group_by(type) %>%
  dplyr::mutate(mean_CNV_inHDR = mean(CNV_inHDR),
                mean_PAV_inHDR = mean(PAV_inHDR),
                mean_HDR_ogCount = mean(HDR_OG_count),
                mean_nonHDR_ogCount = mean(nonHDR_OG_count)) %>%
  dplyr::ungroup()

plt_all_stats <- all_stats %>% dplyr::select(mean_CNV_inHDR, mean_PAV_inHDR, mean_HDR_ogCount, mean_nonHDR_ogCount, type) %>%
  dplyr::distinct() %>%
  dplyr::mutate(across(-last_col(), round)) %>%
  tidyr::pivot_longer(
    cols = -type,
    names_to = "metric",
    values_to = "value") %>%
  dplyr::mutate(metric = factor(metric, levels = c("mean_CNV_inHDR","mean_PAV_inHDR", "mean_HDR_ogCount","mean_nonHDR_ogCount")),
                type = factor(type, levels = c("ALL","Cytochrome_P450s","C_type_lectins","FBOX","GPCRs")))


ggplot(data = plt_all_stats) +
  geom_col(aes(x = metric, y = value, fill = type), position = "dodge") +
  scale_fill_manual(values = c("ALL" = "black", "GPCRs" = "olivedrab", "FBOX" = "firebrick", "C_type_lectins" = "steelblue", "Cytochrome_P450s" = "orange")) +
  scale_y_log10(expand = c(0,0)) +
  theme_bw() +
  theme(
    panel.border = element_rect(color = 'black', fill = NA),
    legend.title = element_blank(),
    legend.text = element_text(size = 14, color = 'black'),
    legend.box.background = element_rect(color = 'black', fill = NA),
    legend.position = 'inside',
    legend.position.inside = c(0.9,0.9),
    axis.text.x = element_text(size = 14, color = 'black', face = 'bold'),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14, color = 'black', face = 'bold'),
    axis.text.y = element_text(size = 12, color = 'black')
  ) +
  labs(y = "OG count")

# The sum of mean_HDR_ogCount for these three gene classes is 32% of the mean for ALL mean_HDR_ogCount
