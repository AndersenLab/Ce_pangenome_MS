library(tidyr)
library(dplyr)
library(ggplot2)
library(readr)
library(cowplot)

# Load in ortholog matrix with counts of paralogs
ortho_matrix <- readr::read_tsv("../../processed_data/orthology/OG_relations_matrix_count.tsv")

# Adding resolution of OG occupancy
ortho_count <- ortho_matrix %>%
  dplyr::mutate(sum = rowSums(!is.na(dplyr::across(-1)))) %>%
  dplyr::mutate(freq = (sum / 142)) 

# Determine an occupancy cutoff:
occupancy <- data.frame()
for (i in seq(0, 1, 0.01)) {
  ortho_occ <- ortho_count %>%
    dplyr::filter(freq >= i) %>%
    dplyr::summarise(occupancy_threshold = i, ortho_count = dplyr::n())

  occupancy <- dplyr::bind_rows(occupancy, ortho_occ)
}

ggplot(occupancy) +
  geom_point(aes(x = occupancy_threshold, y = ortho_count), size = 7, color = 'black') +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(color = 'black', fill = NA),
    axis.text = element_text(size = 18, color = 'black'),
    axis.title = element_text(size = 20, color = 'black')
  ) +
  labs(x = 'Occupancy threshold', y = 'Orthogroup count')
### no clear decrease at a certain occupancy... arbitrarily setting threshold at 50%


ortho_count_filt <- ortho_count %>%
  dplyr::filter(freq >= 0.5) %>%
  dplyr::mutate(
    average_paralog_count = rowMeans(dplyr::across(2:(ncol(.) - 2), ~ as.numeric(.x)), na.rm = TRUE),
    sd_paralog_count = apply(dplyr::across(2:(ncol(ortho_count) - 2), ~ as.numeric(.x)), 1, sd, na.rm = TRUE)) %>%
  dplyr::filter(sd_paralog_count != 0) %>%
  dplyr::arrange(desc(sd_paralog_count)) %>%
  dplyr::slice_head(n = 20) %>%
  dplyr::select(Orthogroup, sum, freq, average_paralog_count, sd_paralog_count)

order <- ortho_count_filt %>% dplyr::pull(Orthogroup)

top <- ggplot(ortho_count_filt %>% dplyr::mutate(Orthogroup = factor(Orthogroup, levels = order)), aes(x = Orthogroup, y = average_paralog_count)) +
  geom_errorbar(aes(ymin = average_paralog_count - sd_paralog_count, ymax = average_paralog_count + sd_paralog_count), 
                width = 0.2) +
  geom_point(aes(fill = freq), size = 2, shape = 22) +
  scale_fill_gradient(low = 'yellow', high = 'red') +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(color = 'black', fill = NA),
    axis.text.y = element_text(size = 10, color = 'black'),
    axis.text.x = element_text(size = 10, color = 'black', angle = 60, hjust = 1),
    legend.position = "inside",
    legend.position.inside = c(0.92,0.7),
    legend.text = element_text(size = 8, color= 'black'),
    legend.title = element_text(size = 8, color = 'black'),
    plot.margin = margin(b = 0, t = 2, l = 2),
    axis.title = element_text(size = 10, color = 'black')
  ) +
  labs(x = 'Orthogroup', y = 'Average paralog count', fill = "Occupancy")

most_variance_orthos <- ortho_count_filt %>% dplyr::pull(Orthogroup)

### Does average paralog count among orthogroups correlate with standard deviation of paralogs?
ortho_count_filt_all <- ortho_count %>%
  dplyr::filter(freq >= 0.5) %>%
  dplyr::mutate(
    average_paralog_count = rowMeans(dplyr::across(2:(ncol(.) - 2), ~ as.numeric(.x)), na.rm = TRUE),
    sd_paralog_count = apply(dplyr::across(2:(ncol(ortho_count) - 2), ~ as.numeric(.x)), 1, sd, na.rm = TRUE)) %>%
  dplyr::filter(sd_paralog_count != 0) %>%
  dplyr::arrange(desc(sd_paralog_count))

ggplot(ortho_count_filt_all) +
  geom_point(aes(x = average_paralog_count, y = sd_paralog_count, fill = freq), size = 4.5, shape = 21) +
  scale_fill_gradient(low = 'yellow', high = 'red') +
  theme(panel.background = element_blank(),
        panel.border = element_rect(color = 'black', fill = NA),
        legend.text = element_text(size = 18, color = 'black'),
        legend.title = element_text(size = 20, color = 'black'),
        axis.text = element_text(size = 22, color = 'black'),
        axis.title = element_text(size = 24, color = 'black')) +
  labs(x = "Average paralog count", y = "Standard deviation", fill = "Occupancy") +
  scale_x_continuous(breaks = seq(0,20,1))



### Adding IPR resolution:
genes_ogs <- readr::read_tsv("../../processed_data/orthology/all_genes_class_OGs.tsv")

ipr_annotations <- readr::read_tsv("../../tables/IPR_annotation_142strains.tsv", col_names = c("strain_transcript", "MD5_digest",	"Seq_length",	"Application",	"Signature_accession",	"Signature_description",	"Amino_acid_start_position", "Amino_acid_end_position", "Score", "status", "date", "IPR_accession", "IPR_description", "GO_ID", "pathway")) %>%
  dplyr::filter(IPR_description != "-") %>%
  dplyr::select(strain_transcript, IPR_accession, IPR_description) %>%
  dplyr::mutate(strain_transcript = sub("\\.[^.]*$", "", strain_transcript)) %>%
  dplyr::mutate(strain_transcript = sub("[a-z]+$", "", strain_transcript))

genes_ogs_ipr <- genes_ogs %>%
  dplyr::mutate(gene = gsub("transcript_","",gene)) %>%
  dplyr::mutate(strain_transcript = paste0(strain,"_",gene)) %>%
  dplyr::left_join(ipr_annotations, by = "strain_transcript") %>%
  dplyr::distinct(strain, Orthogroup, gene, IPR_description, IPR_accession)

ogs_IPR_count <- genes_ogs_ipr %>% dplyr::group_by(Orthogroup, IPR_description) %>%
  dplyr::mutate(ipr_per_og_count = n()) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(Orthogroup, IPR_description, ipr_per_og_count) %>%
  dplyr::group_by(Orthogroup) %>%
  dplyr::arrange(desc(ipr_per_og_count)) %>%
  dplyr::slice_head(n = 1) %>%
  dplyr::mutate(total_ipr_perOG = sum(ipr_per_og_count),
                prop_IPR_perOG = ipr_per_og_count / total_ipr_perOG) %>%
  dplyr::filter(Orthogroup %in% most_variance_orthos) %>%
  dplyr::mutate(Orthogroup = factor(Orthogroup, levels = most_variance_orthos)) %>%
  dplyr::arrange(Orthogroup) %>%
  dplyr::mutate(IPR_description = factor(IPR_description, levels = unique(IPR_description)))

my_colors <- c(
  "#E41A1C",
  "#377EB8",
  "#8DA0CB",
  "skyblue", 
  "plum1",   
  "hotpink", 
  "firebrick",  
  "forestgreen",  
  "#FFD92F", 
  "gray90",  
  "magenta3",
  "#CAB2D6")

bottom <- ggplot(ogs_IPR_count) + 
  geom_col(aes(x = Orthogroup, y = prop_IPR_perOG, fill = IPR_description)) +
  scale_fill_manual(values = my_colors, na.value = "black") +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(color = 'black', fill = NA),
    axis.text.x = element_text(size = 3, color = 'black', angle = 60, hjust = 1),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    plot.margin = margin(t = 0, b = 2, r =2),
    legend.key.size = unit(0.4, "cm"),
    legend.text = element_text(size = 5.5, color= 'black'),
    legend.title = element_blank(),
    axis.title.x = element_blank()) +
  guides(fill = guide_legend(nrow = 13)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(fill = "IPR description")

legend <- cowplot::get_legend(bottom)

combined <- cowplot::plot_grid(
  top + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank()), bottom + theme(legend.position = "none"),  
  align = "v",
  nrow = 2,
  rel_heights = c(1, 0.2))

# Create final plot with IPR annotations added
combined_final <- cowplot::ggdraw(combined) +
  cowplot::draw_grob(
    legend,
    x = 0.43,
    y = 0.39,
    width = 0.22,
    height = 0.65)


# ========================================================================================================================================================================================== #
# LOOKING AT ONLY CORE ORTHOGROUPS
# ========================================================================================================================================================================================== #
ortho_count_filt_core <- ortho_count %>%
  dplyr::filter(freq == 1) %>%
  dplyr::mutate(
    average_paralog_count = rowMeans(dplyr::across(2:(ncol(.) - 2), ~ as.numeric(.x)), na.rm = TRUE),
    sd_paralog_count = apply(dplyr::across(2:(ncol(ortho_count) - 2), ~ as.numeric(.x)), 1, sd, na.rm = TRUE)) %>%
  dplyr::filter(sd_paralog_count != 0) %>%
  dplyr::arrange(desc(sd_paralog_count)) %>%
  dplyr::slice_head(n = 20) %>%
  dplyr::select(Orthogroup, sum, freq, average_paralog_count, sd_paralog_count)

order <- ortho_count_filt_core %>% dplyr::pull(Orthogroup)

top_core <- ggplot(ortho_count_filt_core %>% dplyr::mutate(Orthogroup = factor(Orthogroup, levels = order)), aes(x = Orthogroup, y = average_paralog_count)) +
  geom_errorbar(aes(ymin = average_paralog_count - sd_paralog_count, ymax = average_paralog_count + sd_paralog_count), 
                width = 0.2) +
  geom_point(fill = "red", size = 2, shape = 22) +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(color = 'black', fill = NA),
    axis.text.y = element_text(size = 10, color = 'black'),
    plot.margin = margin(b = 0, l =2 ),
    axis.text.x = element_text(size = 10, color = 'black', angle = 60, hjust = 1),
    legend.text = element_text(size = 9, color= 'black'),
    legend.title = element_text(size = 9, color = 'black'),
    axis.title = element_text(size = 10, color = 'black')
  ) +
  labs(x = 'Orthogroup', y = 'Average paralog count', fill = "Occupancy")

most_variance_orthos_core <- ortho_count_filt_core %>% dplyr::pull(Orthogroup)

ogs_IPR_count_core <- genes_ogs_ipr %>% dplyr::group_by(Orthogroup, IPR_description) %>%
  dplyr::mutate(ipr_per_og_count = n()) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(Orthogroup, IPR_description, ipr_per_og_count) %>%
  dplyr::group_by(Orthogroup) %>%
  dplyr::arrange(desc(ipr_per_og_count)) %>%
  dplyr::slice_head(n = 1) %>%
  dplyr::mutate(total_ipr_perOG = sum(ipr_per_og_count),
                prop_IPR_perOG = ipr_per_og_count / total_ipr_perOG) %>%
  dplyr::filter(Orthogroup %in% most_variance_orthos_core) %>%
  dplyr::mutate(Orthogroup = factor(Orthogroup, levels = most_variance_orthos_core)) %>%
  dplyr::arrange(Orthogroup) %>%
  dplyr::mutate(IPR_description = factor(IPR_description, levels = unique(IPR_description)))

my_colors_core <- c(
  "#E41A1C",  
  "#377EB8",  
  "#4DAF4A",  
  "#984EA3",  
  "#FF7F00",  
  "#A65628",  
  "#F781BF",  
  "#00A6A6",
  "cyan",  
  "#666666", 
  "#6A3D9A", 
  "#1B9E77", 
  "#D95F02")

bottom_core <- ggplot(ogs_IPR_count_core) + 
  geom_col(aes(x = Orthogroup, y = prop_IPR_perOG, fill = IPR_description)) +
  scale_fill_manual(values = my_colors_core, na.value = "black") +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(color = 'black', fill = NA),
    axis.text.x = element_text(size = 3, color = 'black', angle = 60, hjust = 1),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    plot.margin = margin(t = 0, b = 2, r =2),
    legend.key.size = unit(0.45, "cm"),
    legend.text = element_text(size = 5.5, color= 'black'),
    legend.title = element_blank(),
    axis.title.x = element_text(size = 10, color = 'black')) +
  guides(fill = guide_legend(nrow = 14)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(fill = "IPR description")

legend_core <- cowplot::get_legend(bottom_core)

combined_core <- cowplot::plot_grid(
  top_core + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank()), bottom_core + theme(legend.position = "none"),  
  align = "v",
  nrow = 2,
  rel_heights = c(1, 0.2))

combined_final_core <- cowplot::ggdraw(combined_core) +
  cowplot::draw_grob(
    legend_core,
    x = 0.25,
    y = 0.41,
    width = 0.22,
    height = 0.65)


# Combinging paralog variation plots with 50% occupancy threshold and 100% occupancy threshold (core orthogroups)
final_plt <- cowplot::plot_grid(
  combined_final, combined_final_core,
  nrow = 2,
  rel_heights = c(0.8,1),
  align = "v",
  labels = c("a","b"))

# Save the final plot
ggsave("../../figures/supplementary/paralog_variation.png", final_plt, width = 7.5, height = 9, dpi = 600)
