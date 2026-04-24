library(cowplot)
library(readr)
library(dplyr)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)

####################################################################################################
####################################################################################################

# RELATEDNESS MATRIX

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






