library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
library(cowplot)

# Read in the strains isolation world map

# Read in the relatedness matrix plot


# Reading in genome stats
stats <- readr::read_tsv("../../processed_data/genome_resources/wild_strain_genome_stats.tsv") %>%
  dplyr::select(Strain,`Number of contigs`,`Genome size`, `Contig N50`, `Contig L90`, `Contig N90`, `Coverage`, `BUSCO completeness (genome)`) %>%
  dplyr::rename(strain = Strain, n_contigs = `Number of contigs`, contig_bp = `Genome size`, ctg_N50 = `Contig N50`, ctg_L90 = `Contig L90`, ctg_N90 = `Contig N90`) %>%
  dplyr::filter(!strain %in% bad_strains)

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
  # stat_summary(fun=mean, geom="crossbar", width=0.5, color="black") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    panel.grid = element_blank(),
    # panel.grid.major.y = element_line(color = 'gray'),
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
  # stat_summary(fun=mean, geom="crossbar", width=0.5, color="black") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    panel.grid = element_blank(),
    # panel.grid.major.y = element_line(color = 'gray'),
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

final_stats <- cowplot::ggdraw(base) +
  # Upper break marks
  draw_line(x = c(0.145, 0.16), y = c(0.535, 0.540), size = 0.6) +
  draw_line(x = c(0.145, 0.16), y = c(0.47, 0.475), size = 0.6)
final_stats


# Concatenate to make final plot


# Good to go!
# ggsave("../../figures/strain_selection_genome_stats.png", final_plot, width = 7.5, height = 9, dpi = 600)






