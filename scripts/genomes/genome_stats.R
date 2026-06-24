library(readr)
library(ggplot2)
library(dplyr)

data <- readr::read_tsv("../../tables/wild_strain_genome_stats.tsv")

data1 <- data %>%
  dplyr::select(strain = Strain, ctg_N50 = `Contig N50`, ctg_L90 = `Contig L90`, fold_cov = `Coverage`, hifi_mean_readlen = `Mean read length`) %>%
  dplyr::mutate(ctg_N50 = as.numeric(ctg_N50), ctg_L90 = as.numeric(ctg_L90), fold_cov = as.numeric(fold_cov), hifi_mean_readlen = as.numeric(hifi_mean_readlen)) %>%
  dplyr::mutate(fold_cov = ifelse(strain == "ECA741", 107.68495, fold_cov))

mean_N50 <- mean(data1$ctg_N50)

color_limits <- range(data1$ctg_N50 / 1e6, na.rm = TRUE)

final_plt <- ggplot(data1) + 
  geom_vline(xintercept = mean_N50 / 1e6, linetype = 'solid', color = 'black', linewidth = 0.5) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_point(aes(x = ctg_N50 / 1e6, y = ctg_L90, 
                 color = fold_cov, size = hifi_mean_readlen)) + 
  scale_color_gradient(name = "Coverage", low = "deepskyblue", high = "red") +
  annotate("text", x = (mean_N50 / 1e6) + 5.5, y = 40, label = sprintf("Mean N50 (Mb): %.2f", mean_N50 / 1e6), hjust = 1.09, vjust = -0.3, size = 4, color = "black") +
  annotate("text", x = 5.5, y = 80, label = "1 Mb threshold", hjust = 1.09, vjust = -0.3, size = 4, color = "black") +
  scale_size_continuous(name = "Mean read length", range = c(0.2, 4)) + 
  theme(
    legend.position = c(0.97, 0.97),
    legend.justification = c("right", "top"),
    legend.background = element_rect(fill = "white"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    panel.background = element_blank(),
    axis.line = element_line(),
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 13, color = "black")
  ) +
  ylab("Contig L90") + 
  xlab("Contig N50 (Mb)") +
  coord_cartesian(ylim = c(0,115))
final_plt

# Save the plot!
ggsave("../../figures/supplementary/genome_stats.png", final_plt, width = 7.5, height = 6, dpi = 600)
