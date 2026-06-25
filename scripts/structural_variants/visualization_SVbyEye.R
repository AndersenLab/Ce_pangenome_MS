library(readr)
library(ggplot2)
library(dplyr)
library(cowplot)
library(dplyr)
library(SVbyEye)

#############################################################################
# The largest Inversion manually discovered with genome-genome alignment
#############################################################################
# Read in data
paf.eca_table <- readPaf(paf.file = "../../processed_data/structural_variants/ECA3088_N2.paf", include.paf.tags = T, restrict.paf.tags = "cg")
paf.filtered <- filterPaf(paf.table = paf.eca_table) %>% dplyr::filter(t.name == "IV" & q.name == "ptg000002l")
SV_EYE <- plotMiro(paf.table = paf.filtered, color.by = "strand", add.alignment.arrows = FALSE) + 
  ggplot2::scale_color_manual(values = c("+" = "black", "-" = "red")) +
  ggplot2::scale_fill_manual(values = c("+" = "black", "-" = "red")) +
  ggplot2::theme(
    legend.position = 'none',
    axis.text = ggplot2::element_text(size = 10, color = 'black'),
    axis.title = ggplot2::element_text(size = 12, color = 'black'),
    axis.text.y = ggplot2::element_text(size = 11, color = 'black'))
SV_EYE

plotMiro(paf.table = paf.filtered, binsize = 1000) # looking at percent identity


#############################################################################
# The largest Inversion manually discovered with genome-genome alignment
#############################################################################
nucmer <- readr::read_tsv("../../processed_data/genome_resources/genome_data/141_nucmer_ECA741CGC1.tsv", col_names = c("N2S","N2E","WSS","WSE","L1","L2","IDY","LENR","LENQ","N2_chr","contig","strain")) %>%
  dplyr::filter(strain == "ECA3088", N2_chr == "IV", contig == "ptg000002l") %>%
  dplyr::mutate(inv = ifelse(WSE < WSS, TRUE, FALSE))

largest_INV <- ggplot2::ggplot(nucmer) +
  ggplot2::geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6, color = inv), linewidth = 0.7) +
  ggplot2::scale_color_manual(values = c("TRUE" = "red", "FALSE"= "black")) +
  ggplot2::theme_bw() +
  ggplot2::facet_wrap(~N2_chr) +
  ggplot2::theme(
  legend.position = "none",
  strip.text = ggplot2::element_text(size = 14, color = 'black'),
  axis.text = ggplot2::element_text(size = 10, color = 'black'),
  axis.ticks = ggplot2::element_blank(),
  axis.title = ggplot2::element_text(size = 12, color = 'black'),
  panel.background = ggplot2::element_blank(),
  panel.grid = ggplot2::element_blank(),
  panel.border = ggplot2::element_rect(fill = NA)) +
  ggplot2::labs(x = "N2 genome position (Mb)", y = "ECA3088 contig position (Mb)") +
  scale_y_continuous(expand = c(0.008,0.008)) +
  scale_x_continuous(expand = c(0.008,0.008))
largest_INV


final_plt <- cowplot::plot_grid(
  largest_INV, SV_EYE,
  axis = "l",
  align = "v",
  rel_heights = c(1.3,1),
  labels = c("a","b"),
  nrow = 2)
final_plt


# Save the plot
ggsave("../../figures/supplementary/largest_INV_ECA3088.png", final_plt, width = 7.5, height = 8.5, dpi = 600)

