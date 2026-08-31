library(dplyr)
library(ggplot2)
library(readr)
library(ggplot2)
library(tidyr)
library(dplyr)
library(ape)
library(circlize)
library(data.table)
library(grid) 
library(ggrepel)


####################################################################################################
####################################################################################################

# SV cumulative lengths

####################################################################################################
####################################################################################################
allcalls <- readr::read_tsv("../../processed_data/structural_variants/141_over50_PASS_variants.tsv", col_names = c("chrom", "pos", "ref", "alt", "filter", "sv_type","sv_length","strain")) %>% dplyr::select(-filter) %>%
  dplyr::filter(chrom != "MtDNA")

filt_calls <- allcalls %>% 
  dplyr::mutate(sv_length = abs(sv_length)) %>%
  dplyr::group_by(sv_type, strain) %>%
  dplyr::mutate(type_count = n()) %>%
  dplyr::ungroup()

levels3 <- filt_calls %>%
  dplyr::group_by(strain) %>%
  dplyr::mutate(total_SV_length = sum(sv_length)) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(strain,total_SV_length) %>%
  dplyr::arrange(total_SV_length) %>%
  dplyr::distinct(strain) %>%
  dplyr::pull()

total_size2 <- filt_calls %>%
  dplyr::select(sv_type,sv_length,strain) %>%
  dplyr::group_by(strain,sv_type) %>%
  dplyr::mutate(type_total_length = sum(sv_length)) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(strain,sv_type, type_total_length) %>%
  dplyr::mutate(strain = factor(strain, levels = levels3)) %>%
  dplyr::mutate(sv_type = factor(sv_type, levels = c("INV","DEL","INS"))) 

SV_LEN <- ggplot() +
  geom_bar(data = total_size2, aes(x = type_total_length / 1e6, y = strain, fill = sv_type), stat = "identity") +
  scale_fill_manual(values = c("INS" = "blue", "DEL" = "red", "INV" = "gold")) +
  labs(x = "Length of SVs (Mb)", fill = "SV type", y = "Strains") +
  theme(panel.border = element_rect(color = 'black', fill = NA),
        panel.background = element_blank(),
        axis.text = element_text(size = 10, color = 'black'),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.title = element_text(size = 10, color = 'black'),
        legend.box.background = element_rect(color = 'black', fill = NA),
        legend.position = "inside",
        legend.position.inside = c(0.8,0.111),
        legend.text = element_text(size = 10, color = 'black'),
        axis.title.x = element_text(size = 10, color = 'black')) +
  coord_cartesian(xlim = c(0,13)) +
  scale_x_continuous(breaks = seq(0,14, by =2), expand = c(0,0))
SV_LEN

# Write table of cummulative SV legnth of each type for each strain
# write.table(total_size2,"../../processed_data/structural_variants/strain_SVs_cum_length.tsv", sep = '\t', row.names = F, col.names = T, quote = F)



####################################################################################################
####################################################################################################

# PCA OF SVs

####################################################################################################
####################################################################################################
geo_initial <- readr::read_tsv("../../processed_data/genome_resources/isotypes/elegans_isotypes_sampling_geo.tsv")
hawaii_islands <- readr::read_tsv("../../processed_data/genome_resources/isotypes/elegans_isotypes_sampling_geo_hawaii_islands.tsv") %>% dplyr::select(isotype,collection_island_Hawaii)
WSs <- readr::read_tsv("../../tables/wild_strain_genome_stats.tsv") %>% dplyr::select(Strain) %>% dplyr::rename(strain = Strain) %>% dplyr::pull()
merged_SV <- readr::read_tsv("../../processed_data/structural_variants/Jasmine_merged_SVs.tsv")

# Adding Hawaiian island resolution
geo <- geo_initial %>%
  dplyr::left_join(hawaii_islands, by = "isotype") %>%
  dplyr::mutate(geo = ifelse(geo == "Hawaii",collection_island_Hawaii,geo)) %>%
  dplyr::select(isotype, lat, long, geo) %>%
  dplyr::filter(isotype %in% WSs)

# Filter for MAF > 0.05
n <- 141
maf <- 0.05
least <- ceiling(maf * n)      #8
most  <- floor((1 - maf) * n)  #133

common_vcf <- merged_SV %>% 
  dplyr::filter(number_svs_merged >= least & number_svs_merged <= most) %>%
  dplyr::select(-chrom, -pos, -ref, -alt, -sv_type, -sv_length, -number_svs_merged) %>%
  as.matrix()

common_vcf[common_vcf == "./."] <- 0
sv_mat <- apply(common_vcf, 2, as.numeric)

sv_mat_t <- t(sv_mat)

sv_mat_t[is.na(sv_mat_t)] <- colMeans(sv_mat_t, na.rm = TRUE)

sv_pca <- prcomp(sv_mat_t, center = TRUE, scale. = TRUE) 


pca_df <- as.data.frame(sv_pca$x)
pca_df$strain <- rownames(pca_df)

strain_geo <- geo %>% dplyr::rename(strain = isotype) %>% dplyr::select(strain,geo)

pca_df <- pca_df %>%
  dplyr::left_join(strain_geo, by = "strain") %>%
  dplyr::mutate(geo = ifelse(strain == "CGC1", "CGC1",geo)) 

geo.colors <- c("Big Island"="black", "Molokai" = "#66C2A5", "Maui" = "yellow", "Oahu" = "brown", "Kauai" = "purple", "Africa"="green", "North America" = "pink", "Europe" = "#E41A1C", "Atlantic" = "blue", 
                "Oceania" ="cyan", "unknown" = 'gray', "CGC1" = "#DB6333")

pca_df <- pca_df %>%
  dplyr::mutate(label = ifelse(PC2 > 60, strain, NA))

PCA <- ggplot(pca_df, aes(PC1, PC2, color = geo)) +
  geom_text_repel(aes(label = label), size = 2, max.overlaps = Inf, show.legend = FALSE) +
  geom_point(size = 1, alpha = 0.8) +
  scale_color_manual(values = geo.colors) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 10, color = 'black'),
    axis.title = element_text(size = 11, color = 'black'),
    legend.text = element_text(size = 9, color = 'black'),
    legend.title = element_text(size = 9, color = 'black'),
    legend.key.width = unit(0.39, "cm"),
    legend.position = 'top') +
  labs(color = "Collection\nlocation", x = paste0("PC1 (", round(100 * summary(sv_pca)$importance[2,1], 1), "%)"), y = paste0("PC2 (", round(100 * summary(sv_pca)$importance[2,2], 1), "%)"))+
  guides(color = guide_legend(override.aes = list(size = 2), keyheight = unit(0.2, "cm"))) 
PCA


# Create final plot
final_plt <- cowplot::plot_grid(
  SV_LEN, PCA,
  rel_widths = c(0.5,1),
  nrow = 1,
  labels = c("a","b"))
final_plt

# Save the plot
ggsave("../../figures/structural_variants.png", final_plt, width = 7.5, height = 6, dpi = 600)

