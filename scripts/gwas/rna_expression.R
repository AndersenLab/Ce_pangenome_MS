library(readr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggrepel)
library(ape)

# Wild strain glc-1 orthologs
glc_orthos <- readr::read_tsv("../../processed_data/gwas/glc1_orthos.tsv")

# Load in fold-coverage data for every strain
cov <- readr::read_tsv("../../processed_data/gwas/rna_aln/coverage/ALL_strains_cov.txt", col_names = c("contig","start","end","gene","total_cov","mean_cov","strain"))

n2_wb_geneID <- ape::read.gff("../..//processed_data/genome_resources/annotation/c_elegans.PRJNA13758.WS283.csq.PCfeaturesOnly.longest.gff3") %>%
  dplyr::filter(type == "gene") %>%
  dplyr::select(attributes) %>%
  dplyr::mutate(
    wbgene = str_extract(attributes, "(?<=ID=gene:)[^;]+"),
    gene = str_extract(attributes, "(?<=sequence_name=)[^;]+")) %>%
  dplyr::select(wbgene, gene)

cov <- cov %>% dplyr::left_join(n2_wb_geneID, by = c("gene" = "wbgene")) %>%
  dplyr::mutate(gene_final = ifelse(is.na(gene.y) & strain != "N2", gene, gene.y)) %>%
  dplyr::filter(!is.na(gene_final)) %>%
  dplyr::select(contig, start, end, gene = gene_final, total_cov, mean_cov, strain)

# Left-join to get fold-coverage
glc_cov <- glc_orthos %>% dplyr::left_join(cov, by = c("strain","gene")) %>% dplyr::filter(!is.na(mean_cov))

alt_strains <- c("WN2001","NIC199","MY23","MY16","JU2526","CX11314","JU775","MY10","JU258","CX11264","NIC195","CB4856","JU782", "EG4725", "DL238", "ECA36")

pheno_strains <- readr::read_tsv("../..//processed_data/gwas/strains.tsv", col_names = "strain") %>%
  dplyr::mutate(geno = ifelse(strain %in% alt_strains, "ALT", "REF")) 

pheno_strains_vec <- pheno_strains %>% dplyr::pull(strain)

strain_average_cov <- cov %>% dplyr::filter(strain %in% pheno_strains_vec) %>% 
  dplyr::group_by(strain) %>%
  dplyr::mutate(mean_exome_cov = mean(mean_cov)) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(strain, mean_exome_cov)
 
rna_cov_strains <- glc_cov %>% dplyr::left_join(pheno_strains, by = "strain") %>%
  dplyr::select(strain, gene, alias, total_cov, mean_cov, geno) %>%
  dplyr::left_join(strain_average_cov, by = "strain") %>%
  dplyr::mutate(norm_cov = mean_cov / mean_exome_cov) %>%
  dplyr::mutate(geno = factor(geno, levels = c("REF", "ALT"))) 

cov_plt <- ggplot(rna_cov_strains) +
  geom_boxplot(aes(x = geno, y = norm_cov, fill = geno), alpha = 0.5, outliers = FALSE) +
  geom_jitter(aes(x = geno, y = norm_cov), width = 0.05) +
  scale_fill_manual(values = c("REF" = "black", "ALT" = "red")) +
  ggrepel::geom_text_repel(aes(x = geno, y = norm_cov, label = strain), size = 2, max.overlaps = Inf, segment.color = NA) +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(color = 'black', fill = NA),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 10, color = 'black'),
    axis.title.x = element_blank(),
    legend.position = 'none',
    axis.title.y = element_text(size = 12, color ='black')) +
  labs(y = "Normalized RNA coverage")

rna_cov_glc_strain <- rna_cov_strains %>%
  dplyr::group_by(strain, geno) %>%
  dplyr::summarise(normalized_glc1_cov = median(norm_cov, na.rm = TRUE), .groups = "drop")

wilcox_stats <- rna_cov_glc_strain %>%
  rstatix::wilcox_test(normalized_glc1_cov ~ geno, paired = FALSE) %>%
  rstatix::add_significance("p") %>%
  rstatix::add_xy_position(x = "geno", step.increase = 0.2)

cov_plt <- cov_plt +
  ggpubr::stat_pvalue_manual(wilcox_stats, label = "p.signif", size = 4, vjust = 0.1, tip.length = 0.01)
cov_plt

# Look at fold-coverage differences in F11A5.5
f11 <- readr::read_tsv("../../processed_data/gwas/F11A5.5_orthos.tsv") 

f11_cov <- f11 %>% 
  dplyr::left_join(cov, by = c("strain","gene")) %>% dplyr::filter(!is.na(mean_cov)) %>%
  dplyr::left_join(pheno_strains, by = "strain") %>%
  dplyr::select(strain, gene, alias, total_cov, mean_cov, geno) %>%
  dplyr::left_join(strain_average_cov, by = "strain") %>%
  dplyr::mutate(norm_cov = mean_cov / mean_exome_cov) %>%
  dplyr::mutate(geno = factor(geno, levels = c("REF", "ALT"))) 

f11_cov_plt <- ggplot(f11_cov) +
  geom_boxplot(aes(x = geno, y = norm_cov, fill = geno), alpha = 0.5, outliers = FALSE) +
  geom_jitter(aes(x = geno, y = norm_cov), width = 0.05) +
  scale_fill_manual(values = c("REF" = "black", "ALT" = "red")) +
  ggrepel::geom_text_repel(aes(x = geno, y = norm_cov, label = strain), size = 2, max.overlaps = Inf, segment.color = NA) +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(color = 'black', fill = NA),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 10, color = 'black'),
    axis.title.x = element_blank(),
    legend.position = 'none',
    axis.title.y = element_text(size = 12, color ='black')) +
  labs(y = NULL)

rna_cov_F11_strain <- f11_cov %>%
  dplyr::group_by(strain, geno) %>%
  dplyr::summarise(normalized_F11A5.5_cov = median(norm_cov, na.rm = TRUE), .groups = "drop")

wilcox_stats_F11 <- rna_cov_F11_strain %>%
  rstatix::wilcox_test(normalized_F11A5.5_cov ~ geno, paired = FALSE) %>%
  rstatix::add_significance("p") %>%
  rstatix::add_xy_position(x = "geno", step.increase = 0.2) 

f11_cov_plt <- f11_cov_plt +
  ggpubr::stat_pvalue_manual(wilcox_stats_F11, label = "p.signif", size = 4, vjust = 0.1, tip.length = 0.01)
f11_cov_plt


# Make a cowplot:
final_plt <- cowplot::plot_grid(
  cov_plt, f11_cov_plt,
  nrow = 1,
  align = "h",
  labels = c("a","b")
)
final_plt

ggsave("../../figures/supplementary/fold_cov_genesInInterval.png", final_plt, width = 7.5, height = 6, dpi = 600)



