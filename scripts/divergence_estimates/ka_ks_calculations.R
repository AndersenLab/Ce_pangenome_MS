library(dplyr)
library(ggplot2)
library(seqinr)
library(tidyr)
library(rstatix)
library(ggpubr)

# ============================================================================================================================================================================================ #
# For all 12,315 orthogroups
# ============================================================================================================================================================================================ #
# Setting variables and reading in necessary data:
aln_dir <- "../../processed_data/divergence_estimates/codon_alignments"
suffix <- ".aligned.codons.gapped.fa"
og_id_list <- readr::read_tsv("../../processed_data/orthology/orthofinder/orthofinder_output/Orthogroups_SingleCopyOrthologues.txt", col_names = "OGs") %>% dplyr::pull(OGs)
hdr_sc_ogs <- readr::read_tsv("../../processed_data/hdr_liftover/hdr_sc_ogs.tsv") %>%
  dplyr::mutate(inHDR = "HDR") %>% dplyr::rename(OG = sc_hdr_ogs)
all_ws_genes_class_og <- readr::read_tsv("../../processed_data/orthology/all_genes_class_OGs.tsv")

# Function to read in file and filter for occupancy
read_and_filter_alignment <- function(aln_file) {
  aln <- seqinr::read.alignment(aln_file, format = "fasta")
  
  seqs <- toupper(aln$seq)
  
  occupancy <- vapply(seqs, function(x) {
    chars <- strsplit(x, "", fixed = TRUE)[[1]]
    sum(chars %in% c("A", "C", "G", "T")) / length(chars) # total number of nucleotides / total number of positions (nucleotides and gaps "-")
  }, numeric(1))
  
  cds_len <- vapply(seqs, function(x) {
    chars <- strsplit(x, "", fixed = TRUE)[[1]]
    sum(chars %in% c("A", "C", "G", "T")) / 3 # the number of codons present
  }, numeric(1))
  
  keep <- occupancy >= 0.7 & cds_len >= 100
  
  if (sum(keep) < 127) return(NULL) # must have at least four orthologs to calculate Ks and Ka
  
  aln$nam <- aln$nam[keep]
  aln$seq <- aln$seq[keep]
  aln$nb <- length(aln$seq)
  
  return(aln)
}

# Function to calculate Ka and Ks for each OG and create a nice summary table
get_kaks <- function(og_id) {
  
  aln_file <- file.path(aln_dir, paste0(og_id, suffix))

  if (!file.exists(aln_file)) {
    warning("Skipping ", og_id, ": no alignment found")
    return(NULL) }
  
  aln <- read_and_filter_alignment(aln_file)
  
  if (is.null(aln)) {
    warning("Skipping ", og_id, ": fewer than 4 sequences after filtering")
    return(NULL) }
  
  res <- tryCatch(
    seqinr::kaks(aln, rmgap = TRUE),
    error = function(e) {
      warning("Skipping ", og_id, ": ", e$message)
      return(NULL) } )
  
  if (is.null(res)) {
    return(NULL) }
  
  ka_mat <- as.matrix(res$ka)
  ks_mat <- as.matrix(res$ks)
  
  # Only take the upper triangle
  pair_idx <- upper.tri(ka_mat)
  
  ka <- ka_mat[pair_idx]
  ks <- ks_mat[pair_idx]
  
  kaks <- ifelse(
    ks > 0,
    ka / ks,
    NA_real_
  )
  
  summary <- data.frame(
    OG = og_id,
    mean_pairwise_ks = mean(ks, na.rm = TRUE),
    mean_pairwise_ka = mean(ka, na.rm = TRUE),
    mean_pairwise_kaks = mean(kaks, na.rm = TRUE))
  
  final_stats <- summary %>%
    tidyr::pivot_longer(
      cols = c(
        mean_pairwise_ks,
        mean_pairwise_ka,
        mean_pairwise_kaks), 
      names_to = "stats", 
      values_to = "values") %>%
    dplyr::mutate(stats = factor(stats, levels = c("mean_pairwise_ks", "mean_pairwise_ka", "mean_pairwise_kaks")))

  return(final_stats)
}

# Iterating over all 12,315 OGs to create a master summary file
final_kaks_list <- vector("list", length(og_id_list))
c <- 1
for (i in seq_along(og_id_list)) {
  # print(og_id_list[i])
  print(paste0("On orthogroup ", c, " out of 12,315"))
  c <- c + 1
  final_kaks_list[[i]] <- get_kaks(og_id_list[i])
}

# Collapse list into a dataframe 
final_kaks_df <- dplyr::bind_rows(final_kaks_list) %>%
  dplyr::left_join(hdr_sc_ogs, by = "OG") %>% dplyr::mutate(inHDR = ifelse(is.na(inHDR), "non-HDR", inHDR),
                                                            count_inHDR = ifelse(is.na(count_inHDR), 0, count_inHDR))

# How many of the 12,315 SC OGs passed QC for Ka and Ks calculations?
print(length(unique(final_kaks_df$OG))) # 11,631
stats <- final_kaks_df %>% dplyr::distinct(OG, inHDR) %>% dplyr::group_by(inHDR) %>% dplyr::summarize(number_of_OGs = n()) # 2,124 in HDRs and 9,507 outside HDRs

kaks_plt <- final_kaks_df %>% dplyr::filter(values > 0) %>% dplyr::mutate(stats = ifelse(stats == "mean_pairwise_ks", "Ks",
                                                                                         ifelse(stats == "mean_pairwise_ka", "Ka", "Ka / Ks"))) %>%
  dplyr::mutate(stats = factor(stats, levels = c("Ks", "Ka", "Ka / Ks")))

# Wilcoxon tests for each statistic
wilcox_results <- kaks_plt %>%
  dplyr::group_by(stats) %>%
  rstatix::wilcox_test(values ~ inHDR) %>%
  rstatix::adjust_pvalue(method = "BH") %>%  
  rstatix::add_significance()

effect_sizes <- kaks_plt %>%
  dplyr::group_by(stats) %>%
  rstatix::wilcox_effsize(values ~ inHDR) %>%
  dplyr::select(-.y., -group1, -group2, -n1, -n2)

stats_final <- left_join(wilcox_results, effect_sizes, by = "stats")

annotation_df <- stats_final %>%
  dplyr::mutate(label = case_when(
    p.adj < 0.0001 ~ "****",
    p.adj < 0.001 ~ "***",
    p.adj < 0.01 ~ "**",
    p.adj < 0.05 ~ "*",
    TRUE ~ "ns"), y_pos = 6)

# Difference in Ks, Ka, and Ka / Ks between SC OGs in HDRs and not in HDRs?
kaks_all_sc_ogs <- ggplot(kaks_plt %>% dplyr::filter(count_inHDR == 0 | count_inHDR >= 10), aes(x = stats, y = values, fill = inHDR)) +
  geom_boxplot(outlier.size = 0.6, width = 0.7, position = position_dodge(width = 0.75), outliers = FALSE, alpha = 0.5) +
  geom_point(position = position_jitterdodge(jitter.width = 0.25, dodge.width = 0.7), size = 0.5, alpha = 0.5) +
  geom_text(data = annotation_df, aes(x = stats, y = y_pos, label = label), inherit.aes = FALSE, size = 5) +
  # geom_hline(yintercept = 1, linetype = "dashed") +
  scale_fill_manual(values = c("non-HDR" = "blue", "HDR" = "RED")) +
  scale_y_log10(labels = scales::label_number()) +
  theme(
    axis.text.x = element_text(size = 12, color = 'black'),
    axis.text.y = element_text(size = 10, color = 'black'),
    axis.title.y = element_text(size = 12, color = 'black'),
    panel.border = element_rect(color = 'black', fill = NA),
    legend.position = 'inside',
    legend.position.inside = c(0.9, 0.2),
    legend.text = element_text(size = 10, color = 'black'),
    panel.background = element_blank()
  ) +
  labs(x = NULL, y = expression("Mean pairwise value (log"[10]*")"), fill = NULL)
# kaks_all_sc_ogs

# Save the plot
# ggsave("../../figures/supplementary/KaKs_all_singleCopy_OGs.png", kaks_all_sc_ogs, width = 7.5, height = 6, dpi = 600)

# Does sequence divergence (Ks) and non-synonymous rate (Ka) scale with frequence SC OG is found in a HDR?
ggplot(kaks_plt %>% dplyr::filter(stats == "Ks",
                                  inHDR == "HDR"), aes(x = count_inHDR, y = values)) +
  geom_point(size = 2) +
  theme(
    axis.text = element_text(size = 14, color = 'black'),
    axis.title = element_text(size = 16, color = 'black'),
    panel.border = element_rect(color = 'black', fill = NA),
    legend.position = 'inside',
    legend.position.inside = c(0.9, 0.9),
    legend.text = element_text(size = 16, color = 'black'),
    panel.background = element_blank()
  ) +
  labs(x = "Number of wild strains a single-copy orthogroup is in a HDR", y = "Mean pairwise Ks")

ggplot(kaks_plt %>% dplyr::filter(stats == "Ka",
                                  inHDR == "HDR"), aes(x = count_inHDR, y = values)) +
  geom_point(size = 2) +
  theme(
    axis.text = element_text(size = 14, color = 'black'),
    axis.title = element_text(size = 16, color = 'black'),
    panel.border = element_rect(color = 'black', fill = NA),
    legend.position = 'inside',
    legend.position.inside = c(0.9, 0.9),
    legend.text = element_text(size = 16, color = 'black'),
    panel.background = element_blank()
  ) +
  labs(x = "Number of wild strains a single-copy orthogroup is in a HDR", y = "Mean pairwise Ka")









# ============================================================================================================================================================================================ #
# Which single-copy OGs are undergoing positive selection the most?
# ============================================================================================================================================================================================ #
kaks_highest <- kaks_plt %>% dplyr::filter(stats == "Ka / Ks" & values >= 1) %>% dplyr::arrange(desc(values)) 
kaks_highest_OGs <- kaks_highest %>% dplyr::pull(OG)

kaks_highest_OGs_genes <- all_ws_genes_class_og %>% dplyr::filter(Orthogroup %in% kaks_highest_OGs) %>%
  dplyr::mutate(strain_transcript = paste0(strain,"_", gene))

ipr_annotations <- readr::read_tsv("../../tables/IPR_annotation_142strains.tsv", col_names = c("strain_transcript", "MD5_digest",	"Seq_length",	"Application",	"Signature_accession",	"Signature_description",	"Amino_acid_start_position", "Amino_acid_end_position", "Score", "status", "date", "IPR_accession", "IPR_description", "GO_ID", "pathway")) %>%
  dplyr::filter(IPR_description != "-") %>%
  dplyr::select(strain_transcript, IPR_accession, IPR_description) %>%
  dplyr::mutate(strain_transcript = sub("\\.[^.]*$", "", strain_transcript)) %>%
  dplyr::mutate(strain_transcript = sub("[a-z]+$", "", strain_transcript))

kaks_highest_ipr_annotations <- kaks_highest_OGs_genes %>% dplyr::left_join(ipr_annotations, by = "strain_transcript") %>%
  dplyr::distinct(Orthogroup, strain_transcript, IPR_description) %>%
  dplyr::group_by(Orthogroup, IPR_description) %>% dplyr::mutate(count_IPR_inOG = n()) %>%
  dplyr::arrange(desc(count_IPR_inOG)) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(Orthogroup, IPR_description, count_IPR_inOG) %>%
  dplyr::filter(!is.na(IPR_description) | count_IPR_inOG > 140) %>% # 93 OGs don't have an IPR annotation
  dplyr::left_join(kaks_highest, by = c("Orthogroup" = "OG")) %>%
  dplyr::arrange(desc(values)) 

most_common_ipr <- kaks_highest_ipr_annotations %>% dplyr::group_by(IPR_description) %>%
  dplyr::mutate(count_ipr = n()) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(IPR_description, stats, values, count_ipr) %>%
  dplyr::arrange(desc(count_ipr)) %>%
  dplyr::distinct(IPR_description, count_ipr)




# ============================================================================================================================================================================================ #
# Which single-copy OGs are found in HDRs the most???
# ============================================================================================================================================================================================ #
hdr_most <- kaks_plt %>% dplyr::filter(inHDR == "HDR", stats == "Ka / Ks") %>% 
  dplyr::arrange(desc(count_inHDR)) %>% dplyr::pull(OG)

hdr_most_genes <- all_ws_genes_class_og %>% dplyr::filter(Orthogroup %in% hdr_most) %>%
  dplyr::mutate(strain_transcript = paste0(strain,"_", gene))

hdr_most_ipr_anno <- hdr_most_genes %>% dplyr::left_join(ipr_annotations, by = "strain_transcript") %>%
  dplyr::distinct(Orthogroup, strain_transcript, IPR_description) %>%
  dplyr::group_by(Orthogroup, IPR_description) %>% dplyr::mutate(count_IPR_inOG = n()) %>%
  dplyr::arrange(desc(count_IPR_inOG)) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(Orthogroup, IPR_description, count_IPR_inOG) %>% 
  dplyr::left_join(kaks_plt %>% dplyr::filter(inHDR == "HDR", stats == "Ka / Ks"), by = c("Orthogroup" = "OG"))



# ============================================================================================================================================================================================ #
# For GPCRs
# ============================================================================================================================================================================================ #
ws_gpcrs <- readr::read_tsv("../../processed_data/genome_resources/annotation/140_wild_strains_IPR_gpcrs.tsv") 
sc_gpcrs <- ws_gpcrs %>% dplyr::left_join(all_ws_genes_class_og, by = c("strain","gene")) %>% dplyr::filter(Orthogroup %in% og_id_list) %>%
  dplyr::select(OG = Orthogroup) %>% dplyr::distinct() %>% dplyr::pull()

kaks_plt_gpcr <- kaks_plt %>% dplyr::filter(OG %in% sc_gpcrs)

# Wilcoxon tests for each statistic
wilcox_results_gpcr <- kaks_plt_gpcr %>%
  dplyr::group_by(stats) %>%
  rstatix::wilcox_test(values ~ inHDR) %>%
  rstatix::adjust_pvalue(method = "BH") %>%  
  rstatix::add_significance()

effect_sizes_gpcr <- kaks_plt_gpcr %>%
  dplyr::group_by(stats) %>%
  rstatix::wilcox_effsize(values ~ inHDR) %>%
  dplyr::select(-.y., -group1, -group2, -n1, -n2)

stats_final_gpcr <- left_join(wilcox_results_gpcr, effect_sizes_gpcr, by = "stats")

annotation_df_gpcr <- stats_final_gpcr %>%
  dplyr::mutate(label = case_when(
    p.adj < 0.0001 ~ "****",
    p.adj < 0.001 ~ "***",
    p.adj < 0.01 ~ "**",
    p.adj < 0.05 ~ "*",
    TRUE ~ "ns"), y_pos = 6)

# Difference in Ks, Ka, and Ka / Ks between SC OGs in HDRs and not in HDRs?
ggplot(kaks_plt_gpcr, aes(x = stats, y = values, fill = inHDR)) +
  geom_boxplot(outlier.size = 0.6, width = 0.7, position = position_dodge(width = 0.75), outliers = FALSE, alpha = 0.5) +
  geom_point(position = position_jitterdodge(jitter.width = 0.25, dodge.width = 0.7), size = 1, alpha = 0.5) +
  geom_text(data = annotation_df_gpcr, aes(x = stats, y = y_pos, label = label), inherit.aes = FALSE, size = 5) +
  scale_fill_manual(values = c("non-HDR" = "blue", "HDR" = "RED")) +
  scale_y_log10(labels = scales::label_number()) +
  theme(
    axis.text.x = element_text(size = 16, color = 'black'),
    axis.text.y = element_text(size = 14, color = 'black'),
    axis.title.y = element_text(size = 16, color = 'black'),
    panel.border = element_rect(color = 'black', fill = NA),
    legend.position = 'inside',
    legend.position.inside = c(0.9, 0.2),
    legend.text = element_text(size = 16, color = 'black'),
    plot.title = element_text(size = 22, hjust = 0.5, color = 'black'),
    panel.background = element_blank()
  ) +
  labs(x = NULL, y = expression("Mean pairwise value (log"[10]*")"), fill = NULL, title = "GPCRs")





# ============================================================================================================================================================================================ #
# For F-box genes
# ============================================================================================================================================================================================ #
ws_fbox <- readr::read_tsv("../../processed_data/genome_resources/annotation/140_wild_strains_IPR_fBox.tsv") 
sc_fbox <- ws_fbox %>% dplyr::left_join(all_ws_genes_class_og, by = c("strain","gene")) %>% dplyr::filter(Orthogroup %in% og_id_list) %>%
  dplyr::select(OG = Orthogroup) %>% dplyr::distinct() %>% dplyr::pull()

kaks_plt_fbox <- kaks_plt %>% dplyr::filter(OG %in% sc_fbox)

# Wilcoxon tests for each statistic
wilcox_results_fbox <- kaks_plt_fbox %>%
  dplyr::group_by(stats) %>%
  rstatix::wilcox_test(values ~ inHDR) %>%
  rstatix::adjust_pvalue(method = "BH") %>%  
  rstatix::add_significance()

effect_sizes_fbox <- kaks_plt_fbox %>%
  dplyr::group_by(stats) %>%
  rstatix::wilcox_effsize(values ~ inHDR) %>%
  dplyr::select(-.y., -group1, -group2, -n1, -n2)

stats_final_fbox <- left_join(wilcox_results_fbox, effect_sizes_fbox, by = "stats")

annotation_df_fbox <- stats_final_fbox %>%
  dplyr::mutate(label = case_when(
    p.adj < 0.0001 ~ "****",
    p.adj < 0.001 ~ "***",
    p.adj < 0.01 ~ "**",
    p.adj < 0.05 ~ "*",
    TRUE ~ "ns"), y_pos = 6)

# Difference in Ks, Ka, and Ka / Ks between SC OGs in HDRs and not in HDRs?
ggplot(kaks_plt_fbox, aes(x = stats, y = values, fill = inHDR)) +
  geom_boxplot(outlier.size = 0.6, width = 0.7, position = position_dodge(width = 0.75), outliers = FALSE, alpha = 0.5) +
  geom_point(position = position_jitterdodge(jitter.width = 0.25, dodge.width = 0.7), size = 1, alpha = 0.5) +
  geom_text(data = annotation_df_fbox, aes(x = stats, y = y_pos, label = label), inherit.aes = FALSE, size = 5) +
  scale_fill_manual(values = c("non-HDR" = "blue", "HDR" = "RED")) +
  scale_y_log10(labels = scales::label_number()) +
  theme(
    axis.text.x = element_text(size = 16, color = 'black'),
    axis.text.y = element_text(size = 14, color = 'black'),
    axis.title.y = element_text(size = 16, color = 'black'),
    panel.border = element_rect(color = 'black', fill = NA),
    legend.position = 'inside',
    legend.position.inside = c(0.9, 0.2),
    legend.text = element_text(size = 16, color = 'black'),
    plot.title = element_text(size = 22, hjust = 0.5, color = 'black'),
    panel.background = element_blank()
  ) +
  labs(x = NULL, y = expression("Mean pairwise value (log"[10]*")"), fill = NULL, title = "F-box")


# ============================================================================================================================================================================================ #
# For C-type lectins
# ============================================================================================================================================================================================ #
ws_lectin <- readr::read_tsv("../../processed_data/genome_resources/annotation/140_wild_strains_IPR_CtypeLectins.tsv") 
sc_lectin <- ws_lectin %>% dplyr::left_join(all_ws_genes_class_og, by = c("strain","gene")) %>% dplyr::filter(Orthogroup %in% og_id_list) %>%
  dplyr::select(OG = Orthogroup) %>% dplyr::distinct() %>% dplyr::pull()

kaks_plt_lectin <- kaks_plt %>% dplyr::filter(OG %in% sc_lectin)

# Wilcoxon tests for each statistic
wilcox_results_lectin <- kaks_plt_lectin %>%
  dplyr::group_by(stats) %>%
  rstatix::wilcox_test(values ~ inHDR) %>%
  rstatix::adjust_pvalue(method = "BH") %>%  
  rstatix::add_significance()

effect_sizes_lectin <- kaks_plt_lectin %>%
  dplyr::group_by(stats) %>%
  rstatix::wilcox_effsize(values ~ inHDR) %>%
  dplyr::select(-.y., -group1, -group2, -n1, -n2)

stats_final_lectin <- left_join(wilcox_results_lectin, effect_sizes_lectin, by = "stats")

annotation_df_lectin <- stats_final_lectin %>%
  dplyr::mutate(label = case_when(
    p.adj < 0.0001 ~ "****",
    p.adj < 0.001 ~ "***",
    p.adj < 0.01 ~ "**",
    p.adj < 0.05 ~ "*",
    TRUE ~ "ns"), y_pos = 6)

# Difference in Ks, Ka, and Ka / Ks between SC OGs in HDRs and not in HDRs?
ggplot(kaks_plt_lectin, aes(x = stats, y = values, fill = inHDR)) +
  geom_boxplot(outlier.size = 0.6, width = 0.7, position = position_dodge(width = 0.75), outliers = FALSE, alpha = 0.5) +
  geom_point(position = position_jitterdodge(jitter.width = 0.25, dodge.width = 0.7), size = 1, alpha = 0.5) +
  geom_text(data = annotation_df_lectin, aes(x = stats, y = y_pos, label = label), inherit.aes = FALSE, size = 5) +
  scale_fill_manual(values = c("non-HDR" = "blue", "HDR" = "RED")) +
  scale_y_log10(labels = scales::label_number()) +
  theme(
    axis.text.x = element_text(size = 16, color = 'black'),
    axis.text.y = element_text(size = 14, color = 'black'),
    axis.title.y = element_text(size = 16, color = 'black'),
    panel.border = element_rect(color = 'black', fill = NA),
    legend.position = 'inside',
    legend.position.inside = c(0.9, 0.2),
    legend.text = element_text(size = 16, color = 'black'),
    plot.title = element_text(size = 22, hjust = 0.5, color = 'black'),
    panel.background = element_blank()
  ) +
  labs(x = NULL, y = expression("Mean pairwise value (log"[10]*")"), fill = NULL, title = "C-type lectins")


# ============================================================================================================================================================================================ #
# For NHRs
# ============================================================================================================================================================================================ #
ws_nhr <- readr::read_tsv("../../processed_data/genome_resources/annotation/140_wild_strains_IPR_nhr.tsv") 
sc_nhr <- ws_nhr %>% dplyr::left_join(all_ws_genes_class_og, by = c("strain","gene")) %>% dplyr::filter(Orthogroup %in% og_id_list) %>%
  dplyr::select(OG = Orthogroup) %>% dplyr::distinct() %>% dplyr::pull()

kaks_plt_nhr <- kaks_plt %>% dplyr::filter(OG %in% sc_nhr)

# Wilcoxon tests for each statistic
wilcox_results_nhr <- kaks_plt_nhr %>%
  dplyr::group_by(stats) %>%
  rstatix::wilcox_test(values ~ inHDR) %>%
  rstatix::adjust_pvalue(method = "BH") %>%  
  rstatix::add_significance()

effect_sizes_nhr <- kaks_plt_nhr %>%
  dplyr::group_by(stats) %>%
  rstatix::wilcox_effsize(values ~ inHDR) %>%
  dplyr::select(-.y., -group1, -group2, -n1, -n2)

stats_final_nhr <- left_join(wilcox_results_nhr, effect_sizes_nhr, by = "stats")

annotation_df_nhr <- stats_final_nhr %>%
  dplyr::mutate(label = case_when(
    p.adj < 0.0001 ~ "****",
    p.adj < 0.001 ~ "***",
    p.adj < 0.01 ~ "**",
    p.adj < 0.05 ~ "*",
    TRUE ~ "ns"), y_pos = 6)

# Difference in Ks, Ka, and Ka / Ks between SC OGs in HDRs and not in HDRs?
ggplot(kaks_plt_nhr, aes(x = stats, y = values, fill = inHDR)) +
  geom_boxplot(outlier.size = 0.6, width = 0.7, position = position_dodge(width = 0.75), outliers = FALSE, alpha = 0.5) +
  geom_point(position = position_jitterdodge(jitter.width = 0.25, dodge.width = 0.7), size = 1, alpha = 0.5) +
  geom_text(data = annotation_df_nhr, aes(x = stats, y = y_pos, label = label), inherit.aes = FALSE, size = 5) +
  scale_fill_manual(values = c("non-HDR" = "blue", "HDR" = "RED")) +
  scale_y_log10(labels = scales::label_number()) +
  theme(
    axis.text.x = element_text(size = 16, color = 'black'),
    axis.text.y = element_text(size = 14, color = 'black'),
    axis.title.y = element_text(size = 16, color = 'black'),
    panel.border = element_rect(color = 'black', fill = NA),
    legend.position = 'inside',
    legend.position.inside = c(0.9, 0.2),
    legend.text = element_text(size = 16, color = 'black'),
    plot.title = element_text(size = 22, hjust = 0.5, color = 'black'),
    panel.background = element_blank()
  ) +
  labs(x = NULL, y = expression("Mean pairwise value (log"[10]*")"), fill = NULL, title = "NHRs")


# ============================================================================================================================================================================================ #
# For cytochrome P450s
# ============================================================================================================================================================================================ #
ws_cyto <- readr::read_tsv("../../processed_data/genome_resources/annotation/140_wild_strains_IPR_cytochromeP450.tsv")
sc_cyto <- ws_nhr %>% dplyr::left_join(all_ws_genes_class_og, by = c("strain","gene")) %>% dplyr::filter(Orthogroup %in% og_id_list) %>%
  dplyr::select(OG = Orthogroup) %>% dplyr::distinct() %>% dplyr::pull()

kaks_plt_cyto <- kaks_plt %>% dplyr::filter(OG %in% sc_cyto)

# Wilcoxon tests for each statistic
wilcox_results_cyto <- kaks_plt_cyto %>%
  dplyr::group_by(stats) %>%
  rstatix::wilcox_test(values ~ inHDR) %>%
  rstatix::adjust_pvalue(method = "BH") %>%  
  rstatix::add_significance()

effect_sizes_cyto <- kaks_plt_cyto %>%
  dplyr::group_by(stats) %>%
  rstatix::wilcox_effsize(values ~ inHDR) %>%
  dplyr::select(-.y., -group1, -group2, -n1, -n2)

stats_final_cyto <- left_join(wilcox_results_cyto, effect_sizes_cyto, by = "stats")

annotation_df_cyto <- stats_final_cyto %>%
  dplyr::mutate(label = case_when(
    p.adj < 0.0001 ~ "****",
    p.adj < 0.001 ~ "***",
    p.adj < 0.01 ~ "**",
    p.adj < 0.05 ~ "*",
    TRUE ~ "ns"), y_pos = 6)

# Difference in Ks, Ka, and Ka / Ks between SC OGs in HDRs and not in HDRs?
ggplot(kaks_plt_cyto, aes(x = stats, y = values, fill = inHDR)) +
  geom_boxplot(outlier.size = 0.6, width = 0.7, position = position_dodge(width = 0.75), outliers = FALSE, alpha = 0.5) +
  geom_point(position = position_jitterdodge(jitter.width = 0.25, dodge.width = 0.7), size = 1, alpha = 0.5) +
  geom_text(data = annotation_df_cyto, aes(x = stats, y = y_pos, label = label), inherit.aes = FALSE, size = 5) +
  scale_fill_manual(values = c("non-HDR" = "blue", "HDR" = "RED")) +
  scale_y_log10(labels = scales::label_number()) +
  theme(
    axis.text.x = element_text(size = 16, color = 'black'),
    axis.text.y = element_text(size = 14, color = 'black'),
    axis.title.y = element_text(size = 16, color = 'black'),
    panel.border = element_rect(color = 'black', fill = NA),
    legend.position = 'inside',
    legend.position.inside = c(0.9, 0.2),
    legend.text = element_text(size = 16, color = 'black'),
    plot.title = element_text(size = 22, hjust = 0.5, color = 'black'),
    panel.background = element_blank()
  ) +
  labs(x = NULL, y = expression("Mean pairwise value (log"[10]*")"), fill = NULL, title = "Cytochrome P450s")



# ============================================================================================================================================================================================ #
# Plotting Ka and Ks for all enriched gene classes (no log scale y-axis) 
# ============================================================================================================================================================================================ #
### Ks ###
master_ks <- kaks_plt_cyto %>% dplyr::mutate(class = 'Cytochrome P450s') %>% dplyr::bind_rows(kaks_plt_fbox %>% dplyr::mutate(class = 'F-box genes'), 
                                                                                                kaks_plt_nhr %>% dplyr::mutate(class = 'NHRs'), 
                                                                                                kaks_plt_gpcr %>% dplyr::mutate(class = 'GPCRs'), 
                                                                                                kaks_plt_lectin %>% dplyr::mutate(class = 'C-type lectins')) %>%
  dplyr::filter(stats == "Ks")

master_stats_ks <- annotation_df_cyto %>% dplyr::mutate(class = 'Cytochrome P450s') %>% dplyr::bind_rows(annotation_df_fbox %>% dplyr::mutate(class = 'F-box genes'), 
                                                                                                      annotation_df_nhr %>% dplyr::mutate(class = 'NHRs'), 
                                                                                                      annotation_df_gpcr %>% dplyr::mutate(class = 'GPCRs'), 
                                                                                                      annotation_df_lectin %>% dplyr::mutate(class = 'C-type lectins')) %>%
  dplyr::filter(stats == "Ks") %>%
  dplyr::mutate(y_pos = 0.5)

ggplot(master_ks, aes(x = stats, y = values, fill = inHDR)) +
  geom_boxplot(outlier.size = 0.6, width = 0.7, position = position_dodge(width = 0.75), outliers = FALSE, alpha = 0.5) +
  geom_point(position = position_jitterdodge(jitter.width = 0.25, dodge.width = 0.7), size = 1, alpha = 0.5) +
  geom_text(data = master_stats_ks, aes(x = stats, y = y_pos, label = label), inherit.aes = FALSE, size = 5) +
  scale_fill_manual(values = c("non-HDR" = "blue", "HDR" = "RED")) +
  facet_wrap(~class, nrow = 1) +
  theme(
    axis.text.x = element_text(size = 20, color = 'black'),
    axis.text.y = element_text(size = 18, color = 'black'),
    axis.title.y = element_text(size = 20, color = 'black'),
    panel.border = element_rect(color = 'black', fill = NA),
    strip.text = element_text(size = 22, color = 'black'),
    legend.position = 'inside',
    legend.position.inside = c(0.9, 0.5),
    legend.text = element_text(size = 24, color = 'black'),
    panel.background = element_blank()
  ) +
  labs(x = NULL, y = "Mean pairwise value", fill = NULL)

### Ka ###
master_ka <- kaks_plt_cyto %>% dplyr::mutate(class = 'Cytochrome P450s') %>% dplyr::bind_rows(kaks_plt_fbox %>% dplyr::mutate(class = 'F-box genes'), 
                                                                                              kaks_plt_nhr %>% dplyr::mutate(class = 'NHRs'), 
                                                                                              kaks_plt_gpcr %>% dplyr::mutate(class = 'GPCRs'), 
                                                                                              kaks_plt_lectin %>% dplyr::mutate(class = 'C-type lectins')) %>%
  dplyr::filter(stats == "Ka")

master_stats_ka <- annotation_df_cyto %>% dplyr::mutate(class = 'Cytochrome P450s') %>% dplyr::bind_rows(annotation_df_fbox %>% dplyr::mutate(class = 'F-box genes'), 
                                                                                                         annotation_df_nhr %>% dplyr::mutate(class = 'NHRs'), 
                                                                                                         annotation_df_gpcr %>% dplyr::mutate(class = 'GPCRs'), 
                                                                                                         annotation_df_lectin %>% dplyr::mutate(class = 'C-type lectins')) %>%
  dplyr::filter(stats == "Ka") %>%
  dplyr::mutate(y_pos = 0.33)

ggplot(master_ka, aes(x = stats, y = values, fill = inHDR)) +
  geom_boxplot(outlier.size = 0.6, width = 0.7, position = position_dodge(width = 0.75), outliers = FALSE, alpha = 0.5) +
  geom_point(position = position_jitterdodge(jitter.width = 0.25, dodge.width = 0.7), size = 1, alpha = 0.5) +
  geom_text(data = master_stats_ka, aes(x = stats, y = y_pos, label = label), inherit.aes = FALSE, size = 5) +
  scale_fill_manual(values = c("non-HDR" = "blue", "HDR" = "RED")) +
  facet_wrap(~class, nrow = 1) +
  theme(
    axis.text.x = element_text(size = 20, color = 'black'),
    axis.text.y = element_text(size = 18, color = 'black'),
    axis.title.y = element_text(size = 20, color = 'black'),
    panel.border = element_rect(color = 'black', fill = NA),
    strip.text = element_text(size = 22, color = 'black'),
    legend.position = 'inside',
    legend.position.inside = c(0.9, 0.5),
    legend.text = element_text(size = 24, color = 'black'),
    panel.background = element_blank()
  ) +
  labs(x = NULL, y = "Mean pairwise value", fill = NULL)






# ============================================================================================================================================================================================ #
# Plotting Ka / Ks for all enriched gene classes (no log scale y-axis)
# ============================================================================================================================================================================================ #
master_kaks <- kaks_plt_cyto %>% dplyr::mutate(class = 'Cytochrome P450s') %>% dplyr::bind_rows(kaks_plt_fbox %>% dplyr::mutate(class = 'F-box genes'), 
                                                                                         kaks_plt_nhr %>% dplyr::mutate(class = 'NHRs'), 
                                                                                         kaks_plt_gpcr %>% dplyr::mutate(class = 'GPCRs'), 
                                                                                         kaks_plt_lectin %>% dplyr::mutate(class = 'C-type lectins')) %>%
  dplyr::filter(stats == "Ka / Ks")

master_stats_kaks <- annotation_df_cyto %>% dplyr::mutate(class = 'Cytochrome P450s') %>% dplyr::bind_rows(annotation_df_fbox %>% dplyr::mutate(class = 'F-box genes'), 
                                                                                           annotation_df_nhr %>% dplyr::mutate(class = 'NHRs'), 
                                                                                           annotation_df_gpcr %>% dplyr::mutate(class = 'GPCRs'), 
                                                                                           annotation_df_lectin %>% dplyr::mutate(class = 'C-type lectins')) %>%
  dplyr::filter(stats == "Ka / Ks")

ggplot(master_kaks, aes(x = stats, y = values, fill = inHDR)) +
  geom_boxplot(outlier.size = 0.6, width = 0.7, position = position_dodge(width = 0.75), outliers = FALSE, alpha = 0.5) +
  geom_point(position = position_jitterdodge(jitter.width = 0.25, dodge.width = 0.7), size = 1, alpha = 0.5) +
  geom_text(data = master_stats_kaks, aes(x = stats, y = y_pos, label = label), inherit.aes = FALSE, size = 5) +
  scale_fill_manual(values = c("non-HDR" = "blue", "HDR" = "RED")) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  facet_wrap(~class, nrow = 1) +
  theme(
    axis.text.x = element_text(size = 20, color = 'black'),
    axis.text.y = element_text(size = 18, color = 'black'),
    axis.title.y = element_text(size = 20, color = 'black'),
    panel.border = element_rect(color = 'black', fill = NA),
    strip.text = element_text(size = 22, color = 'black'),
    legend.position = 'inside',
    legend.position.inside = c(0.9, 0.5),
    legend.text = element_text(size = 24, color = 'black'),
    panel.background = element_blank()
  ) +
  labs(x = NULL, y = "Mean pairwise value", fill = NULL)







# ============================================================================================================================================================================================ #
# What proportion of single-copy OGs in HDRs are these five gene classes?
# ============================================================================================================================================================================================ #
five_geneClasses <- kaks_plt_gpcr %>% dplyr::bind_rows(kaks_plt_cyto, kaks_plt_lectin, kaks_plt_fbox, kaks_plt_nhr) %>%
  dplyr::filter(inHDR == "HDR") %>% dplyr::distinct(OG) # 394.....













