library(dplyr)
library(ggplot2)
library(seqinr)
library(tidyr)

aln <- read.alignment("../../processed_data/divergence_estimates/codon_alignments/OG0003089.aligned.codons.gapped.fa", format = "fasta")

# Calculate Ka/Ks
res <- seqinr::kaks(aln, rmgap = TRUE)

# View results
print(res)

ka_mat  <- as.matrix(res$ka)
ks_mat  <- as.matrix(res$ks)
vka_mat <- as.matrix(res$vka)
vks_mat <- as.matrix(res$vks)

make_pair_table <- function(mat, value_name) {
  df <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
  colnames(df) <- c("seq1", "seq2", value_name)
  df$seq1 <- as.character(df$seq1)
  df$seq2 <- as.character(df$seq2)
  df
}

ka_table  <- make_pair_table(ka_mat, "Ka")
ks_table  <- make_pair_table(ks_mat, "Ks")
vka_table <- make_pair_table(vka_mat, "vKa")
vks_table <- make_pair_table(vks_mat, "vKs")

final_table <- merge(ka_table, ks_table, by = c("seq1", "seq2"))
final_table <- merge(final_table, vka_table, by = c("seq1", "seq2"))
final_table <- merge(final_table, vks_table, by = c("seq1", "seq2"))

final_table <- final_table[final_table$seq1 < final_table$seq2, ] %>% dplyr::mutate(kaks = dplyr::if_else(Ks > 0, Ka / Ks, NA_real_)) %>% dplyr::select(-vKs, -vKa) # %>%
  # ADD OG ID, and then left-merge to add HDR resolution and the number of sc orthologs found in a HDR for a wild strain

summary <- final_table %>% dplyr::summarize(OG = "OGXX",
                                            HDR = TRUE,
                                            mean_pairwise_ks = mean(Ks),
                                            mean_pairwise_ka = mean(Ka),
                                            mean_pairwise_kaks = mean(kaks, na.rm = T)) %>%
  tidyr::pivot_longer(c(-OG,-HDR), names_to = 'stats', values_to = "values") %>%
  dplyr::mutate(stats = factor(stats, levels = c("mean_pairwise_ks", "mean_pairwise_ka", "mean_pairwise_kaks")))




ggplot(summary) +
  geom_boxplot(aes(x = stats, y = values))













# ============================================================================================================================================================================================ #
# For all 13,315 orthogroups
# ============================================================================================================================================================================================ #

# Setting variables and reading in OG list:
aln_dir <- "../../processed_data/divergence_estimates/codon_alignments"
suffix <- ".aligned.codons.gapped.fa"
og_id_list <- readr::read_tsv("../../processed_data/orthology/orthofinder/orthofinder_output/Orthogroups_SingleCopyOrthologues.txt", col_names = "OGs") %>% dplyr::pull(OGs)

# Function to read in file and filter for occupancy
read_and_filter_alignment <- function(aln_file) {
  aln <- seqinr::read.alignment(aln_file, format = "fasta")
  
  seqs <- toupper(aln$seq)
  
  occupancy <- vapply(seqs, function(x) {
    chars <- strsplit(x, "", fixed = TRUE)[[1]]
    sum(chars %in% c("A", "C", "G", "T")) / length(chars)
  }, numeric(1))
  
  lenest <- vapply(seqs, function(x) {
    chars <- strsplit(x, "", fixed = TRUE)[[1]]
    sum(chars %in% c("A", "C", "G", "T")) / 3
  }, numeric(1))
  
  keep <- occupancy >= 0.7 & lenest >= 100
  
  if (sum(keep) < 4) return(NULL)
  
  aln$nam <- aln$nam[keep]
  aln$seq <- aln$seq[keep]
  aln$nb <- length(aln$seq)
  
  return(aln)
}

# Function to calculate Ka and Ks for each OG and create a nice summary table
get_kaks <- function(og_id) {
  
  aln_file <- file.path(aln_dir, paste0(og_id, suffix))
  # aln_file <- aln_file[file.exists(aln_file)][1]
  
  if (is.na(aln_file)) {
    warning("Skipping ", og_id, ": no alignment found")
    return(NULL) }
  
  aln <- read_and_filter_alignment(aln_file)
  
  if (is.null(aln)) {
    warning("Skipping ", og_id, ": fewer than 4 sequences after filtering")
    return(NULL) }
  
  res <- seqinr::kaks(aln, rmgap = TRUE)
  
  # View results
  print(res)
  
  ka_mat  <- as.matrix(res$ka)
  ks_mat  <- as.matrix(res$ks)
  vka_mat <- as.matrix(res$vka)
  vks_mat <- as.matrix(res$vks)
  
  make_pair_table <- function(mat, value_name) {
    df <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
    colnames(df) <- c("seq1", "seq2", value_name)
    df$seq1 <- as.character(df$seq1)
    df$seq2 <- as.character(df$seq2)
    df
  }
  
  ka_table  <- make_pair_table(ka_mat, "Ka")
  ks_table  <- make_pair_table(ks_mat, "Ks")
  vka_table <- make_pair_table(vka_mat, "vKa")
  vks_table <- make_pair_table(vks_mat, "vKs")
  
  final_table <- merge(ka_table, ks_table, by = c("seq1", "seq2"))
  final_table <- merge(final_table, vka_table, by = c("seq1", "seq2"))
  final_table <- merge(final_table, vks_table, by = c("seq1", "seq2"))
  
  final_table <- final_table[final_table$seq1 < final_table$seq2, ] %>% dplyr::mutate(kaks = dplyr::if_else(Ks > 0, Ka / Ks, NA_real_)) %>% dplyr::select(-vKs, -vKa) #%>%
    # ADD OG ID, and then left-merge to add HDR resolution and the number of sc orthologs found in a HDR for a wild strain
    
  summary <- final_table %>% dplyr::summarize(OG = og_id,
                                              HDR = TRUE,
                                              mean_pairwise_ks = mean(Ks),
                                              mean_pairwise_ka = mean(Ka),
                                              mean_pairwise_kaks = mean(kaks, na.rm = T)) %>%
    tidyr::pivot_longer(c(-OG,-HDR), names_to = 'stats', values_to = "values") %>%
    dplyr::mutate(stats = factor(stats, levels = c("mean_pairwise_ks", "mean_pairwise_ka", "mean_pairwise_kaks")))
  
  return(summary)
}


# Iterating over all 12,315 OGs to create a master summary file
final_kaks_df <- data.frame()
for (i in og_id_list) {
  kaks_stats_i <- get_kaks(i)
  final_kaks_df <- final_kaks_df %>% dplyr::bind_row(kaks_stats_i)
}














