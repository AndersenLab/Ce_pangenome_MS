library(dplyr)
library(ggplot2)
library(seqinr)
library(tidyr)

aln <- read.alignment("../../processed_data/divergence_estimates/codon_alignments/OG0003089.aligned.codons.gapped.fa", format = "fasta")

# Calculate Ka/Ks
res <- seqinr::kaks(aln)

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
