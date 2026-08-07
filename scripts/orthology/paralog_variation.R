library(tidyr)
library(dplyr)
library(ggplot2)
library(readr)

# Load in ortholog matrix with counts of orthologs
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
  geom_point(aes(x = occupancy_threshold, y = ortho_count), size = 3.5, color = 'orange') +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(color = 'black', fill = NA),
    axis.text = element_text(size = 14, color = 'black'),
    axis.title = element_text(size = 17, color = 'black')
  ) +
  labs(x = 'Occupancy threshold', y = 'Orthogroup count')
# no clear decrease at a certain occupancy... arbitrarily setting threshold at 50%


ortho_count_filt <- ortho_count %>%
  dplyr::filter(freq >= 0.5) %>%
  dplyr::mutate(
    average_ortholog_count = rowMeans(dplyr::across(2:(ncol(.) - 2), ~ as.numeric(.x)), na.rm = TRUE),
    sd_ortholog_count = apply(dplyr::across(2:(ncol(ortho_count) - 2), ~ as.numeric(.x)), 1, sd, na.rm = TRUE)) %>%
  dplyr::filter(sd_ortholog_count != 0) %>%
  dplyr::arrange(desc(sd_ortholog_count)) %>%
  dplyr::slice_head(n = 50) %>%
  dplyr::select(Orthogroup, sum, freq, average_ortholog_count, sd_ortholog_count)

order <- ortho_count_filt %>% dplyr::pull(Orthogroup)

ggplot(ortho_count_filt %>% dplyr::mutate(Orthogroup = factor(Orthogroup, levels = order)), aes(x = Orthogroup, y = average_ortholog_count)) +
  geom_errorbar(aes(ymin = average_ortholog_count - sd_ortholog_count, ymax = average_ortholog_count + sd_ortholog_count), 
                width = 0.5) +
  geom_point(aes(fill = freq), size = 6, shape = 22) +
  scale_fill_gradient(low = 'yellow', high = 'red') +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(color = 'black', fill = NA),
    axis.text.y = element_text(size = 14, color = 'black'),
    axis.text.x = element_text(size = 14, color = 'black', angle = 60, hjust = 1),
    legend.text = element_text(size = 14, color= 'black'),
    legend.title = element_text(size = 16, color = 'black'),
    axis.title = element_text(size = 17, color = 'black')
  ) +
  labs(x = 'Orthogroup', y = 'Average ortholog count', fill = "Frequency")




# Adding IPR resolution:
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
  dplyr::slice_head(n = 1)



