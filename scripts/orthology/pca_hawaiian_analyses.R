library(ggplot2)
library(dplyr)
library(ape)
library(readr)
library(tidyr)


# Load in OG matrix
og_mat <- readr::read_tsv("../../tables/all_OGs_matrix.tsv")
og_mat_count <- readr::read_tsv("../../processed_data/orthology/OG_relations_matrix_count.tsv")
all_ogs_class <- readr::read_tsv("../../processed_data/orthology/all_genes_class_OGs.tsv") %>% dplyr::distinct(Orthogroup, class)


# Load in Hawaiian strains that cluster in PCA, BUSCO-inferred relatendess, and accessory orthogroup clustering
pca_strain <- readr::read_tsv("../../processed_data/structural_variants/PCA_hawaii_cluster.tsv") %>% dplyr::pull()


# Load in SV calls
allcalls <- readr::read_tsv("../../processed_data/structural_variants/141_over50_PASS_variants.tsv", col_names = c("chrom", "pos", "ref", "alt", "filter", "sv_type","sv_length","strain")) %>% dplyr::select(-filter) %>% 
  dplyr::filter(chrom != "MtDNA")

filt_calls <- allcalls %>% 
  dplyr::mutate(sv_length = abs(sv_length)) %>%
  dplyr::group_by(sv_type, strain) %>%
  dplyr::mutate(type_count = n()) %>%
  dplyr::ungroup() %>% 
  dplyr::filter(strain %in% pca_strain)

merged_SV <- readr::read_tsv("../../processed_data/structural_variants/Jasmine_merged_SVs.tsv") %>%
  dplyr::mutate(sv_length = abs(sv_length)) %>%
  dplyr::select(1:7, dplyr::all_of(pca_strain)) %>%
  dplyr::filter(if_all(8:last_col(), ~ .x != "./.")) %>%
  dplyr::filter(number_svs_merged == 19)


ggplot(merged_SV) +
  geom_rect(aes(xmin = (pos - 150) / 1e6, xmax = (pos + sv_length + 150) / 1e6, ymin = 0, ymax = 1, fill = sv_type)) +
  facet_wrap(~chrom, scales = "free") +
  scale_fill_manual(values = c("INS" = "blue", "DEL" = "red", "INV" = "gold")) +
  theme(
    panel.background = element_blank(),
    axis.text.x = element_text(size = 16, color = 'black'),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = 'none',
    axis.title.y  = element_blank(),
    strip.text = element_text(size = 24, color = 'black'),
    axis.title.x = element_text(size = 20, color = 'black'),
    panel.border = element_rect(color = 'black', fill = NA)
  ) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "N2 genomic position (Mb)", fill = NULL)
 
# Load in N2 genes
n2_names <- readr::read_tsv("../../processed_data/genome_resources/annotation/N2_geneName_seqName.tsv", col_names = c("name","gene"))

n2_genes <- ape::read.gff("../../processed_data/genome_resources/annotation/c_elegans.PRJNA13758.WS283.csq.PCfeaturesOnly.longest.gff3") %>%
  dplyr::filter(type == "gene", 
                seqid == "I" | seqid == "II" | seqid == "III" | seqid == "IV") %>%
  tidyr::separate(attributes, into = c("junk", "gene"), sep = ';') %>%
  dplyr::mutate(gene = gsub("Alias=","",gene)) %>%
  dplyr::select(-junk,-phase,-score) %>%
  tidyr::separate(gene, into = c("name","seqName"), sep = ",") %>%
  dplyr::select(-seqName) %>%
  dplyr::mutate(name_position = ((start + end) / 2) / 1e6)



n2_exons <- ape::read.gff("../../processed_data/genome_resources/annotation/c_elegans.PRJNA13758.WS283.csq.PCfeaturesOnly.longest.gff3") %>%
  dplyr::filter(type == "exon", 
                seqid == "I" | seqid == "II" | seqid == "III" | seqid == "IV") %>%
  tidyr::separate(attributes, into = c("ID","gene"), sep = ";") %>%
  dplyr::mutate(gene = gsub("Parent=transcript:","", gene)) %>%
  dplyr::mutate(gene = sub("\\.[^.]*$", "", gene)) %>%  # removing trailing isoform numbers
  dplyr::mutate(gene = sub("[A-Za-z]+$", "", gene)) %>% # removing trailing letters
  dplyr::left_join(n2_names, by = "gene") %>%
  dplyr::select(-ID,-gene, -score, -phase)




# Left SV on chromosome I
ggplot(merged_SV %>% dplyr::filter(chrom == "I")) +
  geom_rect(aes(xmin = (pos - 150) / 1e6, xmax = (pos + sv_length + 150) / 1e6, ymin = 0, ymax = 1, fill = sv_type)) +
  geom_rect(data = n2_genes %>% dplyr::filter(name == "tol-1"), aes(xmin = start / 1e6, xmax = end / 1e6, ymin = 0.52, ymax = 0.54), fill = "black", color = "black") +
  geom_rect(data = n2_exons %>%  dplyr::filter(name == "tol-1"), aes(xmin = start / 1e6, xmax = end / 1e6, ymin = 0.4, ymax = 0.5), fill = "grey", color = "black") +
  geom_text(data = n2_genes %>% dplyr::filter(name == "tol-1"), aes(x = name_position, y = 0.55, label = name), hjust = 0, size = 4, inherit.aes = FALSE) +
  facet_wrap(~chrom, scales = "free") +
  scale_fill_manual(values = c("INS" = "blue", "DEL" = "red", "INV" = "gold")) +
  theme(
    panel.background = element_blank(),
    axis.text.x = element_text(size = 16, color = 'black'),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = 'none',
    axis.title.y  = element_blank(),
    strip.text = element_text(size = 24, color = 'black'),
    axis.title.x = element_text(size = 20, color = 'black'),
    panel.border = element_rect(color = 'black', fill = NA)
  ) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "N2 genomic position (Mb)", fill = NULL) +
  coord_cartesian(xlim = c(0.453424 - 0.025, 0.453424 + 0.055))
  
  
# Right SV on chromosome I
ggplot(merged_SV %>% dplyr::filter(chrom == "I")) +
    geom_rect(aes(xmin = (pos - 150) / 1e6, xmax = (pos + sv_length + 150) / 1e6, ymin = 0, ymax = 1, fill = sv_type)) +
    geom_rect(data = n2_genes, aes(xmin = start / 1e6, xmax = end / 1e6, ymin = 0.52, ymax = 0.54), fill = "black", color = "black") +
    geom_rect(data = n2_exons, aes(xmin = start / 1e6, xmax = end / 1e6, ymin = 0.4, ymax = 0.5), fill = "grey", color = "black") +
    geom_text(data = n2_genes %>% dplyr::filter(name_position > 0.72 & name_position < 0.78, seqid == "I"), aes(x = name_position, y = 0.55, label = name), hjust = 0, size = 4, inherit.aes = FALSE) +
    facet_wrap(~chrom, scales = "free") +
    scale_fill_manual(values = c("INS" = "blue", "DEL" = "red", "INV" = "gold")) +
    theme(
      panel.background = element_blank(),
      axis.text.x = element_text(size = 16, color = 'black'),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      legend.position = 'none',
      axis.title.y  = element_blank(),
      strip.text = element_text(size = 24, color = 'black'),
      axis.title.x = element_text(size = 20, color = 'black'),
      panel.border = element_rect(color = 'black', fill = NA)
    ) +
    scale_y_continuous(expand = c(0,0)) +
    labs(x = "N2 genomic position (Mb)", fill = NULL) +
    coord_cartesian(xlim = c(0.728330 - 0.043, 0.728330 + 0.08))


# Largest SV on II
ggplot(merged_SV %>% dplyr::filter(chrom == "II")) +
  geom_rect(aes(xmin = (pos - 150) / 1e6, xmax = (pos + sv_length + 150) / 1e6, ymin = 0, ymax = 1, fill = sv_type)) +
  geom_rect(data = n2_genes %>% dplyr::filter(name == "sra-2" | name == "sra-9"), aes(xmin = start / 1e6, xmax = end / 1e6, ymin = 0.52, ymax = 0.54), fill = "black", color = "black") +
  geom_rect(data = n2_exons %>% dplyr::filter(name == "sra-2" | name == "sra-9"), aes(xmin = start / 1e6, xmax = end / 1e6, ymin = 0.4, ymax = 0.5), fill = "grey", color = "black") +
  geom_text(data = n2_genes %>% dplyr::filter(name == "sra-2" | name == "sra-9") %>% dplyr::filter(name_position > 9.52 & name_position < 9.58, seqid == "II"), aes(x = name_position, y = 0.57, label = name), hjust = 0.5, size = 8, inherit.aes = FALSE) +
  facet_wrap(~chrom, scales = "free") +
  scale_fill_manual(values = c("INS" = "blue", "DEL" = "red", "INV" = "gold")) +
  theme(
    panel.background = element_blank(),
    axis.text.x = element_text(size = 16, color = 'black'),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = 'none',
    axis.title.y  = element_blank(),
    strip.text = element_text(size = 24, color = 'black'),
    axis.title.x = element_text(size = 20, color = 'black'),
    panel.border = element_rect(color = 'black', fill = NA)
  ) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "N2 genomic position (Mb)", fill = NULL) +
  coord_cartesian(xlim = c(9.527080 - 0.003, 9.527080 + 0.009))


og_mat_del_II <- og_mat %>% dplyr::filter(grepl("AH6.14", N2) | grepl("AH6.6", N2)) %>% dplyr::select(N2, dplyr::all_of(pca_strain))

# Lack of alignment???
nucmer <- readr::read_tsv("../../processed_data/genome_resources/genome_data/141_nucmer_ECA741CGC1.tsv", col_names = c("N2S","N2E","WSS","WSE","L1","L2","IDY","LENR","LENQ","N2_chr","contig","strain")) %>%
  dplyr::filter(strain == "ECA730" | strain == "ECA1202" | strain == "ECA722" | strain == "ECA723" | strain == "ECA2948" | strain == "ECA1409", N2_chr == "II") %>%
  dplyr::mutate(inv = ifelse(WSE < WSS, TRUE, FALSE))

sra_2_9_locus <- ggplot2::ggplot(nucmer) +
  geom_rect(data = merged_SV, aes(xmin = (pos - 150) / 1e6, xmax = (pos + sv_length + 150) / 1e6, ymin = -Inf, ymax = Inf), fill = "red", alpha = 0.5) +
  ggplot2::geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6), linewidth = 0.7) +
  geom_rect(data = n2_genes %>% dplyr::filter(name == "sra-2" | name == "sra-9"), aes(xmin = start / 1e6, xmax = end / 1e6, ymin = -Inf, ymax = Inf), fill = NA, color = "black", linetype = "dashed") +
  geom_text(data = n2_genes %>% dplyr::filter(name == "sra-2" | name == "sra-9") %>% dplyr::filter(name_position > 9.52 & name_position < 9.58, seqid == "II"), aes(x = name_position, y = 10, label = name), hjust = 0.5, size = 6, inherit.aes = FALSE) +
  facet_wrap(~chrom, scales = "free") +
  ggplot2::theme_bw() +
  ggplot2::facet_wrap(~strain, nrow = 3) +
  ggplot2::theme(
    legend.position = "none",
    strip.text = ggplot2::element_text(size = 14, color = 'black'),
    axis.text = ggplot2::element_text(size = 10, color = 'black'),
    axis.ticks = ggplot2::element_blank(),
    axis.title = ggplot2::element_text(size = 12, color = 'black'),
    panel.background = ggplot2::element_blank(),
    panel.grid = ggplot2::element_blank(),
    panel.border = ggplot2::element_rect(fill = NA)) +
  ggplot2::labs(x = "N2 chromosome II position (Mb)", y = "Wild strain contig position (Mb)") +
  scale_y_continuous(expand = c(0.008,0.008)) +
  scale_x_continuous(expand = c(0.008,0.008)) +
  coord_cartesian(xlim = c(9.527080 - 0.003, 9.527080 + 0.009))
sra_2_9_locus


# Largest and only INV on II 
ggplot(merged_SV %>% dplyr::filter(chrom == "II", pos == 5346585)) +
  geom_rect(aes(xmin = (pos - 150) / 1e6, xmax = (pos + sv_length + 150) / 1e6, ymin = 0, ymax = 1, fill = sv_type)) +
  geom_rect(data = n2_genes, aes(xmin = start / 1e6, xmax = end / 1e6, ymin = 0.52, ymax = 0.54), fill = "black", color = "black") +
  geom_rect(data = n2_exons, aes(xmin = start / 1e6, xmax = end / 1e6, ymin = 0.4, ymax = 0.5), fill = "grey", color = "black") +
  geom_text(data = n2_genes %>% dplyr::filter(name_position > 5.345 & name_position < 5.355, seqid == "II"), aes(x = name_position, y = 0.55, label = name), hjust = 0, size = 4, inherit.aes = FALSE) +
  facet_wrap(~chrom, scales = "free") +
  scale_fill_manual(values = c("INS" = "blue", "DEL" = "red", "INV" = "gold")) +
  theme(
    panel.background = element_blank(),
    axis.text.x = element_text(size = 16, color = 'black'),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = 'none',
    axis.title.y  = element_blank(),
    strip.text = element_text(size = 24, color = 'black'),
    axis.title.x = element_text(size = 20, color = 'black'),
    panel.border = element_rect(color = 'black', fill = NA)
  ) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "N2 genomic position (Mb)", fill = NULL) +
  coord_cartesian(xlim = c(5.346585 - 0.01, 5.346585 + 0.01))


# Largest INS on II
ggplot(merged_SV %>% dplyr::filter(chrom == "II", pos == 7305372)) +
  geom_rect(aes(xmin = (pos - 150) / 1e6, xmax = (pos + sv_length + 150) / 1e6, ymin = 0, ymax = 1, fill = sv_type)) +
  geom_rect(data = n2_genes, aes(xmin = start / 1e6, xmax = end / 1e6, ymin = 0.52, ymax = 0.54), fill = "black", color = "black") +
  geom_rect(data = n2_exons, aes(xmin = start / 1e6, xmax = end / 1e6, ymin = 0.4, ymax = 0.5), fill = "grey", color = "black") +
  geom_text(data = n2_genes %>% dplyr::filter(name_position > 7.3 & name_position < 7.310, seqid == "II"), aes(x = name_position, y = 0.55, label = name), hjust = 0, size = 4, inherit.aes = FALSE) +
  facet_wrap(~chrom, scales = "free") +
  scale_fill_manual(values = c("INS" = "blue", "DEL" = "red", "INV" = "gold")) +
  theme(
    panel.background = element_blank(),
    axis.text.x = element_text(size = 16, color = 'black'),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = 'none',
    axis.title.y  = element_blank(),
    strip.text = element_text(size = 24, color = 'black'),
    axis.title.x = element_text(size = 20, color = 'black'),
    panel.border = element_rect(color = 'black', fill = NA)
  ) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "N2 genomic position (Mb)", fill = NULL) +
  coord_cartesian(xlim = c(7.305372 - 0.01, 7.305372 + 0.01))

og_mat_ins_II <- og_mat %>% dplyr::filter(grepl("F45E12.3", N2)) %>% dplyr::select(N2, dplyr::all_of(pca_strain))



# Only one on IV
ggplot(merged_SV %>% dplyr::filter(chrom == "IV", pos == 15689620)) +
  geom_rect(aes(xmin = (pos - 150) / 1e6, xmax = (pos + sv_length + 150) / 1e6, ymin = 0, ymax = 1, fill = sv_type)) +
  geom_rect(data = n2_genes, aes(xmin = start / 1e6, xmax = end / 1e6, ymin = 0.52, ymax = 0.54), fill = "black", color = "black") +
  geom_rect(data = n2_exons, aes(xmin = start / 1e6, xmax = end / 1e6, ymin = 0.4, ymax = 0.5), fill = "grey", color = "black") +
  geom_text(data = n2_genes %>% dplyr::filter(name_position > 15.68 & name_position < 15.725, seqid == "IV"), aes(x = name_position, y = 0.55, label = name), hjust = 0, size = 4, inherit.aes = FALSE) +
  facet_wrap(~chrom, scales = "free") +
  scale_fill_manual(values = c("INS" = "blue", "DEL" = "red", "INV" = "gold")) +
  theme(
    panel.background = element_blank(),
    axis.text.x = element_text(size = 16, color = 'black'),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = 'none',
    axis.title.y  = element_blank(),
    strip.text = element_text(size = 24, color = 'black'),
    axis.title.x = element_text(size = 20, color = 'black'),
    panel.border = element_rect(color = 'black', fill = NA)
  ) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "N2 genomic position (Mb)", fill = NULL) +
  coord_cartesian(xlim = c(15.689620 - 0.01, 15.689620 + 0.05))



# Looking for PAVs specific to these strains
pav_find <- og_mat_count %>% dplyr::select(Orthogroup, N2, dplyr::all_of(pca_strain)) %>%
  dplyr::filter(!if_all(-Orthogroup, is.na)) %>%
  dplyr::left_join(all_ogs_class, by = "Orthogroup") %>%
  dplyr::filter(class != "private") %>%
  dplyr::filter((is.na(N2) & if_all(-c(Orthogroup, N2, class), ~ !is.na(.))) | (!is.na(N2) & if_all(-c(Orthogroup, N2, class), is.na))) %>%
  dplyr::mutate(pav_class = ifelse(
    (is.na(N2) & if_all(-c(Orthogroup, N2, class), ~ !is.na(.))), "WS_presence", "N2_presence"))
pav_find_ogs <- pav_find %>% dplyr::pull(Orthogroup)

pres_n2 <- pav_find %>% dplyr::filter(pav_class == "N2_presence") %>% dplyr::select(Orthogroup, N2) %>% dplyr::pull(Orthogroup)

# The four genes that are only absent in the 19 Hawaiian strains
only_absent_19_hawaiian <- og_mat %>% dplyr::filter(Orthogroup %in% pres_n2) %>%
  dplyr::filter(!is.na(N2), if_all(dplyr::all_of(pca_strain), ~ is.na(.)), if_all(-c(Orthogroup, N2, dplyr::all_of(pca_strain)), ~ !is.na(.))) %>%
  dplyr::select(N2) %>% dplyr::mutate(N2 = gsub("transcript_","", N2)) %>%
  dplyr::mutate(N2 = sub("\\.[^.]*$", "", N2)) %>%  # removing trailing isoform numbers
  dplyr::mutate(N2 = sub("[A-Za-z]+$", "", N2)) %>% # removing trailing letters
  dplyr::left_join(n2_names, by = c("N2" = "gene")) 


pres_wss <- pav_find %>% dplyr::filter(pav_class == "WS_presence") %>% dplyr::pull(Orthogroup)

only_present_19_hawaiian <- og_mat %>% dplyr::filter(Orthogroup %in% pres_wss) %>%
  dplyr::filter(is.na(N2), if_all(dplyr::all_of(pca_strain), ~ !is.na(.)), if_all(-c(Orthogroup, N2, dplyr::all_of(pca_strain)), ~ is.na(.))) # only TWO genes.. might be a kinesin protein and SLC12A transporter



pres_n2_genes <- og_mat %>% dplyr::filter(Orthogroup %in% pres_n2) %>% dplyr::select(Orthogroup, N2)# %>%
  tidyr::separate_rows(N2, sep = ", ") %>% dplyr::mutate(N2 = gsub("transcript_","", N2)) %>%
  dplyr::mutate(N2 = sub("\\.[^.]*$", "", N2)) %>%  # removing trailing isoform numbers
  dplyr::mutate(N2 = sub("[A-Za-z]+$", "", N2)) %>% # removing trailing letters
  dplyr::left_join(n2_names, by = c("N2" = "gene")) %>%
  dplyr::arrange(N2) %>%
  dplyr::select(-N2, N2_genes = name)

write.table(pres_n2_genes, "/vast/eande106/projects/Lance/THESIS_WORK/misc/N2_genes_absent_19_hawaiian.tsv", quote = F, row.names = F, col.names = T, sep = "\t")
  
pres_ws <- pav_find %>% dplyr::filter(pav_class == "WS_presence") %>% dplyr::select(-N2, -Orthogroup, -class, -pav_class) %>%
  dplyr::mutate(og_gene_count = rowSums(across(where(is.numeric)))) %>%
  dplyr::arrange(desc(og_gene_count)) %>%
  dplyr::mutate(average_ortholog_count = og_gene_count / 19)

ggplot(pres_ws %>% dplyr::mutate(PAV_class = "WS_presence")) +
  geom_histogram(aes(x = average_ortholog_count)) +
  theme_bw() +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0))

