library(ggplot2)
library(dplyr)
library(readr)

# ========================================================================================================================================================================================================= #
# Load in abamectin and ivermectin trait data
# ========================================================================================================================================================================================================= #
traits <- readr::read_tsv("../../tables/mapping_2018GWAS_MLs.tsv") %>%
  dplyr::select(strain, `Abamectin_q90.TOF_ctrl-regressed`, `Ivermectin_median.TOF_ctrl-regressed`)

geno_matrix <- readr::read_tsv("../../tables/genotype_matrix.tsv")

wild_strains_140 <- readr::read_tsv("../../tables/wild_strain_genome_stats.tsv") %>% dplyr::select(Strain) %>% dplyr::pull()

# ========================================================================================================================================================================================================= #
# Visualizing abamectin haplotype for ALL strains we have assemblies for 
# ========================================================================================================================================================================================================= #
# Looking at REF ALT split
abam_geno <- geno_matrix %>% dplyr::filter(CHROM == "V" & POS == 16198034) %>%
  tidyr::pivot_longer(-c(CHROM, POS, REF, ALT), names_to = "strain", values_to = 'abam_gt') %>%
  dplyr::select(-CHROM, -POS, -REF, -ALT) %>%
  dplyr::mutate(abam_gt = ifelse(abam_gt == "-1", "REF", "ALT"))

# Looking at strains we have assemblies for
abam_traits <- traits %>% dplyr::select(strain,  `Abamectin_q90.TOF_ctrl-regressed`) %>% dplyr::rename(abam_q90_TOF =  `Abamectin_q90.TOF_ctrl-regressed`) %>%
  dplyr::left_join(abam_geno, by = "strain") %>%
  dplyr::filter(!is.na(abam_q90_TOF) & !is.na(abam_gt)) %>%
  dplyr::mutate(we_have_asm = ifelse(strain %in% wild_strains_140, "TRUE", "FALSE")) %>%
  dplyr::filter(we_have_asm == TRUE) %>%
  dplyr::group_by(abam_gt) %>%
  dplyr::arrange(desc(abam_q90_TOF)) %>%
  dplyr::mutate(
    extreme = "YES") %>%
  dplyr::ungroup()

want <- abam_traits %>% dplyr::filter(extreme == "YES") %>% dplyr::arrange(desc(abam_gt), abam_q90_TOF) %>% dplyr::pull(strain) %>% c("N2")

alt_strains <- tail(want,17) 

n2_genes_interval <- readr::read_tsv("../../processed_data/gwas/n2_genes_inInterval.txt", col_names = "gene") %>% dplyr::pull()

n2_gff <- ape::read.gff("../../processed_data/genome_resources/annotation/c_elegans.PRJNA13758.WS283.csq.PCfeaturesOnly.longest.w-mRNA_featuresgff3") %>% dplyr::mutate(STRAIN="N2")


N2_genes <- n2_gff %>%
  dplyr::filter(type=="gene") %>%
  # dplyr::mutate(refStart=N2Start,refEnd=N2End) %>%
  dplyr::filter(grepl("biotype=protein_coding",attributes)) %>%
  tidyr::separate(attributes,into=c("pre","post"),sep=';sequence_name=') %>%
  tidyr::separate(post,into=c("seqname","post2"),sep=';biotype=') %>%
  tidyr::separate(pre,into=c("ID","Name","rest2"),sep=";") %>%
  dplyr::mutate(Name=gsub("Name=","",Name)) %>% 
  # dplyr::select(seqid,start,end,strand,Name,rest2,seqname,STRAIN,refStart,refEnd,seqname) %>%
  dplyr::mutate(rest2=ifelse(grepl("locus",rest2),gsub("locus=","",rest2),seqname)) %>%
  dplyr::rename(alias=rest2)

N2_tran <- n2_gff %>%
  dplyr::filter(type=="mRNA") %>%
  tidyr::separate(attributes, into=c("ID","Parent","Name","wormpep","locus"),sep=';') %>%
  dplyr::mutate(ID=gsub("ID=Transcript:","",ID)) %>%
  tidyr::separate(ID,into = c("fosmid","tseqID",'tranum'),sep="\\.",remove = F) %>%
  dplyr::mutate(tseqname=paste0(fosmid,".",tseqID,".",tranum)) %>%
  dplyr::mutate(Parent=gsub("Parent=Gene:","",Parent)) %>%
  dplyr::select(tseqname,Parent) %>%
  dplyr::rename(tranname=tseqname) %>%
  dplyr::mutate(tranname=paste0("transcript_",tranname)) %>%
  dplyr::left_join(N2_genes, by = c("Parent" = "Name"))

# Tile plot of SNVs
interval_snvs <- readr::read_tsv("../../processed_data/gwas/SNVs_in_Interval.tsv") %>%
  dplyr::filter(POS >= (16198034 - 22000) & POS <= (16198034 + 25000)) %>%
  dplyr::mutate(dplyr::across(4:last_col(), ~ dplyr::case_when(
    .x == "0/0" ~ "0",
    .x == "1/1" ~ "1",
    .x == "./." ~ NA_character_,
    TRUE ~ .x))) %>%
  dplyr::select(-CHROM,-REF,-ALT, -N2) 

interval_snv_nonMissing <- interval_snvs %>% dplyr::filter(dplyr::if_all(4:last_col(), ~ !is.na(.x))) %>%
  dplyr::mutate(alt_specific = dplyr::if_all(
    dplyr::all_of(alt_strains[alt_strains != "N2"]), ~ .x == 1) &
      dplyr::if_all(dplyr::all_of(setdiff(want, alt_strains)), ~ .x == 0))

alt_specific <- interval_snv_nonMissing %>% dplyr::select(POS, alt_specific) %>% dplyr::filter(alt_specific == "TRUE") %>%
  dplyr::mutate(seqid = "V")

pos <- interval_snv_nonMissing$POS

snv_nonMissing_long <- interval_snv_nonMissing %>%
  dplyr::select(-POS, -alt_specific) %>%
  t() %>%
  as.data.frame() %>%
  setNames(pos) %>%
  tibble::rownames_to_column("strain") %>%
  tidyr::pivot_longer(
    cols = -strain,
    names_to = "POS",
    values_to = "genotype") %>%
  dplyr::mutate(POS = as.numeric(POS)) %>%
  dplyr::mutate(strain = factor(strain, levels = rev(want[want != "N2"]))) %>%
  dplyr::mutate(genotype = ifelse(genotype == "0", "REF", "ALT")) 

n2_genes_ROI <- n2_gff %>% dplyr::filter(seqid == "V", start >= (16198034 - 22000) & end <= (16198034 + 25000), type == "gene") 

genes_interval <- ggplot(n2_genes_ROI) +
  geom_hline(yintercept = 1) +
  geom_rect(aes(xmin = start / 1e6, xmax = end / 1e6, ymin = 0.9, ymax = 1.1)) +
  geom_vline(xintercept = 16198034 / 1e6, color = 'red', linetype = 'solid') +
  geom_rect(data = alt_specific, aes(xmin = (POS - 10) / 1e6, xmax = (POS + 10) / 1e6, ymin = 0.8, ymax = 1.2), fill = "magenta3") +
  scale_y_continuous(expand = c(0,0)) +
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  coord_cartesian(xlim = c((16.198034 - 0.022000), (16.198034 + 0.025000)))
genes_interval 

# Adding gene labels
n2_gene_name_pos <- N2_tran %>% dplyr::left_join(n2_genes_ROI, by = c("seqid","start","end")) %>% dplyr::filter(!is.na(STRAIN.y)) %>%
  dplyr::filter(tranname %in% n2_genes_interval) %>%
  dplyr::mutate(middle = (start + end) / 2) %>%
  dplyr::select(alias, middle)

label_plot_2 <- ggplot(n2_gene_name_pos) +
  geom_text(aes(x = middle / 1e6, y = 0, label = alias), angle = 75, size = 2, hjust = 0.5, vjust = 0.5) +
  # scale_x_continuous(expand = c(0.12,0.12)) +
  geom_hline(yintercept = -0.01, color = 'black') +
  # scale_y_continuous(expand = c(0,0)) +
  theme_void() +
  theme(plot.margin = margin(l = 5, r = 5, t = 5, b = -700)) +
  coord_cartesian(xlim = c((16.198034 - 0.022000), (16.198034 + 0.025000)))
label_plot_2

snvs_interval_plt <- ggplot(snv_nonMissing_long %>% dplyr::filter(POS >= min(n2_genes_ROI$start) & POS <= max(n2_genes_ROI$end)), 
                            aes(x = POS / 1e6, y = strain, fill = genotype)) +
  geom_tile() +
  scale_fill_manual(values = c("REF" = "white", "ALT" = "red")) +
  scale_x_continuous(labels = scales::comma) +
  theme_minimal() +
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 11, color = 'black'),
    axis.text.x = element_text(size = 10, color = 'black'),
    axis.text.y = element_text(size = 9, color = 'black'),
    panel.border = element_rect(color = "black", fill = NA),
    panel.grid = element_blank()) +
  labs(x = "N2 genomic position (Mb)", fill = "Genotype") +
  coord_cartesian(xlim = c((16.198034 - 0.022000), (16.198034 + 0.025000)))
# snvs_interval_plt

# Creating the final, labeled plot
top <- cowplot::plot_grid(
  label_plot_2, genes_interval,
  align = "v",
  nrow = 2)

final_plt <- cowplot::plot_grid(
  top, snvs_interval_plt,
  align = "v",
  nrow = 2,
  rel_heights = c(0.3, 1)) + theme(plot.background = element_rect(fill = "white", color = NA))
final_plt

# Save the plot!
# ggsave("../../figures/supplementary/snvs_QTL.png", final_plt, width = 7.5, height = 7.5, dpi = 600)


# Which gene has the most ALT-specific SNVs?
snvs_per_gene <- n2_genes_ROI %>%
  dplyr::mutate(gene_id = dplyr::row_number()) %>%
  dplyr::inner_join(
    alt_specific,
    by = "seqid",
    relationship = "many-to-many") %>%
  dplyr::filter(POS >= start & POS <= end)  %>%
  dplyr::count(gene_id, seqid, start, end, name = "n_snvs") %>%
  dplyr::arrange(desc(n_snvs)) %>% dplyr::left_join(N2_tran, by = c('seqid','start','end')) %>%
  dplyr::distinct(seqid, start, end, n_snvs, alias)
# WHICH OF THESE 68 VARIANTS HAVE A PREDICTED DELETERIOUS EFFECT???


glc_1_variants <- readr::read_csv("../../processed_data/gwas/CSQ_glc-1.csv") %>%
  dplyr::select(-id, -target_consequence, -gene_id) %>%
  dplyr::filter(!is.na(consequence) & consequence != "intron" & consequence != "synonymous") %>%
  tidyr::separate_rows(strains, sep = " ") %>%
  dplyr::filter(strains %in% alt_strains) %>%
  dplyr::group_by(chrom, pos, ref_seq, alt_seq, consequence) %>%
  dplyr::summarise(dplyr::across(dplyr::everything(), ~ paste(unique(.x), collapse = ", ")), .groups = "drop") %>%
  dplyr::left_join(alt_specific, by = c("pos" = "POS")) %>% 
  dplyr::filter(alt_specific == TRUE) %>%
  dplyr::filter(blosum <= 0)


f11_18_variants <- readr::read_csv("../../processed_data/gwas/CSQ_F11A5.18.csv") %>%
  dplyr::select(-id, -target_consequence, -gene_id, -dna_change, -transcript) %>%
  dplyr::filter(!is.na(consequence) & consequence != "intron" & consequence != "synonymous") %>%
  tidyr::separate_rows(strains, sep = " ") %>%
  dplyr::filter(strains %in% c("ECA36","ECA369")) %>%
  dplyr::group_by(chrom, pos, ref_seq, alt_seq, consequence) %>%
  dplyr::summarise(dplyr::across(dplyr::everything(), ~ paste(unique(.x), collapse = ", ")), .groups = "drop") #%>%
  # dplyr::left_join(alt_specific, by = c("pos" = "POS")) %>% 
  # dplyr::filter(alt_specific == TRUE) %>%
  # dplyr::filter(blosum <= 0)


f11_5_variants <- readr::read_csv("../../processed_data/gwas/CSQ_F11A5.5.csv") %>%
  dplyr::select(-id, -target_consequence, -gene_id, -dna_change) %>%
  dplyr::filter(transcript == "F11A5.5b.1") %>%
  dplyr::filter(!is.na(consequence) & consequence != "intron" & consequence != "synonymous") %>%
  tidyr::separate_rows(strains, sep = " ") %>%
  dplyr::filter(strains %in% alt_strains) %>%
  dplyr::group_by(chrom, pos, ref_seq, alt_seq, consequence, amino_acid_change) %>%
  dplyr::summarise(dplyr::across(dplyr::everything(), ~ paste(unique(.x), collapse = ", ")), .groups = "drop") %>%
  dplyr::left_join(alt_specific, by = c("pos" = "POS")) %>%
  dplyr::filter(alt_specific == TRUE) 
  # dplyr::filter(blosum <= 0)


irld26_variants <- readr::read_csv("../../processed_data/gwas/CSQ_irld26.csv") %>%
  dplyr::select(-id, -target_consequence, -gene_id, -dna_change, -transcript) %>%
  dplyr::filter(!is.na(consequence) & consequence != "intron" & consequence != "synonymous") %>%
  tidyr::separate_rows(strains, sep = " ") %>%
  dplyr::filter(strains %in% c("ECA36","ECA369")) %>%
  dplyr::group_by(chrom, pos, ref_seq, alt_seq, consequence) %>%
  dplyr::summarise(dplyr::across(dplyr::everything(), ~ paste(unique(.x), collapse = ", ")), .groups = "drop") %>%
  dplyr::filter(strains == "ECA36", grepl("missense",consequence)) %>%
  dplyr::filter(blosum < 0)



oac15_variants <- readr::read_csv("../../processed_data/gwas/CSQ_oac15.csv") %>%
  dplyr::select(-id, -target_consequence, -gene_id, -dna_change, -transcript) %>%
  dplyr::filter(!is.na(consequence) & consequence != "intron" & consequence != "synonymous") %>%
  tidyr::separate_rows(strains, sep = " ") %>%
  dplyr::filter(strains %in% c("ECA36","ECA369")) %>%
  dplyr::group_by(chrom, pos, ref_seq, alt_seq, consequence) %>%
  dplyr::summarise(dplyr::across(dplyr::everything(), ~ paste(unique(.x), collapse = ", ")), .groups = "drop") %>%
  dplyr::filter(blosum < 0)









