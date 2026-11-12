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
  geom_rect(aes(xmin = start / 1e6, xmax = end / 1e6, ymin = 0.9, ymax = 1.1), fill = "gray45") +
  geom_rect(data = alt_specific, aes(xmin = (POS - 10) / 1e6, xmax = (POS + 10) / 1e6, ymin = 0.8, ymax = 1.2), fill = "black") +
  geom_vline(xintercept = 16198034 / 1e6, color = 'red', linetype = 'solid') +
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
  scale_fill_manual(values = c("REF" = "white", "ALT" = "black")) +
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


f11a515_variants <- readr::read_csv("../../processed_data/gwas/variant_annotation/CSQ_F11A5.15.csv") %>%
  dplyr::select(-id, -target_consequence, -gene_id) %>%
  dplyr::mutate(gene = "F11A5.15") %>%
  dplyr::filter(!is.na(consequence) & consequence != "intron" & consequence != "synonymous") %>%
  tidyr::separate_rows(strains, sep = " ") %>%
  dplyr::filter(strains %in% alt_strains) %>%
  dplyr::group_by(chrom, pos, ref_seq, alt_seq, consequence) %>%
  dplyr::summarise(dplyr::across(dplyr::everything(), ~ paste(unique(.x), collapse = ", ")), .groups = "drop") %>%
  dplyr::left_join(alt_specific, by = c("pos" = "POS")) %>% 
  dplyr::filter(alt_specific == TRUE) #%>%
  # dplyr::filter(blosum <= 0)

srh279_variants <- readr::read_csv("../../processed_data/gwas/variant_annotation/CSQ_srh279.csv") %>%
  dplyr::select(-id, -target_consequence, -gene_id) %>%
  dplyr::filter(!is.na(consequence) & consequence != "intron" & consequence != "synonymous") %>%
  tidyr::separate_rows(strains, sep = " ") %>%
  dplyr::filter(strains %in% alt_strains) %>%
  dplyr::group_by(chrom, pos, ref_seq, alt_seq, consequence) %>%
  dplyr::summarise(dplyr::across(dplyr::everything(), ~ paste(unique(.x), collapse = ", ")), .groups = "drop") %>%
  dplyr::left_join(alt_specific, by = c("pos" = "POS")) %>% 
  dplyr::filter(alt_specific == TRUE) %>%
  dplyr::filter(consequence != "5_prime_utr")

f11a53_variants <- readr::read_csv("../../processed_data/gwas/variant_annotation/CSQ_F11A5.3.csv") %>%
  dplyr::select(-id, -target_consequence, -gene_id) %>%
  dplyr::filter(!is.na(consequence) & consequence != "intron" & consequence != "synonymous") %>%
  tidyr::separate_rows(strains, sep = " ") %>%
  dplyr::filter(strains %in% alt_strains) %>%
  dplyr::group_by(chrom, pos, ref_seq, alt_seq, consequence) %>%
  dplyr::summarise(dplyr::across(dplyr::everything(), ~ paste(unique(.x), collapse = ", ")), .groups = "drop") %>%
  dplyr::left_join(alt_specific, by = c("pos" = "POS")) %>% 
  dplyr::filter(alt_specific == TRUE) 


f11a517_variants <- readr::read_csv("../../processed_data/gwas/variant_annotation/CSQ_F11A5.17.csv") %>%
  dplyr::select(-id, -target_consequence, -gene_id) %>%
  dplyr::mutate(gene = "F11A5.17") %>%
  dplyr::filter(!is.na(consequence) & consequence != "intron" & consequence != "synonymous") %>%
  tidyr::separate_rows(strains, sep = " ") %>%
  dplyr::filter(strains %in% alt_strains) %>%
  dplyr::group_by(chrom, pos, ref_seq, alt_seq, consequence) %>%
  dplyr::summarise(dplyr::across(dplyr::everything(), ~ paste(unique(.x), collapse = ", ")), .groups = "drop") %>%
  dplyr::left_join(alt_specific, by = c("pos" = "POS")) %>% 
  dplyr::filter(alt_specific == TRUE) 

f11a54_variants <- readr::read_csv("../../processed_data/gwas/variant_annotation/CSQ_F11A5.4.csv") %>%
  dplyr::select(-id, -target_consequence, -gene_id) %>%
  dplyr::mutate(gene = "F11A5.4") %>%
  dplyr::filter(!is.na(consequence) & consequence != "intron" & consequence != "synonymous") %>%
  tidyr::separate_rows(strains, sep = " ") %>%
  dplyr::filter(strains %in% alt_strains) %>%
  dplyr::group_by(chrom, pos, ref_seq, alt_seq, consequence) %>%
  dplyr::summarise(dplyr::across(dplyr::everything(), ~ paste(unique(.x), collapse = ", ")), .groups = "drop") %>%
  dplyr::left_join(alt_specific, by = c("pos" = "POS")) %>% 
  dplyr::filter(alt_specific == TRUE) 

f11a59_variants <- readr::read_csv("../../processed_data/gwas/variant_annotation/CSQ_F11A5.9.csv") %>%
  dplyr::select(-id, -target_consequence, -gene_id) %>%
  dplyr::mutate(gene = "F11A5.9") %>%
  dplyr::filter(!is.na(consequence) & consequence != "intron" & consequence != "synonymous") %>%
  tidyr::separate_rows(strains, sep = " ") %>%
  dplyr::filter(strains %in% alt_strains) %>%
  dplyr::group_by(chrom, pos, ref_seq, alt_seq, consequence) %>%
  dplyr::summarise(dplyr::across(dplyr::everything(), ~ paste(unique(.x), collapse = ", ")), .groups = "drop") %>%
  dplyr::left_join(alt_specific, by = c("pos" = "POS")) %>% 
  dplyr::filter(alt_specific == TRUE) 


glc_1_variants <- readr::read_csv("../../processed_data/gwas/variant_annotation/CSQ_glc-1.csv") %>%
  dplyr::select(-id, -target_consequence, -gene_id) %>%
  dplyr::filter(!is.na(consequence) & consequence != "intron" & consequence != "synonymous") %>%
  tidyr::separate_rows(strains, sep = " ") %>%
  dplyr::filter(strains %in% alt_strains) %>%
  dplyr::group_by(chrom, pos, ref_seq, alt_seq, consequence) %>%
  dplyr::summarise(dplyr::across(dplyr::everything(), ~ paste(unique(.x), collapse = ", ")), .groups = "drop") %>%
  dplyr::left_join(alt_specific, by = c("pos" = "POS")) %>% 
  dplyr::filter(alt_specific == TRUE) %>%
  dplyr::filter(blosum <= 0)

f11_18_variants <- readr::read_csv("../../processed_data/gwas/variant_annotation/CSQ_F11A5.18.csv") %>%
  dplyr::select(-id, -target_consequence, -gene_id, -dna_change, -transcript) %>%
  dplyr::filter(!is.na(consequence) & consequence != "intron" & consequence != "synonymous") %>%
  dplyr::mutate(gene = "F11A5.18") %>%
  tidyr::separate_rows(strains, sep = " ") %>%
  dplyr::filter(strains %in% c("ECA36","ECA369")) %>%
  dplyr::group_by(chrom, pos, ref_seq, alt_seq, consequence) %>%
  dplyr::summarise(dplyr::across(dplyr::everything(), ~ paste(unique(.x), collapse = ", ")), .groups = "drop") %>%
  dplyr::filter(consequence != "*synonymous")
  # dplyr::left_join(alt_specific, by = c("pos" = "POS")) %>% 
  # dplyr::filter(alt_specific == TRUE) %>%
  # dplyr::filter(blosum <= 0)

f11_5_variants <- readr::read_csv("../../processed_data/gwas/variant_annotation/CSQ_F11A5.5.csv") %>%
  dplyr::select(-id, -target_consequence, -gene_id, -dna_change) %>%
  dplyr::filter(transcript == "F11A5.5b.1") %>%
  dplyr::mutate(gene = "F11A5.5") %>%
  dplyr::filter(!is.na(consequence) & consequence != "intron" & consequence != "synonymous") %>%
  tidyr::separate_rows(strains, sep = " ") %>%
  dplyr::filter(strains %in% alt_strains) %>%
  dplyr::group_by(chrom, pos, ref_seq, alt_seq, consequence, amino_acid_change) %>%
  dplyr::summarise(dplyr::across(dplyr::everything(), ~ paste(unique(.x), collapse = ", ")), .groups = "drop") %>%
  dplyr::left_join(alt_specific, by = c("pos" = "POS")) %>%
  dplyr::filter(alt_specific == TRUE) 
  # dplyr::filter(blosum <= 0)


irld26_variants <- readr::read_csv("../../processed_data/gwas/variant_annotation/CSQ_irld26.csv") %>%
  dplyr::select(-id, -target_consequence, -gene_id, -dna_change, -transcript) %>%
  dplyr::filter(!is.na(consequence) & consequence != "intron" & consequence != "synonymous") %>%
  tidyr::separate_rows(strains, sep = " ") %>%
  dplyr::filter(strains %in% c("ECA36","ECA369")) %>%
  dplyr::group_by(chrom, pos, ref_seq, alt_seq, consequence) %>%
  dplyr::summarise(dplyr::across(dplyr::everything(), ~ paste(unique(.x), collapse = ", ")), .groups = "drop") %>%
  dplyr::filter(strains == "ECA36", grepl("missense",consequence)) #%>%
  # dplyr::filter(blosum < 0)

oac15_variants <- readr::read_csv("../../processed_data/gwas/variant_annotation/CSQ_oac15.csv") %>%
  dplyr::select(-id, -target_consequence, -gene_id, -dna_change, -transcript) %>%
  dplyr::filter(!is.na(consequence) & consequence != "intron" & consequence != "synonymous") %>%
  tidyr::separate_rows(strains, sep = " ") %>%
  dplyr::filter(strains %in% c("ECA36","ECA369")) %>%
  dplyr::group_by(chrom, pos, ref_seq, alt_seq, consequence) %>%
  dplyr::summarise(dplyr::across(dplyr::everything(), ~ paste(unique(.x), collapse = ", ")), .groups = "drop") %>%
  dplyr::filter(strains != "ECA369") %>%
  dplyr::filter(!grepl("utr",consequence))
  # dplyr::filter(blosum < 0)


avr14_variants <- readr::read_csv("../../processed_data/gwas/variant_annotation/CSQ_avr14.csv") %>%
  dplyr::select(-id, -target_consequence, -gene_id) %>%
  dplyr::filter(!is.na(consequence) & consequence != "intron" & consequence != "synonymous") %>%
  tidyr::separate_rows(strains, sep = " ") %>%
  dplyr::filter(strains %in% alt_strains) %>%
  dplyr::group_by(chrom, pos, ref_seq, alt_seq, consequence) %>%
  dplyr::summarise(dplyr::across(dplyr::everything(), ~ paste(unique(.x), collapse = ", ")), .groups = "drop") %>%
  dplyr::left_join(alt_specific, by = c("pos" = "POS")) %>% 
  dplyr::filter(alt_specific == TRUE)

avr15_variants <- readr::read_csv("../../processed_data/gwas/variant_annotation/CSQ_avr15.csv") %>%
  dplyr::select(-id, -target_consequence, -gene_id) %>%
  dplyr::filter(!is.na(consequence) & consequence != "intron" & consequence != "synonymous") %>%
  tidyr::separate_rows(strains, sep = " ") %>%
  dplyr::filter(strains %in% alt_strains) %>%
  dplyr::group_by(chrom, pos, ref_seq, alt_seq, consequence) %>%
  dplyr::summarise(dplyr::across(dplyr::everything(), ~ paste(unique(.x), collapse = ", ")), .groups = "drop") %>%
  dplyr::left_join(alt_specific, by = c("pos" = "POS")) %>% 
  dplyr::filter(alt_specific == TRUE)


# Concatenating together information from all genes in QTL with shared SNVs and avr-14/15
resistant_snvs_per_gene <- snvs_per_gene %>% dplyr::mutate(n_snvs = as.character(n_snvs)) %>% 
  dplyr::add_row(alias = "F11A5.18", n_snvs = "ECA36") %>%
  dplyr::add_row(alias = "irld-26", n_snvs = "ECA36") %>%
  dplyr::add_row(alias = "oac-15", n_snvs = "ECA36") %>%
  dplyr::add_row(alias = "avr-15", n_snvs = "0") %>%
  dplyr::add_row(alias = "avr-14", n_snvs = "0") %>%
  dplyr::select(Gene = alias, "Number of SNVs specific to resistant strains" = n_snvs) %>% 
  dplyr::left_join(glc_1_variants %>% dplyr::select(gene, Chromosome = chrom, Position = pos, REF = ref_seq, ALT = alt_seq, consequence, amino_acid_change, blosum, grantham), by = c("Gene" = "gene")) %>%
  # dplyr::mutate(amino_acid_change = ifelse(amino_acid_change == "68A>64T, 68A>66T, 68A>68T", "A68T", 
                                           # ifelse(amino_acid_change == "163T>159M, 163T>162M, 163T>163M", "T163M",
                                                  # ifelse(amino_acid_change == "367V>363A, 367V>367A", "V367A", amino_acid_change)))) %>%
  dplyr::left_join(f11_18_variants %>% dplyr::select(gene,  Chromosome = chrom, Position = pos, REF = ref_seq, ALT = alt_seq, consequence, amino_acid_change, blosum, grantham), by = c("Gene" = "gene")) %>%
  dplyr::mutate(consequence = coalesce(consequence.x, consequence.y),
                Chromosome = coalesce(Chromosome.x, Chromosome.y),
                Position = coalesce(Position.x, Position.y),
                REF = coalesce(REF.x, REF.y),
                ALT = coalesce(ALT.x, ALT.y),
                amino_acid_change = coalesce(amino_acid_change.x, amino_acid_change.y),
                BLOSUM62 = coalesce(blosum.x, blosum.y),
                Grantham = coalesce(grantham.x, grantham.y)) %>% dplyr::select(Gene, Chromosome, Position, REF, ALT, `Number of SNVs specific to resistant strains`, consequence, amino_acid_change, BLOSUM62, Grantham) %>%
  dplyr::left_join(f11_5_variants %>% dplyr::select(gene,  Chromosome = chrom, Position = pos, REF = ref_seq, ALT = alt_seq, consequence, amino_acid_change, blosum, grantham), by = c("Gene" = "gene")) %>%
  dplyr::mutate(consequence = coalesce(consequence.x, consequence.y),
                amino_acid_change = coalesce(amino_acid_change.x, amino_acid_change.y),
                Chromosome = coalesce(Chromosome.x, Chromosome.y),
                Position = coalesce(Position.x, Position.y),
                REF = coalesce(REF.x, REF.y),
                ALT = coalesce(ALT.x, ALT.y),
                BLOSUM62 = coalesce(BLOSUM62, blosum),
                Grantham = coalesce(Grantham, grantham)) %>% dplyr::select(Gene, Chromosome, Position, REF, ALT, `Number of SNVs specific to resistant strains`, consequence, amino_acid_change, BLOSUM62, Grantham) %>%
  dplyr::left_join(irld26_variants %>% dplyr::select(gene,  Chromosome = chrom, Position = pos, REF = ref_seq, ALT = alt_seq, consequence, amino_acid_change, blosum, grantham), by = c("Gene" = "gene")) %>%
  dplyr::mutate(consequence = coalesce(consequence.x, consequence.y),
                amino_acid_change = coalesce(amino_acid_change.x, amino_acid_change.y),
                Chromosome = coalesce(Chromosome.x, Chromosome.y),
                Position = coalesce(Position.x, Position.y),
                REF = coalesce(REF.x, REF.y),
                ALT = coalesce(ALT.x, ALT.y),
                BLOSUM62 = coalesce(BLOSUM62, blosum),
                Grantham = coalesce(Grantham, grantham)) %>% dplyr::select(Gene, Chromosome, Position, REF, ALT, `Number of SNVs specific to resistant strains`, consequence, amino_acid_change, BLOSUM62, Grantham)%>%
  dplyr::left_join(oac15_variants %>% dplyr::select(gene,  Chromosome = chrom, Position = pos, REF = ref_seq, ALT = alt_seq, consequence, amino_acid_change, blosum, grantham), by = c("Gene" = "gene")) %>%
  dplyr::mutate(consequence = coalesce(consequence.x, consequence.y),
                amino_acid_change = coalesce(amino_acid_change.x, amino_acid_change.y),
                Chromosome = coalesce(Chromosome.x, Chromosome.y),
                Position = coalesce(Position.x, Position.y),
                REF = coalesce(REF.x, REF.y),
                ALT = coalesce(ALT.x, ALT.y),
                BLOSUM62 = coalesce(BLOSUM62, blosum),
                Grantham = coalesce(Grantham, grantham)) %>% dplyr::select(Gene, Chromosome, Position, REF, ALT, `Number of SNVs specific to resistant strains`, consequence, amino_acid_change, BLOSUM62, Grantham) %>%
  dplyr::left_join(f11a515_variants %>% dplyr::select(gene,  Chromosome = chrom, Position = pos, REF = ref_seq, ALT = alt_seq, consequence, amino_acid_change, blosum, grantham), by = c("Gene" = "gene")) %>%
  dplyr::mutate(consequence = coalesce(consequence.x, consequence.y),
                amino_acid_change = coalesce(amino_acid_change.x, amino_acid_change.y),
                Chromosome = coalesce(Chromosome.x, Chromosome.y),
                Position = coalesce(Position.x, Position.y),
                REF = coalesce(REF.x, REF.y),
                ALT = coalesce(ALT.x, ALT.y),
                BLOSUM62 = coalesce(BLOSUM62, blosum),
                Grantham = coalesce(Grantham, grantham)) %>% dplyr::select(Gene, Chromosome, Position, REF, ALT, `Number of SNVs specific to resistant strains`, consequence, amino_acid_change, BLOSUM62, Grantham) %>%
  dplyr::left_join(f11a517_variants %>%dplyr::select(gene,  Chromosome = chrom, Position = pos, REF = ref_seq, ALT = alt_seq, consequence, amino_acid_change, blosum, grantham), by = c("Gene" = "gene")) %>%
  dplyr::mutate(consequence = coalesce(consequence.x, consequence.y),
                amino_acid_change = coalesce(amino_acid_change.x, amino_acid_change.y),
                Chromosome = coalesce(Chromosome.x, Chromosome.y),
                Position = coalesce(Position.x, Position.y),
                REF = coalesce(REF.x, REF.y),
                ALT = coalesce(ALT.x, ALT.y),
                BLOSUM62 = coalesce(BLOSUM62, blosum),
                Grantham = coalesce(Grantham, grantham)) %>% dplyr::select(Gene, Chromosome, Position, REF, ALT, `Number of SNVs specific to resistant strains`, consequence, amino_acid_change, BLOSUM62, Grantham) %>%
  dplyr::mutate(Chromosome = ifelse(Gene == "avr-14", "I", 
                                    ifelse(is.na(Chromosome), "V", Chromosome)))
  
# ADDING DOMAIN ANNTOATIONS
glc_1_domains <- readr::read_tsv("../../processed_data/gwas/variant_annotation/N2_GLC-1_IPR_annotations.tsv", col_names = c("tran","start","end","IPR_acc","IPR_description")) %>%
  dplyr::mutate(gene = "glc-1")

glc_1_variants_anno <- glc_1_variants %>%
  dplyr::select(chrom, pos, ref_seq, alt_seq, consequence, gene, blosum, grantham, amino_acid_change) %>%
  dplyr::mutate(extra_aa = amino_acid_change) %>%
  tidyr::separate(extra_aa, into = c("aa_pos","rest",sep = ">")) %>%
  dplyr::select(chrom, pos, ref_seq, alt_seq, consequence, gene, blosum, grantham, amino_acid_change, aa_pos) %>%
  dplyr::mutate(aa_pos = as.numeric(sub("[A-Z]$", "", aa_pos))) %>%
  dplyr::left_join(glc_1_domains, by = "gene") %>% dplyr::select(-tran) %>%
  dplyr::mutate(in_domain = ifelse(aa_pos >= start & aa_pos <= end, TRUE, FALSE)) %>%
  dplyr::filter(in_domain == TRUE, IPR_description != "-") %>%
  dplyr::select(-start, -end) %>%
  dplyr::distinct() %>%
  dplyr::group_by(chrom, pos, ref_seq, alt_seq, consequence, gene, blosum, grantham, amino_acid_change, aa_pos) %>%
  dplyr::summarise(dplyr::across(dplyr::everything(), ~ paste(unique(.x), collapse = ", ")), .groups = "drop") %>%
  dplyr::select(Gene = gene, Position = pos, REF = ref_seq, ALT = alt_seq, consequence, amino_acid_change, BLOSUM62 = blosum, Grantham = grantham, IPR_acc, IPR_description)


f11a55_domains <- readr::read_tsv("../../processed_data/gwas/variant_annotation/N2_F11A5.5_IPR_annotations.tsv", col_names = c("tran","start","end","IPR_acc","IPR_description")) %>%
  dplyr::mutate(gene = "F11A5.5")

f11a55_variants_anno <- f11_5_variants %>%
  dplyr::select(chrom, pos, ref_seq, alt_seq, consequence, gene, blosum, grantham, amino_acid_change) %>%
  dplyr::filter(grepl("missense", consequence)) %>%
  dplyr::mutate(extra_aa = amino_acid_change) %>%
  tidyr::separate(extra_aa, into = c("aa_pos","rest",sep = ">")) %>%
  dplyr::select(chrom, pos, ref_seq, alt_seq, consequence, gene, blosum, grantham, amino_acid_change, aa_pos) %>%
  dplyr::mutate(aa_pos = as.numeric(sub("[A-Z]$", "", aa_pos))) %>%
  dplyr::left_join(f11a55_domains, by = "gene") %>% dplyr::select(-tran) %>%
  dplyr::mutate(in_domain = ifelse(aa_pos >= start & aa_pos <= end, TRUE, FALSE)) %>%
  dplyr::filter(in_domain == TRUE, IPR_description != "-") %>%
  dplyr::select(-start, -end) %>%
  dplyr::distinct() %>%
  dplyr::group_by(chrom, pos, ref_seq, alt_seq, consequence, gene, blosum, grantham, amino_acid_change, aa_pos) %>%
  dplyr::summarise(dplyr::across(dplyr::everything(), ~ paste(unique(.x), collapse = ", ")), .groups = "drop") %>%
  dplyr::select(Gene = gene, Position = pos, REF = ref_seq, ALT = alt_seq, consequence, amino_acid_change, BLOSUM62 = blosum, Grantham = grantham, IPR_acc, IPR_description)



irld26_domains <- readr::read_tsv("../../processed_data/gwas/variant_annotation/N2_IRLD-26_IPR_annotations.tsv", col_names = c("tran","start","end","IPR_acc","IPR_description")) %>%
  dplyr::mutate(gene = "irld-26")

irld26_variants_anno <- irld26_variants %>%
  dplyr::select(chrom, pos, ref_seq, alt_seq, consequence, gene, blosum, grantham, amino_acid_change) %>%
  dplyr::filter(grepl("missense", consequence)) %>%
  dplyr::mutate(extra_aa = amino_acid_change) %>%
  tidyr::separate(extra_aa, into = c("aa_pos","rest",sep = ">")) %>%
  dplyr::select(chrom, pos, ref_seq, alt_seq, consequence, gene, blosum, grantham, amino_acid_change, aa_pos) %>%
  dplyr::mutate(aa_pos = as.numeric(sub("[A-Z]$", "", aa_pos))) %>%
  dplyr::left_join(irld26_domains, by = "gene") %>% dplyr::select(-tran) %>%
  dplyr::mutate(in_domain = ifelse(aa_pos >= start & aa_pos <= end, TRUE, FALSE)) %>%
  dplyr::filter(in_domain == TRUE, IPR_description != "-") %>%
  dplyr::select(-start, -end) %>%
  dplyr::distinct() %>%
  dplyr::group_by(chrom, pos, ref_seq, alt_seq, consequence, gene, blosum, grantham, amino_acid_change, aa_pos) %>%
  dplyr::summarise(dplyr::across(dplyr::everything(), ~ paste(unique(.x), collapse = ", ")), .groups = "drop") %>%
  dplyr::select(Gene = gene, Position = pos, REF = ref_seq, ALT = alt_seq, consequence, amino_acid_change, BLOSUM62 = blosum, Grantham = grantham, IPR_acc, IPR_description)



oac15_domains <- readr::read_tsv("../../processed_data/gwas/variant_annotation/N2_OAC-15_IPR_annotations.tsv", col_names = c("tran","start","end","IPR_acc","IPR_description")) %>%
  dplyr::mutate(gene = "oac-15")

oac15_variants_anno <- oac15_variants %>%
  dplyr::select(chrom, pos, ref_seq, alt_seq, consequence, gene, blosum, grantham, amino_acid_change) %>%
  dplyr::filter(grepl("missense", consequence)) %>%
  dplyr::mutate(extra_aa = amino_acid_change) %>%
  tidyr::separate(extra_aa, into = c("aa_pos","rest",sep = ">")) %>%
  dplyr::select(chrom, pos, ref_seq, alt_seq, consequence, gene, blosum, grantham, amino_acid_change, aa_pos) %>%
  dplyr::mutate(aa_pos = as.numeric(sub("[A-Z]$", "", aa_pos))) %>%
  dplyr::left_join(oac15_domains, by = "gene") %>% dplyr::select(-tran) %>%
  dplyr::mutate(in_domain = ifelse(aa_pos >= start & aa_pos <= end, TRUE, FALSE)) %>%
  dplyr::filter(in_domain == TRUE, IPR_description != "-") %>%
  dplyr::select(-start, -end) %>%
  dplyr::distinct() %>%
  dplyr::group_by(chrom, pos, ref_seq, alt_seq, consequence, gene, blosum, grantham, amino_acid_change, aa_pos) %>%
  dplyr::summarise(dplyr::across(dplyr::everything(), ~ paste(unique(.x), collapse = ", ")), .groups = "drop") %>%
  dplyr::select(Gene = gene, Position = pos, REF = ref_seq, ALT = alt_seq, consequence, amino_acid_change, BLOSUM62 = blosum, Grantham = grantham, IPR_acc, IPR_description)

# F11A5.15, F11A5.17, and F11A5.18 don't have any IPR annotations


# Concatenating all domain annotations together
all_domains <- glc_1_variants_anno %>% dplyr::bind_rows(f11a55_variants_anno, irld26_variants_anno, oac15_variants_anno)

# Create final table with domain annotations added
resistant_snvs_per_gene_wDomains <- resistant_snvs_per_gene %>% dplyr::left_join(all_domains, by = c("Gene","Position","REF","ALT","consequence","amino_acid_change","BLOSUM62","Grantham")) %>%
  dplyr::rename(Consequence = consequence, `Amino acid change` = amino_acid_change, `IPR accession` = IPR_acc, `IPR description` = IPR_description)

write.table(resistant_snvs_per_gene_wDomains, "../../tables/QTL_annotated_SNVs.tsv", quote = F, row.names = F, col.names = T, sep = "\t")

