library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(stringr)
library(ape)

# CE and HC orthogroups
orthos <- readr::read_tsv("../../processed_data/orthology/orthofinder/proteomes_HC_CE/OrthoFinder/Results_Nov03/Orthogroups/Orthogroups.tsv")
sc_ogs <- readr::read_tsv("../../processed_data/orthology/orthofinder/proteomes_HC_CE/OrthoFinder/Results_Nov03/Orthogroups/Orthogroups_SingleCopyOrthologues.txt", col_names = "og_ids") %>%
  dplyr::pull()

cols <- colnames(orthos)
name <- gsub(".PRJNA13758.WS283.csq.PCfeaturesOnly.longest.protein","", cols)
name1 <- gsub("haemonchus_contortus.PRJEB506.WBPS19.longest_isoform.protein","HC", name)

colnames(orthos) <- name1

# Entire QTL: V:15658959-18046079 
n2_qtl_start= 15658959
n2_qtl_end = 18046079
# Restricted interval in figure 6: V:16180650 (str-106) -  16221917 (glc-1)
n2_restricted_start = 16180650
n2_restricted_end = 16221917

# N2 transcripts, gene names, and coordinates
N2_gff <- ape::read.gff("../../processed_data/genome_resources/annotation/c_elegans.PRJNA13758.WS283.csq.PCfeaturesOnly.longest.gff3")
n2_gene_info <- N2_gff %>%
  dplyr::filter(type == "gene") %>%
  dplyr::mutate(
    alias = str_extract(attributes, "(?<=Alias=)[^;]+"),
    sequence_name = str_extract(attributes, "(?<=sequence_name=)[^;]+")) %>%
  tidyr::separate(attributes, into = c("wbgene","other"), sep = ";") %>%
  dplyr::mutate(wbgene = gsub("ID=gene:","", wbgene)) %>%
  tidyr::separate(alias, into = c("alias","other"), sep = ",") %>%
  dplyr::select(chrom = seqid,start, end, wbgene, alias, sequence_name) %>%
  dplyr::filter(chrom == "V") %>%
  dplyr::filter(chrom != "MtDNA") %>%
  dplyr::mutate(region = dplyr::case_when(
    chrom == "V" & end >= n2_restricted_start & start <= n2_restricted_end ~ "restricted_QTL",
    chrom == "V" & end >= n2_qtl_start & start <= n2_qtl_end ~ "QTL",
    TRUE ~ "outside_QTL")) 

# HCON Doyle QTL: V:37250000–37550000
hc_qtl_start = 37250000
hc_qtl_end = 37550000
# HCON McIntyre QTL: V:36000000-45000000
hc_qtl_large_start = 36000000
hc_qtl_large_end = 45000000

# HCON transcripts, gene names, and coordinates
hc_tran_gene <- ape::read.gff("../../processed_data/genome_resources/annotation/haemonchus_contortus.PRJEB506.WBPS19.annotations.gene_only.longest_isoform.gff3") %>%
  dplyr::filter(type == "gene") %>%
  dplyr::mutate(
    gene = str_extract(attributes, "(?<=ID=gene:)[^;]+")) %>%
  dplyr::select(chrom = seqid, start, end, gene) %>%
  dplyr::mutate(chrom = dplyr::case_when(
    grepl("chr5", chrom) ~ "V",
    grepl("chr1", chrom) ~ "I",
    grepl("chr2", chrom) ~ "II",
    grepl("chr3", chrom) ~ "III",
    grepl("chr4", chrom) ~ "IV",
    grepl("chrX", chrom) ~ "X")) %>%
  dplyr::mutate(region_doyle = ifelse(chrom == "V" & start <= hc_qtl_end & end >= hc_qtl_start, "inQTL", "outside_QTL")) %>%
  dplyr::mutate(region_mcintyre = ifelse(chrom =="V" & start <= hc_qtl_large_end & end >= hc_qtl_large_start, "inQTL","outside_QTL")) %>%
  dplyr::filter(!is.na(gene))

# Load in orthogroup matrix
ogs <- orthos %>%
  tidyr::separate_rows(N2, sep = ", ") %>% # expanding orthogroups so each row has one N2 gene
  dplyr::mutate(N2 = sub("\\.[^.]*$", "", N2)) %>% # removing the last period in the transcript name and the isoform number
  dplyr::mutate(N2 = sub("[A-Za-z]$", "", N2),
                N2 = gsub("transcript_","", N2)) %>% # removing any letters that are in the last position in a string
  dplyr::left_join(n2_gene_info, by = c("N2" = "sequence_name")) %>% # changing N2 transcript names to gene names
  dplyr::select(Orthogroup, N2_tran = N2, N2_gene = alias, HC) #%>%
  # dplyr::group_by(Orthogroup) %>%
  # dplyr::summarise(dplyr::across(dplyr::everything(), ~ paste(unique(.x), collapse = ", ")), .groups = "drop")  # re-collapsing orthogroups

# None of the 21 N2 genes in the rectricted QTL (peak marker +- ~22 kb) have orthologs in the H. contortus narrowed QTL
restricted_n2_genes <- n2_gene_info %>% dplyr::filter(region == "restricted_QTL") %>% dplyr::pull(alias)
ogs_restricted_interval <- ogs %>% dplyr::filter(N2_gene %in% restricted_n2_genes) %>%
  tidyr::separate_rows(HC, sep = ", ") %>%
  dplyr::mutate(HC = gsub("transcript_","", HC),
                HC = sub("-.*", "", HC)) %>%
  dplyr::left_join(hc_tran_gene, by = c("HC" = "gene")) %>%
  dplyr::select(Orthogroup, N2_gene, HC, chrom, start, end, region_doyle, region_mcintyre)

# Looking at the full N2 QTL
n2_qtl_genes <- n2_gene_info %>% dplyr::filter(region == "QTL") %>% dplyr::pull(alias)

ogs_n2_qtl <- ogs %>% dplyr::filter(N2_gene %in% n2_qtl_genes) %>%
  # dplyr::filter(Orthogroup %in% sc_ogs) %>%
  tidyr::separate_rows(HC, sep = ", ") %>%
  dplyr::mutate(HC = gsub("transcript_","", HC),
                HC = sub("-.*", "", HC)) %>%
  dplyr::left_join(hc_tran_gene, by = c("HC" = "gene")) %>%
  dplyr::select(Orthogroup, N2_gene, HC, chrom, start, end, region_doyle, region_mcintyre)

