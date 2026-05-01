library(plyr) # ALWAYS LOAD BEFORE DPLYR
library(readr)
library(dplyr)
library(ggplot2)
library(stringr)
library(ape)
library(tidyr)
library(data.table)



####################################################################################################
####################################################################################################

# GENE SET CLASSIFICATION AND VISUALIZATION

####################################################################################################
####################################################################################################

# ======================================================================================================================================================================================== #
# Pulling all genes, coordinates, and alignments for all WSs and N2 #
# ======================================================================================================================================================================================== #
# genes_strain <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/raw_data/assemblies/elegans/gff/116_genesOnly_strainRes.tsv", col_names = c("contig","type", "start", "end", "strand", "attributes", "strain")) 
genes_strain <- readr::read_tsv("../../processed_data/genome_resources/genomes/140Ws_CGC1_longestIsoGenes_BRAKER.tsv", col_names = c("seqid","source", "type", "start", "end", "score", "strand", "phase", "attributes", "strain")) %>% dplyr::filter(strain != "ECA396")
N2_gff <- ape::read.gff("../../processed_data/genome_resources/genomes/c_elegans.PRJNA13758.WS283.csq.PCfeaturesOnly.longest.gff3") %>% dplyr::mutate(strain="N2")
genes_strain <- rbind(genes_strain,N2_gff)
all_genes_strain <- genes_strain %>%
  dplyr::mutate(attributes = gsub("ID=gene:","",attributes)) %>%
  dplyr::mutate(attributes = gsub("ID=","",attributes)) %>%
  dplyr::mutate(attributes = sub(";.*", "", attributes)) %>%
  dplyr::filter(type == "gene")

N2_tranGene <- N2_gff %>%
  dplyr::filter(type == "mRNA") %>%
  dplyr::mutate(attributes = gsub("ID=transcript:","",attributes), attributes = gsub("Parent=gene:","",attributes)) %>%
  tidyr::separate_wider_delim(attributes, delim = ";",names = c("tran", "gene", "rest"), too_many = "merge") %>%
  dplyr::select(tran,gene, -rest)

# dplyr::mutate(gene = str_extract(attributes,"(?<=\\bID=gene:)WBGene\\d+|^WBGene\\d+|(?<=\\bName=)WBGene\\d+"),tran = str_extract(attributes, "(?<=\\bsequence_name=)[^;]+")) %>%
# dplyr::select(tran,gene)

# test1 <- all_genes_strain %>%
#   dplyr::filter(strain == "N2") # ~19,000

# write.table(all_genes_strain,"/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/orthofinder/116strain_genes.tsv", quote = F, row.names = F, col.names = T, sep = '\t')

nucmer <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/synteny_vis/elegans/nucmer_aln_WSs/142_nucmer_ECA741CGC1.tsv", col_names = c("N2S","N2E","WSS","WSE","L1","L2","IDY","LENR","LENQ","N2_chr","contig","strain")) %>%
  dplyr::select(-L1,-L2,-IDY,-LENR,-LENQ) %>% dplyr::filter(strain != "ECA396")
# write.table(nucmer,"/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/orthofinder/115WS_nucmer_clean.tsv", quote = F, row.names = F, col.names = T, sep = '\t')





##### Verifying that transcript to gene conversion works correctly ############################################################################################################
# orthogroups_tran <- readr::read_tsv("/vast/eande106/projects/Nicolas/WI_PacBio_genomes/orthology/elegans/prot_115/OrthoFinder/Results_Apr22/Phylogenetic_Hierarchical_Orthogroups/N0.tsv")
# strCol <- colnames(orthogroups_tran)
# strCol_c1 <- gsub(".braker.protein","",strCol)
# strCol_c2 <- gsub("_WS283.protein","",strCol_c1)
# colnames(orthogroups_tran) <- strCol_c2
# 
# # print(nrow(orthogroups_tran)) # 64377
# 
# #### Counting number of transcripts per orthogroup for each strain #### 
# ortho_count_tran <- orthogroups_tran
# 
# strCol_c2_u <- strCol_c2[!strCol_c2 %in% c("OG", "HOG", "Gene Tree Parent Clade")]
# 
# for (i in 1:length(strCol_c2_u)) {
#   print(paste0(i,"out of", length(strCol_c2_u)))
#   temp_colname = paste0(strCol_c2_u[i], "_count")
# 
#   ortho_count_tran <- ortho_count_tran %>%
#     dplyr::mutate(!!sym(temp_colname) := stringr::str_count(!!sym(strCol_c2_u[i]),", ") + 1)
# }
# 
# 
# all_relations_tran <- ortho_count_tran %>%
#   dplyr::select(HOG, dplyr::contains("_count"))
# 
# 
# #### Transcripts converted to genes #### see script: /vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/scripts/Ce_geneAnno-sh/asm_gene_set/tran_gene.sh
# ortho_genes <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/orthofinder/N0_0422_genes.tsv")
# 
# stCol <- colnames(ortho_genes)
# stCol_c1 <- gsub(".braker.protein","",stCol)
# stCol_c2 <- gsub("_WS283.protein","",stCol_c1)
# colnames(ortho_genes) <- stCol_c2
# 
# ortho_count_dup <- ortho_genes
# 
# stCol_c2_u <- stCol_c2[!stCol_c2 %in% c("OG", "HOG", "Gene Tree Parent Clade")]
# 
# for (i in 1:length(stCol_c2_u)) {
#   print(paste0(i,"out of", length(stCol_c2_u)))
#   temp_colname = paste0(stCol_c2_u[i], "_count")
#   
#   ortho_count_dup <- ortho_count_dup %>%
#     dplyr::mutate(!!sym(temp_colname) := stringr::str_count(!!sym(stCol_c2_u[i]),", ") + 1)
# }
# 
# all_relations_duplicated_genes <- ortho_count_dup %>%
#   dplyr::select(HOG, dplyr::contains("_count"))
# 
# ## all_relations and all_relations_g should be identical
# identical(all_relations_tran, all_relations_duplicated_genes) # TRUE
#################################################################################################################################################################################








# 
# 
# N2_braker <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/orthofinder/N0_N2_braker_genes.tsv")
# 
# 
# strainCol <- colnames(N2_braker)
# strainCol_c1 <- gsub(".braker.longest.protein","",strainCol)
# strainCol_c2 <- gsub(".longest.protein","",strainCol_c1)
# colnames(N2_braker) <- strainCol_c2
# 
# ortho_count <- N2_braker
# 
# strainCol_c2_u <- strainCol_c2[!strainCol_c2 %in% c("OG", "HOG", "Gene Tree Parent Clade")]
# 
# for (i in 1:length(strainCol_c2_u)) {
#   print(paste0(i,"out of", length(strainCol_c2_u)))
#   temp_colname = paste0(strainCol_c2_u[i], "_count")
#   
#   ortho_count <- ortho_count %>%
#     dplyr::mutate(!!sym(temp_colname) := stringr::str_count(!!sym(strainCol_c2_u[i]),", ") + 1)
# }
# 
# all_relations <- ortho_count %>%
#   dplyr::select(HOG, dplyr::contains("_count"))
# 
# 
# #### Plotting classification based on all HOGs ####
# private_freq = (0.5)
# 
# classification <- all_relations %>%
#   dplyr::mutate(across(2:(ncol(.)), ~ ifelse(. >= 1, 1, .))) %>%
#   dplyr::mutate(sum = rowSums(across(-1, ~ ., .names = NULL), na.rm = TRUE)) %>%
#   dplyr::mutate(freq = (sum / length(strainCol_c2_u))) %>%
#   dplyr::mutate(
#     class = case_when(
#       freq == 1 ~ "core",
#       freq == private_freq ~ "private",
#       TRUE ~ "undefined"
#     )
#   ) %>%
#   dplyr::count(freq, class) %>%
#   dplyr::mutate(percent = (n / sum(n)) * 100) 
# 
# gs_allOrtho <- ggplot(data = classification, aes(x = freq * 100, y = percent, fill = class)) + 
#   geom_bar(stat = "identity", color = "black", alpha = 0.5) + 
#   scale_fill_manual(values = c(
#     "core" = "green4",
#     "private" = "magenta3"
#   ), 
#   limits = c("core", "private"),  # Manually ordering legend items
#   guide = guide_legend(title = NULL) 
#   ) +
#   ylab("Percent of HOGs - Longest Isoform") +
#   xlab("Frequency") +
#   # ggtitle("All orthogroups") +
#   scale_y_continuous(labels = scales::percent_format(scale = 1)) + 
#   scale_x_continuous(labels = scales::percent_format(scale = 1)) +
#   theme_classic() +
#   theme(
#     axis.title = element_text(size = 16),
#     legend.position = c(0.85, 0.8),
#     plot.title = element_text(size=18, face = 'bold', hjust=0.5),
#     legend.text = element_text(size=13, color = 'black'),
#     axis.text = element_text(size=12, color = 'black')
#   )
# gs_allOrtho
# 
# HOG_class_count <- classification %>%
#   dplyr::group_by(class) %>%
#   dplyr::summarise(n_HOG = sum(n)) %>%
#   dplyr::ungroup()
# 
# N2_braker <- all_relations %>%
#   dplyr::select(N2_count) %>%
#   dplyr::filter(is.na(N2_count)) # 25 - number of private HOGs belonging to the N2.WS283
# 
# N2_WS283 <- all_relations %>%
#   dplyr::select(N2.WS283_count) %>%
#   dplyr::filter(is.na(N2.WS283_count)) # 96 - number of private HOGs belonging to the BRAKER predictions
# 
# n2_braker_genes <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/misc/N2_BRAKER_error/N2.BRAKER.longest.tsv", col_names = c("contig", "type", "start", "end", "strand", "attributes")) 
# ws283_genes <- all_genes_strain %>%
#   dplyr::filter(strain == "N2")
# 
# n2_braker_WS283_genes <- n2_braker_genes %>%
#   dplyr::mutate(attributes = gsub("ID=gene:","",attributes)) %>%
#   dplyr::mutate(attributes = gsub("ID=","",attributes)) %>%
#   dplyr::mutate(attributes = sub(";.*", "", attributes)) %>%
#   dplyr::mutate(strain = "N2_BRAKER") %>%
#   dplyr::bind_rows(ws283_genes) 
# 
# n2_gene_count <- n2_braker_WS283_genes %>%
#   dplyr::count(strain, name = "n_genes")
# 
# ggplot(n2_gene_count) + 
#   geom_bar(aes(x = strain, y = n_genes, fill = strain), stat = "identity") +
#   geom_text(aes(x = strain, y = n_genes, label = n_genes), vjust = -0.3, size = 5) +
#   ylab("Number of PC genes") +
#   theme(
#     legend.position = 'none',
#     panel.background = element_blank(),
#     panel.border = element_rect(fill = NA),
#     axis.title.x = element_blank(),
#     axis.title.y = element_text(size = 14, face = 'bold'),
#     axis.text.y = element_text(size = 13),
#     axis.text.x = element_text(face = 'bold', size=14, color = 'black'))








# ======================================================================================================================================================================================== #
# HOG matrix manipulation and plotting #
# ======================================================================================================================================================================================== #

# Converting transcripts to genes and removing all duplicate genes (not needed anymore, using longest isoform) in a dataframe cell:
# see script - /vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/scripts/Ce_geneAnno-sh/asm_gene_set/tran_gene.sh 
# ortho_genes_dd <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/orthofinder/May_115_longestIso/N0_0504_genes.tsv")
ortho_genes_dd <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/orthology/elegans/orthofinder/64_core/OrthoFinder/Results_Dec07/Orthogroups/Orthogroups.tsv") %>%
  dplyr::filter(!grepl("MTCE",c_elegans.PRJNA13758.WS283.csq.PCfeaturesOnly.longest.protein))
# whatever <- ortho_genes_dd %>%
# dplyr::select(JU1581.braker.longest.protein, N2.longest.protein)

strainCol <- colnames(ortho_genes_dd)
ugh <- gsub(".20251012.inbred.blobFiltered.softMasked.braker.longestIso.protein","", strainCol)
ugh2 <- gsub(".20251014.inbred.blobFiltered.softMasked.braker.longestIso.protein","", ugh)
ugh3 <- gsub(".20251124.inbred.blobFiltered.softMasked.braker.longestIso.protein","", ugh2)
ugh4 <- gsub(".20251012.inbred.onlyONT.blobFiltered.softMasked.braker.longestIso.protein","", ugh3)
ugh5 <- gsub(".Nov2025.softMasked.braker.longest.protein","", ugh4)
ugh6 <- gsub(".20251012.inbred.withONT.blobFiltered.softMasked.braker.longestIso.protein","", ugh5)
strainCol_c2 <- gsub("c_elegans.PRJNA13758.WS283.csq.PCfeaturesOnly.longest.protein","N2", ugh6)
colnames(ortho_genes_dd) <- strainCol_c2

ortho_count <- ortho_genes_dd

# write.table("") ################################## WRITE CLEANED OG TSV TO ../../PROCESSED_DATA/ORHTOLOGY/ AND ADD AS A SUPPLEMENTARY TABLE

strainCol_c2_u <- strainCol_c2[!strainCol_c2 %in% c("Orthogroup")]

for (i in 1:length(strainCol_c2_u)) {
  print(paste0(i,"out of", length(strainCol_c2_u)))
  temp_colname = paste0(strainCol_c2_u[i], "_count")
  
  ortho_count <- ortho_count %>%
    dplyr::mutate(!!sym(temp_colname) := stringr::str_count(!!sym(strainCol_c2_u[i]),", ") + 1)
}

all_relations_pre <- ortho_count %>%
  dplyr::select(Orthogroup, dplyr::contains("_count"))


private_OGs <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/orthology/elegans/orthofinder/64_core/OrthoFinder/Results_Dec07/Orthogroups/Orthogroups_UnassignedGenes.tsv") %>%
  dplyr::filter(!grepl("MTCE",c_elegans.PRJNA13758.WS283.csq.PCfeaturesOnly.longest.protein))

colnames(private_OGs) <- strainCol_c2

private_cols <- strainCol_c2[!strainCol_c2 %in% c("Orthogroup")]

private_ortho_count <- private_OGs
for (i in 1:length(private_cols)) {
  print(paste0(i, " out of ", length(private_cols)))
  temp_colname <- paste0(private_cols[i], "_count")
  
  private_ortho_count <- private_ortho_count %>%
    dplyr::mutate(!!sym(temp_colname) := ifelse(is.na(!!sym(private_cols[i])), NA, 1))
}

all_relations_private <- private_ortho_count %>%
  dplyr::select(Orthogroup, dplyr::contains("_count"))

all_relations <- all_relations_pre %>%
  dplyr::bind_rows(all_relations_private)


#### Plotting ####
private_freq = (1/(length(strainCol_c2_u)))

classification <- all_relations %>%
  dplyr::mutate(across(2:(ncol(.)), ~ ifelse(. >= 1, 1, .))) %>%
  dplyr::mutate(sum = rowSums(across(-1, ~ ., .names = NULL), na.rm = TRUE)) %>%
  dplyr::mutate(freq = (sum / length(strainCol_c2_u))) %>%
  dplyr::mutate(
    class = case_when(
      freq == 1 ~ "core",
      freq > private_freq & freq < 1 ~ "accessory",
      freq == private_freq ~ "private",
      TRUE ~ "undefined"
    )
  ) %>%
  dplyr::count(freq, sum, class) %>%
  dplyr::mutate(percent = (n / sum(n)) * 100) 



gs_allOrtho <- ggplot(data = classification, aes(x = sum, y = n, fill = class)) + 
  geom_bar(stat = "identity", color = "black", alpha = 0.5) + 
  scale_fill_manual(values = c(
    "core" = "green4",
    "accessory" = "#DB6333",
    "private" = "magenta3"
  ), 
  limits = c("core", "accessory", "private"),  # Manually ordering legend items
  guide = guide_legend(title = NULL) 
  ) +
  ylab("Orthogroups") + # longest isoform
  xlab("Genomes") +
  # scale_y_continuous(labels = scales::percent_format(scale = 1)) + 
  # scale_x_continuous(labels = scales::percent_format(scale = 1)) +
  theme_classic() +
  scale_x_continuous(breaks = seq(0, max(classification$sum), by = 25)) +
  theme(
    axis.title = element_text(size = 24, color = 'black', face = 'bold'),
    legend.position = c(0.85, 0.8),
    plot.margin = margin(l = 20, r = 20, t = 20),
    # plot.title = element_text(size=26, face = 'bold', hjust=0.5),
    legend.text = element_text(size=22, color = 'black'),
    axis.text = element_text(size=18, color = 'black')
  )
gs_allOrtho

OG_class_count <- classification %>%
  dplyr::group_by(class) %>%
  dplyr::summarise(n_OG = sum(n)) %>%
  dplyr::ungroup()

# ggsave("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/plots/gene_set_allOrtho_142.png", gs_allOrtho, height = 5, width = 11, dpi = 600)


### PLOTTING BASED ON GENES ###
classification_genes <- all_relations %>%
  dplyr::mutate(across(2:(ncol(.)), ~ ifelse(. >= 1, 1, .))) %>%
  dplyr::mutate(sum = rowSums(across(-1, ~ ., .names = NULL), na.rm = TRUE)) %>%
  dplyr::mutate(freq = (sum / length(strainCol_c2_u))) %>%
  dplyr::mutate(
    class = case_when(
      freq == 1 ~ "core",
      freq > private_freq & freq < 1 ~ "accessory",
      freq == private_freq ~ "private",
      TRUE ~ "undefined"
    )
  ) 

sum_genes <- all_relations %>%
  dplyr::mutate(all_genes = rowSums(dplyr::across(2:ncol(.)), na.rm = TRUE)) %>%
  dplyr::select(all_genes) %>%
  dplyr::bind_cols(classification_genes) %>%
  dplyr::group_by(freq, sum, class) %>%
  dplyr::summarise(total_genes = sum(all_genes, na.rm = TRUE)) %>%
  dplyr::ungroup()


genes_allHOGs <- ggplot(data = sum_genes, aes(x = sum, y = total_genes / 1000, fill = class)) + 
  geom_bar(stat = "identity", color = "black", alpha = 0.5) + 
  scale_fill_manual(values = c(
    "core" = "green4",
    "accessory" = "#DB6333",
    "private" = "magenta3"
  ), 
  limits = c("core", "accessory", "private"),  # Manually ordering legend items
  guide = guide_legend(title = NULL) 
  ) +
  ylab("Genes (1e3)") + # longest isoform
  xlab("Strains") +
  # scale_y_continuous(labels = scales::percent_format(scale = 1)) + 
  # scale_x_continuous(labels = scales::percent_format(scale = 1)) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 16, color = 'black', face = 'bold'),
    legend.position = c(0.85, 0.8),
    # plot.margin = margin(l = 20),
    plot.title = element_text(size=18, face = 'bold', hjust=0.5),
    legend.text = element_text(size=16, color = 'black'),
    axis.text = element_text(size=14, color = 'black')
  )
genes_allHOGs

all_genes_class_count <- sum_genes %>%
  dplyr::group_by(class) %>%
  dplyr::summarise(n_genes = sum(total_genes)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(total = sum(n_genes)) %>%
  dplyr::mutate(prop = (n_genes / total) * 100) %>%
  dplyr::mutate(class = ifelse(class == 'accessory',"Accessory",
                               ifelse(class == "private","Private",
                                      ifelse(class == "core", "Core", class))))


### Pie chart of gene set proportion in pangenome
ggplot(all_genes_class_count, aes(x = 1, y = prop, fill = class)) +
  geom_col(width = 1, color = "black", alpha = 0.5) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = c(
    "Core" = "green4",
    "Accessory" = "#DB6333",
    "Private" = "magenta3")) +
  theme_void() +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 40, color = 'black', hjust = 0.5, vjust = -8)) +
  geom_label(
    data = subset(all_genes_class_count, class == "Core"),
    aes(x = 0.98, y = 32, label = paste0(round(prop,1), "% (", scales::comma(n_genes), ")")), size = 10) +
  geom_label(
    data = subset(all_genes_class_count, class == "Accessory"),
    aes(x = 1.0, y = 80, label = paste0(round(prop,1), "% (", scales::comma(n_genes), ")")), size = 10) +
  geom_label(
    data = subset(all_genes_class_count, class == "Private"),
    aes(x = 1.2, y = 0, label = paste0(round(prop,1), "% (", scales::comma(n_genes), ")")), size = 10) +
  labs(title = "Proportion of genes in pangenome")




####################################################################################################
####################################################################################################

# GENE SET RAREFACTION

####################################################################################################
####################################################################################################














