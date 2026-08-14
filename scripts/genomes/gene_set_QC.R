library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
library(ggrepel)
library(cowplot)

#######################################################################################################################################
# How many non-N2 genes do WSs have on average????
#######################################################################################################################################
# ======================================================================================================================================================================================== #
# Load all genes for N2 and wild strains
# ======================================================================================================================================================================================== #
genes_strain <- readr::read_tsv("../../processed_data/genome_resources/annotation/140Ws_CGC1_longestIsoGenes_BRAKER.tsv", col_names = c("seqid","source", "type", "start", "end", "score", "strand", "phase", "attributes", "strain"))
N2_gff <- ape::read.gff("../../processed_data/genome_resources/annotation/c_elegans.PRJNA13758.WS283.csq.PCfeaturesOnly.longest.gff3") %>% dplyr::mutate(strain="N2")
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


# ======================================================================================================================================================================================== #
# Load in orthology matrix
# ======================================================================================================================================================================================== #
all_relations <- readr::read_tsv("../../processed_data/orthology/OG_relations_matrix_count_withN2_BRAKER.tsv")

#######################################################################################################################################
### Looking at non-N2 ###
#######################################################################################################################################
nonRefGenes = as.data.frame(matrix(ncol = 3, nrow = 141))
colnames(nonRefGenes) <- c("strain","nonN2_genes",'nonN2_Orthogroups')

names_all <- colnames(all_relations %>% dplyr::select(-N2, -Orthogroup))

OG_list <- list()

for (i in 1:length(names_all)) {
  soi <- names_all[i]
  print(paste0("On strain: ", soi, ". ", i, "/141."))
  
  nonRefGenes[i,1] = soi
  
  non_ref <- all_relations %>% dplyr::select(N2, .data[[soi]]) %>% dplyr::filter(is.na(N2) & !is.na(.data[[soi]])) %>% dplyr::select(.data[[soi]]) %>% sum(na.rm = TRUE) 
  
  nonRefGenes[i,2] = non_ref
  
  non_ref_og_count <- all_relations %>% dplyr::select(N2, .data[[soi]]) %>% dplyr::filter(is.na(N2) & !is.na(.data[[soi]])) %>% nrow()
  non_ref_og <- all_relations %>% dplyr::select(Orthogroup, N2, .data[[soi]]) %>% dplyr::filter(is.na(N2) & !is.na(.data[[soi]])) %>% dplyr::pull(Orthogroup)
  OG_list[[i]] <- non_ref_og
  
  nonRefGenes[i,3] = non_ref_og_count
  
}

nonRefGenes_long_n2 <- nonRefGenes %>%
  pivot_longer(
    cols = c(nonN2_genes, nonN2_Orthogroups),
    names_to = "metric",
    values_to = "count"
  )

label_df_top_n2 <- nonRefGenes_long_n2 %>%
  dplyr::filter(metric == "nonN2_genes")  %>% 
  dplyr::arrange(desc(count)) %>% dplyr::slice_head(n = 10)

label_df_bottom_n2 <- nonRefGenes_long_n2 %>%
  dplyr::filter(metric == "nonN2_genes") %>% 
  dplyr::arrange(count) %>% dplyr::slice_head(n = 10)

labels_df_n2 <- label_df_top_n2 %>% dplyr::bind_rows(label_df_bottom_n2)

n2 <- ggplot(nonRefGenes_long_n2, aes(x = metric, y = count)) +
  geom_boxplot(aes(fill = metric), outlier.size = 0.6, width = 0.7, outlier.shape = NA, alpha = 0.5) +
  geom_line(aes(group = strain), alpha = 0.3) +
  geom_point(aes(group = strain),size = 1.5, alpha = 0.6) +
  geom_text_repel(data = labels_df_n2, aes(label = strain), size = 2, max.overlaps = 100) +
  labs(y = "Count") +
  scale_fill_manual(values = c("nonN2_genes" = "gray70", "nonN2_Orthogroups" = "gray30")) +
  theme(
    panel.background = element_blank(),
    legend.position = "none",
    panel.border = element_rect(color = 'black', fill = NA),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 10, color = 'black'),
    axis.text.y = element_text(size = 9, color = 'black'),
    axis.text.x = element_text(size = 9, color = 'black')) +
  coord_cartesian(ylim = c(800,5500))

# How many non-N2 genes on average?
av_non_N2 <- nonRefGenes_long_n2 %>% dplyr::filter(metric == "nonN2_genes") %>% dplyr::summarize(mean_nonN2 = mean(count))

# Looking at the proportion in each gene set
OG_vector_N2 <- unique(unlist(OG_list))

private_freq <- 1 /143

OG_classes <- all_relations %>%
  dplyr::mutate(across(2:(ncol(.)), ~ ifelse(. >= 1, 1, .))) %>%
  dplyr::mutate(sum = rowSums(across(-1, ~ ., .names = NULL), na.rm = TRUE)) %>%
  dplyr::mutate(freq = (sum / 143)) %>%
  dplyr::mutate(
    class = case_when(
      freq == 1 ~ "core",
      freq > private_freq & freq < 1 ~ "accessory",
      freq == private_freq ~ "private",
      TRUE ~ "undefined"
    )
  ) %>% dplyr::select(Orthogroup, class) 

OG_propN2 <- OG_classes %>%
  dplyr::filter(Orthogroup %in% OG_vector_N2) %>%
  dplyr::count(class, name = "class_count") %>%
  dplyr::mutate(prop = class_count / sum(class_count) * 100) %>%
  dplyr::mutate(class = factor(class, levels = c("accessory", "private"))) %>%
  dplyr::mutate(source = "N2")


#######################################################################################################################################
### Looking at non-N2 BRAKER
#######################################################################################################################################
names <- colnames(all_relations %>% dplyr::select(-N2_BRAKER,-Orthogroup))

nonRefGenes = as.data.frame(matrix(ncol = 3, nrow = 141))
colnames(nonRefGenes) <- c("strain","nonN2BRAKER_genes",'nonN2BRAKER_Orthogroups')

OG_listBRAKER_ <- list()

for (i in 1:length(names)) {
  soi <- names[i]
  print(paste0("On strain: ", soi, ". ", i, "/141."))
  
  nonRefGenes[i,1] = soi
  
  non_ref <- all_relations %>% dplyr::select(N2_BRAKER, .data[[soi]]) %>% dplyr::filter(is.na(N2_BRAKER) & !is.na(.data[[soi]])) %>% dplyr::select(.data[[soi]]) %>% sum(na.rm = TRUE) 
  
  nonRefGenes[i,2] = non_ref
  
  non_ref_og_count <- all_relations %>% dplyr::select(N2_BRAKER, .data[[soi]]) %>% dplyr::filter(is.na(N2_BRAKER) & !is.na(.data[[soi]])) %>% nrow()
  non_ref_og <- all_relations %>% dplyr::select(Orthogroup, N2_BRAKER, .data[[soi]]) %>% dplyr::filter(is.na(N2_BRAKER) & !is.na(.data[[soi]])) %>% dplyr::pull(Orthogroup)
  OG_listBRAKER_[[i]] <- non_ref_og
  
  nonRefGenes[i,3] = non_ref_og_count
  
}

nonRefGenes_long <- nonRefGenes %>%
  pivot_longer(
    cols = c(nonN2BRAKER_genes, nonN2BRAKER_Orthogroups),
    names_to = "metric",
    values_to = "count"
  )

label_df_top <- nonRefGenes_long %>%
  dplyr::filter(metric == "nonN2BRAKER_genes")  %>% 
  dplyr::arrange(desc(count)) %>% dplyr::slice_head(n = 10)

label_df_bottom <- nonRefGenes_long %>%
  dplyr::filter(metric == "nonN2BRAKER_genes") %>% 
  dplyr::arrange(count) %>% dplyr::slice_head(n = 10)

labels_df <- label_df_top %>% dplyr::bind_rows(label_df_bottom)

n2BRAKER_ <- ggplot(nonRefGenes_long, aes(x = metric, y = count)) +
  geom_boxplot(aes(fill = metric), outlier.size = 0.6, width = 0.7, outlier.shape = NA, alpha = 0.5) +
  geom_line(aes(group = strain), alpha = 0.3) +
  geom_point(aes(group = strain),size = 1.5, alpha = 0.6) +
  geom_text_repel(data = labels_df, aes(label = strain), size = 2, max.overlaps = 100) +
  labs(y = "Count") +
  scale_fill_manual(values = c("nonN2BRAKER_genes" = "gray70", "nonN2BRAKER_Orthogroups" = "gray30")) +
  theme(
    panel.background = element_blank(),
    legend.position = "none",
    panel.border = element_rect(color = 'black', fill = NA),
    axis.title.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 9, color = 'black')) +
  coord_cartesian(ylim = c(800,5500))


n2_n2Braker <- cowplot::plot_grid(
  n2, n2BRAKER_,
  nrow = 1)

# Looking at the proportion in each gene set?
OG_vector_N2BRAKER_ <- unique(unlist(OG_listBRAKER_))

OG_prop <- OG_classes %>%
  dplyr::filter(Orthogroup %in% OG_vector_N2BRAKER_) %>%
  dplyr::count(class, name = "class_count") %>%
  dplyr::mutate(prop = class_count / sum(class_count) * 100) %>%
  dplyr::mutate(class = factor(class, levels = c("accessory", "private"))) %>%
  dplyr::mutate(source = "N2BRAKER_")

both <- OG_prop %>% dplyr::bind_rows(OG_propN2) %>% dplyr::mutate(source = factor(source, levels = c("N2BRAKER_", "N2")))

both_prop <- ggplot(both, aes(x = prop, y = source, fill = class)) +
  geom_col(position = "stack", width = 0.5, color = "black") +
  scale_fill_manual(values = c("accessory" = "#DB6333", "private" = "magenta3")) +
  labs(x = "Proportion (%)", y = NULL, fill = "Gene set") +
  theme(
    panel.background = element_blank(),
    legend.position = "none",
    panel.border = element_rect(color = 'black', fill = NA),
    axis.title.x = element_text(size = 9, color = 'black'),
    axis.ticks.y = element_blank(),
    axis.title.y = element_text(size = 10, color = 'black'),
    axis.text.y = element_text(size = 9, color = 'black'),
    axis.text.x = element_text(size = 9, color = 'black'))
both_prop


all_plt <- cowplot::plot_grid(
  n2_n2Braker, both_prop,
  nrow = 2,
  rel_heights = c(5,1),
  labels = c("a","b"))
all_plt

ggsave("../../figures/supplementary/gene_model_QC.png", all_plt, width = 7.5, height = 7.5, dpi = 600)


#######################################################################################################################################
### Looking at N2-specific genes that are not found in any wild strains
#######################################################################################################################################
names_all <- colnames(all_relations %>% dplyr::select(-N2, -Orthogroup))

refgenes = as.data.frame(matrix(ncol = 3, nrow = 141))
colnames(refgenes) <- c("strain","N2_specific_genes",'N2_specific_Orthogroups')

OG_list <- list()

for (i in 1:length(names_all)) {
  soi <- names_all[i]
  print(paste0("On strain: ", soi, ". ", i, "/142."))
  
  refgenes[i,1] = soi
  
  non_ref <- all_relations %>% dplyr::select(N2, .data[[soi]]) %>% dplyr::filter(!is.na(N2) & is.na(.data[[soi]])) %>% dplyr::select(N2) %>% sum(na.rm = TRUE) 
  
  refgenes[i,2] = non_ref
  
  non_ref_og_count <- all_relations %>% dplyr::select(N2, .data[[soi]]) %>% dplyr::filter(!is.na(N2) & is.na(.data[[soi]])) %>% nrow()
  non_ref_og <- all_relations %>% dplyr::select(Orthogroup, N2, .data[[soi]]) %>% dplyr::filter(!is.na(N2) & is.na(.data[[soi]])) %>% dplyr::pull(Orthogroup)
  OG_list[[i]] <- non_ref_og
  
  refgenes[i,3] = non_ref_og_count
  
}

nonRefGenes_long <- refgenes %>%
  pivot_longer(
    cols = c(N2_specific_genes, N2_specific_Orthogroups),
    names_to = "metric",
    values_to = "count")

label_df_top <- nonRefGenes_long %>%
  dplyr::filter(metric == "N2_specific_genes")  %>% 
  dplyr::arrange(desc(count)) %>% dplyr::slice_head(n = 10)

label_df_bottom <- nonRefGenes_long %>%
  dplyr::filter(metric == "N2_specific_genes") %>% 
  dplyr::arrange(count) %>% dplyr::slice_head(n = 10)

labels_df <- label_df_top %>% dplyr::bind_rows(label_df_bottom)

n2_spec <- ggplot(nonRefGenes_long, aes(x = metric, y = count)) +
  geom_boxplot(aes(fill = metric), outlier.size = 0.6, width = 0.7, outlier.shape = NA, alpha = 0.5) +
  geom_line(aes(group = strain), alpha = 0.3) +
  geom_point(aes(group = strain),size = 1.5, alpha = 0.6) +
  geom_text_repel(data = labels_df, aes(label = strain), size = 3, max.overlaps = 100) +
  labs(y = "Count") +
  scale_fill_manual(values = c("N2_specific_genes" = "orange", "N2_specific_Orthogroups" = "#DB6333")) +
  theme(
    panel.background = element_blank(),
    legend.position = "none",
    panel.border = element_rect(color = 'black', fill = NA),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 12, color = 'black'),
    axis.text.y = element_text(size = 10, color = 'black'),
    axis.text.x = element_text(size = 12, color = 'black')) 
n2_spec

# Looking at proprotion in each gene set
OG_vector_N2_spec <- unique(unlist(OG_list))

OG_classes_spec <- all_relations %>%
  dplyr::mutate(across(2:(ncol(.)), ~ ifelse(. >= 1, 1, .))) %>%
  dplyr::mutate(sum = rowSums(across(-1, ~ ., .names = NULL), na.rm = TRUE)) %>%
  dplyr::mutate(freq = (sum / 143)) %>%
  dplyr::mutate(
    class = case_when(
      freq == 1 ~ "core",
      freq > private_freq & freq < 1 ~ "accessory",
      freq == private_freq ~ "private",
      TRUE ~ "undefined"
    )
  ) %>% dplyr::select(Orthogroup, class) %>%
  dplyr::filter(Orthogroup %in% OG_vector_N2_spec)

OG_propN2_spec <- OG_classes_spec %>%
  dplyr::count(class, name = "class_count") %>%
  dplyr::mutate(prop = class_count / sum(class_count) * 100) %>%
  dplyr::mutate(class = factor(class, levels = c("accessory", "core", "private"))) %>%
  dplyr::mutate(source = "N2")

n2_spec_prop <- ggplot(OG_propN2_spec, aes(x = prop, y = "", fill = class)) +
  geom_col(position = "stack", width = 0.5, color = "black") +
  scale_fill_manual(values = c("accessory" = "#DB6333", "private" = "magenta3", "core" = "green4")) +
  labs(x = "Proportion (%)", y = NULL, fill = "Gene set") +
  theme(
    panel.background = element_blank(),
    legend.position = "none",
    panel.border = element_rect(color = 'black', fill = NA),
    axis.title.x = element_text(size = 11, color = 'black'),
    axis.ticks.y = element_blank(),
    axis.title.y = element_text(size = 12, color = 'black'),
    axis.text.y = element_text(size = 10, color = 'black'),
    axis.text.x = element_text(size = 10, color = 'black')) 

n2_specific_genes <- cowplot::plot_grid(
  n2_spec, n2_spec_prop,
  nrow = 2,
  rel_heights = c(6,1))
n2_specific_genes

# On average, how many N2-specific genes does N2 have?
av_N2_spec <- nonRefGenes_long %>% dplyr::filter(metric == "N2_specific_genes") %>% dplyr::summarize(mean_n2_spec = mean(count))

ggsave("../../figures/supplementary/N2_specific_gene_model_QC.png", n2_specific_genes, width = 7.5, height = 7.5, dpi = 600)




