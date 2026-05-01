library(readr)
library(dplyr)
library(ggplot2)
library(stringr)
library(ape)
library(tidyr)
library(data.table)
library(cowplot)
library(ggforce)

####################################################################################################
####################################################################################################

# GENE SET CLASSIFICATION AND VISUALIZATION

####################################################################################################
####################################################################################################
# ======================================================================================================================================================================================== #
# Pulling genes for all WSs and N2 
# ======================================================================================================================================================================================== #
genes_strain <- readr::read_tsv("../../processed_data/genome_resources/annotation/140Ws_CGC1_longestIsoGenes_BRAKER.tsv", col_names = c("seqid","source", "type", "start", "end", "score", "strand", "phase", "attributes", "strain")) %>% dplyr::filter(strain != "ECA396")
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
# OF matrix manipulation and plotting #
# ======================================================================================================================================================================================== #
# OrthoFinder-assigned orthogroups among 142 strains of C. elegans
ortho_genes_dd <- readr::read_tsv("../../processed_data/orthology/orthofinder/orthofinder_output/Orthogroups.tsv") %>%
  dplyr::filter(!grepl("MTCE",c_elegans.PRJNA13758.WS283.csq.PCfeaturesOnly.longest.protein))

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


private_OGs <- readr::read_tsv("../../processed_data/orthology/orthofinder/orthofinder_output/Orthogroups_UnassignedGenes.tsv") %>%
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


# Plotting OG gnee set histogram
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
      TRUE ~ "undefined")) %>%
  dplyr::count(freq, sum, class) %>%
  dplyr::mutate(percent = (n / sum(n)) * 100) 

geneset_hist <- ggplot(data = classification, aes(x = sum, y = n, fill = class)) + 
  geom_bar(stat = "identity", color = "black", alpha = 0.5) + 
  scale_fill_manual(values = c("core" = "green4", "accessory" = "#DB6333", "private" = "magenta3"), limits = c("core", "accessory", "private"), guide = guide_legend(title = NULL)) +
  ylab("Orthogroups") + 
  xlab("Genomes") +
  theme_classic() +
  scale_y_continuous(breaks = seq(0, 17000, by = 4000), expand = c(0,0)) +
  coord_cartesian(ylim = c(0,17000)) +
  scale_x_continuous(breaks = seq(0, max(classification$sum), by = 25), expand = c(0,0)) +
  theme(
    axis.title = element_text(size = 10, color = 'black'),
    legend.position = "none",
    # legend.position = c(0.85, 0.8),
    # legend.text = element_text(size=10, color = 'black'),
    axis.text = element_text(size=10, color = 'black')) 
geneset_hist

# OG_class_count <- classification %>%
#   dplyr::group_by(class) %>%
#   dplyr::summarise(n_OG = sum(n)) %>%
#   dplyr::ungroup()

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
      TRUE ~ "undefined")) 

sum_genes <- all_relations %>%
  dplyr::mutate(all_genes = rowSums(dplyr::across(2:ncol(.)), na.rm = TRUE)) %>%
  dplyr::select(all_genes) %>%
  dplyr::bind_cols(classification_genes) %>%
  dplyr::group_by(freq, sum, class) %>%
  dplyr::summarise(total_genes = sum(all_genes, na.rm = TRUE)) %>%
  dplyr::ungroup()

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
all_genes_class_count <- all_genes_class_count %>%
  dplyr::arrange(desc(n_genes)) %>%
  dplyr::mutate(class = factor(class, levels = c("Core","Accessory","Private"))) %>%
  dplyr::mutate(end = cumsum(prop) / 100 * 2 * pi, start = lag(end, default = 0))

pie_base <- ggplot(all_genes_class_count) +
  geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0, r = 0.9, start = start, end = end, fill = class), color = "black", alpha = 0.5) +
  coord_fixed() +
  scale_fill_manual(values = c(
    "Core" = "green4",
    "Accessory" = "#DB6333",
    "Private" = "magenta3")) +
  theme_void() +
  theme(legend.position = "none")

# Adding labels to pie chart
pie_gene_set <- ggdraw(pie_base) +
  draw_label(paste0("Core\n", scales::comma(all_genes_class_count$n_genes[all_genes_class_count$class == "Core"]), " (", 
                    round(all_genes_class_count$prop[all_genes_class_count$class == "Core"], 1), "%)"), 
             x = 0.89, y = 0.2, size = 10, color = 'green4') +
  draw_label(paste0("Accessory\n", scales::comma(all_genes_class_count$n_genes[all_genes_class_count$class == "Accessory"]), " (", 
                    round(all_genes_class_count$prop[all_genes_class_count$class == "Accessory"], 1), "%)"), 
             x = 0.08, y = 0.79, size = 10, color = '#DB6333') +
  draw_label(paste0("Private\n", scales::comma(all_genes_class_count$n_genes[all_genes_class_count$class == "Private"]), " (", 
                    round(all_genes_class_count$prop[all_genes_class_count$class == "Private"], 1), "%)"), 
             x = 0.5, y = 1.05, size = 10, color = 'magenta3')
pie_gene_set


HIST <- ggdraw() +
  draw_plot(geneset_hist) +
  draw_plot(pie_gene_set, 
            x = 0.26,      # horizontal position (0-1)
            y = 0.4,      # vertical position (0-1)
            width = 0.5,  # width of inset (0-1)
            height = 0.5)  # height of inset (0-1)


####################################################################################################
####################################################################################################

# GENE SET RAREFACTION

####################################################################################################
####################################################################################################
all <- all_relations %>% dplyr::select(-Orthogroup)

##################### REMOVE THESE BEFORE PUBLISHING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
pan_final <- readRDS("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/orthology/elegans/plots/pan_iterativeOGcount.rds")
core_final <- readRDS("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/orthology/elegans/plots/core_iterativeOGcount.rds")
accessory_final <- readRDS("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/orthology/elegans/plots/acc_iterativeOGcount.rds")
priv_final <- readRDS("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/orthology/elegans/plots/priv_iterativeOGcount.rds")

set.seed(42)

n_strains_total <- ncol(all)
n_perms <- 100

# # For the pangenome
# pan_list <- vector("list", length = n_strains_total - 1) # -1 because we iterate 2 - 142
# iteration_pan <- 1
#
# for (i in 2:n_strains_total) {
#   for (it_i in 1:n_perms) {
#     # pick k random strains
#     cols <- sample(colnames(all), size = i, replace = FALSE)
#     subset <- all[, cols, drop = FALSE] # subset all df to k strains
# 
#     subset <- subset %>% dplyr::filter(!if_all(everything(), is.na))
#     all_OGs <- nrow(subset)
#     # print(all_OGs)
#     print(paste0("On strain subset: ",i,", and iteration: ", it_i))
# 
#     pan_list[[iteration_pan]] <- data.frame(
#       n_strains = i,
#       replicate = it_i,
#       n_core_ogs = all_OGs)
#     iteration_pan <- iteration_pan + 1
#   }
# }
# 
# pan_final <- dplyr::bind_rows(pan_list)

pan_summary <- pan_final %>%
  dplyr::group_by(n_strains) %>%
  dplyr::summarise(
    median_core = median(n_core_ogs),
    mean_core   = mean(n_core_ogs),
    sd_core     = sd(n_core_ogs),
    q05         = quantile(n_core_ogs, 0.05),
    q95         = quantile(n_core_ogs, 0.95)) %>%
  dplyr::ungroup()

# # For the core pangenome
# res_list <- vector("list", length = n_strains_total - 1) # -1 because we iterate 2 - 141
# iteration <- 1
# 
# for (i in 2:n_strains_total) {
#   for (it_i in 1:n_perms) {
#     # pick k random strains
#     cols <- sample(colnames(all), size = i, replace = FALSE)
#     subset <- all[, cols, drop = FALSE] # subset all df to k strains
# 
#     # binarize and count core
#     core_calc <- subset %>%
#       dplyr::mutate(across(everything(), ~ ifelse(is.na(.),0, ifelse(. >= 1, 1, .)))) %>%
#       dplyr::mutate(sum = rowSums(across(everything()))) %>%
#       dplyr::mutate(freq = (sum / i)) %>%
#       dplyr::mutate(class = case_when(freq == 1 ~ "core")) %>%
#       dplyr::filter(class == "core")
# 
#     # print(head(core_calc))
#     print(paste0("On strain subset: ", i,", and iteration: ", it_i))
# 
#     core_count <- nrow(core_calc)
#     # print(core_count)
# 
#     res_list[[iteration]] <- data.frame(
#       n_strains = i,
#       replicate = it_i,
#       n_core_ogs = core_count)
#     iteration <- iteration + 1
#   }
# }
# 
# core_final <- dplyr::bind_rows(res_list)

core_summary <- core_final %>%
  dplyr::group_by(n_strains) %>%
  dplyr::summarise(
    median_core = median(n_core_ogs),
    mean_core   = mean(n_core_ogs),
    sd_core     = sd(n_core_ogs),
    q05         = quantile(n_core_ogs, 0.05),
    q95         = quantile(n_core_ogs, 0.95)) %>%
  dplyr::ungroup()


# For the accessory pangenome
# res_list <- vector("list", length = n_strains_total - 1) # -1 because we iterate 2 - 141
# iteration <- 1
# 
# for (i in 2:n_strains_total) {
#   for (it_i in 1:n_perms) {
#     # pick k random strains
#     cols <- sample(colnames(all), size = i, replace = FALSE)
#     subset <- all[, cols, drop = FALSE] # subset all df to k strains
# 
#     # binarize and count accessory
#     accessory_calc <- subset %>%
#       dplyr::mutate(across(everything(), ~ ifelse(is.na(.),0, ifelse(. >= 1, 1, .)))) %>%
#       dplyr::mutate(sum = rowSums(across(everything()))) %>%
#       dplyr::mutate(freq = (sum / i)) %>%
#       dplyr::mutate(
#         class = case_when(
#           freq == 1 ~ "core",
#           freq > 1/i & freq < 1 ~ "accessory",
#           freq == 1/i ~ "private",
#           TRUE ~ "undefined")) %>%
#       dplyr::filter(class == "accessory")
# 
#     # print(head(accessory_calc))
#     print(paste0("On strain subset: ", i,", and iteration: ", it_i))
# 
#     accessory_count <- nrow(accessory_calc)
#     # print(accessory_count)
# 
#     res_list[[iteration]] <- data.frame(
#       n_strains = i,
#       replicate = it_i,
#       n_accessory_ogs = accessory_count)
#     iteration <- iteration + 1
#   }
# }
# 
# accessory_final <- dplyr::bind_rows(res_list)

accessory_summary <- accessory_final %>%
  dplyr::group_by(n_strains) %>%
  dplyr::summarise(
    median_accessory = median(n_accessory_ogs),
    mean_accessory   = mean(n_accessory_ogs),
    sd_accessory     = sd(n_accessory_ogs),
    q05         = quantile(n_accessory_ogs, 0.05),
    q95         = quantile(n_accessory_ogs, 0.95)) %>%
  dplyr::ungroup()

# For the private pangenome
# res_list <- vector("list", length = n_strains_total -1)
# iteration <- 1
# 
# for (i in 2:n_strains_total) {              # NEED TO BEGIN ITERATION ON A SINLE STRAIN WHEN CALCULATING PRIVATE
#   for (it_i in 1:n_perms) {
#     # pick k random strains
#     cols <- sample(colnames(all), size = i, replace = FALSE)
#     subset <- all[, cols, drop = FALSE] # subset all df to k strains
# 
#     # binarize and count accessory
#     priv_calc <- subset %>%
#       dplyr::mutate(across(everything(), ~ ifelse(is.na(.),0, ifelse(. >= 1, 1, .)))) %>%
#       dplyr::mutate(sum = rowSums(across(everything()))) %>%
#       dplyr::mutate(freq = (sum / i)) %>%
#       dplyr::mutate(
#         class = case_when(
#           freq == 1 ~ "core",
#           freq > 1/i & freq < 1 ~ "accessory",
#           freq == 1/i ~ "private",
#           TRUE ~ "undefined")) %>%
#       dplyr::filter(class == "private")
# 
#     # print(head(priv_calc))
#     print(paste0("On strain subset: ", i,", and iteration: ", it_i))
# 
#     priv_count <- nrow(priv_calc)
#     # print(priv_count)
# 
#     res_list[[iteration]] <- data.frame(
#       n_strains = i,
#       replicate = it_i,
#       n_priv_ogs = priv_count)
#     iteration <- iteration + 1
#   }
# }
# 
# priv_final <- dplyr::bind_rows(res_list)

priv_summary <- priv_final %>%
  dplyr::group_by(n_strains) %>%
  dplyr::summarise(
    median_priv = median(n_priv_ogs),
    mean_priv   = mean(n_priv_ogs),
    sd_priv     = sd(n_priv_ogs),
    q05         = quantile(n_priv_ogs, 0.05),
    q95         = quantile(n_priv_ogs, 0.95)) %>%
  dplyr::ungroup()

# Plotting the final rarefaction curves for the pangenome and core gene set
RAREFACT <- ggplot() +
  # Pangenome
  geom_errorbar(data = pan_summary, aes(x = n_strains, ymin = mean_core - sd_core, ymax = mean_core + sd_core), width = 0.5, alpha = 0.5) +
  geom_point(data = pan_summary, aes(x = n_strains, y = mean_core, color = "Pangenome"), size = 0.5) +
  # Core
  geom_errorbar(data = core_summary, aes(x = n_strains, ymin = mean_core - sd_core, ymax = mean_core + sd_core), width = 0.5, alpha = 0.5) +
  geom_point(data = core_summary, aes(x = n_strains, y = mean_core, color = "Core"), size = 0.5) +
  scale_color_manual(
    values = c("Pangenome" = "blue", "Core" = "green4"),
    limits = c("Pangenome", "Core")) +
  guides(color = guide_legend(override.aes = list(size = 4))) +  labs(x = "Genomes", y = "Orthogroups") +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA, color = "black"),
    axis.title = element_text(size = 10),
    legend.position = 'inside',
    plot.margin = margin(l = 2, t = 30, r = 2, b = 2),
    legend.position.inside = c(0.25,0.87),
    legend.title = element_blank(),
    legend.text = element_text(size = 10, color = 'black'),
    axis.text = element_text(size = 10, color = 'black')) +
  coord_cartesian(ylim = c(0,60000)) 


# Plotting all gene set rarefaction curve
rarefact_supplement_all <- ggplot() +
  # Pangenome
  geom_errorbar(data = pan_summary, aes(x = n_strains, ymin = mean_core - sd_core, ymax = mean_core + sd_core), width = 0.5, alpha = 0.5) +
  geom_point(data = pan_summary, aes(x = n_strains, y = mean_core, color = "Pangenome"), size = 1) +
  # Core
  geom_errorbar(data = core_summary, aes(x = n_strains, ymin = mean_core - sd_core, ymax = mean_core + sd_core), width = 0.5, alpha = 0.5) +
  geom_point(data = core_summary, aes(x = n_strains, y = mean_core, color = "Core"), size = 1) +
  # Accessory
  geom_errorbar(data = accessory_summary, aes(x = n_strains, ymin = mean_accessory - sd_accessory, ymax = mean_accessory + sd_accessory), width = 0.5, alpha = 0.5) +
  geom_point(data = accessory_summary, aes(x = n_strains, y = mean_accessory, color = "Accessory"), size = 1) +
  # Private
  geom_errorbar(data = priv_summary, aes(x = n_strains, ymin = mean_priv - sd_priv, ymax = mean_priv + sd_priv), width = 0.5, alpha = 0.5) +
  geom_point(data = priv_summary, aes(x = n_strains, y = mean_priv, color = "Private"), size = 1) +
  scale_color_manual(
    values = c("Pangenome" = "blue", "Core" = "green4", "Accessory" = "#DB6333", "Private" = "magenta3"),
    limits = c("Pangenome", "Core", "Accessory", "Private")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +  labs(x = "Genomes", y = "Orthogroups") +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA, color = "black"),
    axis.title = element_text(size = 14),
    legend.position = 'inside',
    legend.position.inside = c(0.15,0.9),
    legend.title = element_blank(),
    legend.text = element_text(size = 14, color = 'black'),
    axis.text = element_text(size = 12, color = 'black')) +
  coord_cartesian(ylim = c(0, 60000)) 

# ggsave("../../figures/supplementary/pangenome_rarefaction.png", rarefact_supplement_all, height = 7, width = 7.5, dpi = 600)

####################################################################################################
####################################################################################################

# N2 gene set density

####################################################################################################
####################################################################################################
count <- all_relations %>%
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

N2_coords <- N2_gff %>%
  dplyr::filter(type == "mRNA") %>%
  dplyr::mutate(attributes = gsub("ID=transcript:","",attributes), attributes = gsub("Parent=gene:","",attributes)) %>%
  tidyr::separate_wider_delim(attributes, delim = ";",names = c("tran", "gene", "rest"), too_many = "merge") %>%
  dplyr::select(seqid, start, end, tran, gene, -rest)

n2_gene <- N2_coords %>%
  dplyr::mutate(tran = paste0("transcript_",tran))

n2_table <- ortho_genes_dd %>%
  dplyr::bind_rows(private_OGs) %>% 
  dplyr::select(Orthogroup,N2)

ortho_count_wCoord <- count %>%
  dplyr::left_join(n2_table, by = "Orthogroup") %>%
  dplyr::select(freq, class, N2) %>%
  dplyr::filter(!is.na(N2)) %>%
  tidyr::separate_rows(N2, sep = ",\\s*") %>% # splitting rows so each gene is on a row and it retains is gene set classification
  dplyr::left_join(n2_gene, by = c("N2" = "tran")) %>%
  dplyr::select(freq,class,N2,gene,seqid,start,end)


plt_data <- ortho_count_wCoord %>%
  dplyr::filter(seqid != "MtDNA") %>%
  dplyr::mutate(mid_mb = (start + end) / 2 / 1e6)

ggplot(data = plt_data) +
  geom_rect(aes(xmin = start / 1e6, xmax = end / 1e6, fill = class), ymin = -Inf, ymax = Inf, alpha = 0.5) +
  geom_density(aes(x = mid_mb, y = after_stat(scaled), color = class), adjust = 0.5, linewidth = 0.9, position = "identity") +
  scale_color_manual(values = c(accessory="#DB6333", private="magenta3", core="green4")) +
  scale_fill_manual(values = c("accessory" = "#DB6333", "private" = "magenta3", "core" = "green4")) +
  facet_wrap(~seqid, scales = "free") +
  theme(axis.title.x = element_text(size = 16, color = 'black', face = 'bold'),
        axis.text.x = element_text(size = 13, color = 'black'),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = 'black')) +
  scale_x_continuous(expand = c(0,0)) +
  xlab("N2 genome position (Mb)") +
  scale_y_continuous(NULL, breaks = NULL) 


plt_data <- plt_data %>%
  mutate(seqid = factor(seqid, levels = c("I","II","III","IV","V","X"))) %>%
  dplyr::rename(Class = class)
ncol_facets <- 3
lev <- levels(plt_data$seqid)
left_facets <- lev[seq(1, length(lev), by = ncol_facets)] 

class_labels <- plt_data %>% dplyr::distinct(Class) %>% dplyr::mutate(y = ifelse(Class == "core", 1.75,
                                                                                 ifelse(Class == "accessory", 1, 0.25)), x = -Inf)
class_labels <- tidyr::crossing(class_labels, seqid = left_facets)

ggplot() + 
  geom_rect(data = plt_data %>% dplyr::filter(Class == "core"), aes(xmin = start / 1e6, xmax = end / 1e6, ymin = 1.5, ymax = 2), fill = "green4", alpha = 0.5) +
  geom_rect(data = plt_data %>% dplyr::filter(Class == "accessory"), aes(xmin = start / 1e6, xmax = end / 1e6, ymin = 0.75, ymax = 1.25), fill = "#DB6333", alpha = 0.5) +
  geom_rect(data = plt_data %>% dplyr::filter(Class == "private"), aes(xmin = start / 1e6, xmax = end / 1e6, ymin = 0, ymax = 0.5), fill = "magenta3", alpha = 0.5) +
  geom_text(data = class_labels,aes(x = x, y = y, label = Class), hjust = 1.1, fontface = "bold", size = 4, inherit.aes = FALSE) +
  facet_wrap(~seqid, scales = "free") +
  coord_cartesian(clip = "off") +
  theme(axis.title.x = element_text(size = 16, color = 'black', face = 'bold'),
        axis.text.x = element_text(size = 13, color = 'black'),
        panel.background = element_blank(),
        panel.border = element_rect(fill = "NA", color = 'black'),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        plot.margin = margin(l = 70)) +
  scale_x_continuous(expand = c(0,0)) +
  xlab("N2 genome position (Mb)") +
  scale_y_continuous(NULL, breaks = NULL)


ggplot(data = plt_data) +
  geom_density(aes(x = mid_mb, y = after_stat(count), color = Class), adjust = 1, linewidth = 0.9, position = "identity") + 
  scale_color_manual(values = c(accessory="#DB6333", private="magenta3", core="green4")) +
  facet_wrap(~seqid, scales = 'free_x') +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        panel.background = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = 'none',
        panel.border = element_rect(fill = "NA", color = 'black'),
        # strip.background = element_blank(),
        axis.text.y = element_text(size = 12, color = 'black'),
        axis.title.y = element_text(size =14, color = 'black', face = 'bold'),
        # strip.text.x = element_blank(),
        plot.margin = margin(l = 70)) +
  ylab("Kernel Density * Count")


# Final plot
lane_h   <- 14        # lane height (y units)
gap_h    <- 6         # gap between lanes
y_core   <- -lane_h                      # [-12, 0)
y_acc    <- -(2*lane_h + gap_h)          # [-30, -18)
y_priv   <- -(3*lane_h + 2*gap_h)        # [-48, -36)

rects <- plt_data %>%
  dplyr::mutate(Class = ifelse(Class == "core", "Core",
                               ifelse(Class == "accessory","Acc.",
                                      ifelse(Class == "private","Priv.", Class)))) %>%
  dplyr::mutate(xmin = start/1e6, xmax = end/1e6, ymin = dplyr::case_when(
    Class == "Core" ~ y_core - lane_h, 
    Class == "Acc." ~ y_acc - lane_h, 
    Class == "Priv." ~ y_priv - lane_h), 
    ymax = dplyr::case_when(
      Class == "Core" ~ y_core, 
      Class == "Acc." ~ y_acc, 
      Class == "Priv." ~ y_priv)) %>% 
  dplyr::rename(`Gene set` = Class)

class_labels_final <- plt_data %>% 
  dplyr::mutate(Class = ifelse(Class == "core", "Core",
                               ifelse(Class == "accessory","Acc.",
                                      ifelse(Class == "private","Priv.", Class)))) %>%
  dplyr::distinct(Class) %>% dplyr::mutate(y = ifelse(Class == "Core", -21,
                                                      ifelse(Class == "Acc.", -41, -61)), x = -Inf)

class_labels_final <- tidyr::crossing(class_labels_final, seqid = left_facets)

density_label <- plt_data %>% dplyr::mutate(label = "Gene density") %>% dplyr::mutate(x = -Inf, y = 100) %>% dplyr::select(label,x,y)
density_label_final <- tidyr::crossing(density_label, seqid = left_facets)

DENSITY <- ggplot() +
  # draw rect "tracks" first so they sit under the density
  geom_rect(data = rects, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = `Gene set`), color = NA) +
  scale_fill_manual(values = c("Acc."="#DB6333", "Priv."="magenta3", "Core"="green4")) +
  # density scaled by counts (absolute abundance)
  geom_density(data = plt_data, aes(x = mid_mb, y = after_stat(count), color = Class), adjust = 0.5, linewidth = 0.4, position = "identity", show.legend = FALSE, alpha = 0.5) +
  geom_text(data = class_labels_final, aes(x = x, y = y, label = Class), hjust = 1.1, size = 2.5, inherit.aes = FALSE) +
  geom_text(data = density_label_final, aes(x = x, y = y, label = label), vjust = -1, size = 3, inherit.aes = FALSE, angle = 90) +
  scale_color_manual(values = c(accessory="#DB6333", private="magenta3", core="green4")) +
  facet_wrap(~ seqid, scales = "free_x", ncol = 3) +
  coord_cartesian(ylim = c(y_priv - 5, NA), clip = "off") +
  scale_x_continuous(expand = c(0.01,0), breaks = seq(0, 25, by = 5)) +
  labs(x = "N2 genome position (Mb)") +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA, color = "black"),
    strip.text = element_text(size = 10, color = 'black'),
    strip.background = element_blank(),
    axis.title.x = element_text(size = 10),
    axis.text = element_blank(),
    axis.ticks.y = element_blank(),
    # panel.spacing=unit(2, "lines"),
    axis.title.y = element_blank(),
    legend.position = 'none',
    plot.margin  = margin(l = 30, r = 5, t = 5, b = 5))
DENSITY

# Looking at bandwith per gene set
# bw_by_class <- plt_data %>%
#   dplyr::group_by(Class) %>%
#   dplyr::summarise(bw = density(mid_mb)$bw)
# bw_by_class
# 
# class_count <- rects %>% dplyr::group_by(`Gene set`) %>% dplyr::mutate(number = n()) %>% dplyr::distinct(`Gene set`,number)


# Concatenating plots for bottom row of final plot
BOTTOM <- cowplot::plot_grid(
  RAREFACT, DENSITY,
  nrow = 1,
  rel_widths = c(1,1.2),
  labels = c("b","c"))
BOTTOM

# Creating final plot
final_plot <- cowplot::plot_grid(
  HIST, BOTTOM,
  nrow = 2,
  rel_heights = c(1.25,1),
  labels = "a")
final_plot

# Save the plot
ggsave("../../figures/gene_set_visualization.png", final_plot, height = 8, width = 7.5, dpi = 600)


