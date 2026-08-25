library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(enrichplot)
library(cowplot)
library(stringr)
library(data.table)
library(rnaturalearth)
library(rnaturalearthdata)

# ======================================================================================================================================================================================== #
# Loading in all genes in the pangenome with resolution on the gene sets they contribute to
# ======================================================================================================================================================================================== #
long_class <- readr::read_tsv("../../processed_data/orthology/all_genes_class_OGs.tsv") %>%
  dplyr::mutate(gene = paste0(strain,"_",gene))

WSs <- long_class %>% 
  dplyr::distinct(strain) %>%
  dplyr::filter(strain != "N2" & strain != "CGC1") %>%
  dplyr::pull()

# ======================================================================================================================================================================================== #
# IPR results
# ======================================================================================================================================================================================== #
all_ipr <- readr::read_tsv("../../tables/IPR_annotation_142strains.tsv", col_names = c("tran", "MD5_digest", "seq_length", "app", "signature_accession", "signature_description", "start", "end", "score", "status", "date", "IPR_accession","IPR_description","GO", "pathways")) %>%
  dplyr::filter(!grepl("CGC1_", tran), !grepl("N2_", tran)) %>% # using IPR results for wild strains only 
  dplyr::select(tran, IPR_accession, IPR_description, GO)

ws_hdr_genes <- readr::read_tsv("../../processed_data/hdr_liftover/wild_strain_genes_inHDRs.tsv", col_names = c("strain", "gene", "class")) %>%
  dplyr::mutate(gene = paste0(strain, "_", gene))

ws_hdr_priv_genes <- ws_hdr_genes %>% dplyr::filter(class == "private") %>% dplyr::pull(gene)
ws_hdr_acc_genes <- ws_hdr_genes %>% dplyr::filter(class == "accessory") %>% dplyr::pull(gene)
ws_hdr_core_genes <- ws_hdr_genes %>% dplyr::filter(class == "core") %>% dplyr::pull(gene)
all_ws_hdr_genes <- ws_hdr_genes %>% dplyr::pull(gene)

all_ipr_background <- all_ipr %>% 
  dplyr::filter(!is.na(IPR_description) & IPR_accession != "-") %>%
  dplyr::distinct(tran, IPR_accession, IPR_description) %>% # only one type of IPR annotation per gene
  dplyr::group_by(IPR_accession) %>%
  dplyr::mutate(n_IPR_acc_background = n()) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(desc(n_IPR_acc_background)) %>%
  dplyr::mutate(tran = sub("\\.[^.]*$", "", tran)) # 2,360,448 wild strain pangenome genes have an IPR annotation

ipr_background_genes <- all_ipr_background %>% dplyr::distinct(tran) %>% dplyr::pull()

# Counting the number of different gene classes in WS HDRs based on IPR classification
ws_hdr_gene_class <- long_class %>%
  dplyr::filter(gene %in% all_ws_hdr_genes) %>%
  dplyr::left_join(all_ipr_background %>% dplyr::select(-n_IPR_acc_background), by = c("gene" = "tran")) %>% # a total of 206,225 genes are in WS HDRs
  dplyr::filter(!is.na(IPR_description)) %>% # 159,191 of the 206,255 HDR genes have an IPR annotation
  dplyr::mutate(rough_gene_class = ifelse(grepl("7TM", IPR_description), "GPCR", 
                                          ifelse(grepl("C-type lectin", IPR_description), "C-type lectin", 
                                                 ifelse(grepl("Cytochrome P450", IPR_description), "Cytochrome P450", 
                                                        ifelse(grepl("F-box", IPR_description), "F-box", IPR_description))))) # assigning gene class based on IPR annotation

# Total number of each (roughly estimated) gene class
ws_hdr_gene_class_count <- ws_hdr_gene_class %>%
  dplyr::select(strain, gene, class, rough_gene_class) %>%
  dplyr::filter(rough_gene_class == "GPCR" | rough_gene_class == "C-type lectin" | rough_gene_class == "Cytochrome P450" | rough_gene_class == "F-box") %>%
  dplyr::group_by(rough_gene_class) %>%
  dplyr::mutate(class_count = n()) %>%
  dplyr::distinct(rough_gene_class, class_count) %>%
  dplyr::mutate(average_per_strain_hdrs = class_count / 140)

# Most abundant IPR terms 
ipr_desc_most <- ws_hdr_gene_class %>%
  dplyr::select(strain, gene, class, IPR_description) %>%
  dplyr::group_by(IPR_description) %>%
  dplyr::mutate(count_ipr = n()) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(IPR_description, count_ipr) %>%
  dplyr::arrange(desc(count_ipr))
  
# Most abundant IPR term for each strains hyper-divergent genome!
ipr_desc_most_strain <- ws_hdr_gene_class %>%
  dplyr::select(strain, gene, class, IPR_description) %>%
  dplyr::group_by(strain, IPR_description) %>%
  dplyr::mutate(count_ipr = n()) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(strain, IPR_description, count_ipr) %>%
  dplyr::group_by(strain) %>%
  dplyr::arrange(desc(count_ipr)) %>%
  dplyr::slice_head(n = 1) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(IPR_description, desc(count_ipr)) %>%
  dplyr::mutate(strain = factor(strain, levels = strain)) %>%
  dplyr::mutate(IPR_description = ifelse(IPR_description == "Domain of unknown function DUF38/FTH, Caenorhabditis species", "Domain of unknown function DUF38/FTH,\nCaenorhabditis species", IPR_description))

geo.colors <- c("7TM GPCR, serpentine receptor class h (Srh)" = "blue", "7TM GPCR, serpentine receptor class w (Srw)" = "cyan",
                "C-type lectin fold" = "forestgreen", "Domain of unknown function DUF38/FTH,\nCaenorhabditis species" = "gray40", 
                "E3 ubiquitin-protein ligase RING-type" = "chocolate", "F-box domain" = "Orange", "GPCR, rhodopsin-like, 7TM" = "skyblue", "MATH/TRAF domain" = "red3", 
                "Nuclear hormone receptor-like domain superfamily" = "red", "SKP1/BTB/POZ domain superfamily" = "purple", "TRA-1 regulated" = 'violet', "Zinc finger, RING/FYVE/PHD-type" = "pink")

# Most abundant IPR term
ipr_abundance_hdr <- ggplot(data = ipr_desc_most_strain) + 
  geom_col(aes(x = strain, y = count_ipr, fill = IPR_description)) +
  scale_fill_manual(values = geo.colors) +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(color = 'black', fill = NA),
    axis.text.y = element_text(size = 10, color = 'black'),
    axis.text.x = element_text(size = 3, color = 'black', angle = 75, hjust = 1),
    axis.title.y = element_text(size = 10, color = 'black'),
    axis.title.x = element_blank(),
    legend.text = element_text(size = 4, color = 'black'),
    legend.title = element_text(size = 5, color = 'black')
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, .05))) +
  labs(y = "IPR description count (number of genes)", fill = "Most abundant IPR description in HDRs")

ggsave("../../figures/supplementary/most_abundant_IPR_inHDRs.png", ipr_abundance_hdr, width = 7.5, height = 5, dpi = 600)

# Assessing if IPR term is enriched in specific geo locations
geo_initial <- readr::read_tsv("../../processed_data/genome_resources/isotypes/elegans_isotypes_sampling_geo.tsv")
hawaii_islands <- readr::read_tsv("../../processed_data/genome_resources/isotypes/elegans_isotypes_sampling_geo_hawaii_islands.tsv") %>% dplyr::select(isotype,collection_island_Hawaii)

# Isolation site of each wild strain
geo <- geo_initial %>%
  dplyr::left_join(hawaii_islands, by = "isotype") %>%
  dplyr::mutate(geo = ifelse(geo == "Hawaii",collection_island_Hawaii,geo)) %>%
  dplyr::select(isotype, lat, long, geo) %>%
  dplyr::filter(isotype %in% WSs)

ipr_most_geo <- ipr_desc_most_strain %>%
  dplyr::left_join(geo, by = c("strain" = "isotype")) %>%
  dplyr::mutate(lat = as.numeric(lat),
                long = as.numeric(long))

# world polygons
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

ggplot() +
  geom_sf(data = world, fill = "gray95", color = "gray70", linewidth = 0.2) +
  geom_point(data = ipr_most_geo,
             aes(x = long, y = lat, color = IPR_description), 
             size = 3, alpha = 0.85) +
  scale_color_manual(values = geo.colors) +
  coord_sf(expand = FALSE) +
  theme_minimal() +
  labs(x = NULL, y = NULL, color = "Most abundant IPR description")

# just the Hawaiian islands
ggplot() +
  geom_sf(data = world, fill = "gray95", color = "gray70", linewidth = 0.2) +
  geom_point(
    data = ipr_most_geo, 
    aes(x = long, y = lat, color = IPR_description),
    size = 3, position = position_jitter(width = 0.1, height = 0.1), alpha = 0.85) +
  scale_color_manual(values = geo.colors) +
  coord_sf(xlim = c(-161, -154), ylim = c(18.5, 23), expand = FALSE) +
  theme_minimal() +
  labs(x = NULL, y = NULL, color = "Most abundant IPR description")


#==============================================================================================================================================================================================================================#

# INTERPROSCAN enrichment test - entire pangenome as background

#==============================================================================================================================================================================================================================#
################################################# PRIVATE HDR genes ##############################################################

# Define universe & HDR membership (annotated-only universe) 
univ_genes <- ipr_background_genes
hdr_genes  <- intersect(ws_hdr_priv_genes, ipr_background_genes) 

N <- length(univ_genes)
n <- length(hdr_genes)

# Counts per IPR (k = in universe, x = in HDR subset) 
k_tbl <- all_ipr_background %>%
  dplyr::count(IPR_accession, name = "k")

x_tbl <- all_ipr_background %>%
  dplyr::filter(tran %in% ws_hdr_priv_genes) %>%
  dplyr::count(IPR_accession, name = "x")

desc_tbl <- all_ipr_background %>%
  dplyr::distinct(IPR_accession, IPR_description)

# Hypergeometric enrichment (one-sided)
ipr_enrichment <- k_tbl %>%
  dplyr::left_join(x_tbl, by = "IPR_accession") %>%
  dplyr::mutate(x = tidyr::replace_na(x, 0L)) %>%
  dplyr::mutate(
    pval = stats::phyper(q = x - 1, m = k, n = N - k, k = n, lower.tail = FALSE),
    expected = (n * k) / N, # if IPR genes are randomly distributed, you’d expect this many HDR genes to carry the IPR.
    enrich_ratio = dplyr::if_else(expected > 0, x / expected, NA_real_), # (x HDR with IPR / k background with IPR) / (n HDR genes / N background genes)
    # odds ratio with Haldane–Anscombe correction (adding 0.5 to each cell to avoid infinities)
    OR = {
      a <- x + 0.5                                # HDR & has IPR
      b <- (n - x) + 0.5                          # HDR & no IPR
      c <- (k - x) + 0.5                          # non-HDR & has IPR
      d <- (N - n - (k - x)) + 0.5                # non-HDR & no IPR
      (a / b) / (c / d)
    },
    FDR_p.adjust = stats::p.adjust(pval, method = "BH")
  ) %>%
  dplyr::left_join(desc_tbl, by = "IPR_accession") %>%
  dplyr::mutate(N = N, n = n) %>%
  dplyr::select(IPR_accession, IPR_description, x, k, n, N, expected, enrich_ratio, OR, pval, FDR_p.adjust) %>%
  dplyr::arrange(FDR_p.adjust, dplyr::desc(enrich_ratio))

ipr_sig <- ipr_enrichment %>%
  dplyr::filter(FDR_p.adjust < 0.05)

ipr_sig %>% dplyr::slice_head(n = 20)

ipr_sig_gene_collapsed <- all_ipr_background %>%
  dplyr::filter(IPR_accession %in% ipr_sig$IPR_accession) %>%
  dplyr::group_by(IPR_accession) %>%
  dplyr::summarise(
    IPR_description = dplyr::first(stats::na.omit(IPR_description)),
    n_genes_HDR = dplyr::n_distinct(tran[tran %in% ws_hdr_priv_genes]),
    genes_HDR   = paste(sort(unique(tran[tran %in% ws_hdr_priv_genes])), collapse = ", "),
    n_genes_all = dplyr::n_distinct(tran),
    genes_all   = paste(sort(unique(tran)), collapse = ", "),
    .groups = "drop") %>%
  dplyr::left_join(ipr_sig %>% dplyr::select(IPR_accession, x, k, n, N, expected, enrich_ratio, OR, pval, FDR_p.adjust), by = "IPR_accession") %>%
  dplyr::mutate(Region = "Private pangenome")

data_plt_priv <- ipr_sig_gene_collapsed %>% 
  dplyr::filter(IPR_accession != "IPR008164") %>% # Excluding annotation of "Repeat of unknown function XGLTT"
  dplyr::arrange(FDR_p.adjust) %>% 
  dplyr::filter(IPR_accession != "") %>% 
  dplyr::slice_head(n = 20) %>% dplyr::arrange(desc(FDR_p.adjust)) %>% 
  dplyr::mutate(plotpoint = dplyr::row_number())



############################################################################### ACCESSORY HDR genes #############################################################
univ_genes <- ipr_background_genes
hdr_genes  <- intersect(ws_hdr_acc_genes, ipr_background_genes) 

N <- length(univ_genes)
n <- length(hdr_genes)

# Counts per IPR (k = in universe, x = in HDR subset) 
k_tbl <- all_ipr_background %>%
  dplyr::count(IPR_accession, name = "k") %>%
  dplyr::mutate(k = as.numeric(k))

x_tbl <- all_ipr_background %>%
  dplyr::filter(tran %in% ws_hdr_acc_genes) %>%
  dplyr::count(IPR_accession, name = "x")

desc_tbl <- all_ipr_background %>%
  dplyr::distinct(IPR_accession, IPR_description)

# Hypergeometric enrichment (one-sided)
ipr_enrichment <- k_tbl %>%
  dplyr::left_join(x_tbl, by = "IPR_accession") %>%
  dplyr::mutate(x = tidyr::replace_na(x, 0L)) %>%
  dplyr::mutate(
    pval = stats::phyper(q = x - 1, m = k, n = N - k, k = n, lower.tail = FALSE),
    expected = (n * k) / N, # if IPR genes are randomly distributed, you’d expect this many HDR genes to carry the IPR.
    enrich_ratio = dplyr::if_else(expected > 0, x / expected, NA_real_), # (x HDR with IPR / k background with IPR) / (n HDR genes / N background genes)
    # odds ratio with Haldane–Anscombe correction (adding 0.5 to each cell to avoid infinities)
    OR = {
      a <- x + 0.5                                # HDR & has IPR
      b <- (n - x) + 0.5                          # HDR & no IPR
      c <- (k - x) + 0.5                          # non-HDR & has IPR
      d <- (N - n - (k - x)) + 0.5                # non-HDR & no IPR
      (a / b) / (c / d)
    },
    FDR_p.adjust = stats::p.adjust(pval, method = "BH")
  ) %>%
  dplyr::left_join(desc_tbl, by = "IPR_accession") %>%
  dplyr::mutate(N = N, n = n) %>%
  dplyr::select(IPR_accession, IPR_description, x, k, n, N, expected, enrich_ratio, OR, pval, FDR_p.adjust) %>%
  dplyr::arrange(FDR_p.adjust, dplyr::desc(enrich_ratio))

ipr_sig <- ipr_enrichment %>%
  dplyr::filter(FDR_p.adjust < 0.05)

ipr_sig %>% dplyr::slice_head(n = 20)

ipr_sig_gene_collapsed <- all_ipr_background %>%
  dplyr::filter(IPR_accession %in% ipr_sig$IPR_accession) %>%
  dplyr::group_by(IPR_accession) %>%
  dplyr::summarise(
    IPR_description = dplyr::first(stats::na.omit(IPR_description)),
    n_genes_HDR = dplyr::n_distinct(tran[tran %in% ws_hdr_acc_genes]),
    genes_HDR   = paste(sort(unique(tran[tran %in% ws_hdr_acc_genes])), collapse = ", "),
    n_genes_all = dplyr::n_distinct(tran),
    genes_all   = paste(sort(unique(tran)), collapse = ", "),
    .groups = "drop") %>%
  dplyr::left_join(ipr_sig %>% dplyr::select(IPR_accession, x, k, n, N, expected, enrich_ratio, OR, pval, FDR_p.adjust), by = "IPR_accession") %>%
  dplyr::mutate(Region = "Accessory pangenome")


data_plt_acc <- ipr_sig_gene_collapsed %>% 
  dplyr::arrange(FDR_p.adjust) %>% 
  dplyr::slice_head(n = 20) %>% 
  dplyr::arrange(desc(FDR_p.adjust)) %>% 
  dplyr::mutate(plotpoint = dplyr::row_number())


############################################################################### CORE HDR genes #############################################################
univ_genes <- ipr_background_genes
hdr_genes  <- intersect(ws_hdr_core_genes, ipr_background_genes) 

N <- length(univ_genes)
n <- length(hdr_genes)

# Counts per IPR (k = in universe, x = in HDR subset) 
k_tbl <- all_ipr_background %>%
  dplyr::count(IPR_accession, name = "k") %>%
  dplyr::mutate(k = as.numeric(k))

x_tbl <- all_ipr_background %>%
  dplyr::filter(tran %in% ws_hdr_core_genes) %>%
  dplyr::count(IPR_accession, name = "x")

desc_tbl <- all_ipr_background %>%
  dplyr::distinct(IPR_accession, IPR_description)

# Hypergeometric enrichment (one-sided)
ipr_enrichment <- k_tbl %>%
  dplyr::left_join(x_tbl, by = "IPR_accession") %>%
  dplyr::mutate(x = tidyr::replace_na(x, 0L)) %>%
  dplyr::mutate(
    pval = stats::phyper(q = x - 1, m = k, n = N - k, k = n, lower.tail = FALSE),
    expected = (n * k) / N, # if IPR genes are randomly distributed, you’d expect this many HDR genes to carry the IPR.
    enrich_ratio = dplyr::if_else(expected > 0, x / expected, NA_real_), # (x HDR with IPR / k background with IPR) / (n HDR genes / N background genes)
    # odds ratio with Haldane–Anscombe correction (adding 0.5 to each cell to avoid infinities)
    OR = {
      a <- x + 0.5                                # HDR & has IPR
      b <- (n - x) + 0.5                          # HDR & no IPR
      c <- (k - x) + 0.5                          # non-HDR & has IPR
      d <- (N - n - (k - x)) + 0.5                # non-HDR & no IPR
      (a / b) / (c / d)
    },
    FDR_p.adjust = stats::p.adjust(pval, method = "BH")
  ) %>%
  dplyr::left_join(desc_tbl, by = "IPR_accession") %>%
  dplyr::mutate(N = N, n = n) %>%
  dplyr::select(IPR_accession, IPR_description, x, k, n, N, expected, enrich_ratio, OR, pval, FDR_p.adjust) %>%
  dplyr::arrange(FDR_p.adjust, dplyr::desc(enrich_ratio))

ipr_sig <- ipr_enrichment %>%
  dplyr::filter(FDR_p.adjust < 0.05)

ipr_sig %>% dplyr::slice_head(n = 20)

ipr_sig_gene_collapsed <- all_ipr_background %>%
  dplyr::filter(IPR_accession %in% ipr_sig$IPR_accession) %>%
  dplyr::group_by(IPR_accession) %>%
  dplyr::summarise(
    IPR_description = dplyr::first(stats::na.omit(IPR_description)),
    n_genes_HDR = dplyr::n_distinct(tran[tran %in% ws_hdr_core_genes]),
    genes_HDR   = paste(sort(unique(tran[tran %in% ws_hdr_core_genes])), collapse = ", "),
    n_genes_all = dplyr::n_distinct(tran),
    genes_all   = paste(sort(unique(tran)), collapse = ", "),
    .groups = "drop") %>%
  dplyr::left_join(ipr_sig %>% dplyr::select(IPR_accession, x, k, n, N, expected, enrich_ratio, OR, pval, FDR_p.adjust), by = "IPR_accession") %>%
  dplyr::mutate(Region = "Core pangenome")


data_plt_core <- ipr_sig_gene_collapsed %>% 
  dplyr::arrange(FDR_p.adjust) %>% 
  dplyr::slice_head(n = 20) %>% 
  dplyr::arrange(desc(FDR_p.adjust)) %>% 
  dplyr::mutate(plotpoint = dplyr::row_number())


# Concatenating enrichment of all three gene sets and pulling the most enriched IPR terms among all
concat_IPR_enrich <- data_plt_priv %>% dplyr::bind_rows(data_plt_acc, data_plt_core) %>%
  dplyr::select(-plotpoint) %>%
  dplyr::arrange(FDR_p.adjust) %>%
  dplyr::slice_head(n = 40) %>%
  dplyr::arrange(desc(FDR_p.adjust)) %>%
  dplyr::rename(`Gene set` = Region) %>%
  dplyr::mutate(`Gene set` = ifelse(`Gene set` == "Accessory pangenome", "Accessory", 
                                    ifelse(`Gene set` == "Core pangenome", "Core", "Private"))) %>%
  dplyr::mutate(FDR_p.adjust = pmax(FDR_p.adjust, 1e-320))  # p values that are "00000e+00" are clipped so that they aren't displayed as infinity


diff <- concat_IPR_enrich %>%
  dplyr::arrange(n_genes_HDR) %>%
  dplyr::mutate(plotpoint = dplyr::row_number())

plot_IPR_all_diff <- ggplot(diff) +
  geom_vline(xintercept = -log10(0.05), color='blue', linewidth=0.4) +
  geom_point(aes(x = n_genes_HDR, y = plotpoint, size = enrich_ratio, fill = -log10(FDR_p.adjust), shape = `Gene set`)) +
  scale_y_continuous(breaks = diff$plotpoint, labels = diff$IPR_description, name = "", expand = c(0.02,0.02)) +
  scale_shape_manual(values = c("Core" = 21, "Accessory" = 22, "Private" = 24)) +
  scale_fill_gradient(low = "yellow", high = "red", breaks = c(round(min(-log10(diff$FDR_p.adjust))), round((max(-log10(diff$FDR_p.adjust)) + min(-log10(diff$FDR_p.adjust))) / 2), round(max(-log10(diff$FDR_p.adjust))))) +
  scale_size_continuous(range = c(1, 5), name = "Fold enrichment", breaks = pretty(diff$enrich_ratio, n = 4)) +
  coord_cartesian(xlim = c(0, 12000)) +
  theme(axis.text.x = element_text(size=10, color='black'),
        axis.text.y = element_text(size=10, color='black'),
        axis.title = element_text(size=11, color='black'),
        plot.title = element_blank(),
        legend.title = element_text(size = 10, color='black', hjust = 1),
        legend.text = element_text(size = 10, color='black', hjust = 1),
        legend.position = "inside",
        legend.position.inside = c(0.8, 0.2),
        legend.direction = "horizontal", legend.box = "vertical",
        legend.spacing.y = unit(0.0001, 'cm'),
        legend.key.height = unit(0.01, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.box.just = "right",
        text = element_text(family="Helvetica"),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        plot.margin = margin(b = 5, t = 10, r = 25, l = 2, unit = "pt")) +
  guides(
    fill = guide_colourbar(nrow=1, order = 1, title.position = "top", force = TRUE, barwidth = 5, barheight = 1),
    size = guide_legend(nrow=1, order = 2, title.position = "top", title.hjust = 1, force = TRUE),
    shape = guide_legend(nrow = 3, order = 3, title.position = "top", title.hjust = 1, override.aes = list(size = 0.5))) +
  labs(title = "Enriched IPR terms for wild strain genes in HDRs",  x = "Gene count", size = "Fold enrichment", fill = expression(-log[10]~"(corrected p-value)"))
plot_IPR_all_diff

write.table(diff, "../../tables/enriched_IPR_inHDRs.tsv", row.names = F, col.names = T, sep = "\t")





##########################################################################################
### Look at rarefaction of enriched gene families ###
##########################################################################################
pan_ipr_cleaned <-readr::read_tsv("../../processed_data/genome_resources/annotation/IPR_annotation_142strains.tsv", col_names = c("tran", "MD5_digest", "seq_length", "app", "signature_accession", "signature_description", "start", "end", "score", "status", "date", "IPR_accession","IPR_description","GO", "pathways")) %>% 
  tidyr::separate(tran, into = c("strain","gene"), sep = "_") %>% 
  dplyr::filter(IPR_description != "-") %>% 
  dplyr::select(strain,gene,IPR_description)

pan_gpcrs <- pan_ipr_cleaned %>% 
  dplyr::filter(grepl("7TM", IPR_description)) %>%
  dplyr::distinct(strain,gene) %>%
  dplyr::group_by(strain) %>%
  dplyr::mutate(n_gpcrs = n()) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(desc(n_gpcrs)) 


cleaned_pan_gpcrs <- pan_gpcrs %>% dplyr::mutate(gene = sub("\\.[^.]*$", "", gene)) %>% dplyr::select(-n_gpcrs) %>% dplyr::filter(strain != "N2", strain != "CGC1")
# write.table(cleaned_pan_gpcrs,"../../processed_data/processed_data/genome_resources/annotation/140_wild_strains_IPR_gpcrs.tsv", quote = F, row.names = F, col.names = T, sep = '\t')

all_relations <- readr::read_tsv("../../processed_data/orthology/OG_relations_matrix_count.tsv")

private_freq <- 1/142

priv <- all_relations %>%
  dplyr::mutate(across(2:(ncol(.)), ~ ifelse(. >= 1, 1, .))) %>%
  dplyr::mutate(sum = rowSums(across(-1, ~ ., .names = NULL), na.rm = TRUE)) %>%
  dplyr::mutate(freq = (sum / 142)) %>%
  dplyr::mutate(
    class = case_when(
      freq == 1 ~ "core",
      freq > private_freq & freq < 1 ~ "accessory",
      freq == private_freq ~ "private",
      TRUE ~ "undefined"
    )
  ) %>%
  dplyr::select(Orthogroup, class) %>%
  dplyr::filter(class == "private") %>%
  dplyr::pull(Orthogroup)


all_OGs <- readr::read_tsv("../../tables/all_OGs_matrix.tsv")

priv_expanded <- all_OGs %>% dplyr::filter(Orthogroup %in% priv) %>% tidyr::separate_rows(everything(), sep = ', ') 

priv_wide <- priv_expanded %>% tidyr::pivot_longer(cols = -Orthogroup, names_to = "strain", values_to = "gene") %>% dplyr::filter(!is.na(gene)) %>% dplyr::select(-Orthogroup)

merged <- priv_wide %>% dplyr::left_join(pan_gpcrs, by = "strain") %>% dplyr::mutate(priv_gpcrs = ifelse(gene.x == gene.y, T, F)) 

priv_gpcrs <- merged %>% dplyr::filter(priv_gpcrs == T) %>% dplyr::group_by(strain) %>% dplyr::mutate(number_priv_gpcrs = n()) %>% dplyr::ungroup() %>% dplyr::arrange(desc(number_priv_gpcrs))

rarefact <- priv_gpcrs %>% dplyr::select(strain,number_priv_gpcrs) %>% dplyr::distinct() %>%
  dplyr::mutate(iterative_sum = cumsum(number_priv_gpcrs)) %>%
  dplyr::mutate(number_genomes = row_number()) %>%
  dplyr::mutate(number_genomes = factor(number_genomes, levels = number_genomes))


#################### Looking at F-box genes ##############################################################################################################################
pan_fbox <- pan_ipr_cleaned %>% 
  dplyr::filter(grepl("F-box", IPR_description)) %>%
  dplyr::distinct(strain,gene) %>%
  dplyr::group_by(strain) %>%
  dplyr::mutate(n_fbox = n()) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(desc(n_fbox)) 


cleaned_pan_fbox <- pan_fbox %>% dplyr::mutate(gene = sub("\\.[^.]*$", "", gene)) %>% dplyr::select(-n_fbox) %>% dplyr::filter(strain != "N2", strain != "CGC1")
# write.table(cleaned_pan_fbox,"../../processed_data/processed_data/genome_resources/annotation/140_wild_strains_IPR_fBox.tsv", quote = F, row.names = F, col.names = T, sep = '\t')


# rarefaction
merged_fbox <- priv_wide %>% dplyr::left_join(pan_fbox, by = "strain") %>% dplyr::mutate(priv_fbox = ifelse(gene.x == gene.y, T, F)) 

private_fbox <- merged_fbox %>% dplyr::filter(priv_fbox == T) %>% dplyr::group_by(strain) %>% dplyr::mutate(number_priv_fbox = n()) %>% dplyr::ungroup() %>% dplyr::arrange(desc(number_priv_fbox))

rarefact_fbox <- private_fbox %>% dplyr::select(strain,number_priv_fbox) %>% dplyr::distinct() %>%
  dplyr::mutate(iterative_sum = cumsum(number_priv_fbox)) %>%
  dplyr::mutate(number_genomes = row_number()) %>%
  dplyr::mutate(number_genomes = factor(number_genomes, levels = number_genomes))



#################### Looking at C-type lectin genes ##############################################################################################################################
pan_lectin <- pan_ipr_cleaned %>% 
  dplyr::filter(grepl("C-type lectin", IPR_description)) %>%
  dplyr::distinct(strain,gene) %>%
  dplyr::group_by(strain) %>%
  dplyr::mutate(n_lectin = n()) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(desc(n_lectin)) 


cleaned_pan_lectin <- pan_lectin %>% dplyr::mutate(gene = sub("\\.[^.]*$", "", gene)) %>% dplyr::select(-n_lectin) %>% dplyr::filter(strain != "N2", strain != "CGC1")
# write.table(cleaned_pan_lectin,"../../processed_data/processed_data/genome_resources/annotation/140_wild_strains_IPR_CtypeLectins.tsv", quote = F, row.names = F, col.names = T, sep = '\t')


# rarefaction
merged_lectin <- priv_wide %>% dplyr::left_join(pan_lectin, by = "strain") %>% dplyr::mutate(priv_lectin = ifelse(gene.x == gene.y, T, F)) 

private_lectin <- merged_lectin %>% dplyr::filter(priv_lectin == T) %>% dplyr::group_by(strain) %>% dplyr::mutate(number_priv_lectin = n()) %>% dplyr::ungroup() %>% dplyr::arrange(desc(number_priv_lectin))

rarefact_lectin <- private_lectin %>% dplyr::select(strain,number_priv_lectin) %>% dplyr::distinct() %>%
  dplyr::mutate(iterative_sum = cumsum(number_priv_lectin)) %>%
  dplyr::mutate(number_genomes = row_number()) %>%
  dplyr::mutate(number_genomes = factor(number_genomes, levels = number_genomes))





#################### Looking at Cytochrome P450 genes ##############################################################################################################################
pan_cyto <- pan_ipr_cleaned %>% 
  dplyr::filter(grepl("Cytochrome P450", IPR_description)) %>% 
  dplyr::distinct(strain,gene) %>%
  dplyr::group_by(strain) %>%
  dplyr::mutate(n_cyto = n()) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(desc(n_cyto)) 


cleaned_pan_cyto <- pan_cyto %>% dplyr::mutate(gene = sub("\\.[^.]*$", "", gene)) %>% dplyr::select(-n_cyto) %>% dplyr::filter(strain != "N2", strain != "CGC1")
# write.table(cleaned_pan_cyto,"../../processed_data/processed_data/genome_resources/annotation/140_wild_strains_IPR_cytochromeP450.tsv", quote = F, row.names = F, col.names = T, sep = '\t')


# rarefaction
merged_cyto <- priv_wide %>% dplyr::left_join(pan_cyto, by = "strain") %>% dplyr::mutate(priv_cyto = ifelse(gene.x == gene.y, T, F)) 

private_cyto <- merged_cyto %>% dplyr::filter(priv_cyto == T) %>% dplyr::group_by(strain) %>% dplyr::mutate(number_priv_cyto = n()) %>% dplyr::ungroup() %>% dplyr::arrange(desc(number_priv_cyto))

rarefact_cyto <- private_cyto %>% dplyr::select(strain, number_priv_cyto) %>% dplyr::distinct() %>%
  dplyr::mutate(iterative_sum = cumsum(number_priv_cyto)) %>%
  dplyr::mutate(number_genomes = row_number()) %>%
  dplyr::mutate(number_genomes = factor(number_genomes, levels = number_genomes))



#################### Looking at nuclear hormone receptors ##############################################################################################################################
pan_nhr <- pan_ipr_cleaned %>% 
  dplyr::filter(grepl("uclear hormone receptor", IPR_description)) %>%
  dplyr::distinct(strain,gene) %>%
  dplyr::group_by(strain) %>%
  dplyr::mutate(n_nhr = n()) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(desc(n_nhr)) 


cleaned_pan_nhr <- pan_nhr %>% dplyr::mutate(gene = sub("\\.[^.]*$", "", gene)) %>% dplyr::select(-n_nhr) %>% dplyr::filter(strain != "N2", strain != "CGC1")
# write.table(cleaned_pan_lectin,"../../processed_data/processed_data/genome_resources/annotation/140_wild_strains_IPR_nhr.tsv", quote = F, row.names = F, col.names = T, sep = '\t')


# rarefaction
merged_nhr <- priv_wide %>% dplyr::left_join(pan_nhr, by = "strain") %>% dplyr::mutate(priv_nhr = ifelse(gene.x == gene.y, T, F)) 

private_nhr <- merged_nhr %>% dplyr::filter(priv_nhr == T) %>% dplyr::group_by(strain) %>% dplyr::mutate(number_priv_nhr = n()) %>% dplyr::ungroup() %>% dplyr::arrange(desc(number_priv_nhr))

rarefact_nhr <- private_nhr %>% dplyr::select(strain,number_priv_nhr) %>% dplyr::distinct() %>%
  dplyr::mutate(iterative_sum = cumsum(number_priv_nhr)) %>%
  dplyr::mutate(number_genomes = row_number()) %>%
  dplyr::mutate(number_genomes = factor(number_genomes, levels = number_genomes))



# Prepping data for rarefaction plotting
rarefact_extended_gp <- rarefact %>% dplyr::bind_rows(tibble(strain = paste0("genome_", (max(as.numeric(.$number_genomes)) + 1):142),
                                                             number_priv_gpcrs = 0, iterative_sum = max(.$iterative_sum), number_genomes = factor((max(as.numeric(.$number_genomes)) + 1):142)))
rarefact_extended_fbox <- rarefact_fbox %>% dplyr::bind_rows(tibble(strain = paste0("genome_", (max(as.numeric(.$number_genomes)) + 1):142),
                                                                    number_priv_gpcrs = 0, iterative_sum = max(.$iterative_sum), number_genomes = factor((max(as.numeric(.$number_genomes)) + 1):142)))
rarefact_extended_lectin <- rarefact_lectin %>% dplyr::bind_rows(tibble(strain = paste0("genome_", (max(as.numeric(.$number_genomes)) + 1):142),
                                                                        number_priv_gpcrs = 0, iterative_sum = max(.$iterative_sum), number_genomes = factor((max(as.numeric(.$number_genomes)) + 1):142)))
rarefact_extended_cyto <- rarefact_cyto %>% dplyr::bind_rows(tibble(strain = paste0("genome_", (max(as.numeric(.$number_genomes)) + 1):142),
                                                                    number_priv_gpcrs = 0, iterative_sum = max(.$iterative_sum), number_genomes = factor((max(as.numeric(.$number_genomes)) + 1):142)))
rarefact_extended_nhr <- rarefact_nhr %>% dplyr::bind_rows(tibble(strain = paste0("genome_", (max(as.numeric(.$number_genomes)) + 1):142),
                                                                  number_priv_gpcrs = 0, iterative_sum = max(.$iterative_sum), number_genomes = factor((max(as.numeric(.$number_genomes)) + 1):142)))


# Rarefaction of all three gene families
rarefaction_gene_families <- ggplot() +
  geom_point(data = rarefact_extended_gp, aes(x = number_genomes, y = iterative_sum, color = "GPCRs"), size = 1) +
  geom_point(data = rarefact_extended_fbox, aes(x = number_genomes, y = iterative_sum, color = "F-box"), size = 1) +
  geom_point(data = rarefact_extended_lectin, aes(x = number_genomes, y = iterative_sum, color = "C-type lectins"), size = 1) +
  geom_point(data = rarefact_extended_cyto, aes(x = number_genomes, y = iterative_sum, color = "Cytochrome P450s"), size = 1) + # these are enriched in core
  geom_point(data = rarefact_extended_nhr, aes(x = number_genomes, y = iterative_sum, color = "Nuclear hormone receptors"), size = 1) + # these are enriched in core 
  scale_color_manual(values = c("GPCRs" = "olivedrab", "F-box" = "firebrick", 
                                "C-type lectins" = "steelblue", "Cytochrome P450s" = "orange", "Nuclear hormone receptors" = "violet")) +
  theme(
    panel.background = element_blank(),
    legend.title = element_blank(),
    legend.box.background = element_rect(color = 'black', fill = NA),
    legend.position = "inside",
    legend.position.inside = c(0.8,0.4),
    legend.text = element_text(size = 11, color = 'black'),
    panel.border = element_rect(fill = NA, color = "black"),
    plot.margin = margin(l = 5, r = 10, t = 5, b = 5),
    axis.title = element_text(size = 11),
    axis.text = element_text(size =10, color = 'black'),
    axis.text.x = element_text(size = 10, color = 'black')) +
  scale_x_discrete(breaks = c(seq(0, 125, 25), 142)) +
  labs(y = "Number of private genes", x = "Number of genomes")
rarefaction_gene_families

# Save rarefaction plot
# ggsave("../../figures/supplementary/gene_family_rarefaction.png",rarefaction_gene_families, width = 7.5, height = 6, dpi = 600)

