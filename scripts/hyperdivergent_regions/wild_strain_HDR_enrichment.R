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
all_ipr <- readr::read_tsv("../../processed_data/genome_resources/annotation/IPR_annotation_142strains.tsv", col_names = c("tran", "MD5_digest", "seq_length", "app", "signature_accession", "signature_description", "start", "end", "score", "status", "date", "IPR_accession","IPR_description","GO", "pathways")) %>%
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

