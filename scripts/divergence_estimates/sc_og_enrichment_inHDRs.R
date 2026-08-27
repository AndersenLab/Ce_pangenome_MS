library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(enrichplot)
library(cowplot)
library(stringr)

# ======================================================================================================================================================================================== #
# Loading in all genes in the pangenome with resolution on the gene sets they contribute to
# ======================================================================================================================================================================================== #
long_class <- readr::read_tsv("../../processed_data/orthology/all_genes_class_OGs.tsv") %>%
  dplyr::mutate(gene = gsub('transcript_','',gene)) %>%
  dplyr::mutate(gene = paste0(strain,"_",gene))

# ======================================================================================================================================================================================== #
# IPR results
# ======================================================================================================================================================================================== #
all_ipr <- readr::read_tsv("../../tables/IPR_annotation_142strains.tsv", col_names = c("tran", "MD5_digest", "seq_length", "app", "signature_accession", "signature_description", "start", "end", "score", "status", "date", "IPR_accession","IPR_description","GO", "pathways")) %>%
  dplyr::select(tran, IPR_accession, IPR_description, GO)

all_sc_ogs <- readr::read_tsv("../../processed_data/divergence_estimates/filtered_sc_OGs.tsv") %>% dplyr::pull(Orthogroup)
all_sc_genes <- long_class %>% dplyr::filter(Orthogroup %in% all_sc_ogs) %>% dplyr::pull(gene)

hdr_sc_ogs <- readr::read_tsv("../../processed_data/divergence_estimates/filtered_sc_OGs_inHDRs.tsv") %>% dplyr::pull(Orthogroup)
hdr_sc_genes <- long_class %>% dplyr::filter(Orthogroup %in% hdr_sc_ogs) %>% dplyr::pull(gene)

all_ipr_background <- all_ipr %>% 
  dplyr::filter(!is.na(IPR_description) & IPR_accession != "-") %>%
  dplyr::distinct(tran, IPR_accession, IPR_description) %>% # only one type of IPR annotation per gene
  dplyr::group_by(IPR_accession) %>%
  dplyr::mutate(n_IPR_acc_background = n()) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(desc(n_IPR_acc_background)) %>%
  dplyr::mutate(tran = sub("\\.[^.]*$", "", tran)) # 2,360,448 wild strain pangenome genes have an IPR annotation

ipr_background_genes <- all_ipr_background %>% dplyr::distinct(tran) %>% dplyr::filter(tran %in% all_sc_genes) %>% dplyr::pull() # only pulling single-copy orthologs for the background set

#==============================================================================================================================================================================================================================#

# INTERPROSCAN enrichment test - only single-copy orthologs as the background

#==============================================================================================================================================================================================================================#
# Define universe & HDR membership (annotated-only universe) 
univ_genes <- ipr_background_genes
hdr_genes  <- intersect(hdr_sc_genes, ipr_background_genes) 

N <- length(univ_genes)
n <- length(hdr_genes)

# Counts per IPR (k = in universe, x = in HDR subset) 
k_tbl <- all_ipr_background %>%
  dplyr::count(IPR_accession, name = "k")

x_tbl <- all_ipr_background %>%
  dplyr::filter(tran %in% hdr_sc_genes) %>%
  dplyr::count(IPR_accession, name = "x")

desc_tbl <- all_ipr_background %>%
  dplyr::distinct(IPR_accession, IPR_description)

# Hypergeometric enrichment (one-sided)
ipr_enrichment <- k_tbl %>%
  dplyr::left_join(x_tbl, by = "IPR_accession") %>%
  dplyr::mutate(x = tidyr::replace_na(x, 0L)) %>%
  dplyr::mutate(
    pval = stats::phyper(q = x - 1, m = k, n = N - k, k = n, lower.tail = FALSE),
    expected = (n * 1.0 * k) / N, # if IPR genes are randomly distributed, you’d expect this many HDR genes to carry the IPR.
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

# ipr_sig %>% dplyr::slice_head(n = 30)

ipr_sig_gene_collapsed <- all_ipr_background %>%
  dplyr::filter(IPR_accession %in% ipr_sig$IPR_accession) %>%
  dplyr::group_by(IPR_accession) %>%
  dplyr::summarise(
    IPR_description = dplyr::first(stats::na.omit(IPR_description)),
    n_genes_HDR = dplyr::n_distinct(tran[tran %in% hdr_sc_genes]),
    genes_HDR   = paste(sort(unique(tran[tran %in% hdr_sc_genes])), collapse = ", "),
    n_genes_all = dplyr::n_distinct(tran),
    genes_all   = paste(sort(unique(tran)), collapse = ", "),
    .groups = "drop") %>%
  dplyr::left_join(ipr_sig %>% dplyr::select(IPR_accession, x, k, n, N, expected, enrich_ratio, OR, pval, FDR_p.adjust), by = "IPR_accession") %>%
  dplyr::mutate(Region = "HDR single-copy orthologs") %>%
  dplyr::filter(x != k, n_genes_HDR > 284) # remove single-copy orthogroups that are ONLY found in HDRs - likely represent one or two orthogroups

data_plt <- ipr_sig_gene_collapsed %>% 
  dplyr::arrange(FDR_p.adjust) %>% 
  dplyr::filter(IPR_accession != "") %>% 
  dplyr::slice_head(n = 30) %>% dplyr::arrange(n_genes_HDR) %>%
  dplyr::mutate(plotpoint = dplyr::row_number()) %>%
  dplyr::mutate(FDR_p.adjust = pmax(FDR_p.adjust, 1e-320)) %>% # p values that are "00000e+00" are clipped so that they aren't displayed as infinity
  dplyr::mutate(orthologRatioInHDRs = x / k)


# Concatenating enrichment of all three gene sets and pulling the most enriched IPR terms among all
# concat_IPR_enrich <- data_plt_priv %>% dplyr::bind_rows(data_plt_acc, data_plt_core) %>%
#   dplyr::select(-plotpoint) %>%
#   dplyr::arrange(FDR_p.adjust) %>%
#   dplyr::slice_head(n = 40) %>%
#   dplyr::arrange(desc(FDR_p.adjust)) %>%
#   dplyr::rename(`Gene set` = Region) %>%
#   dplyr::mutate(`Gene set` = ifelse(`Gene set` == "Accessory pangenome", "Accessory", 
#                                     ifelse(`Gene set` == "Core pangenome", "Core", "Private"))) %>%
#   dplyr::mutate(FDR_p.adjust = pmax(FDR_p.adjust, 1e-320))  # p values that are "00000e+00" are clipped so that they aren't displayed as infinity
# 
# 
# diff <- concat_IPR_enrich %>%
#   dplyr::arrange(n_genes_HDR) %>%
#   dplyr::mutate(plotpoint = dplyr::row_number())

plot_IPR <- ggplot(data_plt) +
  geom_vline(xintercept = -log10(0.05), color='blue', linewidth=0.4) +
  geom_point(aes(x = n_genes_HDR, y = plotpoint, size = enrich_ratio, fill = -log10(FDR_p.adjust), shape = Region)) +
  scale_y_continuous(breaks = data_plt$plotpoint, labels = data_plt$IPR_description, name = "", expand = c(0.02,0.02)) +
  scale_shape_manual(values = c("HDR single-copy orthologs" = 21)) +
  scale_fill_gradient(low = "yellow", high = "red", breaks = c(round(min(-log10(data_plt$FDR_p.adjust))), round((max(-log10(data_plt$FDR_p.adjust)) + min(-log10(data_plt$FDR_p.adjust))) / 2), round(max(-log10(data_plt$FDR_p.adjust))))) +
  scale_size_continuous(range = c(1, 5), name = "Fold enrichment", breaks = pretty(data_plt$enrich_ratio, n = 4)) +
  coord_cartesian(xlim = c(0, 1200)) +
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
plot_IPR


