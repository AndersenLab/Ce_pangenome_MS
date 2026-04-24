library(dplyr)
library(ggplot2)
library(readr)
library(ggplot2)
library(tidyr)
library(dplyr)
library(ape)
library(circlize)
library(data.table)
library(grid)      

####################################################################################################
####################################################################################################

# N2 GENE EXPRESSION THAT OVERLAP WITH SVs

####################################################################################################
####################################################################################################
merged_SV <- readr::read_tsv("../../processed_data/structural_variants/Jasmine_merged_SVs.tsv")
N2_gff <- ape::read.gff("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/raw_data/assemblies/elegans/gff/longest_isoform/c_elegans.PRJNA13758.WS283.csq.PCfeaturesOnly.longest.gff3") 
n2_genes_plt <- N2_gff %>%
  dplyr::filter(type == "gene") %>%
  dplyr::mutate(attributes = gsub("ID=gene:","",attributes)) %>%
  dplyr::mutate(attributes = sub(";.*", "", attributes)) %>%
  dplyr::select(seqid,start,end,attributes) %>%
  dplyr::rename(chrom = seqid) %>% 
  dplyr::select(chrom, start, end) %>% 
  dplyr::filter(chrom != "MtDNA")

MAF_thresh <- round(0.05 * 141)

maf_filt <- merged_SV %>% 
  dplyr::select(chrom,pos,sv_length,sv_type,number_svs_merged) %>%
  dplyr::filter(number_svs_merged > MAF_thresh) %>%
  dplyr::mutate(sv_length = abs(sv_length)) %>%
  dplyr::mutate(end = pos + sv_length) %>% 
  dplyr::rename(start = pos) %>%
  dplyr::select(chrom, start, end, sv_type) %>%
  dplyr::mutate(overlap = F)

svs_dt <- as.data.table(maf_filt)
n2_genes_dt <- as.data.table(n2_genes_plt)

setkey(svs_dt, chrom, start, end)
setkey(n2_genes_dt, chrom, start, end)

svs_inCodingRegions <- data.table::foverlaps(x = svs_dt, y = n2_genes_dt, type = "any") %>% dplyr::filter(!is.na(start)) %>% dplyr::mutate(overlap = T)


# Ensuring that the foverlaps command worked correctly
test <- svs_inCodingRegions %>% dplyr::filter(start > 1600000 & end < 1700000) %>% dplyr::select(chrom,start,end) %>% dplyr::distinct(chrom,start,end)

check <- ggplot(svs_inCodingRegions %>% dplyr::filter(start > 1600000 & end < 1700000)) +
  geom_rect(data = test, aes(xmin = start / 1e6, xmax = end /1e6, ymin = 0, ymax = 0.99, fill = "N2_genes")) +
  geom_rect(aes(xmin = i.start / 1e6, xmax = i.end /1e6, ymin = 1.01, ymax = 2, fill = "SVs")) +
  scale_fill_manual(values = c("N2_genes" = "forestgreen", "SVs" = "purple")) + 
  facet_wrap(~chrom, nrow = 2, scales = "free_x") +
  theme(
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text.x = element_text(size = 12, color = 'black'),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    strip.text = element_text(size = 16, color = "black")) 
check


overlap <- svs_inCodingRegions %>% dplyr::select(chrom, i.start,i.end, sv_type, overlap) %>% dplyr::rename(start = i.start, end = i.end) %>% dplyr::distinct(chrom,start,end,sv_type, .keep_all = T)

final_stats <- maf_filt %>% 
  dplyr::left_join(overlap, by = c("chrom", "start", "end", "sv_type")) %>% 
  dplyr::mutate(region = ifelse(is.na(overlap.y),'non-PC_region','overlaps_PCgene')) %>%
  dplyr::select(chrom,start,end,sv_type,region) %>%
  dplyr::group_by(sv_type) %>%
  dplyr::mutate(total_sv_type = n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(sv_type,region) %>%
  dplyr::mutate(region_count = n()) %>%
  dplyr::ungroup()

plt_stats <- final_stats %>% dplyr::select(sv_type,region,total_sv_type,region_count) %>%
  dplyr::mutate(proportion = (region_count / total_sv_type) * 100) %>%
  dplyr::distinct() %>%
  dplyr::mutate(sv_type = factor(sv_type, levels = c("INS","DEL","INV")))

# PLACEHOLDER FOR RIGHT NOW
N2_EXP <- ggplot() +
  geom_bar(data = plt_stats, aes(x = sv_type, y = proportion, fill = region), stat = "identity") +
  geom_text(data = plt_stats, aes(x = sv_type, y = proportion, label = region_count, group = region), position = position_stack(vjust = 0.5),color = "white", size = 4) +
  scale_fill_manual(values = c("overlaps_PCgene" = "red", "non-PC_region" = "gray40")) +
  labs(y = "Proportion (%)", fill = "Region") +
  theme(panel.border = element_rect(color = 'black', fill = NA),
        panel.background = element_blank(),
        panel.grid.major= element_line(color = 'gray80'),
        panel.grid.major.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 10, color = 'black'),
        axis.text.x = element_text(color = 'black', size = 10),
        axis.title.y = element_text(size = 10, color = 'black')) +
  guides(color = "none") + # to get rid of legend for the horizontal lines
  scale_y_continuous(expand = c(0,0))
N2_EXP


####################################################################################################
####################################################################################################

# SV cumulative lengths

####################################################################################################
####################################################################################################
allcalls <- readr::read_tsv("../../processed_data/structural_variants/141_over50_PASS_variants.tsv", col_names = c("chrom", "pos", "ref", "alt", "filter", "sv_type","sv_length","strain")) %>% dplyr::select(-filter)

filt_calls <- allcalls %>% 
  dplyr::mutate(sv_length = abs(sv_length)) %>%
  dplyr::group_by(sv_type, strain) %>%
  dplyr::mutate(type_count = n()) %>%
  dplyr::ungroup()

levels3 <- filt_calls %>%
  dplyr::group_by(strain) %>%
  dplyr::mutate(total_SV_length = sum(sv_length)) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(strain,total_SV_length) %>%
  dplyr::arrange(total_SV_length) %>%
  dplyr::distinct(strain) %>%
  dplyr::pull()

total_size2 <- filt_calls %>%
  dplyr::select(sv_type,sv_length,strain) %>%
  dplyr::group_by(strain,sv_type) %>%
  dplyr::mutate(type_total_length = sum(sv_length)) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(strain,sv_type, type_total_length) %>%
  dplyr::mutate(strain = factor(strain, levels = levels3)) %>%
  dplyr::mutate(sv_type = factor(sv_type, levels = c("INV","DEL","INS"))) 

SV_LEN <- ggplot() +
  geom_bar(data = total_size2, aes(x = type_total_length / 1e6, y = strain, fill = sv_type), stat = "identity") +
  scale_fill_manual(values = c("INS" = "blue", "DEL" = "red", "INV" = "gold")) +
  labs(x = "Length of SVs (Mb)", fill = "SV type", y = "Strains") +
  theme(panel.border = element_rect(color = 'black', fill = NA),
        panel.background = element_blank(),
        panel.grid.major= element_line(color = 'gray80'),
        panel.grid.major.y = element_blank(),
        axis.text = element_text(size = 10, color = 'black'),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        # plot.margin = margin(l = 20, r = 5, t = 5, b = 5),
        legend.title = element_text(size = 10, color = 'black'),
        legend.box.background = element_rect(color = 'black', fill = NA),
        legend.position = "inside",
        legend.position.inside = c(0.8,0.2),
        legend.text = element_text(size = 10, color = 'black'),
        axis.title.x = element_text(size = 10, color = 'black')) +
  coord_cartesian(xlim = c(0,13)) +
  scale_x_continuous(breaks = seq(0,14, by =2), expand = c(0,0))
SV_LEN

# Write table of cummulative SV legnth of each type for each strain
# write.table(total_size2,"../../processed_data/structural_variants/strain_SVs_cum_length.tsv", sep = '\t', row.names = F, col.names = T, quote = F)



####################################################################################################
####################################################################################################

# PCA OF SVs

####################################################################################################
####################################################################################################
geo_initial <- readr::read_tsv("../../processed_data/genome_resources/elegans_isotypes_sampling_geo.tsv")
hawaii_islands <- readr::read_tsv("../../processed_data/genome_resources/elegans_isotypes_sampling_geo_hawaii_islands.tsv") %>% dplyr::select(isotype,collection_island_Hawaii)
WSs <- readr::read_tsv("../../processed_data/genome_resources/wild_strain_genome_stats.tsv") %>% dplyr::select(Strain) %>% dplyr::rename(strain = Strain) %>% dplyr::pull()

# Adding Hawaiian island resolution
geo <- geo_initial %>%
  dplyr::left_join(hawaii_islands, by = "isotype") %>%
  dplyr::mutate(geo = ifelse(geo == "Hawaii",collection_island_Hawaii,geo)) %>%
  dplyr::select(isotype, lat, long, geo) %>%
  dplyr::filter(isotype %in% WSs)

# Filter for MAF > 0.05
n <- 141
maf <- 0.05
least <- ceiling(maf * n)      #8
most  <- floor((1 - maf) * n)  #133

common_vcf <- merged_SV %>% 
  dplyr::filter(number_svs_merged >= least & number_svs_merged <= most) %>%
  dplyr::select(-chrom, -pos, -ref, -alt, -sv_type, -sv_length, -number_svs_merged) %>%
  as.matrix()

common_vcf[common_vcf == "./."] <- 0
sv_mat <- apply(common_vcf, 2, as.numeric)

sv_mat_t <- t(sv_mat)

sv_mat_t[is.na(sv_mat_t)] <- colMeans(sv_mat_t, na.rm = TRUE)

sv_pca <- prcomp(sv_mat_t, center = TRUE, scale. = TRUE) 


pca_df <- as.data.frame(sv_pca$x)
pca_df$strain <- rownames(pca_df)

strain_geo <- geo %>% dplyr::rename(strain = isotype) %>% dplyr::select(strain,geo)

pca_df <- pca_df %>%
  dplyr::left_join(strain_geo, by = "strain") %>%
  dplyr::mutate(geo = ifelse(strain == "CGC1", "CGC1",geo)) 

geo.colors <- c("Big Island"="black", "Molokai" = "#66C2A5", "Maui" = "yellow", "Oahu" = "brown", "Kauai" = "purple", "Africa"="green", "North America" = "pink", "Europe" = "#E41A1C", "Atlantic" = "blue", 
                "Oceania" ="cyan", "unknown" = 'gray', "CGC1" = "#DB6333")

pca_df <- pca_df %>%
  mutate(label = ifelse(PC2 > 50, strain, NA))

PCA <- ggplot(pca_df, aes(PC1, PC2, color = geo)) +
  geom_text_repel(aes(label = label), size = 2, max.overlaps = Inf, show.legend = FALSE) +
  geom_point(size = 1, alpha = 0.8) +
  scale_color_manual(values = geo.colors) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 10, color = 'black'),
    axis.title = element_text(size = 11, color = 'black'),
    legend.text = element_text(size = 9, color = 'black'),
    legend.title = element_text(size = 9, color = 'black'),
    legend.key.width = unit(0.1, "cm"),
    legend.key.spacing.x = unit(0, "pt"),
    legend.position = 'right') +
  labs(color = "Collection location", x = paste0("PC1 (", round(100 * summary(sv_pca)$importance[2,1], 1), "%)"),y = paste0("PC2 (", round(100 * summary(sv_pca)$importance[2,2], 1), "%)"))+
  guides(color = guide_legend(override.aes = list(size = 1), keyheight = unit(0.2, "cm")), keywidth = unit(0.2, "cm")) 
PCA


####################################################################################################
####################################################################################################

# CIRCOS PLOT

####################################################################################################
####################################################################################################
snps <- readr::read_tsv("../../processed_data/genome_resources/140WSs_biallelicSNPs.tsv", col_names = c("chrom","pos","ref","alt")) 
merged_SV <- readr::read_tsv("../../processed_data/structural_variants/Jasmine_merged_SVs.tsv")

# Outer ring of chromosomes (gene map of gene models represented with black rectangles)
## Chromosome IDs and sizes (start is always equal to zero)
chr_order <- c("I","II","III","IV","V","X")
chrom_sizes <- readr::read_tsv("../../processed_data/genome_resources/N2.WS283.cleaned.fa.fai", col_names = c("chrom","start","end")) %>%
  dplyr::mutate(chrom = factor(chrom, levels = chr_order)) %>%
  dplyr::mutate(chrom = as.character(chrom))

## N2 gene modesl in BED format
n2_genes_bed <- n2_genes_plt %>% dplyr::filter(chrom != "MtDNA") %>%
  transmute(chrom = as.character(chrom), start = as.numeric(start), end = as.numeric(end))

gene_bin_size <- 50000L

# one row per chr with lengths
chr_len_df <- chrom_sizes %>%
  dplyr::filter(chrom %in% chr_order) %>%
  dplyr::transmute(chrom = as.character(chrom), chr_len = as.numeric(end)) %>%
  dplyr::distinct()

# bins that cover the entire chromosome, including the last partial bin
gene_bins <- chr_len_df %>%
  dplyr::group_by(chrom) %>%
  dplyr::do({
    L <- .$chr_len[1]
    starts <- seq(0, L, by = gene_bin_size)  # include L so last bin is created
    tibble(
      chrom = .$chrom[1],
      start = starts[-length(starts)],
      end   = pmin(starts[-1], L)
    )
  }) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(mid = (start + end)/2)

# overlap-count genes per bin
bins_dt  <- as.data.table(gene_bins)
genes_dt <- as.data.table(n2_genes_bed)  # chrom/start/end already numeric in your code

setkey(bins_dt,  chrom, start, end)
setkey(genes_dt, chrom, start, end)

ov <- foverlaps(genes_dt, bins_dt, nomatch = 0)

gene_counts <- ov[, .(gene_count = .N), by = .(chrom, start, end, mid)]

bins_gene <- merge(bins_dt, gene_counts, by = c("chrom","start","end","mid"), all.x = TRUE)
bins_gene[is.na(gene_count), gene_count := 0]

gene_bins_50kb <- as.data.frame(bins_gene) %>%
  dplyr::rename(pos = mid, value = gene_count)

gene_ylim <- c(0, max(gene_bins_50kb$value, na.rm = TRUE))

# Next ring in will have SNV count per kb (plot fitted LOESS line?) 
## SNV count in table: "chrom", "bin", "variant_count"
bin_size <- 1000L

bins <- chrom_sizes %>%
  dplyr::group_by(chrom) %>%
  dplyr::do({
    chr_len <- .$end[1]
    starts <- seq(0, chr_len - 1, by = bin_size)
    tibble(
      chrom = .$chrom[1],
      start = starts,
      end   = pmin(starts + bin_size, chr_len)
    )
  }) %>%
  dplyr::ungroup()

bins_dt <- as.data.table(bins)
setkey(bins_dt, chrom, start, end)
bins_dt[, id := .I]

snps_counts <- snps %>%
  dplyr::select(-ref,-alt) %>%
  dplyr::filter(chrom %in% chr_order) %>%
  dplyr::mutate(
    chrom = as.character(chrom),
    start = (pos %/% bin_size) * bin_size,
    end   = start + bin_size) %>%
  dplyr::count(chrom, start, end, name = "variant_count")

snps_dt <- as.data.table(snps_counts)
setkey(snps_dt, chrom, start, end)

bins_snp <- merge(bins_dt, snps_dt, by = c("chrom","start","end"), all.x = TRUE)
bins_snp[is.na(variant_count), variant_count := 0]
snps_per_bin <- as.data.frame(bins_snp) %>%
  dplyr::mutate(pos = (start + end)/2) %>%
  dplyr::select(chrom,pos,variant_count) %>% dplyr::rename(value = variant_count)

# Next ring in will plot DEL frequncy (plotted as LOESS line) - make sure CGC1 is filtered out
## DEL calls with "chrom", "middle" (of the 1 kb bin), and "freq"
# Calculating DEL frequency
# bins <- snps %>% dplyr::select(chrom, bin)%>% dplyr::group_by(chrom) %>% dplyr::mutate(binEnd = lead(bin)) %>% dplyr::slice(-dplyr::n()) %>% dplyr::ungroup() %>% dplyr::rename(start = bin, end = binEnd)
bins_dt <- as.data.table(bins)
bins_dt[, id := .I]  # optional: keep track of bins

del_calls <- filt_calls %>% dplyr::filter(sv_type == "DEL", strain != "CGC1") %>% dplyr::mutate(end = pos + sv_length) %>% dplyr::rename(start = pos) %>% dplyr::select(chrom,start,end,strain)
del_calls_dt <- as.data.table(del_calls)

setkey(bins_dt, chrom, start, end)
setkey(del_calls_dt, chrom, start, end)

overlaps <- data.table::foverlaps(del_calls_dt, bins_dt, nomatch = 0)

# Count unique strains per bin
counts_PB <- overlaps[, .(n_strains = uniqueN(strain)), by = .(chrom, start, end)]

# If you want to merge with the full bin list (including 0s):
bins_wCounts <- merge(bins_dt, counts_PB, by = c("chrom", "start", "end"), all.x = TRUE)
bins_wCounts[is.na(n_strains), n_strains := 0]

bins_wFreq <- as.data.frame(bins_wCounts) %>%
  dplyr::mutate(freq = n_strains/141) #change me to number of isotypes

del_bin_plt <- bins_wFreq %>% dplyr::mutate(middle = (end + start) / 2)

del_bin_freq <- del_bin_plt %>% dplyr::select(chrom, middle, freq)

# Next ring in will plot INS frequency (plotted as LOESS line) - make sure CGC1 is filtered out
## INS calls with "chrom", "middle" (of the 1 kb bin), and "freq"
bins_dt <- as.data.table(bins)
bins_dt[, id := .I]  # optional: keep track of bins

ins_calls <- filt_calls %>% dplyr::filter(sv_type == "INS", strain != "CGC1") %>% dplyr::mutate(end = pos + sv_length) %>% dplyr::rename(start = pos) %>% dplyr::select(chrom,start,end,strain)
ins_calls_dt <- as.data.table(ins_calls)

setkey(bins_dt, chrom, start, end)
setkey(ins_calls_dt, chrom, start, end)

overlaps <- data.table::foverlaps(ins_calls_dt, bins_dt, nomatch = 0)

# Count unique strains per bin
counts_PB <- overlaps[, .(n_strains = uniqueN(strain)), by = .(chrom, start, end)]

# If you want to merge with the full bin list (including 0s):
bins_wCounts <- merge(bins_dt, counts_PB, by = c("chrom", "start", "end"), all.x = TRUE)
bins_wCounts[is.na(n_strains), n_strains := 0]

bins_wFreq <- as.data.frame(bins_wCounts) %>%
  dplyr::mutate(freq = n_strains/141) #change me to number of isotypes

ins_bin_plt <- bins_wFreq %>% dplyr::mutate(middle = (end + start) / 2)

ins_bin_freq <- ins_bin_plt %>% dplyr::select(chrom, middle, freq)

# Then the final, more inner ring will have INV frequency (plotted as LOESS line) - make sure CGC1 is filtered out
## INV calls with "chrom", "middle" (of the 1 kb bin), and "freq"
bins_dt <- as.data.table(bins)
bins_dt[, id := .I]  # optional: keep track of bins

inv_calls <- filt_calls %>% dplyr::filter(sv_type == "INV", strain != "CGC1") %>% dplyr::mutate(end = pos + sv_length) %>% dplyr::rename(start = pos) %>% dplyr::select(chrom,start,end,strain)
inv_calls_dt <- as.data.table(inv_calls)

setkey(bins_dt, chrom, start, end)
setkey(inv_calls_dt, chrom, start, end)

overlaps <- data.table::foverlaps(inv_calls_dt, bins_dt, nomatch = 0)

# Count unique strains per bin
counts_PB <- overlaps[, .(n_strains = uniqueN(strain)), by = .(chrom, start, end)]

# If you want to merge with the full bin list (including 0s):
bins_wCounts <- merge(bins_dt, counts_PB, by = c("chrom", "start", "end"), all.x = TRUE)
bins_wCounts[is.na(n_strains), n_strains := 0]

bins_wFreq <- as.data.frame(bins_wCounts) %>%
  dplyr::mutate(freq = n_strains/141) #change me to number of isotypes

inv_bin_plt <- bins_wFreq %>% dplyr::mutate(middle = (end + start) / 2)

inv_bin_freq <- inv_bin_plt %>% dplyr::select(chrom, middle, freq)


## =========================
##  Plotting!
## ========================= 
# SV frequency tracks must have: chrom, pos, value
del_bin_freq <- del_bin_freq %>%
  transmute(chrom = as.character(chrom), pos = as.numeric(middle), value = as.numeric(freq))

ins_bin_freq <- ins_bin_freq %>%
  transmute(chrom = as.character(chrom), pos = as.numeric(middle), value = as.numeric(freq))

inv_bin_freq <- inv_bin_freq %>%
  transmute(chrom = as.character(chrom), pos = as.numeric(middle), value = as.numeric(freq))

add_rect_track <- function(bed, col, track_height = 0.08, label = NULL) {
  circos.trackPlotRegion(
    ylim = c(0, 1),
    track.height = track_height,
    bg.border = NA,
    panel.fun = function(x, y) {
      chr <- CELL_META$sector.index
      d <- bed[bed$chrom == chr, , drop = FALSE]
      if (nrow(d) == 0) return()
      circos.rect(d$start, 0, d$end, 1, col = col, border = NA)
    }
  )
  
  if (!is.null(label)) {
    circos.text(
      x = 0, y = 0.5, labels = label,
      sector.index = chr_order[1],
      track.index  = get.current.track.index(),
      facing = "inside", adj = c(1, 0.5), cex = 0.7
    )
  }
}

add_value_track <- function(df, col, ylim, track_height = 0.10, label = NULL,
                            type = c("line","points"), add_loess = FALSE, span = 0.2) {
  type <- match.arg(type)
  
  circos.trackPlotRegion(
    ylim = ylim,
    track.height = track_height,
    bg.border = NA,
    panel.fun = function(x, y) {
      chr <- CELL_META$sector.index
      d <- df[df$chrom == chr, , drop = FALSE]
      if (nrow(d) == 0) return()
      
      d <- d[order(d$pos), , drop = FALSE]
      
      if (type == "points") {
        circos.points(d$pos, d$value, pch = 16, cex = 0.25, col = col)
      } else {
        circos.lines(d$pos, d$value, col = col, lwd = 1)
      }
      
      if (add_loess && nrow(d) >= 50) {
        fit <- stats::loess(value ~ pos, data = d, span = span)
        xs  <- d$pos
        ys  <- stats::predict(fit, newdata = data.frame(pos = xs))
        ok  <- is.finite(ys)
        if (any(ok)) circos.lines(xs[ok], ys[ok], col = col, lwd = 2)
      }
    }
  )
  
  if (!is.null(label)) {
    circos.text(
      x = 0, y = mean(ylim), labels = label,
      sector.index = chr_order[1],
      track.index  = get.current.track.index(),
      facing = "inside", adj = c(1, 0.5), cex = 0.7
    )
  }
}

add_filled_area_track <- function(df, col_fill = "#00BFC480", col_line = "#00BFC4",
                                  ylim, track_height = 0.10, label = NULL) {
  circos.trackPlotRegion(
    ylim = ylim,
    track.height = track_height,
    bg.border = NA,
    panel.fun = function(x, y) {
      chr <- CELL_META$sector.index
      d <- df[df$chrom == chr, , drop = FALSE]
      if (nrow(d) < 2) return()
      d <- d[order(d$pos), , drop = FALSE]
      
      # build polygon (baseline at ylim[1])
      x_poly <- c(d$pos, rev(d$pos))
      y_poly <- c(d$value, rep(ylim[1], nrow(d)))
      
      circos.polygon(x_poly, y_poly, col = col_fill, border = NA)
      circos.lines(d$pos, d$value, col = col_line, lwd = 1.2)
    }
  )
  
  if (!is.null(label)) {
    circos.text(
      x = 0, y = mean(ylim), labels = label,
      sector.index = chr_order[1],
      track.index  = get.current.track.index(),
      facing = "inside", adj = c(1, 0.5), cex = 0.7
    )
  }
}


png("circos_plot.png", width = 3.75, height = 3.75, units = "in", res = 600, bg = "white")

circos.clear()
circos.par(track.margin = c(0.002, 0.002), cell.padding = c(0, 0, 0, 0), canvas.xlim = c(-0.9, 0.9), canvas.ylim = c(-0.99,0.99)) # start.degree = 86, gap.after = c(rep(2, length(chr_order)-1), 8)
circos.initialize(factors = as.character(chrom_sizes$chrom), xlim = cbind(rep(0, nrow(chrom_sizes)), chrom_sizes$end))

# Outer ideogram-like ring + gene models inside it
circos.trackPlotRegion(
  ylim = c(0, 1),
  track.height = 0.02,
  bg.border = NA,
  panel.fun = function(x, y) {
    chr <- CELL_META$sector.index
    xlim <- CELL_META$xlim
    
    # background chromosome band
    circos.rect(xlim[1], 0, xlim[2], 1, col = "grey90", border = "black")
    # gene models as black rectangles (thin band)
    # g <- n2_genes_bed[n2_genes_bed$chrom == chr, , drop = FALSE]
    # if (nrow(g) > 0) {
    #   circos.rect(g$start, 0, g$end, 1, col = "black", border = NA)}
    
    # gene density (50 kb bins) as cyan-blue line
    # gd <- gene_density_50kb[gene_density_50kb$chrom == chr, , drop = FALSE]
    # if (nrow(gd) > 1) {
    #   gd <- gd[order(gd$pos), , drop = FALSE]
    #   
    #   # rescale density to fit nicely within the outer track band
    #   y <- (gd$value - gene_ylim[1]) / (gene_ylim[2] - gene_ylim[1] + 1e-9)
    #   y <- 0.10 + y * 0.80   # occupy 10%..90% of the track height
    #   
    #   circos.lines(gd$pos, y, col = "#00BFC4", lwd = 2)
    # }
    
    
    
    # chromosome label
    circos.text(CELL_META$xcenter, 3.3, chr, facing = "bending.outside", niceFacing = T, cex = 1.5, font = 2)
    # axis ticks
    circos.axis(h = "top", major.at = seq(0, xlim[2], by = 5e6),labels = seq(0, xlim[2], by = 5e6) / 1e6, labels.cex = 1, major.tick.length = 0.001)
  }
)


snp_ylim <- c(0, max(snps_per_bin$value, na.rm = TRUE))
add_value_track(snps_per_bin, col = "#DB6333", ylim = snp_ylim, track_height = 0.2, type = "line", add_loess = FALSE, span = 0.2)
add_value_track(del_bin_freq, col = "red",  ylim = c(0, 1), track_height = 0.2, type = "line", add_loess = FALSE, span = 0.2)
add_value_track(ins_bin_freq, col = "blue", ylim = c(0, 1), track_height = 0.2, type = "line", add_loess = FALSE, span = 0.2)
add_value_track(inv_bin_freq, col = "gold3", ylim = c(0, 1), track_height = 0.2, type = "line", add_loess = FALSE, span = 0.2)

# Legend
lgd <- Legend(labels = c("SNPs/kb","DEL freq.","INS freq.","INV freq."),
              type = "lines", legend_gp = gpar(col = c("#DB6333","red","blue","gold3"),
                                              lwd = c(3, 3, 3, 3)),
              labels_gp = gpar(fontsize = 10),
              title_gp  = gpar(fontsize = 10))

draw(lgd, x = unit(0.99, "npc"), y = unit(0.99, "npc"), just = c("right", "top"))

circos.clear()

dev.off()

CIRC <- ggdraw() + draw_image("circos_plot.png")



# Create final plot
top_plt <- cowplot::plot_grid(
  N2_EXP, SV_LEN,
  nrow = 1,
  labels = c("a","b"))

bot_plt <- cowplot::plot_grid(
  PCA, CIRC,
  nrow = 1,
  rel_widths = c(1.5,1),
  labels = c("c","d"))

final_plt <- cowplot::plot_grid(
  top_plt, bot_plt,
  nrow = 2) + theme(plot.background = element_rect(fill = "white", color = NA))
final_plt

# Save the plot
ggsave("../../figures/structural_variants.png", final_plt, width = 7.5, height = 7.5, dpi = 600)
