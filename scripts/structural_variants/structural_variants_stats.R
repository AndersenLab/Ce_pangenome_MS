library(readr)
library(ggplot2)
library(dplyr)
library(data.table)
library(GenomicRanges)

# Load in all SV calls
allcalls <- readr::read_tsv("../../processed_data/structural_variants/141_over50_PASS_variants.tsv", col_names = c("chrom", "pos", "ref", "alt", "filter", "sv_type","sv_length","strain")) %>% dplyr::select(-filter)

# Looking at the count of each type of SV across all wild strains
filt_calls <- allcalls %>% 
  dplyr::mutate(sv_length = abs(sv_length)) %>%
  dplyr::group_by(sv_type, strain) %>%
  dplyr::mutate(type_count = n()) %>%
  dplyr::ungroup()

levels <- summary <- filt_calls %>%
  dplyr::select(sv_type,strain,type_count) %>%
  dplyr::distinct() %>%
  dplyr::group_by(strain) %>%
  dplyr::arrange(desc(type_count)) %>%
  dplyr::ungroup() %>% 
  dplyr::distinct(strain) %>%
  dplyr::pull()

summary <- filt_calls %>%
  dplyr::select(sv_type,strain,type_count) %>%
  dplyr::distinct() %>%
  dplyr::group_by(strain) %>%
  dplyr::arrange(desc(type_count)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(sv_type = factor(sv_type, levels = c("INS","DEL","INV"))) %>%
  dplyr::mutate(strain = factor(strain, levels = levels))

count_stats <- summary %>%
  dplyr::group_by(sv_type) %>%
  dplyr::summarise(average_count = mean(type_count)) %>%
  dplyr::mutate(label = paste0(sv_type, " mean: ", round(average_count, 0))) %>%
  dplyr::mutate(ypos = max(summary$type_count) * c(0.95, 0.91, 0.87)) # stats of the everage number of SV calls among all wild strains

sv_count <- ggplot() +
  geom_col(data = summary, aes(x = strain, y = type_count, fill = sv_type), position = position_dodge(width = 0.8)) +
  geom_hline(data = count_stats, aes(yintercept = average_count, color = sv_type), linetype = "dashed", linewidth = 0.5) +
  scale_fill_manual(values = c("INS" = "blue", "DEL" = "red", "INV" = "gold")) +
  scale_color_manual(values = c("INS" = "blue", "DEL" = "red", "INV" = "gold")) +
  labs(y = "Count of SV type for wild strains", fill = "SV type") +
  theme(panel.border = element_rect(color = 'black', fill = NA),
        panel.background = element_blank(),
        panel.grid.major= element_line(color = 'gray80'),
        panel.grid.major.x = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.85,0.815),
        legend.title = element_text(size = 8, color = 'black'),
        legend.text = element_text(size = 8, color = 'black'),
        legend.box.background = element_rect(color = 'black', fill = NA),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 9, color = 'black'),
        axis.text.x = element_text(color = 'black', size = 4, angle = 70, hjust = 1),
        axis.title.y = element_text(size = 12, color = 'black')) +
  guides(color = "none") + # to get rid of legend for the horizontal lines
  scale_y_continuous(expand = c(0,0))

# Looking at the proportion of each type of SV for every wild strain
prop <- summary %>%
  dplyr::group_by(strain) %>%
  dplyr::mutate(total = sum(type_count)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(sv_type,strain) %>%
  dplyr::mutate(prop = type_count / total * 100) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(sv_type = factor(sv_type, levels = c("INV","DEL","INS"))) 

prop_plt <- ggplot() +
  geom_bar(data = prop, aes(x = strain, y = prop, fill = sv_type), stat = "identity") +
  scale_fill_manual(values = c("INS" = "blue", "DEL" = "red", "INV" = "gold")) +
  labs(y = "Proportion of SVs", fill = "SV type") +
  theme(panel.border = element_rect(color = 'black', fill = NA),
        panel.background = element_blank(),
        panel.grid.major= element_line(color = 'gray80'),
        panel.grid.major.x = element_blank(),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.text = element_text(size = 9, color = 'black'),
        axis.text.x = element_text(color = 'black', size = 4, angle = 70, hjust = 1),
        axis.title.y = element_text(size = 12, color = 'black')) +
  guides(color = "none") + # to get rid of legend for the horizontal lines
  scale_y_continuous(expand = c(0,0))
prop_plt

final_plt <- cowplot::plot_grid(
  sv_count, prop_plt,
  nrow = 2,
  align = "v",
  labels = c("a","b"))
final_plt


# Save the plot
ggsave("../../figures/supplementary/structural_variant_countProportion.png", final_plt, width = 7.5, height = 8, dpi = 600)




#############################################################################
# How does SV count correlate with SNP count?
#############################################################################
sv_strain_count <- filt_calls %>% dplyr::distinct(strain,sv_type,type_count) %>% dplyr::group_by(strain) %>% dplyr::mutate(all_sv_count = sum(type_count)) %>% dplyr::distinct(strain,all_sv_count) %>%
  dplyr::filter(strain != "CGC1") # we do not have SNP calls for CGC1

sv_strain_type_count <- filt_calls %>% dplyr::distinct(strain,sv_type,type_count) %>%
  dplyr::filter(strain != "CGC1")

snps_per_strain <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/misc/henrique_lab/140_strains_ALT_SNP_count.tsv", col_names = c("strain","snp_count"))

sv_snp_strains <- sv_strain_count %>% dplyr::left_join(snps_per_strain, by = "strain") 

# Correlation plot of SV count with SNP count with founder strains colored differently than the rest
sv_snp_strains <- sv_snp_strains %>% dplyr::arrange(snp_count) %>% dplyr::mutate(strain = factor(strain, levels = strain)) 

model <- lm(snp_count ~ all_sv_count, data = sv_snp_strains)

# View summary (slope, intercept, R-squared)
summary(model)

# Extract slope and intercept
coef(model)

slope <- round(coef(model)[2], 1)
r_squared <- round(summary(model)$r.squared, 3)

snp_sv_corr <- ggplot(sv_snp_strains, aes(x = all_sv_count, y = snp_count)) +
  geom_point(size = 1) +
  geom_smooth(method = "lm", se = TRUE, color = "blue", linetype = "dashed", alpha = 0.2, linewidth = 0.5) + # minimizes the sum of squared residuals, or ordinary least squares
  annotate("text", x = 6000, y = 3.5e5, label = paste0("slope = ", slope, "\nR² = ", r_squared), size = 6, hjust = 0) +
  labs(x = "SV count", y = "SNP count") +
  theme(
    panel.border = element_rect(color = 'black', fill = NA),
    panel.background = element_blank(),
    panel.grid.major = element_line(color = 'gray80'),
    panel.grid.major.x = element_blank(),
    legend.position = "none",
    axis.text = element_text(color = 'black', size = 10),
    axis.title = element_text(size = 12, color = 'black'))
snp_sv_corr

# Save the plot
ggsave("../../figures/supplementary/sv_snp_correlation.png", snp_sv_corr, width = 7.5, height = 7.5, dpi = 600)





#############################################################################
# The largest PAV-identified inversion
#############################################################################
largest <- filt_calls %>%
  dplyr::select(chrom,pos,sv_type, sv_length, strain) %>%
  dplyr::group_by(sv_type) %>%
  dplyr::arrange(desc(sv_length)) %>%
  dplyr::slice_head(n = 4) 

soi <- largest %>% 
  dplyr::distinct(strain) %>%
  dplyr::pull()

nucmer <- readr::read_tsv("../../processed_data/genome_resources/genome_data/141_nucmer_ECA741CGC1.tsv", col_names = c("N2S","N2E","WSS","WSE","L1","L2","IDY","LENR","LENQ","N2_chr","contig","strain")) %>%
  dplyr::filter(strain %in% soi)

inv <- nucmer %>% dplyr::filter(strain == "ECA1413" | strain == "ECA2968" | strain == "XZ1515") %>% dplyr::filter(N2_chr == "IV") %>% 
  dplyr::filter(contig == "ptg000001l" | contig == "ptg000008l" | contig == "ptg000003l") # retaining the main contigs for each strain

inv_rect <- largest %>% dplyr::filter(sv_type == "INV") %>% dplyr::filter(strain != "JU3226") # isolating the largest INV that is shared among three strains from Kaua'i

largest_PAV_INV <- ggplot(inv) +
  geom_rect(data = inv_rect, aes(xmin = (pos + sv_length) / 1e6, xmax = pos / 1e6, ymin = -Inf, ymax = Inf), fill = "gold", alpha = 0.5) +
  geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6), linewidth = 0.7, color = 'black') +
  theme_bw() +
  facet_wrap(~strain) +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 14, color = 'black'),
    axis.text = element_text(size = 10, color = 'black'),
    axis.ticks = element_blank(),
    axis.title = element_text(size = 12, color = 'black'),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA)) +
  coord_cartesian(xlim = c(6, 8)) +
  labs(x = "N2 genome position (Mb)", y = "Wild strain contig position (Mb)")
largest_PAV_INV

# Save the plot
ggsave("../../figures/supplementary/largest_PAV_INV.png", largest_PAV_INV, width = 7.5, height = 6, dpi = 600)



#############################################################################
# Are SVs enriched in HDRs??
#############################################################################
strains <- allcalls %>% dplyr::select(strain) %>% dplyr::distinct() %>% dplyr::pull()
hdrs <- readr::read_tsv("../../processed_data/genome_resources/annotation/20250625_c_elegans_divergent_regions_strain.bed", col_names = c("chrom", "start", "end", "strain")) %>% 
  dplyr::filter(strain %in% strains)

# COLLAPSING HDRS AMONG 140 WSs
all_regions <- hdrs %>%
  dplyr::rename(START = start, CHROM = chrom, END = end) %>%
  dplyr::arrange(CHROM,START) %>%
  dplyr::group_split(CHROM)

strain_count <- hdrs %>% dplyr::distinct(strain, .keep_all = T)
print(nrow(strain_count)) # SHOULD BE 140

# Collapsing all HDRs
getRegFreq <- function(all_regions) {
  all_collapsed <- list()
  for (i in 1:length(all_regions)) {
    temp <- all_regions[[i]]
    k=1
    j=1
    while (k==1) {
      print(paste0("chrom:",i,"/iteration:",j))
      checkIntersect <- temp %>%
        dplyr::arrange(CHROM,START) %>%
        dplyr::mutate(check=ifelse(lead(START) <= END,T,F)) %>%
        dplyr::mutate(check=ifelse(is.na(check),F,check))
      
      #print(nrow(checkIntersect %>% dplyr::filter(check==T)))
      
      if(nrow(checkIntersect %>% dplyr::filter(check==T)) == 0) {
        print("NO MORE INTERSECTS")
        k=0
      } else {
        
        temp <- checkIntersect %>%
          dplyr::mutate(gid=data.table::rleid(check)) %>%
          dplyr::mutate(gid=ifelse((check==F| is.na(check)) & lag(check)==T,lag(gid),gid))
        
        collapse <- temp %>%
          dplyr::filter(check==T | (check==F & lag(check)==T)) %>%
          dplyr::group_by(gid) %>%
          dplyr::mutate(newStart=min(START)) %>%
          dplyr::mutate(newEnd=max(END)) %>%
          dplyr::ungroup() %>%
          dplyr::distinct(gid,.keep_all = T)  %>%
          dplyr::mutate(START=newStart,END=newEnd) %>%
          dplyr::select(-newEnd,-newStart)
        
        retain <- temp %>%
          dplyr::filter(check==F & lag(check)==F)
        
        temp <- rbind(collapse,retain) %>%
          dplyr::select(-gid,-check)
        
        j=j+1
      }
    }
    print(head(temp))
    all_collapsed[[i]] <- temp
  }
  return(all_collapsed)
}

HDR_collapse_master <- getRegFreq(all_regions)

all_collapsed <- plyr::ldply(HDR_collapse_master, data.frame) %>%
  dplyr::select(-strain)

colnames(all_collapsed) <- c("chrom","start","end")



# Calculating DEL frequency
bin_size <- 1000L

chr_order <- c("I","II","III","IV","V","X")

chrom_sizes <- readr::read_tsv("../../data/N2.WS283.cleaned.fa.fai", col_names = c("chrom","start","end")) %>%
  dplyr::mutate(chrom = factor(chrom, levels = chr_order)) %>%
  dplyr::mutate(chrom = as.character(chrom))

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

del_freq <- ggplot(del_bin_plt) +
  geom_rect(data = all_collapsed, aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf), fill = 'gray', alpha = 0.7) +
  geom_point(aes(x = middle, y = freq), color = 'red', size = 0.3) +
  facet_wrap(~chrom, nrow = 1, scales = "free_x") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text.y = element_text(size = 10, color = 'black'),
    # panel.grid = element_blank(),
    panel.border = element_rect(color = 'black', fill = NA),
    panel.background = element_blank(),
    strip.text = element_text(size = 14, color = "black")) +
  coord_cartesian(ylim = c(0,1.01)) +
  scale_y_continuous(expand = c(0,0))


# Calculating INS frequency
# bins <- snps %>% dplyr::select(chrom, bin)%>% dplyr::group_by(chrom) %>% dplyr::mutate(binEnd = lead(bin)) %>% dplyr::slice(-dplyr::n()) %>% dplyr::ungroup() %>% dplyr::rename(start = bin, end = binEnd)
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

ins_freq <- ggplot(ins_bin_plt) +
  geom_rect(data = all_collapsed, aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf), fill = 'gray', alpha = 0.7) +
  geom_point(aes(x = middle, y = freq), color = 'blue', size = 0.3) +
  facet_wrap(~chrom, nrow = 1, scales = "free_x") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text.y = element_text(size = 10, color = 'black'),
    panel.border = element_rect(color = 'black', fill = NA),
    panel.background = element_blank(),
    strip.text = element_blank()) +
  coord_cartesian(ylim = c(0,1.01)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 1), breaks = seq(0, 0.75, 0.25))



# Calculating INV frequency
# bins <- snps %>% dplyr::select(chrom, bin)%>% dplyr::group_by(chrom) %>% dplyr::mutate(binEnd = lead(bin)) %>% dplyr::slice(-dplyr::n()) %>% dplyr::ungroup() %>% dplyr::rename(start = bin, end = binEnd)
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

inv_freq <- ggplot(inv_bin_plt) +
  geom_rect(data = all_collapsed, aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf), fill = 'gray', alpha = 0.7) +
  geom_point(aes(x = middle, y = freq), color = 'gold1', size = 0.3) +
  facet_wrap(~chrom, nrow = 1, scales = "free_x") +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_text(size = 14, color = 'black'),
    axis.text.y = element_text(size = 10, color = 'black'),
    strip.text = element_blank(),
    panel.border = element_rect(color = 'black', fill = NA),
    panel.background = element_blank()) +
  xlab("N2 genomic position (Mb)") +
  coord_cartesian(ylim = c(0,1.01)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 1), breaks = seq(0, 0.75, 0.25))
inv_freq
  
# Concatenate together all three SV types!
all_three_SVs <- cowplot::plot_grid(
  del_freq, ins_freq, inv_freq,
  nrow = 3,
  align = "v",
  rel_heights = c(1,1,1))

# Save the figure
ggsave("../../figures/supplementary/SV_HDR_enrichment.png", all_three_SVs, width = 7.5, height = 6, dpi = 600)


# Looking at enrichment of SVs in HDRs
hdr_bins <- all_collapsed
del_calls <- filt_calls %>% dplyr::filter(sv_type == "DEL", strain != "CGC1") %>% dplyr::mutate(end = pos + sv_length) %>% dplyr::rename(start = pos) %>% dplyr::select(chrom,start,end,strain)

# helper: tibble (BED-like, 0-based) -> GRanges (1-based)
to_gr <- function(df) {
  GRanges(
    seqnames = df$chrom,
    ranges   = IRanges(start = df$start + 1, end = df$end)
  )
}

genome_gr <- to_gr(chrom_sizes)          # whole genome space
hdr_gr    <- to_gr(hdr_bins)  
del_gr    <- to_gr(del_calls)            # DEL intervals

# de-duplicate identical DELs (non-rare HDRs) because we are assessing genomic loci affected by HDRs in or outside of HDRs, not the frequency that a genomic
# locus is affected by a DEL
# Keeping duplicated DEL calls - "Among all DEL calls across 140 strain, are calls more frequent per bp in HDRs than outside?"
# De-duplicating DEL calls - "Are genomic loci that harbor deletions enriched in HDRs vs outside?"
del_gr_uniq <- unique(del_gr)

# count DEL events overlapping HDR
in_hdr <- countOverlaps(del_gr_uniq, hdr_gr) > 0
n_hdr <- sum(in_hdr) # number of DELs that overlap with HDRs
n_out <- length(del_gr_uniq) - n_hdr # number of DELs that do NOT overlap with HDRs

# bp denominators
hdr_bp    <- sum(width(hdr_gr)) # summed genome size of HDRs
genome_bp <- sum(width(genome_gr)) # genome size
nonhdr_bp <- genome_bp - hdr_bp # summed genome size outside of HDRs

# densities and fold enrichment
density_hdr <- n_hdr / hdr_bp # normalizing by bps - # of DELs in HDRs / # bps in HDRs
density_out <- n_out / nonhdr_bp

fold_enrichment <- density_hdr / density_out

# report
list(
  total_DEL = length(del_gr_uniq),
  DEL_in_HDR = n_hdr,
  DEL_out_HDR = n_out,
  HDR_bp = hdr_bp,
  nonHDR_bp = nonhdr_bp,
  density_HDR = density_hdr,
  density_nonHDR = density_out,
  fold_enrichment = fold_enrichment)


# INS enrichment in HDRs
hdr_bins <- all_collapsed
ins_calls <- filt_calls %>% dplyr::filter(sv_type == "INS", strain != "CGC1") %>% dplyr::mutate(end = pos + sv_length) %>% dplyr::rename(start = pos) %>% dplyr::select(chrom,start,end,strain)

to_gr <- function(df) {
  GRanges(
    seqnames = df$chrom,
    ranges   = IRanges(start = df$start + 1, end = df$end)
  )
}

genome_gr <- to_gr(chrom_sizes)          # whole genome space
hdr_gr    <- to_gr(hdr_bins)  
ins_gr    <- to_gr(ins_calls)            # ins intervals


ins_gr_uniq <- unique(ins_gr)

in_hdr <- countOverlaps(ins_gr_uniq, hdr_gr) > 0
n_hdr <- sum(in_hdr) # number of INSs that overlap with HDRs
n_out <- length(ins_gr_uniq) - n_hdr # number of INSs that do NOT overlap with HDRs

hdr_bp    <- sum(width(hdr_gr)) # summed genome size of HDRs
genome_bp <- sum(width(genome_gr)) # genome size
nonhdr_bp <- genome_bp - hdr_bp # summed genome size outside of HDRs

density_hdr <- n_hdr / hdr_bp # normalizing by bps - # of INSs in HDRs / # bps in HDRs
density_out <- n_out / nonhdr_bp

fold_enrichment <- density_hdr / density_out

# report
list(
  total_ins = length(ins_gr_uniq),
  ins_in_HDR = n_hdr,
  ins_out_HDR = n_out,
  HDR_bp = hdr_bp,
  nonHDR_bp = nonhdr_bp,
  density_HDR = density_hdr,
  density_nonHDR = density_out,
  fold_enrichment = fold_enrichment)


# INV enrichment in HDRs
hdr_bins <- all_collapsed
INV_calls <- filt_calls %>% dplyr::filter(sv_type == "INV", strain != "CGC1") %>% dplyr::mutate(end = pos + sv_length) %>% dplyr::rename(start = pos) %>% dplyr::select(chrom,start,end,strain)

to_gr <- function(df) {
  GRanges(
    seqnames = df$chrom,
    ranges   = IRanges(start = df$start + 1, end = df$end)
  )
}

genome_gr <- to_gr(chrom_sizes)          # whole genome space
hdr_gr    <- to_gr(hdr_bins)  
INV_gr    <- to_gr(INV_calls)            # ins intervals


INV_gr_uniq <- unique(INV_gr)

in_hdr <- countOverlaps(INV_gr_uniq, hdr_gr) > 0
n_hdr <- sum(in_hdr) # number of INSs that overlap with HDRs
n_out <- length(INV_gr_uniq) - n_hdr # number of INSs that do NOT overlap with HDRs

hdr_bp    <- sum(width(hdr_gr)) # summed genome size of HDRs
genome_bp <- sum(width(genome_gr)) # genome size
nonhdr_bp <- genome_bp - hdr_bp # summed genome size outside of HDRs

density_hdr <- n_hdr / hdr_bp # normalizing by bps - # of INSs in HDRs / # bps in HDRs
density_out <- n_out / nonhdr_bp

fold_enrichment <- density_hdr / density_out

# report
list(
  total_ins = length(INV_gr_uniq),
  INV_in_HDR = n_hdr,
  INV_out_HDR = n_out,
  HDR_bp = hdr_bp,
  nonHDR_bp = nonhdr_bp,
  density_HDR = density_hdr,
  density_nonHDR = density_out,
  fold_enrichment = fold_enrichment)






















#############################################################################
# SV overlap with protein-coding genes
#############################################################################
# Loading in N2 genes and merged SV calls from PAV and Jasmine
merged_SV <- readr::read_tsv("../../processed_data/structural_variants/Jasmine_merged_SVs.tsv")
N2_gff <- ape::read.gff("../../processed_data/genome_resources/annotation/c_elegans.PRJNA13758.WS283.csq.PCfeaturesOnly.longest.gff3") 
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

sv_gene_overlap <- ggplot() +
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
sv_gene_overlap

# Save the plot:
ggsave("")




### Overlap with an N2 gene or 2kb upstream (encompases promoter region)
twokb_n2_genes_plt <- n2_genes_plt %>% dplyr::mutate(start = start - 2000)
svs_dt <- as.data.table(maf_filt)
n2_genes_dt <- as.data.table(twokb_n2_genes_plt)

setkey(svs_dt, chrom, start, end)
setkey(n2_genes_dt, chrom, start, end)

svs_inCodingRegions <- data.table::foverlaps(x = svs_dt, y = n2_genes_dt, type = "any") %>% dplyr::filter(!is.na(start)) %>% dplyr::mutate(overlap = T)

overlap <- svs_inCodingRegions %>% dplyr::select(chrom, i.start,i.end, sv_type, overlap) %>% dplyr::rename(start = i.start, end = i.end) %>% dplyr::distinct(chrom,start,end,sv_type, .keep_all = T)

final_stats <- maf_filt %>% 
  dplyr::left_join(overlap, by = c("chrom", "start", "end", "sv_type")) %>% 
  dplyr::mutate(region = ifelse(is.na(overlap.y),'non-PC_region','overlaps_PCgene_(plus2kbupstream)')) %>%
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


ggplot() +
  geom_bar(data = plt_stats, aes(x = sv_type, y = proportion, fill = region), stat = "identity") +
  geom_text(data = plt_stats, aes(x = sv_type, y = proportion, label = region_count, group = region), position = position_stack(vjust = 0.5), color = "white", size = 4, fontface = "bold") +
  scale_fill_manual(values = c("overlaps_PCgene_(plus2kbupstream)" = "firebrick", "non-PC_region" = "gray60")) +
  labs(y = "Proportion (%)", fill = "Region") +
  theme(panel.border = element_rect(color = 'black', fill = NA),
        panel.background = element_blank(),
        panel.grid.major= element_line(color = 'gray80'),
        panel.grid.major.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 12, color = 'black'),
        axis.text.x = element_text(color = 'black', size = 12, face = 'bold'),
        axis.title.y = element_text(size = 14, color = 'black', face = 'bold')) +
  guides(color = "none") + # to get rid of legend for the horizontal lines
  scale_y_continuous(expand = c(0,0))

















# ========================================================================================= #
# PCA of SVs
# ========================================================================================= #
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

geo.colors <- c("Big Island"="#66C2A5", "Molokai" = "black", "Maui" = "yellow", "Oahu" = "brown", "Kauai" = "pink", "Africa"="green", "North America" = "purple", "Europe" = "#E41A1C", "Atlantic" = "blue", 
                "Oceania" ="cyan", "unknown" = 'gray', "CGC1" = "#DB6333")

pca_df <- pca_df %>%
  mutate(label = ifelse(PC2 > 50, strain, NA))

ggplot(pca_df, aes(PC1, PC2, color = geo)) +
  geom_text_repel(aes(label = label), size = 4, max.overlaps = Inf, show.legend = FALSE) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = geo.colors) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 12, color = 'black'),
    axis.title = element_text(size = 12, color = 'black', face = 'bold'),
    legend.text = element_text(size = 14, color = 'black'),
    legend.title = element_text(size = 14, color = 'black'),
    plot.title = element_text(size = 16, color = 'black', face = 'bold', hjust = 0.5)) +
  labs(color = "Collection location",title = "All SVs", x = paste0("PC1 (", round(100 * summary(sv_pca)$importance[2,1], 1), "%)"),y = paste0("PC2 (", round(100 * summary(sv_pca)$importance[2,2], 1), "%)"))+
  guides(color = guide_legend(override.aes = list(size = 7))) 


dim(sv_mat_t)
summary(rowSums(sv_mat_t))
cor(rowSums(sv_mat_t), sv_pca$x[,1])  # PC1 is 93% correlated with number of SVs per strain


# PC3/4
ggplot(pca_df, aes(PC3, PC4, color = geo)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = geo.colors) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 12, color = 'black'),
    axis.title = element_text(size = 12, color = 'black', face = 'bold'),
    legend.text = element_text(size = 14, color = 'black'),
    legend.title = element_text(size = 14, color = 'black'),
    plot.title = element_text(size = 16, color = 'black', face = 'bold', hjust = 0.5)) +
  labs(color = "Collection location",title = "All SVs", x = paste0("PC3 (", round(100 * summary(sv_pca)$importance[2,3], 1), "%)"),y = paste0("PC4 (", round(100 * summary(sv_pca)$importance[2,4], 1), "%)"))+
  guides(color = guide_legend(override.aes = list(size = 7))) 

# Scree plot:
scree <- data.frame(PC = seq_along(sv_pca$sdev),variance = sv_pca$sdev^2 / sum(sv_pca$sdev^2))

ggplot(scree[1:10,], aes(PC, variance * 100)) +
  geom_col(fill = "black") +
  theme(
    axis.text = element_text(size = 16, color = 'black'),
    panel.background = element_blank(),
    plot.margin = margin(10,10,10,10),
    panel.border = element_rect(color = 'black', fill = NA),
    axis.title = element_text(size = 16, color = 'black', face = 'bold')) +
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(breaks = 1:10) +
  coord_cartesian(ylim = c(0,16)) +
  labs(y = "Proportion of variance explained (%)", x = "Principal component")


### DELETIONS ###
common_vcf <- merged_SV %>% 
  dplyr::filter(number_svs_merged >= least & number_svs_merged <= most & sv_type == "DEL") %>%
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

geo.colors <- c("Big Island"="#66C2A5", "Molokai" = "black", "Maui" = "yellow", "Oahu" = "brown", "Kauai" = "pink", "Africa"="green", "North America" = "purple", "Europe" = "#E41A1C", "Atlantic" = "blue", 
                "Oceania" ="cyan", "unknown" = 'gray', "CGC1" = "#DB6333")

pca_df <- pca_df %>%
  mutate(label = ifelse(PC2 > 50, strain, NA))

ggplot(pca_df, aes(PC1, PC2, color = geo)) +
  geom_text_repel(aes(label = label), size = 4, max.overlaps = Inf, show.legend = FALSE) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = geo.colors) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 12, color = 'black'),
    axis.title = element_text(size = 12, color = 'black', face = 'bold'),
    legend.text = element_text(size = 14, color = 'black'),
    legend.title = element_text(size = 14, color = 'black'),
    plot.title = element_text(size = 16, color = 'black', face = 'bold', hjust = 0.5)) +
  labs(color = "Collection location",title = "Deletions", x = paste0("PC1 (", round(100 * summary(sv_pca)$importance[2,1], 1), "%)"),y = paste0("PC2 (", round(100 * summary(sv_pca)$importance[2,2], 1), "%)"))+
  guides(color = guide_legend(override.aes = list(size = 7))) 


# PC3/4
ggplot(pca_df, aes(PC3, PC4, color = geo)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = geo.colors) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 12, color = 'black'),
    axis.title = element_text(size = 12, color = 'black', face = 'bold'),
    legend.text = element_text(size = 14, color = 'black'),
    legend.title = element_text(size = 14, color = 'black'),
    plot.title = element_text(size = 16, color = 'black', face = 'bold', hjust = 0.5)) +
  labs(color = "Collection location",title = "Deletions", x = paste0("PC3 (", round(100 * summary(sv_pca)$importance[2,3], 1), "%)"),y = paste0("PC4 (", round(100 * summary(sv_pca)$importance[2,4], 1), "%)"))+
  guides(color = guide_legend(override.aes = list(size = 7))) 

# Scree plot:
scree <- data.frame(PC = seq_along(sv_pca$sdev),variance = sv_pca$sdev^2 / sum(sv_pca$sdev^2))

ggplot(scree[1:10,], aes(PC, variance * 100)) +
  geom_col(fill = "red") +
  theme(
    axis.text = element_text(size = 16, color = 'black'),
    panel.background = element_blank(),
    plot.margin = margin(10,10,10,10),
    panel.border = element_rect(color = 'black', fill = NA),
    axis.title = element_text(size = 16, color = 'black', face = 'bold')) +
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(breaks = 1:10) +
  coord_cartesian(ylim = c(0,16)) +
  labs(y = "Proportion of variance explained (%)", x = "Principal component")





### INSERTIONS ###
common_vcf <- merged_SV %>% 
  dplyr::select(-ref,-alt) %>%
  dplyr::filter(number_svs_merged >= least & number_svs_merged <= most & sv_type == "INS") %>%
  dplyr::select(-chrom, -pos, -sv_type, -sv_length, -number_svs_merged) %>%
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

geo.colors <- c("Big Island"="#66C2A5", "Molokai" = "black", "Maui" = "yellow", "Oahu" = "brown", "Kauai" = "pink", "Africa"="green", "North America" = "purple", "Europe" = "#E41A1C", "Atlantic" = "blue", 
                "Oceania" ="cyan", "unknown" = 'gray', "CGC1" = "#DB6333")

pca_df <- pca_df %>%
  mutate(label = ifelse(PC2 > 50, strain, NA))

ggplot(pca_df, aes(PC1, PC2, color = geo)) +
  geom_text_repel(aes(label = label), size = 4, max.overlaps = Inf, show.legend = FALSE) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = geo.colors) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 12, color = 'black'),
    axis.title = element_text(size = 12, color = 'black', face = 'bold'),
    legend.text = element_text(size = 14, color = 'black'),
    legend.title = element_text(size = 14, color = 'black'),
    plot.title = element_text(size = 16, color = 'black', face = 'bold', hjust = 0.5)) +
  labs(color = "Collection location", title = "Insertions", x = paste0("PC1 (", round(100 * summary(sv_pca)$importance[2,1], 1), "%)"),y = paste0("PC2 (", round(100 * summary(sv_pca)$importance[2,2], 1), "%)"))+
  guides(color = guide_legend(override.aes = list(size = 7))) 


# PC3/4
ggplot(pca_df, aes(PC3, PC4, color = geo)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = geo.colors) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 12, color = 'black'),
    axis.title = element_text(size = 12, color = 'black', face = 'bold'),
    legend.text = element_text(size = 14, color = 'black'),
    legend.title = element_text(size = 14, color = 'black'),
    plot.title = element_text(size = 16, color = 'black', face = 'bold', hjust = 0.5)) +
  labs(color = "Collection location",title = "Insertions", x = paste0("PC3 (", round(100 * summary(sv_pca)$importance[2,3], 1), "%)"),y = paste0("PC4 (", round(100 * summary(sv_pca)$importance[2,4], 1), "%)"))+
  guides(color = guide_legend(override.aes = list(size = 7))) 

# Scree plot:
scree <- data.frame(PC = seq_along(sv_pca$sdev),variance = sv_pca$sdev^2 / sum(sv_pca$sdev^2))

ggplot(scree[1:10,], aes(PC, variance * 100)) +
  geom_col(fill = "blue") +
  theme(
    axis.text = element_text(size = 16, color = 'black'),
    panel.background = element_blank(),
    plot.margin = margin(10,10,10,10),
    panel.border = element_rect(color = 'black', fill = NA),
    axis.title = element_text(size = 16, color = 'black', face = 'bold')) +
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(breaks = 1:10) +
  coord_cartesian(ylim = c(0,17)) +
  labs(y = "Proportion of variance explained (%)", x = "Principal component")



### INVERSIONS ###
common_vcf <- merged_SV %>% 
  dplyr::select(-ref,-alt) %>%
  dplyr::filter(number_svs_merged >= least & number_svs_merged <= most & sv_type == "INV") %>%
  dplyr::select(-chrom, -pos, -sv_type, -sv_length, -number_svs_merged) %>%
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

geo.colors <- c("Big Island"="#66C2A5", "Molokai" = "black", "Maui" = "yellow", "Oahu" = "brown", "Kauai" = "pink", "Africa"="green", "North America" = "purple", "Europe" = "#E41A1C", "Atlantic" = "blue", 
                "Oceania" ="cyan", "unknown" = 'gray', "CGC1" = "#DB6333")

pca_df <- pca_df %>%
  mutate(label = ifelse(PC2 > 5, strain, NA))

ggplot(pca_df, aes(PC1, PC2, color = geo)) +
  geom_text_repel(aes(label = label), size = 4, max.overlaps = Inf, show.legend = FALSE) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = geo.colors) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 12, color = 'black'),
    axis.title = element_text(size = 12, color = 'black', face = 'bold'),
    legend.text = element_text(size = 14, color = 'black'),
    legend.title = element_text(size = 14, color = 'black'),
    plot.title = element_text(size = 16, color = 'black', face = 'bold', hjust = 0.5)) +
  labs(color = "Collection location",title = "Inversions", x = paste0("PC1 (", round(100 * summary(sv_pca)$importance[2,1], 1), "%)"),y = paste0("PC2 (", round(100 * summary(sv_pca)$importance[2,2], 1), "%)")) +
  guides(color = guide_legend(override.aes = list(size = 7))) 


# PC3/4
ggplot(pca_df, aes(PC3, PC4, color = geo)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = geo.colors) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 12, color = 'black'),
    axis.title = element_text(size = 12, color = 'black', face = 'bold'),
    legend.text = element_text(size = 14, color = 'black'),
    legend.title = element_text(size = 14, color = 'black'),
    plot.title = element_text(size = 16, color = 'black', face = 'bold', hjust = 0.5)) +
  labs(color = "Collection location",title = "Inversions", x = paste0("PC3 (", round(100 * summary(sv_pca)$importance[2,3], 1), "%)"),y = paste0("PC4 (", round(100 * summary(sv_pca)$importance[2,4], 1), "%)"))+
  guides(color = guide_legend(override.aes = list(size = 7))) 

# Scree plot:
scree <- data.frame(PC = seq_along(sv_pca$sdev),variance = sv_pca$sdev^2 / sum(sv_pca$sdev^2))

ggplot(scree[1:10,], aes(PC, variance * 100)) +
  geom_col(fill = "gold") +
  theme(
    axis.text = element_text(size = 16, color = 'black'),
    panel.background = element_blank(),
    plot.margin = margin(10,10,10,10),
    panel.border = element_rect(color = 'black', fill = NA),
    axis.title = element_text(size = 16, color = 'black', face = 'bold')) +
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(breaks = 1:10) +
  coord_cartesian(ylim = c(0,16)) +
  labs(y = "Proportion of variance explained (%)", x = "Principal component")


















# 
# ####################################################################################################
# ####################################################################################################
# 
# # CIRCOS PLOT
# 
# ####################################################################################################
# ####################################################################################################
# snps <- readr::read_tsv("../../processed_data/genome_resources/genome_data/140WSs_biallelicSNPs.tsv", col_names = c("chrom","pos","ref","alt")) 
# merged_SV <- readr::read_tsv("../../processed_data/structural_variants/Jasmine_merged_SVs.tsv")
# 
# # Outer ring of chromosomes (gene map of gene models represented with black rectangles)
# ## Chromosome IDs and sizes (start is always equal to zero)
# chr_order <- c("I","II","III","IV","V","X")
# chrom_sizes <- readr::read_tsv("../../data/N2.WS283.cleaned.fa.fai", col_names = c("chrom","start","end")) %>%
#   dplyr::mutate(chrom = factor(chrom, levels = chr_order)) %>%
#   dplyr::mutate(chrom = as.character(chrom))
# 
# ## N2 gene modesl in BED format
# n2_genes_bed <- n2_genes_plt %>% dplyr::filter(chrom != "MtDNA") %>%
#   transmute(chrom = as.character(chrom), start = as.numeric(start), end = as.numeric(end))
# 
# gene_bin_size <- 50000L
# 
# # one row per chr with lengths
# chr_len_df <- chrom_sizes %>%
#   dplyr::filter(chrom %in% chr_order) %>%
#   dplyr::transmute(chrom = as.character(chrom), chr_len = as.numeric(end)) %>%
#   dplyr::distinct()
# 
# # bins that cover the entire chromosome, including the last partial bin
# gene_bins <- chr_len_df %>%
#   dplyr::group_by(chrom) %>%
#   dplyr::do({
#     L <- .$chr_len[1]
#     starts <- seq(0, L, by = gene_bin_size)  # include L so last bin is created
#     tibble(
#       chrom = .$chrom[1],
#       start = starts[-length(starts)],
#       end   = pmin(starts[-1], L)
#     )
#   }) %>%
#   dplyr::ungroup() %>%
#   dplyr::mutate(mid = (start + end)/2)
# 
# # overlap-count genes per bin
# bins_dt  <- as.data.table(gene_bins)
# genes_dt <- as.data.table(n2_genes_bed)  # chrom/start/end already numeric in your code
# 
# setkey(bins_dt,  chrom, start, end)
# setkey(genes_dt, chrom, start, end)
# 
# ov <- foverlaps(genes_dt, bins_dt, nomatch = 0)
# 
# gene_counts <- ov[, .(gene_count = .N), by = .(chrom, start, end, mid)]
# 
# bins_gene <- merge(bins_dt, gene_counts, by = c("chrom","start","end","mid"), all.x = TRUE)
# bins_gene[is.na(gene_count), gene_count := 0]
# 
# gene_bins_50kb <- as.data.frame(bins_gene) %>%
#   dplyr::rename(pos = mid, value = gene_count)
# 
# gene_ylim <- c(0, max(gene_bins_50kb$value, na.rm = TRUE))
# 
# # Next ring in will have SNV count per kb (plot fitted LOESS line?) 
# ## SNV count in table: "chrom", "bin", "variant_count"
# bin_size <- 1000L
# 
# bins <- chrom_sizes %>%
#   dplyr::group_by(chrom) %>%
#   dplyr::do({
#     chr_len <- .$end[1]
#     starts <- seq(0, chr_len - 1, by = bin_size)
#     tibble(
#       chrom = .$chrom[1],
#       start = starts,
#       end   = pmin(starts + bin_size, chr_len)
#     )
#   }) %>%
#   dplyr::ungroup()
# 
# bins_dt <- as.data.table(bins)
# setkey(bins_dt, chrom, start, end)
# bins_dt[, id := .I]
# 
# snps_counts <- snps %>%
#   dplyr::select(-ref,-alt) %>%
#   dplyr::filter(chrom %in% chr_order) %>%
#   dplyr::mutate(
#     chrom = as.character(chrom),
#     start = (pos %/% bin_size) * bin_size,
#     end   = start + bin_size) %>%
#   dplyr::count(chrom, start, end, name = "variant_count")
# 
# snps_dt <- as.data.table(snps_counts)
# setkey(snps_dt, chrom, start, end)
# 
# bins_snp <- merge(bins_dt, snps_dt, by = c("chrom","start","end"), all.x = TRUE)
# bins_snp[is.na(variant_count), variant_count := 0]
# snps_per_bin <- as.data.frame(bins_snp) %>%
#   dplyr::mutate(pos = (start + end)/2) %>%
#   dplyr::select(chrom,pos,variant_count) %>% dplyr::rename(value = variant_count)
# 
# # Next ring in will plot DEL frequncy (plotted as LOESS line) - make sure CGC1 is filtered out
# ## DEL calls with "chrom", "middle" (of the 1 kb bin), and "freq"
# # Calculating DEL frequency
# # bins <- snps %>% dplyr::select(chrom, bin)%>% dplyr::group_by(chrom) %>% dplyr::mutate(binEnd = lead(bin)) %>% dplyr::slice(-dplyr::n()) %>% dplyr::ungroup() %>% dplyr::rename(start = bin, end = binEnd)
# bins_dt <- as.data.table(bins)
# bins_dt[, id := .I]  # optional: keep track of bins
# 
# del_calls <- filt_calls %>% dplyr::filter(sv_type == "DEL", strain != "CGC1") %>% dplyr::mutate(end = pos + sv_length) %>% dplyr::rename(start = pos) %>% dplyr::select(chrom,start,end,strain)
# del_calls_dt <- as.data.table(del_calls)
# 
# setkey(bins_dt, chrom, start, end)
# setkey(del_calls_dt, chrom, start, end)
# 
# overlaps <- data.table::foverlaps(del_calls_dt, bins_dt, nomatch = 0)
# 
# # Count unique strains per bin
# counts_PB <- overlaps[, .(n_strains = uniqueN(strain)), by = .(chrom, start, end)]
# 
# # If you want to merge with the full bin list (including 0s):
# bins_wCounts <- merge(bins_dt, counts_PB, by = c("chrom", "start", "end"), all.x = TRUE)
# bins_wCounts[is.na(n_strains), n_strains := 0]
# 
# bins_wFreq <- as.data.frame(bins_wCounts) %>%
#   dplyr::mutate(freq = n_strains/141) #change me to number of isotypes
# 
# del_bin_plt <- bins_wFreq %>% dplyr::mutate(middle = (end + start) / 2)
# 
# del_bin_freq <- del_bin_plt %>% dplyr::select(chrom, middle, freq)
# 
# # Next ring in will plot INS frequency (plotted as LOESS line) - make sure CGC1 is filtered out
# ## INS calls with "chrom", "middle" (of the 1 kb bin), and "freq"
# bins_dt <- as.data.table(bins)
# bins_dt[, id := .I]  # optional: keep track of bins
# 
# ins_calls <- filt_calls %>% dplyr::filter(sv_type == "INS", strain != "CGC1") %>% dplyr::mutate(end = pos + sv_length) %>% dplyr::rename(start = pos) %>% dplyr::select(chrom,start,end,strain)
# ins_calls_dt <- as.data.table(ins_calls)
# 
# setkey(bins_dt, chrom, start, end)
# setkey(ins_calls_dt, chrom, start, end)
# 
# overlaps <- data.table::foverlaps(ins_calls_dt, bins_dt, nomatch = 0)
# 
# # Count unique strains per bin
# counts_PB <- overlaps[, .(n_strains = uniqueN(strain)), by = .(chrom, start, end)]
# 
# # If you want to merge with the full bin list (including 0s):
# bins_wCounts <- merge(bins_dt, counts_PB, by = c("chrom", "start", "end"), all.x = TRUE)
# bins_wCounts[is.na(n_strains), n_strains := 0]
# 
# bins_wFreq <- as.data.frame(bins_wCounts) %>%
#   dplyr::mutate(freq = n_strains/141) #change me to number of isotypes
# 
# ins_bin_plt <- bins_wFreq %>% dplyr::mutate(middle = (end + start) / 2)
# 
# ins_bin_freq <- ins_bin_plt %>% dplyr::select(chrom, middle, freq)
# 
# # Then the final, more inner ring will have INV frequency (plotted as LOESS line) - make sure CGC1 is filtered out
# ## INV calls with "chrom", "middle" (of the 1 kb bin), and "freq"
# bins_dt <- as.data.table(bins)
# bins_dt[, id := .I]  # optional: keep track of bins
# 
# inv_calls <- filt_calls %>% dplyr::filter(sv_type == "INV", strain != "CGC1") %>% dplyr::mutate(end = pos + sv_length) %>% dplyr::rename(start = pos) %>% dplyr::select(chrom,start,end,strain)
# inv_calls_dt <- as.data.table(inv_calls)
# 
# setkey(bins_dt, chrom, start, end)
# setkey(inv_calls_dt, chrom, start, end)
# 
# overlaps <- data.table::foverlaps(inv_calls_dt, bins_dt, nomatch = 0)
# 
# # Count unique strains per bin
# counts_PB <- overlaps[, .(n_strains = uniqueN(strain)), by = .(chrom, start, end)]
# 
# # If you want to merge with the full bin list (including 0s):
# bins_wCounts <- merge(bins_dt, counts_PB, by = c("chrom", "start", "end"), all.x = TRUE)
# bins_wCounts[is.na(n_strains), n_strains := 0]
# 
# bins_wFreq <- as.data.frame(bins_wCounts) %>%
#   dplyr::mutate(freq = n_strains/141) #change me to number of isotypes
# 
# inv_bin_plt <- bins_wFreq %>% dplyr::mutate(middle = (end + start) / 2)
# 
# inv_bin_freq <- inv_bin_plt %>% dplyr::select(chrom, middle, freq)
# 
# 
# ## =========================
# ##  Plotting!
# ## ========================= 
# # SV frequency tracks must have: chrom, pos, value
# del_bin_freq <- del_bin_freq %>%
#   transmute(chrom = as.character(chrom), pos = as.numeric(middle), value = as.numeric(freq))
# 
# ins_bin_freq <- ins_bin_freq %>%
#   transmute(chrom = as.character(chrom), pos = as.numeric(middle), value = as.numeric(freq))
# 
# inv_bin_freq <- inv_bin_freq %>%
#   transmute(chrom = as.character(chrom), pos = as.numeric(middle), value = as.numeric(freq))
# 
# add_rect_track <- function(bed, col, track_height = 0.08, label = NULL) {
#   circos.trackPlotRegion(
#     ylim = c(0, 1),
#     track.height = track_height,
#     bg.border = NA,
#     panel.fun = function(x, y) {
#       chr <- CELL_META$sector.index
#       d <- bed[bed$chrom == chr, , drop = FALSE]
#       if (nrow(d) == 0) return()
#       circos.rect(d$start, 0, d$end, 1, col = col, border = NA)
#     }
#   )
#   
#   if (!is.null(label)) {
#     circos.text(
#       x = 0, y = 0.5, labels = label,
#       sector.index = chr_order[1],
#       track.index  = get.current.track.index(),
#       facing = "inside", adj = c(1, 0.5), cex = 0.7
#     )
#   }
# }
# 
# add_value_track <- function(df, col, ylim, track_height = 0.10, label = NULL,
#                             type = c("line","points"), add_loess = FALSE, span = 0.2) {
#   type <- match.arg(type)
#   
#   circos.trackPlotRegion(
#     ylim = ylim,
#     track.height = track_height,
#     bg.border = NA,
#     panel.fun = function(x, y) {
#       chr <- CELL_META$sector.index
#       d <- df[df$chrom == chr, , drop = FALSE]
#       if (nrow(d) == 0) return()
#       
#       d <- d[order(d$pos), , drop = FALSE]
#       
#       if (type == "points") {
#         circos.points(d$pos, d$value, pch = 16, cex = 0.25, col = col)
#       } else {
#         circos.lines(d$pos, d$value, col = col, lwd = 1)
#       }
#       
#       if (add_loess && nrow(d) >= 50) {
#         fit <- stats::loess(value ~ pos, data = d, span = span)
#         xs  <- d$pos
#         ys  <- stats::predict(fit, newdata = data.frame(pos = xs))
#         ok  <- is.finite(ys)
#         if (any(ok)) circos.lines(xs[ok], ys[ok], col = col, lwd = 2)
#       }
#     }
#   )
#   
#   if (!is.null(label)) {
#     circos.text(
#       x = 0, y = mean(ylim), labels = label,
#       sector.index = chr_order[1],
#       track.index  = get.current.track.index(),
#       facing = "inside", adj = c(1, 0.5), cex = 0.7
#     )
#   }
# }
# 
# add_filled_area_track <- function(df, col_fill = "#00BFC480", col_line = "#00BFC4",
#                                   ylim, track_height = 0.10, label = NULL) {
#   circos.trackPlotRegion(
#     ylim = ylim,
#     track.height = track_height,
#     bg.border = NA,
#     panel.fun = function(x, y) {
#       chr <- CELL_META$sector.index
#       d <- df[df$chrom == chr, , drop = FALSE]
#       if (nrow(d) < 2) return()
#       d <- d[order(d$pos), , drop = FALSE]
#       
#       # build polygon (baseline at ylim[1])
#       x_poly <- c(d$pos, rev(d$pos))
#       y_poly <- c(d$value, rep(ylim[1], nrow(d)))
#       
#       circos.polygon(x_poly, y_poly, col = col_fill, border = NA)
#       circos.lines(d$pos, d$value, col = col_line, lwd = 1.2)
#     }
#   )
#   
#   if (!is.null(label)) {
#     circos.text(
#       x = 0, y = mean(ylim), labels = label,
#       sector.index = chr_order[1],
#       track.index  = get.current.track.index(),
#       facing = "inside", adj = c(1, 0.5), cex = 0.7
#     )
#   }
# }
# 
# 
# png("circos_plot.png", width = 3.75, height = 3.75, units = "in", res = 600, bg = "white")
# 
# circos.clear()
# circos.par(track.margin = c(0.002, 0.002), cell.padding = c(0, 0, 0, 0), canvas.xlim = c(-0.9, 0.9), canvas.ylim = c(-0.99,0.99)) # start.degree = 86, gap.after = c(rep(2, length(chr_order)-1), 8)
# circos.initialize(factors = as.character(chrom_sizes$chrom), xlim = cbind(rep(0, nrow(chrom_sizes)), chrom_sizes$end))
# 
# # Outer ideogram-like ring + gene models inside it
# circos.trackPlotRegion(
#   ylim = c(0, 1),
#   track.height = 0.02,
#   bg.border = NA,
#   panel.fun = function(x, y) {
#     chr <- CELL_META$sector.index
#     xlim <- CELL_META$xlim
#     
#     # background chromosome band
#     circos.rect(xlim[1], 0, xlim[2], 1, col = "grey90", border = "black")
#     # gene models as black rectangles (thin band)
#     # g <- n2_genes_bed[n2_genes_bed$chrom == chr, , drop = FALSE]
#     # if (nrow(g) > 0) {
#     #   circos.rect(g$start, 0, g$end, 1, col = "black", border = NA)}
#     
#     # gene density (50 kb bins) as cyan-blue line
#     # gd <- gene_density_50kb[gene_density_50kb$chrom == chr, , drop = FALSE]
#     # if (nrow(gd) > 1) {
#     #   gd <- gd[order(gd$pos), , drop = FALSE]
#     #   
#     #   # rescale density to fit nicely within the outer track band
#     #   y <- (gd$value - gene_ylim[1]) / (gene_ylim[2] - gene_ylim[1] + 1e-9)
#     #   y <- 0.10 + y * 0.80   # occupy 10%..90% of the track height
#     #   
#     #   circos.lines(gd$pos, y, col = "#00BFC4", lwd = 2)
#     # }
#     
#     
#     
#     # chromosome label
#     circos.text(CELL_META$xcenter, 3.3, chr, facing = "bending.outside", niceFacing = T, cex = 1.5, font = 2)
#     # axis ticks
#     circos.axis(h = "top", major.at = seq(0, xlim[2], by = 5e6),labels = seq(0, xlim[2], by = 5e6) / 1e6, labels.cex = 1, major.tick.length = 0.001)
#   }
# )
# 
# 
# snp_ylim <- c(0, max(snps_per_bin$value, na.rm = TRUE))
# add_value_track(snps_per_bin, col = "#DB6333", ylim = snp_ylim, track_height = 0.2, type = "line", add_loess = FALSE, span = 0.2)
# add_value_track(del_bin_freq, col = "red",  ylim = c(0, 1), track_height = 0.2, type = "line", add_loess = FALSE, span = 0.2)
# add_value_track(ins_bin_freq, col = "blue", ylim = c(0, 1), track_height = 0.2, type = "line", add_loess = FALSE, span = 0.2)
# add_value_track(inv_bin_freq, col = "gold3", ylim = c(0, 1), track_height = 0.2, type = "line", add_loess = FALSE, span = 0.2)
# 
# # Legend
# lgd <- Legend(labels = c("SNPs/kb","DEL freq.","INS freq.","INV freq."),
#               type = "lines", legend_gp = gpar(col = c("#DB6333","red","blue","gold3"),
#                                                lwd = c(3, 3, 3, 3)),
#               labels_gp = gpar(fontsize = 10),
#               title_gp  = gpar(fontsize = 10))
# 
# draw(lgd, x = unit(0.99, "npc"), y = unit(0.99, "npc"), just = c("right", "top"))
# 
# circos.clear()
# 
# dev.off()
# 
# CIRC <- ggdraw() + draw_image("circos_plot.png")
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
