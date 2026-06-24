library(readr)
library(ggplot2)
library(dplyr)
library(data.table)

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
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text.y = element_text(size = 10, color = 'black'),
    strip.text = element_blank(),
    panel.border = element_rect(color = 'black', fill = NA),
    panel.background = element_blank()) +
  coord_cartesian(ylim = c(0,1.01)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 1), breaks = seq(0, 0.75, 0.25))

# Concatenate together all three SV types!
all_three_SVs <- cowplot::plot_grid(
  del_freq, ins_freq, inv_freq,
  nrow = 3,
  align = "v",
  rel_heights = c(1,1,1))
all_three_SVs

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








