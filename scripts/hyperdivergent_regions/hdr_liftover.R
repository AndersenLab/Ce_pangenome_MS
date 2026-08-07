library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(data.table)
library(cowplot)

# ======================================================================================================================================================================================== #
# Lifting over N2 HDRs to wild strains #
# ======================================================================================================================================================================================== #
# Loading in HDRs and wild strain list....
hdrs <- readr::read_tsv("../../data/20250625_c_elegans_divergent_regions_strain.bed", col_names = c("chrom", "start", "end", "strain"))
WSs <- readr::read_tsv("../../tables/strain_SVs_cum_length.tsv") %>% dplyr::distinct(strain) %>% dplyr::filter(strain != "CGC1") %>% dplyr::pull(strain)

# Filtering for HDRs in pangenome strain set
strain_hdr <- hdrs %>%
  dplyr::filter(strain %in% WSs) %>%
  dplyr::arrange(chrom,start) %>%
  dplyr::rename(og_hdr_start = start, og_hdr_end = end) %>% 
  dplyr::mutate(start = ifelse(og_hdr_start >= 5000, og_hdr_start - 0, og_hdr_start), end = og_hdr_end + 0) %>%
  data.table::as.data.table()

# Load alignments
nucmer <- readr::read_tsv("../../processed_data/genome_resources/genome_data/141_nucmer_ECA741CGC1.tsv", col_names = c("N2S","N2E","WSS","WSE","L1","L2","IDY","LENR","LENQ","N2_chr","contig","strain")) %>% 
  dplyr::filter(strain != "CGC1") %>% # removing alignments for CGC1, there are no HDRs for CGC1 because it is N2
  dplyr::mutate(inv = ifelse((WSS > WSE), T, F)) 

# Rename variables, arrange by WS coordinates, get lag and leading coordinates of alignments in WS space
nucmer_ranges <- nucmer %>%
  dplyr::rename(start = N2S, end = N2E, chrom = N2_chr) %>%
  dplyr::select(chrom, start, end, L1, WSS, WSE, contig, L2, LENQ, inv, strain) %>%
  dplyr::group_by(strain,chrom) %>%
  dplyr::arrange(strain,contig,WSS) %>%
  dplyr::mutate(leadS=lead(start),leadE=lead(end),leadWSS=lead(WSS),leadWSE=lead(WSE),lagS=lag(start),lagE=lag(end),lagWSS=lag(WSS),lagWSE=lag(WSE)) %>%
  dplyr::ungroup() %>%
  data.table::as.data.table()

# Setting keys for finding overlaps of HDRs with WS contig alignments
data.table::setkey(nucmer_ranges, strain, chrom, start, end)
data.table::setkey(strain_hdr, strain, chrom, start, end)

# Find overlaps
hdr_aln <- data.table::foverlaps(
  x = strain_hdr,
  y = nucmer_ranges,
  type = "any",
  nomatch = NA) %>%
  dplyr::rename(hdr_start_extended = i.start, hdr_end_extended = i.end, N2S = start, N2E = end)  %>%
  dplyr::mutate(HDRid=paste0(strain, chrom, og_hdr_start,og_hdr_end))

# Select longest contig among alignments
nucmer_longest <- hdr_aln %>% # equivalent of tigFilt from haplotypePlotter.R
  dplyr::group_by(HDRid) %>%
  dplyr::mutate(nalign = n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(HDRid,contig) %>%
  dplyr::mutate(ntig = n()) %>%
  dplyr::mutate(alnsize = sum(L2)) %>% # summing the number of alignments that overlap with an N2 gene... some contigs align to a single gene many times
  dplyr::ungroup() %>% 
  dplyr::group_by(HDRid) %>%
  dplyr::filter(alnsize == max(alnsize)) %>%
  dplyr::rename(longest_contig = contig) %>%
  dplyr::filter(LENQ == max(LENQ)) %>% # to filter out alignments that are the same size, but from different contigs
  dplyr::ungroup()

# Remove jumps by HDR alignment group
nucmer_longest_jumpRemoved <- nucmer_longest %>%
  dplyr::mutate(St2=ifelse(inv==T,WSE,WSS), Et2=ifelse(inv==T, WSS, WSE)) %>%
  dplyr::arrange(HDRid,St2) %>%
  dplyr::group_by(HDRid) %>%
  dplyr::mutate(leadDiff=lead(St2)-Et2) %>%
  dplyr::mutate(jump=ifelse(leadDiff > 7.5e4, 1 ,0)) %>% # Modified to 75kb - reduced situations where HDRs are called for repetitive alignment with large, interspersed jumps
  dplyr::mutate(leadDiff=ifelse(is.na(leadDiff),0,leadDiff)) %>%
  dplyr::mutate(run_id = cumsum(c(1, head(jump, -1)))) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(HDRid,run_id) %>%
  dplyr::mutate(gsize=n()) %>%
  dplyr::mutate(len=abs(Et2-St2)) %>%
  dplyr::mutate(sumlen=sum(len)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(HDRid) %>%
  dplyr::filter(sumlen == max(sumlen)) %>%
  dplyr::mutate(mark_spurious = ifelse(length(unique(run_id)) > 1, TRUE, FALSE)) %>%
  dplyr::filter(!mark_spurious | dplyr::row_number() == 1) %>%
  dplyr::select(-gsize) %>%
  dplyr::ungroup() 

# Trim boundary alignments
trim_spacer = 1e3
# Trims long alignments to the focal region (i.e. hap_start to hap_end, but transformed to the other genome)
nucmer_longest_trimmed <- nucmer_longest_jumpRemoved %>%
  dplyr::rowwise() %>%
  dplyr::mutate(scale_distortion = ((L2 - L1)/L1)) %>%
  dplyr::mutate(lboundDist=hdr_start_extended-N2S) %>% 
  dplyr::mutate(rboundDist=N2E-hdr_end_extended) %>%
  dplyr::mutate(N2E=ifelse(rboundDist>trim_spacer,(N2E-(rboundDist-trim_spacer)),N2E)) %>%
  dplyr::mutate(N2S=ifelse(lboundDist>trim_spacer,(N2S+(lboundDist-trim_spacer)),N2S)) %>%
  dplyr::mutate(WSE=ifelse(rboundDist>trim_spacer & inv==T,(WSE+(rboundDist-trim_spacer)+(rboundDist*scale_distortion)),WSE)) %>%
  dplyr::mutate(WSS=ifelse(lboundDist>trim_spacer & inv==T,(WSS-(lboundDist-trim_spacer)-(rboundDist*scale_distortion)),WSS)) %>%
  dplyr::mutate(WSE=ifelse(rboundDist>trim_spacer & inv==F,(WSE-(rboundDist-trim_spacer)-(rboundDist*scale_distortion)),WSE)) %>%
  dplyr::mutate(WSS=ifelse(lboundDist>trim_spacer & inv==F,(WSS+(lboundDist-trim_spacer)+(rboundDist*scale_distortion)),WSS)) %>%
  dplyr::ungroup()

# Mark potential extensions
nucmer_mark_extend <- nucmer_longest_trimmed %>%
  dplyr::group_by(HDRid) %>%
  dplyr::mutate(WSE_extend=ifelse(inv==F & WSE == max(WSE) & N2E < hdr_end_extended, T, F),
                WSS_extend=ifelse(inv==F & WSS == min(WSS) & N2S > hdr_start_extended, T,F),
                iWSE_extend=ifelse(inv==T & WSE == min(WSE) & N2E < hdr_end_extended,T,F),
                iWSS_extend=ifelse(inv==T & WSS == max(WSS) & N2S > hdr_start_extended, T,F)) %>%
  dplyr::mutate(any_extend=ifelse(WSE_extend == T | WSS_extend == T | iWSE_extend==T | WSS_extend ==T,T,F)) %>%
  dplyr::ungroup()

# Estimate extension distances
nucmer_to_extend <- nucmer_mark_extend %>% 
  dplyr::filter(any_extend==T) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(extend_length_WI_lead=ifelse(WSE_extend==T & min(leadS,leadE) > hdr_end_extended & WSE < min(leadWSS,leadWSE), min(leadWSS,leadWSE)-WSE,
                                             ifelse(iWSS_extend==T & max(leadS,leadE) < hdr_start_extended & WSS < min(leadWSS,leadWSE), min(leadWSS,leadWSE)-WSS,NA)),
                extend_length_WI_lag=ifelse(WSS_extend==T & max(lagS,lagE) < hdr_start_extended & WSS > max(lagWSS,lagWSE), WSS-max(lagWSS,lagWSE),
                                            ifelse(iWSE_extend==T & min(lagS,lagE) > hdr_end_extended & WSE > max(lagWSS,lagWSE), WSE-max(lagWSS,lagWSE),NA)),
                extend_length_REF_lead=ifelse(WSE_extend==T & min(leadS,leadE) > hdr_end_extended & WSE < min(leadWSS,leadWSE), leadS-N2E,
                                              ifelse(iWSS_extend==T & max(leadS,leadE) < hdr_start_extended & WSS < min(leadWSS,leadWSE), N2S-leadE,NA)),
                extend_length_REF_lag=ifelse(WSS_extend==T & max(lagS,lagE) < hdr_start_extended & WSS > max(lagWSS,lagWSE), N2S-lagE,
                                             ifelse(iWSE_extend==T & min(lagS,lagE) > hdr_end_extended & WSE > max(lagWSS,lagWSE), lagS-N2E,NA))) %>%
  dplyr::ungroup()

# Bind unmarked extensions with aligments-to-extend
nucmer_extensions <- rbind(nucmer_to_extend,nucmer_mark_extend %>% dplyr::filter(any_extend==F) %>% dplyr::mutate(extend_length_WI_lead=NA,extend_length_REF_lead=NA,extend_length_WI_lag=NA,extend_length_REF_lag=NA))

# Extend alignments when conditions are met
ext_max= 5e4
nucmer_extended <- nucmer_extensions %>% 
  dplyr::rowwise() %>%
  dplyr::mutate(N2E=ifelse(WSE_extend==T & min(leadS,leadE) > hdr_end_extended & WSE < min(leadWSS,leadWSE) & !is.na(extend_length_WI_lead) & !is.na(extend_length_REF_lead) & extend_length_WI_lead < ext_max & extend_length_REF_lead < ext_max, leadS,N2E)) %>%
  dplyr::mutate(N2S=ifelse(iWSS_extend==T & max(leadS,leadE) < hdr_start_extended & WSS < min(leadWSS,leadWSE) & !is.na(extend_length_WI_lag) & !is.na(extend_length_REF_lag) & extend_length_WI_lag < ext_max & extend_length_REF_lag < ext_max,leadE ,N2S)) %>%
  dplyr::mutate(N2S=ifelse(WSS_extend==T & max(lagS,lagE) < hdr_start_extended & WSS > max(lagWSS,lagWSE) & !is.na(extend_length_WI_lag) & !is.na(extend_length_REF_lag) & extend_length_WI_lag < ext_max & extend_length_REF_lag < ext_max, lagE,N2S)) %>%
  dplyr::mutate(N2E=ifelse(iWSE_extend==T & min(lagS,lagE) > hdr_end_extended & WSE > max(lagWSS,lagWSE) & !is.na(extend_length_WI_lead) & !is.na(extend_length_REF_lead) & extend_length_WI_lead < ext_max & extend_length_REF_lead < ext_max, lagS,N2E)) %>%
  dplyr::mutate(WSE=ifelse(WSE_extend==T & min(leadS,leadE) > hdr_end_extended & WSE < min(leadWSS,leadWSE) & !is.na(extend_length_WI_lead) & !is.na(extend_length_REF_lead) & extend_length_WI_lead < ext_max & extend_length_REF_lead < ext_max, min(leadWSS,leadWSE),WSE)) %>%
  dplyr::mutate(WSS=ifelse(iWSS_extend==T & max(leadS,leadE) < hdr_start_extended & WSS < min(leadWSS,leadWSE) & !is.na(extend_length_WI_lag) & !is.na(extend_length_REF_lag) & extend_length_WI_lag < ext_max & extend_length_REF_lag < ext_max,min(leadWSS,leadWSE) ,WSS)) %>%
  dplyr::mutate(WSS=ifelse(WSS_extend==T & max(lagS,lagE) < hdr_start_extended & WSS > max(lagWSS,lagWSE) & !is.na(extend_length_WI_lag) & !is.na(extend_length_REF_lag) & extend_length_WI_lag < ext_max & extend_length_REF_lag < ext_max, max(lagWSS,lagWSE),WSS)) %>%
  dplyr::mutate(WSE=ifelse(iWSE_extend==T & min(lagS,lagE) > hdr_end_extended & WSE > max(lagWSS,lagWSE) & !is.na(extend_length_WI_lead) & !is.na(extend_length_REF_lead) & extend_length_WI_lead < ext_max & extend_length_REF_lead < ext_max, max(lagWSS,lagWSE),WSE)) %>% 
  dplyr::ungroup()

# Lift over HDRs from N2 to wild strains
WS_HDRs <- nucmer_extended %>% 
  dplyr::select(strain, longest_contig, WSS, WSE, any_extend,
                HDRid, chrom, hdr_start_extended, hdr_end_extended) %>%
  dplyr::mutate(minStart = pmin(WSS, WSE, na.rm = TRUE),
                maxEnd   = pmax(WSS, WSE, na.rm = TRUE)) %>%
  dplyr::group_by(HDRid) %>%
  dplyr::summarise(strain             = dplyr::first(strain),
                   longest_contig     = dplyr::first(longest_contig),
                   minStart           = min(minStart, na.rm = TRUE),
                   maxEnd             = max(maxEnd,   na.rm = TRUE),
                   any_extend         = dplyr::first(any_extend),
                   chrom              = dplyr::first(chrom),
                   hdr_start_extended = dplyr::first(hdr_start_extended),
                   hdr_end_extended   = dplyr::first(hdr_end_extended),
                   .groups = "drop") %>%
  dplyr::mutate(divSize=maxEnd-minStart,og_divSize=hdr_end_extended-hdr_start_extended,sizeDiff=abs(og_divSize-divSize))


##### DIAGNOSTIC PLOTTING FUNCTION ######
plot_hdr_workflow <- function(HDOI) {
  
  # visualize selected contig
  ctg_select_plt <- ggplot(data = nucmer_longest %>% dplyr::filter(HDRid == HDOI[2])) +
    geom_rect(aes(xmin = hdr_start_extended / 1e6, xmax = hdr_end_extended / 1e6, ymin = -Inf, ymax = Inf), fill = "gray", alpha = 0.1) +
    geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6, color = longest_contig), linewidth = 1) +
    facet_wrap(~chrom) +
    theme(panel.border = element_rect(color = 'black', fill = NA),
          panel.background = element_blank()) +
    labs(x = "N2 genome position (Mb)", y = "WS contig position (Mb)", title = paste0(HDOI[1]))
  
  # visualize pre and post jump removal
  ctg_jumprm_plt <- ggplot(data = nucmer_longest_jumpRemoved %>% dplyr::filter(HDRid == HDOI[2])) +
    geom_rect(aes(xmin = hdr_start_extended / 1e6, xmax = hdr_end_extended / 1e6, ymin = -Inf, ymax = Inf), fill = "gray", alpha = 0.3) +
    geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6, color = longest_contig), linewidth = 1) +
    facet_wrap(~chrom) +
    theme(panel.border = element_rect(color = 'black', fill = NA),
          panel.background = element_blank()) +
    labs(x = "N2 genome position (Mb)", y = "WS contig position (Mb)", title = paste0(HDOI[1], " - POST JUMP RM"))
  
  # visualize trimming after jump removal
  ctg_trim_plt <- ggplot(data = nucmer_longest_trimmed %>% dplyr::filter(HDRid == HDOI[2])) +
    geom_rect(aes(xmin = hdr_start_extended / 1e6, xmax = hdr_end_extended / 1e6, ymin = -Inf, ymax = Inf), fill = "gray", alpha = 0.3) +
    geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6, color = longest_contig), linewidth = 1) +
    facet_wrap(~chrom) +
    theme(panel.border = element_rect(color = 'black', fill = NA),
          panel.background = element_blank()) +
    labs(x = "N2 genome position (Mb)", y = "WS contig position (Mb)", title = paste0(HDOI[1], " - POST TRIM"))
  
  # visualize mark for extending
  ctg_mark_plt <- ggplot(data = nucmer_mark_extend %>% dplyr::filter(HDRid == HDOI[2])) +
    geom_rect(aes(xmin = hdr_start_extended / 1e6, xmax = hdr_end_extended / 1e6, ymin = -Inf, ymax = Inf), fill = "gray", alpha = 0.3) +
    geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6, color = any_extend), linewidth = 1) +
    facet_wrap(~chrom) +
    theme(panel.border = element_rect(color = 'black', fill = NA),
          panel.background = element_blank()) +
    labs(x = "N2 genome position (Mb)", y = "WS contig position (Mb)", title = paste0(HDOI[1], " - POST MARK"))
  
  # visualize potential extensions
  hdr_leadlag_plt <- ggplot(data = nucmer_mark_extend %>% dplyr::filter(HDRid == HDOI[2])) +
    geom_rect(aes(xmin = hdr_start_extended / 1e6, xmax = hdr_end_extended / 1e6, ymin = -Inf, ymax = Inf), fill = "gray", alpha = 0.3) +
    geom_segment(data = nucmer_mark_extend %>%
                   dplyr::filter(any_extend == T) %>%
                   dplyr::filter(HDRid == HDOI[2]) %>%
                   dplyr::filter(N2E < hdr_end_extended),
                 aes(x = leadS / 1e6, xend = leadE / 1e6, y = leadWSS / 1e6, yend = leadWSE / 1e6, color = "leading"), linewidth = 1) +
    geom_segment(data = nucmer_mark_extend %>%
                   dplyr::filter(any_extend == T) %>%
                   dplyr::filter(HDRid == HDOI[2]) %>%
                   dplyr::filter(N2S > hdr_start_extended),
                 aes(x = lagS / 1e6, xend = lagE / 1e6, y = lagWSS / 1e6, yend = lagWSE / 1e6, color = "lagging"), linewidth = 1) +
    geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6, color = "overlapping"), linewidth = 1) +
    facet_wrap(~chrom) +
    theme(panel.border = element_rect(color = 'black', fill = NA),
          panel.background = element_blank(),
          legend.title = element_blank()) +
    labs(x = "N2 genome position (Mb)", y = "WS contig position (Mb)", title = paste0(HDOI[1], " - POTENTIAL EXTENSION"))
  
  # visualize extensions
  hdr_extended_plt <- ggplot(data = nucmer_extended %>% dplyr::filter(HDRid == HDOI[2])) +
    geom_rect(aes(xmin = hdr_start_extended / 1e6, xmax = hdr_end_extended / 1e6, ymin = -Inf, ymax = Inf), fill = "gray", alpha = 0.3) +
    geom_segment(data = nucmer_mark_extend %>%
                   dplyr::filter(any_extend == T) %>%
                   dplyr::filter(HDRid == HDOI[2]) %>%
                   dplyr::filter(N2E < hdr_end_extended),
                 aes(x = leadS / 1e6, xend = leadE / 1e6, y = leadWSS / 1e6, yend = leadWSE / 1e6, color = "leading"), linewidth = 1) +
    geom_segment(data = nucmer_mark_extend %>%
                   dplyr::filter(any_extend == T) %>%
                   dplyr::filter(HDRid == HDOI[2]) %>%
                   dplyr::filter(N2S > hdr_start_extended),
                 aes(x = lagS / 1e6, xend = lagE / 1e6, y = lagWSS / 1e6, yend = lagWSE / 1e6, color = "lagging"), linewidth = 1) +
    geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6, color = "w/Extension"), linewidth = 1) +
    geom_segment(data = nucmer_mark_extend %>%
                   dplyr::filter(HDRid == HDOI[2]),
                 aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6, color = "overlapping"), linewidth = 1) +
    facet_wrap(~chrom) +
    theme(panel.border = element_rect(color = 'black', fill = NA),
          panel.background = element_blank(),
          legend.title = element_blank()) +
    labs(x = "N2 genome position (Mb)", y = "WS contig position (Mb)", title = paste0(HDOI[1], " - POST EXTENSION"))
  
  
  nS <- (WS_HDRs %>% dplyr::filter(HDRid == HDOI[2]))$minStart
  nE <- (WS_HDRs %>% dplyr::filter(HDRid == HDOI[2]))$maxEnd
  
  hdr_transformed_plt <- ggplot() +
    geom_rect(data = nucmer_extended %>% dplyr::filter(HDRid == HDOI[2]), mapping = aes(xmin = hdr_start_extended / 1e6, xmax = hdr_end_extended / 1e6, ymin = -Inf, ymax = Inf), fill = "gray", alpha = 0.3) +
    annotate("rect", ymin = nS / 1e6, ymax = nE   / 1e6, xmin = -Inf, xmax = Inf, fill = "gray", alpha = 1) +
    geom_segment(data = nucmer_mark_extend %>%
                   dplyr::filter(any_extend == T) %>%
                   dplyr::filter(HDRid == HDOI[2]) %>%
                   dplyr::filter(N2E < hdr_end_extended),
                 aes(x = leadS / 1e6, xend = leadE / 1e6, y = leadWSS / 1e6, yend = leadWSE / 1e6, color = "leading"), linewidth = 1) +
    geom_segment(data = nucmer_mark_extend %>%
                   dplyr::filter(any_extend == T) %>%
                   dplyr::filter(HDRid == HDOI[2]) %>%
                   dplyr::filter(N2S > hdr_start_extended),
                 aes(x = lagS / 1e6, xend = lagE / 1e6, y = lagWSS / 1e6, yend = lagWSE / 1e6, color = "lagging"), linewidth = 1) +
    geom_segment(data = nucmer_extended %>%
                   dplyr::filter(HDRid == HDOI[2]),
                 aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6, color = "w/Extension"), linewidth = 1) +
    geom_segment(data = nucmer_mark_extend %>%
                   dplyr::filter(HDRid == HDOI[2]),
                 aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6, color = "overlapping"), linewidth = 1) +
    facet_wrap(~chrom) +
    theme(panel.border = element_rect(color = 'black', fill = NA),
          panel.background = element_blank(),
          legend.title = element_blank()) +
    labs(x = "N2 genome position (Mb)", y = "WS contig position (Mb)", title = paste0(HDOI[1], " - POST EXTENSION"))
  
  # Return plot list
  list(ctg_select = ctg_select_plt,
       ctg_jumprm = ctg_jumprm_plt,
       ctg_trim   = ctg_trim_plt,
       ctg_mark   = ctg_mark_plt,
       lead_lag   = hdr_leadlag_plt,
       extended   = hdr_extended_plt,
       hdr        = hdr_transformed_plt)
}

# Hyper-divergent of interest
HDOI <- c("ECA3005","ECA3005V1842400019264000") 

# Get plot list for HDOI
diag_list <- plot_hdr_workflow(HDOI)

# jump correction errors visualized here
cowplot::plot_grid(diag_list[[1]]+theme(legend.position = 'none'),diag_list[[2]],nrow=1,align = 'h',axis = 'tb',rel_widths = c(0.8,1))
 
# trimm edge errors here
cowplot::plot_grid(diag_list[[2]]+theme(legend.position = 'none'),diag_list[[3]],nrow=1,align = 'h',axis = 'tb',rel_widths = c(0.8,1))
 
# mark extension errors here
cowplot::plot_grid(diag_list[[3]]+theme(legend.position = 'none'),diag_list[[4]]+theme(axis.title.y=element_blank()),diag_list[[5]],nrow=1,align = 'h',axis = 'tb',rel_widths = c(0.6,0.8,1))

# sprevious steps all at once - overview
cowplot::plot_grid(diag_list[[3]]+theme(legend.position = 'none'), diag_list[[4]]+theme(axis.title.y=element_blank()), diag_list[[5]], diag_list[[7]],nrow=1, align = 'h', axis = 'tb',rel_widths = c(0.6,0.8,1,1))



# Plotting examples 
# Example #1
HDOI1 <- c("ECA3005","ECA3005V1842400019264000") 

ws_hdr_start1 <- (WS_HDRs %>% dplyr::filter(HDRid == HDOI1[2]))$minStart
ws_hdr_end1 <- (WS_HDRs %>% dplyr::filter(HDRid == HDOI1[2]))$maxEnd

ex1 <- ggplot() +
  geom_rect(data = nucmer_extended %>% dplyr::filter(HDRid == HDOI1[2]) %>% dplyr::distinct(hdr_start_extended, .keep_all = T), 
            mapping = aes(xmin = hdr_start_extended / 1e6, xmax = hdr_end_extended / 1e6, ymin = -Inf, ymax = Inf), fill = "#DB6333", alpha = 0.3) +
  annotate("rect", ymin = ws_hdr_start1 / 1e6, ymax = ws_hdr_end1 / 1e6, xmin = -Inf, xmax = Inf, fill = "blue", alpha = 0.3) +
  
  geom_segment(data = nucmer_extended %>% dplyr::filter(HDRid == HDOI1[2]) %>% dplyr::distinct(hdr_start_extended, .keep_all = T), 
               aes(x = -Inf, xend = Inf, y = ws_hdr_start1 / 1e6, yend = ws_hdr_start1 / 1e6), color = 'gray18') +
  geom_segment(data = nucmer_extended %>% dplyr::filter(HDRid == HDOI1[2]) %>% dplyr::distinct(hdr_start_extended, .keep_all = T), 
               aes(x = -Inf, xend = Inf, y = ws_hdr_end1 / 1e6, yend = ws_hdr_end1 / 1e6), color = 'gray18') +
  geom_segment(data = nucmer_extended %>% dplyr::filter(HDRid == HDOI1[2]) %>% dplyr::distinct(hdr_start_extended, .keep_all = T), 
               aes(x = hdr_start_extended / 1e6, xend = hdr_start_extended / 1e6, y = -Inf, yend = Inf), color = 'gray18') +
  geom_segment(data = nucmer_extended %>% dplyr::filter(HDRid == HDOI1[2]) %>% dplyr::distinct(hdr_start_extended, .keep_all = T), 
               aes(x = hdr_end_extended / 1e6, xend = hdr_end_extended / 1e6, y = -Inf, yend = Inf), color = 'gray18') +
  
  geom_segment(data = nucmer_mark_extend %>%
                 dplyr::filter(any_extend == T) %>%
                 dplyr::filter(HDRid == HDOI1[2]) %>%
                 dplyr::filter(N2E < hdr_end_extended),
               aes(x = leadS / 1e6, xend = leadE / 1e6, y = leadWSS / 1e6, yend = leadWSE / 1e6, color = "leading"), linewidth = 0.75) +
  geom_segment(data = nucmer_mark_extend %>%
                 dplyr::filter(any_extend == T) %>%
                 dplyr::filter(HDRid == HDOI1[2]) %>%
                 dplyr::filter(N2S > hdr_start_extended),
               aes(x = lagS / 1e6, xend = lagE / 1e6, y = lagWSS / 1e6, yend = lagWSE / 1e6, color = "lagging"), linewidth = 0.75) +
  geom_segment(data = nucmer_extended %>%
                 dplyr::filter(HDRid == HDOI1[2]),
               aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6, color = "w/Extension"), linewidth = 0.75) +
  geom_segment(data = nucmer_mark_extend %>%
                 dplyr::filter(HDRid == HDOI1[2]),
               aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6, color = "overlapping"), linewidth = 0.75) +
  facet_wrap(~chrom) +
  scale_color_manual(values = c("leading" = 'orange', "overlapping" = "blue", "lagging" = "green3", "w/Extension" = "purple")) +
  theme(panel.border = element_rect(color = 'black', fill = NA),
        panel.background = element_blank(),
        axis.text = element_text(size = 10, color = 'black'),
        axis.title = element_text(size = 11, color = 'black'),
        strip.text = element_text(size = 14, color = 'black'),
        axis.title.x = element_blank(),
        legend.position = 'none') +
  labs(x = "N2 genome position (Mb)", y = "ECA3005 contig position (Mb)")
ex1


# Example #2
HDOI2 <- c("ECA2187","ECA2187V26820003414000")

ws_hdr_start2 <- (WS_HDRs %>% dplyr::filter(HDRid == HDOI2[2]))$minStart
ws_hdr_end2 <- (WS_HDRs %>% dplyr::filter(HDRid == HDOI2[2]))$maxEnd

ex2 <- ggplot() +
  geom_rect(data = nucmer_extended %>% dplyr::filter(HDRid == HDOI2[2]) %>% dplyr::distinct(hdr_start_extended, .keep_all = T), 
            mapping = aes(xmin = hdr_start_extended / 1e6, xmax = hdr_end_extended / 1e6, ymin = -Inf, ymax = Inf), fill = "#DB6333", alpha = 0.3) +
  annotate("rect", ymin = ws_hdr_start2 / 1e6, ymax = ws_hdr_end2 / 1e6, xmin = -Inf, xmax = Inf, fill = "blue", alpha = 0.3) +
  
  geom_segment(data = nucmer_extended %>% dplyr::filter(HDRid == HDOI2[2]) %>% dplyr::distinct(hdr_start_extended, .keep_all = T), 
               aes(x = -Inf, xend = Inf, y = ws_hdr_start2 / 1e6, yend = ws_hdr_start2 / 1e6), color = 'gray18') +
  geom_segment(data = nucmer_extended %>% dplyr::filter(HDRid == HDOI2[2]) %>% dplyr::distinct(hdr_start_extended, .keep_all = T), 
               aes(x = -Inf, xend = Inf, y = ws_hdr_end2 / 1e6, yend = ws_hdr_end2 / 1e6), color = 'gray18') +
  geom_segment(data = nucmer_extended %>% dplyr::filter(HDRid == HDOI2[2]) %>% dplyr::distinct(hdr_start_extended, .keep_all = T), 
               aes(x = hdr_start_extended / 1e6, xend = hdr_start_extended / 1e6, y = -Inf, yend = Inf), color = 'gray18') +
  geom_segment(data = nucmer_extended %>% dplyr::filter(HDRid == HDOI2[2]) %>% dplyr::distinct(hdr_start_extended, .keep_all = T), 
               aes(x = hdr_end_extended / 1e6, xend = hdr_end_extended / 1e6, y = -Inf, yend = Inf), color = 'gray18') +
  
  geom_segment(data = nucmer_mark_extend %>%
                 dplyr::filter(any_extend == T) %>%
                 dplyr::filter(HDRid == HDOI2[2]) %>%
                 dplyr::filter(N2E < hdr_end_extended),
               aes(x = leadS / 1e6, xend = leadE / 1e6, y = leadWSS / 1e6, yend = leadWSE / 1e6, color = "leading"), linewidth = 0.75) +
  geom_segment(data = nucmer_mark_extend %>%
                 dplyr::filter(any_extend == T) %>%
                 dplyr::filter(HDRid == HDOI2[2]) %>%
                 dplyr::filter(N2S > hdr_start_extended),
               aes(x = lagS / 1e6, xend = lagE / 1e6, y = lagWSS / 1e6, yend = lagWSE / 1e6, color = "lagging"), linewidth = 0.75) +
  geom_segment(data = nucmer_extended %>%
                 dplyr::filter(HDRid == HDOI2[2]),
               aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6, color = "w/Extension"), linewidth = 0.75) +
  geom_segment(data = nucmer_mark_extend %>%
                 dplyr::filter(HDRid == HDOI2[2]),
               aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6, color = "overlapping"), linewidth = 0.75) +
  facet_wrap(~chrom) +
  scale_color_manual(values = c("leading" = 'orange', "overlapping" = "blue", "lagging" = "green3", "w/Extension" = "purple")) +
  theme(panel.border = element_rect(color = 'black', fill = NA),
        panel.background = element_blank(),
        axis.text = element_text(size = 10, color = 'black'),
        axis.title = element_text(size = 11, color = 'black'),
        plot.margin = margin(t = 5, b = 5, l = 5, r = 5),
        strip.text = element_text(size = 14, color = 'black'),
        axis.title.x = element_blank(),
        legend.position = 'none') +
  labs(x = "N2 genome position (Mb)", y = "ECA2187 contig position (Mb)")
ex2


# Example #3
HDOI3 <- c("ECA3088", "ECA3088II33090003360000") 

ws_hdr_start3 <- (WS_HDRs %>% dplyr::filter(HDRid == HDOI3[2]))$minStart
ws_hdr_end3 <- (WS_HDRs %>% dplyr::filter(HDRid == HDOI3[2]))$maxEnd

ex3 <- ggplot() +
  geom_rect(data = nucmer_extended %>% dplyr::filter(HDRid == HDOI3[2]) %>% dplyr::distinct(hdr_start_extended, .keep_all = T), 
            mapping = aes(xmin = hdr_start_extended / 1e6, xmax = hdr_end_extended / 1e6, ymin = -Inf, ymax = Inf), fill = "#DB6333", alpha = 0.3) +
  annotate("rect", ymin = ws_hdr_start3 / 1e6, ymax = ws_hdr_end3 / 1e6, xmin = -Inf, xmax = Inf, fill = "blue", alpha = 0.3) +
  
  geom_segment(data = nucmer_extended %>% dplyr::filter(HDRid == HDOI3[2]) %>% dplyr::distinct(hdr_start_extended, .keep_all = T), 
               aes(x = -Inf, xend = Inf, y = ws_hdr_start3 / 1e6, yend = ws_hdr_start3 / 1e6), color = 'gray18') +
  geom_segment(data = nucmer_extended %>% dplyr::filter(HDRid == HDOI3[2]) %>% dplyr::distinct(hdr_start_extended, .keep_all = T), 
               aes(x = -Inf, xend = Inf, y = ws_hdr_end3 / 1e6, yend = ws_hdr_end3 / 1e6), color = 'gray18') +
  geom_segment(data = nucmer_extended %>% dplyr::filter(HDRid == HDOI3[2]) %>% dplyr::distinct(hdr_start_extended, .keep_all = T), 
               aes(x = hdr_start_extended / 1e6, xend = hdr_start_extended / 1e6, y = -Inf, yend = Inf), color = 'gray18') +
  geom_segment(data = nucmer_extended %>% dplyr::filter(HDRid == HDOI3[2]) %>% dplyr::distinct(hdr_start_extended, .keep_all = T), 
               aes(x = hdr_end_extended / 1e6, xend = hdr_end_extended / 1e6, y = -Inf, yend = Inf), color = 'gray18') +
  
  geom_segment(data = nucmer_mark_extend %>%
                 dplyr::filter(any_extend == T) %>%
                 dplyr::filter(HDRid == HDOI3[2]) %>%
                 dplyr::filter(N2E < hdr_end_extended),
               aes(x = leadS / 1e6, xend = leadE / 1e6, y = leadWSS / 1e6, yend = leadWSE / 1e6, color = "leading"), linewidth = 0.75) +
  geom_segment(data = nucmer_mark_extend %>%
                 dplyr::filter(any_extend == T) %>%
                 dplyr::filter(HDRid == HDOI3[2]) %>%
                 dplyr::filter(N2S > hdr_start_extended),
               aes(x = lagS / 1e6, xend = lagE / 1e6, y = lagWSS / 1e6, yend = lagWSE / 1e6, color = "lagging"), linewidth = 0.75) +
  geom_segment(data = nucmer_extended %>%
                 dplyr::filter(HDRid == HDOI3[2]),
               aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6, color = "w/Extension"), linewidth = 0.75) +
  geom_segment(data = nucmer_mark_extend %>%
                 dplyr::filter(HDRid == HDOI3[2]),
               aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6, color = "overlapping"), linewidth = 0.75) +
  facet_wrap(~chrom) +
  scale_color_manual(values = c("leading" = 'orange', "overlapping" = "blue", "lagging" = "green3", "w/Extension" = "purple")) +
  theme(panel.border = element_rect(color = 'black', fill = NA),
        panel.background = element_blank(),
        axis.text = element_text(size = 10, color = 'black'),
        axis.title = element_text(size = 11, color = 'black'),
        plot.margin = margin(t = 5, b = 5, l = 5, r = 5),
        strip.text = element_text(size = 14, color = 'black'),
        legend.position = 'none') +
  labs(x = "N2 genome position (Mb)", y = "ECA3088 contig position (Mb)")
ex3


# Example #4
HDOI4 <- c("ECA3088", "ECA3088X1435500014378000") # keep legend 

ws_hdr_start4 <- (WS_HDRs %>% dplyr::filter(HDRid == HDOI4[2]))$minStart
ws_hdr_end4 <- (WS_HDRs %>% dplyr::filter(HDRid == HDOI4[2]))$maxEnd

ex4 <- ggplot() +
  geom_rect(data = nucmer_extended %>% dplyr::filter(HDRid == HDOI4[2]) %>% dplyr::distinct(hdr_start_extended, .keep_all = T), 
            mapping = aes(xmin = hdr_start_extended / 1e6, xmax = hdr_end_extended / 1e6, ymin = -Inf, ymax = Inf), fill = "#DB6333", alpha = 0.3) +
  annotate("rect", ymin = ws_hdr_start4 / 1e6, ymax = ws_hdr_end4 / 1e6, xmin = -Inf, xmax = Inf, fill = "blue", alpha = 0.3) +
  
  geom_segment(data = nucmer_extended %>% dplyr::filter(HDRid == HDOI4[2]) %>% dplyr::distinct(hdr_start_extended, .keep_all = T), 
               aes(x = -Inf, xend = Inf, y = ws_hdr_start4 / 1e6, yend = ws_hdr_start4 / 1e6), color = 'gray18') +
  geom_segment(data = nucmer_extended %>% dplyr::filter(HDRid == HDOI4[2]) %>% dplyr::distinct(hdr_start_extended, .keep_all = T), 
               aes(x = -Inf, xend = Inf, y = ws_hdr_end4 / 1e6, yend = ws_hdr_end4 / 1e6), color = 'gray18') +
  geom_segment(data = nucmer_extended %>% dplyr::filter(HDRid == HDOI4[2]) %>% dplyr::distinct(hdr_start_extended, .keep_all = T), 
               aes(x = hdr_start_extended / 1e6, xend = hdr_start_extended / 1e6, y = -Inf, yend = Inf), color = 'gray18') +
  geom_segment(data = nucmer_extended %>% dplyr::filter(HDRid == HDOI4[2]) %>% dplyr::distinct(hdr_start_extended, .keep_all = T), 
               aes(x = hdr_end_extended / 1e6, xend = hdr_end_extended / 1e6, y = -Inf, yend = Inf), color = 'gray18') +
  
  geom_segment(data = nucmer_mark_extend %>%
                 dplyr::filter(any_extend == T) %>%
                 dplyr::filter(HDRid == HDOI4[2]) %>%
                 dplyr::filter(N2E < hdr_end_extended),
               aes(x = leadS / 1e6, xend = leadE / 1e6, y = leadWSS / 1e6, yend = leadWSE / 1e6, color = "leading"), linewidth = 0.75) +
  geom_segment(data = nucmer_mark_extend %>%
                 dplyr::filter(any_extend == T) %>%
                 dplyr::filter(HDRid == HDOI4[2]) %>%
                 dplyr::filter(N2S > hdr_start_extended),
               aes(x = lagS / 1e6, xend = lagE / 1e6, y = lagWSS / 1e6, yend = lagWSE / 1e6, color = "lagging"), linewidth = 0.75) +
  geom_segment(data = nucmer_extended %>%
                 dplyr::filter(HDRid == HDOI4[2]),
               aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6, color = "w/Extension"), linewidth = 0.75) +
  geom_segment(data = nucmer_mark_extend %>%
                 dplyr::filter(HDRid == HDOI4[2]),
               aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6, color = "overlapping"), linewidth = 0.75) +
  facet_wrap(~chrom) +
  scale_color_manual(values = c("leading" = 'orange', "overlapping" = "blue", "lagging" = "green3", "w/Extension" = "purple")) +
  theme(panel.border = element_rect(color = 'black', fill = NA),
        panel.background = element_blank(),
        axis.text = element_text(size = 10, color = 'black'),
        axis.title = element_text(size = 11, color = 'black'),
        strip.text = element_text(size = 14, color = 'black'),
        legend.position = "inside",
        # legend.background = element_rect(fill = "transparent", color = NA),
        legend.box.background = element_rect(fill = "white", color = 'black'),
        legend.position.inside = c(0.79,0.2),
        legend.text = element_text(size = 10, color = 'black')) +
  labs(x = "N2 genome position (Mb)", y = "ECA3005 contig position (Mb)", color = NULL)
ex4


# Concatenating examples to one final plot
hdr_liftover_examples <- cowplot::plot_grid(
  ex1, ex2, ex3, ex4,
  nrow = 2,
  labels = c("a","b","c","d"))
hdr_liftover_examples

# Save supplementary figure
# ggsave("../../figures/supplementary/hdr_liftover_examples.png", hdr_liftover_examples, width = 7.5, height = 7.5, dpi = 600 )

# Plotting the size difference of the WS HDR lift-overs for all HDRs among all strains
hdr_liftover_size <- ggplot(data = WS_HDRs) + 
  geom_point(aes(x = og_divSize / 1e6, y = divSize / 1e6, color = strain), size = 1) +
  geom_line(data = data.frame(x = c(0, 1.2)), aes(x = x, y = x), linetype = "dashed") +
  theme(
    panel.border = element_rect(color = 'black', fill = NA),
    panel.background = element_blank(),
    axis.text = element_text(size = 10, color= 'black'),
    # legend.text = element_text(size = 3, color = 'black'),
    # legend.title = element_text(size = 10, color = 'black'),
    axis.title = element_text(size = 11, color = 'black'),
    legend.position = 'none'
  ) +
  labs(y = "Wild strain HDR size (Mb)", x = "N2 HDR size (Mb)", color = "Strain") +
  scale_y_continuous(expand = c(0.005,0)) +
  scale_x_continuous(expand = c(0.005,0)) +
  guides(color = guide_legend(nrow = 47))
hdr_liftover_size

# Save supplementary figure
# ggsave("../../figures/supplementary/hdr_liftover_size_comparison.png", hdr_liftover_size, width = 7.5, height = 7.5, dpi = 600 )

# What is the distribution in size difference between N2 HDRs and their lifted-over size??
hdr_diff <- WS_HDRs %>%
  dplyr::select(strain, longest_contig, minStart, maxEnd, divSize, chrom, hdr_start_extended, hdr_end_extended, og_divSize, sizeDiff) %>%
  dplyr::mutate(size_ratio = ifelse(divSize <= og_divSize, 1 - (divSize / og_divSize), 1- (og_divSize / divSize))) %>% # 21,306 HDRs
  dplyr::filter(size_ratio <=0.05) # 2,273


ggplot(hdr_diff) +
  geom_histogram(aes(x = size_ratio), fill = 'orange', binwidth = 0.01) +
  theme_classic() +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))


# Find WS HDR clusters separated by less than 5kb
gap_clust_WS_HDRs <- WS_HDRs %>%
  dplyr::select(longest_contig,minStart,maxEnd,strain) %>%
  dplyr::rename(CHROM=longest_contig,STRAIN=strain) %>%
  dplyr::mutate(divSize=maxEnd-minStart) %>%
  dplyr::arrange(STRAIN,CHROM,minStart) %>%
  dplyr::group_by(STRAIN,CHROM) %>%
  dplyr::mutate(forGapSize=lead(minStart)-maxEnd) %>%
  dplyr::mutate(flag3g=ifelse(forGapSize<=5000,"clust","noclust")) %>%
  dplyr::mutate(dec3g=ifelse(flag3g=="clust" ,"join",
                             ifelse(flag3g=="noclust" & lag(flag3g)=="clust","join","nojoin"))) %>%
  dplyr::mutate(dec3g=ifelse(is.na(dec3g),"nojoin",dec3g)) %>%
  dplyr::ungroup()

# Join WS HDR clusters
joinClust_WS_HDRs<- gap_clust_WS_HDRs %>% 
  dplyr::filter(dec3g=="join") %>%
  dplyr::group_by(STRAIN,CHROM) %>%
  dplyr::mutate(segbreak=ifelse(flag3g=="noclust",paste0(dec3g,row_number()),NA)) %>%
  tidyr::fill(segbreak,.direction = 'up') %>%
  dplyr::mutate(gid=data.table::rleid(segbreak)) %>%
  dplyr::ungroup() %>%
  dplyr::rowwise() %>%
  dplyr::mutate(conID=paste0(CHROM,"-",STRAIN,"-",gid)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(conID) %>%
  dplyr::mutate(newStart=min(minStart),newEnd=max(maxEnd)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(conID) %>%
  dplyr::mutate(newDivSize=newEnd-newStart) %>%
  dplyr::mutate(nclust=n()) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(conID,.keep_all = T) %>%
  dplyr::select(-minStart,-maxEnd,-divSize) %>%
  dplyr::rename(minStart=newStart,maxEnd=newEnd,divSize=newDivSize) %>%
  dplyr::select(CHROM,minStart,maxEnd,divSize,STRAIN,nclust)

# Separate the unclustered regions
nojoin_WS_HDRs <- gap_clust_WS_HDRs %>%
  dplyr::group_by(STRAIN,CHROM) %>%
  dplyr::filter(!(dec3g=="join")) %>%
  dplyr::ungroup() %>%
  dplyr::select(CHROM,minStart,maxEnd,divSize,STRAIN) %>%
  dplyr::mutate(nclust=1)

# Join the joined and unclustered regions
# Size filter
# Order by divergence
all_calls_WS_HDRs<- rbind(joinClust_WS_HDRs,nojoin_WS_HDRs) %>%
  dplyr::filter(divSize/1e3 >= 5) %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::mutate(ncalls=n()) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(desc(ncalls),STRAIN,CHROM,minStart) %>%
  dplyr::mutate(sorter=paste0(ncalls,STRAIN)) %>%
  dplyr::mutate(rleID=data.table::rleid(sorter)) %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::mutate(ystrain=cur_group_id()) %>%
  dplyr::ungroup() %>%
  dplyr::rename(contig=CHROM,strain=STRAIN,strain_order=ystrain,ws_hdr_size=divSize, ws_hdr_start = minStart, ws_hdr_end = maxEnd) %>%
  dplyr::select(-nclust,-ncalls,-sorter,-rleID) %>%
  dplyr::arrange(desc(strain_order),contig,ws_hdr_start) %>% # ordering strains from MOST to least divergent based on number of HDR calls 
  dplyr::select(-strain_order) 

# write.table(all_calls_WS_HDRs, "../../tables/liftedover_WS_HDRs.tsv", col.names = T, row.names = F, quote = F, sep = "\t")

# Calculating total sequence classified as divergent in each wild strain
span_ws_hdrs <- all_calls_WS_HDRs %>%
  dplyr::group_by(strain) %>% 
  dplyr::mutate(span_hdrs = sum(ws_hdr_size)) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(desc(span_hdrs)) %>% # This is extremely concordant with SV calls!
  dplyr::distinct(strain, span_hdrs) %>%
  dplyr::mutate(strain = factor(strain, levels = (strain)))

# Plotting total WS HDR span 
ggplot(data = span_ws_hdrs) + 
  geom_col(aes(x = strain, y = span_hdrs / 1e6), fill = 'chocolate4') +
  theme(
    axis.title.x = element_blank(),
    panel.border = element_rect(color = 'black', fill = NA),
    panel.background = element_blank(),
    axis.text.x = element_text(size = 8, color= 'black', angle = 60, hjust = 1),
    axis.title.y = element_text(size = 14, color = 'black'),
    axis.text.y = element_text(size = 12, color = 'black')
  ) +
  labs(y = "Total WS HDR span (Mb)") +
  scale_y_continuous(expand = expansion(mult = c(0, .05)))


# Calculate the number of HDRs per wild strain
hdr_count <- all_calls_WS_HDRs %>% 
  dplyr::select(strain) %>%
  dplyr::group_by(strain) %>%
  dplyr::mutate(ws_hdr_count = n()) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(span_ws_hdrs, by = "strain") %>%
  dplyr::distinct(strain, ws_hdr_count, span_hdrs)

r_val <- cor(hdr_count$ws_hdr_count, hdr_count$span_hdrs, method = "spearman", use = "complete.obs")

# Plotting total WS HDR span by count of HDRs
ggplot(data = hdr_count) + 
  geom_point(aes(x = ws_hdr_count, y = span_hdrs / 1e6, color = strain)) +
  geom_smooth(aes(x = ws_hdr_count, y = span_hdrs / 1e6), method = "lm", se = TRUE, color = "black", linewidth = 1) +
  annotate("text", x = Inf, y = Inf, label = paste0("Spearman ρ = ", round(r_val, 2)), hjust = 1.1, vjust = 1.5, size = 5) +
  theme(
    panel.border = element_rect(color = 'black', fill = NA),
    panel.background = element_blank(),
    axis.title = element_text(size = 14, color = 'black'),
    axis.text = element_text(size = 12, color = 'black')
  ) +
  labs(y = "Total WS HDR span (Mb)", x = "WS HDR count")



# ======================================================================================================================================================================================== #
# Loading in orthogroups and classifying gene sets 
# ======================================================================================================================================================================================== #
# Loading in matrix of ortholog count orthogroup matrix
all_relations <- readr::read_tsv("../../processed_data/orthology/OG_relations_matrix_count.tsv")

# The frequency of private orthogroupss
private_freq = (1 / ((length(colnames(all_relations))) -1))

# Classifying pangenome gene sets
class <- all_relations %>%
  dplyr::mutate(across(2:(ncol(.)), ~ ifelse(. >= 1, 1, .))) %>%
  dplyr::mutate(sum = rowSums(across(-1, ~ ., .names = NULL), na.rm = TRUE)) %>%
  dplyr::mutate(freq = (sum / (length(colnames(all_relations)) - 1))) %>%
  dplyr::mutate(
    class = case_when(
      freq == 1 ~ "core",
      freq > private_freq & freq < 1 ~ "accessory",
      freq == private_freq ~ "private",
      TRUE ~ "undefined")) %>%
  dplyr::select(Orthogroup,class)

# Loading in orthogroup matrix 
all_OGs <- readr::read_tsv("../../tables/all_OGs_matrix.tsv")

# Adding pangenome gene set classification to OG matrix with gene names
all_class <- all_OGs %>% dplyr::left_join(class, by = "Orthogroup") 

# Creating a long table, where every row has a gene in the pangenome and the gene set that it contributes to
long_class <- all_class %>%
  tidyr::pivot_longer(
    cols = -c(Orthogroup, class),
    names_to = "strain",
    values_to = "gene",
    values_drop_na = TRUE) %>%
  tidyr::separate_rows(gene, sep = ",\\s*") %>%
  dplyr::select(strain, Orthogroup, gene, class) %>%
  dplyr::mutate(gene = sub("\\.[^.]*$", "", gene))

# write.table(long_class, "../../processed_data/orthology/all_genes_class_OGs.tsv",  sep = '\t', quote = F, col.names = T, row.names = F)


# Loading in all genes in pangenome
genes_strain <- readr::read_tsv("../../processed_data/genome_resources/annotation/140Ws_CGC1_longestIsoGenes_BRAKER.tsv", col_names = c("seqid","source", "type", "start", "end", "score", "strand", "phase", "attributes", "strain"))
N2_gff <- ape::read.gff("../../processed_data/genome_resources/annotation/c_elegans.PRJNA13758.WS283.csq.PCfeaturesOnly.longest.gff3") %>% dplyr::mutate(strain="N2")
genes_strain <- rbind(genes_strain,N2_gff)
all_genes_strain <- genes_strain %>%
  dplyr::mutate(attributes = gsub("ID=gene:","",attributes)) %>%
  dplyr::mutate(attributes = gsub("ID=","",attributes)) %>%
  dplyr::mutate(attributes = sub(";.*", "", attributes)) %>%
  dplyr::filter(type == "gene") %>%
  dplyr::rename(gene = attributes)

# Adding the coordinates of every gene in the pangenome
genes_class <- all_genes_strain %>%
  dplyr::left_join(long_class, by = c("gene","strain"))


# ======================================================================================================================================================================================== #
# Pulling WS genes that are in lifted-over WS HDRs #
# ======================================================================================================================================================================================== #
# Prepping input of WS HDRs and WS genes
ws_hdrs <- data.table::as.data.table(all_calls_WS_HDRs %>% dplyr::select(strain, contig, ws_hdr_start, ws_hdr_end) %>% dplyr::rename(start = ws_hdr_start, end = ws_hdr_end))
ws_genes <- genes_class %>% dplyr::select(strain, seqid, start, end, gene, Orthogroup, class) %>% dplyr::filter(strain != "N2" & strain != "CGC1") %>% dplyr::rename(contig = seqid) %>% data.table::as.data.table()

# Setting the keys for foverlaps
data.table::setkey(ws_hdrs, strain, contig, start, end)
data.table::setkey(ws_genes, strain, contig, start, end)

ws_genes_hdrs <- foverlaps(
  x = ws_genes,
  y = ws_hdrs,
  type = "any",
  nomatch = NA) 

# Writing a table of all WS genes in HDRs
ws_hdr_genes <- ws_genes_hdrs %>%
  dplyr::filter(!is.na(start)) %>%
  dplyr::select(strain, gene, Orthogroup, class)

# write.table(ws_hdr_genes, "../../tables/wild_strain_genes_inHDRs.tsv", sep = '\t', quote = F, col.names = T, row.names = F)

# Calculating the number of genes in each gene set for every strain
ws_genes_count <- ws_genes %>%
  dplyr::group_by(strain, class) %>%
  dplyr::mutate(ws_class_count = n()) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(strain,class,ws_class_count) 

# Pulling stats of genes in HDRs and proportion of each gene set
ws_genes_hdrs_stats <- ws_genes_hdrs %>%
  dplyr::filter(!is.na(start)) %>%
  dplyr::group_by(strain, class) %>% 
  dplyr::mutate(ws_class_count_inHDR = n()) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(strain, class, ws_class_count_inHDR)

# Adding stats for the proportion of genes that are in HDRs for each strain and their contribution to each pangenome gene set
ws_genes_hdrs_stats <- ws_genes_count %>%
  dplyr::left_join(ws_genes_hdrs_stats, by = c("strain", "class")) %>%
  dplyr::mutate(ws_class_count_inHDR = ifelse(is.na(ws_class_count_inHDR),0, ws_class_count_inHDR)) %>%
  dplyr::group_by(strain) %>%
  dplyr::mutate(prop_genes_inHDRs_class = ws_class_count_inHDR / ws_class_count) %>%
  dplyr::distinct(strain,class,ws_class_count, ws_class_count_inHDR, prop_genes_inHDRs_class) %>%
  dplyr::group_by(strain) %>%
  dplyr::mutate(gene_set_prop_inHDR = sum(prop_genes_inHDRs_class)) %>%
  dplyr::mutate(non_HDR_geneSet = 1 - gene_set_prop_inHDR) %>%
  dplyr::mutate(ws_total_gene_count = sum(ws_class_count),
                ws_total_inHDR_gene_count = sum(ws_class_count_inHDR),
                prop_total_ws_genes_inHDRs = ws_total_inHDR_gene_count / ws_total_gene_count) %>%
  dplyr::ungroup() 

# Proportion of genes in HDRs for each strain 
prop_genes_in_hdrs <- ws_genes_hdrs_stats %>%
  dplyr::distinct(strain, prop_total_ws_genes_inHDRs) %>%
  dplyr::arrange(desc(prop_total_ws_genes_inHDRs)) %>%
  dplyr::mutate(strain = factor(strain, levels = strain))

# Plotting stats of WS genes in HDRs 
prop_size_genes_inHDRs <- ggplot(data = prop_genes_in_hdrs) + 
  geom_col(aes(x = strain, y = prop_total_ws_genes_inHDRs * 100), fill = 'gray60', color = 'black') +
  geom_point(data = span_ws_hdrs, aes(x = strain, y = (span_hdrs / 1e6)), color = "black", size = 1) +
  geom_line( data = span_ws_hdrs, aes(x = strain, y = (span_hdrs / 1e6) , group = 1), color = "black", linewidth = 0.5) +
  theme(
    axis.title.x = element_blank(),
    panel.border = element_rect(color = 'black', fill = NA),
    panel.background = element_blank(),
    axis.text.x = element_text(size = 3.5, color= 'black', angle = 60, hjust = 1),
    axis.title.y = element_text(size = 10, color = 'black'),
    axis.text.y = element_text(size = 10, color = 'black'),
    plot.margin = margin(t = 8, b = 5, l = 5, r = 5)) +
  labs(y = "Proportion of wild strain genes in HDRs (%)") +
  scale_y_continuous(name = "Proportion of wild strain genes in HDRs (%)", 
                     sec.axis = sec_axis(~. / 1 ,name = "Total wild strain HDRs span (Mb)"), expand = expansion(mult = c(0, .05))) 
prop_size_genes_inHDRs

# ggsave("../../figures/supplementary/prop_size_genes_inHDRs.png", prop_size_genes_inHDRs, width = 7.5, height = 4, dpi = 600)


# Proportions inside and outside HDRs in an admixture-like display:
geneSet_prop <- ws_genes_hdrs_stats %>%
  dplyr::select(strain, class, ws_class_count_inHDR, prop_genes_inHDRs_class, prop_total_ws_genes_inHDRs, gene_set_prop_inHDR, ws_class_count) %>%
  dplyr::mutate(ws_class_outsideHDRs = ws_class_count - ws_class_count_inHDR,
                prop_genes_OUTSIDEhdrs_class = ws_class_outsideHDRs / ws_class_count) %>%
  dplyr::group_by(strain) %>%
  dplyr::mutate(scaling_factor = 1 / gene_set_prop_inHDR,
                gene_set_prop_OUTSIDEHDR = sum(prop_genes_OUTSIDEhdrs_class),
                saling_factor_OUTSIDE = 1 / gene_set_prop_OUTSIDEHDR) %>%
  dplyr::mutate(scaled_geneSet_props = prop_genes_inHDRs_class * scaling_factor,
                scaled_geneSet_props_NONHDR = prop_genes_OUTSIDEhdrs_class * saling_factor_OUTSIDE) %>%
  dplyr::ungroup()

# Order strains by lowest to greatest proportion of WS genes in HDRs
strain_order <- geneSet_prop %>%
  dplyr::filter(class == "private") %>%
  dplyr::arrange(desc(scaled_geneSet_props))%>%   # ascending: smallest -> largest
  dplyr::distinct(strain) %>%
  dplyr::pull(strain)

# Cleaning up for plotting
hdr_nonHDR_prop_geneset <- geneSet_prop %>%
  dplyr::rename(`Gene set` = class) %>%
  dplyr::mutate(`Gene set` = ifelse(`Gene set` == "core","Core",
                                    ifelse(`Gene set` == "accessory", "Accessory", "Private"))) %>%
  dplyr::mutate(
    strain = factor(strain, levels = strain_order),
    `Gene set`  = factor(`Gene set`, levels = c("Core", "Accessory", "Private"))) %>%
  dplyr::group_by(strain) %>%
  dplyr::mutate(scaled_geneSet_props_HDR  = scaled_geneSet_props) %>%
  dplyr::ungroup() 

# Plotting relative fraction in HDRs
hdr <- ggplot(data = hdr_nonHDR_prop_geneset) +
  geom_col(aes(x = strain, y = scaled_geneSet_props_HDR, fill = `Gene set`), alpha = 0.5, width = 1, color = 'black', linewidth = 0.1) +
  scale_fill_manual(values = c(
    "Core" = "green4",
    "Accessory" = "#DB6333",
    "Private" = "magenta3"
  )) +
  theme(
    axis.text.x = element_blank(),
    legend.position = 'none',
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title = element_blank(),
    plot.margin = margin(l = 40, r = 40, t = 10),
    axis.text.y = element_text(size = 10, color = 'black')
  ) +
  scale_y_continuous(expand = c(0,0), breaks = c(0.25, 0.5, 0.75, 1)) 

# Plotting relative fraction outside HDRs
nonhdr <- ggplot(data = hdr_nonHDR_prop_geneset) +
  geom_col(aes(x = strain, y = scaled_geneSet_props_NONHDR, fill = `Gene set`), alpha = 0.5, width = 1, color = 'black', linewidth = 0.1) +
  scale_fill_manual(values = c(
    "Core" = "green4",
    "Accessory" = "#DB6333",
    "Private" = "magenta3"
  )) +
  theme(
    axis.text.x = element_blank(),
    legend.position = 'none',
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.margin = margin(l = 40, r = 40, b = 10),
    axis.text.y = element_text(size = 10, color = 'black')
  ) +
  scale_y_continuous(expand = c(0,0))

aligned <- cowplot::align_plots(hdr, nonhdr, align = "v", axis = "lr")

# Each panel should be no more than 2.5in high and 6in wide!!!!!!!!!!!!!!! - into panels A) and B)
fraction_gene_class_HDRs <- cowplot::plot_grid(cowplot::plot_grid(
  aligned[[1]],aligned[[2]],
  nrow = 2) + draw_label("Relative fraction of each gene set", x=0.02, y=0.5, vjust= 1.5, angle=90, size = 11, color = 'black') +
    draw_label("Genes in HDRs", x=0.97, y=0.75, vjust= 1.5, angle=270, size = 11, color = 'black') +
    draw_label("Genes not in HDRs", x=0.97, y=0.25, vjust= 1.5, angle=270, size = 11, color = 'black'))

fraction_gene_class_HDRs

# ggsave("../../figures/supplementary/fraction_gene_class_HDRs.png", fraction_gene_class_HDRs, width = 7.5, height = 5, dpi = 600)


############### Only HDR! #########################################
relative_gene_set_HDR <- ggplot(data = hdr_nonHDR_prop_geneset) +
  geom_col(aes(x = strain, y = scaled_geneSet_props_HDR, fill = `Gene set`), alpha = 0.5, width = 1, color = 'black', linewidth = 0.1) +
  scale_fill_manual(values = c(
    "Core" = "green4",
    "Accessory" = "#DB6333",
    "Private" = "magenta3"
  )) +
  theme(
    axis.text.x = element_blank(),
    legend.position = 'none',
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 11, color = 'black'),
    plot.margin = margin(l = 10, r = 10, t = 20, b = 20),
    axis.text.y = element_text(size = 10, color = 'black')
  ) +
  scale_y_continuous(expand = c(0,0)) +
  labs(y = "Relative fraction of gene set in HDRs")
relative_gene_set_HDR

# write.table(hdr_nonHDR_prop_geneset, "../../processed_data/hdr_liftover/HDR_nonHDR_relativeFract_geneset.tsv", col.names = T, row.names = F, quote = F, sep = "\t")
###################################################################



# Average proportion of WS genes in each gene set among all wild strains
geneSet_prop_average <- ws_genes_hdrs_stats %>%
  dplyr::select(strain,class,ws_class_count,ws_class_count_inHDR) %>%
  dplyr::group_by(class) %>%
  dplyr::mutate(average_genesetGenes = sum(ws_class_count) / 140,
                average_genesetHDRgenes = sum(ws_class_count_inHDR / 140)) %>%
  dplyr::distinct(class, average_genesetGenes, average_genesetHDRgenes) %>%
  dplyr::mutate(average_prop = (average_genesetHDRgenes / average_genesetGenes) * 100)
# On average:
### 3.4% of core genes are in HDRs among wild strains
### 12.6% of accessory genes are in HDRs among wild strains
### 14.8% of private genes are in HDRs among wild strains

# In the entire pangenome:
hdr_genes_pangenome <- ws_genes_hdrs_stats %>%
  dplyr::select(strain,class,ws_class_count,ws_class_count_inHDR, ws_total_gene_count) %>%
  dplyr::summarise(total_genes = sum(ws_class_count),
                   total_hdr_genes = sum(ws_class_count_inHDR),
                   av_genes_WSs = mean(ws_total_gene_count)) %>%
  dplyr::mutate(prop = total_hdr_genes / total_genes * 100)
### 6.50 %
# In the entire pangenome - by class:
hdr_genes_pangenome_class <- ws_genes_hdrs_stats %>%
  dplyr::select(strain,class,ws_class_count,ws_class_count_inHDR) %>%
  dplyr::group_by(class) %>%
  dplyr::summarise(total_genes = sum(ws_class_count),
                   total_hdr_genes = sum(ws_class_count_inHDR)) %>%
  dplyr::mutate(prop = total_hdr_genes / total_genes * 100)
### 70,320 / 2,101,058 core
### 133,445 / 1,055,492 accessory
### 2,459 / 16,560 private


# ======================================================================================================================================================================================== #
# Creating a final stats table #
# ======================================================================================================================================================================================== #
# Creating a stats table to summarize HDR lift-over
n2_hdr_stats <- strain_hdr %>%
  dplyr::distinct(strain, chrom, og_hdr_start, og_hdr_end) %>%
  dplyr::mutate(n2_hdr_size = og_hdr_end - og_hdr_start) %>%
  dplyr::group_by(strain) %>%
  dplyr::mutate(n2_hdr_count = n()) %>%
  dplyr::mutate(n2_hdr_span = sum(n2_hdr_size)) %>%
  dplyr::mutate(largest_n2_hdr = max(n2_hdr_size)) %>%
  dplyr::mutate(mean_n2_hdr_size = n2_hdr_span / n2_hdr_count) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(strain,n2_hdr_count, n2_hdr_span, largest_n2_hdr, mean_n2_hdr_size)

stats <- all_calls_WS_HDRs %>%
  dplyr::group_by(strain) %>%
  dplyr::mutate(largest_ws_hdr = max(ws_hdr_size)) %>%
  dplyr::select(strain, largest_ws_hdr) %>%
  dplyr::group_by(strain) %>%
  dplyr::mutate(ws_hdr_count = n()) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(span_ws_hdrs, by = "strain") %>%
  dplyr::distinct(strain, ws_hdr_count, span_hdrs, largest_ws_hdr) %>%
  dplyr::mutate(mean_ws_hdr_size = span_hdrs / ws_hdr_count) %>%
  dplyr::left_join(n2_hdr_stats, by = "strain") %>%
  dplyr::rename(span_ws_hdrs = span_hdrs)

# Stats on WS genes in HDRs
addition <- ws_genes_hdrs_stats %>%
  dplyr::select(strain, ws_total_gene_count, ws_total_inHDR_gene_count) %>%
  dplyr::distinct()

final_ws_genes_hdr_stats <- ws_genes_hdrs_stats %>%
  dplyr::select(strain,class,ws_class_count, ws_class_count_inHDR) %>%
  tidyr::pivot_wider(
    id_cols = strain,
    names_from = class,
    values_from = c(ws_class_count, ws_class_count_inHDR),
    names_glue = "{.value}_{class}") %>%
  dplyr::left_join(addition, by = "strain")

final_stats <- stats %>% dplyr::left_join(final_ws_genes_hdr_stats, by = "strain")

# write.table(final_stats, "../../tables/WS_HDR_liftover_finalstats.tsv", sep = '\t', quote = F, col.names = T, row.names = F)
