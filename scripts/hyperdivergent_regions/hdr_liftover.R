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
hdrs <- readr::read_tsv("../../processed_data/genome_resources/annotation/20250625_c_elegans_divergent_regions_strain.bed", col_names = c("chrom", "start", "end", "strain"))
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
HDOI1 <- c("ECA3005","ECA3005V1842400019264000") # A great example of a lift-over!

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
               aes(x = leadS / 1e6, xend = leadE / 1e6, y = leadWSS / 1e6, yend = leadWSE / 1e6, color = "leading"), linewidth = 1) +
  geom_segment(data = nucmer_mark_extend %>%
                 dplyr::filter(any_extend == T) %>%
                 dplyr::filter(HDRid == HDOI1[2]) %>%
                 dplyr::filter(N2S > hdr_start_extended),
               aes(x = lagS / 1e6, xend = lagE / 1e6, y = lagWSS / 1e6, yend = lagWSE / 1e6, color = "lagging"), linewidth = 1) +
  geom_segment(data = nucmer_extended %>%
                 dplyr::filter(HDRid == HDOI1[2]),
               aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6, color = "w/Extension"), linewidth = 1) +
  geom_segment(data = nucmer_mark_extend %>%
                 dplyr::filter(HDRid == HDOI1[2]),
               aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6, color = "overlapping"), linewidth = 1) +
  facet_wrap(~chrom) +
  scale_color_manual(values = c("leading" = 'red', "overlapping" = "blue", "lagging" = "green3", "w/Extension" = "purple")) +
  theme(panel.border = element_rect(color = 'black', fill = NA),
        panel.background = element_blank(),
        axis.text = element_text(size = 10, color = 'black'),
        axis.title = element_text(size = 11, color = 'black'),
        strip.text = element_text(size = 14, color = 'black'),
        axis.title.x = element_blank(),
        legend.position = 'none') +
  labs(x = "N2 genome position (Mb)", y = "ECA3005 contig position (Mb)")
ex1



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
               aes(x = leadS / 1e6, xend = leadE / 1e6, y = leadWSS / 1e6, yend = leadWSE / 1e6, color = "leading"), linewidth = 1) +
  geom_segment(data = nucmer_mark_extend %>%
                 dplyr::filter(any_extend == T) %>%
                 dplyr::filter(HDRid == HDOI2[2]) %>%
                 dplyr::filter(N2S > hdr_start_extended),
               aes(x = lagS / 1e6, xend = lagE / 1e6, y = lagWSS / 1e6, yend = lagWSE / 1e6, color = "lagging"), linewidth = 1) +
  geom_segment(data = nucmer_extended %>%
                 dplyr::filter(HDRid == HDOI2[2]),
               aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6, color = "w/Extension"), linewidth = 1) +
  geom_segment(data = nucmer_mark_extend %>%
                 dplyr::filter(HDRid == HDOI2[2]),
               aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6, color = "overlapping"), linewidth = 1) +
  facet_wrap(~chrom) +
  scale_color_manual(values = c("leading" = 'red', "overlapping" = "blue", "lagging" = "green3", "w/Extension" = "purple")) +
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
               aes(x = leadS / 1e6, xend = leadE / 1e6, y = leadWSS / 1e6, yend = leadWSE / 1e6, color = "leading"), linewidth = 1) +
  geom_segment(data = nucmer_mark_extend %>%
                 dplyr::filter(any_extend == T) %>%
                 dplyr::filter(HDRid == HDOI3[2]) %>%
                 dplyr::filter(N2S > hdr_start_extended),
               aes(x = lagS / 1e6, xend = lagE / 1e6, y = lagWSS / 1e6, yend = lagWSE / 1e6, color = "lagging"), linewidth = 1) +
  geom_segment(data = nucmer_extended %>%
                 dplyr::filter(HDRid == HDOI3[2]),
               aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6, color = "w/Extension"), linewidth = 1) +
  geom_segment(data = nucmer_mark_extend %>%
                 dplyr::filter(HDRid == HDOI3[2]),
               aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6, color = "overlapping"), linewidth = 1) +
  facet_wrap(~chrom) +
  scale_color_manual(values = c("leading" = 'red', "overlapping" = "blue", "lagging" = "green3", "w/Extension" = "purple")) +
  theme(panel.border = element_rect(color = 'black', fill = NA),
        panel.background = element_blank(),
        axis.text = element_text(size = 10, color = 'black'),
        axis.title = element_text(size = 11, color = 'black'),
        plot.margin = margin(t = 5, b = 5, l = 5, r = 5),
        strip.text = element_text(size = 14, color = 'black'),
        legend.position = 'none') +
  labs(x = "N2 genome position (Mb)", y = "ECA3088 contig position (Mb)")
ex3



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
               aes(x = leadS / 1e6, xend = leadE / 1e6, y = leadWSS / 1e6, yend = leadWSE / 1e6, color = "leading"), linewidth = 1) +
  geom_segment(data = nucmer_mark_extend %>%
                 dplyr::filter(any_extend == T) %>%
                 dplyr::filter(HDRid == HDOI4[2]) %>%
                 dplyr::filter(N2S > hdr_start_extended),
               aes(x = lagS / 1e6, xend = lagE / 1e6, y = lagWSS / 1e6, yend = lagWSE / 1e6, color = "lagging"), linewidth = 1) +
  geom_segment(data = nucmer_extended %>%
                 dplyr::filter(HDRid == HDOI4[2]),
               aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6, color = "w/Extension"), linewidth = 1) +
  geom_segment(data = nucmer_mark_extend %>%
                 dplyr::filter(HDRid == HDOI4[2]),
               aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6, color = "overlapping"), linewidth = 1) +
  facet_wrap(~chrom) +
  scale_color_manual(values = c("leading" = 'red', "overlapping" = "blue", "lagging" = "green3", "w/Extension" = "purple")) +
  theme(panel.border = element_rect(color = 'black', fill = NA),
        panel.background = element_blank(),
        axis.text = element_text(size = 10, color = 'black'),
        axis.title = element_text(size = 11, color = 'black'),
        strip.text = element_text(size = 14, color = 'black'),
        legend.position = "inside",
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.box.background = element_rect(fill = "white", color = 'black'),
        legend.position.inside = c(0.8,0.2),
        legend.text = element_text(size = 10, color = 'black')) +
  labs(x = "N2 genome position (Mb)", y = "ECA3005 contig position (Mb)", color = NULL)
ex4


# Concatenating examples to one final plot
hdr_liftover_examples <- cowplot::plot_grid(
  ex1, ex2, ex3, ex4,
  nrow = 2,
  labels = c("a","b","c","d"))
hdr_liftover_examples

ggsave("/vast/eande106/projects/Lance/THESIS_WORK/manuscript_repos/pangenome_Ce_MS/figures/supplementary/hdr_liftover_examples.png", hdr_liftover_examples, width = 7.5, height = 7.5, dpi = 600 )

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

ggsave("/vast/eande106/projects/Lance/THESIS_WORK/manuscript_repos/pangenome_Ce_MS/figures/supplementary/hdr_liftover_size_comparison.png", hdr_liftover_size, width = 7.5, height = 7.5, dpi = 600 )


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
ortho_genes_dd <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/orthology/elegans/orthofinder/64_core/OrthoFinder/Results_Dec07/Orthogroups/Orthogroups.tsv") %>%
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

# write.table(all_relations, "/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/ws_HDR_liftover/OG_relations_matrix_count.tsv",  sep = '\t', quote = F, col.names = T, row.names = F)

private_freq = (1/(length(strainCol_c2_u)))

class <- all_relations %>%
  dplyr::mutate(across(2:(ncol(.)), ~ ifelse(. >= 1, 1, .))) %>%
  dplyr::mutate(sum = rowSums(across(-1, ~ ., .names = NULL), na.rm = TRUE)) %>%
  dplyr::mutate(freq = (sum / length(strainCol_c2_u))) %>%
  dplyr::mutate(
    class = case_when(
      freq == 1 ~ "core",
      freq > private_freq & freq < 1 ~ "accessory",
      freq == private_freq ~ "private",
      TRUE ~ "undefined")) %>%
  dplyr::select(Orthogroup,class)

all_class <- ortho_genes_dd %>% dplyr::bind_rows(private_OGs) %>% dplyr::left_join(class, by = "Orthogroup") 

long_class <- all_class %>%
  tidyr::pivot_longer(
    cols = -c(Orthogroup, class),
    names_to = "strain",
    values_to = "gene",
    values_drop_na = TRUE) %>%
  tidyr::separate_rows(gene, sep = ",\\s*") %>%
  dplyr::select(strain, Orthogroup, gene, class) %>%
  dplyr::mutate(gene = sub("\\.[^.]*$", "", gene))

# write.table(long_class, "/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/ws_HDR_liftover/all_genes_class_OGs.tsv",  sep = '\t', quote = F, col.names = T, row.names = F)


# Loading in all genes in pangenome
genes_strain <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/geneAnno-nf/142_140WSs_andCGC1_longestIsoGenes.tsv", col_names = c("seqid","source", "type", "start", "end", "score", "strand", "phase", "attributes", "strain")) %>% dplyr::filter(strain != "ECA396")
N2_gff <- ape::read.gff("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/raw_data/assemblies/elegans/gff/longest_isoform/c_elegans.PRJNA13758.WS283.csq.PCfeaturesOnly.longest.gff3") %>% dplyr::mutate(strain="N2")
all_genes_strain <- rbind(genes_strain,N2_gff) %>%
  dplyr::mutate(attributes = gsub("ID=gene:","",attributes)) %>%
  dplyr::mutate(attributes = gsub("ID=","",attributes)) %>%
  dplyr::mutate(attributes = sub(";.*", "", attributes)) %>%
  dplyr::filter(type == "gene") %>%
  dplyr::rename(gene = attributes)

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


# Identifying OGs that HDR genes are in
# ws_hdr_genes_og_class <- ws_genes_hdrs %>%
  # dplyr::filter(!is.na(start)) %>%
  # dplyr::distinct(strain, gene, Orthogroup, class)

# write.table(ws_hdr_genes_og_class, "/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/ws_HDR_liftover/hdr_genes_OG_class.tsv",  sep = '\t', quote = F, col.names = T, row.names = F)


# Number of genes in each gene set in HDRs compared to entire gene set
# ws_genes_class <- ws_genes_hdrs_stats %>%
#   dplyr::group_by(class) %>%
#   dplyr::mutate(ws_class_count_total = sum(ws_class_count),
#                 ws_class_count_inHDR_total = sum(ws_class_count_inHDR)) %>%
#   dplyr::distinct(class, ws_class_count_total, ws_class_count_inHDR_total)

# Writing a table of all WS genes in HDRs
# ws_hdr_genes <- ws_genes_hdrs %>%
#   dplyr::filter(!is.na(start)) %>%
#   dplyr::select(strain, gene, class)
# 
# write.table(ws_hdr_genes, "/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/ws_HDR_liftover/WS_genes_inHDRs.tsv", sep = '\t', quote = F, col.names = T, row.names = F)


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

test <- ws_genes_hdrs_stats %>% dplyr::filter(strain == "PX179") # PX179 has zero private genes in HDRs

prop_genes_in_hdrs <- ws_genes_hdrs_stats %>%
  dplyr::distinct(strain, prop_total_ws_genes_inHDRs) %>%
  dplyr::arrange(desc(prop_total_ws_genes_inHDRs)) %>%
  dplyr::mutate(strain = factor(strain, levels = strain))

# Plotting stats of WS genes in HDRs 
ggplot(data = prop_genes_in_hdrs) + 
  geom_col(aes(x = strain, y = prop_total_ws_genes_inHDRs * 100), fill = 'gray45', color = 'black') +
  geom_point(data = span_ws_hdrs, aes(x = strain, y = (span_hdrs / 1e6)), color = "black", size = 5) +
  geom_line( data = span_ws_hdrs, aes(x = strain, y = (span_hdrs / 1e6) , group = 1), color = "black", linewidth = 1) +
  theme(
    axis.title.x = element_blank(),
    panel.border = element_rect(color = 'black', fill = NA),
    panel.background = element_blank(),
    axis.text.x = element_text(size = 13, color= 'black', angle = 60, hjust = 1),
    axis.title.y = element_text(size = 30, color = 'black', face = 'bold'),
    axis.text.y = element_text(size = 26, color = 'black'),
    plot.margin = margin(t = 8, b = 5, l = 5, r = 5)) +
  labs(y = "Proportion of wild strain genes in HDRs (%)") +
  scale_y_continuous(name = "Proportion of wild strain genes in HDRs (%)", 
                     sec.axis = sec_axis(~. / 1 ,name = "Total wild strain HDRs span (Mb)"), expand = expansion(mult = c(0, .05))) 

# Plotting proportion of each gene set in WS HDR genes
strain_order0 <- prop_genes_in_hdrs %>% dplyr::pull(strain)
  
final_prop_gs_final <- ws_genes_hdrs_stats %>%
  dplyr::left_join(span_ws_hdrs, by = "strain") %>%
  dplyr::select(strain,class,ws_class_count_inHDR, ws_total_gene_count, span_hdrs, prop_genes_inHDRs_class, prop_total_ws_genes_inHDRs, gene_set_prop_inHDR) %>%
  dplyr::mutate(prop_gene_class_inHDR = (ws_class_count_inHDR / ws_total_gene_count) * 100,
                span_hdrs = span_hdrs / 1e6) %>%
  dplyr::group_by(strain) %>%
  dplyr::mutate(total_prop = sum(prop_gene_class_inHDR),
                scaling_factor = 1 / gene_set_prop_inHDR,
                scaled_geneSet_props = prop_genes_inHDRs_class * scaling_factor,
                scaled_geneSet_props_final = prop_total_ws_genes_inHDRs * scaled_geneSet_props,
                test = sum(scaled_geneSet_props_final)) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(desc(total_prop), desc(class)) %>%
  dplyr::rename(`Gene set` = class) %>%
  dplyr::mutate(`Gene set` = ifelse(`Gene set` == "core","Core",
                                    ifelse(`Gene set` == "accessory", "Accessory", "Private"))) %>%
  dplyr::mutate(
    strain = factor(strain, levels = strain_order0),
    `Gene set`  = factor(`Gene set`, levels = c("Core", "Accessory", "Private"))) 

ggplot(data = final_prop_gs_final) + 
  geom_col(aes(x = strain, y = prop_gene_class_inHDR, fill = `Gene set`)) +
  geom_point(aes(x = strain, y = span_hdrs), color = "black", size = 2) +
  geom_line( data = span_ws_hdrs, aes(x = strain, y = (span_hdrs / 1e6) , group = 1), color = "black", linewidth = 1) +
  scale_fill_manual(values = c(
    "Core" = "green4",
    "Accessory" = "#DB6333",
    "Private" = "magenta3"
  )) +
  theme(
    axis.title.x = element_blank(),
    panel.border = element_rect(color = 'black', fill = NA),
    panel.background = element_blank(),
    axis.text.x = element_text(size = 8, color= 'black', angle = 60, hjust = 1),
    axis.title.y = element_text(size = 14, color = 'black'),
    axis.text.y = element_text(size = 12, color = 'black'),
    plot.margin = margin(t = 5, b = 5, l = 5, r = 5)) +
  labs(y = "Proportion of WS genes in HDRs (%)") +
  scale_y_continuous(name = "Proportion of WS genes in HDRs (%)", 
                     sec.axis = sec_axis(~. / 1, name = "Total WS HDR span (Mb)"), expand = expansion(mult = c(0, .05))) 

ggplot(data = final_prop_gs_final) + 
  geom_col(aes(x = strain, y = scaled_geneSet_props_final * 100, fill = `Gene set`)) +
  geom_point(aes(x = strain, y = span_hdrs), color = "black", size = 2) +
  geom_line( data = span_ws_hdrs, aes(x = strain, y = (span_hdrs / 1e6) , group = 1), color = "black", linewidth = 1) +
  scale_fill_manual(values = c(
    "Core" = "green4",
    "Accessory" = "#DB6333",
    "Private" = "magenta3"
  )) +
  theme(
    axis.title.x = element_blank(),
    panel.border = element_rect(color = 'black', fill = NA),
    panel.background = element_blank(),
    legend.position = "inside",
    legend.position.inside = c(0.85,0.85),
    legend.title = element_text(size = 20, color = 'black'),
    legend.text = element_text(size = 18, color = 'black'),
    axis.text.x = element_text(size = 12, color= 'black', angle = 60, hjust = 1),
    axis.title.y = element_text(size = 24, color = 'black', face = 'bold'),
    axis.text.y = element_text(size = 18, color = 'black'),
    plot.margin = margin(t = 5, b = 5, l = 5, r = 5)) +
  labs(y = "Proportion of WS genes in HDRs (%)") +
  scale_y_continuous(name = "Proportion of WS genes in HDRs (%)",
                     sec.axis = sec_axis(~. / 1, name = "Total wild strain HDRs span (Mb)"), expand = expansion(mult = c(0, .05)))


# Pie charts to display the proportion of genes in HDRs contributing to each gene set (scaled)
geneSet_prop <- ws_genes_hdrs_stats %>%
  dplyr::select(strain, class, prop_genes_inHDRs_class, prop_total_ws_genes_inHDRs, gene_set_prop_inHDR) %>%
  dplyr::group_by(strain) %>%
  dplyr::mutate(scaling_factor = 1 / gene_set_prop_inHDR) %>%
  dplyr::mutate(scaled_geneSet_props = prop_genes_inHDRs_class * scaling_factor) %>%
  dplyr::ungroup()

# Order the pies (facets) by lowest to greatest proportion of WS genes in HDRs
strain_order <- geneSet_prop %>%
  dplyr::arrange(prop_total_ws_genes_inHDRs) %>%   # ascending: smallest -> largest
  dplyr::distinct(strain) %>%
  dplyr::pull(strain)

# Prep data for pies (ensure all 3 classes exist per strain)
pie_df <- geneSet_prop %>%
  dplyr::rename(`Gene set` = class) %>%
  dplyr::mutate(`Gene set` = ifelse(`Gene set` == "core","Core",
                               ifelse(`Gene set` == "accessory", "Accessory", "Private"))) %>%
  dplyr::mutate(
    strain = factor(strain, levels = strain_order),
    `Gene set`  = factor(`Gene set`, levels = c("Core", "Accessory", "Private"))) %>%
  # tidyr::complete(strain, `Gene set`, fill = list(scaled_geneSet_props = 0)) %>%
  dplyr::group_by(strain) %>%
  dplyr::mutate(pie_frac  = scaled_geneSet_props) %>%
  dplyr::ungroup() 

# Plot 140 pies, 14 high x 10 wide
ggplot(pie_df, aes(x = "", y = pie_frac, fill = `Gene set`), alpha = 0.7) +
  geom_col(width = 1, color = "white", linewidth = 0.2) +
  scale_fill_manual(values = c(
    "Core" = "green4",
    "Accessory" = "#DB6333",
    "Private" = "magenta3"
  )) +
  coord_polar(theta = "y") +
  facet_wrap(~ strain, nrow = 7, ncol = 20, as.table = FALSE) +
  theme_void() +
  theme(
    strip.text = element_text(size = 16, color = "black"),
    legend.position = "bottom",
    legend.title = element_text(size = 22, color = "black"),
    legend.text = element_text(size = 20, color = "black"),
    # plot.title = element_text(size = 26, hjust = 0.5, face = "bold", margin = margin(b = 10)),
    # plot.subtitle = element_text(size = 24, hjust = 0.5, margin = margin(b = 10)),
    plot.margin = margin(t = 15)) +
  labs(
    # title = "Contributions of genes in HDRs to each gene set",
    fill = "Gene set:  ") 
    # subtitle = "• ECA1493: 4,344 genes (most)\n• NIC2: 81 genes (fewest)")










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

# Order the pies (facets) by lowest to greatest proportion of WS genes in HDRs
strain_order <- geneSet_prop %>%
  dplyr::filter(class == "private") %>%
  dplyr::arrange(desc(scaled_geneSet_props))%>%   # ascending: smallest -> largest
  dplyr::distinct(strain) %>%
  dplyr::pull(strain)

# Prep data for pies (ensure all 3 classes exist per strain)
hdr_nonHDR_prop_geneset <- geneSet_prop %>%
  dplyr::rename(`Gene set` = class) %>%
  dplyr::mutate(`Gene set` = ifelse(`Gene set` == "core","Core",
                                    ifelse(`Gene set` == "accessory", "Accessory", "Private"))) %>%
  dplyr::mutate(
    strain = factor(strain, levels = strain_order),
    `Gene set`  = factor(`Gene set`, levels = c("Core", "Accessory", "Private"))) %>%
  dplyr::group_by(strain) %>%
  dplyr::mutate(pie_frac  = scaled_geneSet_props) %>%
  dplyr::ungroup() 

hdr <- ggplot(data = hdr_nonHDR_prop_geneset) +
  geom_col(aes(x = strain, y = pie_frac, fill = `Gene set`), width = 1, color = 'black') +
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
    axis.text.y = element_text(size = 22, color = 'black')
  ) +
  scale_y_continuous(expand = c(0,0), breaks = c(0.25, 0.5, 0.75, 1)) 
hdr


nonhdr <- ggplot(data = hdr_nonHDR_prop_geneset) +
  geom_col(aes(x = strain, y = scaled_geneSet_props_NONHDR, fill = `Gene set`), width = 1, color = 'black') +
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
    axis.text.y = element_text(size = 22, color = 'black')
  ) +
  scale_y_continuous(expand = c(0,0))
nonhdr

aligned <- cowplot::align_plots(hdr, nonhdr, align = "v", axis = "lr")

# Each panel should be no more than 2.5in high and 6in wide!!!!!!!!!!!!!!! - into panels A) and B)
cowplot::plot_grid(cowplot::plot_grid(
  aligned[[1]],aligned[[2]],
  nrow = 2) + draw_label("Relative contribution to each gene set", x=-0.005, y=0.5, vjust= 1.5, angle=90, size = 32, color = 'black', fontface = 'bold') +
    draw_label("Genes in HDRs", x=0.9999, y=0.75, vjust= 1.5, angle=270, size = 32, color = 'black', fontface = 'bold') +
    draw_label("Genes not in HDRs", x=0.9999, y=0.25, vjust= 1.5, angle=270, size = 32, color = 'black', fontface = 'bold'))


############### Only HDR! #########################################
ggplot(data = hdr_nonHDR_prop_geneset) +
  geom_col(aes(x = strain, y = pie_frac, fill = `Gene set`), width = 1, color = 'black') +
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
    axis.title.y = element_text(size = 26, color = 'black', face = 'bold'),
    plot.margin = margin(l = 10, r = 10, t = 20, b = 20),
    axis.text.y = element_text(size = 22, color = 'black')
  ) +
  scale_y_continuous(expand = c(0,0)) +
  labs(y = "Relative fraction of gene set in HDRs")
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

# Look at gene density inside of HDRs compared to genome-wide average
genome_sizes <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/assembly-nf/all_assemblies_sheet/Ce_strains_genome_sizes.tsv", col_names = c("strain", "genome_size")) %>%
  dplyr::filter(strain %in% WSs) 

# Identifying what proportion of WS genomes have gene models predicted
gene_summed_length <- all_genes_strain %>%
  dplyr::select(start,end,strain) %>%
  dplyr::filter(strain %in% WSs) %>%
  dplyr::mutate(gene_size = end - start) %>%
  dplyr::group_by(strain) %>%
  dplyr::mutate(summed_gene_size = sum(gene_size)) %>%
  dplyr::mutate(gene_count = n()) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(strain,gene_count, summed_gene_size) %>%
  dplyr::left_join(genome_sizes, by = "strain") %>%
  dplyr::mutate(prop_coding = summed_gene_size / genome_size)

# Total summed length of genes in HDRs
ws_genes_summed <- ws_genes_hdrs %>%
  dplyr::filter(!is.na(start)) %>%
  dplyr::select(strain,start,i.start,i.end) %>%
  dplyr::mutate(gene_size = i.end - i.start) %>%
  dplyr::group_by(strain) %>%
  dplyr::mutate(summed_genes_inHDRs = sum(gene_size)) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(strain,summed_genes_inHDRs)

# Proportion of HDRs that are coding
prop_hdr_coding <- all_calls_WS_HDRs %>%
  dplyr::select(strain,ws_hdr_size) %>%
  dplyr::group_by(strain) %>%
  dplyr::mutate(summed_hdr_size = sum(ws_hdr_size)) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(strain,summed_hdr_size) %>%
  dplyr::left_join(ws_genes_summed, by = 'strain') %>%
  dplyr::mutate(prop_coding_HDR = summed_genes_inHDRs / summed_hdr_size) %>%
  dplyr::left_join(gene_summed_length, by = "strain") %>%
  dplyr::mutate(fold_enrich_HDR = prop_coding_HDR / prop_coding) # greater than 1 for every wild strain
# All HDRs are more gene dense compared to the genome-wide average
# (summed gene length in HDR / summed length of HDRs) / (summed_genes in entire genome / genome size)

# What is relationship between wild strain genome size and proprtion of genome that is hyper-divergent?
genome_size_hdr <- prop_hdr_coding %>%
  dplyr::select(strain,genome_size,gene_count,summed_hdr_size) %>%
  dplyr::arrange(desc(genome_size)) %>%
  dplyr::mutate(strain = factor(strain, levels = strain))

# Relationship between genome size and span of HDRs
ggplot(data = genome_size_hdr) + 
  geom_col(aes(x = strain, y = genome_size / 1e6, fill = summed_hdr_size / 1e6)) +
  scale_fill_gradient(low = "gold", high = "blue") +
  theme(
    axis.title.x = element_blank(),
    panel.border = element_rect(color = 'black', fill = NA),
    panel.background = element_blank(),
    legend.title = element_text(size = 12, color = 'black'),
    legend.text = element_text(size = 11, color = 'black'),
    axis.text.x = element_text(size = 8, color= 'black', angle = 60, hjust = 1),
    axis.title.y = element_text(size = 14, color = 'black'),
    axis.text.y = element_text(size = 12, color = 'black'),
    plot.margin = margin(t = 5, b = 5, l = 5, r = 5)) +
  labs(y = "Wild strain genome size (Mb)", fill = 'Hyper-divergent genome size (Mb)') +
  scale_y_continuous(expand = c(0,0)) 

# Relationship between genome size and predicted PC genes
ggplot(data = genome_size_hdr) + 
  geom_col(aes(x = strain, y = genome_size / 1e6, fill = gene_count)) +
  scale_fill_gradient(low = "plum1", high = "blue") +
  theme(
    axis.title.x = element_blank(),
    panel.border = element_rect(color = 'black', fill = NA),
    panel.background = element_blank(),
    legend.title = element_text(size = 12, color = 'black'),
    legend.text = element_text(size = 11, color = 'black'),
    axis.text.x = element_text(size = 8, color= 'black', angle = 60, hjust = 1),
    axis.title.y = element_text(size = 14, color = 'black'),
    axis.text.y = element_text(size = 12, color = 'black'),
    plot.margin = margin(t = 5, b = 5, l = 5, r = 5)) +
  labs(y = "Wild strain genome size (Mb)", fill = 'Wild strain gene count') +
  scale_y_continuous(expand = c(0,0)) 

r_val2 <- cor(genome_size_hdr$genome_size, genome_size_hdr$summed_hdr_size, method = "spearman", use = "complete.obs")

hdr_genes <- ws_genes_hdrs_stats %>%
  dplyr::distinct(strain,ws_total_inHDR_gene_count)

# Looking at the relationship among genome size, spans of HDRs, and predicted gene count
ggplot(data = genome_size_hdr %>% dplyr::left_join(hdr_genes, by = 'strain')) +
  geom_point(aes(x = genome_size / 1e6, y = summed_hdr_size / 1e6, fill = gene_count, size = ws_total_inHDR_gene_count), shape = 21) +
  geom_smooth(aes(x = genome_size / 1e6, y = summed_hdr_size / 1e6), method = "lm", se = TRUE, color = "black", linewidth = 1) +
  annotate("text", x = Inf, y = Inf, label = paste0("Spearman ρ = ", round(r_val2, 2)), hjust = 1.1, vjust = 1.5, size = 5) +
  ggrepel::geom_text_repel(data = dplyr::filter(genome_size_hdr, genome_size > 110e6 | summed_hdr_size > 13e6),
                           aes(x = genome_size / 1e6, y = summed_hdr_size / 1e6, label = strain), size = 4, max.overlaps = Inf, point.padding = 12) +
  scale_fill_gradient(low = "yellow", high = "red") +
  scale_size_continuous(name = "Wild strain HDR gene count", breaks = c(1000, 2500, 4000), range = c(1, 8)) +  theme(
    panel.border = element_rect(color = 'black', fill = NA),
    panel.background = element_blank(),
    legend.title = element_text(size = 12, color = 'black'),
    legend.text = element_text(size = 11, color = 'black'),
    axis.title = element_text(size = 16, color = 'black'),
    axis.text = element_text(size = 14, color = 'black'),
    plot.margin = margin(t = 5, b = 5, l = 5, r = 5)) +
  labs(x = "Wild strain genome size (Mb)", y = "Hyper-divergent genome size (Mb)", fill = "Wild strain gene count") +
  coord_cartesian(xlim = c(102,113))
  # scale_x_continuous(expand = c(0,0))

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

# write.table(final_stats, "/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/ws_HDR_liftover/WS_HDR_liftover_finalstats.tsv", sep = '\t', quote = F, col.names = T, row.names = F)












































############### TESTING ################################################################
# OLD METHODS USING SLOPE CALCULATIONS

# 
# SOI <- "ECA3088"
# 
# nucmer_ranges <- nucmer_ranges %>% dplyr::filter(strain == SOI)
# strain_hdr <- all_regions %>% dplyr::filter(strain == SOI) %>% 
#   dplyr::rename(og_hdr_start = start, og_hdr_end = end) %>% 
#   dplyr::mutate(start = ifelse(og_hdr_start >= 5000, og_hdr_start - 5000, og_hdr_start), end = og_hdr_end + 5000) ##################  EXTENDING HDR BOUNDARIES BY 5 KB ON BOTH SIDES TO TRY AND PULL EXTREMELY DIVERGENT REGIONS
# 
# strain_gene_class <- genes_class %>% dplyr::filter(strain == SOI) %>%
#   dplyr::group_by(strain, class) %>%
#   dplyr::mutate(strain_class_count = n()) %>%
#   dplyr::ungroup() %>%
#   dplyr::distinct(seqid,start,end,gene,strain,class,strain_class_count)
# 
# data.table::setkey(nucmer_ranges, chrom, start, end)
# data.table::setkey(strain_hdr, chrom, start, end)
# 
# hdr_aln <- data.table::foverlaps(
#   x = strain_hdr,
#   y = nucmer_ranges,
#   type = "any" # if the start/end of any alignment is within an HDR
# ) %>%
#   dplyr::rename(hdr_start_extended = i.start, hdr_end_extended = i.end, N2S = start, N2E = end) %>%
#   dplyr::select(-i.strain) # 10,333 genes!
# 
# num_hdrs <- nrow(strain_hdr)
# num_hdr_aln <- nrow(hdr_aln %>% dplyr::distinct(hdr_start_extended))
# strain_id <- hdr_aln %>% dplyr::distinct(strain) %>% dplyr::pull()
# print(paste0(num_hdr_aln, " / ", num_hdrs, " HDRs have overlapping ", strain_id, " alignments"))
# 
# # Test plot
# ggplot(data = hdr_aln) +
#   geom_rect(aes(xmin = og_hdr_start / 1e6, xmax = og_hdr_end / 1e6, ymin = -Inf, ymax = Inf), fill = "gray", alpha = 0.1) +
#   geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6,  color = contig), size = 1) +
#   facet_wrap(~chrom) +
#   theme(
#     panel.border = element_rect(color = 'black', fill = NA),
#     panel.background = element_blank()
#   ) +
#   labs(x = "N2 genome position (Mb)", y = "WS contig position (Mb)", title = paste0(hdr_aln$strain))
# 
# # ggplot(data = hdr_aln %>% dplyr::filter(chrom == "IV" & og_hdr_start > 12000000 & og_hdr_end < 14000000)) +
# #   geom_rect(aes(xmin = og_hdr_start / 1e6, xmax = og_hdr_end / 1e6, ymin = -Inf, ymax = Inf), fill = "gray", alpha = 0.1) +
# #   geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6,  color = contig), size = 1) +
# #   facet_wrap(~chrom) +
# #   theme(
# #     panel.border = element_rect(color = 'black', fill = NA),
# #     panel.background = element_blank()
# #   ) +
# #   labs(x = "N2 genome position (Mb)", y = "WS contig position (Mb)", title = paste0(hdr_aln$strain))
# 
# 
# nucmer_longest <- hdr_aln %>% # equivalent of tigFilt from haplotypePlotter.R
#   dplyr::group_by(og_hdr_start) %>%
#   dplyr::mutate(nalign = n()) %>%
#   dplyr::ungroup() %>%
#   dplyr::group_by(og_hdr_start, contig) %>%
#   dplyr::mutate(ntig = n()) %>%
#   dplyr::mutate(tigsize = sum(L2)) %>% # summing the number of alignments that overlap with an N2 gene... some contigs align to a single gene many times
#   dplyr::ungroup() %>% 
#   dplyr::group_by(og_hdr_start) %>%
#   dplyr::filter(tigsize == max(tigsize)) %>%
#   dplyr::ungroup() %>%
#   dplyr::rename(longest_contig = contig) %>%
#   dplyr::group_by(og_hdr_start) %>%
#   dplyr::filter(LENQ == max(LENQ)) %>% # to filter out alignments that are the same size, but from different contigs
#   dplyr::ungroup()
# 
# ggplot(data = nucmer_longest) +
#   geom_rect(aes(xmin = og_hdr_start / 1e6, xmax = og_hdr_end / 1e6, ymin = -Inf, ymax = Inf), fill = "gray", alpha = 0.1) +
#   geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6,  color = longest_contig), size = 1) +
#   facet_wrap(~chrom) +
#   theme(
#     panel.border = element_rect(color = 'black', fill = NA),
#     panel.background = element_blank()
#   ) +
#   labs(x = "N2 genome position (Mb)", y = "WS contig position (Mb)", title = paste0(hdr_aln$strain))
# 
# 
# # Example of farrrrrrr muli-alignments of the same contig:
# jump <- nucmer_longest %>% dplyr::filter(chrom == "IV" & og_hdr_start == "14170000")
# 
# ggplot(data = jump) +
#   geom_rect(data = jump %>% dplyr::distinct(og_hdr_start, .keep_all = T), aes(xmin = og_hdr_start / 1e6, xmax = og_hdr_end / 1e6, ymin = -Inf, ymax = Inf), fill = "#DB6333", alpha = 0.3) +
#   geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6,  color = longest_contig), size = 1) +
#   facet_wrap(~chrom) +
#   theme(
#     panel.border = element_rect(color = 'black', fill = NA),
#     panel.background = element_blank()
#   ) +
#   labs(x = "N2 genome position (Mb)", y = "WS contig position (Mb)", title = paste0(hdr_aln$strain))
# 
# # I need to remove very distanct alignments from the same contig
# nucmer_longest_jumpRemoved <- nucmer_longest %>%
#   dplyr::mutate(St2=ifelse(inv==T,WSE,WSS),Et2=ifelse(inv==T, WSS, WSE)) %>%
#   dplyr::arrange(St2) %>%
#   dplyr::group_by(chrom, og_hdr_start) %>%
#   dplyr::mutate(leadDiff=lead(St2)-Et2) %>%
#   dplyr::mutate(jump=ifelse(leadDiff > 1.5E5, 1 ,0)) %>%
#   dplyr::mutate(leadDiff=ifelse(is.na(leadDiff),0,leadDiff)) %>%
#   dplyr::mutate(run_id = cumsum(c(1, head(jump, -1)))) %>%
#   dplyr::ungroup() %>%
#   dplyr::group_by(chrom,og_hdr_start,run_id) %>%
#   dplyr::mutate(gsize=n()) %>%
#   dplyr::mutate(len=abs(Et2-St2)) %>%
#   dplyr::mutate(sumlen=sum(len)) %>%
#   dplyr::ungroup() %>%
#   dplyr::group_by(chrom,og_hdr_start) %>%
#   dplyr::filter(sumlen==max(sumlen)) %>%
#   dplyr::select(-gsize) %>%
#   dplyr::ungroup()
# 
# jump_rm <- nucmer_longest_jumpRemoved %>% dplyr::filter(chrom == "IV" & og_hdr_start == "14170000")
# 
# ggplot(data = jump_rm) +
#   geom_rect(data = jump_rm %>% dplyr::distinct(og_hdr_start, .keep_all = T), aes(xmin = og_hdr_start / 1e6, xmax = og_hdr_end / 1e6, ymin = -Inf, ymax = Inf), fill = "#DB6333", alpha = 0.3) +
#   geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6,  color = longest_contig), size = 1) +
#   facet_wrap(~chrom) +
#   theme(
#     panel.border = element_rect(color = 'black', fill = NA),
#     panel.background = element_blank()
#   ) +
#   labs(x = "N2 genome position (Mb)", y = "WS contig position (Mb)", title = paste0(hdr_aln$strain))
# 
# 
# 
# # Calculating the WS coordinates that directly overlap with the N2 HDR boundaries
# nucmer_slope <- nucmer_longest_jumpRemoved %>%
#   # dplyr::group_by(og_hdr_start) %>%
#   # dplyr::mutate(longest_L2 = max(L2)) %>% # filtering for the longest alignment for each HDR to calculate the slope
#   # dplyr::ungroup() %>%
#   dplyr::mutate(slope = ((WSE - WSS) / (N2E - N2S))) %>%
#   dplyr::mutate(intercept = WSS - (slope * N2S)) %>% # to find the y-intercept using point-slope form (b = y1 - mx1)
#   dplyr::mutate(WS_hdr_start = ((slope * og_hdr_start) + intercept)) %>%
#   dplyr::mutate(WS_hdr_end = ((slope * og_hdr_end) + intercept)) %>%
#   dplyr::mutate(spans_hdr = (N2S <= og_hdr_start & N2E >= og_hdr_start) | (N2S <= og_hdr_end & N2E >= og_hdr_end)) %>%
#   dplyr::group_by(chrom, og_hdr_start, spans_hdr) %>%
#   # Add resolution on if there is only an alignment for one edge of the HDR
#   dplyr::mutate(
#     WS_hdr_start_min = ifelse(WSS == min(WSS) & inv == F & spans_hdr == T, WS_hdr_start, 
#                               ifelse(WSE == min(WSE) & inv == T & spans_hdr == T, WS_hdr_end, NA))) %>% # conditional based on if alignment is INV and min WSS! - what if there is only an alignment for one HDR boundary???  # (N2S < og_hdr_start & N2E > og_hdr_start | N2S < og_hdr_end & N2E > og_hdr_end)
#   dplyr::mutate(
#     WS_hdr_end_max = ifelse(WSE == max(WSE) & inv == F & spans_hdr == T, WS_hdr_end, 
#                             ifelse(WSS == max(WSS) & inv == T & spans_hdr == T, WS_hdr_start, NA))) %>% # conditional based on if alignment is INV and max WSE! - what if there is only an alignment for one HDR boundary???
#   dplyr::ungroup() %>%
#   dplyr::group_by(chrom, og_hdr_start) %>%
#   tidyr::fill(WS_hdr_start_min, .direction = "updown") %>%
#   tidyr::fill(WS_hdr_end_max, .direction = "updown") %>%
#   dplyr::ungroup() 
# 
# 
# large <- nucmer_slope %>% dplyr::filter(chrom == "IV" & og_hdr_start == "14170000")
# 
# ggplot(data = large) +
#   geom_rect(data = large %>% dplyr::distinct(og_hdr_start, .keep_all = T), aes(xmin = og_hdr_start / 1e6, xmax = og_hdr_end / 1e6, ymin = -Inf, ymax = Inf), fill = "#DB6333", alpha = 0.3) +
#   geom_rect(data = large %>% dplyr::distinct(WS_hdr_start_min, .keep_all = T), aes(xmin = -Inf, xmax = Inf, ymin = WS_hdr_start_min / 1e6, ymax = WS_hdr_end_max / 1e6), fill = 'blue', alpha = 0.3) +
#   geom_segment(aes(x = og_hdr_start / 1e6, xend = og_hdr_start / 1e6, y = -Inf, yend = Inf), color = 'black') +
#   geom_segment(aes(x = og_hdr_end / 1e6, xend = og_hdr_end / 1e6, y = -Inf, yend = Inf), color = 'black') +
#   geom_segment(aes(x = -Inf, xend = Inf, y = WS_hdr_start_min / 1e6, yend = WS_hdr_start_min / 1e6), color = 'black') +
#   geom_segment(aes(x = -Inf, xend = Inf, y = WS_hdr_end_max / 1e6, yend = WS_hdr_end_max / 1e6), color = 'black') +
#   geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6,  color = spans_hdr), size = 1) +
#   scale_color_manual(values = c("TRUE" = "green", "FALSE" = "blue")) +
#   # coord_cartesian(ylim = c(0,5)) +
#   facet_wrap(~chrom) +
#   theme(
#     panel.border = element_rect(color = 'black', fill = NA),
#     panel.background = element_blank()
#   ) +
#   labs(x = "N2 genome position (Mb)", y = "WS contig position (Mb)", title = paste0(hdr_aln$strain))
# 
# # After adding the condition that if an HDR has only one boundary with an alignment, and the other has an alignment withing X distance of the edge of the other alignment, then use the other alignment for calculating the slope
# nucmer_slope_final <- nucmer_slope %>%
#   dplyr::mutate(spans_hdr_start = (N2S <= og_hdr_start & N2E >= og_hdr_start),
#                 spans_hdr_end = (N2S <= og_hdr_end & N2E >= og_hdr_end)) %>%
#   dplyr::mutate(outside_proximity_left = ifelse(spans_hdr_start == FALSE & N2E < og_hdr_start, og_hdr_start - N2E, NA)) %>%
#   dplyr::mutate(outside_proximity_right = ifelse(spans_hdr_end == FALSE & N2S > og_hdr_end, N2S - og_hdr_end, NA)) %>%
#   dplyr::group_by(chrom, og_hdr_start, og_hdr_end) %>%
#   # dplyr::mutate(WS_hdr_start_min_updated = ifelse(inv == T & WS_hdr_end != WS_hdr_start_min & !is.na(outside_proximity_right) & all(spans_hdr_end == F), WS_hdr_end, 
#   #                                                 ifelse(inv == F & WS_hdr_start != WS_hdr_start_min & !is.na(outside_proximity_left) & all(spans_hdr_start == F), WS_hdr_start, NA))) %>%
#   # dplyr::mutate(WS_hdr_end_max_updated = ifelse(inv == T & WS_hdr_start != WS_hdr_end_max & !is.na(outside_proximity_left) & all(spans_hdr_start == F), WS_hdr_start, 
#   #                                                 ifelse(inv == F & WS_hdr_end != WS_hdr_end_max & !is.na(outside_proximity_right) & all(spans_hdr_end == F), WS_hdr_end, NA))) %>%
#   mutate(inside_hdr = !spans_hdr_start & !spans_hdr_end & is.na(outside_proximity_left) & is.na(outside_proximity_right)) %>%
#   dplyr::ungroup() %>%
#   dplyr::group_by(chrom, og_hdr_start, og_hdr_end, inside_hdr) %>%
#   dplyr::mutate(inside_hdr_max = ifelse(inside_hdr == T & Et2 == max(Et2), Et2, NA),
#                 inside_hdr_min = ifelse(inside_hdr == T & St2 == min(St2), St2, NA)) %>%
#   dplyr::mutate(outside_hdr_noSpan_right = ifelse(inside_hdr == F & all(spans_hdr_end == F) & !is.na(outside_proximity_right), T, F),
#                 outside_hdr_noSpan_left = ifelse(inside_hdr == F & all(spans_hdr_start == F) & !is.na(outside_proximity_left), T, F)) %>%
#   dplyr::ungroup() %>%
#   dplyr::group_by(chrom, og_hdr_start, og_hdr_end, inside_hdr, outside_hdr_noSpan_right) %>%
#   dplyr::mutate(outside_hdr_max_right = ifelse(inside_hdr == F & outside_hdr_noSpan_right == T & Et2 == max(Et2), Et2, NA),
#                 outside_hdr_min_right = ifelse(inside_hdr == F & outside_hdr_noSpan_right == T & St2 == min(St2), St2, NA)) %>%
#   dplyr::ungroup() %>%
#   dplyr::group_by(chrom, og_hdr_start, og_hdr_end, inside_hdr, outside_hdr_noSpan_left) %>%
#   dplyr::mutate(outside_hdr_max_left = ifelse(inside_hdr == F & outside_hdr_noSpan_left == T & Et2 == max(Et2), Et2, NA),
#                 outside_hdr_min_left = ifelse(inside_hdr == F & outside_hdr_noSpan_left == T & St2 == min(St2), St2, NA)) %>%
#   dplyr::ungroup() %>%
#   dplyr::group_by(chrom, og_hdr_start, og_hdr_end, inv) %>%
#   dplyr::mutate(len_inv = ifelse(inv == T, sum(L2), NA)) %>%
#   dplyr::mutate(len_noninv = ifelse(inv == F, sum(L2), NA)) %>%
#   dplyr::mutate(len_inv = ifelse(is.na(len_inv),0,len_inv)) %>%
#   dplyr::mutate(len_noninv = ifelse(is.na(len_noninv),0,len_noninv)) %>%
#   dplyr::ungroup() %>%
#   dplyr::group_by(chrom, og_hdr_start, og_hdr_end) %>%
#   tidyr::fill(len_inv, len_noninv, .direction = "updown") %>%
#   dplyr::mutate(mostly_inv = ifelse(len_inv > len_noninv, T, F)) %>%
#   # dplyr::group_by(chrom, og_hdr_start, og_hdr_end) %>%
#   # dplyr::mutate(
#   #   len_inv    = sum(L2[inv == TRUE],  na.rm = TRUE),
#   #   len_noninv = sum(L2[inv == FALSE], na.rm = TRUE),
#   #   mostly_inv = len_inv > len_noninv) %>%
#   # dplyr::ungroup()
#   tidyr::fill(inside_hdr_min, inside_hdr_max, .direction = "updown") %>%
#   dplyr::mutate(WS_hdr_start_min_updated = ifelse(!is.na(outside_hdr_min_left) & outside_hdr_max_right < outside_hdr_min_left & !(inside_hdr_min < outside_hdr_max_right) & mostly_inv == T, outside_hdr_max_right,  # this works correctly for inverted alignments only
#                                                   ifelse(!is.na(outside_hdr_min_left) & outside_hdr_max_right < outside_hdr_min_left & inside_hdr_min < outside_hdr_max_right & mostly_inv == T, inside_hdr_min, 
#                                                          ifelse(is.na(outside_hdr_min_left) & !(inside_hdr_min < outside_hdr_max_right) & all(spans_hdr_end == F) & mostly_inv == T, outside_hdr_max_right, 
#                                                                 ifelse(is.na(outside_hdr_min_left) & inside_hdr_min < outside_hdr_max_right & all(spans_hdr_end == F) & mostly_inv == T, inside_hdr_min, 
#                                                                        ifelse(is.na(outside_hdr_min_left) & inside_hdr_min < outside_hdr_min_right & all(spans_hdr == F) & mostly_inv == F, inside_hdr_min, # non-INV alignment when there is only an outside alignment on the right
#                                                                               ifelse(!is.na(outside_hdr_min_right) & outside_hdr_max_left < outside_hdr_min_right & !(inside_hdr_min < outside_hdr_max_left) & mostly_inv == F, outside_hdr_max_left, # this works correctly for non-INV alignments only
#                                                                                      ifelse(!is.na(outside_hdr_min_right) & outside_hdr_max_left < outside_hdr_min_right & inside_hdr_min < outside_hdr_max_left & mostly_inv == F, inside_hdr_min, 
#                                                                                             ifelse(is.na(outside_hdr_min_right) & !(inside_hdr_min < outside_hdr_max_left) & all(spans_hdr_start == F) & mostly_inv == F, outside_hdr_max_left, 
#                                                                                                    ifelse(is.na(outside_hdr_min_right) & inside_hdr_min < outside_hdr_max_left & all(spans_hdr_start == F) & mostly_inv == F, inside_hdr_min, 
#                                                                                                           ifelse(is.na(outside_hdr_min_right) & inside_hdr_min < outside_hdr_min_left & all(spans_hdr == F) & mostly_inv == T, inside_hdr_min, # for an inverted alignment with only an outside alignment on the left
#                                                                                                                  ifelse(all(is.na(outside_hdr_min_right)) & all(is.na(outside_hdr_max_left)) & all(spans_hdr == F), inside_hdr_min, NA)))))))))))) %>% # this is agnostic to inv or non-inv
#   dplyr::mutate(WS_hdr_end_max_updated = ifelse(!is.na(outside_hdr_max_right) & outside_hdr_min_left > outside_hdr_max_right & !(inside_hdr_max > outside_hdr_min_left) & mostly_inv == T, outside_hdr_min_left,  # this works correctly for inverted alignments only
#                                                 ifelse(!is.na(outside_hdr_max_right) & outside_hdr_min_left > outside_hdr_max_right & inside_hdr_max > outside_hdr_min_left & mostly_inv == T, inside_hdr_max, 
#                                                        ifelse(is.na(outside_hdr_max_right) & !(inside_hdr_max > outside_hdr_min_left) & all(spans_hdr_start == F) & mostly_inv == T, outside_hdr_min_left, 
#                                                               ifelse(is.na(outside_hdr_max_right) & inside_hdr_max > outside_hdr_min_left & all(spans_hdr_start == F) & mostly_inv == T, inside_hdr_max,
#                                                                      ifelse(is.na(outside_hdr_max_right) & inside_hdr_max > outside_hdr_max_left & all(spans_hdr == F) & mostly_inv == F, inside_hdr_max, # for a non-INV alignment with only an outside alignment of the left
#                                                                             ifelse(!is.na(outside_hdr_max_left) & outside_hdr_min_right > outside_hdr_max_left & !(inside_hdr_max > outside_hdr_min_right) & mostly_inv == F, outside_hdr_min_right, # this works correctly for non-INV alignments only
#                                                                                    ifelse(!is.na(outside_hdr_max_left) & outside_hdr_min_right > outside_hdr_max_left & inside_hdr_max > outside_hdr_min_right & mostly_inv == F, inside_hdr_max,
#                                                                                           ifelse(is.na(outside_hdr_max_left) & !(inside_hdr_max > outside_hdr_min_right) & all(spans_hdr_end == F) & mostly_inv == F, outside_hdr_min_right, 
#                                                                                                  ifelse(is.na(outside_hdr_max_left) & inside_hdr_max > outside_hdr_min_right & all(spans_hdr_end == F) & mostly_inv == F, inside_hdr_max,
#                                                                                                         ifelse(is.na(outside_hdr_max_left) & inside_hdr_max > outside_hdr_max_right & all(spans_hdr == F) & mostly_inv == T, inside_hdr_max, # for inverted alignments with only an outside alignment on the right
#                                                                                                                ifelse(all(is.na(outside_hdr_max_left)) & all(is.na(outside_hdr_max_right)) & all(spans_hdr == F), inside_hdr_max, NA)))))))))))) %>% # this is agnostic to inv or non-inv
#   dplyr::mutate(WS_hdr_start_min_updated = ifelse(spans_hdr_end == T & inside_hdr_min < WS_hdr_start_min & inv == T, inside_hdr_min, # when there is an alignment lower inside than what spans the HDR boundary
#                                                   ifelse(spans_hdr_start == T & inside_hdr_min < WS_hdr_start_min & inv == F, inside_hdr_min,
#                                                          ifelse(spans_hdr_end == T & inside_hdr_min == min(WSS) & inv == F & is.na(outside_hdr_max_right) & is.na(outside_hdr_min_left), inside_hdr_min,
#                                                                 # ifelse(spans_hdr_end == T & inside_hdr_min < WS_hdr_start_min & inv == F & is.na(outside_hdr_max_right) & is.na(outside_hdr_min_left), inside_hdr_min, 
#                                                                 ifelse(spans_hdr_start == T & inside_hdr_min == min(WSE) & inv == T & is.na(outside_hdr_max_right) & is.na(outside_hdr_min_left), inside_hdr_min, WS_hdr_start_min_updated)))),
#                 # ifelse(spans_hdr_start == T & inside_hdr_min < WS_hdr_start_min & inv == T & is.na(outside_hdr_max_right) & is.na(outside_hdr_min_left), inside_hdr_min, WS_hdr_start_min_updated)))),
#                 WS_hdr_end_max_updated = ifelse(spans_hdr_end == T & inside_hdr_max > WS_hdr_end_max & inv == F, inside_hdr_max,
#                                                 ifelse(spans_hdr_start == T & inside_hdr_max > WS_hdr_end_max & inv == T, inside_hdr_max, 
#                                                        # ifelse(spans_hdr_start == T & inside_hdr_max > WS_hdr_end_max & inv == F & is.na(outside_hdr_max_right) & is.na(outside_hdr_min_left), inside_hdr_max, 
#                                                        ifelse(spans_hdr_start == T & inside_hdr_max == max(WSE) & inv == F & is.na(outside_hdr_max_right) & is.na(outside_hdr_min_left), inside_hdr_max,
#                                                               ifelse(spans_hdr_end == T & inside_hdr_max == max(WSS) & inv == T & is.na(outside_hdr_max_right) & is.na(outside_hdr_min_left), inside_hdr_max, WS_hdr_end_max_updated))))) %>%
#   tidyr::fill(WS_hdr_start_min_updated, WS_hdr_end_max_updated, .direction = "downup") %>%
#   dplyr::mutate(WS_hdr_start_min_updated = ifelse(is.na(WS_hdr_start_min_updated),WS_hdr_start_min,WS_hdr_start_min_updated),
#                 WS_hdr_end_max_updated = ifelse(is.na(WS_hdr_end_max_updated),WS_hdr_end_max,WS_hdr_end_max_updated)) %>%
#   dplyr::ungroup()
# 
# # TEST <- nucmer_slope_oneside %>% dplyr::filter(chrom == "V" & og_hdr_start == "6430000")
# 
# ################################################ ^^^^^^^^^^^^^^^^^^^^^^^^^^^ NEED TO FIX FOR test_big SITUATION ###############################################################
# # NEED inside_hdr_min TO BE WS_hdr_start_min_updated
# 
# 
# # Account for situations where there are more than one alignment outside of a boundary for a single HDR!!
# ggplot(data = nucmer_slope_final) + 
#   geom_histogram(aes(x = outside_proximity_left)) +
#   theme_classic()
# 
# ggplot(data = nucmer_slope_final) + 
#   geom_histogram(aes(x = outside_proximity_right)) +
#   theme_classic()
# 
# 
# 
# large_updated <- nucmer_slope_final %>% dplyr::filter(chrom == "IV" & og_hdr_start == "14170000") 
# 
# ggplot(data = large_updated) +
#   geom_rect(data = large_updated %>% dplyr::distinct(og_hdr_start, .keep_all = T), aes(xmin = og_hdr_start / 1e6, xmax = og_hdr_end / 1e6, ymin = -Inf, ymax = Inf), fill = "#DB6333", alpha = 0.3) +
#   geom_rect(data = large_updated %>% dplyr::distinct(WS_hdr_start_min_updated, .keep_all = T), aes(xmin = -Inf, xmax = Inf, ymin = WS_hdr_start_min_updated / 1e6, ymax = WS_hdr_end_max / 1e6), fill = 'blue', alpha = 0.3) +
#   geom_segment(aes(x = og_hdr_start / 1e6, xend = og_hdr_start / 1e6, y = -Inf, yend = Inf), color = 'black') +
#   geom_segment(aes(x = og_hdr_end / 1e6, xend = og_hdr_end / 1e6, y = -Inf, yend = Inf), color = 'black') +
#   
#   geom_segment(aes(x = -Inf, xend = Inf, y = outside_hdr_min_left / 1e6, yend = outside_hdr_min_left / 1e6), color = 'gold3', size = 2) +
#   geom_segment(aes(x = -Inf, xend = Inf, y = outside_hdr_max_left / 1e6, yend = outside_hdr_max_left / 1e6), color = 'gold3', size = 2) +
#   geom_segment(aes(x = -Inf, xend = Inf, y = outside_hdr_min_right / 1e6, yend = outside_hdr_min_right / 1e6), color = 'purple', size = 2) +
#   geom_segment(aes(x = -Inf, xend = Inf, y = outside_hdr_max_right / 1e6, yend = outside_hdr_max_right / 1e6), color = 'purple', size = 2) +
#   geom_segment(aes(x = -Inf, xend = Inf, y = inside_hdr_max / 1e6, yend = inside_hdr_max / 1e6), color = 'red', size = 1) +
#   geom_segment(aes(x = -Inf, xend = Inf, y = inside_hdr_min / 1e6, yend = inside_hdr_min / 1e6), color = 'red', size = 1) +
#   
#   geom_segment(aes(x = -Inf, xend = Inf, y = WS_hdr_start_min_updated / 1e6, yend = WS_hdr_start_min_updated / 1e6), color = 'black') +
#   geom_segment(aes(x = -Inf, xend = Inf, y = WS_hdr_end_max / 1e6, yend = WS_hdr_end_max / 1e6), color = 'black') +
#   geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6,  color = spans_hdr), size = 1) +
#   scale_color_manual(values = c("TRUE" = "green", "FALSE" = "blue")) +
#   # coord_cartesian(ylim = c(0,5)) +
#   facet_wrap(~chrom) +
#   theme(
#     panel.border = element_rect(color = 'black', fill = NA),
#     panel.background = element_blank()
#   ) +
#   labs(x = "N2 genome position (Mb)", y = "WS contig position (Mb)", title = paste0(hdr_aln$strain))
# 
# 
# # Now the other side!
# test2 <- nucmer_slope_final %>% dplyr::filter(inv == F & !is.na(outside_proximity_right))
# 
# ugh <- nucmer_slope_final %>% dplyr::filter(chrom == "X" & og_hdr_start == "14355000")
# ugh2 <- nucmer_slope %>% dplyr::filter(chrom == "X" & og_hdr_start == "14355000")
# 
# # ggplot(data = ugh2) +
# #   geom_rect(data = ugh2 %>% dplyr::distinct(og_hdr_start, .keep_all = T), aes(xmin = og_hdr_start / 1e6, xmax = og_hdr_end / 1e6, ymin = -Inf, ymax = Inf), fill = "#DB6333", alpha = 0.3) +
# #   geom_rect(data = ugh2 %>% dplyr::distinct(WS_hdr_start_min, .keep_all = T), aes(xmin = -Inf, xmax = Inf, ymin = WS_hdr_start_min / 1e6, ymax = WS_hdr_end_max / 1e6), fill = 'blue', alpha = 0.3) +
# #   geom_segment(aes(x = og_hdr_start / 1e6, xend = og_hdr_start / 1e6, y = -Inf, yend = Inf), color = 'black') +
# #   geom_segment(aes(x = og_hdr_end / 1e6, xend = og_hdr_end / 1e6, y = -Inf, yend = Inf), color = 'black') +
# #   geom_segment(aes(x = -Inf, xend = Inf, y = WS_hdr_start_min / 1e6, yend = WS_hdr_start_min / 1e6), color = 'black') +
# #   geom_segment(aes(x = -Inf, xend = Inf, y = WS_hdr_end_max / 1e6, yend = WS_hdr_end_max / 1e6), color = 'black') +
# #   geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6,  color = spans_hdr), size = 1) +
# #   scale_color_manual(values = c("TRUE" = "green", "FALSE" = "blue")) +
# #   # coord_cartesian(ylim = c(0,5)) +
# #   facet_wrap(~chrom) +
# #   theme(
# #     panel.border = element_rect(color = 'black', fill = NA),
# #     panel.background = element_blank()
# #   ) +
# #   labs(x = "N2 genome position (Mb)", y = "WS contig position (Mb)", title = paste0(hdr_aln$strain))
# 
# ggplot(data = ugh) +
#   geom_rect(data = ugh %>% dplyr::distinct(og_hdr_start, .keep_all = T), aes(xmin = og_hdr_start / 1e6, xmax = og_hdr_end / 1e6, ymin = -Inf, ymax = Inf), fill = "#DB6333", alpha = 0.3) +
#   geom_rect(data = ugh %>% dplyr::distinct(WS_hdr_start_min_updated, .keep_all = T), aes(xmin = -Inf, xmax = Inf, ymin = WS_hdr_start_min_updated / 1e6, ymax = WS_hdr_end_max_updated / 1e6), fill = 'blue', alpha = 0.3) +
#   geom_segment(aes(x = og_hdr_start / 1e6, xend = og_hdr_start / 1e6, y = -Inf, yend = Inf), color = 'black') +
#   
#   geom_segment(aes(x = -Inf, xend = Inf, y = outside_hdr_min_left / 1e6, yend = outside_hdr_min_left / 1e6), color = 'gold3', size = 2) +
#   geom_segment(aes(x = -Inf, xend = Inf, y = outside_hdr_max_left / 1e6, yend = outside_hdr_max_left / 1e6), color = 'gold3', size = 2) +
#   geom_segment(aes(x = -Inf, xend = Inf, y = outside_hdr_min_right / 1e6, yend = outside_hdr_min_right / 1e6), color = 'purple', size = 2) +
#   geom_segment(aes(x = -Inf, xend = Inf, y = outside_hdr_max_right / 1e6, yend = outside_hdr_max_right / 1e6), color = 'purple', size = 2) +
#   geom_segment(aes(x = -Inf, xend = Inf, y = inside_hdr_max / 1e6, yend = inside_hdr_max / 1e6), color = 'red', size = 1) +
#   geom_segment(aes(x = -Inf, xend = Inf, y = inside_hdr_min / 1e6, yend = inside_hdr_min / 1e6), color = 'red', size = 1) +
#   
#   geom_segment(aes(x = og_hdr_end / 1e6, xend = og_hdr_end / 1e6, y = -Inf, yend = Inf), color = 'black') +
#   geom_segment(aes(x = -Inf, xend = Inf, y = WS_hdr_start_min_updated / 1e6, yend = WS_hdr_start_min_updated / 1e6), color = 'black') +
#   geom_segment(aes(x = -Inf, xend = Inf, y = WS_hdr_end_max_updated / 1e6, yend = WS_hdr_end_max_updated / 1e6), color = 'black') +
#   geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6,  color = spans_hdr), size = 1) +
#   scale_color_manual(values = c("TRUE" = "green", "FALSE" = "blue")) +
#   # coord_cartesian(ylim = c(0,5)) +
#   facet_wrap(~chrom) +
#   theme(
#     panel.border = element_rect(color = 'black', fill = NA),
#     panel.background = element_blank()
#   ) +
#   labs(x = "N2 genome position (Mb)", y = "WS contig position (Mb)", title = paste0(hdr_aln$strain))
# 
# 
# # Plot HDRs that have no alignment spanning the edges!
# # no_aln_boundaries <- nucmer_slope_oneside %>% 
# #   dplyr::group_by(chrom, og_hdr_start, og_hdr_end) %>%
# #   dplyr::filter(all(spans_hdr == F))
# 
# ggplot(nucmer_slope_final %>% dplyr::filter(chrom == "V" & og_hdr_start == "19480000")) + 
#   geom_rect(data = nucmer_slope_final %>% dplyr::filter(chrom == "V" & og_hdr_start == "19480000") %>% dplyr::distinct(og_hdr_start, .keep_all = T), aes(xmin = og_hdr_start / 1e6, xmax = og_hdr_end / 1e6, ymin = -Inf, ymax = Inf), fill = "#DB6333", alpha = 0.3) +
#   # geom_rect(data = no_aln_boundaries %>% dplyr::filter(og_hdr_start == "19480000") %>% dplyr::distinct(WS_hdr_start_min, .keep_all = T), aes(xmin = -Inf, xmax = Inf, ymin = WS_hdr_start_min / 1e6, ymax = WS_hdr_end_max / 1e6), fill = 'blue', alpha = 0.3) +
#   geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6,  color = spans_hdr), size = 1) +
#   geom_rect(data = nucmer_slope_final %>% dplyr::filter(chrom == "V" & og_hdr_start == "19480000") %>% dplyr::distinct(WS_hdr_start_min_updated, .keep_all = T), aes(xmin = -Inf, xmax = Inf, ymin = WS_hdr_start_min_updated / 1e6, ymax = WS_hdr_end_max_updated / 1e6), fill = 'blue', alpha = 0.3) +
#   geom_segment(aes(x = og_hdr_start / 1e6, xend = og_hdr_start / 1e6, y = -Inf, yend = Inf), color = 'black') +
#   geom_segment(aes(x = -Inf, xend = Inf, y = outside_hdr_min_left / 1e6, yend = outside_hdr_min_left / 1e6), color = 'gold3', size = 2) +
#   geom_segment(aes(x = -Inf, xend = Inf, y = outside_hdr_max_left / 1e6, yend = outside_hdr_max_left / 1e6), color = 'gold3', size = 2) +
#   geom_segment(aes(x = -Inf, xend = Inf, y = outside_hdr_min_right / 1e6, yend = outside_hdr_min_right / 1e6), color = 'purple', size = 2) +
#   geom_segment(aes(x = -Inf, xend = Inf, y = outside_hdr_max_right / 1e6, yend = outside_hdr_max_right / 1e6), color = 'purple', size = 2) +
#   geom_segment(aes(x = -Inf, xend = Inf, y = inside_hdr_max / 1e6, yend = inside_hdr_max / 1e6), color = 'red', size = 1) +
#   geom_segment(aes(x = -Inf, xend = Inf, y = inside_hdr_min / 1e6, yend = inside_hdr_min / 1e6), color = 'red', size = 1) +
#   
#   
#   scale_color_manual(values = c("TRUE" = "green", "FALSE" = "blue")) +
#   facet_wrap(~chrom, scales = "free") +
#   theme(
#     panel.border = element_rect(color = 'black', fill = NA),
#     panel.background = element_blank()
#   ) +
#   labs(x = "N2 genome position (Mb)", y = "WS contig position (Mb)", title = paste0(hdr_aln$strain))
# 
# no_aln_bound_inv <- nucmer_slope_final %>% 
#   dplyr::group_by(chrom, og_hdr_start, og_hdr_end) %>%
#   dplyr::filter(all(spans_hdr == F)) %>%
#   dplyr::filter(chrom == "V" & og_hdr_start == "17045000")
# 
# ggplot(no_aln_bound_inv) + 
#   geom_rect(data = no_aln_bound_inv %>% dplyr::distinct(og_hdr_start, .keep_all = T), aes(xmin = og_hdr_start / 1e6, xmax = og_hdr_end / 1e6, ymin = -Inf, ymax = Inf), fill = "#DB6333", alpha = 0.3) +
#   # geom_rect(data = no_aln_bound_inv %>% dplyr::distinct(WS_hdr_start_min, .keep_all = T), aes(xmin = -Inf, xmax = Inf, ymin = WS_hdr_start_min / 1e6, ymax = WS_hdr_end_max / 1e6), fill = 'blue', alpha = 0.3) +
#   geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6,  color = spans_hdr), size = 1) +
#   geom_rect(data = no_aln_bound_inv %>% dplyr::distinct(WS_hdr_start_min_updated, .keep_all = T), aes(xmin = -Inf, xmax = Inf, ymin = WS_hdr_start_min_updated / 1e6, ymax = WS_hdr_end_max_updated / 1e6), fill = 'blue', alpha = 0.3) +
#   geom_segment(aes(x = og_hdr_start / 1e6, xend = og_hdr_start / 1e6, y = -Inf, yend = Inf), color = 'black') +
#   
#   geom_segment(aes(x = -Inf, xend = Inf, y = outside_hdr_min_left / 1e6, yend = outside_hdr_min_left / 1e6), color = 'magenta3', size = 2) +
#   geom_segment(aes(x = -Inf, xend = Inf, y = outside_hdr_max_left / 1e6, yend = outside_hdr_max_left / 1e6), color = 'magenta3', size = 2) +
#   geom_segment(aes(x = -Inf, xend = Inf, y = outside_hdr_min_right / 1e6, yend = outside_hdr_min_right / 1e6), color = 'purple', size = 2) +
#   geom_segment(aes(x = -Inf, xend = Inf, y = outside_hdr_max_right / 1e6, yend = outside_hdr_max_right / 1e6), color = 'purple', size = 2) +
#   geom_segment(aes(x = -Inf, xend = Inf, y = inside_hdr_max / 1e6, yend = inside_hdr_max / 1e6), color = 'red', size = 1) +
#   geom_segment(aes(x = -Inf, xend = Inf, y = inside_hdr_min / 1e6, yend = inside_hdr_min / 1e6), color = 'red', size = 1) +
#   scale_color_manual(values = c("TRUE" = "green", "FALSE" = "blue")) +
#   facet_wrap(~chrom, scales = "free") +
#   theme(
#     panel.border = element_rect(color = 'black', fill = NA),
#     panel.background = element_blank()
#   ) +
#   labs(x = "N2 genome position (Mb)", y = "WS contig position (Mb)", title = paste0(hdr_aln$strain))
# 
# 
# test_big <- nucmer_slope_final %>% dplyr::filter(chrom == "V" & og_hdr_start == "18433000")
# ope <- test_big %>% dplyr::filter(WSE == "2323592")
# op2 <- test_big %>% dplyr::filter(WSS == "1362607")
# 
# ggplot(test_big) + 
#   geom_rect(data = test_big %>% dplyr::distinct(og_hdr_start, .keep_all = T), aes(xmin = og_hdr_start / 1e6, xmax = og_hdr_end / 1e6, ymin = -Inf, ymax = Inf), fill = "#DB6333", alpha = 0.3) +
#   # geom_rect(data = no_aln_bound_inv %>% dplyr::distinct(WS_hdr_start_min, .keep_all = T), aes(xmin = -Inf, xmax = Inf, ymin = WS_hdr_start_min / 1e6, ymax = WS_hdr_end_max / 1e6), fill = 'blue', alpha = 0.3) +
#   geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6,  color = spans_hdr), size = 1) +
#   geom_segment(data = ope, aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6), color = 'cyan', size = 12) + # LOWER BOUNDARY IS BEING SET BASED ON PREDICTED END (START, BUT INVERTED) OF WS HDR FROM SLOPE... (WS_hdr_start_min) 
#   geom_rect(data = test_big %>% dplyr::distinct(WS_hdr_start_min_updated, .keep_all = T), aes(xmin = -Inf, xmax = Inf, ymin = WS_hdr_start_min_updated / 1e6, ymax = WS_hdr_end_max_updated / 1e6), fill = 'blue', alpha = 0.3) +
#   geom_segment(aes(x = og_hdr_start / 1e6, xend = og_hdr_start / 1e6, y = -Inf, yend = Inf), color = 'black') +
#   
#   geom_segment(aes(x = -Inf, xend = Inf, y = outside_hdr_min_left / 1e6, yend = outside_hdr_min_left / 1e6), color = 'magenta3', size = 2) +
#   geom_segment(aes(x = -Inf, xend = Inf, y = outside_hdr_max_left / 1e6, yend = outside_hdr_max_left / 1e6), color = 'magenta3', size = 2) +
#   geom_segment(aes(x = -Inf, xend = Inf, y = outside_hdr_min_right / 1e6, yend = outside_hdr_min_right / 1e6), color = 'purple', size = 2) +
#   geom_segment(aes(x = -Inf, xend = Inf, y = outside_hdr_max_right / 1e6, yend = outside_hdr_max_right / 1e6), color = 'purple', size = 2) +
#   geom_segment(aes(x = -Inf, xend = Inf, y = inside_hdr_max / 1e6, yend = inside_hdr_max / 1e6), color = 'red', size = 1) +
#   geom_segment(aes(x = -Inf, xend = Inf, y = inside_hdr_min / 1e6, yend = inside_hdr_min / 1e6), color = 'red', size = 1) +
#   scale_color_manual(values = c("TRUE" = "green", "FALSE" = "blue")) +
#   facet_wrap(~chrom, scales = "free") +
#   theme(
#     panel.border = element_rect(color = 'black', fill = NA),
#     panel.background = element_blank()
#   ) +
#   labs(x = "N2 genome position (Mb)", y = "WS contig position (Mb)", title = paste0(hdr_aln$strain))
# 
# 
# 
# # nucmer_slope_final <- nucmer_slope_oneside %>%
# #   dplyr::group_by(chrom, og_hdr_start, og_hdr_end) %>%
# #   dplyr::mutate(WS_hdr_start_min_updated = ifelse(is.na(WS_hdr_start_min_updated) & all(spans_hdr == F) & inv == F & WSS == min(WSS) & !is.na(outside_proximity_left) |
# #                                                     is.na(WS_hdr_start_min_updated) & all(spans_hdr == F) & inv == F & WSE == max(WSE) & !is.na(outside_proximity_right), WS_hdr_start, 
# #                                                          ifelse(is.na(WS_hdr_start_min_updated) & all(spans_hdr == F) & inv == T & WSE == min(WSE) & !is.na(outside_proximity_right) | 
# #                                                                   is.na(WS_hdr_start_min_updated) & all(spans_hdr == F) & inv == T & WSS == max(WSS) & !is.na(outside_proximity_left), WS_hdr_end, WS_hdr_start_min_updated))) %>%
# #   tidyr::fill(WS_hdr_start_min_updated, .direction = "updown") %>%
# #   dplyr::mutate(WS_hdr_end_max_updated = ifelse(is.na(WS_hdr_end_max_updated) & all(spans_hdr == F) & inv == F & WSS == min(WSS) & !is.na(outside_proximity_left) | 
# #                                                   is.na(WS_hdr_end_max_updated) & all(spans_hdr == F) & inv == F & WSE == max(WSE) & !is.na(outside_proximity_right), WS_hdr_end, 
# #                                                        ifelse(is.na(WS_hdr_end_max_updated) & all(spans_hdr == F) & inv == T & WSE == min(WSE) & !is.na(outside_proximity_right) |
# #                                                                 is.na(WS_hdr_end_max_updated) & all(spans_hdr == F) & inv == T & WSS == max(WSS) & !is.na(outside_proximity_left), WS_hdr_start, WS_hdr_end_max_updated))) %>%
# #   tidyr::fill(WS_hdr_end_max_updated, .direction = "updown") %>%
# #   dplyr::ungroup()
# 
# 
# ggplot(nucmer_slope_final %>% dplyr::filter(og_hdr_start == "19480000")) + 
#   geom_rect(data = nucmer_slope_final %>% dplyr::filter(og_hdr_start == "19480000") %>% dplyr::distinct(og_hdr_start, .keep_all = T), aes(xmin = og_hdr_start / 1e6, xmax = og_hdr_end / 1e6, ymin = -Inf, ymax = Inf), fill = "#DB6333", alpha = 0.3) +
#   geom_rect(data = nucmer_slope_final %>% dplyr::filter(og_hdr_start == "19480000") %>% dplyr::distinct(WS_hdr_start_min_updated, .keep_all = T), aes(xmin = -Inf, xmax = Inf, ymin = WS_hdr_start_min_updated / 1e6, ymax = WS_hdr_end_max_updated / 1e6), fill = 'blue', alpha = 0.3) +
#   geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6,  color = spans_hdr), size = 1) +
#   scale_color_manual(values = c("TRUE" = "green", "FALSE" = "blue")) +
#   facet_wrap(~chrom, scales = "free") +
#   theme(
#     panel.border = element_rect(color = 'black', fill = NA),
#     panel.background = element_blank()
#   ) +
#   labs(x = "N2 genome position (Mb)", y = "WS contig position (Mb)", title = paste0(hdr_aln$strain))
# 
# 
# ggplot(nucmer_slope_final %>% dplyr::filter(chrom == "V" & og_hdr_start == "17045000")) + 
#   geom_rect(data = nucmer_slope_final %>%   dplyr::filter(chrom == "V" & og_hdr_start == "17045000") %>% dplyr::distinct(og_hdr_start, .keep_all = T), aes(xmin = og_hdr_start / 1e6, xmax = og_hdr_end / 1e6, ymin = -Inf, ymax = Inf), fill = "#DB6333", alpha = 0.3) +
#   geom_rect(data = nucmer_slope_final %>%   dplyr::filter(chrom == "V" & og_hdr_start == "17045000") %>% dplyr::distinct(WS_hdr_start_min_updated, .keep_all = T), aes(xmin = -Inf, xmax = Inf, ymin = WS_hdr_start_min_updated / 1e6, ymax = WS_hdr_end_max_updated / 1e6), fill = 'blue', alpha = 0.3) +
#   geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6,  color = spans_hdr), size = 1) +
#   scale_color_manual(values = c("TRUE" = "green", "FALSE" = "blue")) +
#   facet_wrap(~chrom, scales = "free") +
#   theme(
#     panel.border = element_rect(color = 'black', fill = NA),
#     panel.background = element_blank()
#   ) +
#   labs(x = "N2 genome position (Mb)", y = "WS contig position (Mb)", title = paste0(hdr_aln$strain))
# 
# 
# test <- nucmer_slope_final %>% dplyr::filter(chrom == "IV" & og_hdr_start > 12600000 & og_hdr_end < 13000000)
# ggplot(test) + 
#   geom_rect(data = test %>% dplyr::distinct(og_hdr_start, .keep_all = T), aes(xmin = og_hdr_start / 1e6, xmax = og_hdr_end / 1e6, ymin = -Inf, ymax = Inf), fill = "#DB6333", alpha = 0.3) +
#   geom_rect(data = test %>% dplyr::distinct(WS_hdr_start_min_updated, .keep_all = T), aes(xmin = -Inf, xmax = Inf, ymin = WS_hdr_start_min_updated / 1e6, ymax = WS_hdr_end_max_updated / 1e6), fill = 'blue', alpha = 0.3) +
#   geom_segment(aes(x = og_hdr_start / 1e6, xend = og_hdr_start / 1e6, y = -Inf, yend = Inf), color = 'black') +
#   geom_segment(aes(x = og_hdr_end / 1e6, xend = og_hdr_end / 1e6, y = -Inf, yend = Inf), color = 'black') +
#   geom_segment(aes(x = -Inf, xend = Inf, y = WS_hdr_start_min / 1e6, yend = WS_hdr_start_min / 1e6), color = 'black') +
#   geom_segment(aes(x = -Inf, xend = Inf, y = WS_hdr_end_max / 1e6, yend = WS_hdr_end_max / 1e6), color = 'black') +
#   geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6,  color = spans_hdr), size = 1) +
#   scale_color_manual(values = c("TRUE" = "green", "FALSE" = "blue")) +
#   facet_wrap(~chrom) +
#   theme(
#     panel.border = element_rect(color = 'black', fill = NA),
#     panel.background = element_blank()
#   ) +
#   labs(x = "N2 genome position (Mb)", y = "WS contig position (Mb)", title = paste0(hdr_aln$strain))
# 
# 
# test_inv <- nucmer_slope_final %>% dplyr::filter(chrom == "V" & og_hdr_start > 6200000 & og_hdr_end < 6700000)
# ggplot(test_inv) + 
#   geom_rect(data = test_inv %>% dplyr::distinct(og_hdr_start, .keep_all = T), aes(xmin = og_hdr_start / 1e6, xmax = og_hdr_end / 1e6, ymin = -Inf, ymax = Inf), fill = "#DB6333", alpha = 0.3) +
#   geom_rect(data = test_inv %>% dplyr::distinct(WS_hdr_start_min_updated, .keep_all = T), aes(xmin = -Inf, xmax = Inf, ymin = WS_hdr_start_min_updated / 1e6, ymax = WS_hdr_end_max_updated / 1e6), fill = 'blue', alpha = 0.3) +
#   geom_segment(aes(x = og_hdr_start / 1e6, xend = og_hdr_start / 1e6, y = -Inf, yend = Inf), color = 'black') +
#   geom_segment(aes(x = og_hdr_end / 1e6, xend = og_hdr_end / 1e6, y = -Inf, yend = Inf), color = 'black') +
#   geom_segment(aes(x = -Inf, xend = Inf, y = WS_hdr_start_min_updated / 1e6, yend = WS_hdr_start_min_updated / 1e6), color = 'black') +
#   geom_segment(aes(x = -Inf, xend = Inf, y = WS_hdr_end_max_updated / 1e6, yend = WS_hdr_end_max_updated / 1e6), color = 'black') +
#   geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6,  color = spans_hdr), size = 1) +
#   scale_color_manual(values = c("TRUE" = "green", "FALSE" = "blue")) +
#   facet_wrap(~chrom) +
#   theme(
#     panel.border = element_rect(color = 'black', fill = NA),
#     panel.background = element_blank()
#   ) +
#   labs(x = "N2 genome position (Mb)", y = "WS contig position (Mb)", title = paste0(hdr_aln$strain))
# 
# 
# 
# ggplot(nucmer_slope_final) + 
#   geom_rect(data = nucmer_slope_final %>% dplyr::distinct(og_hdr_start, .keep_all = T), aes(xmin = og_hdr_start / 1e6, xmax = og_hdr_end / 1e6, ymin = -Inf, ymax = Inf), fill = "#DB6333", alpha = 0.3) +
#   geom_rect(data = nucmer_slope_final %>% dplyr::distinct(WS_hdr_start_min_updated, .keep_all = T), aes(xmin = -Inf, xmax = Inf, ymin = WS_hdr_start_min_updated / 1e6, ymax = WS_hdr_end_max_updated / 1e6), fill = 'blue', alpha = 0.3) +
#   # geom_segment(aes(x = og_hdr_start / 1e6, xend = og_hdr_start / 1e6, y = -Inf, yend = Inf), color = 'black') +
#   # geom_segment(aes(x = og_hdr_end / 1e6, xend = og_hdr_end / 1e6, y = -Inf, yend = Inf), color = 'black') +
#   # geom_segment(aes(x = -Inf, xend = Inf, y = WS_hdr_start_min / 1e6, yend = WS_hdr_start_min / 1e6), color = 'black') +
#   # geom_segment(aes(x = -Inf, xend = Inf, y = WS_hdr_end_max / 1e6, yend = WS_hdr_end_max / 1e6), color = 'black') +
#   geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6,  color = spans_hdr), size = 1) +
#   scale_color_manual(values = c("TRUE" = "green", "FALSE" = "blue")) +
#   facet_wrap(~chrom, scales = "free") +
#   theme(
#     panel.border = element_rect(color = 'black', fill = NA),
#     panel.background = element_blank()
#   ) +
#   labs(x = "N2 genome position (Mb)", y = "WS contig position (Mb)", title = paste0(hdr_aln$strain))
# 
# ggplot(nucmer_slope_final %>% dplyr::filter(chrom == "I" & N2S < 2000000)) +
#   geom_rect(data = nucmer_slope_final %>% dplyr::filter(chrom == "I" & N2S < 2000000) %>% dplyr::distinct(og_hdr_start, .keep_all = T), aes(xmin = og_hdr_start / 1e6, xmax = og_hdr_end / 1e6, ymin = -Inf, ymax = Inf), fill = "#DB6333", alpha = 0.3) +
#   geom_rect(data = nucmer_slope_final %>% dplyr::filter(chrom == "I" & N2S < 2000000) %>% dplyr::distinct(WS_hdr_start_min_updated, .keep_all = T), aes(xmin = -Inf, xmax = Inf, ymin = WS_hdr_start_min_updated / 1e6, ymax = WS_hdr_end_max_updated / 1e6), fill = 'blue', alpha = 0.3) +
#   # geom_rect(data = FINAL_hdrs %>% dplyr::filter(chrom == "I" & og_hdr_start < 2000000) %>% dplyr::distinct(WS_hdr_start_min_updated, .keep_all = T), aes(xmin = -Inf, xmax = Inf, ymin = WS_hdr_start_min_updated / 1e6, ymax = WS_hdr_end_max_updated / 1e6), fill = 'red', alpha = 0.3) +
#   # geom_segment(aes(x = og_hdr_start / 1e6, xend = og_hdr_start / 1e6, y = -Inf, yend = Inf), color = 'black') +
#   # geom_segment(aes(x = og_hdr_end / 1e6, xend = og_hdr_end / 1e6, y = -Inf, yend = Inf), color = 'black') +
#   # geom_segment(aes(x = -Inf, xend = Inf, y = WS_hdr_start_min / 1e6, yend = WS_hdr_start_min / 1e6), color = 'black') +
#   # geom_segment(aes(x = -Inf, xend = Inf, y = WS_hdr_end_max / 1e6, yend = WS_hdr_end_max / 1e6), color = 'black') +
#   geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6,  color = spans_hdr), size = 1) +
#   scale_color_manual(values = c("TRUE" = "green", "FALSE" = "blue")) +
#   facet_wrap(~chrom, scales = "free") +
#   theme(
#     panel.border = element_rect(color = 'black', fill = NA),
#     panel.background = element_blank()
#   ) +
#   labs(x = "N2 genome position (Mb)", y = "WS contig position (Mb)", title = paste0(hdr_aln$strain))
# 
# ggplot(nucmer_slope_final %>% dplyr::filter(chrom == "III" & N2S > 11000000)) +
#   geom_rect(data = nucmer_slope_final %>% dplyr::filter(chrom == "III" & N2S > 11000000) %>% dplyr::distinct(og_hdr_start, .keep_all = T), aes(xmin = og_hdr_start / 1e6, xmax = og_hdr_end / 1e6, ymin = -Inf, ymax = Inf), fill = "#DB6333", alpha = 0.3) +
#   geom_rect(data = nucmer_slope_final %>% dplyr::filter(chrom == "III" & N2S > 11000000) %>% dplyr::distinct(WS_hdr_start_min_updated, .keep_all = T), aes(xmin = -Inf, xmax = Inf, ymin = WS_hdr_start_min_updated / 1e6, ymax = WS_hdr_end_max_updated / 1e6), fill = 'blue', alpha = 0.3) +
#   # geom_segment(aes(x = og_hdr_start / 1e6, xend = og_hdr_start / 1e6, y = -Inf, yend = Inf), color = 'black') +
#   # geom_segment(aes(x = og_hdr_end / 1e6, xend = og_hdr_end / 1e6, y = -Inf, yend = Inf), color = 'black') +
#   # geom_segment(aes(x = -Inf, xend = Inf, y = WS_hdr_start_min / 1e6, yend = WS_hdr_start_min / 1e6), color = 'black') +
#   # geom_segment(aes(x = -Inf, xend = Inf, y = WS_hdr_end_max / 1e6, yend = WS_hdr_end_max / 1e6), color = 'black') +
#   geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6,  color = spans_hdr), size = 1) +
#   scale_color_manual(values = c("TRUE" = "green", "FALSE" = "blue")) +
#   facet_wrap(~chrom, scales = "free") +
#   theme(
#     panel.border = element_rect(color = 'black', fill = NA),
#     panel.background = element_blank()
#   ) +
#   labs(x = "N2 genome position (Mb)", y = "WS contig position (Mb)", title = paste0(hdr_aln$strain))
# 
# ggplot(nucmer_slope_final %>% dplyr::filter(chrom == "V" & N2S > 6000000 & N2E < 8000000)) +
#   geom_rect(data = nucmer_slope_final %>% dplyr::filter(chrom == "V" & N2S > 6000000 & N2E < 8000000) %>% dplyr::distinct(og_hdr_start, .keep_all = T), aes(xmin = og_hdr_start / 1e6, xmax = og_hdr_end / 1e6, ymin = -Inf, ymax = Inf), fill = "#DB6333", alpha = 0.3) +
#   geom_rect(data = nucmer_slope_final %>% dplyr::filter(chrom == "V" & N2S > 6000000 & N2E < 8000000) %>% dplyr::distinct(WS_hdr_start_min_updated, .keep_all = T), aes(xmin = -Inf, xmax = Inf, ymin = WS_hdr_start_min_updated / 1e6, ymax = WS_hdr_end_max_updated / 1e6), fill = 'blue', alpha = 0.3) +
#   # geom_segment(aes(x = og_hdr_start / 1e6, xend = og_hdr_start / 1e6, y = -Inf, yend = Inf), color = 'black') +
#   # geom_segment(aes(x = og_hdr_end / 1e6, xend = og_hdr_end / 1e6, y = -Inf, yend = Inf), color = 'black') +
#   # geom_segment(aes(x = -Inf, xend = Inf, y = WS_hdr_start_min / 1e6, yend = WS_hdr_start_min / 1e6), color = 'black') +
#   # geom_segment(aes(x = -Inf, xend = Inf, y = WS_hdr_end_max / 1e6, yend = WS_hdr_end_max / 1e6), color = 'black') +
#   geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6,  color = spans_hdr), size = 1) +
#   scale_color_manual(values = c("TRUE" = "green", "FALSE" = "blue")) +
#   facet_wrap(~chrom, scales = "free") +
#   theme(
#     panel.border = element_rect(color = 'black', fill = NA),
#     panel.background = element_blank()
#   ) +
#   labs(x = "N2 genome position (Mb)", y = "WS contig position (Mb)", title = paste0(hdr_aln$strain))
# 
# ughhhh <- nucmer_slope_final %>% dplyr::filter(chrom == "V" & og_hdr_start == 7234000)
# 
# ggplot(ughhhh) +
#   geom_rect(data = ughhhh  %>% dplyr::distinct(og_hdr_start, .keep_all = T), aes(xmin = og_hdr_start / 1e6, xmax = og_hdr_end / 1e6, ymin = -Inf, ymax = Inf), fill = "#DB6333", alpha = 0.3) +
#   geom_rect(data = ughhhh %>% dplyr::distinct(WS_hdr_start_min_updated, .keep_all = T), aes(xmin = -Inf, xmax = Inf, ymin = WS_hdr_start_min_updated / 1e6, ymax = WS_hdr_end_max_updated / 1e6), fill = 'blue', alpha = 0.3) +
#   # geom_segment(aes(x = og_hdr_start / 1e6, xend = og_hdr_start / 1e6, y = -Inf, yend = Inf), color = 'black') +
#   # geom_segment(aes(x = og_hdr_end / 1e6, xend = og_hdr_end / 1e6, y = -Inf, yend = Inf), color = 'black') +
#   # geom_segment(aes(x = -Inf, xend = Inf, y = WS_hdr_start_min / 1e6, yend = WS_hdr_start_min / 1e6), color = 'black') +
#   # geom_segment(aes(x = -Inf, xend = Inf, y = WS_hdr_end_max / 1e6, yend = WS_hdr_end_max / 1e6), color = 'black') +
#   geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6,  color = spans_hdr), size = 1) +
#   scale_color_manual(values = c("TRUE" = "green", "FALSE" = "blue")) +
#   facet_wrap(~chrom, scales = "free") +
#   theme(
#     panel.border = element_rect(color = 'black', fill = NA),
#     panel.background = element_blank()
#   ) +
#   labs(x = "N2 genome position (Mb)", y = "WS contig position (Mb)", title = paste0(hdr_aln$strain))
# 
# 
# 
# span_one_side <- nucmer_slope_final %>%
#   # dplyr::filter(spans_hdr == T) %>%
#   # dplyr::group_by(chrom,og_hdr_start, og_hdr_end) %>%
#   # dplyr::summarise(
#   #   any_start = any(spans_hdr_start, na.rm = TRUE),
#   #   any_end   = any(spans_hdr_end,   na.rm = TRUE), 
#   #   n_align   = n(), .groups = "drop") #%>%
#   # dplyr::filter(xor(any_start, any_end))
#   dplyr::filter(spans_hdr == T) %>%
#   dplyr::filter(!(spans_hdr_start == T & spans_hdr_end == T)) %>%
#   dplyr::group_by(chrom, og_hdr_start, og_hdr_end) %>%
#   dplyr::mutate(spans_onlyOne = ifelse(any(spans_hdr_start == T) & all(spans_hdr_end == FALSE) | any(spans_hdr_end == TRUE) & all(spans_hdr_start == FALSE), TRUE, FALSE)) %>%
#   dplyr::ungroup() %>%
#   dplyr::select(chrom, og_hdr_start, og_hdr_end, spans_hdr, spans_hdr_start, spans_hdr_end, spans_onlyOne) %>%
#   dplyr::filter(spans_onlyOne == T) %>%
#   dplyr::distinct(chrom, og_hdr_start)
# 
# num_one_side_span <- nrow(span_one_side)
# print(paste0("Number of HDRs that only have an alignment overlapping with either the start or end HDR boundaries, but not both: ", num_one_side_span))
# span_none <- nucmer_slope_final %>%
#   dplyr::group_by(chrom, og_hdr_start, og_hdr_end) %>%
#   dplyr::filter(all(spans_hdr == F)) %>% dplyr::distinct(chrom, og_hdr_start)
# 
# num_noSpan <- nrow(span_none)
# print(paste0("Number of HDRs that don't have any alignments overlapping with HDR boundaries: ", num_noSpan))
# 
# # Finalizing set of WS HDR regionss
# ws_n2_hdr_region <- nucmer_slope_final %>%
#   dplyr::group_by(chrom, og_hdr_start, og_hdr_end) %>%
#   dplyr::mutate(hdr_size = og_hdr_end - og_hdr_start, ws_hdr_size = WS_hdr_end_max_updated - WS_hdr_start_min_updated) %>%
#   dplyr::ungroup() %>%
#   dplyr::distinct(chrom, og_hdr_start, .keep_all = T) %>%
#   dplyr::mutate(ws_hdr_size = as.numeric(ws_hdr_size))
# 
# max_size <- ws_n2_hdr_region %>% dplyr::select(ws_hdr_size) %>% dplyr::arrange(desc(ws_hdr_size)) %>% dplyr::slice_head(n = 1)
# 
# # Shows that HDRs are very structurally different than N2 and not just high SNV density
# ggplot(data = ws_n2_hdr_region) +
#   geom_point(aes(x = hdr_size / 1e6, y = abs(ws_hdr_size / 1e6)), color = 'firebrick') +
#   geom_line(data = data.frame(x = c(0, 1.2)), aes(x = x, y = x), linetype = "dashed") +
#   theme_bw() +
#   labs(x = "N2 HDR size (Mb)", y = paste0(strain_id, " HDR size (Mb)"))
# 
# print(paste0("The largest HDR lift-over in ",SOI," is ", round(max_size)))
# 
# # Distinct wild strain HDR lift-overs!
# ws_hdr_coords <- nucmer_slope_final %>% 
#   dplyr::distinct(chrom,N2S,N2E,og_hdr_start,og_hdr_end,longest_contig,WS_hdr_start_min_updated,WS_hdr_end_max_updated) %>%
#   dplyr::filter(!is.na(WS_hdr_start_min_updated) & !is.na(WS_hdr_end_max_updated))
# 
# number_hdrs <- ws_hdr_coords %>% 
#   dplyr::distinct(chrom,og_hdr_start,og_hdr_end, longest_contig, WS_hdr_start_min_updated, WS_hdr_end_max_updated) 
# num <- nrow(number_hdrs)
# 
# two_starts <- nucmer_slope_final %>% dplyr::filter(chrom == "II" & og_hdr_start == "1040000") %>%
#   dplyr::select(chrom, og_hdr_start, og_hdr_end, longest_contig, N2S, N2E, WSS, WSE, WS_hdr_start_min_updated, WS_hdr_end_max_updated, spans_hdr, L2)
# 
# ggplot(data = two_starts) +
#   geom_rect(data = two_starts %>% dplyr::distinct(og_hdr_start, .keep_all = T), aes(xmin = og_hdr_start / 1e6, xmax = og_hdr_end / 1e6, ymin = -Inf, ymax = Inf), fill = "#DB6333", alpha = 0.3) +
#   geom_rect(data = two_starts %>% dplyr::distinct(WS_hdr_start_min_updated, .keep_all = T), aes(xmin = -Inf, xmax = Inf, ymin = WS_hdr_start_min_updated / 1e6, ymax = WS_hdr_end_max_updated / 1e6), fill = 'blue', alpha = 0.3) +
#   geom_segment(aes(x = og_hdr_start / 1e6, xend = og_hdr_start / 1e6, y = -Inf, yend = Inf), color = 'black') +
#   geom_segment(aes(x = og_hdr_end / 1e6, xend = og_hdr_end / 1e6, y = -Inf, yend = Inf), color = 'black') +
#   geom_segment(aes(x = -Inf, xend = Inf, y = WS_hdr_start_min_updated / 1e6, yend = WS_hdr_start_min_updated / 1e6), color = 'black') +
#   geom_segment(aes(x = -Inf, xend = Inf, y = WS_hdr_end_max_updated / 1e6, yend = WS_hdr_end_max_updated / 1e6), color = 'black') +
#   geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6,  color = spans_hdr), size = 1) +
#   scale_color_manual(values = c("TRUE" = "green", "FALSE" = "blue")) +
#   # coord_cartesian(ylim = c(0,5)) +
#   facet_wrap(~chrom) +
#   theme(
#     panel.border = element_rect(color = 'black', fill = NA),
#     panel.background = element_blank()
#   ) +
#   labs(x = "N2 genome position (Mb)", y = "WS contig position (Mb)", title = paste0(hdr_aln$strain))
# 
# 
# # aln_lengths <- nucmer_slope_final %>% dplyr::select(chrom, og_hdr_start, og_hdr_end, WS_hdr_start, WS_hdr_end, inv, WS_hdr_start_min_updated, WS_hdr_end_max_updated, L2) %>%
# #   dplyr::filter(!is.na(WS_hdr_start_min_updated) & !is.na(WS_hdr_end_max_updated)) %>%
# #   dplyr::filter((inv == T & WS_hdr_end == WS_hdr_start_min_updated | inv == T & WS_hdr_start == WS_hdr_end_max_updated) | 
# #                   (inv == F & WS_hdr_start == WS_hdr_start_min_updated | inv == F & WS_hdr_end == WS_hdr_end_max_updated))
# # 
# # number_hdrs_correct <- number_hdrs %>%
# #   dplyr::left_join(aln_lengths, by = c("chrom","og_hdr_start","og_hdr_end","WS_hdr_start_min_updated","WS_hdr_end_max_updated")) %>%
# #   dplyr::group_by(chrom, og_hdr_start) %>%
# #   dplyr::mutate(count = n()) %>%
# #   dplyr::ungroup() %>%
# #   dplyr::group_by(chrom, og_hdr_start, og_hdr_end) %>%
# #   dplyr::filter(L2 == max(L2))
# # 
# # num <- nrow(number_hdrs_correct)
# 
# 
# # !!!!!!!!!!!!!!!!!!!!! THIS ONLY WORKS IF THERE ARE TWO ALIGNMENTS OVERLAPPING WITH A BOUNDARY...... EDGE CASES BELOW 
# two_start_count <- two_starts %>% dplyr::group_by(chrom, og_hdr_start, WS_hdr_start_min_updated) %>%
#   dplyr::mutate(hdr_start_count = n()) %>%
#   dplyr::ungroup() %>%
#   dplyr::group_by(chrom, og_hdr_start, WS_hdr_end_max_updated) %>%
#   dplyr::mutate(hdr_end_count = n()) %>%
#   dplyr::ungroup() %>%
#   dplyr::group_by(chrom, og_hdr_start, og_hdr_end) %>%
#   dplyr::filter(hdr_start_count == max(hdr_start_count) & hdr_end_count == max(hdr_end_count)) 
# 
# FINAL_hdrs <- ws_hdr_coords %>% 
#   dplyr::group_by(chrom, og_hdr_start, WS_hdr_start_min_updated) %>%
#   dplyr::mutate(hdr_start_count = n()) %>%
#   dplyr::ungroup() %>%
#   dplyr::group_by(chrom, og_hdr_start, WS_hdr_end_max_updated) %>%
#   dplyr::mutate(hdr_end_count = n()) %>%
#   dplyr::ungroup() %>%
#   dplyr::group_by(chrom, og_hdr_start, og_hdr_end) %>%
#   dplyr::filter(hdr_start_count == max(hdr_start_count) & hdr_end_count == max(hdr_end_count)) %>%
#   dplyr::distinct(chrom,og_hdr_start,og_hdr_end,longest_contig, WS_hdr_start_min_updated, WS_hdr_end_max_updated)
# 
# no_hdr_liftover <- dplyr::anti_join(ws_hdr_coords, FINAL_hdrs, by = c("chrom","og_hdr_start","og_hdr_end", "WS_hdr_start_min_updated", "WS_hdr_end_max_updated"))
# 
# plot_test <- no_hdr_liftover %>% dplyr::left_join(nucmer_slope_final, by = c("chrom","og_hdr_start","og_hdr_end","WS_hdr_start_min_updated","WS_hdr_end_max_updated"))
# plot_test1 <- nucmer_slope_final %>% dplyr::filter(chrom == "II" & og_hdr_start == "3309000")
# plot_test2 <- nucmer_slope_final %>% dplyr::filter(chrom == "X" & og_hdr_start == "17099000")
# plot_test3 <- nucmer_slope_final %>% dplyr::filter(chrom == "I" & og_hdr_start == "12266000")
# plot_test4 <- nucmer_slope_final %>% dplyr::filter(chrom == "II" & og_hdr_start == "1703000")
# plot_test5 <- nucmer_slope_final %>% dplyr::filter(chrom == "II" & og_hdr_start == "1040000")
# 
# 
# ggplot(data = plot_test1) +
#   geom_rect(data = plot_test1 %>% dplyr::distinct(og_hdr_start, .keep_all = T), aes(xmin = og_hdr_start / 1e6, xmax = og_hdr_end / 1e6, ymin = -Inf, ymax = Inf), fill = "#DB6333", alpha = 0.3) +
#   geom_rect(data = plot_test1 %>% dplyr::distinct(WS_hdr_start_min_updated, .keep_all = T), aes(xmin = -Inf, xmax = Inf, ymin = WS_hdr_start_min_updated / 1e6, ymax = WS_hdr_end_max_updated / 1e6), fill = 'blue', alpha = 0.3) +
#   geom_segment(aes(x = og_hdr_start / 1e6, xend = og_hdr_start / 1e6, y = -Inf, yend = Inf), color = 'black') +
#   geom_segment(aes(x = og_hdr_end / 1e6, xend = og_hdr_end / 1e6, y = -Inf, yend = Inf), color = 'black') +
#   geom_segment(aes(x = -Inf, xend = Inf, y = WS_hdr_start_min_updated / 1e6, yend = WS_hdr_start_min_updated / 1e6), color = 'black') +
#   geom_segment(aes(x = -Inf, xend = Inf, y = WS_hdr_end_max_updated / 1e6, yend = WS_hdr_end_max_updated / 1e6), color = 'black') +
#   
#   # geom_segment(aes(x = -Inf, xend = Inf, y = outside_hdr_min_left / 1e6, yend = outside_hdr_min_left / 1e6), color = 'magenta3', size = 2) +
#   # geom_segment(aes(x = -Inf, xend = Inf, y = outside_hdr_max_left / 1e6, yend = outside_hdr_max_left / 1e6), color = 'magenta3', size = 2) +
#   # geom_segment(aes(x = -Inf, xend = Inf, y = outside_hdr_min_right / 1e6, yend = outside_hdr_min_right / 1e6), color = 'purple', size = 2) +
#   # geom_segment(aes(x = -Inf, xend = Inf, y = outside_hdr_max_right / 1e6, yend = outside_hdr_max_right / 1e6), color = 'purple', size = 2) +
#   # geom_segment(aes(x = -Inf, xend = Inf, y = inside_hdr_max / 1e6, yend = inside_hdr_max / 1e6), color = 'red', size = 1) +
#   # geom_segment(aes(x = -Inf, xend = Inf, y = inside_hdr_min / 1e6, yend = inside_hdr_min / 1e6), color = 'red', size = 1) +
#   
#   geom_segment(aes(x = N2S / 1e6, xend = N2E/ 1e6, y = WSS / 1e6, yend = WSE / 1e6,  color = spans_hdr), size = 1) +
#   scale_color_manual(values = c("TRUE" = "green", "FALSE" = "blue")) +
#   # coord_cartesian(ylim = c(0,5)) +
#   facet_wrap(~chrom) +
#   theme(
#     panel.border = element_rect(color = 'black', fill = NA),
#     panel.background = element_blank()
#   ) +
#   labs(x = "N2 genome position (Mb)", y = "WS contig position (Mb)", title = paste0(hdr_aln$strain))
# 
# ggplot(data = plot_test2) +
#   geom_rect(data = plot_test2 %>% dplyr::distinct(og_hdr_start, .keep_all = T), aes(xmin = og_hdr_start / 1e6, xmax = og_hdr_end / 1e6, ymin = -Inf, ymax = Inf), fill = "#DB6333", alpha = 0.3) +
#   geom_rect(data = plot_test2 %>% dplyr::distinct(WS_hdr_start_min_updated, .keep_all = T), aes(xmin = -Inf, xmax = Inf, ymin = WS_hdr_start_min_updated / 1e6, ymax = WS_hdr_end_max_updated / 1e6), fill = 'blue', alpha = 0.3) +
#   geom_segment(aes(x = og_hdr_start / 1e6, xend = og_hdr_start / 1e6, y = -Inf, yend = Inf), color = 'black') +
#   geom_segment(aes(x = og_hdr_end / 1e6, xend = og_hdr_end / 1e6, y = -Inf, yend = Inf), color = 'black') +
#   geom_segment(aes(x = -Inf, xend = Inf, y = WS_hdr_start_min_updated / 1e6, yend = WS_hdr_start_min_updated / 1e6), color = 'black') +
#   geom_segment(aes(x = -Inf, xend = Inf, y = WS_hdr_end_max_updated / 1e6, yend = WS_hdr_end_max_updated / 1e6), color = 'black') +
#   
#   geom_segment(aes(x = -Inf, xend = Inf, y = outside_hdr_min_left / 1e6, yend = outside_hdr_min_left / 1e6), color = 'magenta3', size = 2) +
#   geom_segment(aes(x = -Inf, xend = Inf, y = outside_hdr_max_left / 1e6, yend = outside_hdr_max_left / 1e6), color = 'magenta3', size = 2) +
#   geom_segment(aes(x = -Inf, xend = Inf, y = outside_hdr_min_right / 1e6, yend = outside_hdr_min_right / 1e6), color = 'purple', size = 2) +
#   geom_segment(aes(x = -Inf, xend = Inf, y = outside_hdr_max_right / 1e6, yend = outside_hdr_max_right / 1e6), color = 'purple', size = 2) +
#   geom_segment(aes(x = -Inf, xend = Inf, y = inside_hdr_max / 1e6, yend = inside_hdr_max / 1e6), color = 'red', size = 1) +
#   geom_segment(aes(x = -Inf, xend = Inf, y = inside_hdr_min / 1e6, yend = inside_hdr_min / 1e6), color = 'red', size = 1) +
#   
#   geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6,  color = spans_hdr), size = 1) +
#   scale_color_manual(values = c("TRUE" = "green", "FALSE" = "blue")) +
#   # coord_cartesian(ylim = c(0,5)) +
#   facet_wrap(~chrom) +
#   theme(
#     panel.border = element_rect(color = 'black', fill = NA),
#     panel.background = element_blank()
#   ) +
#   labs(x = "N2 genome position (Mb)", y = "WS contig position (Mb)", title = paste0(hdr_aln$strain))
# 
# ggplot(data = plot_test3) +
#   geom_rect(data = plot_test3 %>% dplyr::distinct(og_hdr_start, .keep_all = T), aes(xmin = og_hdr_start / 1e6, xmax = og_hdr_end / 1e6, ymin = -Inf, ymax = Inf), fill = "#DB6333", alpha = 0.3) +
#   geom_rect(data = plot_test3 %>% dplyr::distinct(WS_hdr_start_min_updated, .keep_all = T), aes(xmin = -Inf, xmax = Inf, ymin = WS_hdr_start_min_updated / 1e6, ymax = WS_hdr_end_max_updated / 1e6), fill = 'blue', alpha = 0.3) +
#   geom_segment(aes(x = og_hdr_start / 1e6, xend = og_hdr_start / 1e6, y = -Inf, yend = Inf), color = 'black') +
#   geom_segment(aes(x = og_hdr_end / 1e6, xend = og_hdr_end / 1e6, y = -Inf, yend = Inf), color = 'black') +
#   geom_segment(aes(x = -Inf, xend = Inf, y = WS_hdr_start_min_updated / 1e6, yend = WS_hdr_start_min_updated / 1e6), color = 'black') +
#   geom_segment(aes(x = -Inf, xend = Inf, y = WS_hdr_end_max_updated / 1e6, yend = WS_hdr_end_max_updated / 1e6), color = 'black') +
#   
#   geom_segment(aes(x = -Inf, xend = Inf, y = outside_hdr_min_left / 1e6, yend = outside_hdr_min_left / 1e6), color = 'magenta3', size = 2) +
#   geom_segment(aes(x = -Inf, xend = Inf, y = outside_hdr_max_left / 1e6, yend = outside_hdr_max_left / 1e6), color = 'magenta3', size = 2) +
#   geom_segment(aes(x = -Inf, xend = Inf, y = outside_hdr_min_right / 1e6, yend = outside_hdr_min_right / 1e6), color = 'purple', size = 2) +
#   geom_segment(aes(x = -Inf, xend = Inf, y = outside_hdr_max_right / 1e6, yend = outside_hdr_max_right / 1e6), color = 'purple', size = 2) +
#   geom_segment(aes(x = -Inf, xend = Inf, y = inside_hdr_max / 1e6, yend = inside_hdr_max / 1e6), color = 'red', size = 1) +
#   geom_segment(aes(x = -Inf, xend = Inf, y = inside_hdr_min / 1e6, yend = inside_hdr_min / 1e6), color = 'red', size = 1) +
#   
#   geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6,  color = spans_hdr), size = 1) +
#   scale_color_manual(values = c("TRUE" = "green", "FALSE" = "blue")) +
#   # coord_cartesian(ylim = c(0,5)) +
#   facet_wrap(~chrom) +
#   theme(
#     panel.border = element_rect(color = 'black', fill = NA),
#     panel.background = element_blank()
#   ) +
#   labs(x = "N2 genome position (Mb)", y = "WS contig position (Mb)", title = paste0(hdr_aln$strain))
# 
# 
# ggplot(data = plot_test4) +
#   geom_rect(data = plot_test4 %>% dplyr::distinct(og_hdr_start, .keep_all = T), aes(xmin = og_hdr_start / 1e6, xmax = og_hdr_end / 1e6, ymin = -Inf, ymax = Inf), fill = "#DB6333", alpha = 0.3) +
#   geom_rect(data = plot_test4 %>% dplyr::distinct(WS_hdr_start_min_updated, .keep_all = T), aes(xmin = -Inf, xmax = Inf, ymin = WS_hdr_start_min_updated / 1e6, ymax = WS_hdr_end_max_updated / 1e6), fill = 'blue', alpha = 0.3) +
#   geom_segment(aes(x = og_hdr_start / 1e6, xend = og_hdr_start / 1e6, y = -Inf, yend = Inf), color = 'black') +
#   geom_segment(aes(x = og_hdr_end / 1e6, xend = og_hdr_end / 1e6, y = -Inf, yend = Inf), color = 'black') +
#   geom_segment(aes(x = -Inf, xend = Inf, y = WS_hdr_start_min_updated / 1e6, yend = WS_hdr_start_min_updated / 1e6), color = 'black') +
#   geom_segment(aes(x = -Inf, xend = Inf, y = WS_hdr_end_max_updated / 1e6, yend = WS_hdr_end_max_updated / 1e6), color = 'black') +
#   
#   geom_segment(aes(x = -Inf, xend = Inf, y = outside_hdr_min_left / 1e6, yend = outside_hdr_min_left / 1e6), color = 'magenta3', size = 2) +
#   geom_segment(aes(x = -Inf, xend = Inf, y = outside_hdr_max_left / 1e6, yend = outside_hdr_max_left / 1e6), color = 'magenta3', size = 2) +
#   geom_segment(aes(x = -Inf, xend = Inf, y = outside_hdr_min_right / 1e6, yend = outside_hdr_min_right / 1e6), color = 'purple', size = 2) +
#   geom_segment(aes(x = -Inf, xend = Inf, y = outside_hdr_max_right / 1e6, yend = outside_hdr_max_right / 1e6), color = 'purple', size = 2) +
#   geom_segment(aes(x = -Inf, xend = Inf, y = inside_hdr_max / 1e6, yend = inside_hdr_max / 1e6), color = 'red', size = 1) +
#   geom_segment(aes(x = -Inf, xend = Inf, y = inside_hdr_min / 1e6, yend = inside_hdr_min / 1e6), color = 'red', size = 1) +
#   
#   geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6,  color = spans_hdr), size = 1) +
#   scale_color_manual(values = c("TRUE" = "green", "FALSE" = "blue")) +
#   # coord_cartesian(ylim = c(0,5)) +
#   facet_wrap(~chrom) +
#   theme(
#     panel.border = element_rect(color = 'black', fill = NA),
#     panel.background = element_blank()
#   ) +
#   labs(x = "N2 genome position (Mb)", y = "WS contig position (Mb)", title = paste0(hdr_aln$strain))
# 
# 
# 
# ggplot(data = plot_test5) +
#   geom_rect(data = plot_test5 %>% dplyr::distinct(og_hdr_start, .keep_all = T), aes(xmin = og_hdr_start / 1e6, xmax = og_hdr_end / 1e6, ymin = -Inf, ymax = Inf), fill = "#DB6333", alpha = 0.3) +
#   geom_rect(data = plot_test5 %>% dplyr::distinct(WS_hdr_start_min_updated, .keep_all = T), aes(xmin = -Inf, xmax = Inf, ymin = WS_hdr_start_min_updated / 1e6, ymax = WS_hdr_end_max_updated / 1e6), fill = 'blue', alpha = 0.3) +
#   geom_segment(aes(x = og_hdr_start / 1e6, xend = og_hdr_start / 1e6, y = -Inf, yend = Inf), color = 'black') +
#   geom_segment(aes(x = og_hdr_end / 1e6, xend = og_hdr_end / 1e6, y = -Inf, yend = Inf), color = 'black') +
#   geom_segment(aes(x = -Inf, xend = Inf, y = WS_hdr_start_min_updated / 1e6, yend = WS_hdr_start_min_updated / 1e6), color = 'black') +
#   geom_segment(aes(x = -Inf, xend = Inf, y = WS_hdr_end_max_updated / 1e6, yend = WS_hdr_end_max_updated / 1e6), color = 'black') +
#   
#   geom_segment(aes(x = -Inf, xend = Inf, y = outside_hdr_min_left / 1e6, yend = outside_hdr_min_left / 1e6), color = 'magenta3', size = 2) +
#   geom_segment(aes(x = -Inf, xend = Inf, y = outside_hdr_max_left / 1e6, yend = outside_hdr_max_left / 1e6), color = 'magenta3', size = 2) +
#   geom_segment(aes(x = -Inf, xend = Inf, y = outside_hdr_min_right / 1e6, yend = outside_hdr_min_right / 1e6), color = 'purple', size = 2) +
#   geom_segment(aes(x = -Inf, xend = Inf, y = outside_hdr_max_right / 1e6, yend = outside_hdr_max_right / 1e6), color = 'purple', size = 2) +
#   geom_segment(aes(x = -Inf, xend = Inf, y = inside_hdr_max / 1e6, yend = inside_hdr_max / 1e6), color = 'red', size = 1) +
#   geom_segment(aes(x = -Inf, xend = Inf, y = inside_hdr_min / 1e6, yend = inside_hdr_min / 1e6), color = 'red', size = 1) +
#   
#   geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6,  color = spans_hdr), size = 1) +
#   scale_color_manual(values = c("TRUE" = "green", "FALSE" = "blue")) +
#   # coord_cartesian(ylim = c(0,5)) +
#   facet_wrap(~chrom) +
#   theme(
#     panel.border = element_rect(color = 'black', fill = NA),
#     panel.background = element_blank()
#   ) +
#   labs(x = "N2 genome position (Mb)", y = "WS contig position (Mb)", title = paste0(hdr_aln$strain))
# 
# 
# 
# 
# print(paste0("Total number of HDRs lifted over to wild strains: ", nrow(FINAL_hdrs), " / ", num_hdrs))
# 
# # two_starts_after <- number_hdrs_correct %>% dplyr::filter(chrom == "II" & og_hdr_start == "1040000") %>%
# # dplyr::left_join(two_starts, by = c("chrom","og_hdr_start","og_hdr_end", "WS_hdr_start_min_updated", "WS_hdr_end_max_updated"))
# # dplyr::filter(!is.na(L2.y))
# 
# ggplot(data = two_start_count) +
#   geom_rect(data = two_start_count %>% dplyr::distinct(og_hdr_start, .keep_all = T), aes(xmin = og_hdr_start / 1e6, xmax = og_hdr_end / 1e6, ymin = -Inf, ymax = Inf), fill = "#DB6333", alpha = 0.3) +
#   geom_rect(data = two_start_count %>% dplyr::distinct(WS_hdr_start_min_updated, .keep_all = T), aes(xmin = -Inf, xmax = Inf, ymin = WS_hdr_start_min_updated / 1e6, ymax = WS_hdr_end_max_updated / 1e6), fill = 'blue', alpha = 0.3) +
#   geom_segment(aes(x = og_hdr_start / 1e6, xend = og_hdr_start / 1e6, y = -Inf, yend = Inf), color = 'black') +
#   geom_segment(aes(x = og_hdr_end / 1e6, xend = og_hdr_end / 1e6, y = -Inf, yend = Inf), color = 'black') +
#   geom_segment(aes(x = -Inf, xend = Inf, y = WS_hdr_start_min_updated / 1e6, yend = WS_hdr_start_min_updated / 1e6), color = 'black') +
#   geom_segment(aes(x = -Inf, xend = Inf, y = WS_hdr_end_max_updated / 1e6, yend = WS_hdr_end_max_updated / 1e6), color = 'black') +
#   geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6,  color = spans_hdr), size = 1) +
#   scale_color_manual(values = c("TRUE" = "green", "FALSE" = "blue")) +
#   # coord_cartesian(ylim = c(0,5)) +
#   facet_wrap(~chrom) +
#   theme(
#     panel.border = element_rect(color = 'black', fill = NA),
#     panel.background = element_blank()
#   ) +
#   labs(x = "N2 genome position (Mb)", y = "WS contig position (Mb)", title = paste0(hdr_aln$strain))
# 
# 
# # FINAL_hdrs <- number_hdrs_correct %>%
# #   dplyr::left_join(nucmer_slope_final, by = c("chrom","og_hdr_start","og_hdr_end", "WS_hdr_start_min_updated", "WS_hdr_end_max_updated")) %>%
# #   dplyr::mutate(WS_hdr_start_min_updated = ifelse(WS_hdr_start_min_updated < 0, 0, WS_hdr_start_min_updated))
# 
# ggplot(FINAL_hdrs) + 
#   geom_rect(data = FINAL_hdrs %>% dplyr::distinct(og_hdr_start, .keep_all = T), aes(xmin = og_hdr_start / 1e6, xmax = og_hdr_end / 1e6, ymin = -Inf, ymax = Inf), fill = "#DB6333", alpha = 0.3) +
#   geom_rect(data = FINAL_hdrs %>% dplyr::distinct(WS_hdr_start_min_updated, .keep_all = T), aes(xmin = -Inf, xmax = Inf, ymin = WS_hdr_start_min_updated / 1e6, ymax = WS_hdr_end_max_updated / 1e6), fill = 'blue', alpha = 0.3) +
#   # geom_segment(aes(x = og_hdr_start / 1e6, xend = og_hdr_start / 1e6, y = -Inf, yend = Inf), color = 'black') +
#   # geom_segment(aes(x = og_hdr_end / 1e6, xend = og_hdr_end / 1e6, y = -Inf, yend = Inf), color = 'black') +
#   # geom_segment(aes(x = -Inf, xend = Inf, y = WS_hdr_start_min / 1e6, yend = WS_hdr_start_min / 1e6), color = 'black') +
#   # geom_segment(aes(x = -Inf, xend = Inf, y = WS_hdr_end_max / 1e6, yend = WS_hdr_end_max / 1e6), color = 'black') +
#   # geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6,  color = spans_hdr), size = 1) +
#   scale_color_manual(values = c("TRUE" = "green", "FALSE" = "blue")) +
#   facet_wrap(~chrom, scales = "free") +
#   theme(
#     panel.border = element_rect(color = 'black', fill = NA),
#     panel.background = element_blank()
#   ) +
#   labs(x = "N2 genome position (Mb)", y = "WS contig position (Mb)", title = paste0(hdr_aln$strain))
# 
# 
# 
# 
# 
# 
# 
# ################## Look at HDRs that have "overlapping alignments" but then do not have HDR liftover ################################
# final_count <- FINAL_hdrs %>% dplyr::distinct(chrom,og_hdr_start) # 308
# pre_count <- nucmer_slope_final %>% dplyr::distinct(chrom,og_hdr_start) # 315 - duplicated for HDRs that have multiple starts and stops calculated
# 
# no_hdr_liftover <- dplyr::anti_join(nucmer_slope, FINAL_hdrs, by = c("chrom","og_hdr_start","og_hdr_end"))
# 
# ggplot(data = no_hdr_liftover) +
#   geom_rect(aes(xmin = og_hdr_start / 1e6, xmax = og_hdr_end / 1e6, ymin = -Inf, ymax = Inf), fill = "gray", alpha = 0.1) +
#   geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6,  color = longest_contig), size = 1) +
#   facet_wrap(~chrom) +
#   theme(
#     panel.border = element_rect(color = 'black', fill = NA),
#     panel.background = element_blank()) +
#   labs(x = "N2 genome position (Mb)", y = "WS contig position (Mb)", title = paste0(hdr_aln$strain))
# 
# 
# ggplot(data = nucmer_slope_final %>% dplyr::filter(og_hdr_start == "1269000")) +
#   geom_rect(aes(xmin = og_hdr_start / 1e6, xmax = og_hdr_end / 1e6, ymin = -Inf, ymax = Inf), fill = "gray", alpha = 0.1) +
#   geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6,  color = longest_contig), size = 1) +
#   facet_wrap(~chrom) +
#   theme(
#     panel.border = element_rect(color = 'black', fill = NA),
#     panel.background = element_blank()) +
#   # coord_cartesian(xlim = c(19.4,19.8)) +
#   labs(x = "N2 genome position (Mb)", y = "WS contig position (Mb)", title = paste0(hdr_aln$strain))
# 
# 
# ggplot(data = nucmer_slope_final %>% dplyr::filter(og_hdr_start == "13517000" & chrom == "II")) +
#   geom_rect(aes(xmin = og_hdr_start / 1e6, xmax = og_hdr_end / 1e6, ymin = -Inf, ymax = Inf), fill = "gray", alpha = 0.1) +
#   geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6,  color = longest_contig), size = 1) +
#   facet_wrap(~chrom) +
#   theme(
#     panel.border = element_rect(color = 'black', fill = NA),
#     panel.background = element_blank()) +
#   # coord_cartesian(xlim = c(19.4,19.8)) +
#   labs(x = "N2 genome position (Mb)", y = "WS contig position (Mb)", title = paste0(hdr_aln$strain))
# 
# 
# test_final <- FINAL_hdrs %>% dplyr::filter(chrom == "V" & og_hdr_start > 17000000 & og_hdr_end < 17250000)
# test_finalPRE <- nucmer_slope_final %>% dplyr::filter(chrom == "V" & og_hdr_start > 17000000 & og_hdr_end < 17250000)
# 
# # THESE SHOULD LOOK THE SAME......
# cowplot::plot_grid( ggplot(test_final) + 
#                       geom_rect(data = test_final %>% dplyr::distinct(og_hdr_start, .keep_all = T), aes(xmin = og_hdr_start / 1e6, xmax = og_hdr_end / 1e6, ymin = -Inf, ymax = Inf), fill = "#DB6333", alpha = 0.3) +
#                       # geom_rect(data = no_aln_bound_inv %>% dplyr::distinct(WS_hdr_start_min, .keep_all = T), aes(xmin = -Inf, xmax = Inf, ymin = WS_hdr_start_min / 1e6, ymax = WS_hdr_end_max / 1e6), fill = 'blue', alpha = 0.3) +
#                       geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6,  color = spans_hdr), size = 1) +
#                       geom_rect(data = test_final %>% dplyr::distinct(WS_hdr_start_min_updated, .keep_all = T), aes(xmin = -Inf, xmax = Inf, ymin = WS_hdr_start_min_updated / 1e6, ymax = WS_hdr_end_max_updated / 1e6), fill = 'blue', alpha = 0.3) +
#                       geom_segment(aes(x = og_hdr_start / 1e6, xend = og_hdr_start / 1e6, y = -Inf, yend = Inf), color = 'black') +
#                       
#                       geom_segment(aes(x = -Inf, xend = Inf, y = outside_hdr_min_left / 1e6, yend = outside_hdr_min_left / 1e6), color = 'magenta3', size = 2) +
#                       geom_segment(aes(x = -Inf, xend = Inf, y = outside_hdr_max_left / 1e6, yend = outside_hdr_max_left / 1e6), color = 'magenta3', size = 2) +
#                       geom_segment(aes(x = -Inf, xend = Inf, y = outside_hdr_min_right / 1e6, yend = outside_hdr_min_right / 1e6), color = 'purple', size = 2) +
#                       geom_segment(aes(x = -Inf, xend = Inf, y = outside_hdr_max_right / 1e6, yend = outside_hdr_max_right / 1e6), color = 'purple', size = 2) +
#                       geom_segment(aes(x = -Inf, xend = Inf, y = inside_hdr_max / 1e6, yend = inside_hdr_max / 1e6), color = 'red', size = 1) +
#                       geom_segment(aes(x = -Inf, xend = Inf, y = inside_hdr_min / 1e6, yend = inside_hdr_min / 1e6), color = 'red', size = 1) +
#                       scale_color_manual(values = c("TRUE" = "green", "FALSE" = "blue")) +
#                       facet_wrap(~chrom, scales = "free") +
#                       theme(
#                         panel.border = element_rect(color = 'black', fill = NA),
#                         panel.background = element_blank()
#                       ) +
#                       labs(x = "N2 genome position (Mb)", y = "WS contig position (Mb)", title = paste0(hdr_aln$strain)),
#                     
#                     ggplot(test_finalPRE) + 
#                       geom_rect(data = test_finalPRE %>% dplyr::distinct(og_hdr_start, .keep_all = T), aes(xmin = og_hdr_start / 1e6, xmax = og_hdr_end / 1e6, ymin = -Inf, ymax = Inf), fill = "#DB6333", alpha = 0.3) +
#                       # geom_rect(data = no_aln_bound_inv %>% dplyr::distinct(WS_hdr_start_min, .keep_all = T), aes(xmin = -Inf, xmax = Inf, ymin = WS_hdr_start_min / 1e6, ymax = WS_hdr_end_max / 1e6), fill = 'blue', alpha = 0.3) +
#                       geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6,  color = spans_hdr), size = 1) +
#                       geom_rect(data = test_finalPRE %>% dplyr::distinct(WS_hdr_start_min_updated, .keep_all = T), aes(xmin = -Inf, xmax = Inf, ymin = WS_hdr_start_min_updated / 1e6, ymax = WS_hdr_end_max_updated / 1e6), fill = 'blue', alpha = 0.3) +
#                       geom_segment(aes(x = og_hdr_start / 1e6, xend = og_hdr_start / 1e6, y = -Inf, yend = Inf), color = 'black') +
#                       
#                       geom_segment(aes(x = -Inf, xend = Inf, y = outside_hdr_min_left / 1e6, yend = outside_hdr_min_left / 1e6), color = 'magenta3', size = 2) +
#                       geom_segment(aes(x = -Inf, xend = Inf, y = outside_hdr_max_left / 1e6, yend = outside_hdr_max_left / 1e6), color = 'magenta3', size = 2) +
#                       geom_segment(aes(x = -Inf, xend = Inf, y = outside_hdr_min_right / 1e6, yend = outside_hdr_min_right / 1e6), color = 'purple', size = 2) +
#                       geom_segment(aes(x = -Inf, xend = Inf, y = outside_hdr_max_right / 1e6, yend = outside_hdr_max_right / 1e6), color = 'purple', size = 2) +
#                       geom_segment(aes(x = -Inf, xend = Inf, y = inside_hdr_max / 1e6, yend = inside_hdr_max / 1e6), color = 'red', size = 1) +
#                       geom_segment(aes(x = -Inf, xend = Inf, y = inside_hdr_min / 1e6, yend = inside_hdr_min / 1e6), color = 'red', size = 1) +
#                       scale_color_manual(values = c("TRUE" = "green", "FALSE" = "blue")) +
#                       facet_wrap(~chrom, scales = "free") +
#                       theme(
#                         panel.border = element_rect(color = 'black', fill = NA),
#                         panel.background = element_blank()
#                       ) +
#                       labs(x = "N2 genome position (Mb)", y = "WS contig position (Mb)", title = paste0(hdr_aln$strain)))
# 
# 
# 
# 
# 
# # chaos <- nucmer_slope_final %>% dplyr::filter(og_hdr_start == 17126000)
# # ooo <- chaos %>% dplyr::filter(St2 == 3961588)
# # 
# # ggplot(chaos) + 
# #   geom_rect(data = chaos %>% dplyr::distinct(og_hdr_start, .keep_all = T), aes(xmin = og_hdr_start / 1e6, xmax = og_hdr_end / 1e6, ymin = -Inf, ymax = Inf), fill = "#DB6333", alpha = 0.3) +
# #   # geom_rect(data = no_aln_bound_inv %>% dplyr::distinct(WS_hdr_start_min, .keep_all = T), aes(xmin = -Inf, xmax = Inf, ymin = WS_hdr_start_min / 1e6, ymax = WS_hdr_end_max / 1e6), fill = 'blue', alpha = 0.3) +
# #   geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6,  color = spans_hdr), size = 1) +
# #   geom_segment(data = ooo, aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6), color = 'cyan', size = 4) + # LOWER BOUNDARY IS BEING SET BASED ON PREDICTED END (START, BUT INVERTED) OF WS HDR FROM SLOPE... (WS_hdr_start_min) 
# #   geom_rect(data = chaos %>% dplyr::distinct(WS_hdr_start_min_updated, .keep_all = T), aes(xmin = -Inf, xmax = Inf, ymin = WS_hdr_start_min_updated / 1e6, ymax = WS_hdr_end_max_updated / 1e6), fill = 'blue', alpha = 0.3) +
# #   geom_segment(aes(x = og_hdr_start / 1e6, xend = og_hdr_start / 1e6, y = -Inf, yend = Inf), color = 'black') +
# #   
# #   geom_segment(aes(x = -Inf, xend = Inf, y = outside_hdr_min_left / 1e6, yend = outside_hdr_min_left / 1e6), color = 'magenta3', size = 2) +
# #   geom_segment(aes(x = -Inf, xend = Inf, y = outside_hdr_max_left / 1e6, yend = outside_hdr_max_left / 1e6), color = 'magenta3', size = 2) +
# #   geom_segment(aes(x = -Inf, xend = Inf, y = outside_hdr_min_right / 1e6, yend = outside_hdr_min_right / 1e6), color = 'purple', size = 2) +
# #   geom_segment(aes(x = -Inf, xend = Inf, y = outside_hdr_max_right / 1e6, yend = outside_hdr_max_right / 1e6), color = 'purple', size = 2) +
# #   geom_segment(aes(x = -Inf, xend = Inf, y = inside_hdr_max / 1e6, yend = inside_hdr_max / 1e6), color = 'red', size = 1) +
# #   geom_segment(aes(x = -Inf, xend = Inf, y = inside_hdr_min / 1e6, yend = inside_hdr_min / 1e6), color = 'red', size = 1) +
# #   scale_color_manual(values = c("TRUE" = "green", "FALSE" = "blue")) +
# #   facet_wrap(~chrom, scales = "free") +
# #   theme(
# #     panel.border = element_rect(color = 'black', fill = NA),
# #     panel.background = element_blank()
# #   ) +
# #   labs(x = "N2 genome position (Mb)", y = "WS contig position (Mb)", title = paste0(hdr_aln$strain))