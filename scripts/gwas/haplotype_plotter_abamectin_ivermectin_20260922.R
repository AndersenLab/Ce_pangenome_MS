library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(cowplot)
library(ape)
library(data.table)
library(stringr)
library(ggh4x)
library(Biostrings)

# ========================================================================================================================================================================================================= #
# Load in abamectin and ivermectin trait data
# ========================================================================================================================================================================================================= #
traits <- readr::read_tsv("../../tables/mapping_2018GWAS_MLs.tsv") %>%
  dplyr::select(strain, `Abamectin_q90.TOF_ctrl-regressed`, `Ivermectin_median.TOF_ctrl-regressed`)

geno_matrix <- readr::read_tsv("../../tables/genotype_matrix.tsv")

wild_strains_140 <- readr::read_tsv("../../tables/wild_strain_genome_stats.tsv") %>% dplyr::select(Strain) %>% dplyr::pull()


# ========================================================================================================================================================================================================= #
# First trait to look at is q90 TOF for abamectin (peak marker of V:16198034 log10 7.66)
# ========================================================================================================================================================================================================= #
hdr_chrom = "V"
hdr_start_pos = 16198034 - 22000
hdr_end_pos = 16198034 + 22000
# hdr_start_pos = 15658959 # qtl start
# hdr_end_pos = 18046079 # qtl end


# ========================================================================================================================================================================================================= #
# Visualizing abamectin haplotype for ALL strains we have assemblies for 
# ========================================================================================================================================================================================================= #
# Looking at REF ALT split
abam_geno <- geno_matrix %>% dplyr::filter(CHROM == "V" & POS == 16198034) %>%
  tidyr::pivot_longer(-c(CHROM, POS, REF, ALT), names_to = "strain", values_to = 'abam_gt') %>%
  dplyr::select(-CHROM, -POS, -REF, -ALT) %>%
  dplyr::mutate(abam_gt = ifelse(abam_gt == "-1", "REF", "ALT"))

# Looking at strains we have assemblies for
abam_traits <- traits %>% dplyr::select(strain,  `Abamectin_q90.TOF_ctrl-regressed`) %>% dplyr::rename(abam_q90_TOF =  `Abamectin_q90.TOF_ctrl-regressed`) %>%
  dplyr::left_join(abam_geno, by = "strain") %>%
  dplyr::filter(!is.na(abam_q90_TOF) & !is.na(abam_gt)) %>%
  dplyr::mutate(we_have_asm = ifelse(strain %in% wild_strains_140, "TRUE", "FALSE")) %>%
  dplyr::filter(we_have_asm == TRUE) %>%
  dplyr::group_by(abam_gt) %>%
  dplyr::arrange(desc(abam_q90_TOF)) %>%
  dplyr::mutate(
    extreme = "YES") %>%
  dplyr::ungroup()

# Look at phenotype distribution for strains we have
ggplot(abam_traits %>% dplyr::mutate(abam_gt = factor(abam_gt, levels = c("REF","ALT")))) +
  geom_boxplot(aes(x = abam_gt, y = abam_q90_TOF), outliers = FALSE) +
  geom_jitter(aes(x = abam_gt, y = abam_q90_TOF, color = extreme, shape = we_have_asm), width = 0.2, size = 3.5) +
  scale_color_manual(values = c("YES" = "red", "NO" = "black")) +
  theme(
    legend.position = 'none',
    panel.background = element_blank(),
    panel.border = element_rect(color = 'black', fill = NA),
    axis.title.y = element_text(size = 16, color = 'black'),
    axis.text.y = element_text(size = 14, color = 'black'),
    axis.text.x = element_text(size = 16, color = 'black')) +
  labs(y = "Abamectin q90 TOF", x = NULL, color = NULL)

# Look at strains we have assemblies for in relation to all strains phenotyped
all_pheno_asm <- traits %>% dplyr::select(strain,  `Abamectin_q90.TOF_ctrl-regressed`) %>% dplyr::rename(abam_q90_TOF =  `Abamectin_q90.TOF_ctrl-regressed`) %>%
  dplyr::left_join(abam_geno, by = "strain") %>%
  dplyr::filter(!is.na(abam_q90_TOF) & !is.na(abam_gt)) %>%
  dplyr::mutate(we_have_asm = ifelse(strain %in% wild_strains_140, "TRUE", "FALSE")) %>%
  dplyr::group_by(abam_gt) %>%
  dplyr::arrange(desc(abam_q90_TOF)) %>%
  dplyr::ungroup() 

pxg <- ggplot(all_pheno_asm %>% dplyr::mutate(abam_gt = factor(abam_gt, levels = c("REF","ALT")))) +
  geom_boxplot(aes(x = abam_gt, y = abam_q90_TOF), outliers = FALSE) +
  geom_jitter(aes(x = abam_gt, y = abam_q90_TOF, fill = we_have_asm), shape = 21, color = 'black', width = 0.2, size = 1) +
  scale_fill_manual(values = c("TRUE" = "black", "FALSE" = "white")) +
  theme(
    legend.position = 'none',
    panel.background = element_blank(),
    panel.border = element_rect(color = 'black', fill = NA),
    axis.title.y = element_text(size = 10, color = 'black'),
    axis.text.y = element_text(size = 10, color = 'black'),
    axis.text.x = element_text(size = 10, color = 'black')) +
  labs(y = "Abamectin q90 TOF", x = NULL, color = NULL)
pxg 

want <- abam_traits %>% dplyr::filter(extreme == "YES") %>% dplyr::arrange(desc(abam_gt), abam_q90_TOF) %>% dplyr::pull(strain) %>% c("N2")


# ========================================================================================================================================================================================================= #
# Visualizing abamectin haplotype
# ========================================================================================================================================================================================================= #
# Load in gneome-genome alignments
transformed_coords <- readr::read_tsv("../../processed_data/genome_resources/genome_data/141_nucmer_ECA741CGC1.tsv",col_names = F) 
colnames(transformed_coords) <- c("S1","E1","S2","E2","L1","L2","IDY","LENR","LENQ","REF","HIFI","STRAIN") 
transformed_coords <- transformed_coords %>% dplyr::filter(STRAIN != "ECA396") %>%
  dplyr::filter(STRAIN %in% want)

# Read in gene models
gffCat1 <- readr::read_tsv("../../processed_data/genome_resources/annotation/140WSs_CGC1_longestIso_mRNAandGenes_BRAKER.tsv", col_names = F)
colnames(gffCat1) <- c("seqid","source","type","start","end","score","strand","phase","attributes","STRAIN")
gffCat2 <- ape::read.gff("../../processed_data/genome_resources/annotation/c_elegans.PRJNA13758.WS283.csq.PCfeaturesOnly.longest.w-mRNA_featuresgff3") %>% dplyr::mutate(STRAIN="N2")
gffCat <- rbind(gffCat1 %>% dplyr::filter(STRAIN != "ECA396"), gffCat2) %>% 
  dplyr::filter(STRAIN %in% want)


# Read in orthogroups
orthos <- readr::read_tsv("../../processed_data/orthology/orthofinder/orthofinder_output/Orthogroups.tsv")
strainCol <- colnames(orthos)
ugh <- gsub(".20251012.inbred.blobFiltered.softMasked.braker.longestIso.protein","", strainCol)
ugh2 <- gsub(".20251014.inbred.blobFiltered.softMasked.braker.longestIso.protein","", ugh)
ugh3 <- gsub(".20251124.inbred.blobFiltered.softMasked.braker.longestIso.protein","", ugh2)
ugh4 <- gsub(".20251012.inbred.onlyONT.blobFiltered.softMasked.braker.longestIso.protein","", ugh3)
ugh5 <- gsub(".Nov2025.softMasked.braker.longest.protein","", ugh4)
ugh6 <- gsub(".20251012.inbred.withONT.blobFiltered.softMasked.braker.longestIso.protein","", ugh5)
strainCol_c2 <- gsub("c_elegans.PRJNA13758.WS283.csq.PCfeaturesOnly.longest.protein","N2", ugh6)
colnames(orthos) <- strainCol_c2


# offset lets you explore adjacent regions
offset = 0
hap_chrom = hdr_chrom
hap_start = hdr_start_pos - offset
hap_end = hdr_end_pos + offset 


#use reference coordinates from g2g alginments to pull the contigs that contain the alt haplotypes for the HDR
hap_coords <- transformed_coords %>%
  dplyr::filter((REF == hap_chrom & hap_start >= S1 & hap_start <= E1 ) | 
                  (REF == hap_chrom & hap_end >= S1 & hap_end <= E1) | 
                  (REF == hap_chrom & S1 >= hap_start & E1 <= hap_end)) %>%
  dplyr::mutate(inv=ifelse(S2>E2,T,F))  %>%
  dplyr::mutate(St2=ifelse(inv==T,E2,S2),Et2=ifelse(inv==T,S2,E2))

# naive visualization of g2g alignments for the target region
# multiple contigs may map to the REF region, we need to filter those!
ggplot(hap_coords) + 
  geom_rect(xmin=hap_start/1e6,xmax=hap_end/1e6,ymin=-Inf,ymax=Inf,fill="lightgrey")+
  geom_segment(aes(x=S1/1e6,xend=E1/1e6,y=S2/1e6,yend=E2/1e6,color=HIFI)) +
  facet_wrap(~STRAIN,scales = 'free') +
  xlab("REF genome position (Mb)") +
  ylab("WILD contig position (Mb)") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill=NA),
        # axis.text = element_blank(),
        # axis.ticks = element_blank(),
        legend.position = 'none') 


#keep only the contig with the largest extent of alignment with the REF HDR
tigFilt <- hap_coords %>% 
  dplyr::group_by(STRAIN) %>%
  dplyr::mutate(nalign = n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(STRAIN,HIFI) %>%
  dplyr::mutate(ntig= n()) %>%
  dplyr::mutate(tigsize=sum(L1)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::filter(tigsize == max(tigsize)) %>%
  dplyr::mutate(rangeDiff=(max(S1,E1)-min(S1,E1))-(hap_end-hap_start)) %>%
  dplyr::ungroup()

#important diagnostic plot
#after the optimal contig is selected, we can see if the HDR boundaries are within the alignment
#for proper visualization, the HDR (grey box in plot) needs to be encompassed by the selected contig
#otherwise some haplotypes will be truncated 
#a step to drop genomes with incomplete coverage of the HDR could be added
#have in mind that the alignments look fragmented because of sequence divergence, but the contig is linear in the genome file
ggplot(tigFilt) +
  geom_rect(xmin=hap_start/1e6,xmax=hap_end/1e6,ymin=-Inf,ymax=Inf,fill="lightgrey")+
  geom_segment(aes(x=S1/1e6,xend=E1/1e6,y=S2/1e6,yend=E2/1e6,color=HIFI)) +
  facet_wrap(~STRAIN,scales = 'free') +
  xlab("N2 genome position (Mb)") +
  ylab("WILD contig position (Mb)") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill=NA),
        legend.position = 'none') 

#keep the set of alignments with the largest span (i.e. removes small distant alignments ("jumps"))
tigFilt2 <- tigFilt %>%
  dplyr::arrange(St2) %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::mutate(leadDiff=lead(St2)-Et2) %>%
  dplyr::mutate(jump=ifelse(leadDiff > 1.5E5,1,0)) %>%
  dplyr::mutate(leadDiff=ifelse(is.na(leadDiff),0,leadDiff)) %>%
  dplyr::mutate(run_id = cumsum(c(1, head(jump, -1)))) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(STRAIN,run_id) %>%
  dplyr::mutate(gsize=n()) %>%
  dplyr::mutate(len=abs(Et2-St2)) %>%
  dplyr::mutate(sumlen=sum(len)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::filter(sumlen==max(sumlen)) %>%
  dplyr::select(-gsize) %>%
  dplyr::ungroup()

ggplot(tigFilt2) +
  geom_rect(xmin=hap_start/1e6,xmax=hap_end/1e6,ymin=-Inf,ymax=Inf,fill="lightgrey")+
  geom_segment(aes(x=S1/1e6,xend=E1/1e6,y=S2/1e6,yend=E2/1e6,color=HIFI)) +
  facet_wrap(~STRAIN,scales = 'free') +
  xlab("N2 genome position (Mb)") +
  ylab("WILD contig position (Mb)") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill=NA),
        legend.position = 'none') 

trim_spacer = 2e4
#trims long alignments to the focal region (i.e. hap_start to hap_end, but transformed to the other genome)
tigTrim <- tigFilt2 %>%
  dplyr::arrange(STRAIN,S1) %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(rboundDist=max(S1,E1)-hap_end) %>%
  dplyr::mutate(E2=ifelse(rboundDist>trim_spacer & inv==F,(E2-(rboundDist-trim_spacer)),E2)) %>%
  dplyr::mutate(E2=ifelse(rboundDist>trim_spacer & inv==T,(E2+(rboundDist-trim_spacer)),E2)) %>%
  dplyr::mutate(E1=ifelse(rboundDist>trim_spacer,(E1-(rboundDist-trim_spacer)),E1)) %>%
  dplyr::mutate(lboundDist=hap_start-min(S1,E1)) %>%
  dplyr::mutate(S2=ifelse(lboundDist>trim_spacer & inv==F,(S2+(lboundDist-trim_spacer)),S2)) %>%
  dplyr::mutate(S2=ifelse(lboundDist>trim_spacer & inv==T,(S2-(lboundDist-trim_spacer)),S2)) %>%
  dplyr::mutate(S1=ifelse(lboundDist>trim_spacer,(S1+(lboundDist-trim_spacer)),S1))

ggplot(tigTrim) +
  geom_rect(xmin=hap_start/1e6,xmax=hap_end/1e6,ymin=-Inf,ymax=Inf,fill="lightgrey")+
  geom_segment(aes(x=S1/1e6,xend=E1/1e6,y=S2/1e6,yend=E2/1e6,color=HIFI), size = 2) +
  facet_wrap(~STRAIN,scales = 'free') +
  xlab("N2 genome position (Mb)") +
  ylab("WILD contig position (Mb)") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill=NA),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = 'none') 

# Coloring by REF and ALT
gt_key <- abam_traits %>% dplyr::filter(extreme == "YES") %>% dplyr::select(STRAIN = strain, abam_gt)

tigTrim_TEST <- tigTrim %>% dplyr::left_join(gt_key, by = "STRAIN") %>%
  dplyr::mutate(S2_plot = ifelse(inv == T, tigsize - E2, S2),
                E2_plot = ifelse(inv == T, tigsize - S2, E2),
                St2_plot = ifelse(inv == T, E2_plot, S2_plot),
                Et2_plot = ifelse(inv == T, S2_plot, E2_plot))  
# re-orienting contigs so that they display as non -inverted alignments
strain_order <- tigTrim_TEST %>% dplyr::arrange(desc(abam_gt)) %>% dplyr::distinct(STRAIN) %>% dplyr::pull(STRAIN)

alns <- ggplot(tigTrim_TEST %>% dplyr::mutate(STRAIN = factor(STRAIN, levels = strain_order))) +
  geom_rect(xmin=hap_start/1e6,xmax=hap_end/1e6,ymin=-Inf,ymax=Inf,fill="lightgrey") +
  geom_vline(xintercept = 16198034 / 1e6, color = 'black', linetype = 'dashed', size = 0.4) +
  geom_segment(aes(x=S1/1e6,xend=E1/1e6,y=St2_plot/1e6,yend=Et2_plot/1e6,color=abam_gt), size = 0.75) +
  scale_color_manual(values = c("ALT" = "red", "REF" = 'black')) +
  facet_wrap(~STRAIN, scales = 'free', nrow = ) +
  xlab("N2 genome position (Mb)") +
  ylab("Wild strain contig position (Mb)") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill=NA),
        axis.title = element_text(size = 10, color = 'black'),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = 'none') 
alns 

# Concatenating pheno x geno plot and dotplots 
pheno_dotplot <- cowplot::plot_grid(
  pxg, alns,
  nrow = 1,
  rel_widths = c(0.5,1),
  labels = c("a","b"))

# Save the plot
# ggsave("../../figures/supplementary/pheno_spread_dotplot_GWAS.png", pheno_dotplot, width = 7.5, height = 6, dpi = 600)


#get the minimum and maximum boundary of the WILD genome alignments that contain the HDR
HV_boundary <- tigTrim %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::select(-leadDiff,-inv) %>%
  dplyr::mutate(inv=ifelse(sum(E2)-sum(S2) >0,F,T)) %>%
  dplyr::mutate(refStart=min(S1,E1),refEnd=max(S1,E1)) %>%
  dplyr::mutate(boundStart=min(S2,E2), boundEnd=max(S2,E2)) %>%
  dplyr::distinct(STRAIN, .keep_all = T) %>%
  dplyr::select(HIFI,boundStart,boundEnd,STRAIN,REF,refStart,refEnd,inv) %>%
  dplyr::ungroup() %>%
  dplyr::rename(boundChrom=HIFI)


#filter the concatenated GFF to extract the gene models of each WILD genome contig boundary
wild_genes <- gffCat %>%
  dplyr::filter(type=="gene" & !(STRAIN=="N2")) %>%
  dplyr::mutate(attributes=gsub(";","",attributes)) %>% 
  dplyr::mutate(attributes=gsub("ID=","",attributes)) %>%
  dplyr::select(attributes,seqid,start,end,strand,STRAIN) %>%
  dplyr::rename(Name=attributes)

wild_tran <-  gffCat  %>%
  dplyr::filter(type=="mRNA" & !(STRAIN=="N2")) %>%
  tidyr::separate(attributes,into=c("tranname","Parent"),sep=";Parent=") %>%
  dplyr::mutate(Parent=gsub(";","",Parent)) %>%
  dplyr::mutate(tranname=gsub("ID=","",tranname)) %>%
  dplyr::select(tranname,Parent,STRAIN) %>%
  dplyr::left_join(wild_genes,by=c("Parent"="Name","STRAIN")) %>%
  dplyr::left_join(HV_boundary,by="STRAIN")

#extract the REF genes 
N2Start = min(HV_boundary$refStart)
N2End = max(HV_boundary$refEnd)
N2_genes <- gffCat %>%
  dplyr::filter(type=="gene" & STRAIN=="N2") %>%
  dplyr::mutate(refStart=N2Start,refEnd=N2End) %>%
  dplyr::filter(grepl("biotype=protein_coding",attributes)) %>%
  tidyr::separate(attributes,into=c("pre","post"),sep=';sequence_name=') %>%
  tidyr::separate(post,into=c("seqname","post2"),sep=';biotype=') %>%
  tidyr::separate(pre,into=c("ID","Name","rest2"),sep=";") %>%
  dplyr::mutate(Name=gsub("Name=","",Name)) %>% 
  dplyr::select(seqid,start,end,strand,Name,rest2,seqname,STRAIN,refStart,refEnd,seqname) %>%
  dplyr::mutate(rest2=ifelse(grepl("locus",rest2),gsub("locus=","",rest2),seqname)) %>%
  dplyr::rename(alias=rest2)

#extract the REF protein-coding transcripts
N2_tran <- gffCat %>%
  dplyr::filter(type=="mRNA" & STRAIN=="N2") %>%
  tidyr::separate(attributes, into=c("ID","Parent","Name","wormpep","locus"),sep=';') %>%
  dplyr::mutate(ID=gsub("ID=Transcript:","",ID)) %>%
  tidyr::separate(ID,into = c("fosmid","tseqID",'tranum'),sep="\\.",remove = F) %>%
  dplyr::mutate(tseqname=paste0(fosmid,".",tseqID,".",tranum)) %>%
  dplyr::mutate(Parent=gsub("Parent=Gene:","",Parent)) %>%
  dplyr::filter(Parent %in% N2_genes$Name) %>%
  dplyr::select(tseqname,Parent) %>%
  dplyr::rename(tranname=tseqname) %>%
  dplyr::mutate(tranname=paste0("transcript_",tranname)) %>%
  dplyr::left_join(N2_genes,by=c('Parent'='Name'))


N2_tran_reg <- N2_tran %>%
  dplyr::filter((start >= hap_start & start <= hap_end) | (end >= hap_start & end <= hap_end))  %>%
  dplyr::filter(seqid==hap_chrom)


#get gene list
HV_genelist <- N2_tran_reg$tranname 
#get alt gene names/aliases
aliases <- N2_tran %>% dplyr::select(seqname,tranname,alias)

#minor diagnostic plot to visualize the REF loci captured by the HDR
#this is your REF haplotype
ggplot(N2_tran_reg) + geom_rect(aes(xmin=start,xmax=end,ymin=1,ymax=2))


#filter orthologous groups using REF genes
#this will establish your orthology relationships between REF and WILD haplotypes
all_orthos_unnest <- orthos %>%
  dplyr::mutate(N2 = strsplit(as.character(N2), ",")) %>%
  tidyr::unnest(N2) %>%
  dplyr::mutate(N2=trimws(N2)) %>%
  dplyr::mutate(na_count = rowSums(is.na(.))) %>%
  dplyr::filter(na_count < length(strainCol_c2) - 2) %>%
  dplyr::left_join(N2_tran %>% dplyr::select(tranname,seqid,seqname,start,end),by=c("N2"="tranname"))  %>%
  dplyr::select(-na_count)

filtOrthos <- orthos %>%
  dplyr::filter(grepl(paste(HV_genelist,collapse="|"),N2)) %>%
  dplyr::mutate(N2 = strsplit(as.character(N2), ",")) %>%
  tidyr::unnest(N2) %>%
  dplyr::mutate(N2=trimws(N2)) %>%
  dplyr::left_join(N2_tran_reg %>% dplyr::select(tranname,seqid,seqname,start,end) %>% dplyr::mutate(og_loc="in_region"),by=c("N2"="tranname")) 

inreg_orthos <- filtOrthos %>% dplyr::filter(!is.na(seqid)) 
outreg_orthos <- filtOrthos %>% dplyr::filter(is.na(seqid)) %>% 
  dplyr::select(-seqid,-seqname,-start,-end,-og_loc) %>%
  dplyr::left_join(N2_tran %>% dplyr::select(tranname,seqid,seqname,start,end,refStart,refEnd) %>% 
                     dplyr::mutate(refChrom=hap_chrom) %>%
                     dplyr::mutate(og_loc="out_region"),by=c("N2"="tranname")) %>%
  dplyr::filter(seqid==refChrom) %>%
  dplyr::mutate(start_dist=abs(refStart-end),end_dist=abs(start-refEnd)) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(min_dist_bases=min(start_dist,end_dist)) %>%
  dplyr::ungroup() %>% 
  dplyr::mutate(updown=ifelse(end < refStart,"upstream","downstream")) %>%
  dplyr::mutate(status=ifelse(min_dist_bases <2e4,"out_expand","outside"))

if (nrow(outreg_orthos  %>% dplyr::filter(min_dist_bases < 2e4)) > 0) {
  print("WARNING: There is at least one paralog that is within 10 kb of a gene within your defined boundary in N2. Your boundary will be automatically expanded to include:")
  print(outreg_orthos %>% dplyr::select(seqid,seqname,start,end,N2,min_dist_bases) %>% dplyr::filter(min_dist_bases < 2e4))
}

#generate a lookup table (all_ortho_pairs) which contains all pairwise gene orthologs between REF and WILD
orthoList <- list()
orthoList_bound <- list()
orthoList_raw <- list()
strainCol_iter <- strainCol_c2[!strainCol_c2 %in% c("Orthogroup","N2")]
strainCol_iter <- strainCol_iter[strainCol_iter %in% want] ############################# For when you are only looking at a subset of strains

for (i in 1:length(strainCol_iter)) {
  
  id=strainCol_iter[i]
  raw_tmp <- orthos %>%
    dplyr::select(Orthogroup,strainCol_iter[i],N2) %>%
    dplyr::mutate(str=!!sym(strainCol_iter[i])) %>%
    dplyr::mutate(str = strsplit(as.character(str), ",")) %>%
    tidyr::unnest(str) %>%
    dplyr::mutate(str=trimws(str)) %>%
    dplyr::filter(!is.na(N2)) %>%
    dplyr::select(-strainCol_iter[i]) %>%
    dplyr::mutate(STRAIN=strainCol_iter[i]) %>%
    dplyr::mutate(has_any_ortho=T) %>%
    dplyr::left_join(wild_tran,by=c("STRAIN","str"="tranname"))
  
  orthoList_raw[[i]] <- raw_tmp
  
  print(paste0("Mapped orthologs for ",i,"/",length(strainCol_iter)," strains."))
  tmp <- rbind(inreg_orthos %>% dplyr::mutate(status="within") %>% dplyr::select(Orthogroup,strainCol_iter[i],N2,seqid,seqname,start,end,og_loc,status),outreg_orthos %>% 
                 dplyr::select(Orthogroup,strainCol_iter[i],N2,seqid,seqname,start,end,og_loc,status)) %>%
    dplyr::select(Orthogroup,N2,strainCol_iter[i],og_loc,status) %>%
    dplyr::rename(tmpSel=strainCol_iter[i]) %>%
    dplyr::mutate(newSel = strsplit(as.character(tmpSel), ",")) %>%
    tidyr::unnest(newSel) %>%
    dplyr::mutate(newSel=trimws(newSel)) %>%
    dplyr::select(-tmpSel) %>%
    dplyr::mutate(STRAIN=strainCol_iter[i]) %>%
    tidyr::separate(newSel,into=c("Name","tnum"),sep="\\.",remove = F) %>%
    dplyr::select(Orthogroup,newSel,Name,STRAIN,N2,-tnum,og_loc,status) %>%
    dplyr::rename(tranname=newSel,Parent=Name) %>%
    dplyr::left_join(wild_tran,by=c("tranname","Parent","STRAIN")) 
  
  orthoList[[i]] <- tmp
  boundg <- tmp %>% 
    dplyr::filter(og_loc=="in_region" | status=="out_expand") %>%
    dplyr::select(-og_loc,-status) %>%
    dplyr::filter(seqid==boundChrom) %>%
    dplyr::mutate(og_loc=ifelse(((start >= boundStart & start <= boundEnd) | (end >= boundStart & end <= boundEnd)),"in_region","out_region")) %>%
    dplyr::mutate(start_dist=ifelse(og_loc=="out_region",abs(boundStart-end),NA),end_dist=ifelse(og_loc=="out_region",abs(start-boundEnd),NA)) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(min_dist_bases=min(start_dist,end_dist)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(status=ifelse(og_loc=="in_region","within",ifelse(min_dist_bases < 10000 & !is.na(min_dist_bases),"out_expand","outside")))
  check <- boundg %>% dplyr::filter(status=="out_expand")
  
  if (nrow(check) > 0) {
    print(paste0("WARNING: There is at least one paralog that is within 10 kb of a gene within your derived boundary in ",strainCol_iter[[i]],". Your boundary will be automatically expanded to include:"))
    print(check %>% dplyr::select(seqid,tranname,start,end,strand,STRAIN,min_dist_bases) %>% dplyr::distinct(tranname,.keep_all = T))
    
    sorter <- check %>% dplyr::mutate(updown=ifelse(min_dist_bases==start_dist,"upstream","downstream"))
    upstream <- sorter %>% dplyr::filter(updown=="upstream")
    downstream <- sorter %>% dplyr::filter(updown=="downstream")
    
    if (nrow(upstream) > 0) {
      outer_lim <- max(upstream$end)
      inner <- boundg %>% dplyr::filter(status=="within") %>% dplyr::arrange(start) %>% dplyr::filter(start==min(start))
      inner_lim <- min(inner$start)
      seqid_match <- as.character(unique(inner$seqid))
      extension <- raw_tmp %>%
        dplyr::filter(seqid==seqid_match & start > outer_lim & end <inner_lim) %>%
        dplyr::rename(tranname=str) %>%
        dplyr::select(Orthogroup,tranname,Parent,STRAIN,N2,everything(),-has_any_ortho) %>%
        dplyr::mutate(og_loc="out_region",status="out_extend") %>%
        dplyr::mutate(N2 = strsplit(as.character(N2), ",")) %>%
        tidyr::unnest(N2) %>%
        dplyr::mutate(N2=trimws(N2))
      
      boundg_inc <- rbind(extension, boundg %>% dplyr::select(-start_dist,-end_dist,-min_dist_bases))
      orthoList_bound[[i]]  <- boundg_inc %>% dplyr::arrange(start)
    } 
    
    if(nrow(downstream) > 0) {
      outer_lim <- min(downstream$end)
      inner <- boundg %>% dplyr::filter(status=="within") %>% dplyr::arrange(start) %>% dplyr::filter(end==max(end))
      inner_lim <- max(inner$end)
      seqid_match <- as.character(unique(inner$seqid))
      extension <- raw_tmp %>%
        dplyr::filter(seqid==seqid_match & start > inner_lim & end < outer_lim) %>%
        dplyr::rename(tranname=str) %>%
        dplyr::select(Orthogroup,tranname,Parent,STRAIN,N2,everything(),-has_any_ortho) %>%
        dplyr::mutate(og_loc="out_region",status="out_extend") %>%
        dplyr::mutate(N2 = strsplit(as.character(N2), ",")) %>%
        tidyr::unnest(N2) %>%
        dplyr::mutate(N2=trimws(N2))
      
      boundg_inc <- rbind(extension, boundg %>% dplyr::select(-start_dist,-end_dist,-min_dist_bases))
      orthoList_bound[[i]]  <- boundg_inc %>% dplyr::arrange(start)
    } 
    
    
  } else {
    orthoList_bound[[i]] <- boundg %>% 
      dplyr::select(-start_dist,-end_dist,-min_dist_bases) %>% dplyr::arrange(start)
  }
}

all_ortho_pairs  <- ldply(orthoList,data.frame) 
all_ortho_pairs_bound_pre <-ldply(orthoList_bound,data.frame) %>% 
  dplyr::filter(!status=="outside") 

corr_jumps <- all_ortho_pairs_bound_pre %>%
  dplyr::distinct(STRAIN,Parent,.keep_all = T) %>%
  dplyr::arrange(STRAIN,start) %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::mutate(leadDist=lead(start)-start) %>%
  dplyr::mutate(leadDist=ifelse(is.na(leadDist),0,leadDist)) %>%
  dplyr::mutate(jump=ifelse(lag(leadDist)>5e4 & lag(status)=="within","JUMP","NOJUMP")) %>%
  dplyr::mutate(jump=ifelse(is.na(jump),"NOJUMP",jump)) %>%
  dplyr::mutate(jumpID=rleid(jump)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(STRAIN,jumpID) %>%
  dplyr::mutate(jgroup_size=n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::filter(jgroup_size==max(jgroup_size)) %>%
  dplyr::mutate(keep=T)


all_ortho_pairs_bound <- all_ortho_pairs_bound_pre %>%
  dplyr::arrange(STRAIN,start)

all_ortho_pairs_raw <- ldply(orthoList_raw,data.frame) %>% dplyr::select(STRAIN,str,has_any_ortho) %>% dplyr::rename(tranname=str) 

new_boundaries_WI <-  all_ortho_pairs_bound %>%
  dplyr::select(seqid,start,end,STRAIN,tranname,Parent) %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::mutate(minStart=min(start), maxEnd=max(end)) %>%
  dplyr::distinct(tranname,.keep_all = T) %>%
  dplyr::filter(start==minStart | end==maxEnd) %>%
  dplyr::mutate(gene2gene=paste(Parent,collapse="-")) %>%
  dplyr::distinct(minStart,.keep_all = T) %>%
  dplyr::ungroup() %>%
  dplyr::select(seqid,minStart,maxEnd,STRAIN,gene2gene)

N2_expand <- rbind(inreg_orthos %>% 
                     dplyr::mutate(status="within"),outreg_orthos %>% 
                     dplyr::select(-refStart,-refEnd,-refChrom,-start_dist,-end_dist,-min_dist_bases,-updown)) %>%
  dplyr::filter(!status=="outside")

new_boundaries_N2 <- N2_expand %>%
  dplyr::mutate(minStart=min(start), maxEnd=max(end)) %>%
  dplyr::distinct(N2,.keep_all = T) %>% 
  dplyr::filter(start==minStart | end==maxEnd) %>%
  dplyr::mutate(gene2gene=paste(seqname,collapse="-")) %>%
  dplyr::mutate(STRAIN="N2") %>%
  dplyr::distinct(minStart,.keep_all = T) %>%
  dplyr::select(seqid,minStart,maxEnd,STRAIN,gene2gene)

new_boundaries <- rbind(new_boundaries_WI,new_boundaries_N2) %>%
  dplyr::rename(boundStart=minStart,boundEnd=maxEnd)

#find the bound genes for each strain that are not orthologous
boundGenes <- rbind(wild_tran %>% 
                      dplyr::select(-boundChrom,-boundStart,-boundEnd,-REF,-refStart,-refEnd,-inv) %>% 
                      dplyr::mutate(alias=NA),
                    N2_tran %>% dplyr::select(tranname,seqname,STRAIN,seqid,start,end,strand,alias) %>%
                      dplyr::rename(Parent=seqname)) %>%
  dplyr::left_join(new_boundaries,by=c("STRAIN","seqid")) %>%
  dplyr::filter(!is.na(boundStart)) %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::filter(start >= boundStart & start <= boundEnd) %>%
  dplyr::ungroup()


N2_ad <- boundGenes %>% 
  dplyr::filter(STRAIN=="N2") %>%
  dplyr::mutate(tr_has_any_ortho=ifelse(tranname %in% all_orthos_unnest$N2,T,F)) %>%
  dplyr::mutate(tr_has_bound_ortho=ifelse(tranname %in% all_ortho_pairs_bound$N2,T,F)) %>%
  dplyr::arrange(start) %>%
  dplyr::group_by(Parent) %>%
  dplyr::mutate(has_any_ortho = any(tr_has_any_ortho)) %>% 
  dplyr::mutate(has_bound_ortho = any(tr_has_bound_ortho)) %>%
  dplyr::select(-tr_has_any_ortho,-tr_has_bound_ortho) %>%
  dplyr::ungroup()

reassess_distal <- N2_ad %>% dplyr::filter(has_any_ortho==T & has_bound_ortho==F) 

orthos_tran <- orthos %>%
  dplyr::mutate(N2 = strsplit(as.character(N2), ",")) %>%
  tidyr::unnest(N2) %>%
  dplyr::mutate(N2=trimws(N2)) %>%
  dplyr::left_join(N2_tran_reg %>% dplyr::select(tranname,seqid,seqname,start,end) %>% dplyr::mutate(og_loc="in_region"),by=c("N2"="tranname"))


distal_ortho <- orthos_tran %>%  # inreg_orthos
  dplyr::filter(N2 %in% reassess_distal$tranname) %>%
  dplyr::select(any_of(strainCol_iter), N2)%>%
  # dplyr::mutate(comma_count = stringr::str_count(CB4856, ",")+1) %>%
  # dplyr::group_by(CB4856) %>%
  # dplyr::mutate(comma_count=sum(comma_count)) %>%
  # dplyr::filter(comma_count > 1)
  dplyr::mutate(has_any_in_want = dplyr::if_any(dplyr::all_of(strainCol_iter), ~ !is.na(.) & . != "")) %>%
  dplyr::filter(!has_any_in_want) %>%     # keep rows where NONE of the strains has an ortholog
  dplyr::distinct(N2)

N2_ad_corr <- N2_ad %>%
  dplyr::mutate(has_any_ortho=ifelse(tranname %in% distal_ortho$N2,F,has_any_ortho)) 

g_count <- length(unique(N2_ad_corr$Parent))


# Ordering the strains by N2, REF strains, ALT strains
desired_order <- rev(want)

WI_ad <- boundGenes %>% 
  dplyr::filter(!STRAIN=="N2") %>%
  dplyr::left_join(all_ortho_pairs_raw,by=c("STRAIN","tranname")) %>%
  dplyr::mutate(tr_has_any_ortho=ifelse(is.na(has_any_ortho),F,has_any_ortho)) %>%
  dplyr::left_join(all_ortho_pairs_bound %>% dplyr::select(tranname,STRAIN,N2,status) %>% dplyr::filter(N2 %in% N2_ad$tranname),by=c("STRAIN","tranname")) %>%
  dplyr::mutate(tr_has_bound_ortho=ifelse(!is.na(status),T,F)) %>%
  dplyr::select(-status,-has_any_ortho) %>% 
  dplyr::rename(N2_name=N2) %>%
  dplyr::left_join(aliases %>% dplyr::select(-seqname),by=c("N2_name"="tranname")) %>%
  dplyr::mutate(alias.x=alias.y) %>%
  dplyr::select(-alias.y,-N2_name) %>%
  dplyr::rename(alias=alias.x) %>%
  dplyr::group_by(STRAIN,Parent) %>%
  dplyr::mutate(has_any_ortho = any(tr_has_any_ortho)) %>% 
  dplyr::mutate(has_bound_ortho = any(tr_has_bound_ortho)) %>%
  dplyr::select(-tr_has_any_ortho,-tr_has_bound_ortho) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(HV_boundary %>% dplyr::select(STRAIN,inv),by="STRAIN") %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::mutate(n_gene=n_distinct(Parent)) %>%
  dplyr::mutate(start_sort = if_else(rep(dplyr::first(inv), dplyr::n()), -start, start)) %>%
  dplyr::arrange(start_sort, .by_group = TRUE) %>%
  dplyr::mutate(first_gene=dplyr::first(alias)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(STRAIN = factor(STRAIN, levels = desired_order)) %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::arrange(first_gene, n_gene, .by_group = TRUE) %>%
  dplyr::mutate(order_gene = dplyr::first(first_gene), order_num = dplyr::first(n_gene)) %>%
  dplyr::ungroup() %>%
  # dplyr::arrange(desc(order_gene), desc(order_num)) %>%
  # dplyr::arrange(STRAIN, desc(order_gene), order_num) %>%
  # dplyr::arrange(desc(STRAIN)) %>%
  dplyr::select(-order_gene, -order_num) %>%
  dplyr::mutate(g_diff = abs(n_gene-g_count)) %>%
  dplyr::mutate(y_pos=rleid(STRAIN)) %>%
  dplyr::select(-g_diff,-start_sort,-inv,-first_gene,-n_gene) 

N2_ad_corr <- N2_ad_corr %>%
  dplyr::mutate(y_pos=max(WI_ad$y_pos)+1)


all_ad <- rbind(N2_ad_corr,WI_ad) %>% 
  dplyr::arrange(STRAIN,start) %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::mutate(gen_pos=rleid(Parent)) %>%
  dplyr::mutate(shift=min(start)) %>%
  dplyr::mutate(end=end-min(start),start=start-min(start)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(col=ifelse(has_any_ortho==T & has_bound_ortho ==T,2,ifelse(has_any_ortho==T,1,0))) %>%
  dplyr::left_join(HV_boundary %>% dplyr::select(STRAIN,inv),by="STRAIN") %>%
  dplyr::mutate(inv=ifelse(is.na(inv),F,inv)) %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::mutate(bound_corr=max(end)) %>%
  dplyr::mutate(start=ifelse(inv==T,abs(start-bound_corr),start)) %>%
  dplyr::mutate(end=ifelse(inv==T,abs(end-bound_corr),end))


hlines <- new_boundaries %>% 
  dplyr::left_join(all_ad %>% dplyr::select(STRAIN,y_pos) %>% dplyr::distinct(STRAIN,.keep_all = T),by="STRAIN") %>%
  dplyr::left_join(all_ad %>% dplyr::select(STRAIN,shift) %>% dplyr::distinct(STRAIN,.keep_all = T),by="STRAIN") #%>%
  # dplyr::mutate(geno = ifelse(STRAIN %in% alt_strains[alt_strains != "N2"], "ALT", "REF")) %>%
  # dplyr::mutate(geno = factor(geno, levels = c("REF","ALT"))) 

segments <- all_ortho_pairs_bound %>%
  dplyr::select(STRAIN,Parent,start,end,N2,strand) %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::distinct(Parent,N2,.keep_all = T) %>%
  dplyr::left_join(N2_tran %>% 
                     dplyr::rename(start_N2=start,end_N2=end,strand_N2=strand,chrom_N2=seqid,N2id=STRAIN) %>% 
                     dplyr::select(tranname,chrom_N2,start_N2,end_N2,strand_N2,alias,seqname,N2id),
                   by=c("N2"="tranname")) %>% 
  dplyr::filter(N2 %in% N2_ad$tranname) %>%
  dplyr::left_join(all_ad %>% dplyr::select(STRAIN,y_pos) %>% dplyr::distinct(STRAIN,.keep_all = T),by="STRAIN") %>%
  dplyr::rename(WI_y_pos=y_pos) %>%
  dplyr::left_join(all_ad %>% dplyr::select(STRAIN,y_pos) %>% dplyr::distinct(STRAIN,.keep_all = T),by=c("N2id"="STRAIN")) %>%
  dplyr::rename(N2_y_pos=y_pos) %>%
  dplyr::mutate(N2_shift=min(start_N2),WI_shift=min(start)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(WI_x_pos=(start+((end-start)/2))-WI_shift, N2_x_pos=(start_N2+((end_N2-start_N2)/2)-N2_shift)) %>%
  dplyr::mutate(WI_y_pos=WI_y_pos+0.2,N2_y_pos=N2_y_pos-0.2) %>%
  dplyr::distinct(STRAIN,Parent,seqname,.keep_all = T) %>%
  dplyr::group_by(STRAIN,Parent) %>%
  dplyr::mutate(n1=n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(STRAIN,seqname) %>%
  dplyr::mutate(n2=n()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(col=ifelse(n1>1 | n2>1,"multi_copy","single_copy")) 

plot_ad <- all_ad %>%
  dplyr::group_by(STRAIN,Parent) %>%
  dplyr::filter(col==max(col)) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(STRAIN,Parent,.keep_all = T) %>%
  dplyr::mutate(class=ifelse(col==0,"no_known_ortho",ifelse(col==1,"has_distal_ortho","has_local_ortho"))) 

all_hap <- ggplot() +
  geom_segment(data=hlines,aes(x=boundStart-shift,xend=boundEnd-shift,y=y_pos,yend=y_pos))+
  geom_rect(data=plot_ad, aes(xmin=start,xmax=end,ymin=y_pos+0.2,ymax=y_pos-0.2,fill=class),color="black") +
  theme(legend.title = element_blank(),
        panel.background = element_blank(),
        axis.title = element_blank(),
        axis.text = element_text(size = 12, color = 'black'),
        axis.ticks = element_blank()) +
  scale_fill_manual(values=c("has_local_ortho"="grey","has_distal_ortho"="black","no_known_ortho"="red","no_known_allelic_CB"="blue")) +
  scale_x_continuous(expand = c(0.01,0)) +
  scale_y_continuous(expand = c(0.01, 0), breaks = hlines$y_pos, labels = hlines$STRAIN)
all_hap

# Exclude "non-ortho" from the trapezium joining
plot_ad_filtered <- plot_ad %>%
  dplyr::mutate(alias=ifelse(is.na(alias),"non-ortho",alias)) %>%
  filter(alias != "non-ortho")

# Join filtered data frames for many-to-many connections
trapeziums <- dplyr::inner_join(
  plot_ad_filtered, plot_ad_filtered,
  by = "alias",
  suffix = c("_upper", "_lower"),
  relationship = "many-to-many") %>% 
  dplyr::filter(y_pos_upper - y_pos_lower == 1)

# Create trapezium polygons using min/max for x-coordinates so that start/end orientation is corrected.
trapezium_polys <- trapeziums %>% 
  rowwise() %>%
  do({
    # Calculate corrected x coordinates for the upper rectangle
    x_left_upper <- min(.$start_upper, .$end_upper)
    x_right_upper <- max(.$start_upper, .$end_upper)
    
    # Calculate corrected x coordinates for the lower rectangle
    x_left_lower <- min(.$start_lower, .$end_lower)
    x_right_lower <- max(.$start_lower, .$end_lower)
    
    data.frame(
      alias = .$alias,
      strain_lower = .$STRAIN_lower,
      group = paste(.$alias, .$y_pos_upper, sep = "_"),
      x = c(x_left_upper, x_right_upper, x_right_lower, x_left_lower),
      y = c(.$y_pos_upper - 0.2,  # bottom edge of the upper rectangle
            .$y_pos_upper - 0.2,
            .$y_pos_lower + 0.2,  # top edge of the lower rectangle
            .$y_pos_lower + 0.2)
    )
  }) %>%
  dplyr::ungroup()# %>%
  # dplyr::mutate(geno = ifelse(strain_lower %in% alt_strains[alt_strains != "N2"], "ALT", "REF")) %>%
  # dplyr::mutate(geno = factor(geno, levels = c("REF","ALT")))

# Extract unique aliases at y_pos 77 in order of increasing start position
ordered_aliases <- plot_ad %>%
  dplyr::filter(y_pos == max(plot_ad$y_pos)) %>%
  dplyr::arrange(start) %>%
  dplyr::pull(alias) %>%
  unique()

# Reorder the factor levels so that the legend follows the ordered aliases
plot_ad <- plot_ad %>%
  dplyr::mutate(alias = factor(alias, levels = ordered_aliases)) #%>%
  # dplyr::mutate(geno = ifelse(STRAIN %in% alt_strains[alt_strains != "N2"], "ALT", "REF")) %>%
  # dplyr::mutate(geno = factor(geno, levels = c("REF","ALT"))) 

# Also update any other data frames with alias info, e.g. trapezium_polys:
trapezium_polys <- trapezium_polys %>%
  dplyr::mutate(alias = factor(alias, levels = ordered_aliases))

# Shuffle the assignment of colors to the ordered aliases
set.seed(9) 
shuffled_aliases <- sample(ordered_aliases)

# Generate colors using hcl.colors() for the shuffled aliases
default_colors <- setNames(hcl.colors(length(shuffled_aliases), "Dark 3"), shuffled_aliases)

# But to keep the legend order as ordered_aliases, we re-map these colors back:
final_colors <- default_colors[ordered_aliases]

# Optionally, if you have the "non-ortho" alias (or any other), add it explicitly:
final_colors <- c(final_colors, "non-N2 gene" = "darkgrey")

# plot_ad_segments <- plot_ad %>%
#   dplyr::mutate(
#     # Adjust strand logic if inverted
#     strand_logic = case_when(
#       strand == "+" & !inv ~ "+",
#       strand == "-" & !inv ~ "-",
#       strand == "+" & inv  ~ "-",
#       strand == "-" & inv  ~ "+"),
#     seg_color = ifelse(strand_logic == "+", "black", "red"),
#     
#     x_start = start,
#     x_end   = end,
#     y_seg   = y_pos - 0.25  # just under the geom_rect (geom_rect is y_pos ± 0.2))

# Create final plot
all_hap_bg <- ggplot() +
  geom_segment(data = hlines, 
               aes(x = boundStart - shift, xend = boundEnd - shift, y = y_pos, yend = y_pos)) +
  geom_polygon(data = trapezium_polys,
               aes(x = x, y = y, group = group, fill = alias)) +
  geom_rect(data = plot_ad %>% dplyr::mutate(alias=ifelse(is.na(alias),"non-N2 gene", as.character(alias))),
            aes(xmin = start, xmax = end, ymin = y_pos + 0.2, ymax = y_pos - 0.2, fill = alias), color = "black") +
  annotate("rect", xmin = -1000, xmax = -100, ymin = 0.7, ymax = 16.3, fill = 'red') +
  annotate("rect", xmin = -1000, xmax = -100, ymin = 16.7, ymax = 50.3, fill = 'black') +
  scale_y_continuous(expand = c(0.01, 0), breaks = hlines$y_pos, labels = hlines$STRAIN) +
  # annotate("text", x = 14180, y = 49.8, label = "*", size = 5, color = "black") +
  scale_x_continuous(expand = c(0.01, 0), labels = function(x) x / 1000) +
  scale_fill_manual(values = final_colors, breaks = names(final_colors)) +
  scale_color_identity()  +
  # ggh4x::facet_grid2(geno ~ ., scales = "free_y", space = "free_y", switch = "y",
                     # strip = ggh4x::strip_themed(
                       # background_y = list(
                         # element_rect(fill = "black", color = "black"),
                         # element_rect(fill = "red", color = "black")),
                       # text_y = list(
                         # element_text(angle = 90, size = 14, color = "white", face = "bold"),
                         # element_text(angle = 90, size = 14, color = "white", face = "bold")))) +
  labs(fill="Reference\ngene") +
  xlab("Physical distance (kb)") +
  theme(
    panel.background = element_blank(),
    axis.title = element_blank(),
    axis.text = element_text(size = 11, color = 'black'), 
    axis.ticks.y = element_blank(),
    axis.line.x = element_line(),
    panel.spacing.y = unit(0, "pt"),
    axis.title.x = element_text(color = 'black', size  = 14),
    # legend.position = 'none',
    legend.position = "right",
    legend.direction = "horizontal",
    legend.key.size = unit(0.4, "lines"),
    legend.text = element_text(size = 8, color = 'black', face = 'italic'),
    legend.title = element_text(size = 10, color = 'black')) +
  guides(fill = guide_legend(title.position = "top", nrow = 27, byrow = TRUE, override.aes = list(size = 2)))
all_hap_bg

# Save plot:
ggsave("../../figures/supplementary/abamectin_REF_ALT_haplotypes.png", all_hap_bg, width = 7.5, height = 7, dpi = 600)


# Ensure that none of the 5 missing N2 genes are in any of the ALT strains
alt_strains <- tail(want,17) 
n2_pav <- c("F11A5.18", "F11A5.7", "F11A5.8", "F11A5.16", "F11A5.15")

alt_orthos <- orthos %>% dplyr::select(Orthogroup, all_of(alt_strains)) %>%
  dplyr::mutate(N2 = gsub("transcript_","",N2)) %>%
  tidyr::separate_rows(N2, sep = ", ") %>%
  dplyr::mutate(N2 = sub("\\.[^.]*$", "", N2)) %>%
  dplyr::filter(N2 %in% n2_pav) # CONFIRMED! All five N2 genes are missing orthologs in ALT strains (except for ECA36)



# # Visualizing gene trees for the 5 genes present in ECA36 but absent in all other ALT strains
# # F11A5.18 
# tree <- ape::read.tree("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/orthology/elegans/orthofinder/64_core/OrthoFinder/Results_Dec07/Resolved_Gene_Trees/F11A5.18_geneTree.nwk")
# 
# tree$tip.label <- stringr::str_remove(tree$tip.label, "_20251012.*|_Nov2025.*|_PRJNA13758.*|_20251014.*|_20251124.*")
# 
# tree$tip.label <- ifelse(tree$tip.label == "c_elegans", "N2", tree$tip.label)
# 
# ref_strains <- head(want,33) %>% c("N2")  
# 
# tip_colors <- ifelse(tree$tip.label %in% "ECA36", "red", 
#                      ifelse(tree$tip.label %in% ref_strains, "black", "grey"))
# 
# plot(tree, cex = 0.7, font = 1, tip.color = tip_colors, label.offset = 0.01)
# title(main = "Gene Tree for F11A5.18", cex.main = 1.5, cex.sub = 1.0, col.main = "black", font.main = 2)
# length(tree$tip.label) # only 104 of the 142 have this gene
# 
# 
# # F11A5.7
# tree <- ape::read.tree("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/orthology/elegans/orthofinder/64_core/OrthoFinder/Results_Dec07/Resolved_Gene_Trees/F11A5.7_geneTree.nwk")
# 
# tree$tip.label <- stringr::str_remove(tree$tip.label, "_20251012.*|_Nov2025.*|_PRJNA13758.*|_20251014.*|_20251124.*")
# 
# tree$tip.label <- ifelse(tree$tip.label == "c_elegans", "N2", tree$tip.label)
# 
# ref_strains <- head(want,33) %>% c("N2")  
# 
# tip_colors <- ifelse(tree$tip.label %in% "ECA36", "red", 
#                      ifelse(tree$tip.label %in% ref_strains, "black", "grey"))
# 
# plot(tree, cex = 0.7, font = 1, tip.color = tip_colors, label.offset = 0.01)
# title(main = "Gene Tree for irld-26 (F11A5.7)", cex.main = 1.5, cex.sub = 1.0, col.main = "black", font.main = 2)
# length(tree$tip.label) # only 104 of the 142 have this gene
# 
# 
# # F11A5.8
# tree <- ape::read.tree("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/orthology/elegans/orthofinder/64_core/OrthoFinder/Results_Dec07/Resolved_Gene_Trees/F11A5.8_geneTree.nwk")
# 
# tree$tip.label <- stringr::str_remove(tree$tip.label, "_20251012.*|_Nov2025.*|_PRJNA13758.*|_20251014.*|_20251124.*")
# 
# tree$tip.label <- ifelse(tree$tip.label == "c_elegans", "N2", tree$tip.label)
# 
# ref_strains <- head(want,33) %>% c("N2")  
# 
# tip_colors <- ifelse(tree$tip.label %in% "ECA36", "red", 
#                      ifelse(tree$tip.label %in% ref_strains, "black", "grey"))
# 
# plot(tree, cex = 0.7, font = 1, tip.color = tip_colors, label.offset = 0.01)
# title(main = "Gene Tree for oac-15 (F11A5.8)", cex.main = 1.5, cex.sub = 1.0, col.main = "black", font.main = 2)
# 
# length(tree$tip.label) # only 104 of the 142 have this gene
# 
# 
# # F11A5.16
# tree <- ape::read.tree("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/orthology/elegans/orthofinder/64_core/OrthoFinder/Results_Dec07/Resolved_Gene_Trees/F11A5.16_geneTree.nwk")
# 
# tree$tip.label <- stringr::str_remove(tree$tip.label, "_20251012.*|_Nov2025.*|_PRJNA13758.*|_20251014.*|_20251124.*")
# 
# tree$tip.label <- ifelse(tree$tip.label == "c_elegans", "N2", tree$tip.label)
# 
# ref_strains <- head(want,33) %>% c("N2")  
# 
# tip_colors <- ifelse(tree$tip.label %in% "ECA36", "red", 
#                      ifelse(tree$tip.label %in% ref_strains, "black", "grey"))
# 
# plot(tree, cex = 0.7, font = 1, tip.color = tip_colors, label.offset = 0.01)
# title(main = "Gene Tree for F11A5.16", cex.main = 1.5, cex.sub = 1.0, col.main = "black", font.main = 2)
# 
# length(tree$tip.label) # only 102 of the 142 have this gene
# 
# 
# 
# # F11A5.15
# tree <- ape::read.tree("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/orthology/elegans/orthofinder/64_core/OrthoFinder/Results_Dec07/Resolved_Gene_Trees/F11A5.15_geneTree.nwk")
# 
# tree$tip.label <- stringr::str_remove(tree$tip.label, "_20251012.*|_Nov2025.*|_PRJNA13758.*|_20251014.*|_20251124.*")
# 
# tree$tip.label <- ifelse(tree$tip.label == "c_elegans", "N2", tree$tip.label)
# 
# ref_strains <- head(want,33) %>% c("N2")  
# 
# tip_colors <- ifelse(tree$tip.label %in% "ECA36", "red", 
#                      ifelse(tree$tip.label %in% ref_strains, "black", "grey"))
# 
# plot(tree, cex = 0.7, font = 1, tip.color = tip_colors, label.offset = 0.01)
# title(main = "Gene Tree for F11A5.15", cex.main = 1.5, cex.sub = 1.0, col.main = "black", font.main = 2)
# 
# length(tree$tip.label) # only 99 of the 142 have this gene


























################ ADD PSEUDOGENES FOR WILD STRAINS (CHECK TO SEE IF THEY ARE EXPRESSED) AND ADD OVERLAPPING N2 GENES AS DIFFERENT COLUMNS - EXPRESSION EVIDENCE IN WILD STRAINS???
# Make all genes the same size and in the exact same position for orthologs across all strains
updated_gene_sizes <- plot_ad %>% dplyr::filter(STRAIN == "N2") %>% dplyr::select(alias, start, end) %>% 
  dplyr::filter(!is.na(alias)) %>% # ADD BACK SRZ-13 AND SRZ-14 - they both have distal orthologs in only some strains...
  dplyr::arrange(start) %>%
  dplyr::mutate(
    new_start = cumsum(c(0, rep(1000, n() - 1))),
    new_end = new_start + 500) # making genes evenly spaces

plot_ad_new <- plot_ad %>% 
  dplyr::filter(!is.na(alias)) %>% # For visualization purposes
  dplyr::left_join(updated_gene_sizes, by = "alias") %>%
  dplyr::mutate(geno = ifelse(STRAIN %in% alt_strains[alt_strains != "N2"], "ALT", "REF")) %>%
  dplyr::mutate(geno = factor(geno, levels = c("REF","ALT"))) %>%
  dplyr::select(STRAIN, alias, start = start.x, end = end.x, new_start, new_end, geno, y_pos) 

hlines_new <- plot_ad_new %>% dplyr::select(STRAIN, new_start, new_end, y_pos) %>% 
  dplyr::mutate(start = min(new_start, na.rm = TRUE), 
                end = max(new_end, na.rm = TRUE)) %>%
  dplyr::distinct(STRAIN, start, end, y_pos) %>%
  dplyr::mutate(geno = ifelse(STRAIN %in% alt_strains[alt_strains != "N2"], "ALT", "REF")) %>%
  dplyr::mutate(geno = factor(geno, levels = c("REF","ALT")))

gene_labels <- plot_ad_new %>% dplyr::mutate(middle = (new_start + new_end) / 2) %>%
  dplyr::distinct(alias, middle) %>%
  dplyr::mutate(geno = factor("REF", levels = c("REF", "ALT")))



# Filling genes by %IDY with N2 
n2_genes_interval <- plot_ad %>% dplyr::filter(STRAIN == "N2") %>% dplyr::distinct(tranname) %>% dplyr::filter(tranname != "transcript_K05D4.6b.1" & tranname != "transcript_F11A5.5a.1") %>% 
  dplyr::pull() %>% c("transcript_K05D4.6a.1") %>% c("transcript_F11A5.5b.1")

orthos_interval <- orthos %>% dplyr::select(Orthogroup, all_of(want)) %>%
  tidyr::separate_rows(N2, sep = ", ") %>%
  dplyr::filter(N2 %in% n2_genes_interval) %>%
  dplyr::distinct(Orthogroup)

# write.table(orthos_interval, "../../processed_data/gwas/interval_orthogroups.tsv", quote = F, col.names = F, row.names = F)

# write.table(data.frame(strains = want), "../../processed_data/gwas/strains.tsv", sep = "\t", row.names = FALSE, col.names = F, quote = FALSE)



# OG MSA files
ogs <- list.files("../../processed_data/gwas", pattern = "^OG.*\\.MSA.fa$", full.names = TRUE)

# Calculate percent IDY between each wild strain ortholog and N2
msa_idy <- lapply(ogs, function(file) {
  
  og <- tools::file_path_sans_ext(basename(file))
  
  aln <- readAAStringSet(file)
  
  # Extract strain name from sequence IDs
  strain <- sub("_.*$", "", names(aln))
  
  # Keep desired strains
  aln <- aln[strain %in% want]
  strain <- strain[strain %in% want]
  
  # Find N2
  n2 <- aln[strain == "N2"]
  
  identity <- sapply(setdiff(unique(strain), "N2"), function(s) {
    
    seq <- aln[strain == s][1]
    
    n2_chars <- strsplit(as.character(n2), "")[[1]]
    seq_chars <- strsplit(as.character(seq), "")[[1]]
    
    # Check alignment length
    if (length(n2_chars) != length(seq_chars)) {
      warning(
        og, ": ", s,
        " has alignment length ", length(seq_chars),
        " vs N2 length ", length(n2_chars)
      )
      return(NA_real_)
    }
    
    # Only compare positions where neither sequence has a gap
    keep <- n2_chars != "-" & seq_chars != "-"
    
    mean(n2_chars[keep] == seq_chars[keep]) * 100
  })
  
  # Add N2 = 100
  identity <- c(N2 = 100, identity)
  
  data.frame(orthogroup = og, t(identity), check.names = FALSE)
}) %>% dplyr::bind_rows()


# Pivot for left merging and plotting
og_idy_long <- msa_idy %>% tidyr::pivot_longer(-orthogroup, names_to = "strain", values_to = "identity") %>%
  dplyr::mutate(orthogroup = gsub("^[^_]+_", "", orthogroup),
                orthogroup = gsub(".MSA", "", orthogroup)) %>%
  dplyr::rename(tranname = orthogroup) %>%
  dplyr::left_join(N2_tran_reg, by = "tranname") %>%
  dplyr::distinct(STRAIN = strain, identity, alias)


plot_ad_new_idy <- plot_ad_new %>% dplyr::left_join(og_idy_long, by = c("STRAIN", "alias")) %>% dplyr::mutate(identity = as.numeric(as.character(identity)))


# Create the final aligned plot!
all_hap_bg_new <- ggplot() +
  geom_segment(data = hlines_new,
               aes(x = start, xend = end, y = y_pos, yend = y_pos)) +
  geom_rect(data = plot_ad_new_idy %>% dplyr::mutate(alias=ifelse(is.na(alias),"Unknown gene",as.character(alias))),
            aes(xmin = new_start, xmax = new_end, ymin = y_pos + 0.4, ymax = y_pos - 0.4, fill = identity),color = "black") +
  scale_y_continuous(expand = c(0.01, 0), breaks = hlines$y_pos, labels = hlines$STRAIN) +
  scale_x_continuous(expand = c(0.01, 0),labels = function(x) x / 1000) +
  scale_fill_gradientn(colors = c("yellow", "blue")) +
  # scale_color_identity()  +
  labs(fill="Percent\nidentity (%)") +
  ggh4x::facet_grid2(geno ~ ., scales = "free_y", space = "free_y", switch = "y",
    strip = ggh4x::strip_themed(
      background_y = list(
        element_rect(fill = "black", color = "black"),
        element_rect(fill = "red", color = "black")),
      text_y = list(
        element_text(angle = 90, size = 10, color = "white", face = "bold"),
        element_text(angle = 90, size = 10, color = "white", face = "bold")))) +
  # geom_text(data = gene_labels, aes(x = middle, y = 52, label = alias), angle = 0, size = 6, hjust = 0.5, vjust = 0.5) +
  theme(
    panel.background = element_blank(),
    axis.title = element_blank(),
    axis.text.y = element_text(size = 8, color = 'black'), 
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks = element_blank(),
    axis.line.x = element_blank(),
    plot.margin = margin(l = 5,0,0,0),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.box.spacing = unit(0, "pt"),
    legend.title = element_text(size = 8, color = 'black'),
    legend.text = element_text(size = 8, color = 'black'),
    # legend.position = 'none',
    # strip.placement = "outside",           
    strip.background = element_rect(fill = "grey85", color = "black"), 
    strip.text.y.left = element_text(angle = 90, size = 8)) +
  guides(fill = guide_colorbar(title.position = "left", title.theme = element_text(margin = margin(r = 10)), byrow = TRUE, override.aes = list(size = 9)))
all_hap_bg_new


label_plot <- ggplot(gene_labels) +
  geom_text(aes(x = middle, y = 0, label = alias), angle = 50, size = 2.5, hjust = 0.4, vjust = 0, fontface = "italic") +
  scale_x_continuous(expand = c(0.01, 0), limits = range(c(hlines_new$start, hlines_new$end))) +
  scale_y_continuous(expand = c(0, 1), limits = c(0, 1)) +
  theme_void()

# Combine using cowplot
final_labeled_plt <- cowplot::plot_grid(
  label_plot,
  all_hap_bg_new,
  ncol = 1,
  rel_heights = c(0.05, 0.85),
  align = "v",
  axis = "lr")  + theme(plot.background = element_rect(fill = "white", color = NA))
final_labeled_plt


# Save the plot:
ggsave("../../figures/abamectin_haplotype.png", width = 7.5, height = 7.5, dpi = 600)





# Tile plot of SNVs
interval_snvs <- readr::read_tsv("../../processed_data/gwas/SNVs_in_Interval.tsv") %>%
  dplyr::mutate(dplyr::across(4:last_col(), ~ dplyr::case_when(
    .x == "0/0" ~ "0",
    .x == "1/1" ~ "1",
    .x == "./." ~ NA_character_,
    TRUE ~ .x))) %>%
  dplyr::select(-CHROM,-REF,-ALT, -N2) 

interval_snv_nonMissing <- interval_snvs %>% dplyr::filter(dplyr::if_all(4:last_col(), ~ !is.na(.x))) %>%
  dplyr::select(-CHROM,-REF,-ALT, -N2) %>%
  dplyr::mutate(alt_specific = dplyr::if_all(
      dplyr::all_of(alt_strains[alt_strains != "N2"]), ~ .x == 1) &
      dplyr::if_all(dplyr::all_of(setdiff(want, alt_strains)), ~ .x == 0))

alt_specific <- interval_snv_nonMissing %>% dplyr::select(POS, alt_specific) %>% dplyr::filter(alt_specific == "TRUE") %>%
  dplyr::mutate(seqid = "V")

pos <- interval_snv_nonMissing$POS

snv_nonMissing_long <- interval_snv_nonMissing %>%
  dplyr::select(-POS, -alt_specific) %>%
  t() %>%
  as.data.frame() %>%
  setNames(pos) %>%
  tibble::rownames_to_column("strain") %>%
  tidyr::pivot_longer(
    cols = -strain,
    names_to = "POS",
    values_to = "genotype") %>%
  dplyr::mutate(POS = as.numeric(POS)) %>%
  dplyr::mutate(strain = factor(strain, levels = rev(want[want != "N2"]))) %>%
  dplyr::mutate(genotype = ifelse(genotype == "0", "REF", "ALT"))

n2_genes_ROI <- gffCat2 %>% dplyr::filter(seqid == "V", start >= (16198034 - 22000) & end <=(16198034 + 32000), type == "gene") 

genes_interval <- ggplot(n2_genes_ROI) +
  geom_hline(yintercept = 1) +
  geom_rect(aes(xmin = start / 1e6, xmax = end / 1e6, ymin = 0.9, ymax = 1.1)) +
  geom_vline(xintercept = 16198034 / 1e6, color = 'red', linetype = 'solid') +
  geom_rect(data = alt_specific, aes(xmin = (POS - 10) / 1e6, xmax = (POS + 10) / 1e6, ymin = 0.8, ymax = 1.2), fill = "magenta3") +
  scale_y_continuous(expand = c(0,0)) +
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())
genes_interval 

# Adding gene labels
n2_gene_name_pos <- N2_tran %>% dplyr::left_join(n2_genes_ROI, by = c("seqid","start","end")) %>% dplyr::filter(!is.na(STRAIN.y)) %>%
  dplyr::filter(tranname %in% n2_genes_interval) %>%
  dplyr::mutate(middle = (start + end) / 2) %>%
  dplyr::select(alias, middle)

label_plot_2 <- ggplot(n2_gene_name_pos) +
  geom_text(aes(x = middle, y = 0, label = alias), angle = 75, size = 2, hjust = 0.5, vjust = 0.5) +
  scale_x_continuous(expand = c(0.08,0.08)) +
  scale_y_continuous(expand = c(0,0)) +
  theme_void() +
  theme(plot.margin = margin(l = 5, r = 5, t = 5, b = -27))
# label_plot_2

snvs_interval_plt <- ggplot(snv_nonMissing_long %>% dplyr::filter(POS >= min(n2_genes_ROI$start) & POS <= max(n2_genes_ROI$end)), 
                            aes(x = POS / 1e6, y = strain, fill = genotype)) +
  geom_tile() +
  scale_fill_manual(values = c("REF" = "white", "ALT" = "red")) +
  scale_x_continuous(labels = scales::comma) +
  theme_minimal() +
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 11, color = 'black'),
    axis.text.x = element_text(size = 10, color = 'black'),
    axis.text.y = element_text(size = 9, color = 'black'),
    panel.border = element_rect(color = "black", fill = NA),
    panel.grid = element_blank()) +
  labs(x = "N2 genomic position (Mb)", fill = "Genotype")
# snvs_interval_plt

# Creating the final, labeled plot
top <- cowplot::plot_grid(
  label_plot_2, genes_interval,
  align = "v",
  nrow = 2)

final_plt <- cowplot::plot_grid(
  top, snvs_interval_plt,
  align = "v",
  nrow = 2,
  rel_heights = c(0.3, 1)) + theme(plot.background = element_rect(fill = "white", color = NA))

# ggsave("/vast/eande106/projects/Lance/THESIS_WORK/misc/TESTER.png", final_plt, width = 7.5, height = 7.5, dpi = 600)


# Which gene has the most ALT-specific SNVs?
snvs_per_gene <- n2_genes_ROI %>%
  dplyr::mutate(gene_id = dplyr::row_number()) %>%
  dplyr::inner_join(
    alt_specific,
    by = "seqid",
    relationship = "many-to-many"
  ) %>%
  dplyr::filter(
    POS >= start,
    POS <= end
  ) %>%
  dplyr::count(gene_id, seqid, start, end, name = "n_snvs") %>%
  dplyr::arrange(desc(n_snvs)) %>% dplyr::left_join(N2_tran, by = c('seqid','start','end')) %>%
  dplyr::distinct(seqid, start, end, n_snvs, alias)
# WHICH OF THESE 115 VARIANTS HAVE A PREDICTED DELETERIOUS EFFECT???

















