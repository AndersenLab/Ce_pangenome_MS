library(tidyverse)
library(data.table)
library(patchwork)
library(ggbeeswarm)
library(ggrepel)
library(glue)
library(circlize)
library(ggnewscale)
library(cowplot)

# ============================================================
# USER SETTINGS + DATA
# ============================================================
mapping_file <- "../../tables/mapping_2018GWAS_MLs.tsv"
genotype_matrix_file <- "../../tables/gwas/genotype_matrix.tsv"
mapping_id_col <- "strain"
independent_tests <- 9462.00000000002
LOCO_gwa_file <- "../../tables/abamectin_q90.TOF_ctrl-regressed_loco.gwa"
hdrs <- readr::read_tsv("../../data/genomic_resources/20250625_c_elegans_divergent_regions_strain.bed", col_names = c("chrom", "start", "end", "strain"))
glc1 <- data.frame(CHR = "V", xmin = (16219634 / 1e6) - 0.02, xmax = (16221917 / 1e6) + 0.02)
strains_phenotyped <- readr::read_tsv("../../tables/mapping_2018GWAS_MLs.tsv") %>% dplyr::pull(strain)

# ============================================================
# HELPERS
# ============================================================
safe_read_tsv <- function(path) {
  if (is.null(path) || !file.exists(path)) return(NULL)
  x <- tryCatch(
    fread(path, data.table = FALSE, sep = "\t", header = TRUE, fill = TRUE),
    error = function(e) NULL
  )
  if (is.null(x) || nrow(x) == 0) return(NULL)
  names(x) <- trimws(names(x))
  x
}

safe_read_csv <- function(path) {
  if (is.null(path) || !file.exists(path)) return(NULL)
  x <- tryCatch(
    fread(path, data.table = FALSE, sep = ",", header = TRUE, fill = TRUE),
    error = function(e) NULL
  )
  if (is.null(x) || nrow(x) == 0) return(NULL)
  names(x) <- trimws(names(x))
  x
}

read_genotype_matrix <- function(path) {
  x <- safe_read_tsv(path)
  if (is.null(x)) return(NULL)
  names(x) <- trimws(names(x))
  if (!all(c("CHROM", "POS", "REF", "ALT") %in% names(x))) return(NULL)
  x %>%
    mutate(
      CHROM = as.character(CHROM),
      POS = suppressWarnings(as.numeric(POS))
    )
}

standardize_gwa <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(NULL)
  names(df) <- trimws(names(df))

  # harmonize chromosome column
  if (!("CHR" %in% names(df))) {
    if ("CHROM" %in% names(df)) df <- dplyr::rename(df, CHR = CHROM)
    if ("Chr"   %in% names(df)) df <- dplyr::rename(df, CHR = Chr)
    if ("chr"   %in% names(df)) df <- dplyr::rename(df, CHR = chr)
  }

  # harmonize position column
  if (!("POS" %in% names(df))) {
    if ("BP" %in% names(df)) df <- dplyr::rename(df, POS = BP)
    if ("bp" %in% names(df)) df <- dplyr::rename(df, POS = bp)
    if ("pos" %in% names(df)) df <- dplyr::rename(df, POS = pos)
  }

  # harmonize p-value column
  if (!("P" %in% names(df))) {
    if ("p" %in% names(df)) df <- dplyr::rename(df, P = p)
    if ("PVAL" %in% names(df)) df <- dplyr::rename(df, P = PVAL)
    if ("pval" %in% names(df)) df <- dplyr::rename(df, P = pval)
  }

  if (!all(c("CHR", "POS", "P") %in% names(df))) return(NULL)

  df %>%
    mutate(
      CHR = as.character(CHR),
      POS = suppressWarnings(as.numeric(POS)),
      P   = suppressWarnings(as.numeric(P))
    ) %>%
    filter(!is.na(CHR), !is.na(POS), !is.na(P), P > 0)
}

find_trait_dirs <- function(root_dir) {
  dirs <- list.dirs(root_dir, recursive = FALSE, full.names = TRUE)
  dirs <- dirs[!basename(dirs) %in% c("Reports", "summary_figures")]
  dirs
}

find_gwa_file <- function(trait_dir, method) {
  exact <- file.path(trait_dir, "gwa", method, paste0(basename(trait_dir), "_", method, ".gwa"))
  if (file.exists(exact)) return(exact)

  alt1 <- Sys.glob(file.path(trait_dir, "gwa", method, paste0("*_", method, ".gwa")))
  if (length(alt1) > 0) return(alt1[[1]])

  alt2 <- Sys.glob(file.path(trait_dir, "gwa", method, "*.gwa"))
  if (length(alt2) > 0) return(alt2[[1]])

  NULL
}

find_qtl_file <- function(trait_dir, method) {
  candidates <- c(
    file.path(trait_dir, "gwa", method, paste0(basename(trait_dir), "_", method, "_qtl.tsv")),
    file.path(trait_dir, "gwa", method, paste0("*_", method, "_qtl.tsv")),
    file.path(trait_dir, "gwa", method, "*.qtl.tsv")
  )
  for (pat in candidates) {
    files <- Sys.glob(pat)
    if (length(files) > 0) return(files[[1]])
  }
  return(NULL)
}

read_qtl_and_sort <- function(qtl_file, method) {
  qtl <- safe_read_tsv(qtl_file)
  if (is.null(qtl) || nrow(qtl) == 0) return(NULL)

  qtl <- as_tibble(qtl)
  # keep only rows with marker / peakPOS
  qtl <- qtl %>%
    mutate(
      method = method,
      CHROM = as.character(CHROM),
      peakPOS = as.numeric(peakPOS),
      startPOS = as.numeric(startPOS),
      endPOS = as.numeric(endPOS)
    ) %>%
    arrange(
      factor(CHROM, levels = roman_chr_levels),
      peakPOS
    )
  qtl
}

find_finemap_folder <- function(trait_dir, method, qtl_row) {
  if (!all(c("CHROM", "startPOS", "endPOS") %in% names(qtl_row))) return(NULL)
  region_name <- paste0(qtl_row$CHROM, "_", qtl_row$startPOS, "_", qtl_row$endPOS)
  folder <- file.path(trait_dir, "fine_mapping", method, region_name)
  if (!dir.exists(folder)) return(NULL)
  folder
}

find_annotated_gwa <- function(finemap_folder) {
  if (is.null(finemap_folder) || !dir.exists(finemap_folder)) return(NULL)
  files <- Sys.glob(file.path(finemap_folder, "*.annotated.gwa"))
  if (length(files) == 0) return(NULL)
  files[[1]]
}

roman_chr_levels <- c("I", "II", "III", "IV", "V", "X", "MtDNA")


# ============================================================
# Collapse HDRs
# ============================================================
all_regions <- hdrs %>%
  dplyr::filter(strain %in% strains_phenotyped) %>%
  dplyr::rename(START = start, CHROM = chrom, END = end) %>%
  dplyr::arrange(CHROM,START) %>%
  dplyr::group_split(CHROM)

# Collapsing all HDRs
getRegFreq <- function(all_regions) {
  all_collapsed <- list()
  for (i in 1:length(all_regions)) {
    temp <- all_regions[[i]]
    k=1
    j=1
    while (k==1) {
      checkIntersect <- temp %>%
        dplyr::arrange(CHROM,START) %>%
        dplyr::mutate(check=ifelse(lead(START) <= END,T,F)) %>%
        dplyr::mutate(check=ifelse(is.na(check),F,check))
      
      print(nrow(checkIntersect %>% dplyr::filter(check==T)))
      
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
          dplyr::filter(check == FALSE & dplyr::coalesce(dplyr::lag(check), FALSE) == FALSE)
        # dplyr::filter(check==F & lag(check)==F)
        
        temp <- rbind(collapse,retain) %>%
          dplyr::select(-gid,-check)
        
        j=j+1
      }
    }
    print("WRITING TO LIST")
    print(head(temp))
    all_collapsed[[i]] <- temp
  }
  return(all_collapsed)
}

HDR_collapse_master <- getRegFreq(all_regions)

all_collapsed <- plyr::ldply(HDR_collapse_master, data.frame) %>%
  dplyr::select(-strain)

colnames(all_collapsed) <- c("CHR","start","end")


# ============================================================
# MANHATTAN
# ============================================================
make_manhattan_df <- function(gwa, bf_threshold) {
  gwa <- standardize_gwa(gwa)
  if (is.null(gwa) || nrow(gwa) == 0) return(NULL)

  gwa <- gwa %>%
    mutate(
      neglogp = -log10(P),
      sig = case_when(
        neglogp >= bf_threshold ~ "BF",
        TRUE ~ "NONSIG"
      ),
      CHR = factor(CHR, levels = c("I", "II", "III", "IV", "V", "X"))  # (1) Removed MtDNA
    ) %>%
    filter(!is.na(CHR))

  gwa
}

plot_manhattan <- function(gwa, bf_threshold) {
  df <- make_manhattan_df(gwa, bf_threshold)
  if (is.null(df) || nrow(df) == 0) {
    return(ggplot() + theme_void())
  }

  sig.colors <- c("BF" = "red", "NONSIG" = "black")
  sig.alpha  <- c("BF" = 1, "NONSIG" = 0.25)

  # per-chromosome horizontal threshold lines
  line_df <- df %>%
    dplyr::distinct(CHR) %>%
    dplyr::mutate(bf = bf_threshold)

  ggplot(df) +
    geom_rect(data = all_collapsed, aes(xmin = start / 1e6, xmax = end / 1e6, ymin = 0, ymax = 8), fill = 'black', alpha = 0.2) +
    geom_rect(data = glc1, aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = 8), fill = 'purple') +
    geom_point(aes(x = POS / 1e6, y = neglogp, colour = sig, alpha = sig), size = 0.6) +
    scale_colour_manual(values = sig.colors) +
    scale_alpha_manual(values = sig.alpha) +
    geom_hline(data = line_df, aes(yintercept = bf), linetype = 1) +
    facet_grid(. ~ CHR, scales = "free_x", space = "free_x", drop = FALSE) +
    labs(x = "N2 genomic position (Mb)", y = expression(-log[10](italic(p)))) +
    scale_y_continuous(expand = c(0,0)) +
    # scale_x_continuous(expand = c(0,0)) +
    theme_bw(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_text(face = "plain"))
  }

# ============================================================
# LOAD DATA AND PLOT
# ============================================================
bf_threshold <- -log10(0.05 / independent_tests)

# LOCO GWA
gwa <- safe_read_tsv(LOCO_gwa_file)
head(gwa)

# manhattan plot
p <- plot_manhattan(gwa, bf_threshold)
print(p)

# Save the plot
ggsave("../../figures/supplementary/abamectin_q90.TOF_ctrl-regressed_manhattan.png", p, width = 7.5, height = 4, dpi = 600)
