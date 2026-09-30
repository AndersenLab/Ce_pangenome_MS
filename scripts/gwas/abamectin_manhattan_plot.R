library(tidyverse)
library(data.table)
library(patchwork)
library(ggbeeswarm)
library(ggrepel)
library(glue)
library(circlize)
library(ggnewscale)

# ============================================================
# USER SETTINGS + DATA
# ============================================================
mapping_file <- "../../processed_data/gwas/mapping_2018GWAS_MLs.tsv"
genotype_matrix_file <- "../../processed_data/gwas/genotype_matrix.tsv"
mapping_id_col <- "strain"
independent_tests <- 9462.00000000002
LOCO_gwa_file <- "../../processed_data/gwas/abamectin_q90.TOF_ctrl-regressed_loco.gwa"

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

  ggplot(df, aes(x = POS / 1e6, y = neglogp, colour = sig, alpha = sig)) +
    geom_point(size = 0.6) +
    scale_colour_manual(values = sig.colors) +
    scale_alpha_manual(values = sig.alpha) +
    geom_hline(data = line_df, aes(yintercept = bf), linetype = 1) +
    facet_grid(. ~ CHR, scales = "free_x", space = "free_x", drop = FALSE) +
    labs(x = "N2 genomic position (Mb)", y = expression(-log[10](italic(p)))) +
    scale_y_continuous(expand = expansion(mult = c(0,0.01))) +
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
