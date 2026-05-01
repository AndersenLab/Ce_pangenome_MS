library(ggplot2)
library(ape)
library(cluster)
library(dplyr)

set.seed(13)

tree <- ape::read.tree("../../processed_data/genome_resources/trees/Ce_isotypes_LD_0.9.phy.contree")

# Patristic distance between all strains - how evolutionarily far apart are the strains? pairwise
D <- ape::cophenetic.phylo(tree)

# 140 wild strains
strains <- readr::read_tsv("../../processed_data/genome_resources/wild_strain_genome_stats.tsv") %>%
  dplyr::select(Strain) %>%
  dplyr::rename(strain = Strain) %>% 
  dplyr::pull()

tree$tip.label <- stringr::str_trim(tree$tip.label)

have_vec <- strains$strain
tips     <- tree$tip.label

in_both  <- intersect(tips, have_vec)
missing_in_tree <- setdiff(have_vec, tips)       # strains in your TSV but not in tree
extra_in_tree   <- setdiff(tips, have_vec)       # tips in tree you don't have

message("# strains you have: ", length(have_vec))
message("# tips in tree:     ", length(tips))
message("# overlap:           ", length(in_both))

# 4) subset the tree to only strains you actually have
subtree <- keep.tip(tree, in_both)

# 5) patristic distances on the subset only
D <- ape::cophenetic.phylo(subtree)
dist_obj <- as.dist(D)

# 6) run PAM on the subset; cap k if overlap < 64
k <- min(64, attr(dist_obj, "Size"))
pam_fit <- pam(dist_obj, k = k, diss = TRUE)

# 7) map medoid indices back to strain names
core_set <- attr(dist_obj, "Labels")[pam_fit$id.med]

present_core_set <- strains %>% dplyr::filter(strain %in% core_set)

# Final set of 64 strains 
write.table(present_core_set, "../../processed_data/orthology/orthofinder/core64_FINAL.tsv", quote = F, col.names = F, row.names = F)

D <- as.matrix(dist_obj)
core_idx <- match(core_set, attr(dist_obj, "Labels"))

nearest_medoid_dist <- apply(D, 1, function(row) min(row[core_idx]))
names(nearest_medoid_dist) <- attr(dist_obj, "Labels")

# histogram + ECDF
dfd <- data.frame(strain = names(nearest_medoid_dist),
                  d_to_core = as.numeric(nearest_medoid_dist),
                  is_core = names(nearest_medoid_dist) %in% core_set)

p_hist <- ggplot(dfd, aes(d_to_core)) +
  geom_histogram(bins = 40) +
  labs(x = "Patristic distance to nearest core strain",
       y = "Count",
       title = "Coverage of diversity by the 64-core set")

p_hist

p_ecdf <- ggplot(dfd, aes(d_to_core)) +
  stat_ecdf() +
  labs(x = "Patristic distance to nearest core strain",
       y = "ECDF",
       title = "Cumulative coverage by distance") +
  annotate("text", x = quantile(dfd$d_to_core, 0.9),
           y = 0.1,
           label = sprintf("Median=%.4f; 90th=%.4f",
                           median(dfd$d_to_core),
                           quantile(dfd$d_to_core,0.9)))

p_ecdf
