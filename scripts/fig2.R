#!/usr/bin/env Rscript

# Load required packages and functions
library(tidyverse) 
library(ggtree)
library(tidytree)
library(ggplot2)
library(labdsv)
library(data.table)
source("scripts/r_functions.R")

`%nin%` = Negate(`%in%`)

##### ANI MATRIX ####

# Load and root the dated tree (wASTRAL topology)
dated_tree <- read.tree("analyses/phylogenetics/set103/divtime/1_part/mcmc/c1/dated.tree")
dated_tree <- dated_tree %>%
  treeio::root(outgroup = c("Aphanizomenon_flos_aquae_NIES_81.fa",
                            "Anabaena_cylindrica_PCC_7122.fa",
                            "Cylindrospermum_stagnale_PCC_7417.fa"), 
               resolve.root = T)
outgroup = c("Aphanizomenon_flos_aquae_NIES_81.fa",
             "Anabaena_cylindrica_PCC_7122.fa",
             "Cylindrospermum_stagnale_PCC_7417.fa")
# Load FastANI output
fastani_df <- read_delim(file = "analyses/species_delimitation/fastani/set12c/fastani_set12c_ql_out", 
                         col_names = FALSE) %>%
  mutate(X1 = 
           str_remove(X1, "analyses/cyano_genomes/set12c/")) %>%
  mutate(X1 = 
           str_remove(X1, "_chromosome")) %>%
  mutate(X2 = 
           str_remove(X2, "analyses/cyano_genomes/set12c/")) %>%
  mutate(X2 = 
           str_remove(X2, "_chromosome")) %>%
  mutate(alignment_fraction = 
           X4/X5) %>%
  rowwise() %>%
  mutate(genome1 = 
           min(c(X1, X2))) %>%
  mutate(genome2 = 
           max(c(X1, X2))) %>%
  ungroup() %>%
  distinct(genome1, genome2, .keep_all = T) %>%
  filter(genome1 %nin% outgroup & genome2 %nin% outgroup) %>%
  select(!c(X1, X2, X4, X5)) %>%
  rename(ani = X3) %>%
  relocate(genome1, genome2)
#
gap <- function(fastani_df, low_lim, hi_lim) {
  genomes <- fastani_df %>% dplyr::pull(1) %>%
    base::unique()
  gap <- numeric()
  gap_low_lim <- numeric()
  gap_hi_lim <- numeric()
  for (i in 1:length(genomes)) {
    genome <- genomes[i]
    anis <- fastani_df %>% dplyr::filter(genome1 == genome | 
                                              genome2 == genome) %>% 
      pull(ani) %>% 
      base::subset(. > low_lim & . <= hi_lim) %>%
      sort()
    diffs <- diff(anis)
    gap[i] <- max(diffs)
    if (is.finite(gap[i])) {
      gap_low_lim_index <- which(diffs == gap[i]) %>%
        unique()
      gap_low_lim[i] <- anis[gap_low_lim_index]
      gap_hi_lim[i] <- anis[gap_low_lim_index+1]
    } else {
      gap_low_lim[i] <- NA
      gap_hi_lim[i] <- NA
    }
  }
  out <- dplyr::bind_cols(genomes, gap, gap_low_lim, gap_hi_lim) %>%
    mutate(`...2` =
             ifelse(gap == -Inf, NA, gap))
  colnames(out) <- c("genomes", "gap",  "gap_low_lim", "gap_hi_lim")
  out
}

gap_df <- gap(fastani_df, low_lim = 88, hi_lim = 100)

anis <- fastani_df %>% filter(genome1 == "P8231_bin_15.fa" | 
                                genome2 == "P8231_bin_15.fa") %>% 
  pull(ani) %>% 
  base::subset(. > 88) %>%
  sort() %>%
  diff() %>%
  max()



# Converta ANI data to matrix
ani_matrix <- fastani_df %>%
  select(1:3) %>% 
  as.matrix() %>%
  matrify() %>%
  mutate_if(is.character, as.numeric)
ani_matrix[ani_matrix == 0] <- NA


#
ani_hist <- ggplot(data = fastani_df, aes(x = ani)) +
  geom_histogram(fill = "gray40") +
  labs(x = "ANI", y = "frequency") +
  scale_x_continuous(n.breaks = 13) +
  theme(panel.background = NULL, 
        panel.border = element_rect(fill = "transparent", linewidth = 0.75),
        axis.text = element_text(size = 12))
#
ani_vs_aln <- ggplot(data = fastani_df, aes(x = ani, y = alignment_fraction)) +
  geom_point(size = 0.5, alpha = 0.5) +
  labs(x = "ANI", y = "Alignment fraction") +
  scale_x_continuous(n.breaks = 13) +
  theme(panel.background = NULL, 
        panel.border = element_rect(fill = "transparent", linewidth = 0.75),
        axis.text = element_text(size = 12))






ani_matrix <- ani_matrix %>%
  select(all_of(ordered_tips)) %>%
  as.data.frame()
rownames(ani_matrix) <- ordered_tips
# Plot tree with matrix
astral_tree_plot <- ggtree(astral_tree)
gheatmap(astral_tree_plot, ani_matrix, colnames=F) +
  scale_fill_continuous(type = "viridis")
ggsave(h, filename = "test.pdf", device = "pdf")
#

tree <- read.tree(text = "(((A,B),(C,D)),E);")
tree2 <- ladderize(tree, right = FALSE)
tree$tip.label
#> [1] "A" "B" "C" "D" "E"
tree2$tip.label
#> [1] "A" "B" "C" "D" "E"
plot(tree2)
nodelabels()
tiplabels()

astral_tree1 <- ladderize(astral_tree, right = F)
is_tip <- astral_tree1$edge[,2] <= length(astral_tree1$tip.label)
ordered_tips <- astral_tree1$edge[is_tip, 2]
ordered_tips <- astral_tree1$tip.label[ordered_tips] %>% rev()

