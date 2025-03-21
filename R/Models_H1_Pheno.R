# ---------------------------------------------------------------------------------------------------------------------------- #
# Script to model population growth rate in relation to climatic anomalies & range position with a phylogenetic analyses
# Author: [YM]
# Inputs:
#   - Data_analyses.csv available in Zenodo.
#   - Butterflies_Europe_tree.nwk is avalible from Dapporto et al (2019) Doi 10.1111/1755-0998.13059, and Wiemers et al (2020) doi 10.3897/zookeys.938.50878.
# Outputs:
#   - model results, predictions and plots
## --------------------------------------------------------------------------------------------------------------------------- #

# A. Frequentist Option
library(ape)
library(phyr)
library(ggeffects)
library(ggplot2)

# 0. Read and clean tree
phylo_tree <- read.tree("Butterflies_Europe_tree.nwk")
phylo_tree$tip.label <- gsub("_", " ", phylo_tree$tip.label)

# match tree and data
butterfly_species <- unique(Data_Analyses$species)
matched_species <- intersect(butterfly_species, phylo_tree$tip.label)
pruned_tree <- drop.tip(phylo_tree, setdiff(phylo_tree$tip.label, matched_species))
phylo_tree$edge.length <- rep(1, nrow(phylo_tree$edge))  # set edge lengths if missing

# assign phylogeny column
Data_Analyses$phylo <- Data_Analyses$species

# 1. Quadratic PGLMM models to be performed separately for locally and globally adapted species
model <- pglmm(
  n_change ~ logNt + rangeposition * local.anomaly + I(local.anomaly^2) +
    (1 | species__) + (1 | site),
  data = Data_Analyses,
  cov_ranef = list(species = pruned_tree),
  family = "gaussian",
  REML = TRUE
)

# 2. Model predictions
mydf.pheno <- ggpredict(model, terms = c("local.anomaly[all]",
                                         paste0("rangeposition[", paste(vals, collapse = ","), "]")))
# 3. Plot predictions
plot3 <- ggplot(mydf.pheno, aes(x, predicted, colour = group)) +
  geom_line(size = 1) +
  scale_color_manual(values = pal) +
  coord_cartesian(xlim = c(-4, 6), ylim = c(-1, 0.5)) +
  annotate("text", x = -3.7, y = 0.45, label = "(b) locally adapted species", vjust = 0.5, hjust = 0) +
  labs(
    x = "Climatic anomaly of the year",
    y = "Pop change - model with Pheno (Frequentist)"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black")
  )

# B. Bayesian =ption: Results are equivalent
library(brms)
phylo_cov <- ape::vcv(pruned_tree)

# 1. Phylogenetic lmer Bayesian model
model_bayes <- brm(
  n_change ~ logNt + rangeposition * local.anomaly + I(local.anomaly^2) +
    (1 | site) + (1 | species) + (1 | gr(phylo, cov = phylo_cov)),
  data =Data_Analyses,
  family = gaussian(),
  data2 = list(phylo_cov = phylo_cov),
  cores = 4,
  iter = 2000,
  warmup = 1000,
  control = list(adapt_delta = 0.95)
)
