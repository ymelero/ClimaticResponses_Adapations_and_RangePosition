# ---------------------------------------------------------------------------------------------- #
# Script to model population growth rate in relation to climatic anomalies & range position
# Author: [YM]
# Inputs:
#   - Data_analyses.csv available in Zenodo.
# Outputs:
#   - model results, predictions and plots
## ---------------------------------------------------------------------------------------------- #

# 1. Quadratic models to be performed separately for locally and globally adapted species
Data_analyses <- data.frame(Data_analyses %>%
                     group_by(site, species) %>%
                     mutate(
                       n_change = log((lead(Nt)/Nt))
                     ))

gr.model<-lmer(n_change ~logNt +rangeposition*local.anomaly+I(local.anomaly^2)+(1|species) +(1|site),control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)) ,REML = T ,data=butterflydatadens.global)

# 2. Model predictions
vals <- seq(1, -1, by = -0.1)
mydf.predicted <- ggpredict(gr.model, terms = c("local.anomaly[all]", paste0("rangeposition[", paste(vals, collapse = ","), "]")))

# 3. Plot predictions
# Define custom reversed color palette
pal <- rev(c("#FF0000FF", "#FF0000FF", "#FF2A00FF", "#FF5500FF", "#FF8000FF", "#FFAA00FF", "#FFD500FF",
             "#FFFF00FF", "#FFFF40FF", "#FFFFBFFF", "#FFFFBFFF", "#F7FBFF", "#DEEBF7", "#C6DBEF",
             "#9ECAE1", "#6BAED6", "#4292C6", "#2171B5", "#08519C", "#08306B", "#08306B"))

plot1 <- ggplot(mydf.gr, aes(x, predicted, colour = group)) + 
  geom_line(size = 1) +
  scale_color_manual(values = pal) +
  coord_cartesian(xlim = c(-7, 5), ylim = c(-1, 0.5)) +
  annotate("text", x = -3.7, y = 0.45, label = "(a) Locally adapted species", vjust = 0.5, hjust = 0) +
  labs(x = "Climatic anomaly of the year", y = NULL) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black")
  )
