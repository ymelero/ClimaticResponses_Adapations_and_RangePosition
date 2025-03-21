# ---------------------------------------------------------------------------------------------- #
# Script to model abundance over time in relation to climatic anomalies & range position
# Author: [YM, LE]
# Inputs:
#   - Data_analyses.csv available in Zenodo. 
#   - Year data is available via signed license agreement (https://butterfly-monitoring.net/)
# Outputs:
#   - model results, predictions and plots
## ---------------------------------------------------------------------------------------------- #


# 0. Add time series index (year data from eBMS under license)
Data_analyses <- Data_analyses %>%
  group_by(species) %>%
  mutate(Timeseries = year - min(year) + 1)

# 1. Inverse hyperbolic sine transformation
IHST <- function(x) log(x + sqrt(x^2 + 1))

# 2. Fit abundance trends over time model
Abd.model <- lmer(IHST(gam_index) ~ Timeseries * rangeposition + 
                    (Timeseries | species) + (1 | site), data = mydf)

# 3. Model predictions
mydf.predicted2 <- ggpredict(Abd.model, terms = c("Timeseries[all]", 
                                                 paste0("rangeposition[", paste(vals, collapse = ","), "]")))
# Back-transform predicted values (to see abundance's values directly)
mydf.predicted2 <- mydf.predicted2 %>%
  mutate(predicted.exp = (exp(2 * predicted) - 1) / (2 * exp(predicted)))

# 4. Plot predictions
plot2<- ggplot(mydf.predicted2, aes(x, predicted.exp, colour = group)) + 
  geom_line(linewidth = 1) +
  scale_color_manual(values = pal) +
  coord_cartesian(xlim = c(0, 18), ylim = c(1, 5.5)) +
  annotate("text", x = 0, y = 5.5, label = "(a) Locally adapted species", hjust = 0, vjust = 0) +
  labs(x = "Number of years", y = "Population abundance") +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black")
  )