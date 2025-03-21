# ---------------------------------------------------------------------------------------------- #
# Script to calculate bioclimatic range position per species population (site)
# Author: [LE, YM]
# Inputs: 
#   - Data from eBMS license agreement & ECAD website (https://www.ecad.eu/).
# Outputs:
#   - Bioclimatic Range position. Available at Zenodo.
## ---------------------------------------------------------------------------------------------- #

# Function to generate a standardized range position
standardrange <- function(x) {
  x <- na.omit(x)  # Remove missing values
  x1 <- x - min(x)  # Shift values so the minimum becomes zero
  mm <- max(x) - min(x)  # Calculate the total range of the variable
  
  # Scale the values between -1 and 1
  ans <- 2 * x1 / mm - 1
  return(ans)
}

# Compute the mean climatic value at each site and standardize range position for each species
data <- data.ebms %>%
  group_by(site) %>%
  distinct(year, .keep_all = TRUE) %>%  # Ensure only unique year entries per site
  summarise(sitemeanWeatherPeriod = mean(clima.site, na.rm = TRUE)) %>%  # Compute mean climate per site
  right_join(data.ebms, by = "site") %>%  # Merge with original dataset
  group_by(species) %>%  # Process data per species
  mutate(rangeposition = standardrange(sitemeanWeatherPeriod)) %>%
  ungroup()  # Ensure calculations apply independently per species

# clima.site represents the climatic variable of interest (Temperature, Precipitation, or Aridity) 
# for each species at a given site and date, extracted from ECAD.
# Note: Precipitation should be inverted (negative values) to align -1 and 1 with trailing and leading margins.
