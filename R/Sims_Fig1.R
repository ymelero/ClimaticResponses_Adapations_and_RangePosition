# ---------------------------------------------------------------------------------------------- #
# Script to plot sites and map and simulated hypthesis 1
# Author: [LE, YM]
# Inputs:
#   - Climate data is available via ECAD website (https://www.ecad.eu/).
#   - Location data is available via signed license agreement (https://butterfly-monitoring.net/)
# Outputs:
#   - model results, predictions and plots
## ---------------------------------------------------------------------------------------------- #

library(raster)
library(ncdf4)
library(tidyverse)
library(sf)
library(rnaturalearth)
library(scico)  

# 1. Temperature maps

# Load the NetCDF file downloaded from ECAD 
nc_file <- 'tg_0.25deg_reg_v17.0.nc'
euro <- nc_open(nc_file)

# Extract longitude, latitude, and time variables
lon <- ncvar_get(euro, "longitude")
lat <- ncvar_get(euro, "latitude")
time <- ncvar_get(euro, "time")

# Get time units and convert time to actual dates
time_units <- ncatt_get(euro, "time", "units")$value
origin_date <- as.Date(substr(time_units, 12, 21))
dates <- origin_date + time

# Define the period of interest (1995-2017)
study_start <- as.Date("1995-01-01")
study_end <- as.Date("2017-12-31")
start_index <- which.min(abs(dates - study_start))
end_index <- which.min(abs(dates - study_end))

# Extract and calculate mean temperature for the study period
ndims <- euro$ndims
vsize <- euro$var$tg$varsize

# Initialize raster with the first time slice
start <- rep(1, ndims)
start[ndims] <- start_index
count <- vsize
count[ndims] <- 1
baser <- raster(ncvar_get(euro, "tg", start = start, count = count))

# Loop through the time slices and calculate the mean
for (i in (start_index + 1):end_index) {
  start[ndims] <- i
  temp_raster <- raster(ncvar_get(euro, "tg", start = start, count = count))
  baser <- baser + temp_raster
}

# Calculate the average temperature
baser <- baser / (end_index - start_index + 1)

# Flip and correct the extent for visualization
baser <- t(flip(baser, direction = 1))
extent(baser) <- c(min(lon), max(lon), min(lat), max(lat))

# Close the NetCDF file
nc_close(euro)

# Crop the raster to Europe with adjusted extent (removing Iceland)
europe_extent <- extent(-25, 31, 35, 71)  # Adjusted northern limit to remove Iceland
baser_cropped <- crop(baser, europe_extent)

# Convert raster to data frame for ggplot
baser_df <- as.data.frame(rasterToPoints(baser_cropped), xy = TRUE)
colnames(baser_df) <- c("lon", "lat", "temp")

# Load Europe shapefile for overlay
europe <- ne_countries(scale = "large", continent = "Europe", returnclass = "sf")

# Exclude Iceland from the shapefile
europe <- europe %>% filter(!grepl("Iceland", name_long))

map_plot<-ggplot() +
  geom_raster(data = baser_df, aes(x = lon, y = lat, fill = temp)) +
  geom_sf(data = europe, fill = NA, color = "black", size = 0.2) +
  scale_fill_scico(palette = "roma", direction = -1) +  # Invertimos los colores con direction = -1
  coord_sf(xlim = c(-25, 31), ylim = c(35, 71), expand = FALSE) +
  theme_void() +
  theme(
    legend.position = "none"  # Quitamos la leyenda
  )

# 2. Add random selection of data points:
library(sp)
library(ggsci)
library(RColorBrewer)

# Data avalible under lisence agreement via the eBMS
catcord <- readRDS("cat/section_coordR.rds")
fincord <- readRDS("fin/section_coordR.rds")
ukcord <- readRDS("uk/section_coordR.rds")

catcord <- catcord %>% distinct(transect_id, .keep_all = TRUE) %>% filter(!is.na(section_lat))
fincord <- fincord %>% distinct(transect_id, .keep_all = TRUE) %>% filter(!is.na(section_lat))
ukcord <- ukcord %>% distinct(transect_id, .keep_all = TRUE) %>% filter(!is.na(section_lat))

# to WGS84 (EPSG:4326)
catcord <- st_as_sf(catcord, coords = c("section_lon", "section_lat"), crs = 3035) %>%
  st_transform(crs = 4326)
fincord <- st_as_sf(fincord, coords = c("section_lon", "section_lat"), crs = 3035) %>%
  st_transform(crs = 4326)
ukcord <- st_as_sf(ukcord, coords = c("section_lon", "section_lat"), crs = 3035) %>%
  st_transform(crs = 4326)

# Select
set.seed(1234)  # Aseguramos reproducibilidad
scatcord <- catcord %>% slice_sample(n = 10)
sfincord <- fincord %>% slice_sample(n = 10)
sukcord <- ukcord %>% slice_sample(n = 20)

# sf to data.frames for ggplot
scatcord_df <- as.data.frame(st_coordinates(scatcord))
sfincord_df <- as.data.frame(st_coordinates(sfincord))
sukcord_df <- as.data.frame(st_coordinates(sukcord))

map_plot + geom_point(data = scatcord_df, aes(X, Y), pch = 21, color = "black", fill = "yellow", alpha = 0.9, size = 1.7) +
  geom_point(data = sukcord_df, aes(X, Y), pch = 24, color = "black", fill = "magenta", alpha = 0.9, size = 1.5) +
  geom_point(data = sfincord_df, aes(X, Y), pch = 22, color = "black", fill = "black", alpha = 0.9, size = 1.7)

# 3. Sidebar
# Convert to sf objects using the raster CRS
library(raster)
scatcord_sp <- as_Spatial(scatcord)
sfincord_sp <- as_Spatial(sfincord)
sukcord_sp <- as_Spatial(sukcord)

# Extract temperature values for all combined coordinates
temp_cat <- na.omit(data.frame(temp = raster::extract(baser_cropped, scatcord_sp), country = 'CAT'))
temp_fin <- na.omit(data.frame(temp = raster::extract(baser_cropped, sfincord_sp), country = 'FIN'))
temp_uk <- na.omit(data.frame(temp = raster::extract(baser_cropped, sukcord_sp), country = 'UK'))

temp_all_df <- rbind(temp_cat, temp_fin, temp_uk)
temp_all_range <- range(temp_all_df$temp)

# Standardize temperature values from -1 to 1
temp_all_df$temp_standardized <- (temp_all_df$temp - temp_all_range[1]) / (temp_all_range[2] - temp_all_range[1]) * 2 - 1

# Creta a gprah bar - with ranges of values 
library(cowplot)
ranges <- seq(1, -1, by = -0.2)  # Reverse order of gradient
lev <- rep(1, length(ranges))  # Ensure all squares are in the same row
colorbar <- data.frame(ranges, lev)

colorbar1 <- ggplot(colorbar, aes(lev, ranges)) +  # Swap x and y to make it vertical
  geom_point(aes(color = ranges), size = 20, shape = 15) +  
  scale_color_gradient2(midpoint = 0, low = "#1B58A9", high = "#A0522D") +  
  theme_cowplot() +  
  xlim(0.95, 1.01) +  # Reduce space between bar and labels
  ylim(1, -1) +  # Reverse y-axis
  theme(
    axis.title = element_blank(), 
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    axis.text = element_blank(),  
    axis.ticks = element_blank(),  
    legend.position = "none"  
  ) +  
  annotate("text", y = ranges, x = rep(1.004, length(ranges)),   
           label = as.character(ranges), angle = 0, size = 4, hjust = 0) +
  
  geom_point(data = temp_all_df[temp_all_df$country == 'FIN', ], aes(y = temp_standardized, x = 1), pch = 22, color = "black", fill = "black", alpha = 0.9, size = 4)+
  geom_point(data = temp_all_df[temp_all_df$country == 'UK', ], aes(y = temp_standardized, x = 1), pch = 24, color = "black", fill = "magenta", alpha = 0.9, size = 4) +
  geom_point(data = temp_all_df[temp_all_df$country == 'CAT', ], aes(y = temp_standardized, x = 1), pch = 21, color = "black", fill = "yellow", alpha = 0.9, size = 4) 

# 4. Simulated population responses (H1)
library(tidyverse)
library(scico)

# define quadratic function
globalfunc = function(x, a, p1, p2) a + x * p1 + x^2 * p2

# input segments
segmentssim = c(
  seq(-1,-0.9,length=10), seq(-0.9,-0.7,length=10), seq(-0.7,-0.5,length=10),
  seq(-0.5,-0.3,length=10), seq(-0.3,-0.1,length=10), seq(-0.09,0.09,length=10),
  seq(0.1,0.3,length=10), seq(0.3,0.5,length=10), seq(0.5,0.7,length=10),
  seq(0.7,0.9,length=10), seq(0.9,1,length=10)
)

segnum = rep(1:11, each = 10)

# global case
restab = cbind.data.frame(
  pos = segmentssim,
  func = globalfunc(segmentssim, 1, 0, -1),
  segnum = segnum
)

predictedmod = c()
coefs = c()
for(i in unique(restab$segnum)) {
  d = filter(restab, segnum == i)
  lmmod = lm(func ~ pos + I(pos^2), data = d)
  predictedmod = c(predictedmod, predict(lmmod))
  coefs = c(coefs, coef(lmmod)[[2]])
}
restab = cbind.data.frame(restab, predictedmod)

cols1 = rev(scico::scico(11, palette = "roma"))

rlong = restab %>%
  group_by(segnum) %>%
  mutate(
    maxh = max(predictedmod),
    medh = median(predictedmod),
    minh = min(predictedmod),
    y = predictedmod / maxh
  ) %>%
  ungroup()

qmod = lm(func ~ pos + I(pos^2), data = restab)
rlong = cbind.data.frame(rlong, pred2 = predict(qmod))
rlong = rlong %>%
  mutate(dy = 0.5 - medh) %>%
  mutate(y2 = predictedmod + dy)

rlong$pos2 = rep(rlong$pos[51:60], 11)

globalcross = filter(rlong, segnum %in% 2:10) %>%
  ggplot() +
  geom_line(aes(pos2, y2, group = segnum, color = as.factor(segnum)), lwd = 2.5, alpha = 0.8) +
  theme_void() +
  theme(legend.position = "none") +
  scale_color_manual(values = cols1) +
  ylim(0, 1)

# local case
ranges3 = seq(-1, -0.9, length = 10)
funclo = rep(
  globalfunc(ranges3 + 0.9, 0, -90, -900) / max(globalfunc(ranges3 + 0.9, 0, -90, -900)),
  11
)

restablo = cbind.data.frame(
  pos = segmentssim,
  func = funclo,
  segnum = segnum
)

predictedmodlo = c()
coefslo = c()
for(i in unique(restab$segnum)) {
  d = filter(restablo, segnum == i)
  lmmod = lm(func ~ pos + I(pos^2), data = d)
  predictedmodlo = c(predictedmodlo, predict(lmmod))
  coefslo = c(coefslo, coef(lmmod)[[2]])
}
restablo = cbind.data.frame(restablo, predictedmodlo)

rlong_lo = restablo %>%
  group_by(segnum) %>%
  mutate(
    maxh = max(predictedmodlo),
    medh = median(predictedmodlo),
    minh = min(predictedmodlo),
    y = (predictedmodlo - minh) * 0.2
  ) %>%
  ungroup()

qmodlo = lm(func ~ pos + I(pos^2), data = restablo)
rlong_lo$y2 = predict(qmodlo) * 0.2

rlong_lo = rlong_lo %>%
  mutate(
    dy = 0.5 - medh,
    y2 = (predictedmodlo - minh) * 0.2
  )

rlong_lo$pos2 = rep(rlong_lo$pos[51:60], 11)
rlong_lo = rlong_lo %>%
  group_by(segnum) %>%
  mutate(pos2 = pos2 + rnorm(1, 0, 0.003)) %>%
  ungroup()

localcross = filter(rlong_lo, segnum %in% 2:10) %>%
  ggplot() +
  geom_line(aes(pos2, y2, group = segnum, color = as.factor(segnum)), lwd = 1, alpha = 0.8) +
  theme_void() +
  theme(legend.position = "none") +
  scale_color_manual(values = cols1) +
  ylim(0, 0.4)
