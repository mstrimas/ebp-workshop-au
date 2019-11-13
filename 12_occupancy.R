## ----occupancy-data------------------------------------------------------
library(auk)
library(lubridate)
library(sf)
library(dggridR)
library(unmarked)
library(raster)
library(ebirdst)
library(MuMIn)
library(AICcmodavg)
library(fields)
library(tidyverse)
# resolve namespace conflicts
select <- dplyr::select
projection <- raster::projection

set.seed(1)

# ebird data
ebird <- read_csv("data/ebd_woothr_june_bcr27_zf.csv") %>% 
  mutate(year = year(observation_date),
         # occupancy modeling requires an integer response
         species_observed = as.integer(species_observed))

# modis land cover covariates
habitat <- read_csv("data/pland-elev_location-year.csv") %>% 
  mutate(year = as.integer(year))

# combine ebird and modis data
ebird_habitat <- inner_join(ebird, habitat, by = c("locality_id", "year"))

# prediction surface
pred_surface <- read_csv("data/pland-elev_prediction-surface.csv")
# latest year of landcover data
max_lc_year <- pred_surface$year[1]
r <- raster("data/prediction-surface.tif")

# load gis data for making maps
map_proj <- st_crs(102003)
ne_land <- read_sf("data/gis-data.gpkg", "ne_land") %>% 
  st_transform(crs = map_proj) %>% 
  st_geometry()
bcr <- read_sf("data/gis-data.gpkg", "bcr") %>% 
  st_transform(crs = map_proj) %>% 
  st_geometry()
ne_country_lines <- read_sf("data/gis-data.gpkg", "ne_country_lines") %>% 
  st_transform(crs = map_proj) %>% 
  st_geometry()
ne_state_lines <- read_sf("data/gis-data.gpkg", "ne_state_lines") %>% 
  st_transform(crs = map_proj) %>% 
  st_geometry()


## ----occupancy-prep-filter-----------------------------------------------
# filter to a single year of data
ebird_filtered <- filter(ebird_habitat, 
                         number_observers <= 5,
                         year == max(year))


## ----occupancy-prep-repeats, class.source="livecode"---------------------
# subset for occupancy modeling

### LIVE CODE ###


## ----occupancy-prep-unmarked, class.source="livecode"--------------------
# format for unmarked, select occupancy and detection covariates

### LIVE CODE ###


## ----encounter-prep-sss, results="hide"----------------------------------
# generate hexagonal grid with ~ 5 km betweeen cells
dggs <- dgconstruct(spacing = 5)
# get hexagonal cell id for each site
occ_wide_cell <- occ_wide %>% 
  mutate(cell = dgGEO_to_SEQNUM(dggs, longitude, latitude)$seqnum)
# sample one checklist per grid cell
occ_ss <- occ_wide_cell %>% 
  group_by(cell) %>% 
  sample_n(size = 1) %>% 
  ungroup() %>% 
  select(-cell)


## ----encounter-data-unmarked, class.source="livecode"--------------------
# creat unmarked object

### LIVE CODE ###


## ----occupancy-model-fit, class.source="livecode"------------------------
# fit model

### LIVE CODE ###


## ----occupancy-model-assess, eval=FALSE, echo=1:2, class.source="dontrun"----
occ_gof <- mb.gof.test(occ_model, nsim = 10, plot.hist = FALSE)
print(occ_gof)


## ----occupancy-predict-predict, eval=FALSE, echo=1:10, class.source="livecode"----

### LIVE CODE ###


## ----occupancy-predict-rasterize-----------------------------------------
r_pred <- pred_occ %>% 
  # convert to spatial features
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
  st_transform(crs = projection(r)) %>% 
  # rasterize
  rasterize(r)
r_pred <- r_pred[[c("occ_prob", "occ_se")]]


## ----occupancy-predict-map, fig.asp=1.236--------------------------------
# project predictions
r_pred_proj <- projectRaster(r_pred, crs = map_proj$proj4string, method = "ngb")

par(mfrow = c(2, 1))
for (nm in names(r_pred)) {
  r_plot <- r_pred_proj[[nm]]
  
  par(mar = c(3.5, 0.25, 0.25, 0.25))
  # set up plot area
  plot(bcr, col = NA, border = NA)
  plot(ne_land, col = "#dddddd", border = "#888888", lwd = 0.5, add = TRUE)

  # occupancy probability or standard error
  if (nm == "occ_prob") {
    title <- "Wood Thrush Occupancy Probability"
    brks <- seq(0, 1, length.out = 21)
    lbl_brks <- seq(0, 1, length.out = 11) %>% 
      round(2)
  } else {
    title <- "Wood Thrush Occupancy Uncertainty (SE)"
    mx <- ceiling(1000 * cellStats(r_plot, max)) / 1000
    brks <- seq(0, mx, length.out = 21)
    lbl_brks <- seq(0, mx, length.out = 11) %>% 
      round(2)
  }
  pal <- abundance_palette(length(brks) - 1)
  plot(r_plot, 
       col = pal, breaks = brks, 
       maxpixels = ncell(r_plot),
       legend = FALSE, add = TRUE)
  
  # borders
  plot(bcr, border = "#000000", col = NA, lwd = 1, add = TRUE)
  plot(ne_state_lines, col = "#ffffff", lwd = 0.75, add = TRUE)
  plot(ne_country_lines, col = "#ffffff", lwd = 1.5, add = TRUE)
  box()
  
  # legend
  par(new = TRUE, mar = c(0, 0, 0, 0))
  image.plot(zlim = range(brks), legend.only = TRUE, 
             breaks = brks, col = pal,
             smallplot = c(0.25, 0.75, 0.06, 0.09),
             horizontal = TRUE,
             axis.args = list(at = lbl_brks, labels = lbl_brks,
                              fg = "black", col.axis = "black",
                              cex.axis = 0.75, lwd.ticks = 0.5,
                              padj = -1.5),
             legend.args = list(text = title,
                                side = 3, col = "black",
                                cex = 1, line = 0))
}


## ----occupancy-model-select-dredge, class.source="livecode"--------------

### LIVE CODE ###


## ----occupancy-model-select-average, class.source="livecode"-------------

### LIVE CODE ###


## ----occupancy-model-select-predict, eval=FALSE, class.source="dontrun"----
## #
## occ_pred_avg <- predict(occ_avg,
##                         newdata = as.data.frame(pred_surface),
##                         type = "state")
## 
## # add to prediction surface
## pred_occ_avg <- bind_cols(pred_surface,
##                           occ_prob = occ_pred_avg$fit,
##                           occ_se = occ_pred_avg$se.fit) %>%
##   select(latitude, longitude, occ_prob, occ_se)


## ----occupancy-select-detection-define-----------------------------------
# define occupancy model formula with only effort in detection submodel
det_mod <- ~ time_observations_started + 
  duration_minutes + 
  effort_distance_km + 
  number_observers + 
  protocol_type ~ 
  pland_04 + pland_05 + pland_12 + pland_13

# create new formulae with landcover covariates added to the detection submodel
mods <- list(det_mod_null = det_mod, 
             det_mod_dec = update.formula(det_mod, ~ . + pland_04 ~ .),
             det_mod_mix = update.formula(det_mod, ~ . + pland_05 ~ .),
             global = update.formula(det_mod, 
                                     ~ . + pland_04 + pland_05 ~ .)) %>% 
  # fit candidate models
  map(occu, data = occ_um)


## ----occupancy-select-detection-compare----------------------------------
mod_sel <- fitList(fits = mods) %>% 
  modSel()
mod_sel


## ----occupancy-select-detection-coef-------------------------------------
coef(occ_model) %>% 
  enframe() %>% 
  filter(str_detect(name, "pland_0"))


# Exercises ----

# 1. Try sampling more than a single checklist per grid cell in the
# spatiotemporal sampling. How does that affect model fit and predictions?

# 2. What happens to the size of dataset if you only use stationary counts, or
# reduce the distance traveled to 1 km? How does it impact the results? How does
# the different input data affect your interpretation of the results?

# 3. What happens to the size of the dataset if you allow repeat visits to be by
# multiple observers? How does this impact the results.

# 4. Produce a map based on model averaged predictions. Note that making these
# predictions may take up to an hour.