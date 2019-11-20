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
ebird <- read_csv("data/ebd_zf_sep_tst.csv") %>% 
  filter(common_name == "Australian Ibis") %>% 
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
map_proj <- st_crs(3577)
ne_land <- read_sf("data/gis-data.gpkg", "ne_country") %>% 
  st_transform(crs = map_proj) %>% 
  st_geometry()
cmz <- read_sf("data/gis-data.gpkg", "cmz") %>% 
  filter(cmz_name == "Eastern Australia Temperate and Subtropical forests") %>% 
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

# look at the regression coefficients from the model

### LIVE CODE ###


## ----occupancy-model-assess, eval=FALSE, echo=1:2, class.source="dontrun"----
occ_gof <- mb.gof.test(occ_model, nsim = 10, plot.hist = FALSE)
print(occ_gof)


## ----occupancy-predict-predict, eval=FALSE, echo=1:10, class.source="livecode"----
# make prediction

# add to prediction surface

### LIVE CODE ###


## ----occupancy-predict-rasterize-----------------------------------------
r_pred <- pred_occ %>% 
  # convert to spatial features
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
  st_transform(crs = projection(r)) %>% 
  # rasterize
  rasterize(r)
r_pred <- r_pred[[c("occ_prob", "occ_se")]]


## ----occupancy-predict-map, fig.asp=1.6----------------------------------
# project predictions
r_pred_proj <- projectRaster(r_pred, crs = map_proj$proj4string, method = "ngb")

par(mfrow = c(2, 1))
for (nm in names(r_pred)) {
  r_plot <- r_pred_proj[[nm]]
  
  par(mar = c(0.25, 0.25, 0.25, 5))
  # set up plot area
  plot(cmz, col = NA, border = NA)
  plot(ne_land, col = "#dddddd", border = "#888888", lwd = 0.5, add = TRUE)

  # occupancy probability or standard error
  if (nm == "occ_prob") {
    title <- "Australian Ibis Occupancy Probability"
    brks <- seq(0, 1, length.out = 21)
    lbl_brks <- seq(0, 1, length.out = 11) %>% 
      round(2)
  } else {
    title <- "Australian Ibis Occupancy Uncertainty (SE)"
    mx <- ceiling(1000 * cellStats(r_plot, max)) / 1000
    brks <- seq(0, mx, length.out = 21)
    lbl_brks <- seq(0, mx, length.out = 11) %>% 
      round(3)
  }
  pal <- abundance_palette(length(brks) - 1)
  plot(r_plot, 
       col = pal, breaks = brks, 
       maxpixels = ncell(r_plot),
       legend = FALSE, add = TRUE)
  
  # borders
  plot(ne_state_lines, col = "#ffffff", lwd = 0.75, add = TRUE)
  box()
  
  # legend
  par(new = TRUE, mar = c(0, 0, 0, 0))
  image.plot(zlim = range(brks), legend.only = TRUE, 
             breaks = brks, col = pal,
             smallplot = c(0.89, 0.92, 0.25, 0.75),
             horizontal = FALSE,
             axis.args = list(at = lbl_brks, labels = lbl_brks,
                              fg = "black", col.axis = "black",
                              cex.axis = 0.75, lwd.ticks = 0.5,
                              padj = 0),
             legend.args = list(text = title,
                                side = 2, col = "black",
                                cex = 1, line = 0))
}


## ----occupancy-model-select-dredge, class.source="livecode"--------------
# get list of all possible terms, then subset to those we want to keep

# fit all possible combinations of the occupancy covariates

# model comparison

### LIVE CODE ###


## ----occupancy-model-select-average, class.source="livecode"-------------
# select models with the most suport for model averaging


# average models based on model weights 

### LIVE CODE ###


## ----occupancy-select-detection-define-----------------------------------
# define occupancy model formula with only effort in detection submodel
det_mod <- ~ time_observations_started + 
  duration_minutes + 
  effort_distance_km + 
  number_observers + 
  protocol_type ~ 
  pland_13 + pland_09 + pland_02 + pland_08

# create new formulae with landcover covariates added to the detection submodel
mods <- list(det_mod_null = det_mod, 
             det_mod_for = update.formula(det_mod, ~ . + pland_02 ~ .),
             det_mod_ws = update.formula(det_mod, ~ . + pland_08 ~ .),
             global = update.formula(det_mod, 
                                     ~ . + pland_02 + pland_08 ~ .)) %>% 
  # fit candidate models
  map(occu, data = occ_um)


## ----occupancy-select-detection-compare----------------------------------
mod_sel <- fitList(fits = mods) %>% 
  modSel()
mod_sel


## ----occupancy-select-detection-coef-------------------------------------
coef(occ_model) %>% 
  enframe() %>% 
  filter(str_detect(name, "pland_"))