## ----encounter-----------------------------------------------------------
library(sf)
library(raster)
library(dggridR)
library(lubridate)
library(ranger)
library(scam)
library(PresenceAbsence)
library(verification)
library(edarf)
library(ebirdst)
library(fields)
library(gridExtra)
library(tidyverse)
# resolve namespace conflicts
select <- dplyr::select
projection <- raster::projection
map <- purrr::map

set.seed(1)

# ebird data
ebird <- read_csv("data/ebd_woothr_june_bcr27_zf.csv") %>% 
  # year required to join to habitat data
  mutate(year = year(observation_date))

# modis habitat covariates
habitat <- read_csv("data/pland-elev_location-year.csv") %>% 
  mutate(year = as.integer(year))

# combine ebird and habitat data
ebird_habitat <- inner_join(ebird, habitat, by = c("locality_id", "year"))

# prediction surface
pred_surface <- read_csv("data/pland-elev_prediction-surface.csv")
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


## ----encounter-prep-ss-grid----------------------------------------------
# generate hexagonal grid with ~ 5 km betweeen cells
dggs <- dgconstruct(spacing = 5)
# get hexagonal cell id and week number for each checklist
checklist_cell <- ebird_habitat %>% 
  mutate(cell = dgGEO_to_SEQNUM(dggs, longitude, latitude)$seqnum,
         year = year(observation_date),
         week = week(observation_date))


## ----encounter-prep-ss, class.source="livecode"--------------------------
# sample one checklist per grid cell per week
# sample detection/non-detection independently 

### LIVE CODE ###


## ----encounter-prep-select, class.source="livecode"----------------------
# select covariates for model

### LIVE CODE ###


## ----encounter-prep-tt---------------------------------------------------
# split 80/20
ebird_split <- ebird_ss %>% 
  split(if_else(runif(nrow(.)) <= 0.8, "train", "test"))


## ----encounter-rf-detfreq------------------------------------------------
detection_freq <- mean(ebird_split$train$species_observed)


## ----encounter-rf-fit, class.source="livecode"---------------------------
# ranger requires a factor response to do classification

### LIVE CODE ###

# grow random forest

### LIVE CODE ###


## ----encounter-rf-cal----------------------------------------------------
# make predictions on training data
occ_pred <- rf$predictions[, 2]
# convert the observered response back to a numeric value from factor
occ_obs <- ebird_split$train$species_observed %>% 
  as.logical() %>% 
  as.integer()
rf_pred_train <- tibble(obs = occ_obs, pred = occ_pred) %>% 
  drop_na()


## ----encounter-rf-cal-scam, class.source="livecode"----------------------
# fit gam calibration model
# scam allows us to use constrained shapes for the smooths

### LIVE CODE ###


## ----encounter-rf-cal-plot-----------------------------------------------
# plot calibration curve
cal_pred <- tibble(pred = seq(0, 1, length.out = 100))
cal_pred <- predict(calibration_model, cal_pred, type = "response") %>% 
  bind_cols(cal_pred, calibrated = .)
ggplot(cal_pred) +
  aes(x = pred, y = calibrated) +
  geom_line() +
  labs(x = "RF prediction",
       y = "Calibrated prediction",
       title = "Calibration model") + 
  xlim(0, 1) + ylim(0, 1)


## ----encounter-rf-assess-------------------------------------------------
# predict on test data using calibrated model
p_fitted <- predict(rf, data = ebird_split$test, type = "response")
# extract probability of detection
p_fitted <- p_fitted$predictions[, 2]
p_calibrated <- predict(calibration_model, 
                        newdata = tibble(pred = p_fitted), 
                        type = "response")
rf_pred_test <- data.frame(id = seq_along(p_calibrated),
                           # actual detection/non-detection
                           obs = ebird_split$test$species_observed,
                           # uncalibrated prediction
                           fit = p_fitted,
                           # calibrated prediction
                           cal = p_calibrated) %>%
  # constrain probabilities to 0-1
  mutate(cal = pmin(pmax(cal, 0), 1)) %>% 
  drop_na()

# mean squared error (mse)
mse_fit <- mean((rf_pred_test$obs - rf_pred_test$fit)^2, na.rm = TRUE)
mse_cal <- mean((rf_pred_test$obs - rf_pred_test$cal)^2, na.rm = TRUE)

# pick threshold to maximize kappa
opt_thresh <- optimal.thresholds(rf_pred_test, opt.methods = "MaxKappa")

# calculate accuracy metrics: auc, kappa, sensitivity, specificity, brier
metrics_fit <- rf_pred_test %>% 
  select(id, obs, fit) %>% 
  presence.absence.accuracy(threshold = opt_thresh$fit, 
                            na.rm = TRUE, 
                            st.dev = FALSE)
metrics_cal <- rf_pred_test %>% 
  select(id, obs, cal) %>% 
  presence.absence.accuracy(threshold = opt_thresh$cal, 
                            na.rm = TRUE, 
                            st.dev = FALSE)

# combine various performance metrics together
tibble(
  model = c("RF", "Calibrated RF"),
  mse = c(mse_fit, mse_cal),
  sensitivity = c(metrics_fit$sensitivity, metrics_cal$sensitivity),
  specificity = c(metrics_fit$specificity, metrics_cal$specificity),
  auc = c(metrics_fit$AUC, metrics_cal$AUC),
  kappa = c(metrics_fit$Kappa, metrics_cal$Kappa)
)


## ----encounter-habitat-pi------------------------------------------------
pi <- enframe(rf$variable.importance, "predictor", "importance")
# plots
ggplot(pi) + 
  aes(x = fct_reorder(predictor, importance), y = importance) +
  geom_col() +
  geom_hline(yintercept = 0, size = 2, colour = "#555555") +
  scale_y_continuous(expand = c(0, 0)) +
  coord_flip() +
  labs(x = NULL, 
       y = "Predictor Importance (Gini Index)") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        panel.grid.major.x = element_line(colour = "#cccccc", size = 0.5))


## ----encounter-habitat-pi-pland, echo=FALSE------------------------------
read_csv("data/mcd12q1_classes.csv") %>% 
  select(class, name)


## ----encounter-habitat-pd------------------------------------------------
# top 9 predictors other than date
top_pred <- pi %>% 
  filter(!predictor %in% c("year", "day_of_year")) %>% 
  top_n(n = 9, wt = importance) %>% 
  arrange(desc(importance))

# calculate partial dependence for each predictor
pd <- top_pred %>% 
  mutate(pd = map(predictor, partial_dependence, 
                  fit = rf, data = ebird_split$train),
         pd = map(pd, ~ .[, c(1, 3)]),
         pd = map(pd, set_names, nm = c("value",  "encounter_rate"))) %>% 
  unnest(cols = pd)

# calibrate predictions
pd$encounter_rate <- predict(calibration_model, 
                             newdata = tibble(pred = pd$encounter_rate), 
                             type = "response") %>% 
  as.numeric()
# constrain probabilities to 0-1
pd$encounter_rate <- pmin(pmax(pd$encounter_rate, 0), 1)
  
# plot
ggplot(pd) +
  aes(x = value, y = encounter_rate) +
  geom_line() +
  geom_point() +
  scale_y_continuous(labels = scales::percent) +
  facet_wrap(~ as_factor(predictor), nrow = 3, scales = "free") +
  labs(x = NULL, y = "Encounter Rate") +
  theme_minimal() +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "grey60"),
        axis.ticks  = element_line(color = "grey60"))


## ----encounter-predict-time----------------------------------------------
# find peak time of day from partial dependence
pd_time <- partial_dependence(rf, 
                              vars = "time_observations_started", 
                              # make estimates at 30 minute intervals
                              # use the entire training dataset for estimation
                              n = c(24 * 2, nrow(ebird_split$train)), 
                              data = ebird_split$train) %>% 
  select(time_observations_started, encounter_rate = "TRUE")

# hours with at least 1% of checklists
search_hours <- ebird_split$train %>% 
  mutate(hour = floor(time_observations_started)) %>%
  count(hour) %>% 
  mutate(pct = n / sum(n)) %>% 
  filter(pct >= 0.01)

# constrained peak time
t_peak <- pd_time %>% 
  filter(floor(time_observations_started) %in% search_hours$hour) %>% 
  top_n(1, wt = desc(time_observations_started)) %>% 
  pull(time_observations_started)
t_peak


## ----encounter-predict-readable, echo = FALSE----------------------------
human_time <- str_glue("{h}:{m} {ap}", 
                       h = floor(t_peak),
                       m = str_pad(round((t_peak %% 1) * 60), 2, pad = "0"),
                       ap = ifelse(t_peak > 12, "PM", "AM"))


## ----encounter-predict-effort, class.source="livecode"-------------------
# add effort covariates to prediction 

### LIVE CODE ###

# predict

### LIVE CODE ###

# apply calibration models

### LIVE CODE ###

# add to prediction surface

### LIVE CODE ###


## ----encounter-predict-rasterize, class.source="livecode"----------------
# rasterize predictions

### LIVE CODE ###


## ----encounter-predict-map-----------------------------------------------
# project predictions
r_pred_proj <- projectRaster(r_pred, crs = map_proj$proj4string, method = "ngb")

par(mar = c(3.5, 0.25, 0.25, 0.25))
# set up plot area
plot(bcr, col = NA, border = NA)
plot(ne_land, col = "#dddddd", border = "#888888", lwd = 0.5, add = TRUE)

# encounter rate
r_max <- ceiling(10 * cellStats(r_pred_proj, max)) / 10
brks <- seq(0, r_max, by = 0.025)
lbl_brks <- seq(0, r_max, by = 0.1)
# ebird status and trends color palette
pal <- abundance_palette(length(brks) - 1)
plot(r_pred_proj, 
     col = pal, breaks = brks, 
     maxpixels = ncell(r_pred_proj),
     legend = FALSE, add = TRUE)

# borders
plot(bcr, border = "#000000", col = NA, lwd = 1, add = TRUE)
plot(ne_state_lines, col = "#ffffff", lwd = 0.75, add = TRUE)
plot(ne_country_lines, col = "#ffffff", lwd = 1.5, add = TRUE)
box()

# legend
par(new = TRUE, mar = c(0, 0, 0, 0))
title <- "Wood Thrush Encounter Rate"
image.plot(zlim = range(brks), legend.only = TRUE, 
           col = pal, breaks = brks,
           smallplot = c(0.25, 0.75, 0.06, 0.09),
           horizontal = TRUE,
           axis.args = list(at = lbl_brks, labels = lbl_brks,
                            fg = "black", col.axis = "black",
                            cex.axis = 0.75, lwd.ticks = 0.5,
                            padj = -1.5),
           legend.args = list(text = title,
                              side = 3, col = "black",
                              cex = 1, line = 0))


# Exercises ----

# 1. How does changing the subsampling grid cell size affect the model
# performance?

# 2. What happens to the predictions if you make them for an eBirder traveling
# further than 1 km, or birding for longer than 1 hour?

# 3. Filter the data to only shorter duration checklists or shorter distances
# traveled. How does this affect model performance?

# 4. An alternative approach to dealing with class imbalance, is to grid sample
# only the non-detections, while keeping all the detections. Try this
# subsampling approach and see what the affect is on the predictive performance
# metrics.