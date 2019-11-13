## ----abundance-data------------------------------------------------------
library(lubridate)
library(sf)
library(raster)
library(dggridR)
library(pdp)
library(edarf)
library(mgcv)
library(fitdistrplus)
library(viridis)
library(fields)
library(tidyverse)
# resolve namespace conflicts
select <- dplyr::select
map <- purrr::map
projection <- raster::projection

set.seed(1)

# ebird data
ebird <- read_csv("data/ebd_woothr_june_bcr27_zf.csv") %>% 
  mutate(protocol_type = factor(protocol_type, 
                                levels = c("Stationary" , "Traveling"))) %>%
  # remove observations with no count
  filter(!is.na(observation_count))

# modis habitat covariates
habitat <- read_csv("data/pland-elev_location-year.csv") %>% 
  mutate(year = as.integer(year))

# combine ebird and habitat data
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


## ----abundance-nocount-sol-----------------------------------------------
read_csv("data/ebd_woothr_june_bcr27_zf.csv") %>% 
  summarize(n_total = n(),
            n_nocount = sum(is.na(observation_count)),
            prop = mean(is.na(observation_count)))


## ----abundance-prep-sss, results="hide"----------------------------------
# generate hexagonal grid with ~ 5 km betweeen cells
dggs <- dgconstruct(spacing = 5)
# get hexagonal cell id and week number for each checklist
checklist_cell <- ebird_habitat %>% 
  mutate(cell = dgGEO_to_SEQNUM(dggs, longitude, latitude)$seqnum,
         week = week(observation_date))
# sample one checklist per grid cell per week
# sample detection/non-detection independently 
ebird_ss <- checklist_cell %>% 
  group_by(species_observed, year, week, cell) %>% 
  sample_n(size = 1) %>% 
  ungroup() %>% 
  select(-cell, -week)


## ----abundance-prep-tt---------------------------------------------------
# split 80/20
ebird_split <- ebird_ss %>% 
  split(if_else(runif(nrow(.)) <= 0.8, "train", "test"))


## ----abundance-model-dist, fig.asp=0.65----------------------------------
par(mfrow = c(1, 2))
# counts with zeros
hist(ebird_ss$observation_count, main = "Histogram of counts", 
     xlab = "Observed count")
# counts without zeros
pos_counts <- keep(ebird_ss$observation_count, ~ . > 0)
hist(pos_counts, main = "Histogram of counts > 0", 
     xlab = "Observed non-zero count")
prop_zero <- sum(ebird_ss$observation_count == 0) / nrow(ebird_ss)
prop_zero


## ----abundance-model-formula, class.source="livecode"--------------------
# gam formula

### LIVE CODE ###

# explicitly specify where the knots should occur for time_observations_started
# this ensures that the cyclic spline joins the variable at midnight
# this won't happen by default if there are no data near midnight

### LIVE CODE ###


## ----abundance-model-formula-sol-----------------------------------------
gam_formula_elev <- observation_count ~ s(day_of_year, k = 5) + 
  s(duration_minutes, k = 5) + 
  s(effort_distance_km, k = 5) + 
  s(number_observers, k = 5) + 
  s(pland_04, k = 5) + 
  s(pland_05, k = 5) + 
  s(pland_12, k = 5) + 
  s(pland_13, k = 5) + 
  s(elevation_median, k = 5) + 
  s(elevation_sd, k = 5) + 
  protocol_type + 
  s(time_observations_started, bs = "cc", k = k_time)


## ----abundance-model-gams, class.source="livecode"-----------------------
# zero-inflated poisson
### LIVE CODE ###

# negative binomial
### LIVE CODE ###

# tweedie distribution
### LIVE CODE ###


## ----abundance-assess-pred-----------------------------------------------
obs_count <- select(ebird_split$test, obs = observation_count)

# presence probability is on the complimentary log-log scale
# we can get the inverse link function with
inv_link <- binomial(link = "cloglog")$linkinv
# combine ziplss presence and count predictions
m_ziplss_pred <- predict(m_ziplss, ebird_split$test, type = "link") %>% 
  as.data.frame() %>%
  transmute(family = "Zero-inflated Poisson",
            pred = inv_link(V2) * exp(V1)) %>% 
  bind_cols(obs_count)

m_nb_pred <- predict(m_nb, ebird_split$test, type = "response") %>% 
  tibble(family = "Negative Binomial", pred = .) %>% 
  bind_cols(obs_count)

m_tw_pred <- predict(m_tw, ebird_split$test, type = "response") %>% 
  tibble(family = "Tweedie", pred = .) %>% 
  bind_cols(obs_count)

# combine predictions from all three models
test_pred <- bind_rows(m_ziplss_pred, m_nb_pred, m_tw_pred) %>% 
  mutate(family = as_factor(family))


## ----abundance-assess-metrics, class.source="livecode"-------------------
# spearman’s rank correlation
### LIVE CODE ###


## ----abundance-assess-metrics-sol----------------------------------------
test_pred %>% 
  filter(obs > 0) %>% 
  group_by(family) %>% 
  summarise(n_under = sum(obs / pred < 10),
            pct_under = mean(obs / pred < 10)) %>% 
  ungroup()


## ----abundance-assess-plot-----------------------------------------------
# plot predicted vs. observed
ticks <- c(0, 1, 10, 100, 1000)
mx <- round(max(test_pred$obs))
ggplot(test_pred) +
  aes(x = log10(obs + 1), 
      y = log10(pred + 1)) +
  geom_jitter(alpha = 0.2, height = 0) +
  # y = x line
  geom_abline(slope = 1, intercept = 0, alpha = 0.5) +
  # area where counts off by a factor of 10
  geom_area(data = tibble(x = log10(seq(0, mx - 1) + 1), 
                          y = log10(seq(0, mx - 1) / 10 + 1)),
            mapping = aes(x = x, y = y),
            fill = "red", alpha = 0.2) +
  # loess fit
  geom_smooth(method = "loess", 
              method.args = list(span = 2 / 3, degree = 1)) +
  scale_x_continuous(breaks = log10(ticks + 1), labels = ticks) +
  scale_y_continuous(breaks = log10(ticks + 1), labels = ticks) +
  labs(x = "Observed count",
       y = "Predicted count") +
  facet_wrap(~ family, nrow = 1)


## ----abundance-assess-decision-------------------------------------------
pred_model <- m_nb


## ----abundance-model-cov-plot, fig.asp=1---------------------------------
par(mai = c(0.75, 0.25, 0.2, 0.25))
plot(pred_model, pages = 1, ylab = "")


## ----abundance-predict-peak----------------------------------------------
# create a dataframe of covariates with a range of start times
seq_tod <- seq(0, 24, length.out = 300)
tod_df <- ebird_split$train %>% 
  # find average pland habitat covariates
  select(starts_with("pland")) %>% 
  summarize_all(mean, na.rm = TRUE) %>% 
  ungroup() %>% 
  # use standard checklist
  mutate(day_of_year = yday(ymd(str_glue("{max_lc_year}-06-15"))),
         duration_minutes = 60,
         effort_distance_km = 1,
         number_observers = 1,
         protocol_type = "Traveling") %>% 
  cbind(time_observations_started = seq_tod)

# predict at different start times
pred_tod <- predict(pred_model, newdata = tod_df, 
                    type = "link", 
                    se.fit = TRUE) %>% 
  as_tibble() %>% 
  # calculate backtransformed confidence limits
  transmute(time_observations_started = seq_tod,
            pred = pred_model$family$linkinv(fit),
            pred_lcl = pred_model$family$linkinv(fit - 1.96 * se.fit),
            pred_ucl = pred_model$family$linkinv(fit + 1.96 * se.fit))

# find optimal time of day
t_peak <- pred_tod$time_observations_started[which.max(pred_tod$pred_lcl)]

# plot the partial dependence plot
ggplot(pred_tod) +
  aes(x = time_observations_started, y = pred,
      ymin = pred_lcl, ymax = pred_ucl) +
  geom_ribbon(fill = "grey80", alpha = 0.5) +
  geom_line() +
  geom_vline(xintercept = t_peak, color = "blue", linetype = "dashed") +
  labs(x = "Hours since midnight",
       y = "Predicted relative abundance",
       title = "Effect of observation start time on Wood Thrush reporting",
       subtitle = "Peak detectability shown as dashed blue line")


## ----abundance-predict-readable, echo=FALSE------------------------------
human_time <- str_glue("{h}:{m} {ap}", 
                       h = floor(t_peak),
                       m = str_pad(round((t_peak %% 1) * 60), 2, pad = "0"),
                       ap = ifelse(t_peak > 12, "PM", "AM"))


## ----abundance-predict-effort--------------------------------------------
# add effort covariates to prediction surface
pred_surface_eff <- pred_surface %>% 
  mutate(day_of_year = yday(ymd(str_glue("{max_lc_year}-06-15"))),
         time_observations_started = t_peak,
         duration_minutes = 60,
         effort_distance_km = 1,
         number_observers = 1,
         protocol_type = "Traveling")


## ----abundance-predict-pred, class.source="livecode"---------------------
# predict
### LIVE CODE ###


## ----abundance-predict-rasterize-----------------------------------------
r_pred <- pred %>% 
  # convert to spatial features
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
  select(abd, abd_se) %>% 
  st_transform(crs = projection(r)) %>% 
  # rasterize
  rasterize(r)
r_pred <- r_pred[[-1]]


## ----abundance-predict-map, fig.asp=1.236--------------------------------
# any expected abundances below this threshold are set to zero
zero_threshold <- 0.05

# project predictions
r_pred_proj <- projectRaster(r_pred, crs = map_proj$proj4string, method = "ngb")

par(mfrow = c(2, 1))
for (nm in names(r_pred)) {
  r_plot <- r_pred_proj[[nm]]
  
  par(mar = c(3.5, 0.25, 0.25, 0.25))
  # set up plot area
  plot(bcr, col = NA, border = NA)
  plot(ne_land, col = "#dddddd", border = "#888888", lwd = 0.5, add = TRUE)
  
  # modified plasma palette
  plasma_rev <- rev(plasma(25, end = 0.9))
  gray_int <- colorRampPalette(c("#dddddd", plasma_rev[1]))
  pal <- c(gray_int(4)[2], plasma_rev)
  
  # abundance vs. se
  if (nm == "abd") {
    title <- "Wood Thrush Relative Abundance"
    # set very low values to zero
    r_plot[r_plot <= zero_threshold] <- NA
    # log transform
    r_plot <- log10(r_plot)
    # breaks and legend
    mx <- ceiling(100 * cellStats(r_plot, max)) / 100
    mn <- floor(100 * cellStats(r_plot, min)) / 100
    brks <- seq(mn, mx, length.out = length(pal) + 1)
    lbl_brks <- sort(c(-2:2, mn, mx))
    lbls <- round(10^lbl_brks, 2)
  } else {
    title <- "Wood Thrush Abundance Uncertainty (SE)"
    # breaks and legend
    mx <- ceiling(1000 * cellStats(r_plot, max)) / 1000
    mn <- floor(1000 * cellStats(r_plot, min)) / 1000
    brks <- seq(mn, mx, length.out = length(pal) + 1)
    lbl_brks <- seq(mn, mx, length.out = 5)
    lbls <- round(lbl_brks, 2)
  }
  
  # abundance
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
  image.plot(zlim = range(brks), legend.only = TRUE, col = pal,
             smallplot = c(0.25, 0.75, 0.06, 0.09),
             horizontal = TRUE,
             axis.args = list(at = lbl_brks, 
                              labels = lbls,
                              fg = "black", col.axis = "black",
                              cex.axis = 0.75, lwd.ticks = 0.5,
                              padj = -1.5),
             legend.args = list(text = title,
                                side = 3, col = "black",
                                cex = 1, line = 0))
}

# Exercises ----

# 1. Refit the model without effort variables and see how model performance
# changes.

# 2. Predict from the same model for checklists of 10 minutes duration, instead
# of 1 hour. Compare the results and consider how the interpretation changes.

# 3. Change the degrees of freedom for the covariate smooths and compare the
# fitted relationships.

# 4. Refit the model with a random forest. Compare the predictions and the
# fitted relationships with covariates.

# 5. Compare the encounter rate map to the relative abundance map.

# 6. Fit an encounter rate random forest model then use the predicted encounter
# rate as a new covariate in the abundance model. Compare the model performance.