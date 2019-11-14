---
output: html_document
editor_options: 
  chunk_output_type: console
---

# Encounter Rate {#encounter}

In this lesson, we'll estimate the **encounter rate** of Austrlian Ibis on eBird checklists in September in the Temperate and Subtropical Forest CMZ, where encounter rate is defined as the probability of an eBirder encountering a species on a standard eBird checklist. We'll be using random forests in this lesson, a machine learning technique that uses an ensemble of many decision trees, each of which is fit using a bootstrap sampled of the data. For the purposes of this tutorial, we'll treat the random forest as a black box method. 

Let's start by loading all the packages and data we'll need for this lesson.


```r
library(sf)
```

```
## Linking to GEOS 3.7.2, GDAL 2.4.2, PROJ 5.2.0
```

```r
library(raster)
```

```
## Loading required package: sp
```

```
## 
## Attaching package: 'raster'
```

```
## The following object is masked from 'package:magrittr':
## 
##     extract
```

```r
library(dggridR)
```

```
## Loading required package: rgdal
```

```
## rgdal: version: 1.4-4, (SVN revision 833)
##  Geospatial Data Abstraction Library extensions to R successfully loaded
##  Loaded GDAL runtime: GDAL 2.1.3, released 2017/20/01
##  Path to GDAL shared files: /Users/mes335/Library/R/3.6/library/rgdal/gdal
##  GDAL binary built with GEOS: FALSE 
##  Loaded PROJ.4 runtime: Rel. 4.9.3, 15 August 2016, [PJ_VERSION: 493]
##  Path to PROJ.4 shared files: /Users/mes335/Library/R/3.6/library/rgdal/proj
##  Linking to sp version: 1.3-1
```

```
## Loading required package: ggplot2
```

```
## Loading required package: dplyr
```

```
## 
## Attaching package: 'dplyr'
```

```
## The following objects are masked from 'package:raster':
## 
##     intersect, select, union
```

```
## The following objects are masked from 'package:stats':
## 
##     filter, lag
```

```
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
```

```r
library(lubridate)
```

```
## 
## Attaching package: 'lubridate'
```

```
## The following object is masked from 'package:base':
## 
##     date
```

```r
library(ranger)
library(scam)
```

```
## Loading required package: mgcv
```

```
## Loading required package: nlme
```

```
## 
## Attaching package: 'nlme'
```

```
## The following object is masked from 'package:dplyr':
## 
##     collapse
```

```
## The following object is masked from 'package:raster':
## 
##     getData
```

```
## This is mgcv 1.8-29. For overview type 'help("mgcv-package")'.
```

```
## This is scam 1.2-5.
```

```r
library(PresenceAbsence)
library(verification)
```

```
## Loading required package: fields
```

```
## Loading required package: spam
```

```
## Loading required package: dotCall64
```

```
## Loading required package: grid
```

```
## Spam version 2.3-0 (2019-09-13) is loaded.
## Type 'help( Spam)' or 'demo( spam)' for a short introduction 
## and overview of this package.
## Help for individual functions is also obtained by adding the
## suffix '.spam' to the function name, e.g. 'help( chol.spam)'.
```

```
## 
## Attaching package: 'spam'
```

```
## The following objects are masked from 'package:base':
## 
##     backsolve, forwardsolve
```

```
## Loading required package: maps
```

```
## See https://github.com/NCAR/Fields for
##  an extensive vignette, other supplements and source code
```

```
## Loading required package: boot
```

```
## Loading required package: CircStats
```

```
## Loading required package: MASS
```

```
## 
## Attaching package: 'MASS'
```

```
## The following object is masked from 'package:dplyr':
## 
##     select
```

```
## The following objects are masked from 'package:raster':
## 
##     area, select
```

```
## Loading required package: dtw
```

```
## Loading required package: proxy
```

```
## 
## Attaching package: 'proxy'
```

```
## The following object is masked from 'package:spam':
## 
##     as.matrix
```

```
## The following object is masked from 'package:raster':
## 
##     as.matrix
```

```
## The following objects are masked from 'package:stats':
## 
##     as.dist, dist
```

```
## The following object is masked from 'package:base':
## 
##     as.matrix
```

```
## Loaded dtw v1.21-3. See ?dtw for help, citation("dtw") for use in publication.
```

```r
library(edarf)
library(ebirdst)
library(fields)
library(gridExtra)
```

```
## 
## Attaching package: 'gridExtra'
```

```
## The following object is masked from 'package:dplyr':
## 
##     combine
```

```r
library(tidyverse)
```

```
## ── Attaching packages ────────────────────────────────────────────────────── tidyverse 1.2.1 ──
```

```
## ✔ tibble  2.1.3     ✔ purrr   0.3.3
## ✔ tidyr   1.0.0     ✔ stringr 1.4.0
## ✔ readr   1.3.1     ✔ forcats 0.4.0
```

```
## ── Conflicts ───────────────────────────────────────────────────────── tidyverse_conflicts() ──
## ✖ lubridate::as.difftime() masks base::as.difftime()
## ✖ nlme::collapse()         masks dplyr::collapse()
## ✖ gridExtra::combine()     masks dplyr::combine()
## ✖ lubridate::date()        masks base::date()
## ✖ tidyr::extract()         masks raster::extract(), magrittr::extract()
## ✖ dplyr::filter()          masks stats::filter()
## ✖ lubridate::intersect()   masks raster::intersect(), base::intersect()
## ✖ dplyr::lag()             masks stats::lag()
## ✖ purrr::map()             masks maps::map()
## ✖ MASS::select()           masks dplyr::select(), raster::select()
## ✖ purrr::set_names()       masks magrittr::set_names()
## ✖ lubridate::setdiff()     masks base::setdiff()
## ✖ lubridate::union()       masks raster::union(), base::union()
```

```r
# resolve namespace conflicts
select <- dplyr::select
projection <- raster::projection
map <- purrr::map

set.seed(1)

# ebird data
ebird <- read_csv("data/ebd_zf_sep_tst.csv") %>% 
  filter(common_name == "Australian Ibis") %>% 
  # year required to join to habitat data
  mutate(year = year(observation_date))
```

```
## Parsed with column specification:
## cols(
##   .default = col_character(),
##   observation_count = col_double(),
##   species_observed = col_logical(),
##   latitude = col_double(),
##   longitude = col_double(),
##   all_species_reported = col_logical(),
##   observation_date = col_date(format = ""),
##   year = col_double(),
##   day_of_year = col_double(),
##   time_observations_started = col_double(),
##   duration_minutes = col_double(),
##   effort_distance_km = col_double(),
##   number_observers = col_double()
## )
```

```
## See spec(...) for full column specifications.
```

```r
# modis habitat covariates
habitat <- read_csv("data/pland-elev_location-year.csv") %>% 
  mutate(year = as.integer(year))
```

```
## Parsed with column specification:
## cols(
##   .default = col_double(),
##   year_lc = col_character(),
##   locality_id = col_character()
## )
## See spec(...) for full column specifications.
```

```r
# combine ebird and habitat data
ebird_habitat <- inner_join(ebird, habitat, by = c("locality_id", "year"))

# prediction surface
pred_surface <- read_csv("data/pland-elev_prediction-surface.csv")
```

```
## Parsed with column specification:
## cols(
##   .default = col_double()
## )
## See spec(...) for full column specifications.
```

```r
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
```

## Data preparation {#encounter-prep}

As we learned in [Part I](#subsampling) of this workshop, spatiotemporal subsampling can reduce spatial and temporal bias, and class imbalance, provided we sample detections and non-detections separately. So, we'll apply subsampling prior to fitting the random forest model.

<div class="tip">
  <h2>Tip</h2>
  Sampling detections and non-detections separately will change the prevalence rate of the detections in the data. As a result, the estimated probability of occurrence based on these subsampled data will be larger than the true occurrence rate. When examining the outputs from the models it will be important to recall that we altered the prevalence rate at this stage. 
</div>


```r
# generate hexagonal grid with ~ 5 km betweeen cells
dggs <- dgconstruct(spacing = 5)
```

```
## Resolution: 13, Area (km^2): 31.9926151554038, Spacing (km): 5.58632116604266, CLS (km): 6.38233997895802
```

```r
# get hexagonal cell id and week number for each checklist
checklist_cell <- ebird_habitat %>% 
  mutate(cell = dgGEO_to_SEQNUM(dggs, longitude, latitude)$seqnum,
         year = year(observation_date),
         week = week(observation_date))
```


```{.r .livecode}
# sample one checklist per grid cell per week
# sample detection/non-detection independently 
ebird_ss <- checklist_cell %>% 
  group_by(species_observed, year, week, cell) %>% 
  sample_n(size = 1) %>% 
  ungroup()
```

<div class="exercise">
  <h2>Exercise</h2>
  For very rare species, a more drastic approach to dealing with class imbalance is needed: only subsampling the non-detections and keeping all the detections. How could you modify the above code to do this?
  
  <button class="solution">Solution</button>
<div class="solution-content">

There are several valid ways to code this in R, but here is our approach.


```r
split_det <- split(checklist_cell, checklist_cell$species_observed)
ebird_all_det <- split_det$`FALSE` %>% 
  group_by(year, week, cell) %>% 
  sample_n(size = 1) %>% 
  ungroup() %>% 
  bind_rows(split_det$`TRUE`)
```

This approach leads to many more detections being kept in the data. 


```r
sum(ebird_ss$species_observed)
```

```
## [1] 1533
```

```r
sum(ebird_all_det$species_observed)
```

```
## [1] 3041
```

However, some of the extra detections we have with this approach are in the same 5km cell and the same week, they may not be independent.  There are trade-offs to many of these decisions about post-hoc sampling. 

</div>
</div>

In preparation for modeling, we'll select only the the columns that will be used as predictors in the model. We include both habitat predictors, which we expect to influence whether a species is present at a site, and also effort predictors to help control for variation in detectability.


```{.r .livecode}
# select covariates for model
ebird_ss <- ebird_ss %>% 
  select(species_observed,
         year, day_of_year,
         time_observations_started, duration_minutes,
         effort_distance_km, number_observers, 
         starts_with("pland_"),
         starts_with("elevation_")) %>% 
  drop_na()
```

Finally, we'll hold 20% of the data aside so we have an independent test set, which we can later use to assess the performance of our model.


```r
# split 80/20
ebird_split <- ebird_ss %>% 
  split(if_else(runif(nrow(.)) <= 0.8, "train", "test"))
```

<div class="exercise">
  <h2>Exercise</h2>
  How would you modify the above code to include 25% of data in the test set?
  
  <button class="solution">Solution</button>
<div class="solution-content">


```r
ebird_split_25 <- ebird_ss %>% 
  split(if_else(runif(nrow(.)) <= 0.75, "train", "test"))
```

</div>
</div>

## Random forests {#encounter-rf}

Random forests are an excellent, general purpose machine learning method suitable for modeling encounter rate in a wide variety of scenarios. To address the issue of class imbalance, we'll use a **balanced random forest** approach, a modification of the traditional random forest algorithm specifically designed to handle scenarios in which one class (in our case: species detections) is much more common than the other (non-detections). To implement a balanced random forest, we'll first need to calculate the frequency of detections (the smaller class).


```r
detection_freq <- mean(ebird_split$train$species_observed)
```

Now we can use the [`ranger`](https://github.com/imbs-hl/ranger) package to fit a random forest model to the eBird data.


```{.r .livecode}
# ranger requires a factor response to do classification
ebird_split$train$species_observed <- factor(ebird_split$train$species_observed)
# grow random forest
rf <- ranger(formula =  species_observed ~ ., 
             data = ebird_split$train,
             importance = "impurity",
             probability = TRUE,
             replace = TRUE, 
             sample.fraction = c(detection_freq, detection_freq))
```

### Calibration {#encounter-rf-cal}

For various reasons, the predicted probabilities from models do not always align with the observed frequencies of detections. We'll address this mismatch using **model calibration**, which aligns the estimated probabilities to the observed frequencies. In particular, to calibrate our model results, we predict encounter rate for each checklist in the training set, then fit a binomial Generalized Additive Model (GAM) with the real observations as the response and the predicted encounter rate as the predictor variable.


```r
# make predictions on training data
occ_pred <- rf$predictions[, 2]
# convert the observered response back to a numeric value from factor
occ_obs <- ebird_split$train$species_observed %>% 
  as.logical() %>% 
  as.integer()
rf_pred_train <- tibble(obs = occ_obs, pred = occ_pred) %>% 
  drop_na()
```


```{.r .livecode}
# fit gam calibration model
# scam allows us to use constrained shapes for the smooths
calibration_model <- scam(obs ~ s(pred, k = 6, bs = "mpi"), 
                          gamma = 2,
                          data = rf_pred_train)
```

```
## Warning: partial match of 'dev' to 'deviance'

## Warning: partial match of 'dev' to 'deviance'

## Warning: partial match of 'dev' to 'deviance'

## Warning: partial match of 'dev' to 'deviance'

## Warning: partial match of 'dev' to 'deviance'

## Warning: partial match of 'dev' to 'deviance'

## Warning: partial match of 'dev' to 'deviance'

## Warning: partial match of 'dev' to 'deviance'

## Warning: partial match of 'dev' to 'deviance'

## Warning: partial match of 'dev' to 'deviance'

## Warning: partial match of 'dev' to 'deviance'

## Warning: partial match of 'dev' to 'deviance'

## Warning: partial match of 'dev' to 'deviance'

## Warning: partial match of 'dev' to 'deviance'

## Warning: partial match of 'dev' to 'deviance'

## Warning: partial match of 'dev' to 'deviance'
```


```r
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
```

```
## Warning: Removed 4 rows containing missing values (geom_path).
```

![plot of chunk encounter-rf-cal-plot](figure/encounter-rf-cal-plot-1.png)

### Model assessment {#encounter-rf-assess}

To assess model quality, we'll validate the model’s ability to predict the observed patterns of occurrence using independent validation data (i.e. the 20% test data set that we removed earlier). We'll use a range of predictive performance metrics to compare the predictions to the actual observations. Different performance metrics reveal different aspects of the data and we can choose to emphasise some over others, depending on the specific goals of our analysis. For example, if we want to minimise the number of false negatives (i.e. we don't want to miss any places where the species actually occurs), we would focus on sensitivity. 


```r
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
rf_assessment <- tibble(
  model = c("RF", "Calibrated RF"),
  mse = c(mse_fit, mse_cal),
  sensitivity = c(metrics_fit$sensitivity, metrics_cal$sensitivity),
  specificity = c(metrics_fit$specificity, metrics_cal$specificity),
  auc = c(metrics_fit$AUC, metrics_cal$AUC),
  kappa = c(metrics_fit$Kappa, metrics_cal$Kappa)
)
knitr::kable(rf_assessment, digits = 3)
```



|model         |   mse| sensitivity| specificity|  auc| kappa|
|:-------------|-----:|-----------:|-----------:|----:|-----:|
|RF            | 0.163|       0.606|       0.847| 0.82| 0.434|
|Calibrated RF | 0.134|       0.592|       0.855| 0.82| 0.434|

## Habitat associations {#encounter-habitat}

From the random forest model, we can glean two important sources of information about the association between Austrlian Ibis detection and features of their local environment. First, **predictor importance** is a measure of the predictive power of each covariate, and is calculated as a byproduct of fitting a random forest model. Second, **partial dependence** plots estimate the marginal effect of one predictor holding all other predictors constant.

### Predictor importance {#encounter-habitat-pi}

During the process of fitting a random forest model, some variables are removed at each node of the trees that make up the random forest. **Predictor importance** is based on the mean decrease in accuracy of the model when a given covariate is not used.


```r
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
```

![plot of chunk encounter-habitat-pi](figure/encounter-habitat-pi-1.png)

<div class="tip">
  <h2>Tip</h2>
  Consult the file `data/mcd12q1_classes.csv` for a key to the different `pland_` variables.

```
## Parsed with column specification:
## cols(
##   class = col_double(),
##   name = col_character(),
##   description = col_character()
## )
```



| class|name                                |
|-----:|:-----------------------------------|
|     0|Water bodies                        |
|     1|Evergreen Needleleaf Forests        |
|     2|Evergreen Broadleaf Forests         |
|     3|Deciduous Needleleaf Forests        |
|     4|Deciduous Broadleaf Forests         |
|     5|Mixed Forests                       |
|     6|Closed Shrublands                   |
|     7|Open Shrublands                     |
|     8|Woody Savannas                      |
|     9|Savannas                            |
|    10|Grasslands                          |
|    11|Permanent Wetlands                  |
|    12|Croplands                           |
|    13|Urban and Built-up Lands            |
|    14|Cropland/Natural Vegetation Mosaics |
|    15|Non-Vegetated Lands                 |
|   255|Unclassified                        |
</div>

### Partial dependence {#encounter-habitat-pd}

**Partial dependence** plots show the marginal effect of a given predictor on encounter rate averaged across the other predictors. We'll use the R package `edarf` to construct partial dependence plots for the most important predictors.


```r
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
```

![plot of chunk encounter-habitat-pd](figure/encounter-habitat-pd-1.png)

## Prediction {#encounter-predict}

Finally, we can use the calibrated random forest model to make a map of Austrlian Ibis encounter rate in the CMZ! The data package contains a prediction surface consisting of the PLAND habitat covariates summarized on a regular grid of points across the CMZ. We'll make predictions of encounter rate at these points. However, first we need to bring effort variables into this prediction surface. We'll make predictions for a **standard eBird checklist**: a 1 km, 1 hour traveling count at the peak time of day for detecting this species.

To find the time of day with the highest detection probability, we can look for the peak of the partial dependence plot, constraining the search to times of day for which there are enough data to make reasonable predictions (hours with at least 1% of checklists).


```r
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
```

```
## [1] 5.320922
```



Based on this analysis, the best time for detecting Austrlian Ibis is at 5:19 AM. Now we can use this time to make predictions.


```{.r .livecode}
# add effort covariates to prediction 
pred_surface_eff <- pred_surface %>% 
  mutate(observation_date = ymd("2019-09-15"),
         year = year(observation_date),
         day_of_year = yday(observation_date),
         time_observations_started = t_peak,
         duration_minutes = 60,
         effort_distance_km = 1,
         number_observers = 1)

# predict
pred_rf <- predict(rf, data = pred_surface_eff, type = "response")
pred_rf <- pred_rf$predictions[, 2]
# apply calibration models
pred_rf_cal <- predict(calibration_model, 
                       data.frame(pred = pred_rf), 
                       type = "response")
# add to prediction surface
pred_er <- bind_cols(pred_surface_eff, encounter_rate = pred_rf_cal) %>% 
  select(latitude, longitude, encounter_rate) %>% 
  mutate(encounter_rate = pmin(pmax(encounter_rate, 0), 1))
```

Next, we'll convert this data frame to spatial features using `sf`, then rasterize the points using the prediction surface raster template.


```{.r .livecode}
# rasterize predictions
r_pred <- pred_er %>% 
  # convert to spatial features
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
  st_transform(crs = projection(r)) %>% 
  # rasterize
  rasterize(r)
r_pred <- r_pred[[-1]]
```

Finally, we can map these predictions!


```r
# project predictions
r_pred_proj <- projectRaster(r_pred, crs = map_proj$proj4string, method = "ngb")

par(mar = c(0.25, 0.25, 0.25, 3.5))
# set up plot area
plot(cmz, col = NA, border = NA)
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
plot(ne_state_lines, col = "#ffffff", lwd = 0.75, add = TRUE)
box()

# legend
par(new = TRUE, mar = c(0, 0, 0, 0))
title <- "Austrlian Ibis Encounter Rate"
image.plot(zlim = range(brks), legend.only = TRUE, 
           col = pal, breaks = brks,
           smallplot = c(0.91, 0.94, 0.25, 0.75),
           horizontal = FALSE,
           axis.args = list(at = lbl_brks, labels = lbl_brks,
                            fg = "black", col.axis = "black",
                            cex.axis = 0.75, lwd.ticks = 0.5,
                            padj = -1.5),
           legend.args = list(text = title,
                              side = 2, col = "black",
                              cex = 1, line = 0))
```

```
## Warning: partial match of 'z' to 'zlim'
```

![plot of chunk encounter-predict-map](figure/encounter-predict-map-1.png)

## Exercises {#encouter-exercises}

Now that you've completed this lesson, try modifying your script to complete at least one of the following exercises:

1. How does changing the subsampling grid cell size affect the model performance?

2. What happens to the predictions if you make them for an eBirder traveling further than 1 km, or birding for longer than 1 hour?

3. Filter the data to only shorter duration checklists or shorter distances traveled. How does this affect model performance?

4. An alternative approach to dealing with class imbalance, is to grid sample only the non-detections, while keeping all the detections. Try this subsampling approach and see what the affect is on the predictive performance metrics.

5. Try producing a map of encounter rate using one of the other species in the example data.


