library(raster)
library(dggridR)
library(lubridate)
library(mgcv)
library(DHARMa)
library(PresenceAbsence)
library(tidyverse)
# resolve namespace conflicts
select <- dplyr::select
map <- purrr::map

species_list <- c("Little Corella", "Pied Butcherbird", 
                  "Australian Ibis", "Rainbow Bee-eater", 
                  "Leaden Flycatcher", "Rose Robin",
                  "Rufous Whistler", "Gray Shrikethrush",
                  "Gray Fantail", "Black-chinned Honeyeater", 
                  "Brown Treecreeper", "Jacky-winter",
                  "Channel-billed Cuckoo", "Red-capped Robin", 
                  "Yellow-faced Honeyeater")

#for (species in species_list) {
species <- species_list[7]
  set.seed(1)
  
  sp_code <- filter(auk::ebird_taxonomy, common_name == species) %>% 
    pull(species_code)
  out_dir <- file.path("output", sp_code)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # ebird data
  ebird <- read_csv("data/ebd_zf_sep_tst.csv") %>% 
    filter(common_name == species) %>% 
    # year required to join to habitat data
    mutate(year = year(observation_date))
  
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
    

  # generate hexagonal grid with ~ 5 km betweeen cells
  dggs <- dgconstruct(spacing = 5)
  # get hexagonal cell id and week number for each checklist
  checklist_cell <- ebird_habitat %>% 
    mutate(cell = dgGEO_to_SEQNUM(dggs, longitude, latitude)$seqnum,
           week = week(observation_date))

  # --------------------------------------------------------------------
  # now we aggregate the information within each grid cell. 
  # this takes advantage of all the information within each cell rather 
  # than just selecting one checklist. 
  # this approach of data aggregation can give the model more power to 
  # detect change, especially in areas with few grid cells. 


  # aggregate checklist information within each grid cell and week
  mean.na <- function(x) {mean(x, na.rm = TRUE)}
  sum.na <- function(x) {sum(x, na.rm = TRUE)}
  ebird_ag <- checklist_cell %>% 
    group_by(week, cell, year) %>% 
    summarise(mean_species_observed = mean.na(species_observed),
              mean_count = mean.na(observation_count),
              total_count = sum.na(observation_count),
              mean_longitude = mean.na(longitude),
              mean_latitude = mean.na(latitude),
              sum_duration_minutes = sum.na(duration_minutes), 
              sum_effort_distance_km = sum.na(effort_distance_km),
              mean_time_observations_started = mean.na(time_observations_started),
              wt_mean_time_observations_started = sum(time_observations_started*duration_minutes)/sum(duration_minutes),
              median_number_observers = median(number_observers, na.rm = TRUE),
              day_of_year = mean.na(day_of_year),
              median_pland_00 = median(pland_00),
              median_pland_01 = median(pland_01),             
              median_pland_02 = median(pland_02),             
              median_pland_04 = median(pland_04),             
              median_pland_05 = median(pland_05),             
              median_pland_06 = median(pland_06),             
              median_pland_07 = median(pland_07),             
              median_pland_08 = median(pland_08),             
              median_pland_09 = median(pland_09),             
              median_pland_10 = median(pland_10),             
              median_pland_11 = median(pland_11),             
              median_pland_12 = median(pland_12),             
              median_pland_13 = median(pland_13),             
              median_pland_14 = median(pland_14),
              median_pland_15 = median(pland_15),
              number_checklists = n(),
              pos_checklists = sum.na(species_observed)) %>%
    ungroup()

  # split into train and test by randomly selecting 80% of grid cells for training
  ebird_ag_split <- ebird_ag %>% 
    split(if_else(runif(nrow(.)) <= 0.8, "train", "test"))
  map_int(ebird_ag_split, nrow)
  

  # --------------------------------------------------------------------
  # fit an occurrence rate GAM to these data

  par(mfrow = c(2, 2))

  # the response here is the mean occurrence rate on checklists within each hexagonal grid cell
  # some grid cells have more checklists, so the occurrence rate will be more precise. this makes
  # us think about weighting towards those data points with more precision. however, if we did 
  # this, we would put greater weight on the cells that have been surveyed more often, which may
  # result in bias. 
  # so how do we account for the fact that some occurrence rates are imprecise, without weighting
  # too much towards popular birding hotspots. So there are two goals here: 1) weighting towards places
  # with more information; and 2) not weighting towards places with more information! 
  # there are a number of approaches we could take here, but here we have opted to down-weight 
  # cells with fewer than 5 checklists, but treat the information from all others equally. this is 
  # a compromise between the two goals. 

  # define the weights (there are many different and justifiable ways to do this)
  weights = apply(cbind(ebird_ag_split$train$number_checklists, 5), 1, min)

  # we know that occurrence rate varies across the landscape, and here we are modelling that with 
  # a 2-dimensional spline. 

  # Duchon splines are less prone to edge effects, so we use those here for the year and spatial 
  # smooth components, by specififying bs = "ds"
  # https://www.rdocumentation.org/packages/mgcv/versions/1.8-28/topics/smooth.construct.ds.smooth.spec

  # info on different types of smooth: 
  # https://www.rdocumentation.org/packages/mgcv/versions/1.8-28/topics/smooth.terms

  # We use REML as the method to fit the splines. See here for a brief description about why:
  # https://github.com/DistanceDevelopment/dsm/wiki/Why-is-the-default-smoothing-method-%22REML%22-rather-than-%22GCV.Cp%22%3F

  # We set gamma at 1.4, which reduces overfitting by the GAM function and is recommended by Simon Wood. 

  occ_gam <- gam(mean_species_observed ~ s(year, k=5, m=c(1, 0.5), bs = "ds") + 
              week + s(mean_latitude, mean_longitude, k = 50, bs = "ds", m=c(1, 0.5)) + 
              s(sum_duration_minutes, k = 5) + s(sum_effort_distance_km, k = 5) + 
              s(wt_mean_time_observations_started, k = 10) + s(median_number_observers, k = 4),
              method = "REML",
              gamma = 1.4, 
              weights = weights,
              data = ebird_ag_split$train, family = binomial)

  summary(occ_gam)

  # We see that the deviance explained is fairly low, suggesting there is a lot more variation that cannot be
  # captured by the covariates we included. For now, that's ok. We're here focussed on the general change over time. 

  # There is a handy function that helps us check the residuals of our model. 

  gam.check(occ_gam)

  # we see that the QQ plot (upper left) is close to the line, so this suggests a good fit. 
  # plots 2 and 3 show weird residual distributions, but this is typical of binomial models, so this is
  # not too concerning here. 
  # response vs fitted values in plot 4 shows some positive correlation, but it's hard to see detail, because there
  # are so many points overlaid 

  # we can get better residuals for a binomial model from a different method in a new package DHARMa

  simulationOutput <- simulateResiduals(fittedModel = occ_gam, n = 250)
  plot(simulationOutput)

  # we see good fit to the QQ plot and a close fit of residual vs predicted lines in plot 2. 

  # now we'll use this model to predict to the test dataset and plot the observed against predicted values
  pred_dat <- ebird_ag_split$test

  pred_test <- predict(occ_gam, newdata = pred_dat, type = "response", se.fit = FALSE)

  par(mfrow=c(1,2))
  plot(pred_test, pred_dat$mean_species_observed, cex = sqrt(apply(cbind(pred_dat$number_checklists, 40), 1, min))*0.5, 
    xlim=c(0, 1), ylim = c(0, 1), xlab = "predited occurrence rate", ylab = "observed occurrence rate")
  boxplot(pred_dat$mean_species_observed ~ cut(pred_test, seq(0, 1, by=0.1)))
  abline(-0.1, 0.1, col="red")

  # here we see reasonable correlation, but not great calibration. 

  # we'll use the functions in another package to explore the fit of the test data a little more
  presence.absence.summary(data.frame(id = 1:length(pred_test), obs = pred_dat$mean_species_observed, pred1 = pred_test))

  # to double check things look sensible, let's look at the fitted relationships with the effort covariates
  plot(occ_gam, pages = 1)

  # now that we've verified the model fit, let's explore the modelled change over time

  par(mfrow=c(1,1))
  plot(occ_gam, select = 1, ylim=c(-1, 1))

  # the y-axis here is a centred response on the link scale. so it's hard to interpret the absolute 
  # values. but the pattern is similar to the real scale. 

  # we'll convert this modelled trend to something a little more interpretable

  # --------------------------------------------------------------------
  # predictions for the trend across time

  mean_covariates <- ebird_ag %>%
    group_by() %>%
    summarise_all(mean, na.rm = TRUE) %>%
    ungroup() %>%
    select(-year) %>%
    as.data.frame() %>%
    mutate(id = 1)

  year_seq <- seq(min(ebird_ag_split$train$year), max(ebird_ag_split$train$year))
  nd <- data.frame(year = year_seq, id = 1) %>%
        left_join(mean_covariates, by = "id")

  expit <- function(x) { exp(x)/(1+exp(x)) }
  pred_occ_gam <- predict(occ_gam, newdata=nd, se.fit = TRUE)
  pred_occ <- data.frame(year = nd$year, est = expit(pred_occ_gam$fit),
    lcl = expit(pred_occ_gam$fit - 1.96*pred_occ_gam$se.fit),
    ucl = expit(pred_occ_gam$fit + 1.96*pred_occ_gam$se.fit))

  # plot out the raw reporting rates and the fitted gam over time
  plot(0, 0, 
    xlab = "year", ylab = "Encounter rate", 
    xlim = c(min(pred_occ$year), max(pred_occ$year)), ylim = c(0, max(pred_occ$ucl)*1.1), 
    col = "white")
  polygon(x = c(pred_occ$year, rev(pred_occ$year), pred_occ$year[1]), 
          y = c(pred_occ$lcl, rev(pred_occ$ucl), pred_occ$lcl[1]), 
          border = alpha("white", 0), col = "grey80")
  lines(pred_occ$year, pred_occ$est, pch = 16)

  # we see that we have very high uncertainty for the trend of this species in this region. 
  # we have no precision to understand whether the species is increasing, decreasing, or stable. 
  # we can say that we have no evidence of a decreasing or increasing trend, but we don't know 
  # whether this is due to a stable population, or low power. we could investigate this with 
  # a power analysis or simulation to understand more. 



  # ####################################################################
  # ABUNDANCE TREND
  # ####################################################################

  # Now we'll move onto estimating an abundance trend for the same species in the same region. 
  # We use the same dataset above and we've already calculated the mean abundance in each cell. 
  # For estimating changes in abundance, cells which have never recorded the species do not
  # provide us any information, and they can cause additional challenges because there are often
  # many zeroes. so we'll remove those from our analysis. 

  # first aggregate the counts within each site to find sites that never record the species
  ebird_site_total_positive <- ebird_ag %>%
      dplyr:::select(cell, total_count) %>%
      group_by(cell) %>%
      summarise_all(list(alltime = sum), na.rm = TRUE) %>%
      ungroup() %>%
      filter(alltime>0)

  # filter to only the sites with counts at some point and split into train and test
  ebird_ag_count <- ebird_ag %>% 
    filter(cell %in% ebird_site_total_positive$cell) %>%
    split(if_else(runif(nrow(.)) <= 0.8, "train", "test"))
  map_int(ebird_ag_count, nrow)

  # let's visualise the data
  par(mfrow=c(1,2))
  hist(ebird_ag_count$train$total_count, main = "total count per cell", xlab = "total count per cell")
  hist(ebird_ag_count$train$total_count / (ebird_ag_count$train$sum_duration_minutes/60), main = "mean count per hour", xlab = "count per hour")

  # as above, we'll fit a GAM to these data. as we have counts we need to use a different distribution. 
  # we'll use a negative binomial model 

  weights <- ifelse(ebird_ag_count$train$number_checklists>5, 5, ebird_ag_count$train$number_checklists) 

  abd_gam <- gam(total_count ~ s(year, k=5, m=2, bs = "ds") + 
              week + s(mean_latitude, mean_longitude, k = 50, bs = "ds") + 
              s(sum_duration_minutes, k = 5) + s(sum_effort_distance_km, k = 5) + 
              s(wt_mean_time_observations_started, k = 10) + s(median_number_observers, k = 4),
              weights = weights,
              data = ebird_ag_count$train, 
              method = "REML",
              gamma = 1.4, 
              family = nb)

  # let's look at the model
  summary(abd_gam)
  gam.check(abd_gam)

  plot(abd_gam, pages = 1)


  # to estimate the year trend we'll predict counts for a cell with 
  # average landcover across each year. 
  year_seq <- seq(min(ebird_ag_split$train$year), max(ebird_ag_split$train$year))
  nd <- data.frame(year = year_seq, id = 1) %>%
        left_join(mean_covariates, by = "id")
  pred_abd_gam <- predict(abd_gam, newdata=nd, se.fit = TRUE)
  pred_abd <- data.frame(year = nd$year, est = exp(pred_abd_gam$fit),
    lcl = exp(pred_abd_gam$fit - 1.96*pred_abd_gam$se.fit),
    ucl = exp(pred_abd_gam$fit + 1.96*pred_abd_gam$se.fit))

  # plot out the raw reporting rates and the fitted gam over time
  plot(0, 0, 
    xlab = "year", ylab = "Expected count", 
    xlim = c(min(pred_abd$year), max(pred_abd$year)), ylim = c(0, max(pred_abd$ucl)*1.1), 
    col = "white")
  polygon(x = c(pred_abd$year, rev(pred_abd$year), pred_abd$year[1]), 
          y = c(pred_abd$lcl, rev(pred_abd$ucl), pred_abd$lcl[1]), 
          border = alpha("white", 0), col = "grey80")
  lines(pred_abd$year, pred_abd$est, pch = 16)


  # let's compare occurrence and abundance

  # plot out the raw reporting rates and the fitted gam over time
  par(mfrow=c(1,2))

  plot(0, 0, 
    xlab = "year", ylab = "Encounter rate", main = "Occurrence model",
    xlim = c(min(pred_occ$year), max(pred_occ$year)), ylim = c(0, max(pred_occ$ucl)*1.1), 
    col = "white")
  polygon(x = c(pred_occ$year, rev(pred_occ$year), pred_occ$year[1]), 
          y = c(pred_occ$lcl, rev(pred_occ$ucl), pred_occ$lcl[1]), 
          border = alpha("white", 0), col = "grey80")
  lines(pred_occ$year, pred_occ$est, pch = 16)

  # plot out the raw reporting rates and the fitted gam over time
  plot(0, 0, 
    xlab = "year", ylab = "Expected count", main = "Abundance model",
    xlim = c(min(pred_abd$year), max(pred_abd$year)), ylim = c(0, max(pred_abd$ucl)*1.1), 
    col = "white")
  polygon(x = c(pred_abd$year, rev(pred_abd$year), pred_abd$year[1]), 
          y = c(pred_abd$lcl, rev(pred_abd$ucl), pred_abd$lcl[1]), 
          border = alpha("white", 0), col = "grey80")
  lines(pred_abd$year, pred_abd$est, pch = 16)



  # ####################################################################
  # TWEEDIE WITH ZEROS
  # ####################################################################

  weights <- ifelse(ebird_ag_split$train$number_checklists>5, 5, ebird_ag_split$train$number_checklists)

  abd_gam2 <- gam(total_count ~ s(year, k=5, m=1, bs = "ds") + 
              week + s(mean_latitude, mean_longitude, k = 50, bs = "ds") + 
              s(sum_duration_minutes, k = 5) + s(sum_effort_distance_km, k = 5) + 
              s(wt_mean_time_observations_started, k = 10) + s(median_number_observers, k = 4),
              method = "REML",
              gamma = 1.4, 
              weights = weights,
              data = ebird_ag_split$train, family = tw)

  # let's look at the model
  summary(abd_gam2)
  gam.check(abd_gam2)

  plot(abd_gam2, pages = 1)


  # we'll predict counts for a cell with average landcover across each year. 
  pred_abd_gam2 <- predict(abd_gam2, newdata=nd, se.fit = TRUE)
  pred_abd2 <- data.frame(year = nd$year, est = exp(pred_abd_gam2$fit),
    lcl = exp(pred_abd_gam2$fit - 1.96*pred_abd_gam2$se.fit),
    ucl = exp(pred_abd_gam2$fit + 1.96*pred_abd_gam2$se.fit))

  # plot out the raw reporting rates and the fitted gam over time
  plot(0, 0, 
    xlab = "year", ylab = "Expected count", 
    xlim = c(min(pred_abd2$year), max(pred_abd2$year)), ylim = c(0, max(pred_abd2$ucl)*1.1), 
    col = "white")
  polygon(x = c(pred_abd2$year, rev(pred_abd2$year), pred_abd2$year[1]), 
          y = c(pred_abd2$lcl, rev(pred_abd2$ucl), pred_abd2$lcl[1]), 
          border = alpha("white", 0), col = "grey80")
  lines(pred_abd2$year, pred_abd2$est, pch = 16)


  # let's compare occurrence and abundance

  # plot out the raw reporting rates and the fitted gam over time
  par(mfrow=c(1,2))

  # plot out the raw reporting rates and the fitted gam over time
  plot(0, 0, 
    xlab = "year", ylab = "Expected count", main = "Truncated data NB",
    xlim = c(min(pred_abd$year), max(pred_abd$year)), ylim = c(0, max(pred_abd$ucl)*1.1), 
    col = "white")
  polygon(x = c(pred_abd$year, rev(pred_abd$year), pred_abd$year[1]), 
          y = c(pred_abd$lcl, rev(pred_abd$ucl), pred_abd$lcl[1]), 
          border = alpha("white", 0), col = "grey80")
  lines(pred_abd$year, pred_abd$est, pch = 16)

 # plot out the raw reporting rates and the fitted gam over time
  plot(0, 0, 
    xlab = "year", ylab = "Expected count", main = "All data Tweedie",
    xlim = c(min(pred_abd2$year), max(pred_abd2$year)), ylim = c(0, max(pred_abd2$ucl)*1.1), 
    col = "white")
  polygon(x = c(pred_abd2$year, rev(pred_abd2$year), pred_abd2$year[1]), 
          y = c(pred_abd2$lcl, rev(pred_abd2$ucl), pred_abd2$lcl[1]), 
          border = alpha("white", 0), col = "grey80")
  lines(pred_abd2$year, pred_abd2$est, pch = 16)