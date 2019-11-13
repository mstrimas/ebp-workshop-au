library(auk)
library(here)
library(tidyverse)

ebd_dir <- here("raw-data", "ebd")
data_dir <- here("raw-data", "data")
dir.create(ebd_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

# countries
f_ebd <- here("raw-data", "ebd_relSep-2019_yucatan.txt")
f_sed <- here("raw-data", "ebd_sampling_relSep-2019_yucatan.txt")
if (!file.exists(f_ebd)) {
  auk_ebd("ebd_relApr-2019.txt", "ebd_sampling_relApr-2019.txt") %>% 
    auk_country(c("Mexico", "Guatemala", "Belize")) %>% 
    auk_date(c("2014-01-01", "2015-12-31")) %>% 
    auk_filter(f_ebd, f_sed)
}

# mexican states in yucatan
f_ebd_workshop <- here("raw-data", "ebd_2014-2015_yucatan.txt")
f_sed_workshop <- here("raw-data", "ebd_sampling_2014-2015_yucatan.txt")
auk_ebd(f_ebd, f_sed) %>% 
  auk_state(c("MX-CAM", "MX-CHP", "MX-TAB", "MX-ROO", "MX-YUC")) %>% 
  auk_filter(f_ebd_workshop, f_sed_workshop, overwrite = TRUE)
# other countries
f_ebd_nomx <- here("raw-data", "ebd_2014-2015_nomx.txt")
f_sed_nomx <- here("raw-data", "ebd_sampling_2014-2015_nomx.txt")
auk_ebd(f_ebd, f_sed) %>% 
  auk_country(c("Guatemala", "Belize")) %>% 
  auk_filter(f_ebd_nomx, f_sed_nomx, overwrite = TRUE)

# trim headers and combine
str_glue("sed '1d' {f_ebd_nomx} >> {f_ebd_workshop}") %>% system()
str_glue("sed '1d' {f_sed_nomx} >> {f_sed_workshop}") %>% system()
unlink(f_ebd_nomx)
unlink(f_sed_nomx)

# create download file
# part I
f_ebd_workshop %>% file.copy(., file.path(ebd_dir, basename(.)))
f_sed_workshop %>% file.copy(., file.path(ebd_dir, basename(.)))

# part II files
book_files <- c("ebd_woothr_june_bcr27_zf.csv",
                "gis-data.gpkg",
                "mcd12q1_classes.csv",
                "pland-elev_location-year.csv",
                "pland-elev_prediction-surface.csv",
                "prediction-surface.tif")
file.copy(file.path("~/projects/ebird-best-practices/data/", book_files),
          file.path(data_dir, book_files))

# zip
f_zip <- "ebp-workshop-data.zip"
unlink(f_zip)
str_glue("cd {here('raw-data')}; zip -r {f_zip} data/ ebd/") %>% 
  system()
unlink(here("raw-data", "data"), recursive = TRUE)
unlink(here("raw-data", "ebd"), recursive = TRUE)
unlink(f_ebd_workshop)
unlink(f_sed_workshop)
