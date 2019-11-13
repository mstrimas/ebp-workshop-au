library(raster)
library(sf)
library(dplyr)
library(RandomFields)
library(here)

set.seed(1)
r <- raster(nrow = 20, ncol = 20, 
            xmn = -1, xmx = 1, ymn = -1, ymx = 1, vals = 1) %>% 
  prioritizr::simulate_data(n = 1, 
                            model = RPbernoulli(RMexp(), threshold = 0.1))


extents <- list(extent(-1, 0, 0, 1),
                extent(-1, 0, -1, 0),
                extent(0, 1, 0, 1),
                extent(-0, 1, -1, 0))
png(here("assets", "img", "09_advanced_lc-tiles.png"), 
    width = 600, height = 600)
par(mar = c(1, 1, 1, 1), mfrow = c(2, 2))
for (i in seq_along(extents)) {
  r_crop <- crop(r, extents[[i]])
  plot(st_as_sfc(st_bbox(extent(r_crop))), col = NA, border = NA)
  plot(r_crop, col = c("forestgreen", "grey80"), legend = FALSE, axes = FALSE, 
       add = TRUE)
  #title(paste("tile ", i))
}
dev.off()


png(here("assets", "img", "09_advanced_lc-mosaic.png"), 
    width = 600, height = 600)
par(mar = c(0, 0, 0, 0), mfrow = c(1, 1))
plot(st_as_sfc(st_bbox(extent(r))), col = NA, border = NA)
plot(r, col = c("forestgreen", "grey80"), legend = FALSE, axes = FALSE, 
     add = TRUE)
dev.off()

r_agg <- aggregate(r, 5, fun = mean)
p <- rasterToPoints(r_agg, spatial = TRUE) %>% 
  st_as_sf() %>% 
  slice(10) %>% 
  st_geometry()
buff <- st_buffer(p, 0.25)

png(here("assets", "img", "09_advanced_lc-buffer.png"), 
    width = 800, height = 800)
par(mar = c(0, 0, 0, 0), mfrow = c(1, 1))
plot(st_as_sfc(st_bbox(extent(r))), col = NA, border = NA)
plot(r, col = c("forestgreen", "grey80"), legend = FALSE, axes = FALSE, 
     add = TRUE)
plot(p, pch = 19, cex = 3, col = "black", add = TRUE)
plot(buff, border = "black", lwd = 5, add = TRUE)
dev.off()


png(here("assets", "img", "09_advanced_lc-pland-1.png"), 
    width = 600, height = 600)
par(mar = c(0, 0, 0, 0), mfrow = c(1, 1))
plot(st_as_sfc(st_bbox(extent(r))), col = NA, border = NA)
plot(r_agg, col = RColorBrewer::brewer.pal(5, "Greens"), axes = FALSE, 
     legend = FALSE, add = TRUE)
dev.off()

png(here("assets", "img", "09_advanced_lc-pland-0.png"), 
    width = 600, height = 600)
par(mar = c(0, 0, 0, 0), mfrow = c(1, 1))
plot(st_as_sfc(st_bbox(extent(r))), col = NA, border = NA)
plot(1 - r_agg, col = RColorBrewer::brewer.pal(5, "Greys"), axes = FALSE, 
     legend = FALSE, add = TRUE)
dev.off()
