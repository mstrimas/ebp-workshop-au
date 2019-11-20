library(tidyverse)
set.seed(1)

# generate a 10 x 10 grid
psi <- 0.3
p <- 0.5
p_surveyed <- 0.25
grid <- crossing(x = seq(0.5, 9.5, 1),
                 y = seq(0.5, 9.5, 1)) %>% 
  mutate(occupied = FALSE,
         surveyed = FALSE)
grid$occupied[sample.int(nrow(grid), nrow(grid) * psi)] <- TRUE
grid$surveyed[sample.int(nrow(grid), nrow(grid) * p_surveyed)] <- TRUE

col_pres <- "#4daf4a"
col_abs <- "grey80"

g_grid <- ggplot(grid) +
  aes(x = x, y = y) +
  geom_tile(aes(fill = occupied), color = "#000000") +
  scale_fill_manual("Site occupied",
                    values = c(col_abs, col_pres),
                    labels = c("No", "Yes")) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_equal() +
  theme_void() +
  theme(legend.text = element_text(size = 18),
        legend.title = element_text(size = 22),
        legend.position = "bottom")
ggsave("assets/img/12_occupancy_grid.svg", g_grid, width = 7, height = 7)

# surveys
survey <- grid %>% 
  filter(surveyed) %>% 
  crossing(visit = paste("Visit", 1:3)) %>% 
  mutate(detected = occupied & (runif(nrow(.)) <= p),
         missed = occupied & !detected)
# single survey
g_s1 <- ggplot(grid) +
  aes(x = x, y = y) +
  geom_tile(aes(fill = occupied), color = "#000000", show.legend = FALSE) +
  geom_point(data = survey %>% filter(visit == "Visit 1"), 
             aes(shape = detected), size = 5, stroke = 2) +
  # geom_point(data = survey %>% filter(visit == "Visit 1", missed), 
  #            shape = 1, size = 5, stroke = 2) +
  scale_fill_manual("Site occupied",
                    values = c(col_abs, col_pres),
                    labels = c("No", "Yes")) +
  scale_shape_manual("Species detected",
                     values = c(1, 19),
                     labels = c("No", "Yes")) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_equal() +
  theme_void() +
  theme(legend.text = element_text(size = 18),
        legend.title = element_text(size = 22),
        legend.position = "bottom")
ggsave("assets/img/12_occupancy_1survey.svg", g_s1, width = 7, height = 7)

# single survey
g_s3 <- ggplot(grid) +
  aes(x = x, y = y) +
  geom_tile(aes(fill = occupied), color = "#000000", show.legend = FALSE) +
  geom_point(data = survey, 
             aes(shape = detected), size = 5, stroke = 2) +
  scale_fill_manual("Site occupied",
                    values = c(col_abs, col_pres),
                    labels = c("No", "Yes")) +
  scale_shape_manual("Species detected",
                     values = c(1, 19),
                     labels = c("No", "Yes")) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_equal() +
  facet_wrap(~ visit, nrow = 1) +
  theme_void() +
  theme(strip.text = element_text(size = 22),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 22),
        legend.position = "bottom")
ggsave("assets/img/12_occupancy_3surveys.svg", g_s3, width = 18, height = 6)


# observations only
g_o1 <- ggplot(grid) +
  aes(x = x, y = y) +
  geom_tile(aes(fill = occupied), color = "#000000", show.legend = FALSE) +
  geom_point(data = survey %>% filter(visit == "Visit 1"), 
             aes(shape = detected), size = 5, stroke = 2) +
  scale_fill_manual("Site occupied",
                    values = c("white", "white"),
                    labels = c("No", "Yes")) +
  scale_shape_manual("Species detected",
                     values = c(1, 19),
                     labels = c("No", "Yes")) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_equal() +
  theme_void() +
  theme(legend.text = element_text(size = 18),
        legend.title = element_text(size = 22),
        legend.position = "bottom")
ggsave("assets/img/12_occupancy_obsonly.svg", g_o1, width = 7, height = 7)


# obs only 3 surveys
g_o3 <- ggplot(grid) +
  aes(x = x, y = y) +
  geom_tile(aes(fill = occupied), color = "#000000", show.legend = FALSE) +
  geom_point(data = survey, 
             aes(shape = detected), size = 5, stroke = 2) +
  scale_fill_manual("Site occupied",
                    values = c("white", "white"),
                    labels = c("No", "Yes")) +
  scale_shape_manual("Species detected",
                     values = c(1, 19),
                     labels = c("No", "Yes")) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_equal() +
  facet_wrap(~ visit, nrow = 1) +
  theme_void() +
  theme(strip.text = element_text(size = 22),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 22),
        legend.position = "bottom")
ggsave("assets/img/12_occupancy_3surveys_obsonly.svg", g_o3, width = 18, height = 6)



