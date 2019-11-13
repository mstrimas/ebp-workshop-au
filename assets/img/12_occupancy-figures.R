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
g_grid <- ggplot(grid) +
  aes(x = x, y = y) +
  geom_tile(aes(fill = occupied), color = "#000000") +
  scale_fill_manual("Site occupied",
                    values = c("#999999", "#4daf4a"),
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
  geom_point(data = survey %>% filter(visit == "Visit 1", !missed), 
             aes(shape = detected), size = 4) +
  geom_point(data = survey %>% filter(visit == "Visit 1", missed), 
             shape = 1, color = "#ffffff", size = 4) +
  scale_fill_manual("Site occupied",
                    values = c("#999999", "#4daf4a"),
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
  geom_point(data = survey %>% filter(!missed), 
             aes(shape = detected), size = 4) +
  geom_point(data = survey %>% filter(missed), 
             shape = 1, color = "#ffffff", size = 4) +
  scale_fill_manual("Site occupied",
                    values = c("#999999", "#4daf4a"),
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
