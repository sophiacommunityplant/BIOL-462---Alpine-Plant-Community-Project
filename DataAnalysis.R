#===============================================================================
#Title: Shifts in alpine plant communities due to changing climates
#       in British Columbia Parks
#Author: Sophia Johnson
#Date: March 2026
#Course: BIOL 462, University of Victoria

#Description
#   This script analyses long-term changes in alpine plant community
#   diversity and composition in relation to mean annual temperature (MAT)
#   and mean annual precipitation (MAP) across three BC provincial parks:
#   Mount Assiniboine, Big White, and Garibaldi.

# Data:
#   AlpinePlants_community_matrix.csv
#   Rows = park x site x year observations
#   Columns 1-12 = metadata (Park, Year, MAT, MAP, etc.)
#   Columns 13-141 = species percent cover values
#
# R version: 4.4.2
# Required packages: tidyverse, vegan, lme4, lmerTest, ggplot2, patchwork
#
# Output files:
#   study_area_map.svg          _ Figure 1: map of three study parks
#   alpha_diversity_4panel.svg  — Figure 2: four alpha diversity metrics
#   pca_park_year.svg           — Figure 3: PCA ordination by park and year
#   regression_MAT_MAP_combined.svg — Figure 4: Shannon ~ MAT and MAP
#   focal_species_all.svg       — Figure 5: focal species cover over time
#===============================================================================

#Load packages
rm(list = ls())
remotes::install_github("ropensci/rnaturalearthhires")

library(tidyverse)
library(vegan)
library(lme4)
library(ggplot2)
library(lmerTest)
library(patchwork)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggspatial)
library(rnaturalearthhires)

#Load and prepare data
dat <- read.csv("AlpinePlants_community_matrix.csv", check.names = FALSE)
names(dat) <- make.unique(names(dat))

#Define species columns (Excel M to FK = 13:141)
start_pos <- 13
end_pos   <- 141
species_cols <- start_pos:end_pos

#Convert species to numeric and replace NA with 0
sp <- dat[, species_cols] %>%
  mutate(across(everything(), ~ suppressWarnings(as.numeric(.x)))) %>%
  mutate(across(everything(), ~ replace_na(.x, 0)))

sp_mat <- as.matrix(sp)

#Convert grouping variables
dat$Park <- as.factor(dat$Park)
dat$Sample.Station.Label <- as.factor(dat$Sample.Station.Label)
dat$Year <- as.numeric(as.character(dat$Year))

#Calculate Shannon Diversity per plot
dat$Shannon <- diversity(sp_mat, index = "shannon")

#Park-level summary, MAT and MAP are measured at the park level, not plot level
park_year_summary <- dat %>%
  group_by(Park, Year) %>%
  summarise(
    mean_Shannon = mean(Shannon, na.rm = TRUE),
    se_Shannon   = sd(Shannon, na.rm = TRUE) / sqrt(n()),
    MAT          = mean(MAT, na.rm = TRUE),
    MAP          = mean(MAP, na.rm = TRUE),
    .groups      = "drop")

#Combine white mountain heather columns
dat$CASSMER_combined <- dat$CASSMER + dat$CASSMER1

#===============================================================================
#Figue 1 - Map of Study Areas
#Get Canada and provinces
canada <- ne_states(country = "canada", returnclass = "sf")
bc     <- canada[canada$name == "British Columbia", ]

#Park coordinates (approximate centroids)
parks <- data.frame(
  name = c("Mt. Assiniboine\nProvincial Park",
           "Big White Mountain\nEcological Reserve",
           "Garibaldi\nProvincial Park"),
  lon  = c(-115.65, -118.93, -122.98),
  lat  = c(50.87,    49.72,   49.93),
  colour = c("#C1440E", "#4A7C9E", "#2D6A4F")
)
parks_sf <- st_as_sf(parks, coords = c("lon", "lat"), crs = 4326)

#Main map — zoomed into southern BC
main_map <- ggplot() +
  geom_sf(data = bc, fill = "grey92", colour = "grey50", linewidth = 0.5) +
  geom_sf(data = parks_sf, aes(colour = name), size = 5, shape = 16) +
  geom_sf_label(data = parks_sf, aes(label = name, colour = name),
                size = 2.8, nudge_y = 0.3, fontface = "bold",
                label.size = 0.2, fill = "white") +
  scale_colour_manual(values = c("#C1440E", "#4A7C9E", "#2D6A4F")) +
  coord_sf(xlim = c(-126, -113), ylim = c(48.5, 52.5)) +
  annotation_scale(location = "bl") +
  annotation_north_arrow(location = "tr",
                         style = north_arrow_fancy_orienteering(
                           fill = c("grey40", "white"))) +
  theme_classic(base_size = 11) +
  theme(legend.position = "none",
        axis.title = element_blank()) +
  labs()

#Inset — all of BC
inset_map <- ggplot() +
  geom_sf(data = canada, fill = "grey88", colour = "grey60", linewidth = 0.3) +
  geom_sf(data = bc, fill = "grey70", colour = "grey40", linewidth = 0.4) +
  geom_sf(data = parks_sf, aes(colour = name), size = 2, shape = 16) +
  scale_colour_manual(values = c("#C1440E", "#4A7C9E", "#2D6A4F")) +
  coord_sf(xlim = c(-140, -113), ylim = c(48, 60)) +
  theme_void() +
  theme(legend.position  = "none",
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))

#Combine with patchwork
final_map <- main_map +
  inset_element(inset_map, left = 0.6, bottom = 0.0, right = 1.1, top = 0.4)

print(final_map)
ggsave("study_area_map.svg", final_map, width = 8, height = 7, dpi = 300)

#===============================================================================
#Figure 2 - 4 panel alpha diversity

#Calculate all four metrics per plot:
#   Shannon, Simpson, Pielou's eveness, species richness
alpha_summary <- dat %>%
  mutate(
    richness = apply(sp_mat, 1, function(x) sum(x > 0, na.rm = TRUE)),
    shannon  = diversity(sp_mat, index = "shannon"),
    simpson  = diversity(sp_mat, index = "simpson"),
    evenness = if_else(richness > 1, shannon / log(richness), NA_real_)) %>%
  group_by(Park, Year) %>%
  summarise(
    richness_mean = mean(richness, na.rm = TRUE),
    shannon_mean  = mean(shannon,  na.rm = TRUE),
    simpson_mean  = mean(simpson,  na.rm = TRUE),
    evenness_mean = mean(evenness, na.rm = TRUE),
    richness_se   = sd(richness, na.rm = TRUE) / sqrt(n()),
    shannon_se    = sd(shannon,  na.rm = TRUE) / sqrt(n()),
    simpson_se    = sd(simpson,  na.rm = TRUE) / sqrt(n()),
    evenness_se   = sd(evenness, na.rm = TRUE) / sqrt(n()),
    .groups = "drop")

park_only_alpha_summary <- dat %>%
  mutate(
    richness = apply(sp_mat, 1, function(x) sum(x > 0, na.rm = TRUE)),
    shannon  = diversity(sp_mat, index = "shannon"),
    simpson  = diversity(sp_mat, index = "simpson"),
    evenness = if_else(richness > 1, shannon / log(richness), NA_real_)) %>%
  group_by(Park) %>%
  summarise(
    richness_mean = mean(richness, na.rm = TRUE),
    shannon_mean  = mean(shannon,  na.rm = TRUE),
    simpson_mean  = mean(simpson,  na.rm = TRUE),
    evenness_mean = mean(evenness, na.rm = TRUE),
    richness_se   = sd(richness, na.rm = TRUE) / sqrt(n()),
    shannon_se    = sd(shannon,  na.rm = TRUE) / sqrt(n()),
    simpson_se    = sd(simpson,  na.rm = TRUE) / sqrt(n()),
    evenness_se   = sd(evenness, na.rm = TRUE) / sqrt(n()),
    .groups = "drop")

#Reshape to long format for faceting
plot_df <- alpha_summary %>%
  pivot_longer(cols = c(richness_mean, shannon_mean, simpson_mean, 
                        evenness_mean),
               names_to = "metric", values_to = "value") %>%
  mutate(se = case_when(
      metric == "richness_mean" ~ richness_se,
      metric == "shannon_mean"  ~ shannon_se,
      metric == "simpson_mean"  ~ simpson_se,
      metric == "evenness_mean" ~ evenness_se),
    metric = recode(metric,
                    richness_mean = "Species Richness",
                    shannon_mean  = "Shannon (H')",
                    simpson_mean  = "Simpson",
                    evenness_mean = "Evenness (Pielou)"))
label_df_4panel <- data.frame(
  metric      = c("Evenness (Pielou)", "Shannon (H')", "Simpson", "Species Richness"),
  panel_label = c("(a)", "(b)", "(c)", "(d)"))

fig_4panel <- ggplot(plot_df, aes(x = Year, y = value, group = Park, 
                                  colour = Park)) +
  geom_errorbar(aes(ymin = value - se, ymax = value + se),
                width = 0.4, linewidth = 0.4, alpha = 0.5, na.rm = TRUE) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.5) +
  geom_text(data = label_df_4panel,
            aes(x = -Inf, y = Inf, label = panel_label),
            hjust    = 1.52,
            vjust    = -1.5,
            colour   = "black",
            fontface = "plain",
            size     = 3.8,
            inherit.aes = FALSE) +
  coord_cartesian(clip = "off") +
  scale_colour_manual(values = c("#C1440E", "#4A7C9E", "#2D6A4F")) +
  scale_x_continuous(breaks = sort(unique(dat$Year))) +
  facet_wrap(~ metric, scales = "free_y", ncol = 2) +
  theme_classic(base_size = 12) +
  theme(legend.position  = "bottom",
        legend.title     = element_blank(),
        axis.text.x      = element_text(angle = 45, hjust = 1),
        plot.margin = margin(t = 15, r = 5, b = 5, l = 5),
        strip.background = element_blank(),
        strip.text       = element_text(size = 11),
        panel.spacing    = unit(1.2, "lines")) +
  labs(x = "Year", y = "Mean ± SE")

print(fig_4panel)
ggsave("alpha_diversity_4panel.svg", fig_4panel, width = 10, height = 7,
       dpi = 300)

#===============================================================================
#Figure 3 - Ordination: PCA by Park and Year
#Hellinger-transformed cover data

#Remove any all-zero rows before transformation
nonzero_rows2 <- rowSums(sp_mat) > 0
sp_mat_ord    <- sp_mat[nonzero_rows2, ]
dat_ord       <- dat[nonzero_rows2, ]

#Hellinger-transform before PCA
sp_hel <- decostand(sp_mat_ord, method = "hellinger")
pca <- rda(sp_hel, scale = FALSE)
pca_summary  <- summary(pca)
pct_var      <- round(pca_summary$cont$importance[2, 1:2] * 100, 1)
cat("PC1 explains:", pct_var[1], "% of variance\n")
cat("PC2 explains:", pct_var[2], "% of variance\n")

#extract site scores and attach meta data
scores_df <- as.data.frame(scores(pca, display = "sites"))
scores_df$Park    <- dat_ord$Park
scores_df$Year    <- as.factor(dat_ord$Year)
scores_df$Shannon <- dat_ord$Shannon

#Calculate park-year centroids
centroids <- scores_df %>%
  group_by(Park, Year) %>%
  summarise(PC1 = mean(PC1),
            PC2 = mean(PC2),
            .groups = "drop")

# Axis labels with variance explained
xlab <- paste0("PC1 (", pct_var[1], "% variance explained)")
ylab <- paste0("PC2 (", pct_var[2], "% variance explained)")

fig_pca_park <- ggplot(scores_df, aes(x = PC1, y = PC2, colour = Park)) +
  geom_point(alpha = 0.25, size = 1.8) +
  geom_point(data = centroids, size = 5, stroke = 1.2) +
  stat_ellipse(aes(group = Park), level = 0.95, linewidth = 0.8)+
  geom_text(data = centroids, aes(label = Year),
            colour = "black", size = 2.8,
            vjust = -1, fontface = "bold",  nudge_y = 0.025, nudge_x = 0.01) +
  scale_colour_manual(values = c("#C1440E", "#4A7C9E", "#2D6A4F"))+
  scale_fill_manual(values   = c("#C1440E", "#4A7C9E", "#2D6A4F"))+
  theme_classic(base_size = 12) +
  theme(legend.position  = "right",
        legend.title     = element_blank()) +
  labs(    x = xlab, y = ylab)
print(fig_pca_park)
ggsave("pca_park_year.svg", fig_pca_park, width = 8, height = 6, dpi = 300)

#===============================================================================
#How has community composition changed with temperature?

#Mixed Effects Model
model_RQ1_alpha <- lmer(
  Shannon ~ MAT + Year + (1 | Park/Sample.Station.Label),
  data = dat)

summary(model_RQ1_alpha)
plot(model_RQ1_alpha)          # residuals vs fitted
qqnorm(resid(model_RQ1_alpha)) # Q-Q plot
qqline(resid(model_RQ1_alpha))

#PERMANOVA: tests whether community composition differs with MAT, Year, Park
# 999 permutations, Bray-Curtis dissimilarity
set.seed(42)
adonis_RQ1 <- adonis2(
  sp_mat ~ MAT + Year + Park,
  data = dat,
  method = "bray")

adonis_RQ1
bray_dist <- vegdist(sp_mat, method = "bray")
disp <- betadisper(bray_dist, dat$Park)
permutest(disp, permutations = 999)

#Useful for visualizing the overall MAT-diversity relationship


#Linear Regression: MAT Regression using park-year means
lm_park <- lm(mean_Shannon ~ MAT, data = park_year_summary)
summary(lm_park)

#Figure 4a
fig_reg_park <- ggplot(park_year_summary,
                       aes(x = MAT, y = mean_Shannon,
                           colour = Park, shape = Park)) +
  geom_smooth(method = "lm", se = TRUE,
              colour = "grey40", fill = "grey85", linewidth = 0.8,
              data = park_year_summary,
              aes(x = MAT, y = mean_Shannon),
              inherit.aes = FALSE) +
  geom_errorbar(aes(ymin = mean_Shannon - se_Shannon,
                    ymax = mean_Shannon + se_Shannon),
                width = 0.05, alpha = 0.6) +
  geom_point(size = 4) +
  geom_text(aes(label = Year), size = 2.8, vjust = -1.2,
            show.legend = FALSE) +
  scale_colour_manual(values = c("#C1440E", "#4A7C9E", "#2D6A4F"))+
scale_fill_manual(values   = c("#C1440E", "#4A7C9E", "#2D6A4F"))+
  scale_shape_manual(values = c(16, 17, 15)) +
  theme_classic() +
  theme(legend.position = "none")+
  labs(subtitle = paste0("R² = ", round(summary(lm_park)$r.squared, 3),
                      ",  p = ", round(summary(lm_park)$coefficients[2,4], 4)),
    x = "Mean Annual Temperature (°C)",
    y = "Mean Shannon Diversity (H') ± SE")
print(fig_reg_park)

#===============================================================================
#How has community composition changed with precipitation?
#How have temperature and precipitation interacted?

#Mixed effects model with MAT x MAP interactions
model_RQ2 <- lmer(
  Shannon ~ MAT * MAP + Year + (1 | Park/Sample.Station.Label),
  data = dat)

summary(model_RQ2)
plot(model_RQ2)          # residuals vs fitted
qqnorm(resid(model_RQ2)) # Q-Q plot
qqline(resid(model_RQ2))

#PERMANOVA with interaction term
#Remove rows with NA in predictors
complete_rows <- complete.cases(dat[, c("MAT", "MAP", "Year")]) &
  rowSums(sp_mat) > 0   # also remove empty species rows

sp_mat_clean <- sp_mat[complete_rows, ]
dat_clean    <- dat[complete_rows, ]

# Run PERMANOVA
adonis_RQ2 <- adonis2(
  sp_mat_clean ~ MAT * MAP + Year,
  data = dat_clean,
  method = "bray")

adonis_RQ2

#Regression using park-year means
lm_MAP_park <- lm(mean_Shannon ~ MAP, data = park_year_summary)
summary(lm_MAP_park)
#Figure 4b
fig_MAP_park <- ggplot(park_year_summary,
                       aes(x = MAP, y = mean_Shannon,
                           colour = Park, shape = Park)) +
  geom_smooth(method = "lm", se = TRUE,
              colour = "grey40", fill = "grey85", linewidth = 0.8,
              data = park_year_summary,
              aes(x = MAP, y = mean_Shannon),
              inherit.aes = FALSE) +
  geom_errorbar(aes(ymin = mean_Shannon - se_Shannon,
                    ymax = mean_Shannon + se_Shannon),
                width = 10, alpha = 0.6) +
  geom_point(size = 4) +
  geom_text(aes(label = Year), size = 2.8, vjust = -1.2,
            show.legend = FALSE) +
  scale_colour_manual(values = c("#C1440E", "#4A7C9E", "#2D6A4F"))+
  scale_fill_manual(values   = c("#C1440E", "#4A7C9E", "#2D6A4F"))+
  scale_shape_manual(values = c(16, 17, 15)) +
  theme_classic()+
  labs(subtitle = paste0("R² = ", round(summary(lm_MAP_park)$r.squared, 3),
                      ",  p = ", round(summary(lm_MAP_park)$coefficients[2,4], 4)),
    x = "Mean Annual Precipitation (mm)",
    y = "Mean Shannon Diversity (H') ± SE")
print(fig_MAP_park)

#Figure 4: combined MAP and MAT regressions
fig_reg_combined <- fig_reg_park + fig_MAP_park+
  plot_annotation(tag_levels = "a",
                  tag_prefix  = "(",
                  tag_suffix  = ")")
print(fig_reg_combined)
ggsave("regression_MAT_MAP_combined.svg", fig_reg_combined,
       width = 14, height = 5, dpi = 300)

#===============================================================================
#How have focal species changed in abundance over time?

summary_PHYLEMP <- dat %>%
  group_by(Park, Year) %>%
  summarise(mean_cover = mean(PHYLEMP, na.rm = TRUE),
            se_cover   = sd(PHYLEMP, na.rm = TRUE) / sqrt(n()),
            MAT        = mean(MAT, na.rm = TRUE),
            MAP        = mean(MAP, na.rm = TRUE),
            .groups    = "drop")

summary_VACCSCO <- dat %>%
  group_by(Park, Year) %>%
  summarise(mean_cover = mean(VACCSCO, na.rm = TRUE),
            se_cover   = sd(VACCSCO, na.rm = TRUE) / sqrt(n()),
            MAT        = mean(MAT, na.rm = TRUE),
            MAP        = mean(MAP, na.rm = TRUE),
            .groups    = "drop")

summary_CASSMER <- dat %>%
  group_by(Park, Year) %>%
  summarise(mean_cover = mean(CASSMER_combined, na.rm = TRUE),
            se_cover   = sd(CASSMER_combined, na.rm = TRUE) / sqrt(n()),
            MAT        = mean(MAT, na.rm = TRUE),
            MAP        = mean(MAP, na.rm = TRUE),
            .groups    = "drop")

#Mixed effects models

#Pink Mountain Heather
model_heather_pink <- lmer(
  PHYLEMP ~ MAT + MAP + Year + (1 | Park/Sample.Station.Label),
  data = dat)
cat("=== Pink Mountain Heather (PHYLEMP) ===\n")
print(summary(model_heather_pink))

#Grouseberry
model_grouseberry <- lmer(
  VACCSCO ~ MAT + MAP + Year + (1 | Park/Sample.Station.Label),
  data = dat)
cat("\n=== Grouseberry (VACCSCO) ===\n")
print(summary(model_grouseberry))

# White Mountain Heather
model_heather_white <- lmer(
  CASSMER_combined ~ MAT + MAP + Year + (1 | Park/Sample.Station.Label),
  data = dat)
cat("\n=== White Mountain Heather (CASSMER combined) ===\n")
print(summary(model_heather_white))


#Figure 5: All three species on one panel

all_focal <- bind_rows(
  summary_PHYLEMP %>% mutate(Species = "Pink Mtn Heather\n(Phyllodoce empetriformis)"),
  summary_VACCSCO %>% mutate(Species = "Grouseberry\n(Vaccinium scoparium)"),
  summary_CASSMER %>% mutate(Species = "White Mtn Heather\n(Cassiope mertensiana)")) %>%
  mutate(Species = factor(Species, levels = c(
    "Pink Mtn Heather\n(Phyllodoce empetriformis)",
    "Grouseberry\n(Vaccinium scoparium)",
    "White Mtn Heather\n(Cassiope mertensiana)" )))
label_df <- data.frame(
  Species = factor(
    c("Pink Mtn Heather\n(Phyllodoce empetriformis)",
      "Grouseberry\n(Vaccinium scoparium)",
      "White Mtn Heather\n(Cassiope mertensiana)"),
    levels = levels(all_focal$Species)),
  label = c("(a)", "(b)", "(c)"))

fig_all_focal <- ggplot(all_focal,
                        aes(x = Year, y = mean_cover,
                            colour = Park, group = Park)) +
  geom_ribbon(aes(ymin = pmax(mean_cover - se_cover, 0),
                  ymax = mean_cover + se_cover,
                  fill = Park),
              alpha = 0.15, colour = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  geom_text(data = label_df,
            aes(x = -Inf, y = Inf, label = label),
            hjust    = 1.52,
            vjust    = -1.5,
            colour   = "black",
            size     = 3.87,
            inherit.aes = FALSE) +
  coord_cartesian(clip = "off") +    
   scale_colour_manual(values = c("#C1440E", "#4A7C9E", "#2D6A4F"))+
  scale_fill_manual(values   = c("#C1440E", "#4A7C9E", "#2D6A4F"))+
  scale_x_continuous(breaks = sort(unique(dat$Year))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) + 
  facet_wrap(~ Species, scales = "free_y", ncol = 1) +
  theme_classic(base_size = 12) +
  theme(legend.position  = "bottom",
        legend.title     = element_blank(),
        axis.text.x      = element_text(angle = 45, hjust = 1),
        strip.background = element_blank(),
        strip.text       = element_text(size = 11, family = "sans"),
        panel.spacing    = unit(1.2, "lines"),    
        strip.placement  = "outside",
        plot.margin      = margin(t = 15, r = 5, b = 5, l = 5)) +
  labs(       x = "Year", y = "Mean % Cover")

print(fig_all_focal)
ggsave("focal_species_all.svg", fig_all_focal,
       width = 8, height = 11, dpi = 300)

#===============================================================================
#Package citations
citation()           
citation("vegan")
citation("lme4")
citation("lmerTest")
citation("ggplot2")
citation("tidyverse")
citation("patchwork")
