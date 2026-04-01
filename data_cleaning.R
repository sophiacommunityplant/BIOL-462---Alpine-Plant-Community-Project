## data cleaning for community project
setwd("C:/Users/sophi/OneDrive - University of Victoria/UVIC/UVIC Year 4 Semester 2 2025/BIOL 462/Data")
library(dplyr)
library(tidyr)

AlpinePlants <- read.csv("Alpine Plants.csv") 
BAlpinePlants <- AlpinePlants %>%
  filter(!is.na(Species), Species != "NULL") %>%
  select(Park, Sample.Station.Label, Year,
         Species, Percent.Cover)
CAlpinePlants <- BAlpinePlants %>%
  mutate(site_id = paste(Park, Sample.Station.Label, Year, sep = "_"))

CAlpinePlants <- CAlpinePlants %>%
  mutate(Species = toupper(Species))

summarised_alpine_plants <- CAlpinePlants %>%
  group_by(site_id, Species, Park, Year, Sample.Station.Label) %>%
  summarise(Percent.Cover = sum(Percent.Cover), .groups = "drop")

AlpinePlantsdone <- summarised_alpine_plants %>%
  pivot_wider(id_cols = c(site_id, Park, Year, Sample.Station.Label), names_from = Species, values_from = Percent.Cover, values_fill = 0)

write.csv(AlpinePlantsdone, file = "AlpinePlants_community_matrix.csv", row.names = FALSE)
  