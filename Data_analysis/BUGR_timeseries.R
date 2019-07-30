## Set your datpath!! (file in Data_cleaning)

library(tidyverse)
library(readr)
library(ggplot2)

dat <- read_csv(paste(datpath_clean, "/bugrdat.csv", sep=""))


View(dat)

richness <- dat %>%
  filter(cover != 0, spname != c("Unknown", "Moss")) %>%
  group_by(year, quadratNew) %>%
  summarize(richness = length(unique(spname))) %>%
  group_by(year) %>%
  summarize(mean_rich = mean(richness))

ggplot(richness, aes(year, mean_rich)) +
  geom_point() +
  geom_line()
