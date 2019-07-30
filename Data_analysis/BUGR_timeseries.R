## Set your datpath!! (file in Data_cleaning)

library(tidyverse)
library(readr)
library(ggplot2)

dat <- read_csv(paste(datpath_clean, "/bugrdat.csv", sep=""))

SC <- read_csv(paste(datpath_clean, "/SpeciesCodes.csv", sep=""))


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


tog<-left_join(dat, SC)

functog<-tog%>%
  group_by(quadratNew, year, status, func)%>%
  summarize(sumcov=sum(cover))%>%
  filter(!is.na(status), !is.na(func))

ggplot(functog, aes(year, sumcov))+geom_boxplot()+facet_grid(status~func)

