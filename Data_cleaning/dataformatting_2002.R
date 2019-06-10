library(tidyverse)
library(readxl)
library(stringr)

## GOPHER, COWPIE, BARE, LITTER
datpath <- "C:/Users/eliza/University Of Oregon/O365.hallett-lab - Documents/Tulare-LabPaper/Serpentine_Plants2001_2018_MK_working.xlsx"

dat <- read_excel(datpath, 
                  sheet = "Plants 2002", skip = 3)

datnames <- read_excel(datpath, 
                       sheet = "Plants 2002")[0,]

names(dat) = names(datnames)
names(dat)[1:2] = c("quadrat", "site")


dat2002env <- dat %>%
  select(quadrat, site, GOPHER, BARE, ROCK, LITTER, COWPIE)

dat2002 <- dat %>%
  select(-c(GOPHER, BARE, ROCK, LITTER, COWPIE)) %>%
  gather(species, cover, "X__3":"X__69") %>%
  mutate(dummyspp = substr(species, 1,3)) %>%
  filter(dummyspp != "X__") %>%
  select(-dummyspp) %>%
  #mutate(site = substr(quadrat, 1,2)) %>%
  filter(site == "TH") %>%
  mutate(year = 2002)


key2002 <- dat2002 %>%
  select(quadrat, site, year) %>%
  unique()

rm(dat, datnames)
