library(tidyverse)
library(readxl)
library(stringr)

## GOPHER, COWPIE, BARE, LITTER

dat <- read_excel(datpath, 
                  sheet = "Plants 2006", skip = 1) 

datnames <- read_excel(datpath, 
                       sheet = "Plants 2006")[0,]

names(dat) = names(datnames)
names(dat)[1:6] = c("year", "Page", "OddEven", "quadrat","transect", "site")


dat2006env <- dat %>%
  select(quadrat, site, transect, year, GOPHER, BARE, ROCK, LITTER, COWPIE)

dat2006 <- dat %>%
  select(-c(GOPHER, BARE, ROCK, LITTER, COWPIE, Page, OddEven)) %>%
  gather(species, cover, "X__7":"X__77") %>%
  mutate(dummyspp = substr(species, 1,3)) %>%
  filter(dummyspp != "X__") %>%
  select(-dummyspp) %>%
  # mutate(site = substr(quadrat, 1,2)) %>%
  filter(site == "TH") %>%
  mutate(year = 2006)


key2006 <- dat2006 %>%
  select(quadrat, site, year, transect) %>%
  unique()

rm(dat, datnames)
