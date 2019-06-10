library(tidyverse)
library(readxl)
library(stringr)

## GOPHER, COWPIE, BARE, LITTER

dat <- read_excel(datpath, 
                  sheet = "Plants2013", skip = 1) 

datnames <- read_excel(datpath, 
                       sheet = "Plants2013")[0,]

names(dat) = names(datnames)
names(dat)[1:4] = c("year", "quadrat","transect", "site")


dat2013env <- dat %>%
  select(quadrat, site, transect, year, GOPHER, BARE, ROCK, LITTER, COWPIE)

# note - is PGE for the metcalf energy site? 
dat2013 <- dat %>%
  select(-c(GOPHER, BARE, ROCK, LITTER, COWPIE)) %>%
  gather(species, cover, "X__5":"X__29") %>%
  mutate(dummyspp = substr(species, 1,3)) %>%
  filter(dummyspp != "X__") %>%
  select(-dummyspp) %>%
  # mutate(site = substr(quadrat, 1,2)) %>%
  filter(site == "TH"| site == "TH-BUG") 


key2012 <- dat2012 %>%
  select(quadrat, site, year, transect) %>%
  unique()

rm(dat, datnames)
