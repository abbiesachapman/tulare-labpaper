library(tidyverse)
library(readxl)
library(stringr)

## GOPHER, COWPIE, BARE, LITTER

dat <- read_excel(datpath, 
                  sheet = "Plants2009", skip = 1) 

datnames <- read_excel(datpath, 
                       sheet = "Plants2009")[0,]

names(dat) = names(datnames)
names(dat)[1:6] = c("year", "quadrat","transect", "site", "subsite", "numberedsite")


dat2009env <- dat %>%
  select(quadrat, site, transect, year, GOPHER, BARE, ROCK, LITTER, COWPIE)

# note - is PGE for the metcalf energy site? 
dat2009 <- dat %>%
  select(-c(GOPHER, BARE, ROCK, LITTER, COWPIE, subsite, numberedsite)) %>%
  gather(species, cover, "X__7":"X__28") %>%
  mutate(dummyspp = substr(species, 1,3)) %>%
  filter(dummyspp != "X__") %>%
  select(-dummyspp) %>%
  # mutate(site = substr(quadrat, 1,2)) %>%
  filter(site == "TH"| site == "TH-BUG") 


key2009 <- dat2009 %>%
  select(quadrat, site, year, transect) %>%
  unique()

rm(dat, datnames)
