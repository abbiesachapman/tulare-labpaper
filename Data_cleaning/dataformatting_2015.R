library(tidyverse)
library(readxl)
library(stringr)

## GOPHER, COWPIE, BARE, LITTER
dat <- read_excel(datpath, 
                  sheet = "Plants2015", skip = 2) 

datnames <- read_excel(datpath, 
                       sheet = "Plants2015")[0,]

names(dat) = names(datnames)
names(dat)[4:6] = c("quadrat", "transect", "site")


dat2015env <- dat %>%
  select(quadrat, transect, site, GOPHER, BARE, ROCK, LITTER, COWPIE) %>%
  mutate(year = 2015)
 

dat2015 <- dat %>%
  select(-c(GOPHER, BARE, ROCK, LITTER, COWPIE, X__1, X__2, X__3)) %>%
  gather(species, cover, "X__7":"X__29") %>%
  mutate(dummyspp = substr(species, 1,3)) %>%
  filter(dummyspp != "X__") %>%
  select(-dummyspp) %>%
  #mutate(site = substr(quadrat, 1,2)) %>%
  filter(site =="TH" | site =="TH-MEC" | site == "TH-North" | site == "PGE") %>%
  mutate(year = 2015)


key2015 <- dat2015 %>%
  select(transect, quadrat, site, year) %>%
  unique()

rm(dat, datnames)
