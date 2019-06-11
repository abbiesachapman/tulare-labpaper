library(tidyverse)
library(readxl)
library(stringr)

## GOPHER, COWPIE, BARE, LITTER
dat <- read_excel(datpath, 
                  sheet = "Plants2018", skip = 2) 

datnames <- read_excel(datpath, 
                       sheet = "Plants2018")[0,]

names(dat) = names(datnames)
names(dat)[4:6] = c("quadrat", "transect", "site")


dat2018env <- dat %>%
  select(quadrat, transect, site, GOPHER, BARE, ROCK, LITTER, COWPIE) %>%
  mutate(year = 2018)
 

dat2018 <- dat %>%
  select(-c(GOPHER, BARE, ROCK, LITTER, COWPIE, X__1, X__2, X__3)) %>%
  gather(species, cover, "X__7":"X__24") %>%
  mutate(dummyspp = substr(species, 1,3)) %>%
  filter(dummyspp != "X__") %>%
  select(-dummyspp) %>%
  #mutate(site = substr(quadrat, 1,2)) %>%
  filter(site =="TH-MEC" | site == "TH-North") %>%
  mutate(year = 2018)


key2018 <- dat2018 %>%
  select(transect, quadrat, site, year) %>%
  unique()

rm(dat, datnames)
