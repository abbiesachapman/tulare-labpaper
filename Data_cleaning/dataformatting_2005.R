library(tidyverse)
library(readxl)
library(stringr)

## GOPHER, COWPIE, BARE, LITTER
dat <- read_excel(datpath, 
                  sheet = "Plants2005", skip = 4) 

datnames <- read_excel(datpath, 
                       sheet = "Plants2005")[0,]

names(dat) = names(datnames)
names(dat)[4:6] = c("quadrat", "transect", "site")


dat2005env <- dat %>%
  select(quadrat, transect, site, GOPHER, BARE, ROCK, LITTER, COWPIE) %>%
  mutate(year = 2005)

dat2005 <- dat %>%
  select(-c(GOPHER, BARE, ROCK, LITTER, COWPIE, X__1, X__2, X__3)) %>%
  gather(species, cover, "X__7":"X__64") %>%
  mutate(dummyspp = substr(species, 1,3)) %>%
  filter(dummyspp != "X__") %>%
  select(-dummyspp) %>%
  #mutate(site = substr(quadrat, 1,2)) %>%
  filter(site == "TH") %>%
  mutate(year = 2005)


key2005 <- dat2005 %>%
  select(quadrat, site, year) %>%
  unique()

rm(dat, dat2, datnames)
