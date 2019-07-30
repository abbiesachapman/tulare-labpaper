library(tidyverse)
library(readxl)
library(stringr)

## GOPHER, COWPIE, BARE, LITTER

dat <- read_excel(datpath, 
                  sheet = "Plants2001", skip = 3)

datnames <- read_excel(datpath, 
                       sheet = "Plants2001")[0,]

names(dat) = names(datnames)
names(dat)[1:2] = c("quadrat", "site")


dat2001env <- dat %>%
  select(quadrat, site, GOPHER, BARE, ROCK, LITTER, COWPIE)

dat2001 <- dat %>%
  select(-c(GOPHER, BARE, ROCK, LITTER, COWPIE)) %>%
  gather(species, cover, "X__3":"X__47") %>%
  mutate(dummyspp = substr(species, 1,3)) %>%
  filter(dummyspp != "X__") %>%
  select(-dummyspp) %>%
  mutate(site = substr(quadrat, 1,2)) %>%
  filter(site == "TH") %>%
  mutate(year = 2001)


key2001 <- dat2001 %>%
  select(quadrat, site, year) %>%
  unique()

rm(dat, datnames)

