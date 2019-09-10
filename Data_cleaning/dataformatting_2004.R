library(tidyverse)
library(readxl)
library(stringr)

## GOPHER, COWPIE, BARE, LITTER

dat <- read_excel(datpath, 
                  sheet = "Plants2004", skip = 3) 

datnames <- read_excel(datpath, 
                       sheet = "Plants2004")[0,]

names(dat) = names(datnames)
names(dat)[1:3] = c("quadrat", "quadrat1", "site")


dat2004env <- dat %>%
  select(quadrat, site, GOPHER, BARE, ROCK, LITTER, COWPIE) %>%
  mutate(year = 2004)

dat2004 <- dat %>%
  select(-c(GOPHER, BARE, ROCK, LITTER, COWPIE, quadrat1)) %>%
  gather(species, cover, "X__4":"X__57") %>%
  mutate(dummyspp = substr(species, 1,3)) %>%
  filter(dummyspp != "X__") %>%
  select(-dummyspp) %>%
  # mutate(site = substr(quadrat, 1,2)) %>%
  filter(site == "TH") %>%
  mutate(year = 2004)


key2004 <- dat2004 %>%
  select(quadrat, site, year) %>%
  unique()

rm(dat, datnames)
