library(tidyverse)
library(readxl)
library(stringr)

## GOPHER, COWPIE, BARE, LITTER

dat <- read_excel(datpath, 
                  sheet = "Plants2003", skip = 3) 

datnames <- read_excel(datpath, 
                       sheet = "Plants2003")[0,]

names(dat) = names(datnames)
names(dat)[1:3] = c("quadrat", "quadrat1", "site")


dat2003env <- dat %>%
  select(quadrat, site, GOPHER, BARE, ROCK, LITTER, COWPIE)

dat2003 <- dat %>%
  select(-c(GOPHER, BARE, ROCK, LITTER, COWPIE, quadrat1)) %>%
  gather(species, cover, "X__4":"X__57") %>%
  mutate(dummyspp = substr(species, 1,3)) %>%
  filter(dummyspp != "X__") %>%
  select(-dummyspp) %>%
  mutate(site = substr(quadrat, 1,2)) %>%
  filter(site == "TH") %>%
  mutate(year = 2003)


key2003 <- dat2003 %>%
  select(quadrat, site, year) %>%
  unique()

rm(dat, datnames)
