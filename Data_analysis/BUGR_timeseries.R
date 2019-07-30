## Set your datpath!! (file in Data_cleaning)

library(tidyverse)
library(readr)

dat <- read_csv(paste(datpath_clean, "/bugrdat.csv", sep=""))

