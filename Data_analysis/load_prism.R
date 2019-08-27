library(tidyverse) # a suite of packages for wrangling and tidying data
library(prism)     # package to access and download climate data
library(raster)    # the climate data comes in raster files- this package helps process those
library(stringr)   # character manipulation
library(magrittr)

#set a file path where prism data will be stored
options(prism.path = 'C:\\Users\\Lina\\Dropbox\\Academics\\PRISM\\PRISM_data.path')

#select the type of data and data range
get_prism_monthlys(type = "ppt", years = 2000:2018, mo = 1:12, keepZip = F)
get_prism_monthlys(type = "tmean", years = 2000:2018, mo = 1:12, keepZip = F)
get_prism_monthlys(type = "tmax", years = 2000:2018, mo = 1:12, keepZip = F)
get_prism_monthlys(type = "tmin", years = 2000:2018, mo = 1:12, keepZip = F)

#grab the prism data and compile the files
climate_data <- ls_prism_data() %>%  
  prism_stack(.) 

#extract project coordinates from raster stack
climate_crs <- climate_data@crs@projargs

#set coordiante points
mypoints <- data.frame(id = c(1), lat = c(37.220838), long = c(-121.754843))

#convert points to spatial points data frame
coordinates(mypoints) <- c('long', 'lat')
proj4string(mypoints) <- CRS(climate_crs)

#extract data from raster
climate_extract <- data.frame(coordinates(mypoints), mypoints$id, extract(climate_data, mypoints))

#reshape data - gather dates in one column
climate_extract <- climate_extract %>%  gather(date, value, 4:ncol(climate_extract))

#split date into type, year, month
climate_extract$date <- gsub("PRISM_", "", climate_extract$date) %>% 
                        gsub("_stable_4kmM3_", "", .) %>% 
                        gsub("_stable_4kmM2_", "", .) %>%
                        gsub("ppt", "prcp", .) %>%
                        gsub("tmean", "tavg", .) %>%
                        gsub("_bil", "", .)
climate_extract <- separate(climate_extract, "date", c("type", "Year", "Month"), sep = c(4, 8, 10))
climate_extract <- climate_extract[, -7]
  
#reshape data again - ppt, tmax, tmin have their own columns
climate_extract <- climate_extract %>%  spread(type, value)

#order data by coordinates
climate_extract <- climate_extract[order(climate_extract$mypoints.id),]

#export data as csv
write.csv(climate_extract, file = "Tulare_prism_data.csv")
