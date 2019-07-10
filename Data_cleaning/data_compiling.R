# manually set the datpath using the data_pathway file 

# then source the raw data
source("Data_cleaning/dataformatting_2001.R")
source("Data_cleaning/dataformatting_2002.R")
source("Data_cleaning/dataformatting_2003.R")
source("Data_cleaning/dataformatting_2004.R")
source("Data_cleaning/dataformatting_2005.R")
source("Data_cleaning/dataformatting_2006.R")
source("Data_cleaning/dataformatting_2007.R")
source("Data_cleaning/dataformatting_2008.R")
source("Data_cleaning/dataformatting_2009.R")
source("Data_cleaning/dataformatting_2010.R")
source("Data_cleaning/dataformatting_2011.R")
source("Data_cleaning/dataformatting_2012.R")
source("Data_cleaning/dataformatting_2013.R")
source("Data_cleaning/dataformatting_2014.R")
source("Data_cleaning/dataformatting_2015.R")
source("Data_cleaning/dataformatting_2016.R")
source("Data_cleaning/dataformatting_2017.R")
source("Data_cleaning/dataformatting_2018.R")


unique(key2001$quadrat)
unique(key2002$quadrat)

length(intersect(unique(key2001$quadrat),
          unique(key2005$quadrat)))

length(intersect(unique(key2001$quadrat),
                 unique(key2003$quadrat)))

length(outersect(unique(key2001$quadrat),
                 unique(key2003$quadrat)))

length(intersect(unique(key2001$quadrat),
                 unique(key2006$quadrat)))

length(unique(key2006$quadrat))
