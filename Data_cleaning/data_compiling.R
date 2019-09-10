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

## game plan with this: 
# append all the years of data
# add the thermal column to the quadrat master
# merge in the quadrat master - this will give the treatments and the unique plot code
# select only the burned grazed plots 
# then link to the species key!
## one way to check that you have all the species, is to make a data frame of all the unique species in alldat
## do a full join with the species key code - and look for gaps

# remove column "transect" from 2005-2018 to use rbind
dat_2005 <- dat2005 %>%
  select(-transect)

dat_2006 <- dat2006 %>%
  select(-transect)

dat_2007 <- dat2007 %>%
  select(-transect)

dat_2008 <- dat2008 %>%
  select(-transect)

dat_2009 <- dat2009 %>%
  select(-transect)

dat_2010 <- dat2010 %>%
  select(-transect)

dat_2011 <- dat2011 %>%
  select(-transect)

dat_2012 <- dat2012 %>%
  select(-transect)

dat_2013 <- dat2013 %>%
  select(-transect)

dat_2014 <- dat2014 %>%
  select(-transect)

dat_2015 <- dat2015 %>%
  select(-transect)

dat_2016 <- dat2016 %>%
  select(-transect)

dat_2017 <- dat2017 %>%
  select(-transect)

dat_2018 <- dat2018 %>%
  select(-transect)

# append all the data for species
alldat <- rbind(dat2001, dat2002, dat2003, dat2004, dat_2005, dat_2006, 
                dat_2007, dat_2008, dat_2009, dat_2010, dat_2011, 
                dat_2012, dat_2013, dat_2014, dat_2015, dat_2016, 
                dat_2017, dat_2018)
names(alldat)[3]=c("spcode")
alldat[is.na(alldat)] <- 0 

# species code check
uniquesp <- data.frame(unique(alldat$spcode))
names(uniquesp)[1]=c("spcode")

spkey <- read_excel(file.choose()) 
#choose file called "SpeciesCodes.xlsx" in our OneDrive folder, Tulare-LabPaper
names(spkey)[1:2]=c("spcode","spname")

diffsp <- full_join(spkey, uniquesp, by = "spcode")

# DEL SP is only species code that does not have a full plant name
delcheck <- alldat %>%
  filter(spcode == "DEL SP")

# link species key to alldat
alldatsp <- inner_join(alldat, spkey, by = "spcode")

# join datasets to include thermal treatment 
quadkey <- read_excel(file.choose())
#choose file called "MasterQuadratKey.xlsx" in our OneDrive folder: Tulare-LabPaper
names(quadkey)[3:4]=c("quadratNew", "quadrat")

# join alldatsp and quadkey for all data
alldatsptrt <- inner_join(alldatsp, quadkey, by = "quadrat") %>%
  select(quadratNew, treatment, spcode, spname, cover, year, thermal) %>%
  separate(treatment, sep = " ", c("burn", "graze")) %>%
  filter(graze != "NA", thermal != "?")

#check existence of all data across years
alldatcheck <- alldatsptrt %>%
  select(-cover, -spcode,-spname) %>%
  unique() %>%
  mutate(exist = 1) %>%
  spread(year, exist, fill="")

# print data frames
write.csv(alldatcheck, "alldatcheck.csv")
write.csv(alldatsptrt, "alldatsptrt.csv")

# join alldatsp and quadkey for only burned and grazed data
bugr <- inner_join(alldatsp, quadkey, by = "quadrat") %>%
  select(quadratNew, treatment, spcode, spname, cover, year, thermal) %>%
  filter(treatment == "burned grazed")

#check existence of bugr data across years
bugrcheck <- bugr %>%
  select(-cover, -spcode,-spname) %>%
  unique() %>%
  mutate(exist = 1) %>%
  spread(year, exist, fill="")

# print data frames
write.csv(bugrcheck, "BUGRdatacheck.csv")
write.csv(bugr, "bugrdat.csv")

# remove column "transect" from 2005-2018 to use rbind
dat_2005env <- dat2005env %>%
  select(-transect)

dat_2006env <- dat2006env %>%
  select(-transect)

dat_2007env <- dat2007env %>%
  select(-transect)

dat_2008env <- dat2008env %>%
  select(-transect)

dat_2009env <- dat2009env %>%
  select(-transect)

dat_2010env <- dat2010env %>%
  select(-transect)

dat_2011env <- dat2011env %>%
  select(-transect)

dat_2012env <- dat2012env %>%
  select(-transect)

dat_2013env <- dat2013env %>%
  select(-transect)

dat_2014env <- dat2014env %>%
  select(-transect)

dat_2015env <- dat2015env %>%
  select(-transect)

dat_2016env <- dat2016env %>%
  select(-transect)

dat_2017env <- dat2017env %>%
  select(-transect)

dat_2018env <- dat2018env %>%
  select(-transect)

#append all environmental data
alldatenv <- rbind(dat2001env, dat2002env, dat2003env, dat2004env, dat_2005env, dat_2006env, 
                dat_2007env, dat_2008env, dat_2009env, dat_2010env, dat_2011env, 
                dat_2012env, dat_2013env, dat_2014env, dat_2015env, dat_2016env, 
                dat_2017env, dat_2018env)
alldatenv[is.na(alldatenv)] <- 0 
alldatenvTH <- alldatenv %>%
  filter(site == "TH"| site == "TH-BG"| site == "TH-BUG"| site == "TH-UBUG"| site == "PGE"|
       site == "TH-MEC")

# join alldatenvTH and quadkey for all data
alldatenv <- inner_join(alldatenvTH, quadkey, by = "quadrat") %>%
  select(quadratNew, site, GOPHER, BARE, ROCK, LITTER, COWPIE, year)

write.csv(alldatenv, "envdat.csv")





