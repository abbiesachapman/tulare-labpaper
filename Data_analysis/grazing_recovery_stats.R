library(tidyverse)
library(readr)

#load updated master data
alldat<-read_csv(paste(datpath_clean, "/alldatsptrt.csv", sep="")) %>%
  select(-1)%>%
  filter(transect%in%c("THBUGM1", "THBUGM2", "THM1", "THM2", "THM3", "THM4", "THUBUGM1", "THUBUGM2"))%>%
  filter(thermal=="moderate")%>%
  group_by(year, spname, spcode, quadratNew, status, type, transect, burn, graze)%>%
  summarize(cover=sum(cover))%>%
  filter(cover!=0)

########
#one-way ANOVA analyze within each year 
#richness, cover, evenness, litter as dependent variables
#focus on two functional groups: native forbs and non-native grasses
#nest quadrats within transects
########



########
#one-way ANOVA analyze two blocks of years: 2005-2008 and 2008-2012 
########



########
#linear mixed models
#grazed and burned as fixed effects, year as random
#nest quadrats within transects
########
