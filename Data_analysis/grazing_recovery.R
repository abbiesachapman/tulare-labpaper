#run after BUGR_timeseries.R
# Objective: compare grazed and ungrazed plots as grazing was reintroduced
# 2006-2008 = ungrazed plots are ungrazed, 2009-2012 = ungrazed plots have cattle reintroduced
# all years = grazed plots are grazed

alldat<-read_csv(paste(datpath_clean, "/alldatsp.csv", sep="")) %>%
  select(-1)%>%
  group_by(year, spname, quadrat)%>%
  summarize(cover=sum(cover))%>%
  filter(cover!=0)%>%
  separate(quadrat, into=c("transect", "quadrat"), sep="-")

#treatments is a dataframe of burning, grazing, thermal treatment for each quadrat, and which years have data
treatments<-read_csv(paste(datpath_clean, "/quadrat_trt.csv", sep=""))%>%
  gather(key="year", value="data", -quadratID, -burn, -graze, -thermal)
#treatments612 is only quadrats with data in 2006-2012
treatments_612<-treatments%>%
  filter(year==2006|year==2007|year==2008|year==2009|year==2010|year==2011|year==2012)%>%
  mutate(year=paste("y", year, sep=""))%>%
  spread(year, data)%>%
  filter(!is.na(y2006),!is.na(y2007),!is.na(y2008),!is.na(y2009),!is.na(y2010),!is.na(y2011),!is.na(y2012))%>%
  select(1:4)%>%
  mutate(transect=quadratID)%>%
  select(-1)

grztog<-left_join(treatments_612, alldat)%>%
  mutate(transect.quad=paste(transect, quadrat, sep="_"))%>%
  select(transect.quad, transect, quadrat, year, graze, thermal, spname, cover)
                       