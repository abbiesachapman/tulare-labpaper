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
  select(transect.quad, transect, quadrat, year, graze, thermal, spname, cover)%>%
  filter(year==2006|year==2007|year==2008|year==2009|year==2010|year==2011|year==2012)
grztog2<-left_join(grztog, SC)

#plot timeseries of richness 
grzrich <- grztog2 %>%
  filter(func!="NA", status!="NA")%>%
  mutate(func=paste(func, status))%>%
  filter(cover != 0, spname != c("Unknown", "Moss")) %>%
  group_by(year, transect.quad, graze, func, thermal)%>%
  summarize(richness = length(unique(spname)))

grzrich1<-grzrich%>%
  group_by(year, graze, func) %>%
  summarize(mean_rich = mean(richness), se_rich=calcSE(richness))

ggplot(grzrich1, aes(year, mean_rich)) +
  geom_line(aes(color=as.factor(graze)))+facet_wrap(~func) +
  geom_errorbar(aes(ymin=mean_rich-se_rich, ymax=mean_rich+se_rich, color=as.factor(graze)), width=.2)+
  geom_vline(xintercept=2009)

#plot time series of shannon diversity
dat3<-dat1%>%
  filter(cover!=0)
dat3<-left_join(dat3, SC)%>%
  mutate(sitetrt=paste(sitetrt, status, func, sep="_"))
simp<-community_diversity(dat3, time.var = "year", abundance.var="cover", replicate.var="sitetrt", metric = c("InverseSimpson"))
shandiv<- community_diversity(dat3, time.var = "year", abundance.var="cover", replicate.var="sitetrt", metric = c("Shannon")) 

evenness<-left_join(simp, shandiv)%>%
  separate(sitetrt, into=c("site", "trt", "status", "func"), sep="_")%>%
  group_by(year, trt, status, func)%>%
  summarize(meanShan=mean(Shannon), meanSimp=mean(InverseSimpson))%>%
  filter(!is.na(status), !is.na(func), func!=("NA"), status!="NA")

ggplot(evenness, aes(year, meanShan)) +
  geom_line(aes(color=as.factor(interaction(status, func))))+facet_wrap(~trt)

ggplot(evenness, aes(year, meanSimp)) +
  geom_line(aes(color=as.factor(interaction(status, func))))+facet_wrap(~trt)

