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
grztog2<-left_join(grztog, SC)%>%
  filter(thermal=="moderate")

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

library(codyn)
#plot time series of shannon diversity
grzshan<-grztog2%>%
  filter(func!="NA", status!="NA")%>%
  mutate(func=paste(func, status))%>%
  filter(cover != 0, spname != c("Unknown", "Moss"))%>%
  mutate(alltrt=paste(transect.quad, graze, func, sep="_"))
simp.g<-community_diversity(grzshan, time.var = "year", abundance.var="cover", replicate.var="alltrt", metric = c("InverseSimpson"))
shandiv.g<- community_diversity(grzshan, time.var = "year", abundance.var="cover", replicate.var="alltrt", metric = c("Shannon")) 

grz.even<-left_join(simp.g, shandiv.g)%>%
  separate(alltrt, into=c("transect", "quad", "graze", "func"), sep="_")%>%
  mutate(transect.quad=paste(transect, quad, sep="_"))%>%
  group_by(year, graze, func)%>%
  summarize(meanShan=mean(Shannon), seShan=calcSE(Shannon), meanSimp=mean(InverseSimpson))%>%
  filter(!is.na(graze), !is.na(func), func!=("NA"), graze!="NA")

ggplot(grz.even, aes(year, meanShan)) +
  geom_line(aes(color=as.factor(graze)))+facet_wrap(~func) +
  geom_errorbar(aes(ymin=meanShan-seShan, ymax=meanShan+seShan, color=as.factor(graze)), width=.2)+
  geom_vline(xintercept=2009)+ylab("Shannon Diversity")

#plot timeseries of cover
grzfunc<-grztog2%>%
  group_by(transect.quad, year, graze, status, func, spcode)%>%
  summarize(sumcov=sum(cover))%>%
  filter(!is.na(status), !is.na(func))
grzfuncagg<-grzfunc%>%
  group_by(year, graze, status, func)%>%
  summarize(meancov=mean(sumcov), se_cov=calcSE(sumcov))
 
ggplot(grzfuncagg, aes((year), meancov))+
  geom_line(aes(color=graze))+
  geom_point(aes(color=graze))+
  geom_errorbar(aes(ymin=meancov-se_cov, ymax=meancov+se_cov, color=graze), width=.2)+
  facet_grid(status~func) +geom_vline(xintercept=2009)

##########
# Indicator Species
###########
library(indicspecies)
grztog3<-grztog2%>%
  mutate(prepost=ifelse(year<2009, "pre", "post"))%>%
  mutate(trtgroup=paste(graze, prepost, sep="_"))%>%
  mutate(rep=paste(transect.quad, year))%>%
  select(-transect, -transect.quad,-year, -quadrat, -graze, -prepost, -thermal, -spname, -status, -func)

indic_treatments<-select(grztog3, 3, 4)%>%
  unique()
indic_species<-select(grztog3, 1, 2, 4)%>%
  spread(spcode, cover, fill=0)%>%
  remove_rownames()%>%
  column_to_rownames("rep")
indic_species<-decostand(indic_species, "total")

indic_treatments2<-indic_species%>%
  rownames_to_column("rep")%>%
  select(1)
indic_treatments3<-left_join(indic_treatments2, indic_treatments)

indicators<-multipatt(indic_species, indic_treatments$trtgroup, func="IndVal.g", control=how(nperm=999))

indsum<-indicators$sign%>%
  rownames_to_column("species")%>%
  filter(p.value<.05)
indsum1<-left_join(indsum, SC, by=c("species"="spcode"))%>%
  select(-index, -stat, -p.value)

write.csv(indsum1, file = "Grazing_indicator_sp.csv")
