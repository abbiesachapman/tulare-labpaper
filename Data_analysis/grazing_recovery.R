#run after BUGR_timeseries.R
# Objective: compare grazed and ungrazed plots as grazing was reintroduced
# 2006-2008 = ungrazed plots are ungrazed, 2009-2012 = ungrazed plots have cattle reintroduced
# all years = grazed plots are grazed
library(ggplot2); theme_set(theme_bw())

alldat<-read_csv(paste(datpath_clean, "/alldatsptrt.csv", sep="")) %>%
  select(-1)%>%
  filter(transect%in%c("THBUGM1", "THBUGM2", "THM1", "THM2", "THM3", "THM4", "THUBUGM1", "THUBUGM2"))%>%
  filter(!quadratNew%in%c("THM1-1", "THM3-3", "THM1-10"))%>%
  filter(thermal=="moderate")%>%
  group_by(year, spname, spcode, quadratNew, status, type, transect, burn, graze)%>%
  summarize(cover=sum(cover))%>%
  filter(cover!=0)

#plot timeseries of richness 
rich <- alldat %>%
  filter(type!="NA", status!="NA")%>%
  mutate(func=paste(type, status))%>%
  mutate(trt=paste(graze, burn))%>%
  filter(cover != 0, spname !="Unknown", spname!="Moss") %>%
  group_by(year, quadratNew, trt, func)%>%
  summarize(richness = length(unique(spname)))%>%
  mutate(prepost=ifelse(year<2009, "pre", "post"))

rich1<-rich%>%
  group_by(year, trt, func, prepost) %>%
  summarize(mean_rich = mean(richness), se_rich=calcSE(richness))

ggplot(rich1, aes(year, mean_rich)) +
  geom_line(aes(color=as.factor(trt)))+
  geom_point(aes(color=as.factor(trt)))+
  facet_wrap(~func, scales="free") +
  geom_errorbar(aes(ymin=mean_rich-se_rich, ymax=mean_rich+se_rich, color=as.factor(trt)), width=.2)+
  geom_vline(xintercept=2008.5)+geom_vline(xintercept=2004.5, color="red") +
  labs(x = "Year", y = "Mean Species Richness", color = "Treatment")

library(codyn)
#plot time series of shannon diversity
shan<-alldat%>%
  filter(type!="NA", status!="NA")%>%
  mutate(func=paste(type, status))%>%
  mutate(trt=paste(graze, burn))%>%
  filter(cover != 0, spname !="Unknown", spname!="Moss") %>%
  mutate(alltrt=paste(quadratNew, trt, func, sep="_"))
simp.g<-community_diversity(shan, time.var = "year", abundance.var="cover", replicate.var="alltrt", metric = c("InverseSimpson"))
shandiv.g<- community_diversity(shan, time.var = "year", abundance.var="cover", replicate.var="alltrt", metric = c("Shannon")) 

shan2<-left_join(simp.g, shandiv.g)%>%
  separate(alltrt, into=c("transectNew", "trt", "func"), sep="_")%>%
  mutate(prepost=ifelse(year<2009, "pre", "post"))

shan1<-shan2%>%
  group_by(year, trt, func, prepost)%>%
  summarize(meanShan=mean(Shannon), seShan=calcSE(Shannon), meanSimp=mean(InverseSimpson))%>%
  filter(!is.na(trt), !is.na(func), func!=("NA"), trt!="NA")

ggplot(shan1, aes(year, meanShan)) +
  geom_line(aes(color=as.factor(trt)))+
  geom_point(aes(color=as.factor(trt)))+
  facet_wrap(~func) +
  geom_errorbar(aes(ymin=meanShan-seShan, ymax=meanShan+seShan, color=as.factor(trt)), width=.2)+
  geom_vline(xintercept=2008.5)+geom_vline(xintercept=2004.5, color="red")+ylab("Shannon Diversity") +
  labs(x = "Year", y = "Mean Shannon Diversity", color = "Treatment")

#plot timeseries of cover
cov<-alldat%>%
  filter(type!="NA", status!="NA")%>%
  mutate(func=paste(type, status))%>%
  mutate(trt=paste(graze, burn))%>%
  filter(cover != 0, spname !="Unknown", spname!="Moss") %>%
  group_by(year, quadratNew, trt, func)%>%
  summarize(sumcov=sum(cover))%>%
  filter(!is.na(trt), !is.na(func))%>%
  mutate(prepost=ifelse(year<2009, "pre", "post"))%>%
  group_by(year, prepost, quadratNew, trt)%>%
  mutate(totcov=sum(sumcov))%>%
  mutate(relcov=sumcov/totcov)
cov1<-cov%>%
  group_by(year, trt, func, prepost)%>%
  summarize(meancov=mean(sumcov), se_cov=calcSE(sumcov), meanrelcov=mean(relcov), se_relcov=calcSE(relcov))


ggplot(cov1, aes((year), meancov))+
  geom_line(aes(color=trt))+
  geom_point(aes(color=trt))+
  geom_errorbar(aes(ymin=meancov-se_cov, ymax=meancov+se_cov, color=trt), width=.2)+
  facet_wrap(~func) +geom_vline(xintercept=2008.5)+geom_vline(xintercept=2004.5, color="red") +
  labs(x = "Year", y = "Mean Cover (%)", color = "Treatment")

ggplot(cov1, aes((year), meanrelcov))+
  geom_line(aes(color=trt))+
  geom_point(aes(color=trt))+
  geom_errorbar(aes(ymin=meanrelcov-se_relcov, ymax=meanrelcov+se_relcov, color=trt), width=.2)+
  facet_wrap(~func) +geom_vline(xintercept=2008.5)+geom_vline(xintercept=2004.5, color="red") +
  labs(x = "Year", y = "Mean Relative Cover", color = "Treatment")

##########
#Litter
##########

#load environmental data
envdat <- read_csv(paste(datpath_clean, "/envdat.csv", sep = "")) %>%
  select(-1) 

#join env data with master
joindat <- left_join(alldat, envdat)

calcSE<-function(x){
  x <- x[!is.na(x)]
  sd(x)/sqrt(length(x))
}

#take avg litter and se
litter <- joindat %>%
  filter(type != "NA", status != "NA") %>% #remove NAs
  mutate(trt=paste(graze, burn)) %>%
  group_by(trt, year) %>%
  summarise(mean_litter = mean(litter),
            se_litter = calcSE(litter))

#graph litter time series
ggplot(litter, aes(year, mean_litter)) +
  geom_line(aes(color=as.factor(trt))) +
  geom_point(aes(color=as.factor(trt))) +
  geom_errorbar(aes(ymin=mean_litter-se_litter, ymax=mean_litter+se_litter, color=as.factor(trt)), width=.2) +
  geom_vline(xintercept=2008.5)+geom_vline(xintercept=2004.5, color="red") +
  labs(x = "Year", y = "Mean Litter Cover (%)", color = "Treatment")


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