## Set your datpath!! (file in Data_cleaning)

library(tidyverse)
library(readr)
library(ggplot2)
library(vegan)

##FN for Calculating SE
calcSE<-function(x){
  x <- x[!is.na(x)]
  sd(x)/sqrt(length(x))
}

#load data
dat <- read_csv(paste(datpath_clean, "/bugrdat.csv", sep=""))
SC <- read_csv(paste(datpath_clean, "/SpeciesCodes.csv", sep=""))
clim <- read_csv(paste(datpath_clean, "/san_jose_clim.csv", sep =""), skip = 9) %>%
  mutate(stand_ppt = scale(PRCP, center = T, scale = T)) %>% #standardize to z-scores
  mutate(stand_temp = scale(TAVG, center = T, scale = T)) #standardize to z-scores
dat1<-dat%>%
  select(-1)%>%
  group_by(year, spname, quadratNew, thermal)%>%
  summarize(cover=sum(cover))%>%
  ungroup()%>%
  mutate(sitetrt=as.factor(paste(quadratNew, thermal, sep="_")))%>%
  filter(thermal!="?")
nitro <- read_csv(paste(datpath_clean, "/NTN-CA66-deposition.csv", sep = "")) %>%
  mutate(stand_NH4 = scale(NH4, center = T, scale = T)) %>%
  mutate(stand_NO3 = scale(NO3, center = T, scale = T)) %>%
  mutate(stand_N = scale(totalN, center = T, scale = T)) 
colnames(nitro)[colnames(nitro)=="yr"] <- "year"

#plot timeseries of richness 
dat2<-left_join(dat, SC)
richness0 <- dat2 %>%
  filter(func!="NA", status!="NA")%>%
  mutate(func=paste(func, status))%>%
  filter(cover != 0, spname != c("Unknown", "Moss")) %>%
  group_by(year, quadratNew, func, thermal)%>%
  summarize(richness = length(unique(spname)))
richness<-richness0%>%
  group_by(year,func, thermal) %>%
  filter(thermal!="?")%>%
  summarize(mean_rich = mean(richness), se_rich=calcSE(richness))

ggplot(richness, aes(year, mean_rich)) +
  geom_line(aes(color=as.factor(func)))+facet_wrap(~thermal) #+ geom_errorbar(aes((color=as.factor(func)), ymin=mean_rich-se_rich, ymax=mean_rich+se_rich), width=.2) +
  geom_line(aes(color=as.factor(func))) + ylab("mean species richness +/-se")

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
  

#join data and species key
tog<-left_join(dat, SC)
View(tog)

#plot timeseries of cover by functional group
functog<-tog%>%
  group_by(quadratNew, year, status, func, thermal)%>%
  summarize(sumcov=sum(cover))%>%
  filter(!is.na(status), !is.na(func))
functogagg<-functog%>%
  group_by(year, status, func)%>%
  summarize(meancov=mean(sumcov), se_cov=calcSE(sumcov))

View(functoagg)
View(functog)

ggplot(functog, aes(as.factor(year), sumcov))+
  geom_boxplot()+
  facet_grid(status~func)
ggplot(functogagg, aes((year), meancov))+
  geom_line(aes(color=interaction(status, func)))+
  geom_point(aes(color=interaction(status, func)))+
  geom_errorbar(aes(ymin=meancov-se_cov, ymax=meancov+se_cov, color=interaction(status, func)), width=.2)

#plot timeseries of cover by functional group by thermal
thermal_trend <- functog %>%
  group_by(year, status, func, thermal) %>%
  summarize(meancov = mean(sumcov), se_cov=calcSE(sumcov)) %>%
  filter(thermal != "?")

ggplot(thermal_trend, aes(year, meancov)) +
  geom_line(aes(color = thermal)) +
  geom_point(aes(color = thermal)) +
  facet_grid(status~func)

#plot timeseries of richness by functional group
richtog <- tog %>%
  filter(cover != 0, spname != c("Unknown", "Moss")) %>%
  group_by(quadratNew, year, status, func, thermal) %>%
  summarize(richness = length(unique(spname))) %>%
  filter(status != "NA") 
richtogagg <- richtog %>%
  group_by(year, status, func) %>%
  summarize(meanrich = mean(richness))

ggplot(richtogagg, aes(year, meanrich)) +
  geom_line(aes(color = interaction(status, func))) +
  geom_point(aes(color = interaction(status, func)))

#plot timeseries of richness by functional group and thermal
thermal_rich <- richtog %>%
  filter(thermal != "?") %>%
  group_by(year, status, func, thermal) %>%
  summarize(meanrich = mean(richness))

ggplot(thermal_rich, aes(year, meanrich)) +
  geom_line(aes(color = thermal)) +
  geom_point(aes(color = thermal)) +
  facet_grid(status~func)

#timeseries of ppt
ggplot(clim, aes(DATE,PRCP)) + geom_line() + geom_point()

#cover by thermal transposed on annual precip
ggplot(thermal_trend, aes((year), meancov)) + facet_grid(status~func) +
  geom_bar(data = clim, aes(x = DATE, y = PRCP*3), stat = "identity", fill = "lightgrey") +
  geom_line(aes(color= thermal))+ geom_point(aes(color = thermal))  +
  geom_errorbar(aes(ymin=meancov-se_cov, ymax=meancov+se_cov, color=thermal), width=.2)+
  scale_y_continuous(sec.axis = sec_axis(~./3, name = "Annual Precipitation in inches"))

#cover by thermal transpoesd on annual mean temp
ggplot(thermal_trend, aes((year), meancov)) + facet_grid(status~func) +
  geom_bar(data = clim, aes(x = DATE, y = stand_temp*30), stat = "identity", fill = "lightgrey") +
  geom_line(aes(color= thermal))+ geom_point(aes(color = thermal))  +
  geom_errorbar(aes(ymin=meancov-se_cov, ymax=meancov+se_cov, color=thermal), width=.2)+
  scale_y_continuous(sec.axis = sec_axis(~./30, name = "Annual Mean Temp Deviation (z-scores)"))

#join clim and thermal_trend
clim1 <- clim %>%
  rename(year = DATE)
joined_dat <- left_join(thermal_trend, clim1)

#cover by temp
ggplot(joined_dat, aes((PRCP), meancov)) + facet_grid(status~func) +
   geom_point(aes(color = thermal)) + geom_smooth(aes(color = thermal), method = "lm", se = F)

#lag effect of clim?
clim_lag <- clim1 %>%
  mutate(year = year + 1)
clim_lag <- left_join(thermal_trend, clim_lag)  

ggplot(clim_lag, aes((year), meancov)) + facet_grid(status~func) +
  geom_bar(data = clim_lag, aes(x = year, y = stand_ppt*5), stat = "identity", fill = "lightgrey") +
  geom_line(aes(color= thermal))+ geom_point(aes(color = thermal))  +
  scale_y_continuous(sec.axis = sec_axis(~./20, name = "Previous Year's Precipitation Deviation (z-scores)"))

ggplot(clim_lag, aes((year), meancov)) + facet_grid(status~func) +
  geom_bar(data = clim_lag, aes(x = year, y = stand_temp*5), stat = "identity", fill = "lightgrey") +
  geom_line(aes(color= thermal))+ geom_point(aes(color = thermal))  +
  scale_y_continuous(sec.axis = sec_axis(~./15, name = "Previous Year's Temperature Deviation (z-scores)"))

#join nitro and thermal trend 
nitro_dat <- left_join(thermal_trend, nitro)

#plot cover by thermal transposed on nitrogen
ggplot(thermal_trend, aes(year, meancov)) +
  facet_grid(status~func) +
  geom_bar(data = nitro, aes(x = year, y = stand_NH4*20), stat = "identity", fill = "lightgrey") +
  geom_line(aes(color = thermal)) + 
  geom_point(aes(color = thermal)) +
  scale_y_continuous(sec.axis = sec_axis(~./30, name = "Annual Mean NH4 Deviation (z-scores)"))

ggplot(thermal_trend, aes(year, meancov)) +
  facet_grid(status~func) +
  geom_bar(data = nitro, aes(x = year, y = stand_NO3*20), stat = "identity", fill = "lightgrey") +
  geom_line(aes(color = thermal)) + 
  geom_point(aes(color = thermal)) +
  scale_y_continuous(sec.axis = sec_axis(~./30, name = "Annual Mean NO3 Deviation (z-scores)"))

#cover by nitrogen
ggplot(nitro_dat, aes(NH4, meancov)) +
  facet_grid(status~func) +
  geom_point(aes(color = thermal)) +
  geom_smooth(aes(color = thermal), method = "lm", se = F)

ggplot(nitro_dat, aes(NO3, meancov)) +
  facet_grid(status~func) +
  geom_point(aes(color = thermal)) +
  geom_smooth(aes(color = thermal), method = "lm", se = F)

#lag effect of nitrogen?
nitro_lag <- nitro %>%
  mutate(year = year + 1)
nitro_lag <- left_join(thermal_trend, nitro_lag)  

ggplot(nitro_lag, aes(NH4, meancov)) +
  facet_grid(status~func) +
  geom_point(aes(color = thermal)) +
  geom_smooth(aes(color = thermal), method = "lm", se = F)

ggplot(nitro_lag, aes(NO3, meancov)) +
  facet_grid(status~func) +
  geom_point(aes(color = thermal)) +
  geom_smooth(aes(color = thermal), method = "lm", se = F)

#year by year species turnover
  turnover1<-turnover(dat1,
                      time.var="year",
                      species.var="spname",
                      abundance.var="cover",
                      replicate.var="sitetrt")%>%
  separate(sitetrt, into=c("site", "trt"), sep="_")%>%
  group_by( trt, year)%>%
  summarize(mean.turnover=mean(as.numeric(total)), se.turnover=calcSE(as.numeric(total)))

turnoverg<-turnover(dat1,
                    time.var="year",
                    species.var="spname",
                    abundance.var="cover",
                    replicate.var="sitetrt", metric="appearance")%>%
  separate(sitetrt, into=c("site", "trt"), sep="_")%>%
  group_by( trt, year)%>%
  summarize(mean.app=mean(as.numeric(appearance)), se.app=calcSE(as.numeric(appearance)))

turnoverl<-turnover(dat1,
                    time.var="year",
                    species.var="spname",
                    abundance.var="cover",
                    replicate.var="sitetrt", metric="disappearance")%>%
  separate(sitetrt, into=c("site", "trt"), sep="_")%>%
  group_by( trt, year)%>%
  summarize(mean.dapp=mean(as.numeric(disappearance)), se.dapp=calcSE(as.numeric(disappearance)))

turnover<-left_join(turnover1, turnoverg)
turnover<-left_join(turnover, turnoverl)
rm(turnover1, turnoverg, turnoverl)

ggplot(turnover)+
  geom_point(aes((year), mean.turnover, color=trt))+
  geom_line(aes((year), mean.turnover, color=trt))+
  labs(x="Year", y="Annual Turnover")+
  geom_errorbar(aes(color=trt, (year), ymin=mean.turnover-se.turnover, ymax=mean.turnover+se.turnover ), width=.2)

ggplot(turnover)+
  geom_point(aes((year), mean.app, color=trt))+
  geom_line(aes((year), mean.app, color=trt))+
  labs(x="Year", y="Annual Gains")+
  geom_errorbar(aes(color=trt, (year), ymin=mean.app-se.app, ymax=mean.app+se.app ), width=.2)

ggplot(turnover)+
  geom_point(aes((year), mean.dapp, color=trt))+
  geom_line(aes((year), mean.dapp, color=trt))+
  labs(x="Year", y="Annual Losses")+
  geom_errorbar(aes(color=trt, (year), ymin=mean.dapp-se.dapp, ymax=mean.dapp+se.dapp ), width=.2)

ggplot(turnover, aes(x=trt, y=mean.turnover))+geom_boxplot()


# year by year species rank abundance shift
rankshift <- rank_shift(dat1,
                        time.var="year",
                        species.var="spname",
                        abundance.var="cover",
                        replicate.var="sitetrt")%>%
  separate(sitetrt, into=c("site", "trt"), sep="_")%>%
  group_by( trt, year_pair)%>%
  summarize(mean.MRS=mean(MRS), se.MRS=calcSE(MRS))%>%
  mutate(year=substr(year_pair, 6,9))

ggplot(rankshift)+
  geom_line(aes(as.factor(year), mean.MRS, group=trt, color=trt))+
  geom_errorbar(aes(color=trt, as.factor(year), ymin=mean.MRS-se.MRS, ymax=mean.MRS+se.MRS ), width=.1)+
  labs(x="Year", y="Mean Rank Shift")

meanmeanrank<-rankshift%>%
  group_by(trt)%>%
  summarize(meanmean=mean(mean.MRS), se=calcSE(mean.MRS))

ggplot(meanmeanrank, aes(x=trt, y=meanmean))+geom_bar(stat="identity")

ggplot(rankshift, aes(x=trt, y=mean.MRS))+geom_boxplot()
