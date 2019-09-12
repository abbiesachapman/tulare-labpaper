library(tidyverse)
library(readr)
library(nlme)

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
########

#Richness
rich2 <- alldat %>%
  filter(type!="NA", status!="NA")%>%
  filter(cover != 0, spname !="Unknown", spname!="Moss") %>%
  mutate(func=paste(type, status))%>%
  mutate(trt=paste(graze, burn))%>%
  group_by(year, quadratNew, func, trt, type, status, burn, graze, transect)%>%
  summarize(richness = length(unique(spname)))%>%
  mutate(prepost=ifelse(year<2009, "pre", "post"))

rich2005 <- rich2 %>%
  filter(trt != "ungrazed unburned") %>%
  filter(year == 2005)
anova(lm(richness~graze, data = rich2005%>%filter(func == "forb native")))
anova(lm(richness~graze, data = rich2005%>%filter(func == "grass non-native")))
anova(lm(richness~graze, data = rich2005%>%filter(func == "forb non-native")))
anova(lm(richness~graze, data = rich2005%>%filter(func == "grass native")))

rich2006 <- rich2 %>%
  filter(trt != "ungrazed unburned") %>%
  filter(year == 2006)
anova(lm(richness~graze, data = rich2006%>%filter(func == "forb native")))
anova(lm(richness~graze, data = rich2006%>%filter(func == "grass non-native")))
anova(lm(richness~graze, data = rich2006%>%filter(func == "forb non-native")))
anova(lm(richness~graze, data = rich2006%>%filter(func == "grass native")))

rich2007 <- rich2 %>%
  filter(trt != "ungrazed unburned") %>%
  filter(year == 2007)
anova(lm(richness~graze, data = rich2007%>%filter(func == "forb native")))
anova(lm(richness~graze, data = rich2007%>%filter(func == "grass non-native")))
anova(lm(richness~graze, data = rich2007%>%filter(func == "forb non-native")))
anova(lm(richness~graze, data = rich2007%>%filter(func == "grass native")))

rich2008 <- rich2 %>%
  filter(trt != "ungrazed unburned") %>%
  filter(year == 2008)
anova(lm(richness~graze, data = rich2008%>%filter(func == "forb native")))
anova(lm(richness~graze, data = rich2008%>%filter(func == "grass non-native")))
anova(lm(richness~graze, data = rich2008%>%filter(func == "forb non-native")))
anova(lm(richness~graze, data = rich2008%>%filter(func == "grass native")))

rich2009 <- rich2 %>%
  filter(year == 2009)
anova(lm(richness~trt, data = rich2009%>%filter(func == "forb native")))
anova(lm(richness~trt, data = rich2009%>%filter(func == "grass non-native")))
anova(lm(richness~trt, data = rich2009%>%filter(func == "forb non-native")))
anova(lm(richness~trt, data = rich2009%>%filter(func == "grass native")))

rich2010 <- rich2 %>%
  filter(year == 2010)
anova(lm(richness~trt, data = rich2010%>%filter(func == "forb native")))
anova(lm(richness~trt, data = rich2010%>%filter(func == "grass non-native")))
anova(lm(richness~trt, data = rich2010%>%filter(func == "forb non-native")))
anova(lm(richness~trt, data = rich2010%>%filter(func == "grass native")))

rich2011 <- rich2 %>%
  filter(year == 2011)
anova(lm(richness~trt, data = rich2011%>%filter(func == "forb native")))
anova(lm(richness~trt, data = rich2011%>%filter(func == "grass non-native")))
anova(lm(richness~trt, data = rich2011%>%filter(func == "forb non-native")))
anova(lm(richness~trt, data = rich2011%>%filter(func == "grass native")))

rich2012 <- rich2 %>%
  filter(year == 2012)
anova(lm(richness~trt, data = rich2012%>%filter(func == "forb native")))
anova(lm(richness~trt, data = rich2012%>%filter(func == "grass non-native")))
anova(lm(richness~trt, data = rich2012%>%filter(func == "forb non-native")))
anova(lm(richness~trt, data = rich2012%>%filter(func == "grass native")))

#cover
cov<-alldat%>%
  filter(type!="NA", status!="NA")%>%
  mutate(func=paste(type, status))%>%
  mutate(trt=paste(graze, burn))%>%
  filter(cover != 0, spname !="Unknown", spname!="Moss") %>%
  group_by(year, quadratNew, trt, func, spcode, spname, burn, graze, transect)%>%
  summarize(sumcov=sum(cover))%>%
  filter(!is.na(trt), !is.na(func))%>%
  mutate(prepost=ifelse(year<2009, "pre", "post"))

cov2005 <- cov %>%
  filter(year == 2005) %>%
  filter(trt != "ungrazed unburned")
anova(lm(sumcov~graze, data = cov2005%>%filter(func == "forb native")))
anova(lm(sumcov~graze, data = cov2005%>%filter(func == "grass non-native")))
anova(lm(sumcov~graze, data = cov2005%>%filter(func == "forb non-native")))
anova(lm(sumcov~graze, data = cov2005%>%filter(func == "grass native")))

cov2006 <- cov %>%
  filter(year == 2006) %>%
  filter(trt != "ungrazed unburned")
anova(lm(sumcov~graze, data = cov2006%>%filter(func == "forb native")))
anova(lm(sumcov~graze, data = cov2006%>%filter(func == "grass non-native")))
anova(lm(sumcov~graze, data = cov2006%>%filter(func == "forb non-native")))
anova(lm(sumcov~graze, data = cov2006%>%filter(func == "grass native")))

cov2007 <- cov %>%
  filter(year == 2007) %>%
  filter(trt != "ungrazed unburned")
anova(lm(sumcov~graze, data = cov2007%>%filter(func == "forb native")))
anova(lm(sumcov~graze, data = cov2007%>%filter(func == "grass non-native")))
anova(lm(sumcov~graze, data = cov2007%>%filter(func == "forb non-native")))
anova(lm(sumcov~graze, data = cov2007%>%filter(func == "grass native")))

cov2008 <- cov %>%
  filter(year == 2008) %>%
  filter(trt != "ungrazed unburned")
anova(lm(sumcov~graze, data = cov2008%>%filter(func == "forb native")))
anova(lm(sumcov~graze, data = cov2008%>%filter(func == "grass non-native")))
anova(lm(sumcov~graze, data = cov2008%>%filter(func == "forb non-native")))
anova(lm(sumcov~graze, data = cov2008%>%filter(func == "grass native")))

cov2009 <- cov %>%
  filter(year == 2009)
anova(lm(sumcov~trt, data = cov2009%>%filter(func == "forb native")))
anova(lm(sumcov~trt, data = cov2009%>%filter(func == "grass non-native")))
anova(lm(sumcov~trt, data = cov2009%>%filter(func == "forb non-native")))
anova(lm(sumcov~trt, data = cov2009%>%filter(func == "grass native")))

cov2010 <- cov %>%
  filter(year == 2010)
anova(lm(sumcov~trt, data = cov2010%>%filter(func == "forb native")))
anova(lm(sumcov~trt, data = cov2010%>%filter(func == "grass non-native")))
anova(lm(sumcov~trt, data = cov2010%>%filter(func == "forb non-native")))
anova(lm(sumcov~trt, data = cov2010%>%filter(func == "grass native")))

cov2011 <- cov %>%
  filter(year == 2011)
anova(lm(sumcov~trt, data = cov2011%>%filter(func == "forb native")))
anova(lm(sumcov~trt, data = cov2011%>%filter(func == "grass non-native")))
anova(lm(sumcov~trt, data = cov2011%>%filter(func == "forb non-native")))
anova(lm(sumcov~trt, data = cov2011%>%filter(func == "grass native")))

cov2012 <- cov %>%
  filter(year == 2012)
anova(lm(sumcov~trt, data = cov2012%>%filter(func == "forb native")))
anova(lm(sumcov~trt, data = cov2012%>%filter(func == "grass non-native")))
anova(lm(sumcov~trt, data = cov2012%>%filter(func == "forb non-native")))
anova(lm(sumcov~trt, data = cov2012%>%filter(func == "grass native")))

########
#one-way ANOVA analyze two blocks of years: 2005-2008 and 2008-2012 
########
rich3 <- rich2 %>%
  filter(trt != "ungrazed unburned") %>%
  filter(year%in%c(2005, 2006, 2007, 2008))

anova(lm(richness~graze, data = rich3%>%filter(func == "forb native")))
anova(lm(richness~graze, data = rich3%>%filter(func == "grass non-native")))
anova(lm(richness~graze, data = rich3%>%filter(func == "forb non-native")))
anova(lm(richness~graze, data = rich3%>%filter(func == "grass native")))

rich4 <- rich2 %>%
  filter(year%in%c(2008:2012))

anova(lm(richness~trt, data = rich4%>%filter(func == "forb native")))
anova(lm(richness~trt, data = rich4%>%filter(func == "grass non-native")))
anova(lm(richness~trt, data = rich4%>%filter(func == "forb non-native")))
anova(lm(richness~trt, data = rich4%>%filter(func == "grass native")))

########
#linear mixed models
#grazed and burned as fixed effects, year as random
#nest quadrats within transects
########

fit2005fn<- lme(richness~graze, random = ~1|transect/quadratNew, data = rich2005%>%filter(func == "forb native"))
summary(fit2005fn)
fit2005ge<- lme(richness~graze, random = ~1|transect/quadratNew, data = rich2005%>%filter(func == "grass non-native"))
summary(fit2005ge)

fitpostfn <- lme(richness~graze+burn, random = ~1|transect/quadratNew, data = rich4%>%filter(func == "forb native"))
summary(fitpostfn)
