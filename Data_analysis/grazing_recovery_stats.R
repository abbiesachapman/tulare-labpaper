library(tidyverse)
library(readr)
library(sjstats) #standardized effect size cohen's f
library(nlme) #linear mixed models

#load updated master data
alldat<-read_csv(paste(datpath_clean, "/alldatsptrt.csv", sep="")) %>%
  select(-1)%>%
  filter(transect%in%c("THBUGM1", "THBUGM2", "THM1", "THM2", "THM3", "THM4", "THUBUGM1", "THUBUGM2"))%>%
  filter(!quadratNew%in%c("THM1-1", "THM3-3", "THM1-10"))%>%
  filter(thermal=="moderate")%>%
  group_by(year, spname, spcode, quadratNew, status, type, transect, burn, graze)%>%
  summarize(cover=sum(cover))

########
#one-way ANOVA analyze within each year 
#richness, cover, evenness, litter as dependent variables
#focus on two functional groups: native forbs and non-native grasses
########

#Richness
rich2 <- alldat %>%
  filter(type!="NA", status!="NA")%>%
  filter(spname !="Unknown", spname!="Moss") %>%
  mutate(func=paste(type, status))%>%
  mutate(trt=paste(graze, burn))%>%
  mutate(present=ifelse(cover>0, 1, 0))%>%
  group_by(year, quadratNew, func, trt, type, status, burn, graze, transect, present)%>%
  summarize(richness = sum(present))

TukeyHSD(aov(richness~trt, data = rich2%>%filter(year == 2005, func == "forb native")))
TukeyHSD(aov(richness~trt, data = rich2%>%filter(year == 2006, func == "forb native")))
TukeyHSD(aov(richness~trt, data = rich2%>%filter(year == 2007, func == "forb native")))
TukeyHSD(aov(richness~trt, data = rich2%>%filter(year == 2008, func == "forb native")))
TukeyHSD(aov(richness~trt, data = rich2%>%filter(year == 2009, func == "forb native")))
TukeyHSD(aov(richness~trt, data = rich2%>%filter(year == 2010, func == "forb native")))
TukeyHSD(aov(richness~trt, data = rich2%>%filter(year == 2011, func == "forb native")))
TukeyHSD(aov(richness~trt, data = rich2%>%filter(year == 2012, func == "forb native")))

rich2005 <- rich2 %>%
  filter(trt != "ungrazed unburned") %>%
  filter(year == 2005)
anova_stats(anova(lm(richness~graze, data = rich2005%>%filter(func == "forb native"))))
anova_stats(anova(lm(richness~graze, data = rich2005%>%filter(func == "grass non-native"))))
anova_stats(anova(lm(richness~graze, data = rich2005%>%filter(func == "forb non-native"))))
anova_stats(anova(lm(richness~graze, data = rich2005%>%filter(func == "grass native"))))

rich2006 <- rich2 %>%
  filter(trt != "ungrazed unburned") %>%
  filter(year == 2006)
anova_stats(anova(lm(richness~graze, data = rich2006%>%filter(func == "forb native"))))
anova_stats(anova(lm(richness~graze, data = rich2006%>%filter(func == "grass non-native"))))
anova_stats(anova(lm(richness~graze, data = rich2006%>%filter(func == "forb non-native"))))
anova_stats(anova(lm(richness~graze, data = rich2006%>%filter(func == "grass native"))))

rich2007 <- rich2 %>%
  filter(trt != "ungrazed unburned") %>%
  filter(year == 2007)
anova_stats(anova(lm(richness~graze, data = rich2007%>%filter(func == "forb native"))))
anova_stats(anova(lm(richness~graze, data = rich2007%>%filter(func == "grass non-native"))))
anova_stats(anova(lm(richness~graze, data = rich2007%>%filter(func == "forb non-native"))))
anova_stats(anova(lm(richness~graze, data = rich2007%>%filter(func == "grass native"))))

rich2008 <- rich2 %>%
  filter(trt != "ungrazed unburned") %>%
  filter(year == 2008)
anova_stats(anova(lm(richness~graze, data = rich2008%>%filter(func == "forb native"))))
anova_stats(anova(lm(richness~graze, data = rich2008%>%filter(func == "grass non-native"))))
anova_stats(anova(lm(richness~graze, data = rich2008%>%filter(func == "forb non-native"))))
anova_stats(anova(lm(richness~graze, data = rich2008%>%filter(func == "grass native"))))

rich2009 <- rich2 %>%
  filter(year == 2009)
anova_stats(anova(lm(richness~trt, data = rich2009%>%filter(func == "forb native"))))
anova_stats(anova(lm(richness~trt, data = rich2009%>%filter(func == "grass non-native"))))
anova_stats(anova(lm(richness~trt, data = rich2009%>%filter(func == "forb non-native"))))
anova_stats(anova(lm(richness~trt, data = rich2009%>%filter(func == "grass native"))))

rich2010 <- rich2 %>%
  filter(year == 2010)
anova_stats(anova(lm(richness~trt, data = rich2010%>%filter(func == "forb native"))))
anova_stats(anova(lm(richness~trt, data = rich2010%>%filter(func == "grass non-native"))))
anova_stats(anova(lm(richness~trt, data = rich2010%>%filter(func == "forb non-native"))))
anova_stats(anova(lm(richness~trt, data = rich2010%>%filter(func == "grass native"))))

rich2011 <- rich2 %>%
  filter(year == 2011)
anova_stats(anova(lm(richness~trt, data = rich2011%>%filter(func == "forb native"))))
anova_stats(anova(lm(richness~trt, data = rich2011%>%filter(func == "grass non-native"))))
anova_stats(anova(lm(richness~trt, data = rich2011%>%filter(func == "forb non-native"))))
anova_stats(anova(lm(richness~trt, data = rich2011%>%filter(func == "grass native"))))

rich2012 <- rich2 %>%
  filter(year == 2012)
anova_stats(anova(lm(richness~trt, data = rich2012%>%filter(func == "forb native"))))
anova_stats(anova(lm(richness~trt, data = rich2012%>%filter(func == "grass non-native"))))
anova_stats(anova(lm(richness~trt, data = rich2012%>%filter(func == "forb non-native"))))
anova_stats(anova(lm(richness~trt, data = rich2012%>%filter(func == "grass native"))))

#cover
cov<-alldat%>%
  filter(type!="NA", status!="NA")%>%
  mutate(func=paste(type, status))%>%
  mutate(trt=paste(graze, burn))%>%
  filter(spname !="Unknown", spname!="Moss") %>%
  group_by(year, quadratNew, trt, func, spcode, spname, burn, graze, transect)%>%
  summarize(sumcov=sum(cover))%>%
  filter(!is.na(trt), !is.na(func)) %>%
  mutate(totcov=sum(sumcov))%>%
  mutate(relcov=sumcov/totcov)

cov2005 <- cov %>%
  filter(year == 2005) %>%
  filter(trt != "ungrazed unburned")
anova_stats(anova(lm(sumcov~graze, data = cov2005%>%filter(func == "forb native"))))
anova_stats(anova(lm(sumcov~graze, data = cov2005%>%filter(func == "grass non-native"))))
anova_stats(anova(lm(sumcov~graze, data = cov2005%>%filter(func == "forb non-native"))))
anova_stats(anova(lm(sumcov~graze, data = cov2005%>%filter(func == "grass native"))))

cov2006 <- cov %>%
  filter(year == 2006) %>%
  filter(trt != "ungrazed unburned")
anova_stats(anova(lm(sumcov~graze, data = cov2006%>%filter(func == "forb native"))))
anova_stats(anova(lm(sumcov~graze, data = cov2006%>%filter(func == "grass non-native"))))
anova_stats(anova(lm(sumcov~graze, data = cov2006%>%filter(func == "forb non-native"))))
anova_stats(anova(lm(sumcov~graze, data = cov2006%>%filter(func == "grass native"))))

cov2007 <- cov %>%
  filter(year == 2007) %>%
  filter(trt != "ungrazed unburned")
anova_stats(anova(lm(sumcov~graze, data = cov2007%>%filter(func == "forb native"))))
anova_stats(anova(lm(sumcov~graze, data = cov2007%>%filter(func == "grass non-native"))))
anova_stats(anova(lm(sumcov~graze, data = cov2007%>%filter(func == "forb non-native"))))
anova_stats(anova(lm(sumcov~graze, data = cov2007%>%filter(func == "grass native"))))

cov2008 <- cov %>%
  filter(year == 2008) %>%
  filter(trt != "ungrazed unburned")
anova_stats(anova(lm(sumcov~graze, data = cov2008%>%filter(func == "forb native"))))
anova_stats(anova(lm(sumcov~graze, data = cov2008%>%filter(func == "grass non-native"))))
anova_stats(anova(lm(sumcov~graze, data = cov2008%>%filter(func == "forb non-native"))))
anova_stats(anova(lm(sumcov~graze, data = cov2008%>%filter(func == "grass native"))))

cov2009 <- cov %>%
  filter(year == 2009)
anova_stats(anova(lm(sumcov~trt, data = cov2009%>%filter(func == "forb native"))))
anova_stats(anova(lm(sumcov~trt, data = cov2009%>%filter(func == "grass non-native"))))
anova_stats(anova(lm(sumcov~trt, data = cov2009%>%filter(func == "forb non-native"))))
anova_stats(anova(lm(sumcov~trt, data = cov2009%>%filter(func == "grass native"))))

cov2010 <- cov %>%
  filter(year == 2010)
anova_stats(anova(lm(sumcov~trt, data = cov2010%>%filter(func == "forb native"))))
anova_stats(anova(lm(sumcov~trt, data = cov2010%>%filter(func == "grass non-native"))))
anova_stats(anova(lm(sumcov~trt, data = cov2010%>%filter(func == "forb non-native"))))
anova_stats(anova(lm(sumcov~trt, data = cov2010%>%filter(func == "grass native"))))

cov2011 <- cov %>%
  filter(year == 2011)
anova_stats(anova(lm(sumcov~trt, data = cov2011%>%filter(func == "forb native"))))
anova_stats(anova(lm(sumcov~trt, data = cov2011%>%filter(func == "grass non-native"))))
anova_stats(anova(lm(sumcov~trt, data = cov2011%>%filter(func == "forb non-native"))))
anova_stats(anova(lm(sumcov~trt, data = cov2011%>%filter(func == "grass native"))))

cov2012 <- cov %>%
  filter(year == 2012)
anova_stats(anova(lm(sumcov~trt, data = cov2012%>%filter(func == "forb native"))))
anova_stats(anova(lm(sumcov~trt, data = cov2012%>%filter(func == "grass non-native"))))
anova_stats(anova(lm(sumcov~trt, data = cov2012%>%filter(func == "forb non-native"))))
anova_stats(anova(lm(sumcov~trt, data = cov2012%>%filter(func == "grass native"))))

#evenness
#pull in object "shan2" from grazing recovery
shan3 <- shan2 %>%
  mutate(trt1 = trt) %>%
  separate(trt, into=c("grazed", "burn"), sep=" ")

shan2005 <- shan3%>%
  filter(year == 2005) %>%
  filter(trt1 != "ungrazed unburned")
anova_stats(anova(lm(Shannon~grazed, data = shan2005%>%filter(func == "forb native"))))
anova_stats(anova(lm(Shannon~grazed, data = shan2005%>%filter(func == "grass non-native"))))
anova_stats(anova(lm(Shannon~grazed, data = shan2005%>%filter(func == "forb non-native"))))
anova_stats(anova(lm(Shannon~grazed, data = shan2005%>%filter(func == "grass native"))))

shan2006 <- shan3%>%
  filter(year == 2006) %>%
  filter(trt1 != "ungrazed unburned")
anova_stats(anova(lm(Shannon~grazed, data = shan2006%>%filter(func == "forb native"))))
anova_stats(anova(lm(Shannon~grazed, data = shan2006%>%filter(func == "grass non-native"))))
anova_stats(anova(lm(Shannon~grazed, data = shan2006%>%filter(func == "forb non-native"))))
anova_stats(anova(lm(Shannon~grazed, data = shan2006%>%filter(func == "grass native"))))

shan2007 <- shan3%>%
  filter(year == 2007) %>%
  filter(trt1 != "ungrazed unburned")
anova_stats(anova(lm(Shannon~grazed, data = shan2007%>%filter(func == "forb native"))))
anova_stats(anova(lm(Shannon~grazed, data = shan2007%>%filter(func == "grass non-native"))))
anova_stats(anova(lm(Shannon~grazed, data = shan2007%>%filter(func == "forb non-native"))))
anova_stats(anova(lm(Shannon~grazed, data = shan2007%>%filter(func == "grass native"))))

shan2008 <- shan3%>%
  filter(year == 2008) %>%
  filter(trt1 != "ungrazed unburned")
anova_stats(anova(lm(Shannon~grazed, data = shan2008%>%filter(func == "forb native"))))
anova_stats(anova(lm(Shannon~grazed, data = shan2008%>%filter(func == "grass non-native"))))
anova_stats(anova(lm(Shannon~grazed, data = shan2008%>%filter(func == "forb non-native"))))
anova_stats(anova(lm(Shannon~grazed, data = shan2008%>%filter(func == "grass native"))))

shan2009 <- shan3%>%
  filter(year == 2009)
anova_stats(anova(lm(Shannon~trt1, data = shan2009%>%filter(func == "forb native"))))
anova_stats(anova(lm(Shannon~trt1, data = shan2009%>%filter(func == "grass non-native"))))
anova_stats(anova(lm(Shannon~trt1, data = shan2009%>%filter(func == "forb non-native"))))
anova_stats(anova(lm(Shannon~trt1, data = shan2009%>%filter(func == "grass native"))))

shan2010 <- shan3%>%
  filter(year == 2010)
anova_stats(anova(lm(Shannon~trt1, data = shan2010%>%filter(func == "forb native"))))
anova_stats(anova(lm(Shannon~trt1, data = shan2010%>%filter(func == "grass non-native"))))
anova_stats(anova(lm(Shannon~trt1, data = shan2010%>%filter(func == "forb non-native"))))
anova_stats(anova(lm(Shannon~trt1, data = shan2010%>%filter(func == "grass native"))))

shan2011 <- shan3%>%
  filter(year == 2011)
anova_stats(anova(lm(Shannon~trt1, data = shan2011%>%filter(func == "forb native"))))
anova_stats(anova(lm(Shannon~trt1, data = shan2011%>%filter(func == "grass non-native"))))
anova_stats(anova(lm(Shannon~trt1, data = shan2011%>%filter(func == "forb non-native"))))
anova_stats(anova(lm(Shannon~trt1, data = shan2011%>%filter(func == "grass native"))))

shan2012 <- shan3%>%
  filter(year == 2012)
anova_stats(anova(lm(Shannon~trt1, data = shan2012%>%filter(func == "forb native"))))
anova_stats(anova(lm(Shannon~trt1, data = shan2012%>%filter(func == "grass non-native"))))
anova_stats(anova(lm(Shannon~trt1, data = shan2012%>%filter(func == "forb non-native"))))
anova_stats(anova(lm(Shannon~trt1, data = shan2012%>%filter(func == "grass native"))))

#litter
#pull in object "joindat" from grazingrecovery
lit <- joindat %>%
  mutate(trt=paste(graze, burn)) %>%
  ungroup() %>%
  select(quadratNew, year, trt, transect, burn, graze, litter)

lit2005 <- lit%>%
  filter(year == 2005) %>%
  filter(trt != "ungrazed unburned")
anova_stats(anova(lm(litter~graze, data = lit2005)))

lit2006 <- lit%>%
  filter(year == 2006) %>%
  filter(trt != "ungrazed unburned")
anova_stats(anova(lm(litter~graze, data = lit2006)))

lit2007 <- lit%>%
  filter(year == 2007) %>%
  filter(trt != "ungrazed unburned")
anova_stats(anova(lm(litter~graze, data = lit2007)))

lit2008 <- lit%>%
  filter(year == 2008) %>%
  filter(trt != "ungrazed unburned")
anova_stats(anova(lm(litter~graze, data = lit2008)))

lit2009 <- lit%>%
  filter(year == 2009) 
anova_stats(anova(lm(litter~trt, data = lit2009)))

lit2010 <- lit%>%
  filter(year == 2010) 
anova_stats(anova(lm(litter~trt, data = lit2010)))

lit2011 <- lit%>%
  filter(year == 2011) 
anova_stats(anova(lm(litter~trt, data = lit2011)))

lit2012 <- lit%>%
  filter(year == 2012) 
anova_stats(anova(lm(litter~trt, data = lit2012)))

########
#one-way ANOVA analyze two blocks of years: 2005-2008 and 2008-2012 
########
rich3 <- rich2 %>%
  filter(trt != "ungrazed unburned") %>%
  filter(year%in%c(2005, 2006, 2007, 2008))

anova_stats(anova(lm(richness~graze, data = rich3%>%filter(func == "forb native"))))
anova_stats(anova(lm(richness~graze, data = rich3%>%filter(func == "grass non-native"))))
anova_stats(anova(lm(richness~graze, data = rich3%>%filter(func == "forb non-native"))))
anova_stats(anova(lm(richness~graze, data = rich3%>%filter(func == "grass native"))))

rich4 <- rich2 %>%
  filter(year%in%c(2008:2012))

anova_stats(anova(lm(richness~trt, data = rich4%>%filter(func == "forb native"))))
anova_stats(anova(lm(richness~trt, data = rich4%>%filter(func == "grass non-native"))))
anova_stats(anova(lm(richness~trt, data = rich4%>%filter(func == "forb non-native"))))
anova_stats(anova(lm(richness~trt, data = rich4%>%filter(func == "grass native"))))

cov_pre <- cov %>%
  filter(trt != "ungrazed unburned") %>%
  filter(year%in%c(2005:2008))

anova_stats(anova(lm(sumcov~graze, data = cov_pre%>%filter(func == "forb native"))))
anova_stats(anova(lm(sumcov~graze, data = cov_pre%>%filter(func == "grass non-native"))))
anova_stats(anova(lm(sumcov~graze, data = cov_pre%>%filter(func == "forb non-native"))))
anova_stats(anova(lm(sumcov~graze, data = cov_pre%>%filter(func == "grass native"))))

cov_post <- cov %>%
  filter(year%in%c(2008:2012))
anova_stats(anova(lm(sumcov~trt, data = cov_post%>%filter(func == "forb native"))))
anova_stats(anova(lm(sumcov~trt, data = cov_post%>%filter(func == "grass non-native"))))
anova_stats(anova(lm(sumcov~trt, data = cov_post%>%filter(func == "forb non-native"))))
anova_stats(anova(lm(sumcov~trt, data = cov_post%>%filter(func == "grass native"))))

shan_pre <- shan3 %>%
  filter(trt1 != "ungrazed unburned") %>%
  filter(year%in%c(2005:2008))
anova_stats(anova(lm(Shannon~grazed, data = shan_pre%>%filter(func == "forb native"))))
anova_stats(anova(lm(Shannon~grazed, data = shan_pre%>%filter(func == "grass non-native"))))
anova_stats(anova(lm(Shannon~grazed, data = shan_pre%>%filter(func == "forb non-native"))))
anova_stats(anova(lm(Shannon~grazed, data = shan_pre%>%filter(func == "grass native"))))

shan_post <- shan3 %>%
  filter(year%in%c(2008:2012))
anova_stats(anova(lm(Shannon~trt1, data = shan_post%>%filter(func == "forb native"))))
anova_stats(anova(lm(Shannon~trt1, data = shan_post%>%filter(func == "grass non-native"))))
anova_stats(anova(lm(Shannon~trt1, data = shan_post%>%filter(func == "forb non-native"))))
anova_stats(anova(lm(Shannon~trt1, data = shan_post%>%filter(func == "grass native"))))

lit_pre <- lit %>%
  filter(trt != "ungrazed unburned") %>%
  filter(year%in%c(2005:2008))
anova_stats(anova(lm(litter~graze, data = lit_pre)))

lit_post <- lit %>%
  filter(year%in%c(2008:2012))
anova_stats(anova(lm(litter~trt, data = lit_post)))

########
#linear mixed models
#grazed and burned as fixed effects, year as random
#nest quadrats within transects
########

#richness

fit2005fn<- lme(richness~graze, random = ~1|transect/quadratNew, data = rich2005%>%filter(func == "forb native"))
summary(fit2005fn)
fit2005ge<- lme(richness~graze, random = ~1|transect/quadratNew, data = rich2005%>%filter(func == "grass non-native"))
summary(fit2005ge)

fitpostfn <- lme(richness~graze+burn, random = ~1|transect/quadratNew, data = rich4%>%filter(func == "forb native"))
summary(fitpostfn)


# all years??
allshan<-lme(Shannon~grazed+burn, random = list(~1|transectNew, ~1|year), data = shan3%>%filter(func == "forb native"))
summary(allshan)

allshan<-lme(Shannon~trt1, random = list(~1|transectNew, ~1|year), data = shan3%>%filter(func == "forb native"))
summary(allshan)

##example from onoline
lme(y ~ p_gender*t_gender + part_gen, data=grdata,
    + random = list(~ p_gender | therapist, ~ 1 | group))

