library(tidyverse)
library(readr)
library(sjstats) #standardized effect size cohen's f
library(nlme) #linear mixed models

#load updated master data, "alldat", in grazing_recovery.R

########
#one-way ANOVA analyze within each year 
#richness, cover, evenness, litter as dependent variables
#focus on two functional groups: native forbs and non-native grasses
########

#Richness
#load "rich" from grazing_recovery.R

TukeyHSD(aov(richness~trt, data = rich%>%filter(year == 2005, func == "forb native")))
TukeyHSD(aov(richness~trt, data = rich%>%filter(year == 2006, func == "forb native")))
TukeyHSD(aov(richness~trt, data = rich%>%filter(year == 2007, func == "forb native")))
TukeyHSD(aov(richness~trt, data = rich%>%filter(year == 2008, func == "forb native")))
TukeyHSD(aov(richness~trt, data = rich%>%filter(year == 2009, func == "forb native")))
TukeyHSD(aov(richness~trt, data = rich%>%filter(year == 2010, func == "forb native")))
TukeyHSD(aov(richness~trt, data = rich%>%filter(year == 2011, func == "forb native")))
TukeyHSD(aov(richness~trt, data = rich%>%filter(year == 2012, func == "forb native")))

TukeyHSD(aov(richness~trt, data = rich%>%filter(year == 2005, func == "grass non-native")))
TukeyHSD(aov(richness~trt, data = rich%>%filter(year == 2006, func == "grass non-native")))
TukeyHSD(aov(richness~trt, data = rich%>%filter(year == 2007, func == "grass non-native")))
TukeyHSD(aov(richness~trt, data = rich%>%filter(year == 2008, func == "grass non-native")))
TukeyHSD(aov(richness~trt, data = rich%>%filter(year == 2009, func == "grass non-native")))
TukeyHSD(aov(richness~trt, data = rich%>%filter(year == 2010, func == "grass non-native")))
TukeyHSD(aov(richness~trt, data = rich%>%filter(year == 2011, func == "grass non-native")))
TukeyHSD(aov(richness~trt, data = rich%>%filter(year == 2012, func == "grass non-native")))

rich2005 <- rich %>%
  filter(trt != "ungrazed unburned") %>%
  filter(year == 2005)
anova_stats(anova(lm(richness~graze, data = rich2005%>%filter(func == "forb native"))))
anova_stats(anova(lm(richness~graze, data = rich2005%>%filter(func == "grass non-native"))))
anova_stats(anova(lm(richness~graze, data = rich2005%>%filter(func == "forb non-native"))))
anova_stats(anova(lm(richness~graze, data = rich2005%>%filter(func == "grass native"))))

rich2006 <- rich %>%
  filter(trt != "ungrazed unburned") %>%
  filter(year == 2006)
anova_stats(anova(lm(richness~graze, data = rich2006%>%filter(func == "forb native"))))
anova_stats(anova(lm(richness~graze, data = rich2006%>%filter(func == "grass non-native"))))
anova_stats(anova(lm(richness~graze, data = rich2006%>%filter(func == "forb non-native"))))
anova_stats(anova(lm(richness~graze, data = rich2006%>%filter(func == "grass native"))))

rich2007 <- rich %>%
  filter(trt != "ungrazed unburned") %>%
  filter(year == 2007)
anova_stats(anova(lm(richness~graze, data = rich2007%>%filter(func == "forb native"))))
anova_stats(anova(lm(richness~graze, data = rich2007%>%filter(func == "grass non-native"))))
anova_stats(anova(lm(richness~graze, data = rich2007%>%filter(func == "forb non-native"))))
anova_stats(anova(lm(richness~graze, data = rich2007%>%filter(func == "grass native"))))

rich2008 <- rich %>%
  filter(trt != "ungrazed unburned") %>%
  filter(year == 2008)
anova_stats(anova(lm(richness~graze, data = rich2008%>%filter(func == "forb native"))))
anova_stats(anova(lm(richness~graze, data = rich2008%>%filter(func == "grass non-native"))))
anova_stats(anova(lm(richness~graze, data = rich2008%>%filter(func == "forb non-native"))))
anova_stats(anova(lm(richness~graze, data = rich2008%>%filter(func == "grass native"))))

rich2009 <- rich %>%
  filter(year == 2009)
anova_stats(anova(lm(richness~trt, data = rich2009%>%filter(func == "forb native"))))
anova_stats(anova(lm(richness~trt, data = rich2009%>%filter(func == "grass non-native"))))
anova_stats(anova(lm(richness~trt, data = rich2009%>%filter(func == "forb non-native"))))
anova_stats(anova(lm(richness~trt, data = rich2009%>%filter(func == "grass native"))))

rich2010 <- rich %>%
  filter(year == 2010)
anova_stats(anova(lm(richness~trt, data = rich2010%>%filter(func == "forb native"))))
anova_stats(anova(lm(richness~trt, data = rich2010%>%filter(func == "grass non-native"))))
anova_stats(anova(lm(richness~trt, data = rich2010%>%filter(func == "forb non-native"))))
anova_stats(anova(lm(richness~trt, data = rich2010%>%filter(func == "grass native"))))

rich2011 <- rich %>%
  filter(year == 2011)
anova_stats(anova(lm(richness~trt, data = rich2011%>%filter(func == "forb native"))))
anova_stats(anova(lm(richness~trt, data = rich2011%>%filter(func == "grass non-native"))))
anova_stats(anova(lm(richness~trt, data = rich2011%>%filter(func == "forb non-native"))))
anova_stats(anova(lm(richness~trt, data = rich2011%>%filter(func == "grass native"))))

rich2012 <- rich %>%
  filter(year == 2012)
anova_stats(anova(lm(richness~trt, data = rich2012%>%filter(func == "forb native"))))
anova_stats(anova(lm(richness~trt, data = rich2012%>%filter(func == "grass non-native"))))
anova_stats(anova(lm(richness~trt, data = rich2012%>%filter(func == "forb non-native"))))
anova_stats(anova(lm(richness~trt, data = rich2012%>%filter(func == "grass native"))))

#cover
#load "cov" from grazing_recovery.R

TukeyHSD(aov(relcov~trt, data = cov%>%filter(year == 2005, func == "forb native")))
TukeyHSD(aov(relcov~trt, data = cov%>%filter(year == 2006, func == "forb native")))
TukeyHSD(aov(relcov~trt, data = cov%>%filter(year == 2007, func == "forb native")))
TukeyHSD(aov(relcov~trt, data = cov%>%filter(year == 2008, func == "forb native")))
TukeyHSD(aov(relcov~trt, data = cov%>%filter(year == 2009, func == "forb native")))
TukeyHSD(aov(relcov~trt, data = cov%>%filter(year == 2010, func == "forb native")))
TukeyHSD(aov(relcov~trt, data = cov%>%filter(year == 2011, func == "forb native")))
TukeyHSD(aov(relcov~trt, data = cov%>%filter(year == 2012, func == "forb native")))

TukeyHSD(aov(relcov~trt, data = cov%>%filter(year == 2005, func == "grass non-native")))
TukeyHSD(aov(relcov~trt, data = cov%>%filter(year == 2006, func == "grass non-native")))
TukeyHSD(aov(relcov~trt, data = cov%>%filter(year == 2007, func == "grass non-native")))
TukeyHSD(aov(relcov~trt, data = cov%>%filter(year == 2008, func == "grass non-native")))
TukeyHSD(aov(relcov~trt, data = cov%>%filter(year == 2009, func == "grass non-native")))
TukeyHSD(aov(relcov~trt, data = cov%>%filter(year == 2010, func == "grass non-native")))
TukeyHSD(aov(relcov~trt, data = cov%>%filter(year == 2011, func == "grass non-native")))
TukeyHSD(aov(relcov~trt, data = cov%>%filter(year == 2012, func == "grass non-native")))
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
#load object "shan2" from grazing_recovery.R
shan3 <- shan2 %>%
  mutate(trt1 = trt) %>%
  separate(trt, into=c("grazed", "burn"), sep=" ")

TukeyHSD(aov(Shannon~trt, data = shan2%>%filter(year == 2005, func == "forb native")))
TukeyHSD(aov(Shannon~trt, data = shan2%>%filter(year == 2006, func == "forb native")))
TukeyHSD(aov(Shannon~trt, data = shan2%>%filter(year == 2007, func == "forb native")))
TukeyHSD(aov(Shannon~trt, data = shan2%>%filter(year == 2008, func == "forb native")))
TukeyHSD(aov(Shannon~trt, data = shan2%>%filter(year == 2009, func == "forb native")))
TukeyHSD(aov(Shannon~trt, data = shan2%>%filter(year == 2010, func == "forb native")))
TukeyHSD(aov(Shannon~trt, data = shan2%>%filter(year == 2011, func == "forb native")))
TukeyHSD(aov(Shannon~trt, data = shan2%>%filter(year == 2012, func == "forb native")))

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

TukeyHSD(aov(litter~trt, data = lit%>%filter(year == 2005)))
TukeyHSD(aov(litter~trt, data = lit%>%filter(year == 2006)))
TukeyHSD(aov(litter~trt, data = lit%>%filter(year == 2007)))
TukeyHSD(aov(litter~trt, data = lit%>%filter(year == 2008)))
TukeyHSD(aov(litter~trt, data = lit%>%filter(year == 2009)))
TukeyHSD(aov(litter~trt, data = lit%>%filter(year == 2010)))
TukeyHSD(aov(litter~trt, data = lit%>%filter(year == 2011)))
TukeyHSD(aov(litter~trt, data = lit%>%filter(year == 2012)))

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
rich3 <- rich %>%
  filter(trt != "ungrazed unburned") %>%
  filter(year%in%c(2005, 2006, 2007, 2008))

anova_stats(anova(lm(richness~graze, data = rich3%>%filter(func == "forb native"))))
anova_stats(anova(lm(richness~graze, data = rich3%>%filter(func == "grass non-native"))))
anova_stats(anova(lm(richness~graze, data = rich3%>%filter(func == "forb non-native"))))
anova_stats(anova(lm(richness~graze, data = rich3%>%filter(func == "grass native"))))

rich4 <- rich %>%
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
#grazed and burned as fixed effects, quadrats nested within transects as random
########

# richness by year
richrich<-rich%>%
  separate(quadratNew, into=c("transect2", "quadrat"), sep="-") 

rich_NF2005<-lme(richness~trt, random = ~1|transect2/quadrat, data = subset(richrich, func=="forb native"&year==2005))
rich_NF2006<-lme(richness~trt, random = ~1|transect2/quadrat, data = subset(richrich, func=="forb native"&year==2006))
rich_NF2007<-lme(richness~trt, random = ~1|transect2/quadrat, data = subset(richrich, func=="forb native"&year==2007))
rich_NF2008<-lme(richness~trt, random = ~1|transect2/quadrat, data = subset(richrich, func=="forb native"&year==2008))
rich_NF2009<-lme(richness~trt, random = ~1|transect2/quadrat, data = subset(richrich, func=="forb native"&year==2009))
rich_NF2010<-lme(richness~trt, random = ~1|transect2/quadrat, data = subset(richrich, func=="forb native"&year==2010))
rich_NF2011<-lme(richness~trt, random = ~1|transect2/quadrat, data = subset(richrich, func=="forb native"&year==2011))
rich_NF2012<-lme(richness~trt, random = ~1|transect2/quadrat, data = subset(richrich, func=="forb native"&year==2012))

rich_IG2005<-lme(richness~trt, random = ~1|transect2/quadrat, data = subset(richrich, func=="grass non-native"&year==2005))
rich_IG2006<-lme(richness~trt, random = ~1|transect2/quadrat, data = subset(richrich, func=="grass non-native"&year==2006))
rich_IG2007<-lme(richness~trt, random = ~1|transect2/quadrat, data = subset(richrich, func=="grass non-native"&year==2007))
rich_IG2008<-lme(richness~trt, random = ~1|transect2/quadrat, data = subset(richrich, func=="grass non-native"&year==2008))
rich_IG2009<-lme(richness~trt, random = ~1|transect2/quadrat, data = subset(richrich, func=="grass non-native"&year==2009))
rich_IG2010<-lme(richness~trt, random = ~1|transect2/quadrat, data = subset(richrich, func=="grass non-native"&year==2010))
rich_IG2011<-lme(richness~trt, random = ~1|transect2/quadrat, data = subset(richrich, func=="grass non-native"&year==2011))
rich_IG2012<-lme(richness~trt, random = ~1|transect2/quadrat, data = subset(richrich, func=="grass non-native"&year==2012))

summary(rich_NF2005)
summary(rich_NF2006)
summary(rich_NF2007)
summary(rich_NF2008)
summary(rich_NF2009)
summary(rich_NF2010)
summary(rich_NF2011)
summary(rich_NF2012)

summary(rich_IG2005)
summary(rich_IG2006)
summary(rich_IG2007)
summary(rich_IG2008)
summary(rich_IG2009)
summary(rich_IG2010)
summary(rich_IG2011)
summary(rich_IG2012)

# cover by year
covcov<-cov%>%
  separate(quadratNew, into=c("transect2", "quadrat"), sep="-") 

covcov2<-covcov%>%
  ungroup()%>%
  mutate(trt=as.factor(trt))
covcov2$trt <- factor(covcov2$trt, levels = c("ungrazed burned", "ungrazed unburned", "grazed burned"))

cov_NF2005<-lme(relcov~trt, random = ~1|transect2/quadrat, data = subset(covcov, func=="forb native"&year==2005))
cov_NF2005.2<-lme(relcov~trt, random = ~1|transect2/quadrat, data = subset(covcov2, func=="forb native"&year==2005))
cov_NF2006<-lme(relcov~trt, random = ~1|transect2/quadrat, data = subset(covcov, func=="forb native"&year==2006))
cov_NF2006.2<-lme(relcov~trt, random = ~1|transect2/quadrat, data = subset(covcov2, func=="forb native"&year==2006))
cov_NF2007<-lme(relcov~trt, random = ~1|transect2/quadrat, data = subset(covcov, func=="forb native"&year==2007))
cov_NF2007.2<-lme(relcov~trt, random = ~1|transect2/quadrat, data = subset(covcov2, func=="forb native"&year==2007))
cov_NF2008<-lme(relcov~trt, random = ~1|transect2/quadrat, data = subset(covcov, func=="forb native"&year==2008))
cov_NF2008.2<-lme(relcov~trt, random = ~1|transect2/quadrat, data = subset(covcov2, func=="forb native"&year==2008))
cov_NF2009<-lme(relcov~trt, random = ~1|transect2/quadrat, data = subset(covcov, func=="forb native"&year==2009))
cov_NF2009.2<-lme(relcov~trt, random = ~1|transect2/quadrat, data = subset(covcov2, func=="forb native"&year==2009))
cov_NF2010<-lme(relcov~trt, random = ~1|transect2/quadrat, data = subset(covcov, func=="forb native"&year==2010))
cov_NF2010.2<-lme(relcov~trt, random = ~1|transect2/quadrat, data = subset(covcov2, func=="forb native"&year==2010))
cov_NF2011<-lme(relcov~trt, random = ~1|transect2/quadrat, data = subset(covcov, func=="forb native"&year==2011))
cov_NF2011.2<-lme(relcov~trt, random = ~1|transect2/quadrat, data = subset(covcov2, func=="forb native"&year==2011))
cov_NF2012<-lme(relcov~trt, random = ~1|transect2/quadrat, data = subset(covcov, func=="forb native"&year==2012))
cov_NF2012.2<-lme(relcov~trt, random = ~1|transect2/quadrat, data = subset(covcov2, func=="forb native"&year==2012))

cov_IG2005<-lme(relcov~trt, random = ~1|transect2/quadrat, data = subset(covcov, func=="grass non-native"&year==2005))
cov_IG2005.2<-lme(relcov~trt, random = ~1|transect2/quadrat, data = subset(covcov2, func=="grass non-native"&year==2005))
cov_IG2006<-lme(relcov~trt, random = ~1|transect2/quadrat, data = subset(covcov, func=="grass non-native"&year==2006))
cov_IG2006.2<-lme(relcov~trt, random = ~1|transect2/quadrat, data = subset(covcov2, func=="grass non-native"&year==2006))
cov_IG2007<-lme(relcov~trt, random = ~1|transect2/quadrat, data = subset(covcov, func=="grass non-native"&year==2007))
cov_IG2007.2<-lme(relcov~trt, random = ~1|transect2/quadrat, data = subset(covcov2, func=="grass non-native"&year==2007))
cov_IG2008<-lme(relcov~trt, random = ~1|transect2/quadrat, data = subset(covcov, func=="grass non-native"&year==2008))
cov_IG2008.2<-lme(relcov~trt, random = ~1|transect2/quadrat, data = subset(covcov2, func=="grass non-native"&year==2008))
cov_IG2009<-lme(relcov~trt, random = ~1|transect2/quadrat, data = subset(covcov, func=="grass non-native"&year==2009))
cov_IG2009.2<-lme(relcov~trt, random = ~1|transect2/quadrat, data = subset(covcov2, func=="grass non-native"&year==2009))
cov_IG2010<-lme(relcov~trt, random = ~1|transect2/quadrat, data = subset(covcov, func=="grass non-native"&year==2010))
cov_IG2010.2<-lme(relcov~trt, random = ~1|transect2/quadrat, data = subset(covcov2, func=="grass non-native"&year==2010))
cov_IG2011<-lme(relcov~trt, random = ~1|transect2/quadrat, data = subset(covcov, func=="grass non-native"&year==2011))
cov_IG2011.2<-lme(relcov~trt, random = ~1|transect2/quadrat, data = subset(covcov2, func=="grass non-native"&year==2011))
cov_IG2012<-lme(relcov~trt, random = ~1|transect2/quadrat, data = subset(covcov, func=="grass non-native"&year==2012))
cov_IG2012.2<-lme(relcov~trt, random = ~1|transect2/quadrat, data = subset(covcov2, func=="grass non-native"&year==2012))

summary(cov_NF2005)
summary(cov_NF2005.2)
summary(cov_NF2006)
summary(cov_NF2006.2)
summary(cov_NF2007)
summary(cov_NF2007.2)
summary(cov_NF2008)
summary(cov_NF2008.2)
summary(cov_NF2009)
summary(cov_NF2009.2)
summary(cov_NF2010)
summary(cov_NF2010.2)
summary(cov_NF2011)
summary(cov_NF2011.2)
summary(cov_NF2012)
summary(cov_NF2012.2)

summary(cov_IG2005)
summary(cov_IG2005.2)
summary(cov_IG2006)
summary(cov_IG2006.2)
summary(cov_IG2007)
summary(cov_IG2007.2)
summary(cov_IG2008)
summary(cov_IG2008.2)
summary(cov_IG2009)
summary(cov_IG2009.2)
summary(cov_IG2010)
summary(cov_IG2010.2)
summary(cov_IG2011)
summary(cov_IG2011.2)
summary(cov_IG2012)
summary(cov_IG2012.2)

##example from onoline
lme(y ~ p_gender*t_gender + part_gen, data=grdata,
    + random = list(~ p_gender | therapist, ~ 1 | group))

##plot lme stats

lmestats<-read_csv(paste(datpath_clean, "/Grazing_recovery_stats_lme.csv", sep=""))%>%
  select(-Test)%>%
  rename("p"="p-value", "p1"="p-value_1", "t"="t-value", "t1"="t-value_1", "ugub"="ungrazed unburned", "ugb"="ungrazed burned")

lmestats<-lmestats%>%
  mutate(p=ifelse(p<.05, .05, ifelse(p<.1, .1, NA)))%>%
  mutate(p1=ifelse(p1<.05, .05, ifelse(p1<.1, .1, NA)))
  
lme1<-ggplot(subset(lmestats, func=="native forb"&response=="cover"), aes(year, ugb*100)) +geom_line(color="blue") +geom_point(color="blue") +
  geom_line(aes(year, ugub*100), color="red") +geom_point(aes(year, ugub*100), color="red") +
  geom_point(aes(year-.1, y=-p*100),shape=8, color="blue")+geom_point(aes(year+.1,y=-p1*100),shape=8, color="red")+
  geom_hline(yintercept=0) +
  annotate("text", x= 2004.75, y = .5, label = "fire", size = 3) +
  annotate("text", x= 2009, y = .5, label = "grazing", size = 3)+
  geom_vline(xintercept=2008.5, color = "grey66", lty =2)+
  geom_vline(xintercept=2004.5, color="grey66", lty = 2)+
  xlab("")+ylab("% Cover Effect Size vs Burned-Grazed")

lme2<-ggplot(subset(lmestats, func=="non-native grass"&response=="cover"), aes(year, ugb*100)) +geom_line(color="blue") +geom_point(color="blue") +
  geom_line(aes(year, ugub*100), color="red") +geom_point(aes(year, ugub*100), color="red") +
  geom_point(aes(year-.1, y=p*100),shape=8, color="blue")+geom_point(aes(year+.1,y=p1*100),shape=8, color="red")+
  geom_hline(yintercept=0) +
  annotate("text", x= 2004.75, y = .5, label = "fire", size = 3) +
  annotate("text", x= 2009, y = .5, label = "grazing", size = 3)+
  geom_vline(xintercept=2008.5, color = "grey66", lty =2)+
  geom_vline(xintercept=2004.5, color="grey66", lty = 2)+
  xlab("")+ylab("")


lme3<-ggplot(subset(lmestats, func=="native forb"&response=="richness"), aes(year, ugb)) +geom_line(color="blue") +geom_point(color="blue") +
  geom_line(aes(year, ugub), color="red") +geom_point(aes(year, ugub), color="red") +
  geom_point(aes(year-.1, y=p*10),shape=8, color="blue")+geom_point(aes(year+.1,y=p1*10),shape=8, color="red")+
  geom_hline(yintercept=0) +
  annotate("text", x= 2004.75, y = .5, label = "fire", size = 3) +
  annotate("text", x= 2009, y = .5, label = "grazing", size = 3)+
  geom_vline(xintercept=2008.5, color = "grey66", lty =2)+
  geom_vline(xintercept=2004.5, color="grey66", lty = 2)+
  xlab("")+ylab("Richness Effect Size vs Burned-Grazed")


lme4<-ggplot(subset(lmestats, func=="non-native grass"&response=="richness"), aes(year, ugb)) +geom_line(color="blue") +geom_point(color="blue") +
  geom_line(aes(year, ugub), color="red") +geom_point(aes(year, ugub), color="red") +
  geom_point(aes(year-.1, y=p),shape=8, color="blue")+geom_point(aes(year+.1,y=p1),shape=8, color="red")+
  geom_hline(yintercept=0) +
  annotate("text", x= 2004.75, y = .5, label = "fire", size = 3) +
  annotate("text", x= 2009, y = .5, label = "grazing", size = 3)+
  geom_vline(xintercept=2008.5, color = "grey66", lty =2)+
  geom_vline(xintercept=2004.5, color="grey66", lty = 2)+
  xlab("")+ylab("")


ggarrange(lme1, lme2, lme3, lme4,  ncol = 2, nrow = 2, 
          labels = c("a) Native forb cover", "b) Non-native grass cover",
                     "c) Native forb richness", "d) Non-native grass richness"),
          common.legend = TRUE, legend = "right", 
          font.label = list(size = 10),
          hjust = c(-0.5, -0.35, -0.5, -0.35))
  


                                                     