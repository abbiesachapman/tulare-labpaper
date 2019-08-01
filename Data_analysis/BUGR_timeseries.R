## Set your datpath!! (file in Data_cleaning)

library(tidyverse)
library(readr)
library(ggplot2)

#load data
dat <- read_csv(paste(datpath_clean, "/bugrdat.csv", sep=""))
SC <- read_csv(paste(datpath_clean, "/SpeciesCodes.csv", sep=""))
clim <- read_csv(paste(datpath_clean, "/san_jose_clim.csv", sep =""), skip = 9) %>%
  mutate(stand_ppt = scale(PRCP, center = T, scale = T)) %>% #standardize to z-scores
  mutate(stand_temp = scale(TAVG, center = T, scale = T)) #standardize to z-scores

View(clim)

#plot timeseries of richness 
richness <- dat %>%
  filter(cover != 0, spname != c("Unknown", "Moss")) %>%
  group_by(year, quadratNew) %>%
  summarize(richness = length(unique(spname))) %>%
  group_by(year) %>%
  summarize(mean_rich = mean(richness))

ggplot(richness, aes(year, mean_rich)) +
  geom_point() +
  geom_line()

#join data and species key
tog<-left_join(dat, SC)
View(tog)

#plot timeseries of cover by thermal
functog<-tog%>%
  group_by(quadratNew, year, status, func, thermal)%>%
  summarize(sumcov=sum(cover))%>%
  filter(!is.na(status), !is.na(func))
functogagg<-functog%>%
  group_by(year, status, func)%>%
  summarize(meancov=mean(sumcov))

View(functoagg)
View(functog)

thermal_trend <- functog %>%
  group_by(year, status, func, thermal) %>%
  summarize(meancov = mean(sumcov)) %>%
  filter(thermal != "?")

ggplot(functog, aes(as.factor(year), sumcov))+geom_boxplot()+facet_grid(status~func)
ggplot(functogagg, aes((year), meancov))+geom_line(aes(color=interaction(status, func)))+geom_point(aes(color=interaction(status, func)))

#cover by thermal transposed on annual precip
ggplot(thermal_trend, aes((year), meancov)) + facet_grid(status~func) +
  geom_bar(data = clim, aes(x = DATE, y = PRCP*3), stat = "identity", fill = "lightgrey") +
  geom_line(aes(color= thermal))+ geom_point(aes(color = thermal))  +
  scale_y_continuous(sec.axis = sec_axis(~./3, name = "Annual Precipitation in inches"))

#timeseries of ppt
ggplot(clim, aes(DATE,PRCP)) + geom_line() + geom_point()

#cover by thermal transpoesd on annual mean temp
ggplot(thermal_trend, aes((year), meancov)) + facet_grid(status~func) +
  geom_bar(data = clim, aes(x = DATE, y = stand_temp*30), stat = "identity", fill = "lightgrey") +
  geom_line(aes(color= thermal))+ geom_point(aes(color = thermal))  +
  scale_y_continuous(sec.axis = sec_axis(~./30, name = "Annual Mean Temp Deviation (z-scores)"))

#join clim and thermal_trend
clim1 <- clim %>%
  rename(year = DATE)
joined_dat <- left_join(thermal_trend, clim1)

#cover by temp
ggplot(joined_dat, aes((PRCP), meancov)) + facet_grid(status~func) +
   geom_point(aes(color = thermal)) + geom_smooth(aes(color = thermal), method = "lm", se = F)

