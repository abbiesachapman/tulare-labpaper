#run after BUGR_timeseries.R
# Objective: compare grazed and ungrazed plots as grazing was reintroduced
# 2006-2008 = ungrazed plots are ungrazed, 2009-2012 = ungrazed plots have cattle reintroduced
# all years = grazed plots are grazed
library(ggplot2); theme_set(theme_classic())

alldat<-read_csv(paste(datpath_clean, "/alldatsptrt.csv", sep="")) %>%
  select(-1)%>%
  filter(transect%in%c("THBUGM1", "THBUGM2", "THM1", "THM2", "THM3", "THM4", "THUBUGM1", "THUBUGM2"))%>%
  filter(!quadratNew%in%c("THM1-1", "THM3-3", "THM1-10"))%>%
  filter(thermal=="moderate")%>%
  group_by(year, spname, spcode, quadratNew, status, type, transect, burn, graze)%>%
  summarize(cover=sum(cover))

#plot timeseries of richness 
rich <- alldat %>%
  filter(type!="NA", status!="NA")%>%
  filter(spname !="Unknown", spname!="Moss") %>%
  mutate(func=paste(type, status))%>%
  mutate(trt=paste(graze, burn))%>%
  mutate(present=ifelse(cover>0, 1, 0))%>%
  group_by(year, quadratNew, func, trt, type, status, burn, graze, transect, present)%>%
  summarize(richness = sum(present))%>%
  group_by(year, quadratNew, func, trt, type, status, burn, graze, transect)%>%
  summarize(richness=sum(richness))

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
#shan<-alldat%>%
#  filter(type!="NA", status!="NA")%>%
#  mutate(func=paste(type, status))%>%
#  mutate(trt=paste(graze, burn))%>%
#  filter(cover != 0, spname !="Unknown", spname!="Moss") %>%
#  mutate(alltrt=paste(quadratNew, trt, func, sep="_"))
#simp.g<-community_diversity(shan, time.var = "year", abundance.var="cover", replicate.var="alltrt", metric = c("InverseSimpson"))
#shandiv.g<- community_diversity(shan, time.var = "year", abundance.var="cover", replicate.var="alltrt", metric = c("Shannon")) 

#shan2<-left_join(simp.g, shandiv.g)%>%
#  separate(alltrt, into=c("transectNew", "trt", "func"), sep="_")%>%
#  mutate(prepost=ifelse(year<2009, "pre", "post"))

#shan1<-shan2%>%
#  group_by(year, trt, func, prepost)%>%
#  summarize(meanShan=mean(Shannon), seShan=calcSE(Shannon), meanSimp=mean(InverseSimpson))%>%
#  filter(!is.na(trt), !is.na(func), func!=("NA"), trt!="NA")

#ggplot(shan1, aes(year, meanShan)) +
#  geom_line(aes(color=as.factor(trt)))+
#  geom_point(aes(color=as.factor(trt)))+
#  facet_wrap(~func) +
#  geom_errorbar(aes(ymin=meanShan-seShan, ymax=meanShan+seShan, color=as.factor(trt)), width=.2)+
#  geom_vline(xintercept=2008.5)+geom_vline(xintercept=2004.5, color="red")+ylab("Shannon Diversity") +
#  labs(x = "Year", y = "Mean Shannon Diversity", color = "Treatment")

#plot timeseries of cover
cov<-alldat%>%
  filter(type!="NA", status!="NA")%>%
  mutate(func=paste(type, status))%>%
  mutate(trt=paste(graze, burn))%>%
  filter(spname !="Unknown", spname!="Moss") %>%
  group_by(year, quadratNew, trt, func, burn, graze, transect)%>%
  summarize(sumcov=sum(cover))%>%
  filter(!is.na(trt), !is.na(func), !is.na(sumcov)) %>%
  group_by(year, quadratNew, trt)%>%
  mutate(totcov=sum(sumcov))%>%
  mutate(relcov=sumcov/totcov)

cov1<-cov%>%
  group_by(year, trt, func)%>%
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

##raw data
#ggplot(cov, aes((year), relcov))+
#  geom_point(aes(color=trt))+
#  facet_wrap(~func) +geom_vline(xintercept=2008.5)+geom_vline(xintercept=2004.5, color="red") +
#  labs(x = "Year", y = "Relative Cover", color = "Treatment")

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
  filter(year != "2005") %>%
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
#Publication figure
##########
library(ggpubr)
f1 <- ggplot(rich1%>%filter(func == "forb native", year%in%c(2005:2012)), aes(year, mean_rich)) +
        geom_line(aes(color=as.factor(trt))) +
        geom_point(aes(fill=as.factor(trt)), pch = 21) +
        geom_errorbar(aes(ymin=mean_rich-se_rich, ymax=mean_rich+se_rich, color=as.factor(trt)), width=.2)+
        geom_vline(xintercept=2008.5, color = "grey66", lty=2)+geom_vline(xintercept=2004.5, color="grey66", lty=2) +
        labs(x = NULL, y = "Mean Species Richness", fill = "Treatment")+
        scale_color_manual(values= c("grey0", "grey36", "grey65"), guide = FALSE) +
        scale_fill_manual(values = c("grey0", "grey85", "grey100")) +
        annotate("text", x= 2004.5, y = 1, label = "fire", size = 3) +
        annotate("text", x= 2008.5, y = 1, label = "grazing", size = 3) +
        ggtitle("")

f2 <- ggplot(rich1%>%filter(func == "grass non-native", year%in%c(2005:2012)), aes(year, mean_rich)) +
        geom_line(aes(color=as.factor(trt))) +
        geom_point(aes(fill=as.factor(trt)), pch = 21) +
        geom_errorbar(aes(ymin=mean_rich-se_rich, ymax=mean_rich+se_rich, color=as.factor(trt)), width=.2)+
        geom_vline(xintercept=2008.5, color = "grey66", lty=2)+geom_vline(xintercept=2004.5, color="grey66", lty=2) +
        labs(x = NULL, y = "Mean Species Richness", fill = "Treatment")+
        scale_color_manual(values= c("grey0", "grey36", "grey65"), guide = FALSE) +
        scale_fill_manual(values = c("grey0", "grey85", "grey100")) +
        annotate("text", x= 2004.5, y = 1.6, label = "fire", size = 3) +
        annotate("text", x= 2008.5, y = 1.6, label = "grazing", size = 3) +
        ggtitle("")

f3 <- ggplot(cov1%>%filter(func == "forb native", year%in%c(2005:2012)), aes((year), meanrelcov))+
        geom_line(aes(color=trt))+
        geom_point(aes(fill=trt), pch = 21)+
        geom_errorbar(aes(ymin=meanrelcov-se_relcov, ymax=meanrelcov+se_relcov, color=trt), width=.2)+
        geom_vline(xintercept=2008.5, color = "grey66", lty=2)+geom_vline(xintercept=2004.5, color = "grey66", lty=2) +
        labs(x = NULL, y = "Mean Relative Cover (%)", fill = "Treatment") +
        scale_color_manual(values= c("grey0", "grey36", "grey65"), guide = FALSE) +
        scale_fill_manual(values = c("grey0", "grey85", "grey100")) +
        annotate("text", x= 2004.5, y = 0.05, label = "fire", size = 3) +
        annotate("text", x= 2008.5, y = 0.05, label = "grazing", size = 3) +
        ggtitle("")

f4 <- ggplot(cov1%>%filter(func == "grass non-native", year%in%c(2005:2012)), aes((year), meanrelcov))+
        geom_line(aes(color=trt))+
        geom_point(aes(fill=trt), pch = 21)+
        geom_errorbar(aes(ymin=meanrelcov-se_relcov, ymax=meanrelcov+se_relcov, color=trt), width=.2)+
        geom_vline(xintercept=2008.5, color = "grey66", lty=2)+geom_vline(xintercept=2004.5, color = "grey66", lty=2) +
        labs(x = NULL, y = "Mean Relative Cover (%)", fill = "Treatment") +
        scale_color_manual(values= c("grey0", "grey36", "grey65"), guide = FALSE) +
        scale_fill_manual(values = c("grey0", "grey85", "grey100")) +
        annotate("text", x= 2004.5, y = 0.2, label = "fire", size = 3) +
        annotate("text", x= 2008.5, y = 0.2, label = "grazing", size = 3) +
        ggtitle("")

f5 <- ggplot(litter%>%filter(year%in%c(2005:2012)), aes(year, mean_litter)) +
        geom_line(aes(color=as.factor(trt))) +
        geom_point(aes(fill=as.factor(trt)), pch = 21) +
        geom_errorbar(aes(ymin=mean_litter-se_litter, ymax=mean_litter+se_litter, color=as.factor(trt)), width=.2) +
        geom_vline(xintercept=2008.5, color = "grey66", lty =2)+geom_vline(xintercept=2004.5, color="grey66", lty = 2) +
        labs(x = NULL, y = "Mean Litter Cover (%)", fill = "Treatment") +
        scale_color_manual(values= c("grey0", "grey36", "grey65"), guide = FALSE) +
        scale_fill_manual(values = c("grey0", "grey85", "grey100")) +
        annotate("text", x= 2004.5, y = 0, label = "fire", size = 3) +
        annotate("text", x= 2008.5, y = 0, label = "grazing", size = 3) +
        ggtitle("")

#load "prism grow" from BUGR_timeseries.R
f6 <- ggplot(prism_grow%>%filter(year%in%c(2005:2012)), aes(year, prcp)) +
        geom_bar(stat = "identity", fill = "lightgrey") +
        labs(x = NULL, y = "Mean Annual Precip (mm)") +
        ggtitle("")

ggarrange(f1, f2, f3, f4, f5, f6,  ncol = 2, nrow = 3, 
          labels = c("a) Native forb", "b) Non-native grass",
                     "c) Native forb", "d) Non-native grass", "e) Litter", "f) Precipitation"),
          common.legend = TRUE, legend = "bottom", 
          font.label = list(size = 10),
          hjust = c(-0.5, -0.35, -0.5, -0.35, -0.9, -0.5))

#####
# Publication v2: with significance:
f1b <- ggplot(subset(rich1, func == "forb native"&year%in%c(2005:2008)), aes(year, mean_rich)) +
  geom_line(aes(color=as.factor(trt)), linetype="dashed") +
  geom_errorbar(aes(ymin=mean_rich-se_rich, ymax=mean_rich+se_rich, color=as.factor(trt)), width=.2)+
  geom_point(aes(fill=as.factor(trt)), pch = 21) +
  
  geom_line(data=subset(rich1, func=="forb native"&year%in%c(2008:2012)), aes(color=as.factor(trt))) +
  geom_errorbar(data=subset(rich1, func=="forb native"&year%in%c(2008:2012)), aes(ymin=mean_rich-se_rich, ymax=mean_rich+se_rich, color=as.factor(trt)), width=.2)+
  geom_point(data=subset(rich1, func=="forb native"&year%in%c(2008:2012)), aes(fill=as.factor(trt)), pch = 21) +
  
  geom_line(data=subset(rich1, func=="forb native"&year%in%c(2005:2008)&trt=="grazed burned"), aes(color=as.factor(trt))) +
  
  labs(x = NULL, y = "Mean Species Richness", fill = "Treatment")+
  scale_color_manual(values= c("grey0", "grey36", "grey65"), guide = FALSE) +
  scale_fill_manual(values = c("grey0", "grey85", "grey100")) +
  
  annotate("text", x= 2004.5, y = 1.5, label = "fire", size = 3) +
  annotate(geom = 'text', x= 2008.5, y = 2, 
           label = "atop(cattle, reintroduced)", 
           parse = TRUE, size=3)+
  geom_segment(aes(x=2004.5, y=1, xend=2004.5, yend=.5), arrow=arrow(length = unit(0.03, "npc"))) + 
  geom_segment(aes(x=2008.5, y=1, xend=2008.5, yend=.5), arrow=arrow(length = unit(0.03, "npc"))) + 
  
  ggtitle("")+
  annotate("text", x= 2005, y = 10.5, label = "*", size = 4) +
  annotate("text", x= 2006, y = 10.5, label = "**", size = 4) +
  annotate("text", x= 2007, y = 10.5, label = "**", size = 4) +
  annotate("text", x= 2008, y = 10.5, label = "**", size = 4) +
  annotate("text", x= 2009, y = 10.5, label = "**", size = 4) +
  annotate("text", x= 2010, y = 10.5, label = "*", size = 4) +
  annotate("text", x= 2011, y = 10.5, label = "*", size = 4) +
  annotate("text", x= 2012, y = 10.5, label = "", size = 4)

f2b <- ggplot(subset(rich1, func == "grass non-native"&year%in%c(2005:2008)), aes(year, mean_rich)) +
  geom_line(aes(color=as.factor(trt)), linetype="dashed") +
  geom_errorbar(aes(ymin=mean_rich-se_rich, ymax=mean_rich+se_rich, color=as.factor(trt)), width=.2)+
  geom_point(aes(fill=as.factor(trt)), pch = 21) +
  
  geom_line(data=subset(rich1, func=="grass non-native"&year%in%c(2008:2012)), aes(color=as.factor(trt))) +
  geom_errorbar(data=subset(rich1, func=="grass non-native"&year%in%c(2008:2012)), aes(ymin=mean_rich-se_rich, ymax=mean_rich+se_rich, color=as.factor(trt)), width=.2)+
  geom_point(data=subset(rich1, func=="grass non-native"&year%in%c(2008:2012)), aes(fill=as.factor(trt)), pch = 21) +
  
  geom_line(data=subset(rich1, func=="grass non-native"&year%in%c(2005:2008)&trt=="grazed burned"), aes(color=as.factor(trt))) +
  
  geom_errorbar(aes(ymin=mean_rich-se_rich, ymax=mean_rich+se_rich, color=as.factor(trt)), width=.2)+
  geom_point(aes(fill=as.factor(trt)), pch = 21) +
  labs(x = NULL, y = "Mean Species Richness", fill = "Treatment")+
  scale_color_manual(values= c("grey0", "grey36", "grey65"), guide = FALSE) +
  scale_fill_manual(values = c("grey0", "grey85", "grey100")) +
  ggtitle("")+
  annotate("text", x= 2004.5, y = 1.3, label = "fire", size = 3) +
  annotate(geom = 'text', x= 2008.5, y = 1.4, 
           label = "atop(cattle, reintroduced)", 
           parse = TRUE, size=3)+
  geom_segment(aes(x=2004.5, y=1.2, xend=2004.5, yend=1.1), arrow=arrow(length = unit(0.03, "npc"))) + 
  geom_segment(aes(x=2008.5, y=1.2, xend=2008.5, yend=1.1), arrow=arrow(length = unit(0.03, "npc")))

f3b <- ggplot(subset(cov1,func == "forb native"&year%in%c(2005:2008)), aes((year), meanrelcov*100))+

  geom_line(aes(color=as.factor(trt)), linetype="dashed") +
  geom_errorbar(aes(ymin=meanrelcov*100-se_relcov*100, ymax=meanrelcov*100+se_relcov*100, color=trt), width=.2)+
  geom_point(aes(fill=as.factor(trt)), pch = 21) +
  
  geom_line(data=subset(cov1, func=="forb native"&year%in%c(2008:2012)), aes(color=as.factor(trt))) +
  geom_errorbar(data=subset(cov1, func=="forb native"&year%in%c(2008:2012)), aes(ymin=meanrelcov*100-se_relcov*100, ymax=meanrelcov*100+se_relcov*100, color=trt), width=.2)+
  geom_point(data=subset(cov1, func=="forb native"&year%in%c(2008:2012)), aes(fill=as.factor(trt)), pch = 21) +
  
  geom_line(data=subset(cov1, func=="forb native"&year%in%c(2005:2008)&trt=="grazed burned"), aes(color=as.factor(trt))) +
  
  geom_errorbar(aes(ymin=meanrelcov*100-se_relcov*100, ymax=meanrelcov*100+se_relcov*100, color=trt), width=.2)+
  geom_point(aes(fill=trt), pch = 21)+
  labs(x = NULL, y = "Mean Relative Cover (%)", fill = "Treatment") +
  scale_color_manual(values= c("grey0", "grey36", "grey65"), guide = FALSE) +
  scale_fill_manual(values = c("grey0", "grey85", "grey100")) +
  ggtitle("")+
  annotate("text", x= 2005, y = 75, label = "**", size = 4) +
  annotate("text", x= 2006, y = 75, label = "***", size = 4) +
  annotate("text", x= 2007, y = 75, label = "**", size = 4) +
  annotate("text", x= 2008, y = 75, label = "**", size = 4) +
  annotate("text", x= 2009, y = 75, label = "**", size = 4) +
  annotate("text", x= 2010, y = 75, label = "**", size = 4) +
  annotate("text", x= 2011, y = 75, label = "", size = 4) +
  annotate("text", x= 2012, y = 75, label = "*", size = 4) +
  annotate("text", x= 2004.5, y = 8, label = "fire", size = 3) +
  annotate(geom = 'text', x= 2008.5, y = 12, 
           label = "atop(cattle, reintroduced)", 
           parse = TRUE, size=3)+
  geom_segment(aes(x=2004.5, y=5, xend=2004.5, yend=.5), arrow=arrow(length = unit(0.03, "npc"))) + 
  geom_segment(aes(x=2008.5, y=5, xend=2008.5, yend=.5), arrow=arrow(length = unit(0.03, "npc"))) +   ggtitle("")

f4b <- ggplot(subset(cov1, func == "grass non-native"&year%in%c(2005:2008)), aes((year), meanrelcov*100))+
  
  geom_line(aes(color=as.factor(trt)), linetype="dashed") +
  geom_errorbar(aes(ymin=meanrelcov*100-se_relcov*100, ymax=meanrelcov*100+se_relcov*100, color=trt), width=.2)+
  geom_point(aes(fill=as.factor(trt)), pch = 21) +
  
  geom_line(data=subset(cov1, func=="grass non-native"&year%in%c(2008:2012)), aes(color=as.factor(trt))) +
  geom_errorbar(data=subset(cov1, func=="grass non-native"&year%in%c(2008:2012)), aes(ymin=meanrelcov*100-se_relcov*100, ymax=meanrelcov*100+se_relcov*100, color=trt), width=.2)+
  geom_point(data=subset(cov1, func=="grass non-native"&year%in%c(2008:2012)), aes(fill=as.factor(trt)), pch = 21) +
  
  geom_line(data=subset(cov1, func=="grass non-native"&year%in%c(2005:2008)&trt=="grazed burned"), aes(color=as.factor(trt))) +
  
  labs(x = NULL, y = "Mean Relative Cover (%)", fill = "Treatment") +
  scale_color_manual(values= c("grey0", "grey36", "grey65"), guide = FALSE) +
  scale_fill_manual(values = c("grey0", "grey85", "grey100")) +
  annotate("text", x= 2004.5, y = 9, label = "fire", size = 3) +
  annotate(geom = 'text', x= 2008.5, y = 15, 
           label = "atop(cattle, reintroduced)", 
           parse = TRUE, size=3)+
  geom_segment(aes(x=2004.5, y=6, xend=2004.5, yend=.5), arrow=arrow(length = unit(0.03, "npc"))) + 
  geom_segment(aes(x=2008.5, y=6, xend=2008.5, yend=.5), arrow=arrow(length = unit(0.03, "npc"))) + 
  ggtitle("")

f5b <- ggplot(subset(litter, year%in%c(2005:2008)), aes(year, mean_litter)) +
  geom_line(aes(color=as.factor(trt)), linetype="dashed") +
  geom_errorbar(aes(ymin=mean_litter-se_litter, ymax=mean_litter+se_litter, color=as.factor(trt)), width=.2) +
  geom_point(aes(fill=as.factor(trt)), pch = 21) +
  
  geom_line(data=subset(litter, year%in%c(2008:2012)), aes(color=as.factor(trt))) +
  geom_errorbar(data=subset(litter, year%in%c(2008:2012)), aes(ymin=mean_litter-se_litter, ymax=mean_litter+se_litter, color=as.factor(trt)), width=.2)+
  geom_point(data=subset(litter, year%in%c(2008:2012)), aes(fill=as.factor(trt)), pch = 21) +
  
  geom_line(data=subset(litter, year%in%c(2005:2008)&trt=="grazed burned"), aes(color=as.factor(trt))) +
  

  labs(x = NULL, y = "Mean Litter Cover (%)", fill = "Treatment") +
  scale_color_manual(values= c("grey0", "grey36", "grey65"), guide = FALSE) +
  scale_fill_manual(values = c("grey0", "grey85", "grey100")) +
  ggtitle("")+
  annotate("text", x= 2006, y = 40, label = "**", size = 4) +
  annotate("text", x= 2007, y = 40, label = "***", size = 4) +
  annotate("text", x= 2008, y = 40, label = "**", size = 4) +
  annotate("text", x= 2009, y = 40, label = "***", size = 4) +
  annotate("text", x= 2010, y = 40, label = "**", size = 4) +
  annotate("text", x= 2011, y = 40, label = "***", size = 4) +
  annotate("text", x= 2012, y = 40, label = "***", size = 4) +
  annotate("text", x= 2004.5, y = 5, label = "fire", size = 3) +
  geom_segment(aes(x=2004.5, y=2.5, xend=2004.5, yend=.5), arrow=arrow(length = unit(0.03, "npc"))) + 
  ggtitle("")

#load "prism grow" from BUGR_timeseries.R
f6b <- ggplot(subset(prism_grow, year%in%c(2005:2012)), aes(year, prcp)) +
  geom_bar(stat = "identity", fill = "lightgrey") +
  labs(x = NULL, y = "Mean Annual Precip (mm)") +
  ggtitle("")

ggarrange(f1b, f2b, f3b, f4b, f5b, f6b,  ncol = 2, nrow = 3, 
          labels = c("a) Native forb", "b) Non-native grass",
                     "c) Native forb", "d) Non-native grass", "e) Litter", "f) Precipitation"),
          common.legend = TRUE, legend = "bottom", 
          font.label = list(size = 10),
          hjust = c(-0.5, -0.35, -0.5, -0.35, -0.9, -0.5))
###
#Litter and native forb regression
###
forb_rich <- rich %>%
  filter(func == "forb native")
rich_litter <- left_join(forb_rich, envdat) %>%
  filter(year%in%c(2006:2012)) %>%
  filter(litter != "NA")
ggplot(rich_litter, aes(litter, richness)) +
  geom_point(aes(color=as.factor(trt))) +
  labs(x = "Litter Cover (%)", y = "Native Forb Richness", color = "Treatment") +
  geom_smooth(method = lm, se = FALSE, color = "black") + 
  theme_classic() +
  scale_color_manual(values= c("grey0", "grey36", "grey65"))
##########
# Indicator Species
###########
#library(indicspecies)
#grztog3<-grztog2%>%
#  mutate(prepost=ifelse(year<2009, "pre", "post"))%>%
#  mutate(trtgroup=paste(graze, prepost, sep="_"))%>%
#  mutate(rep=paste(transect.quad, year))%>%
#  select(-transect, -transect.quad,-year, -quadrat, -graze, -prepost, -thermal, -spname, -status, -func)

#indic_treatments<-select(grztog3, 3, 4)%>%
#  unique()
#indic_species<-select(grztog3, 1, 2, 4)%>%
#  spread(spcode, cover, fill=0)%>%
##  remove_rownames()%>%
#  column_to_rownames("rep")
#indic_species<-decostand(indic_species, "total")

#indic_treatments2<-indic_species%>%
#  rownames_to_column("rep")%>%
#  select(1)
#indic_treatments3<-left_join(indic_treatments2, indic_treatments)

#indicators<-multipatt(indic_species, indic_treatments$trtgroup, func="IndVal.g", control=how(nperm=999))

#indsum<-indicators$sign%>%
#  rownames_to_column("species")%>%
#  filter(p.value<.05)
#indsum1<-left_join(indsum, SC, by=c("species"="spcode"))%>%
#  select(-index, -stat, -p.value)

#write.csv(indsum1, file = "Grazing_indicator_sp.csv")