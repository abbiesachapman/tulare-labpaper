## Set your datpath!! (file in Data_cleaning)

library(tidyverse)
library(readr)
library(vegan)
library(MASS)
library(dplyr)
library(RVAideMemoire) #for posthoc tests on permanova

# Import csv file, transform to wide data, call it data
all.dat <- read_csv(paste(datpath_clean, "/alldatsptrt.csv", sep=""))
all.dat<-as.data.frame(all.dat)
#dat2<-dat %>% mutate(cover=cover+0.00001) #a work around for data with lots of zeros
#create wide data, first filter so year is 2005 to 2012
all.data<- all.dat %>% dplyr::select(-X1, -spcode) %>% filter(year>2004 & year < 2013) %>% spread(spname, cover)
all.data[is.na(all.data)] <- 0 #replace NAs with 0 (species not counted in plots have NAs when wide dataset created)
levels(as.factor(all.dat$quadratNew))
levels(as.factor(all.dat$thermal))
levels(as.factor(all.dat$spname)) #any to remove?
levels(as.factor(all.dat$burn))
levels(as.factor(all.dat$graze))

#Import environmental data - NADP pinnacles deposition data and NOAA San Jose temperature data
env <- read_csv(paste(datpath_clean, "/NTN-CA66-deposition.csv", sep=""))
env<-env %>% rename(year=yr) %>% filter(year>2004 & year < 2013)
temp <- read_csv(paste(datpath_clean, "/san_jose_clim.csv", sep=""),skip=9)
temp <- temp %>% rename(year=DATE, temp=TAVG) %>% filter(year>2004 & year <2013) %>% dplyr::select(year, temp)
env<-merge(env,temp)

#merge env to match full data
all.env<-merge(env,all.data) %>% dplyr::select(NH4, NO3, totalN, ppt, temp)

#check count of thermal factor
all.data %>% 
  group_by(thermal, burn, graze) %>%
  summarise(no_rows = length(thermal)) 

#check count of graze and burn factor
all.data %>% 
  group_by(graze, burn) %>%
  summarise(no_rows = length(graze))


str(all.data)

all.data$ID <- seq.int(nrow(all.data))
plotnames<-all.data[,1]
cover.gb<- all.data %>% dplyr::select(-c(1:5)) #wide data with ID columns removed, only species/cover for NMDS
rownames(cover.gb)<-plotnames

#check for empty rows
cover.Biodrop.gb<-cover.gb[rowSums(cover.gb[, (1:157)]) ==0, ] #if no empty rows, next step not needed
#cover.Biodrop.gb<-cover.gb[rowSums(cover.gb[, (1:157)])  >0 ]#remove empty rows

#if needed, relativize by row or column or calculate presence/absence
#cover.rowsums <- rowSums(cover.Bio [1:157])
#cover.relrow <- data.frame(cover.Bio /cover.rowsums)
#cover.colmax<-sapply(cover.Bio ,max)
#cover.relcolmax <- data.frame(sweep(cover.Bio ,2,cover.colmax,'/'))
#cover.pa <- cover.Bio %>% mutate_each(funs(ifelse(.>0,1,0)), 1:57)


######################
#NMDS
######################
#make bray-curtis dissimilarity matrix
gb.bcd <- vegdist(cover.gb)

#quick run to check out NMS ordination
gb.mds0 <-isoMDS(gb.bcd) #runs nms only once
gb.mds0  #by default 2 dimensions returned, stress is 6.4, converged
ordiplot(gb.mds0) #ugly

#prefer to run multiple NMS ordinations
#with different starting configurations, and choose the best
#this function does that, and in addition does several other steps too including: 
#standardizing the data (though fuction call below turns this off with autotransform=F)
#calculating distance matrix (default bray-curtis)
#running NMDS with random starts
#rotation of axes to maximize variance of site scores on axis 1
#calculate species scores based on weighted averaging
#help(metaMDS)
gb.mds<-metaMDS(cover.gb, trace = TRUE, autotransform=T, trymax=100, k=3) #runs several with different starting configurations
#trace= TRUE will give output for step by step what its doing
#default is 2 dimensions, can put k=4 for 4 dimensions
gb.mds #solution did not converge after 100 tries
summary(gb.mds)

#quick plot of results
stressplot(gb.mds, gb.bcd) #stressplot to show fit
ordiplot(gb.mds)

#overlay environmental variables on full data
envvec.nms.gb<-envfit(gb.mds,all.env, na.rm=TRUE)
envvec.nms.gb
plot(gb.mds)
plot(envvec.nms.gb) #add vectors to previous ordination

#store scores in new dataframe
spscores1.gb<-scores(gb.mds,display="sites",choices=1)
spscores2.gb<-scores(gb.mds,display="sites",choices=2)
tplots.gb<-all.data[,2]
tplot_levels_gb<-levels(tplots.gb)
spscoresall.gb<-data.frame(tplots.gb,spscores1.gb,spscores2.gb)

#make nicer plot colored based on thermal, shapes on pre/post fire
#help(ordiplot)
#first, set colors and shapes
cols1<-all.data %>% dplyr::select(burn, graze) %>% mutate(color = "forestgreen", 
                                                   color = ifelse(burn == "burned" & graze=="grazed", "red",
                                                                  ifelse(burn=="burned" & graze =="ungrazed", "orange",
                                                                         ifelse(burn=="unburned" & graze=="graze", "purple", color)))) #colors based on burn trt
Lcols <- rep(c("orange", "Red", "forestgreen","yellow4")) #colors for the legend
shapes <- all.data %>% dplyr::select(year) %>%
  mutate(shape = 0, shape = ifelse(year == "2005", 16, ifelse(year=="2006", 17, ifelse(year=="2007", 18,
                                    ifelse(year=="2008", 15, ifelse(year=="2009", 1, ifelse(year=="2010",2, 
                                          ifelse(year=="2011", 5, shape)))))))) #shapes based on grazing 
Lshapes <-rep(c(16,17,18,15,1,2,5,0))#shapes for legend
#make the plot
gb.plot <- ordiplot(gb.mds, choices=c(1,2), type = "none")   #Set up the plot
points(spscoresall.gb$NMDS1,spscoresall.gb$NMDS2,col=cols1$color,pch=shapes$shape) 
#plot(envvec.nms.gb, col="green")
#text(gb.mds, display = "species", cex=0.5, col="grey30") #label species
#legend("bottomright",legend=levels(as.factor(cols1$burn)), col=Lcols, pch=15, cex=0.9,inset=0.1,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)
legend("topright",legend=levels(as.factor(shapes$year)), col="black", pch=Lshapes, cex=0.9,inset=0.1,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)

#make nicer plot colored based on thermal, shapes on pre/post fire
#help(ordiplot)
#first, set colors and shapes
cols1<- all.data %>% dplyr::select(thermal) %>% mutate(color = "black", 
                                                   color = ifelse(thermal == "cool", "purple", 
                                                                  ifelse(thermal=="moderate", "orange", 
                                                                         ifelse(thermal=="very cool", "blue",
                                                                                ifelse(thermal=="very warm", "red", color))))) #colors based on thermal group
Lcols <- rep(c("purple", "orange", "blue","red","black")) #colors for the legend
shapes <- all.data %>% dplyr::select(year) %>%
  mutate(shape = 0, shape = ifelse(year == "2005", 16, ifelse(year=="2006", 17, ifelse(year=="2007", 18,
                                                                                       ifelse(year=="2008", 15, ifelse(year=="2009", 1, ifelse(year=="2010",2, 
                                                                                                                                               ifelse(year=="2011", 5, shape)))))))) #shapes based on grazing 
Lshapes <-rep(c(16,17,18,15,1,2,5,0))#shapes for legend
#make the plot
gb.plot <- ordiplot(gb.mds, choices=c(1,2), type = "none")   #Set up the plot
points(spscoresall.gb$NMDS1,spscoresall.gb$NMDS2,col=cols1$color,pch=shapes$shape) 
#plot(envvec.nms.gb, col="green")
#text(gb.mds, display = "species", cex=0.5, col="grey30") #label species
legend("bottomright",legend=levels(as.factor(cols1$thermal)), col=Lcols, pch=15, cex=0.9,inset=0.1,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)
legend("topright",legend=levels(as.factor(shapes$year)), col="black", pch=Lshapes, cex=0.9,inset=0.1,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)


###early years####
#create wide data, first filter so year is 2005 to 2008
all.data.early<- all.dat %>% dplyr::select(-X1, -spcode) %>% filter(year>2004 & year < 2009) %>% spread(spname, cover)
all.data.early[is.na(all.data.early)] <- 0 #replace NAs with 0 (species not counted in plots have NAs when wide dataset created)
levels(as.factor(all.dat$quadratNew))
levels(as.factor(all.dat$thermal))
levels(as.factor(all.dat$spname)) #any to remove?
levels(as.factor(all.dat$burn))
levels(as.factor(all.dat$graze))

#Import environmental data - NADP pinnacles deposition data and NOAA San Jose temperature data
env <- read_csv(paste(datpath_clean, "/NTN-CA66-deposition.csv", sep=""))
env<-env %>% rename(year=yr) %>% filter(year>2004 & year < 2009)
temp <- read_csv(paste(datpath_clean, "/san_jose_clim.csv", sep=""),skip=9)
temp <- temp %>% rename(year=DATE, temp=TAVG) %>% filter(year>2004 & year <2009) %>% dplyr::select(year, temp)
env<-merge(env,temp)

#merge env to match full data
all.env.early<-merge(env,all.data.early) %>% dplyr::select(NH4, NO3, totalN, ppt, temp)

#check count of thermal factor
all.data.early %>% 
  group_by(thermal) %>%
  summarise(no_rows = length(thermal)) 

#check count of graze and burn factor
all.data.early %>% 
  group_by(graze, burn) %>%
  summarise(no_rows = length(graze))


str(all.data)

all.data.early$ID <- seq.int(nrow(all.data))
plotnames<-all.data.early[,1]
cover.gb.early<- all.data.early %>% dplyr::select(-c(1:5)) #wide data with ID columns removed, only species/cover for NMDS
rownames(cover.gb.early)<-plotnames

#check for empty rows
cover.Biodrop.gb.early<-cover.gb.early[rowSums(cover.gb.early[, (1:156)]) ==0, ] #if no empty rows, next step not needed
#cover.Biodrop.gb<-cover.gb[rowSums(cover.gb[, (1:157)])  >0 ]#remove empty rows

#if needed, relativize by row or column or calculate presence/absence
#cover.rowsums <- rowSums(cover.Bio [1:157])
#cover.relrow <- data.frame(cover.Bio /cover.rowsums)
#cover.colmax<-sapply(cover.Bio ,max)
#cover.relcolmax <- data.frame(sweep(cover.Bio ,2,cover.colmax,'/'))
#cover.pa <- cover.Bio %>% mutate_each(funs(ifelse(.>0,1,0)), 1:57)


######################
#NMDS
######################
#make bray-curtis dissimilarity matrix
gb.bcd.early <- vegdist(cover.gb.early)

gb.mds.early<-metaMDS(gb.bcd.early , trace = TRUE, autotransform=T, trymax=100, k=3) #runs several with different starting configurations
#trace= TRUE will give output for step by step what its doing
#default is 2 dimensions, can put k=4 for 4 dimensions
gb.mds.early #solution did not converge after 100 tries
summary(gb.mds.early)

#quick plot of results
stressplot(gb.mds.early, gb.bcd.early) #stressplot to show fit
ordiplot(gb.mds.early)

#overlay environmental variables on full data
envvec.nms.gb.early<-envfit(gb.mds.early,all.env.late, na.rm=TRUE)
envvec.nms.gb.early
plot(gb.mds.early)
plot(envvec.nms.gb.early) #add vectors to previous ordination

#store scores in new dataframe
spscores1.gb.e<-scores(gb.mds.early,display="sites",choices=1)
spscores2.gb.e<-scores(gb.mds.early,display="sites",choices=2)
tplots.gb.e<-all.data.early[,2]
tplot_levels_gb_e<-levels(tplots.gb.e)
spscoresall.gb.e<-data.frame(tplots.gb.e,spscores1.gb.e,spscores2.gb.e)

#make nicer plot colored based on thermal, shapes on pre/post fire
#help(ordiplot)
#first, set colors and shapes
cols1.e<-all.data.early %>% dplyr::select(burn, graze) %>% mutate(color = "forestgreen", 
                                                          color = ifelse(burn == "burned" & graze=="grazed", "red",
                                                                         ifelse(burn=="burned" & graze =="ungrazed", "orange",
                                                                                ifelse(burn=="unburned" & graze=="graze", "purple", color)))) #colors based on burn trt
Lcols.e <- rep(c("orange", "Red", "forestgreen","yellow4")) #colors for the legend
shapes.e <- all.data.early %>% dplyr::select(year) %>%
  mutate(shape = 0, shape = ifelse(year == "2005", 16, ifelse(year=="2006", 17, ifelse(year=="2007", 18,
                                                                                       ifelse(year=="2008", 15, shape))))) #shapes based on grazing 
Lshapes.e <-rep(c(16,17,18,15))#shapes for legend
#make the plot
gb.plot.e <- ordiplot(gb.mds.early, choices=c(1,2), type = "none")   #Set up the plot
points(spscoresall.gb.e$NMDS1,spscoresall.gb.e$NMDS2,col=cols1.e$color,pch=shapes.e$shape) 
#plot(envvec.nms.gb, col="green")
#text(gb.mds, display = "species", cex=0.5, col="grey30") #label species
legend("bottomright",legend=levels(as.factor(cols1$thermal)), col=Lcols, pch=15, cex=0.9,inset=0.1,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)
legend("topright",legend=levels(as.factor(shapes.e$year)), col="black", pch=Lshapes.e, cex=0.9,inset=0.1,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)

###late years####
#create wide data, first filter so year is 2005 to 2008
all.data.late<- all.dat %>% dplyr::select(-X1, -spcode) %>% filter(year>2008 & year < 2012) %>% spread(spname, cover)
all.data.late[is.na(all.data.late)] <- 0 #replace NAs with 0 (species not counted in plots have NAs when wide dataset created)

#Import environmental data - NADP pinnacles deposition data and NOAA San Jose temperature data
env <- read_csv(paste(datpath_clean, "/NTN-CA66-deposition.csv", sep=""))
env<-env %>% rename(year=yr) %>% filter(year>2008 & year < 2012)
temp <- read_csv(paste(datpath_clean, "/san_jose_clim.csv", sep=""),skip=9)
temp <- temp %>% rename(year=DATE, temp=TAVG) %>% filter(year>2008 & year <2012) %>% dplyr::select(year, temp)
env<-merge(env,temp)

#merge env to match full data
all.env.late<-merge(env,all.data.late) %>% dplyr::select(NH4, NO3, totalN, ppt, temp)

#check count of thermal factor
all.data.late %>% 
  group_by(thermal) %>%
  summarise(no_rows = length(thermal)) 

#check count of graze and burn factor
all.data.late %>% 
  group_by(graze, burn) %>%
  summarise(no_rows = length(graze))


str(all.data.late)

all.data.late$ID <- seq.int(nrow(all.data.late))
plotnames<-all.data.late[,1]
cover.gb.late<- all.data.late %>% dplyr::select(-c(1:5)) #wide data with ID columns removed, only species/cover for NMDS
rownames(cover.gb.late)<-plotnames

#check for empty rows
cover.Biodrop.gb.late<-cover.gb.late[rowSums(cover.gb.late[, (1:157)]) ==0, ] #if no empty rows, next step not needed
#cover.Biodrop.gb<-cover.gb[rowSums(cover.gb[, (1:157)])  >0 ]#remove empty rows

#if needed, relativize by row or column or calculate presence/absence
#cover.rowsums <- rowSums(cover.Bio [1:157])
#cover.relrow <- data.frame(cover.Bio /cover.rowsums)
#cover.colmax<-sapply(cover.Bio ,max)
#cover.relcolmax <- data.frame(sweep(cover.Bio ,2,cover.colmax,'/'))
#cover.pa <- cover.Bio %>% mutate_each(funs(ifelse(.>0,1,0)), 1:57)


######################
#NMDS
######################
#make bray-curtis dissimilarity matrix
gb.bcd.late <- vegdist(cover.gb.late)

gb.mds.late<-metaMDS(gb.bcd.late , trace = TRUE, autotransform=T, trymax=100, k=3) #runs several with different starting configurations
#trace= TRUE will give output for step by step what its doing
#default is 2 dimensions, can put k=4 for 4 dimensions
gb.mds.late #solution did not converge after 100 tries
summary(gb.mds.late)

#quick plot of results
stressplot(gb.mds.late, gb.bcd.late) #stressplot to show fit
ordiplot(gb.mds.late)

#overlay environmental variables on full data
envvec.nms.gb.late<-envfit(gb.mds.late,all.env.late, na.rm=TRUE)
envvec.nms.gb.late
plot(gb.mds.late)
plot(envvec.nms.gb.late) #add vectors to previous ordination

#store scores in new dataframe
spscores1.gb.l<-scores(gb.mds.late,display="sites",choices=1)
spscores2.gb.l<-scores(gb.mds.late,display="sites",choices=2)
tplots.gb.l<-all.data.late[,2]
tplot_levels_gb_l<-levels(tplots.gb.l)
spscoresall.gb.l<-data.frame(tplots.gb.l,spscores1.gb.l,spscores2.gb.l)

#make nicer plot colored based on thermal, shapes on pre/post fire
#help(ordiplot)
#first, set colors and shapes
cols1.l<-all.data.late %>% dplyr::select(burn, graze) %>% mutate(color = "forestgreen", 
                                                                  color = ifelse(burn == "burned" & graze=="grazed", "red",
                                                                                 ifelse(burn=="burned" & graze =="ungrazed", "orange",
                                                                                        ifelse(burn=="unburned" & graze=="graze", "purple", color)))) #colors based on burn trt
Lcols.l <- rep(c("orange", "Red", "forestgreen","yellow4")) #colors for the legend
shapes.l <- all.data.late %>% dplyr::select(year) %>%
  mutate(shape = 0, shape = ifelse(year=="2009", 1, ifelse(year=="2010",2, 
                                                           ifelse(year=="2011", 5, shape)))) #shapes based on grazing 
Lshapes.l <-rep(c(1,2,5))#shapes for legend
#make the plot
gb.plot.l <- ordiplot(gb.mds.late, choices=c(1,2), type = "none")   #Set up the plot
points(spscoresall.gb.l$NMDS1,spscoresall.gb.l$NMDS2,col=cols1.l$color,pch=shapes.l$shape) 
#plot(envvec.nms.gb, col="green")
#text(gb.mds, display = "species", cex=0.5, col="grey30") #label species
legend("bottomright",legend=levels(as.factor(cols1$thermal)), col=Lcols, pch=15, cex=0.9,inset=0.1,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)
legend("topright",legend=levels(as.factor(shapes.l$year)), col="black", pch=Lshapes.l, cex=0.9,inset=0.1,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)

#make nicer plot colored based on thermal, shapes on pre/post fire
#help(ordiplot)
#first, set colors and shapes
cols2.l<- all.data.late %>% dplyr::select(thermal) %>% mutate(color = "black", 
                                                   color = ifelse(thermal == "cool", "purple", 
                                                                  ifelse(thermal=="moderate", "orange", 
                                                                         ifelse(thermal=="very cool", "blue",
                                                                                ifelse(thermal=="very warm", "red", color))))) #colors based on thermal group
Lcols.l2 <- rep(c("purple", "orange", "blue","red","black")) #colors for the legend
shapes.l <- all.data.late %>% dplyr::select(year) %>%
  mutate(shape = 0, shape = ifelse(year=="2009", 1, ifelse(year=="2010",2, 
                                                           ifelse(year=="2011", 5, shape)))) #shapes based on grazing 
Lshapes.l <-rep(c(1,2,5))#shapes for legend
#make the plot
gb.plot.l <- ordiplot(gb.mds.late, choices=c(1,2), type = "none")   #Set up the plot
points(spscoresall.gb.l$NMDS1,spscoresall.gb.l$NMDS2,col=cols2.l$color,pch=shapes.l$shape) 
#plot(envvec.nms.gb, col="green")
#text(gb.mds, display = "species", cex=0.5, col="grey30") #label species
legend("bottomright",legend=levels(as.factor(cols2.l$thermal)), col=Lcols.l2, pch=15, cex=0.9,inset=0.1,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)
legend("topright",legend=levels(as.factor(shapes.l$year)), col="black", pch=Lshapes.l, cex=0.9,inset=0.1,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)
