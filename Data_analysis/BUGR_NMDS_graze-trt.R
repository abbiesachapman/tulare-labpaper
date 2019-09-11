## Set your datpath!! (file in Data_cleaning)
library(tidyverse)
library(readr)
library(vegan)
library(MASS)
library(dplyr)
library(RVAideMemoire) #for posthoc tests on permanova
library(indicspecies)
library(gridExtra)
library(goeveg)#for scree plot of NMDS to test number of dimensions

# Import csv file, transform to wide data, call it data
all.dat <- read_csv(paste(datpath_clean, "/alldatsptrt.csv", sep=""))
all.dat<-as.data.frame(all.dat)
all.dat<-all.dat %>% dplyr::select(-status, -type)

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

plotnames<-all.data[,1]
cover.gb<- all.data %>% dplyr::select(-c(1:6)) #wide data with ID columns removed, only species/cover for NMDS
rownames(cover.gb)<-plotnames

#check for empty rows
cover.Biodrop.gb<-cover.gb[rowSums(cover.gb[, (1:156)]) ==0, ] #if no empty rows, next step not needed
#cover.Biodrop.gb<-cover.gb[rowSums(cover.gb[, (1:156)])  >0 ]#remove empty rows

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

########MODERATE FOR ALL YEARS##############
#create wide data, first filter so year is 2005 to 2012
mod.data<- all.dat %>% dplyr::select(-X1, -spcode) %>% filter(year>2004 & year < 2013) %>% filter(thermal=="moderate")%>% spread(spname, cover)
mod.data[is.na(mod.data)] <- 0 #replace NAs with 0 (species not counted in plots have NAs when wide dataset created)
levels(as.factor(mod.data$quadratNew))
levels(as.factor(mod.data$thermal))
levels(as.factor(mod.data$burn))
levels(as.factor(mod.data$graze))

#Import environmental data - NADP pinnacles deposition data and NOAA San Jose temperature data
env <- read_csv(paste(datpath_clean, "/NTN-CA66-deposition.csv", sep=""))
env<-env %>% rename(year=yr) %>% filter(year>2004 & year < 2013)
temp <- read_csv(paste(datpath_clean, "/san_jose_clim.csv", sep=""),skip=9)
temp <- temp %>% rename(year=DATE, temp=TAVG) %>% filter(year>2004 & year <2013) %>% dplyr::select(year, temp)
env<-merge(env,temp)

#merge env to match full data
mod.env<-merge(env,mod.data) %>% dplyr::select(NH4, NO3, totalN, ppt, temp)

#check count of graze and burn factor
mod.data %>% 
  group_by(graze, burn) %>%
  summarise(no_rows = length(graze))

plotnames<-mod.data[,1]
cover.mod<- mod.data %>% dplyr::select(-c(1:5)) #wide data with ID columns removed, only species/cover for NMDS
rownames(cover.mod)<-plotnames

#check for empty rows
cover.Biodrop.mod<-cover.gb[rowSums(cover.mod[, (1:156)]) ==0, ] 

######################
#NMDS
######################
#make bray-curtis dissimilarity matrix
mod.bcd <- vegdist(cover.mod)
#run mds
mod.mds<-metaMDS(cover.mod, trace = TRUE, autotransform=T, trymax=100, k=2) #runs several with different starting configurations
#trace= TRUE will give output for step by step what its doing
#default is 2 dimensions, can put k=4 for 4 dimensions
mod.mds #solution did not converge after 100 tries
summary(mod.mds)

#quick plot of results
stressplot(mod.mds, mod.bcd) #stressplot to show fit
ordiplot(mod.mds)

#overlay environmental variables on full data
envvec.nms.mod<-envfit(mod.mds,mod.env, na.rm=TRUE)
envvec.nms.mod
plot(mod.mds)
plot(envvec.nms.mod) #add vectors to previous ordination

#store scores in new dataframe
spscores1.mod<-scores(mod.mds,display="sites",choices=1)
spscores2.mod<-scores(mod.mds,display="sites",choices=2)
tplots.mod<-mod.data[,2]
tplot_levels_mod<-levels(tplots.mod)
spscoresall.mod<-data.frame(tplots.mod,spscores1.mod,spscores2.mod)

#make nicer plot colored based on burning and grazing, shapes on year
#help(ordiplot)
#first, set colors and shapes
cols1<-mod.data %>% dplyr::select(burn, graze) %>% mutate(color = "forestgreen", 
                                                          color = ifelse(burn == "burned" & graze=="grazed", "red",
                                                                         ifelse(burn=="burned" & graze =="ungrazed", "orange",
                                                                                ifelse(burn=="unburned" & graze=="graze", "purple", color)))) #colors based on burn trt
Lcols <- rep(c("orange", "Red", "forestgreen","yellow4")) #colors for the legend
shapes <- mod.data %>% dplyr::select(year) %>%
  mutate(shape = 0, shape = ifelse(year == "2005", 16, ifelse(year=="2006", 17, ifelse(year=="2007", 18,
                                                                                       ifelse(year=="2008", 15, ifelse(year=="2009", 1, ifelse(year=="2010",2, 
                                                                                                                                               ifelse(year=="2011", 5, shape)))))))) #shapes based on grazing 
Lshapes <-rep(c(16,17,18,15,1,2,5,0))#shapes for legend
#make the plot
mod.plot <- ordiplot(mod.mds, choices=c(1,2), xlim=c(-0.05,0.1), ylim=c(-0.05,0.05),type = "none")   #Set up the plot
points(spscoresall.mod$NMDS1,spscoresall.mod$NMDS2,col=cols1$color,pch=shapes$shape) 
#plot(envvec.nms.gb, col="green")
#text(mod.mds, display = "species", cex=0.5, col="grey30") #label species
#legend("bottomright",legend=levels(as.factor(cols1$burn)), col=Lcols, pch=15, cex=0.9,inset=0.1,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)
legend("topright",legend=levels(as.factor(shapes$year)), col="black", pch=Lshapes, cex=0.9,inset=0.1,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)


###early years####
#create wide data, first filter so year is 2005 to 2008
all.data.early<- all.dat %>% dplyr::select(-X1, -spcode) %>% filter(year>2004 & year < 2009) %>% spread(spname, cover)
all.data.early[is.na(all.data.early)] <- 0 #replace NAs with 0 (species not counted in plots have NAs when wide dataset created)

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
  group_by(thermal, graze, burn) %>%
  summarise(no_rows = length(thermal)) 

#check count of graze and burn factor
all.data.early %>% 
  group_by(graze, burn) %>%
  summarise(no_rows = length(graze))

plotnames<-all.data.early[,1]
cover.gb.early<- all.data.early %>% dplyr::select(-c(1:5)) #wide data with ID columns removed, only species/cover for NMDS
rownames(cover.gb.early)<-plotnames

#check for empty rows
cover.Biodrop.gb.early<-cover.gb.early[rowSums(cover.gb.early[, (1:156)]) ==0, ] #if no empty rows, next step not needed
#cover.Biodrop.gb<-cover.gb[rowSums(cover.gb[, (1:157)])  >0 ]#remove empty rows

#if needed, relativize by row or column or calculate presence/absence
cover.rowsums <- rowSums(cover.gb.early [1:156])
cover.relrow <- data.frame(cover.gb.early /cover.rowsums)
#cover.colmax<-sapply(cover.Bio ,max)
#cover.relcolmax <- data.frame(sweep(cover.Bio ,2,cover.colmax,'/'))
#cover.pa <- cover.Bio %>% mutate_each(funs(ifelse(.>0,1,0)), 1:57)


######################
#NMDS
######################
#make bray-curtis dissimilarity matrix
gb.bcd.early <- vegdist(cover.gb.early)
gb.mds.early<-metaMDS(cover.gb.early, dist="bray", trace = TRUE, autotransform=T, trymax=100, k=4)
#gb.mds.early<-metaMDS(gb.bcd.early , trace = TRUE, autotransform=T, trymax=100, k=3) #runs several with different starting configurations
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
#legend("bottomright",legend=levels(as.factor(cols1$thermal)), col=Lcols, pch=15, cex=0.9,inset=0.1,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)
legend("topright",legend=levels(as.factor(shapes.e$year)), col="black", pch=Lshapes.e, cex=0.9,inset=0.1,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)

#make nicer plot colored based on thermal, shapes on pre/post fire
#help(ordiplot)
#first, set colors and shapes
cols2.e<- all.data.early %>% dplyr::select(thermal) %>% mutate(color = "black", 
                                                              color = ifelse(thermal == "cool", "purple", 
                                                                             ifelse(thermal=="moderate", "orange", 
                                                                                    ifelse(thermal=="very cool", "blue",
                                                                                           ifelse(thermal=="very warm", "red", color))))) #colors based on thermal group
Lcols.e2 <- rep(c("purple", "orange", "blue","red","black")) #colors for the legend
shapes.e <- all.data.early %>% dplyr::select(year) %>%
  mutate(shape = 0, shape = ifelse(year == "2005", 16, ifelse(year=="2006", 17, ifelse(year=="2007", 18,
                                                                                       ifelse(year=="2008", 15, shape))))) #shapes based on grazing 
Lshapes.e <-rep(c(16,17,18,15))#shapes for legend
#make the plot
gb.plot.e <- ordiplot(gb.mds.early, choices=c(1,2), type = "none")   #Set up the plot
points(spscoresall.gb.e$NMDS1,spscoresall.gb.e$NMDS2,col=cols2.e$color,pch=shapes.e$shape) 
#plot(envvec.nms.gb, col="green")
#text(gb.mds, display = "species", cex=0.5, col="grey30") #label species
legend("bottomright",legend=levels(as.factor(cols2.e$thermal)), col=Lcols.e2, pch=15, cex=0.9,inset=0.1,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)
legend("topright",legend=levels(as.factor(shapes.e$year)), col="black", pch=Lshapes.e, cex=0.9,inset=0.1,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)

###MODERATE######
###early years####
#create wide data, first filter so year is 2005 to 2008
mod.data.early<- all.dat %>% dplyr::select(-X1, -spcode) %>% filter(year>2004 & year < 2009) %>% filter(thermal=="moderate") %>% spread(spname, cover)
mod.data.early[is.na(mod.data.early)] <- 0 #replace NAs with 0 (species not counted in plots have NAs when wide dataset created)

#Import environmental data - NADP pinnacles deposition data and NOAA San Jose temperature data
env <- read_csv(paste(datpath_clean, "/NTN-CA66-deposition.csv", sep=""))
env<-env %>% rename(year=yr) %>% filter(year>2004 & year < 2009)
temp <- read_csv(paste(datpath_clean, "/san_jose_clim.csv", sep=""),skip=9)
temp <- temp %>% rename(year=DATE, temp=TAVG) %>% filter(year>2004 & year <2009) %>% dplyr::select(year, temp)
env<-merge(env,temp)

#merge env to match full data
mod.env.early<-merge(env,mod.data.early) %>% dplyr::select(NH4, NO3, totalN, ppt, temp)

#check count of thermal factor
mod.data.early %>% 
  group_by(graze, burn, transect) %>%
  summarise(no_rows = length(thermal)) 

#check count of graze and burn factor
mod.data.early %>% 
  group_by(graze, burn) %>%
  summarise(no_rows = length(graze))

#see which species are not present using colmax
cover.colmax<-sapply(mod.data.early ,max)
i<-cover.colmax != 0
mod.data.early.nozero<-mod.data.early[, i]

plotnames<-mod.data.early.nozero[,1]
cover.mod.early<- mod.data.early.nozero %>% dplyr::select(-c(1:6)) #wide data with ID columns removed, only species/cover for NMDS
rownames(cover.mod.early)<-plotnames

#check for empty rows
cover.Biodrop.mod.early<-cover.mod.early[rowSums(cover.mod.early[, (1:69)]) ==0, ] #if no empty rows, next step not needed
#cover.Biodrop.gb<-cover.gb[rowSums(cover.gb[, (1:157)])  >0 ]#remove empty rows

#if needed, relativize by row or column or calculate presence/absence
cover.rowsums.me <- rowSums(cover.mod.early [1:69])
cover.relrow.me <- data.frame(cover.mod.early /cover.rowsums.me)
#cover.colmax<-sapply(cover.mod.early ,max)
#cover.relcolmax <- data.frame(sweep(cover.Bio ,2,cover.colmax,'/'))
#cover.pa <- cover.Bio %>% mutate_each(funs(ifelse(.>0,1,0)), 1:57)
covercheck<-as.data.frame(cover.rowsums.me)

######################
#NMDS
######################
#make bray-curtis dissimilarity matrix
mod.bcd.early <- vegdist(cover.relrow.me)

mod.mds.early<-metaMDS(mod.bcd.early , trace = TRUE, autotransform=F, trymax=100, k=5) #runs several with different starting configurations
mod.mds.early #solution did not converge after 100 tries
summary(mod.mds.early)
dimcheckMDS(cover.relrow.me, distance = "bray", k = 8, trymax = 20, autotransform = F) #check for optimal dimensions - choose when starts to flatten out

#quick plot of results
stressplot(mod.mds.early, mod.bcd.early) #stressplot to show fit
ordiplot(mod.mds.early)

#store scores in new dataframe
mod.data.early <- mod.data.early %>%
  unite(treatment, c(burn, graze), remove=FALSE, sep = " ")
spscores1.mod.e<-scores(mod.mds.early,display="sites",choices=1)
spscores2.mod.e<-scores(mod.mds.early,display="sites",choices=2)
year<-mod.data.early$year
treatment<-mod.data.early$treatment
spscoresall.mod.e<-data.frame(year,treatment,spscores1.mod.e,spscores2.mod.e)


############################
#Indicator species for early MODERATE#
############################
mod.trt.early <- mod.data.early %>%
  unite(treatment, c(burn, graze), remove=FALSE, sep = " ")
mod_trt_isa_early = multipatt(cover.mod.early, mod.trt.early$treatment, control=how(nperm=999))
summary(mod_trt_isa_early)

#make nicer plot colored based on burn/graze, shapes on year
#moderate
species.e<-as.data.frame(mod.mds.early$species) 
species.e$name<-row.names(species.e)
#store custom transparent color 
mycol1<- rgb(0, 0, 255, max = 255, alpha = 0, names = "blue50")
spc.e<- species.e %>% mutate(color = mycol1, 
                         color = ifelse(name == "Trifolium depauperatum", name=="Crassula connata", name=="Calandrinia ciliata", name=="Poa secunda ssp. secunda", name=="Festuca bromoides",name=="Koeleria macrantha", name=="Galium aparine", name == "Rigiopappus leptoclodus", "red", 
                                        ifelse(name=="Festuca myuros", name=="Layia gaillardiodes", "orange",
                                               ifelse(name=="Epilobium sp." | name=="Chlorogalum pomeridianum" |name=="Sanicula bipinnatifida"|name=="Trifolium willdenovii",name=="Triteleia laxa", name=="Allium serra", "forestgreen", 
                                                      ifelse(name=="Plantago erecta"|name=="Lasthenia californica"|name=="Aphanes occidentalis"|name=="Erodium cicutarium"|name=="Gilia tricolor"|name=="Lepidium nitidum", name=="Astragalus gambellianus",name=="Hemizonia congesta", name=="Castilleja densiflora", name=="Microseris douglasii", "magenta", 
                                                             ifelse(name=="Brodiaea spp."|name=="Hordeum murinum ssp. leporinum"|name=="Eschscholzia californica"|name=="Muilla maritima", "blue",
                                                                    ifelse(name=="Festuca perennis", "green", color)))))))

cols1.e<-mod.data.early %>% dplyr::select(burn, graze) %>% mutate(color = "forestgreen", 
                                                                  color = ifelse(burn == "burned" & graze=="grazed", "red",
                                                                                 ifelse(burn=="burned" & graze =="ungrazed", "orange",
                                                                                        ifelse(burn=="unburned" & graze=="graze", "purple", color)))) #colors based on burn trt
Lcols.e <- rep(c("orange", "Red", "forestgreen","yellow4")) #colors for the legend
shapes.e <- mod.data.early %>% dplyr::select(year) %>%
  mutate(shape = 0, shape = ifelse(year == "2005", 16, ifelse(year=="2006", 17, ifelse(year=="2007", 8,
                                                                                       ifelse(year=="2008", 15, shape))))) #shapes based on grazing 
Lshapes.e <-rep(c(16,17,8,15))#shapes for legend
#make the plot
mod.plot.e <- ordiplot(mod.mds.early, choices=c(1,2), ylim=c(-0.4, 0.4), type = "none")   #Set up the plot
points(spscoresall.mod.e$NMDS1,spscoresall.mod.e$NMDS2,col=cols1.e$color,pch=shapes.e$shape) 
#plot(envvec.nms.gb, col="green")
text(species.e$MDS1,species.e$MDS2, cex=0.9, col=spc.e$color, label=species.e$name) #label species
#legend("bottomright",legend=levels(as.factor(cols1$thermal)), col=Lcols, pch=15, cex=0.9,inset=0.1,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)
legend("topright",legend=levels(as.factor(shapes.e$year)), col="black", pch=Lshapes.e, cex=0.9,inset=0.1,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)


#create plot in ggplot 
fig1b<-ggplot(subset(spscoresall.mod.e, year==2005), aes(x=NMDS1, y=NMDS2, col=treatment, shape=as.factor(year)))+
  geom_point(cex=2)+
  ggtitle("b) 2005")+
  xlim(-0.4,0.4)+
  ylim(-0.4,0.4)+
  scale_shape_manual(values=c(16),guide = guide_legend(title = "Year"))+
  scale_color_manual(values=c("red", "orange", "forestgreen"), guide = guide_legend(title = "Treatment"), #change legend title
                     labels=c("Burned & Grazed", "Burned & Ungrazed", "Unburned & Ungrazed"))+ #change labels in the legend)+
  theme_bw()+
  theme(plot.title = element_text(color="black", size=14, face="bold.italic"))+
  theme(legend.position="none")
fig1b

fig1c<-ggplot(subset(spscoresall.mod.e, year==2006), aes(x=NMDS1, y=NMDS2, col=treatment, shape=as.factor(year)))+
  geom_point(cex=2)+
  ggtitle("c) 2006")+
  xlim(-0.4,0.4)+
  ylim(-0.4,0.4)+
  scale_shape_manual(values=c(17),guide = guide_legend(title = "Year"))+
  scale_color_manual(values=c("red", "orange", "forestgreen"), guide = guide_legend(title = "Treatment"), #change legend title
                     labels=c("Burned & Grazed", "Burned & Ungrazed", "Unburned & Ungrazed"))+ #change labels in the legend)+
  theme_bw()+
  theme(plot.title = element_text(color="black", size=14, face="bold.italic"))+
  theme(legend.position="none")
fig1c
fig1c+geom_text(subset(spscoresall.mod.e, year==2006), mapping=aes(x=NMDS1, y=NMDS2, label=mod.2006$transect), cex=4)

fig1d<-ggplot(subset(spscoresall.mod.e, year==2007), aes(x=NMDS1, y=NMDS2, col=treatment, shape=as.factor(year)))+
  geom_point(cex=2)+
  ggtitle("d) 2007")+
  xlim(-0.4,0.4)+
  ylim(-0.4,0.4)+
  scale_shape_manual(values=c(8),guide = guide_legend(title = "Year"))+
  scale_color_manual(values=c("red", "orange", "forestgreen"), guide = guide_legend(title = "Treatment"), #change legend title
                     labels=c("Burned & Grazed", "Burned & Ungrazed", "Unburned & Ungrazed"))+ #change labels in the legend)+
  theme_bw()+
  theme(legend.position="none")+
  theme(plot.title = element_text(color="black", size=14, face="bold.italic"))
fig1d

fig1e<-ggplot(subset(spscoresall.mod.e, year==2008), aes(x=NMDS1, y=NMDS2, col=treatment, shape=as.factor(year)))+
  geom_point(cex=2)+
  ggtitle("e) 2008")+
  xlim(-0.4,0.4)+
  ylim(-0.4,0.4)+
  scale_shape_manual(values=c(15),guide = guide_legend(title = "Year"))+
  scale_color_manual(values=c("red", "orange","forestgreen"), guide = guide_legend(title = "Treatment"), #change legend title
                     labels=c("Burned & Grazed", "Burned & Ungrazed", "Unburned & Ungrazed"))+ #change labels in the legend)+
  theme_bw()+
  theme(legend.position="none")+
  theme(plot.title = element_text(color="black", size=14, face="bold.italic"))
#theme(legend.position="bottom", legend.title=element_text(size=11), legend.text=element_text(size=10), axis.text=element_text(size=8), axis.title=element_text(size=11))+
#theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))
fig1e


fig1f<-ggplot(spscoresall.mod.e, aes(x=NMDS1, y=NMDS2, col=mod.data.early$treatment, shape=as.factor(mod.data.early$year)))+
  geom_point(cex=2)+
  ggtitle("f)")+
  xlim(-0.4,0.4)+
  ylim(-0.4,0.4)+
  scale_shape_manual(values=c(16,17,8,15),guide = guide_legend(title = "Year:"))+
  scale_color_manual(values=c("red", "orange","forestgreen"), guide = guide_legend(title = "Treatment:"), #change legend title
                     labels=c("Burned & Grazed", "Burned & Ungrazed", "Unburned & Ungrazed"))+ #change labels in the legend)+
  theme_bw()+
  #theme(legend.position="none")+
  theme(plot.title = element_text(color="black", size=14, face="bold.italic"))+
  theme(legend.position="bottom", legend.title=element_text(size=10), legend.text=element_text(size=8), axis.text=element_text(size=8), axis.title=element_text(size=11))+
  theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))
fig1f

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(fig1f)

##test for differences in treatment by year#######
mod.2005<-subset(mod.data.early, year==2005)
cover.2005<-mod.2005 %>% dplyr::select(-quadratNew,-treatment,-thermal,-burn,-graze, -year, -transect)
cover.rowsums.05 <- rowSums(cover.2005 [1:156])
cover.relrow.05 <- data.frame(cover.2005/cover.rowsums.05)
mod.bcd.05 <- vegdist(cover.relrow.05)
permanova05<-adonis(cover.relrow.05~mod.2005$treatment,perm=100, method="bray")
permanova05
pairwise.perm.manova(mod.bcd.05,mod.2005$treatment, nperm=100) #all three treatments differ in 2005

mod.2006<-subset(mod.data.early, year==2006)
cover.2006<-mod.2006 %>% dplyr::select(-quadratNew,-treatment,-thermal,-burn,-graze, -year, -transect)
cover.rowsums.06 <- rowSums(cover.2006 [1:156])
cover.relrow.06 <- data.frame(cover.2006/cover.rowsums.06)
mod.bcd.06 <- vegdist(cover.relrow.06)
permanova06<-adonis(cover.relrow.06~mod.2006$treatment, perm=100, method="bray")
permanova06
pairwise.perm.manova(mod.bcd.06,mod.2006$treatment, nperm=100) #all three treatments differ in 2006

mod.2007<-subset(mod.data.early, year==2007)
cover.2007<-mod.2007 %>% dplyr::select(-quadratNew,-treatment,-thermal,-burn,-graze, -year, -transect)
cover.rowsums.07 <- rowSums(cover.2007 [1:156])
cover.relrow.07 <- data.frame(cover.2007/cover.rowsums.07)
mod.bcd.07 <- vegdist(cover.relrow.07)
permanova07<-adonis(cover.relrow.07~mod.2007$treatment, perm=100, method="bray")
permanova07
pairwise.perm.manova(mod.bcd.07,mod.2007$treatment, nperm=100) #ungrazed communities are same, both differ from grazed

mod.2008<-subset(mod.data.early, year==2008)
cover.2008<-mod.2008 %>% dplyr::select(-quadratNew,-treatment,-thermal,-burn,-graze, -year, -transect)
cover.rowsums.08 <- rowSums(cover.2008 [1:156])
cover.relrow.08 <- data.frame(cover.2008/cover.rowsums.08)
mod.bcd.08 <- vegdist(cover.relrow.08)
permanova08<-adonis(cover.2008~mod.2008$treatment, perm=100, method="bray")
permanova08
pairwise.perm.manova(mod.bcd.08,mod.2008$treatment, nperm=100) #ungrazed communities are same, both differ from grazed


#####################
#successional vectors on summarized MODERATE data for burn (2005-2008)
####################
mod_yr_burn<-all.dat %>% group_by(thermal, burn, graze, year, spname) %>% filter(thermal == "moderate", year>2004 & year<2009) %>% 
  summarize(mean=mean(cover))%>% arrange(burn)%>%  arrange(graze)%>% arrange(year)%>%
  spread(spname, mean) 
mod_yr_burn[is.na(mod_yr_burn)] <- 0 

cover.yr <- mod_yr_burn %>% ungroup %>% dplyr::select(-thermal,-burn,-graze, -year)
cover.rowsums <- rowSums(cover.yr [1:156])
cover.relrow <- data.frame(cover.yr /cover.rowsums)

#merge env to match full data
#yr.env<-merge(env,mod_yr) %>% dplyr::select(NH4, NO3, totalN, ppt, temp)

#make bray-curtis dissimilarity matrix
vec.bcd <- vegdist(cover.relrow)
dimcheckMDS(cover.relrow)#check for optimal dimensions

#NMDS 
vec.mds<-metaMDS(cover.relrow, trace = TRUE, autotransform = F, trymax=100, k=4) #runs several with different starting configurations
#trace= TRUE will give output for step by step what its doing
#default is 2 dimensions, can put k=4 for 4 dimensions
vec.mds #solution converged after 20 tries
summary(vec.mds)


#quick plot of results
stressplot(vec.mds, vec.bcd) #stressplot to show fit
ordiplot(vec.mds)

#overlay environmental variables, only showing significant drivers
#envvec.vec<-envfit(vec.mds,yr.env, na.rm=TRUE)
#envvec.vec
#plot(vec.mds)
#plot(envvec.vec, p.max=0.05) #add vectors to previous ordination

#store scores in new dataframe
spscores1.vec<-scores(vec.mds,display="sites",choices=1)
spscores2.vec<-scores(vec.mds,display="sites",choices=2)
year<-mod_yr_burn$year
burn<-mod_yr_burn$burn
spscoresall.vec<-data.frame(burn,year,spscores1.vec,spscores2.vec)

#make plot to show successional vectors
#first, set colors and shapes
mod_yr_burn <- mod_yr_burn %>%
  unite(treatment, c(burn, graze), remove=FALSE, sep = " ")
cols.yr<- mod_yr_burn %>% dplyr::select(burn, graze) %>% mutate(color = "forestgreen", 
                                                                color = ifelse(burn == "burned" & graze=="grazed", "red",
                                                                               ifelse(burn=="burned" & graze =="ungrazed", "orange",
                                                                                      ifelse(burn=="unburned" & graze=="graze", "purple", color)))) #colors based on burn trt
Lcols.yr <- rep(c( "red", "orange","forestgreen")) #colors for the legend
#shapes.yr <- dat_yr%>% dplyr::select(year) %>%
#mutate(shape = 1, shape = ifelse(year == "2004", 8, ifelse(year>="2005" & year<"2014", 16, ifelse(year>="2014", 15,shape)))) #shapes based on year 
#shapes.yr<-shapes.yr %>% mutate(time="Pre-Fire", 
#time= ifelse(year==2004, "Fire",
#ifelse(year>=2005&year<2014, "Post-Fire",
#ifelse(year>=2014, "2014 and after",time)))) 
#Lshapes <-rep(c(15,8,16,1))#shapes for legend
#make the plot
vec.plot <- ordiplot(vec.mds, choices=c(1,2), type = "none")   #Set up the plot
points(spscoresall.vec$NMDS1,spscoresall.vec$NMDS2, col=cols.yr$color,pch=19) 
#plot(envvec.vec, p.max=0.05, col="green")
ordiarrows(vec.mds, groups=mod_yr_burn$treatment, order.by=mod_yr_burn$year, label=F, col=Lcols.yr)
#text(spp.mds, display = "species", cex=0.5, col="grey30") #label species
legend("topleft",legend=levels(as.factor(mod_yr_burn$treatment)), col=Lcols.yr, pch=15, cex=0.9,inset=0.07,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)
#legend("topright",legend=levels(as.factor(shapes$time)), col="black", pch=Lshapes, cex=0.9,inset=0.07,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)

############################
#Indicator species for early MODERATE#
############################
mod_vec_isa_early = multipatt(cover.relrow, mod_yr_burn$treatment, control=how(nperm=999))
summary(mod_vec_isa_early)

############################
###Make vector plot again in ggplot
###########################

#to plot indicator species (from above) on plot
species.e<-as.data.frame(vec.mds$species)
species.e$name<-row.names(species.e)
spc.e<- species.e %>% filter(name == "Trifolium depauperatum"| name=="Crassula connata"| name=="Calandrinia ciliata"| name == "Rigiopappus leptoclodus" |name=="Poa secunda ssp. secunda"| name=="Festuca bromoides"|name=="Koeleria macrantha"| name=="Galium aparine"| name == "Rigiopappus leptoclodus" |name=="Festuca myuros"| name=="Layia gaillardiodes"| name=="Silene gallica"|name=="Athysanus pusilus"|name=="Sisyrinchium bellum"|name=="Epilobium sp."| name=="Chlorogalum pomeridianum"|name=="Sanicula bipinnatifida"|name=="Lessingia micradenia glabratai"|name=="Triteleia laxa"| name=="Allium serra"|name=="Plantago erecta"|name=="Lasthenia californica"|name=="Aphanes occidentalis"|name=="Erodium cicutarium"|name=="Gilia tricolor"|name=="Lepidium nitidum"| name=="Hemizonia congesta"| name=="Castilleja densiflora"| name=="Microseris douglasii"|name=="Agoseris heterophylla"|name=="Brodiaea spp."|name=="Hordeum murinum ssp. leporinum"|name=="Muilla maritima"|name=="Festuca perennis")

vec1<-ggplot(spscoresall.vec, aes(x=NMDS1, y=NMDS2))+
  geom_point(cex=3, aes(shape=as.factor(mod_yr_burn$year), col=mod_yr_burn$treatment))+
  ggtitle("a)")+
  scale_color_manual(values=c("red", "orange","forestgreen"), guide = guide_legend(title = "Treatment"), #change legend title
                     labels=c("Burned & Grazed", "Burned & Ungrazed", "Unburned & Ungrazed"))+ #change labels in the legend)+
  scale_shape_manual(values=c(16,17,8,15,18),guide = guide_legend(title = "Year"))+
  geom_path(arrow=arrow(), aes(col=mod_yr_burn$treatment))+
  theme_bw()+
  theme(legend.position="none")+
  theme(plot.title = element_text(color="black", size=14, face="bold.italic"))
  #theme(legend.position=c(0.8,0.8), legend.title=element_text(size=14), legend.text=element_text(size=12), axis.text=element_text(size=16), axis.title=element_text(size=16))+
  #theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))
vec1

vec1<-vec1+geom_text(data=spc.e, mapping=aes(x=MDS1, y=MDS2, label=name), cex=1.5)
vec1

#create layout for panel
lay <- rbind(c(1,1,1,1),
             c(1,1,1,1),
             c(1,1,1,1),
             c(1,1,1,1),
             c(2,2,3,3),
             c(2,2,3,3),
             c(2,2,3,3),
             c(4,4,5,5),
             c(4,4,5,5),
             c(4,4,5,5),
             c(6,6,6,6))
grid.arrange(vec1, fig1b, fig1c, fig1d, fig1e, mylegend, layout_matrix = lay) #put panel together



###MODERATE#####
###late years####
#create wide data, first filter so year is 2005 to 2008
mod.data.late<- all.dat %>% dplyr::select(-X1, -spcode) %>% filter(year>2007 & year < 2013) %>% filter(thermal=="moderate") %>% spread(spname, cover)
mod.data.late[is.na(mod.data.late)] <- 0 #replace NAs with 0 (species not counted in plots have NAs when wide dataset created)

#Import environmental data - NADP pinnacles deposition data and NOAA San Jose temperature data
env <- read_csv(paste(datpath_clean, "/NTN-CA66-deposition.csv", sep=""))
env<-env %>% rename(year=yr) %>% filter(year>2008 & year < 2013)
temp <- read_csv(paste(datpath_clean, "/san_jose_clim.csv", sep=""),skip=9)
temp <- temp %>% rename(year=DATE, temp=TAVG) %>% filter(year>2008 & year <2012) %>% dplyr::select(year, temp)
env<-merge(env,temp)

#merge env to match full data
mod.env.late<-merge(env,mod.data.late) %>% dplyr::select(NH4, NO3, totalN, ppt, temp)

#check count of graze and burn factor
mod.data.late %>% 
  group_by(graze, burn) %>%
  summarise(no_rows = length(graze))

plotnames<-mod.data.late[,1]
cover.mod.late<- mod.data.late %>% dplyr::select(-c(1:6)) #wide data with ID columns removed, only species/cover for NMDS
rownames(cover.mod.late)<-plotnames
##delete outlier row 122 (nearly empty row: has only 2 species with 1% cover... is this a data recording issue?)
mod.data.late <- mod.data.late[-c(122),] 
cover.mod.late<-cover.mod.late[-c(122),] 

#now relativize by row
cover.rowsums.ml <- rowSums(cover.mod.late [1:156])
cover.relrow.ml <- data.frame(cover.mod.late /cover.rowsums.ml)


#check for empty rows
cover.Biodrop.mod.late<-cover.mod.late[rowSums(cover.mod.late[, (1:156)]) ==0, ] #if no empty rows, next step not needed
#cover.Biodrop.gb<-cover.gb[rowSums(cover.gb[, (1:157)])  >0 ]#remove empty rows

######################
#NMDS
######################
#make bray-curtis dissimilarity matrix
mod.bcd.late <- vegdist(cover.relrow.ml)

#check for optimal number of dimensions (optimum is when scree plot flattens out/break in slope, we'd also like the stress to be under 10% or 0.10)
dimcheckMDS(cover.relrow.ml, distance = "bray", k = 8, trymax = 20, autotransform = F)

mod.mds.late<-metaMDS(cover.relrow.ml, trace = TRUE, autotransform=F, trymax=100, k=4) #runs several with different starting configurations
#trace= TRUE will give output for step by step what its doing
#default is 2 dimensions, k=4 for 4 dimensions
mod.mds.late #no solution reached
summary(mod.mds.late)

#quick plot of results
stressplot(mod.mds.late, mod.bcd.late) #stressplot to show fit
ordiplot(mod.mds.late)
orditorp(mod.mds.late,display="sites",cex=0.8,air=0.01)

#store scores in new dataframe
mod.trt.late <- mod.data.late %>%
  unite(treatment, c(burn, graze), remove=FALSE, sep = " ")
spscores1.mod.l<-scores(mod.mds.late,display="sites",choices=1)
spscores2.mod.l<-scores(mod.mds.late,display="sites",choices=2)
treatment<-mod.trt.late$treatment
year<-mod.data.late$year
spscoresall.mod.l<-data.frame(year,treatment,spscores1.mod.l,spscores2.mod.l)

############################
#Indicator species for late MODERATE#
############################
mod_trt_isa_late = multipatt(cover.relrow.ml, mod.trt.late$treatment, control=how(nperm=999))
summary(mod_trt_isa_late)

#make  plot colored based on burn/graze, shapes on year
cols1.l<-mod.data.late %>% dplyr::select(burn, graze) %>% mutate(color = "forestgreen", 
                                                                 color = ifelse(burn == "burned" & graze=="grazed", "red",
                                                                                ifelse(burn=="burned" & graze =="ungrazed", "orange",
                                                                                       ifelse(burn=="unburned" & graze=="graze", "purple", color)))) #colors based on burn trt
Lcols.l <- rep(c("orange", "Red", "forestgreen","yellow4")) #colors for the legend
shapes.l <- mod.data.late %>% dplyr::select(year) %>%
  mutate(shape = 15, shape = ifelse(year=="2009", 0, ifelse(year=="2010",2, 
                                                           ifelse(year=="2011", 5, 
                                                                  ifelse(year=="2012", 6, shape))))) #shapes based on grazing 
Lshapes.l <-rep(c(15,0,2, 5,6))#shapes for legend
#make the plot
mod.plot.l <- ordiplot(mod.mds.late, choices=c(1,2),  type = "none")   #Set up the plot
points(spscoresall.mod.l$NMDS1,spscoresall.mod.l$NMDS2,col=cols1.l$color,pch=shapes.l$shape) 
#plot(envvec.nms.gb, col="green")
#text(gb.mds, display = "species", cex=0.5, col="grey30") #label species
#legend("bottomright",legend=levels(as.factor(cols1$thermal)), col=Lcols, pch=15, cex=0.9,inset=0.1,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)
legend("topright",legend=levels(as.factor(shapes.l$year)), col="black", pch=Lshapes.l, cex=0.9,inset=0.1,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)


###late years, all thermal####
#create wide data, first filter so year is 2005 to 2008
all.data.late<- all.dat %>% dplyr::select(-X1, -spcode) %>% filter(year>2008 & year < 2013) %>% spread(spname, cover)
all.data.late[is.na(all.data.late)] <- 0 #replace NAs with 0 (species not counted in plots have NAs when wide dataset created)

#Import environmental data - NADP pinnacles deposition data and NOAA San Jose temperature data
env <- read_csv(paste(datpath_clean, "/NTN-CA66-deposition.csv", sep=""))
env<-env %>% rename(year=yr) %>% filter(year>2008 & year < 2013)
temp <- read_csv(paste(datpath_clean, "/san_jose_clim.csv", sep=""),skip=9)
temp <- temp %>% rename(year=DATE, temp=TAVG) %>% filter(year>2008 & year <2012) %>% dplyr::select(year, temp)
env<-merge(env,temp)

#merge env to match full data
all.env.late<-merge(env,all.data.late) %>% dplyr::select(NH4, NO3, totalN, ppt, temp)

#check count of thermal factor
all.data.late %>% 
  group_by(thermal, burn, graze) %>%
  summarise(no_rows = length(thermal)) 

#check count of graze and burn factor
all.data.late %>% 
  group_by(graze, burn) %>%
  summarise(no_rows = length(graze))

plotnames<-all.data.late[,1]
cover.gb.late<- all.data.late %>% dplyr::select(-c(1:5)) #wide data with ID columns removed, only species/cover for NMDS
rownames(cover.gb.late)<-plotnames

#check for empty rows
cover.Biodrop.gb.late<-cover.gb.late[rowSums(cover.gb.late[, (1:156)]) ==0, ] #if no empty rows, next step not needed
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
gb.mds.late #no solution reached
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

#make nicer plot colored based on burn/graZE, shapes on YEAR
#help(ordiplot)
#first, set colors and shapes
cols1.l<-all.data.late %>% dplyr::select(burn, graze) %>% mutate(color = "forestgreen", 
                                                                 color = ifelse(burn == "burned" & graze=="grazed", "red",
                                                                                ifelse(burn=="burned" & graze =="ungrazed", "orange",
                                                                                       ifelse(burn=="unburned" & graze=="graze", "purple", color)))) #colors based on burn trt
Lcols.l <- rep(c("orange", "Red", "forestgreen","yellow4")) #colors for the legend
shapes.l <- all.data.late %>% dplyr::select(year) %>%
  mutate(shape = 0, shape = ifelse(year=="2009", 19, ifelse(year=="2010",2, 
                                                           ifelse(year=="2011", 5, shape)))) #shapes based on grazing 
Lshapes.l <-rep(c(19,2,5, 0))#shapes for legend
#make the plot
gb.plot.l <- ordiplot(gb.mds.late, choices=c(1,2), xlim=c(-0.2, 0.2), type = "none")   #Set up the plot
points(spscoresall.gb.l$NMDS1,spscoresall.gb.l$NMDS2,col=cols1.l$color,pch=shapes.l$shape) 
#plot(envvec.nms.gb, col="green")
#text(gb.mds, display = "species", cex=0.5, col="grey30") #label species
#legend("bottomright",legend=levels(as.factor(cols1$thermal)), col=Lcols, pch=15, cex=0.9,inset=0.1,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)
legend("topright",legend=levels(as.factor(shapes.l$year)), col="black", pch=Lshapes.l, cex=0.9,inset=0.1,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)

#make plot colored based on thermal, shapes on year
cols2.l<- all.data.late %>% dplyr::select(thermal) %>% mutate(color = "black", 
                                                              color = ifelse(thermal == "cool", "purple", 
                                                                             ifelse(thermal=="moderate", "orange", 
                                                                                    ifelse(thermal=="very cool", "blue",
                                                                                           ifelse(thermal=="very warm", "red", color))))) #colors based on thermal group
Lcols.l2 <- rep(c("purple", "orange", "blue","red","black")) #colors for the legend
shapes.l <- all.data.late %>% dplyr::select(year) %>%
  mutate(shape = 0, shape = ifelse(year=="2009", 19, ifelse(year=="2010",2, 
                                                           ifelse(year=="2011", 5, shape)))) #shapes based on grazing 
Lshapes.l <-rep(c(19,2,5,0))#shapes for legend
#make the plot
gb.plot.l <- ordiplot(gb.mds.late, choices=c(1,2),xlim=c(-0.2, 0.2), type = "none")   #Set up the plot
points(spscoresall.gb.l$NMDS1,spscoresall.gb.l$NMDS2,col=cols2.l$color,pch=shapes.l$shape) 
#plot(envvec.nms.gb, col="green")
#text(gb.mds, display = "species", cex=0.5, col="grey30") #label species
legend("bottomright",legend=levels(as.factor(cols2.l$thermal)), col=Lcols.l2, pch=15, cex=0.9,inset=0.08,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)
legend("topright",legend=levels(as.factor(shapes.l$year)), col="black", pch=Lshapes.l, cex=0.9,inset=0.08,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)




########################
#####NMDS for 2006######
########################
# Import csv file, transform to wide data, call it data
all.dat <- read_csv(paste(datpath_clean, "/alldatsptrt.csv", sep=""))
all.dat<-as.data.frame(all.dat)
#dat2<-dat %>% mutate(cover=cover+0.00001) #a work around for data with lots of zeros
#create wide data, first filter so year is 2005 to 2012

#2006 data only
dat.06<- all.dat %>% 
  dplyr::select(-X1, -spcode) %>% 
  filter(year == 2006) %>% 
  spread(spname, cover)
dat.06[is.na(dat.06)] <- 0 #replace NAs with 0 (species not counted in plots have NAs when wide dataset created)
#2006 moderate data only
mod.dat.06 <- dat.06 %>%
  filter(thermal == "moderate")

levels(as.factor(dat.06$quadratNew))
levels(as.factor(dat.06$thermal))
levels(as.factor(dat.06$spname)) #any to remove?
levels(as.factor(dat.06$burn))
levels(as.factor(dat.06$graze))

#Import environmental data - NADP pinnacles deposition data and NOAA San Jose temperature data
env <- read_csv(paste(datpath_clean, "/NTN-CA66-deposition.csv", sep=""))
env<-env %>% 
  rename(year=yr) %>% 
  filter(year == 2006)
temp <- read_csv(paste(datpath_clean, "/san_jose_clim.csv", sep=""),skip=9)
temp <- temp %>% 
  rename(year=DATE, temp=TAVG) %>% 
  filter(year == 2006) %>% 
  dplyr::select(year, temp)
env<-merge(env,temp)

#merge env to match full data
all.env<-merge(env,dat.06) %>% 
  dplyr::select(NH4, NO3, totalN, ppt, temp)

#check count of thermal factor
dat.06 %>% 
  group_by(thermal) %>%
  summarise(no_rows = length(thermal)) 

#check count of graze factor
dat.06 %>% 
  group_by(graze) %>%
  summarise(no_rows = length(graze))

#check count of burn factor
dat.06 %>% 
  group_by(burn) %>%
  summarise(no_rows = length(burn))

str(dat.06)

dat.06$ID <- seq.int(nrow(dat.06))
plotnames<-dat.06[,1]
cover.gb.06<- dat.06 %>% 
  dplyr::select(-c(1:5,161)) #wide data with ID columns removed, only species/cover for NMDS
rownames(cover.gb.06)<-plotnames

mod.dat.06$ID <- seq.int(nrow(mod.dat.06))
plotnames<-mod.dat.06[,1]
cover.gb.mod06<- mod.dat.06 %>% 
  dplyr::select(-c(1:5,161)) #wide data with ID columns removed, only species/cover for NMDS
rownames(cover.gb.mod06)<-plotnames

#check for empty rows
cover.Biodrop.gb.06<-cover.gb.06[rowSums(cover.gb.06[, (1:155)]) ==0, ]
cover.Biodrop.gb.mod06<-cover.gb.mod06[rowSums(cover.gb.mod06[, (1:155)]) ==0, ]#if no empty rows, next step not needed
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
gb.bcd.06 <- vegdist(cover.gb.06)
gb.bcd.mod06 <- vegdist(cover.gb.mod06) 

#quick run to check out NMS ordination
gb.mds0 <-isoMDS(gb.bcd.06) #runs nms only once
gb.mds0  #by default 2 dimensions returned, stress is 6.4, converged
ordiplot(gb.mds0) #ugly

gb.mds0 <-isoMDS(gb.bcd.mod06) #runs nms only once
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
gb.mds.06<-metaMDS(cover.gb.06, trace = TRUE, autotransform=T, trymax=100, k=6)
gb.mds.mod06<-metaMDS(cover.gb.mod06, trace = TRUE, autotransform=T, trymax=100, k=6)#runs several with different starting configurations
#trace= TRUE will give output for step by step what its doing
#default is 2 dimensions, can put k=4 for 4 dimensions
gb.mds.06 #solution did not converge after 100 tries
summary(gb.mds.06)

gb.mds.mod06 #solution did not converge after 100 tries
summary(gb.mds.mod06)

#quick plot of results
stressplot(gb.mds.06, gb.bcd.06) #stressplot to show fit
ordiplot(gb.mds.06)

stressplot(gb.mds.mod06, gb.bcd.mod06) #stressplot to show fit
ordiplot(gb.mds.mod06)

#overlay environmental variables on full data
envvec.nms<-envfit(gb.mds.06,all.env, na.rm=TRUE)
envvec.nms
plot(gb.mds.06)
plot(envvec.nms) #add vectors to previous ordination

#store scores in new dataframe
spscores1<-scores(gb.mds.06,display="sites",choices=1)
spscores2<-scores(gb.mds.06,display="sites",choices=2)
tplots<-dat.06[,4]
tplot_levels<-levels(tplots)
spscoresall<-data.frame(tplots,spscores1,spscores2)

spscores1<-scores(gb.mds.mod06,display="sites",choices=1)
spscores2<-scores(gb.mds.mod06,display="sites",choices=2)
tplots<-mod.dat.06[,4]
tplot_levels<-levels(tplots)
spscoresall<-data.frame(tplots,spscores1,spscores2)
#make nicer plot colored based burn and graze treatment in 2006
#help(ordiplot)
#first, set colors and shapes
cols1<- dat.06 %>% 
  dplyr::select(burn) %>% 
  mutate(color = "black", color = ifelse(burn == "burned", "red", 
                                  ifelse(burn =="unburned", "black", color))) #colors based on thermal group
Lcols <- rep(c("Red", "Black")) #colors for the legend
shapes <- dat.06 %>% 
  dplyr::select(graze) %>%
  mutate(shape = 17, shape = ifelse(graze == "grazed", 8, 
                            ifelse(graze == "ungrazed", 17,shape))) #shapes based on year 
Lshapes <-rep(c(8,17))#shapes for legend

#moderate 2006
cols1<- mod.dat.06 %>% 
  dplyr::select(burn) %>% 
  mutate(color = "black", color = ifelse(burn == "burned", "red", 
                                         ifelse(burn =="unburned", "black", color))) #colors based on thermal group
Lcols <- rep(c("Red", "Black")) #colors for the legend
shapes <- mod.dat.06 %>% 
  dplyr::select(graze) %>%
  mutate(shape = 17, shape = ifelse(graze == "grazed", 8, 
                                    ifelse(graze == "ungrazed", 17,shape))) #shapes based on year 
Lshapes <-rep(c(8,17))#shapes for legend
#shapes<-shapes %>% 
#  mutate(time="Pre-Fire", 
#         time= ifelse(year==2004, "Fire",
#               ifelse(year>=2005&year<2014, "Post-Fire",
#              ifelse(year>=2014, "2014 and after",time)))) 


#make the plot
bio.plot <- ordiplot(gb.mds.mod06, choices=c(1,2), type = "none", xlim = c(-3,3), ylim = c(-1,1))   #Set up the plot
points(spscoresall$NMDS1,spscoresall$NMDS2,col=cols1$color,pch=shapes$shape) 
#plot(envvec.nms, col="green")
#text(gb.mds.06, display = "species", cex=0.5, col="grey30") #label species
legend("topleft",legend=levels(as.factor(cols1$burn)), col=Lcols, pch=15, cex=0.9,inset=0.1,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)
legend("bottomleft",legend=levels(as.factor(shapes$graze)), col="black", pch=Lshapes, cex=0.9,inset=0.1,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)

#make nicer plot colored based on thermal
#help(ordiplot)
#first, set colors and shapes
cols1<- dat.06 %>% 
  dplyr::select(thermal) %>% 
  mutate(color = "black", 
         color = ifelse(thermal == "cool", "purple", 
                 ifelse(thermal=="very cool", "blue",
                 ifelse(thermal=="moderate", "orange", 
                 ifelse(thermal=="warm", "yellow",     
                 ifelse(thermal=="very warm", "red", color))))))
Lcols <- rep(c("Purple", "Orange", "Blue", "Red", "Yellow"))

shapes <- dat.06 %>% 
  dplyr::select(graze) %>%
  mutate(shape = 17, shape = ifelse(graze == "grazed", 8, 
                                    ifelse(graze == "ungrazed", 17,shape))) #shapes based on year 
#shapes<-shapes %>% 
#  mutate(time="Pre-Fire", 
#         time= ifelse(year==2004, "Fire",
#               ifelse(year>=2005&year<2014, "Post-Fire",
#              ifelse(year>=2014, "2014 and after",time)))) 

Lshapes <-rep(c(8,17))#shapes for legend
#make the plot
bio.plot <- ordiplot(gb.mds.06, choices=c(1,2), type = "none")   #Set up the plot
points(spscoresall$NMDS1,spscoresall$NMDS2,col=cols1$color,pch=shapes$shape) 
plot(envvec.nms, col="green")
#text(gb.mds.06, display = "species", cex=0.5, col="grey30") #label species
legend("topleft",legend=levels(as.factor(cols1$thermal)), col=Lcols, pch=15, cex=0.9,inset=0.1,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)
legend("bottomleft",legend=levels(as.factor(shapes$graze)), col="black", pch=Lshapes, cex=0.9,inset=0.1,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)

#moderate
species<-as.data.frame(gb.mds.mod06$species) 
species$name<-row.names(species)
#store custom transparent color 
mycol1<- rgb(0, 0, 255, max = 255, alpha = 0, names = "blue50")
spc<- species %>% mutate(color = mycol1, 
         color = ifelse(name == "Rigiopappus leptoclodus", "red", 
                ifelse(name=="Festuca myuros", "orange",
                ifelse(name=="Festuca perennis" | name=="Sanicula bipinnatifida" |name=="Trifolium willdenovii"|name=="Triteleia laxa", "black", 
                ifelse(name=="Plantago erecta"|name=="Lasthenia californica"|name=="Aphanes occidentalis"|name=="Erodium cicutarium"|name=="Gilia tricolor"|name=="Lepidium nitidum", "magenta", 
                ifelse(name=="Brodiaea spp."|name=="Hordeum murinum ssp. leporinum"|name=="Eschscholzia californica"|name=="Muilla maritima", "blue",color))))))

bio.plot <- ordiplot(gb.mds.mod06, choices=c(1,2), type = "none")   #Set up the plot
points(spscoresall$NMDS1,spscoresall$NMDS2,col=cols1$color,pch=shapes$shape) 
text(species$MDS1,species$MDS2, cex=0.9, col=spc$color, label=species$name) #label species
legend("topleft",legend=levels(as.factor(cols1$burn)), col=Lcols, pch=15, cex=0.9,inset=0.1,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)
legend("bottomleft",legend=levels(as.factor(shapes$graze)), col="black", pch=Lshapes, cex=0.9,inset=0.1,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)

############################
#Indicator species for 2006#
############################
trt.dat.06 <- dat.06 %>%
  unite(treatment, c(burn, graze), remove=FALSE, sep = " ")
library(indicspecies)
trt_isa_06 = multipatt(cover.gb.06, trt.dat.06$treatment, control=how(nperm=999))
summary(trt_isa_06)

############################
#Indicator species for 2006 MODERATE#
############################
mod.trt.dat.06 <- mod.dat.06 %>%
  unite(treatment, c(burn, graze), remove=FALSE, sep = " ")
library(indicspecies)
mod_trt_isa_06 = multipatt(cover.gb.mod06, mod.trt.dat.06$treatment, control=how(nperm=999))
summary(mod_trt_isa_06)

########################
#####NMDS for 2005######
########################
# Import csv file, transform to wide data, call it data
all.dat <- read_csv(paste(datpath_clean, "/alldatsptrt.csv", sep=""))
all.dat<-as.data.frame(all.dat)
#dat2<-dat %>% mutate(cover=cover+0.00001) #a work around for data with lots of zeros
#create wide data, first filter so year is 2005 to 2012
dat.05<- all.dat %>% 
  dplyr::select(-X1, -spcode) %>% 
  filter(year == 2005) %>% 
  spread(spname, cover)
dat.05[is.na(dat.05)] <- 0 #replace NAs with 0 (species not counted in plots have NAs when wide dataset created)

mod.dat.05 <- dat.05 %>%
  filter(thermal == "moderate")

levels(as.factor(dat.05$quadratNew))
levels(as.factor(dat.05$thermal))
levels(as.factor(dat.05$spname)) #any to remove?
levels(as.factor(dat.05$burn))
levels(as.factor(dat.05$graze))

#Import environmental data - NADP pinnacles deposition data and NOAA San Jose temperature data
env <- read_csv(paste(datpath_clean, "/NTN-CA66-deposition.csv", sep=""))
env<-env %>% 
  rename(year=yr) %>% 
  filter(year == 2005)
temp <- read_csv(paste(datpath_clean, "/san_jose_clim.csv", sep=""),skip=9)
temp <- temp %>% 
  rename(year=DATE, temp=TAVG) %>% 
  filter(year == 2005) %>% 
  dplyr::select(year, temp)
env<-merge(env,temp)

#merge env to match full data
all.env<-merge(env,dat.05) %>% 
  dplyr::select(NH4, NO3, totalN, ppt, temp)

#check count of thermal factor
dat.05 %>% 
  group_by(thermal) %>%
  summarise(no_rows = length(thermal)) 

#check count of graze factor
dat.05 %>% 
  group_by(graze) %>%
  summarise(no_rows = length(graze))

#check count of burn factor
dat.05 %>% 
  group_by(burn) %>%
  summarise(no_rows = length(burn))

str(dat.05)

dat.05$ID <- seq.int(nrow(dat.05))
plotnames<-dat.05[,1]
cover.gb.05<- dat.05 %>% 
  dplyr::select(-c(1:5,162)) #wide data with ID columns removed, only species/cover for NMDS
rownames(cover.gb.05)<-plotnames

mod.dat.05$ID <- seq.int(nrow(mod.dat.05))
plotnames<-mod.dat.05[,1]
cover.gb.mod05<- mod.dat.05 %>% 
  dplyr::select(-c(1:5,162)) #wide data with ID columns removed, only species/cover for NMDS
rownames(cover.gb.mod05)<-plotnames

#check for empty rows
cover.Biodrop.gb.05<-cover.gb.05[rowSums(cover.gb.05[, (1:156)]) ==0, ] #if no empty rows, next step not needed
cover.Biodrop.gb.mod05<-cover.gb.mod05[rowSums(cover.gb.mod05[, (1:156)]) ==0, ]
#cover.Biodrop.gb<-cover.gb[rowSums(cover.gb[, (1:157)])  >0 ]#remove empty rows

#if needed, relativize by row or column or calculate presence/absence
#cover.rowsums <- rowSums(cover.Bio [1:157])
#cover.relrow <- data.frame(cover.Bio /cover.rowsums)
#cover.colmax<-sapply(cover.Bio ,max)
#cover.relcolmax <- data.frame(sweep(cover.Bio ,2,cover.colmax,'/'))
#cover.pa <- cover.Bio %>% mutate_each(funs(ifelse(.>0,1,0)), 1:57)

gb.bcd.05 <- vegdist(cover.gb.05)
gb.bcd.mod05 <- vegdist(cover.gb.mod05)

#quick run to check out NMS ordination
gb.mds0 <-isoMDS(gb.bcd.05) #runs nms only once
gb.mds0  #by default 2 dimensions returned, stress is 6.4, converged
ordiplot(gb.mds0) #ugly

gb.mds0 <-isoMDS(gb.bcd.mod05) #runs nms only once
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
gb.mds.05<-metaMDS(cover.gb.05, trace = TRUE, autotransform=T, trymax=100, k=6)
gb.mds.mod05<-metaMDS(cover.gb.mod05, trace = TRUE, autotransform=T, trymax=100, k=6)#runs several with different starting configurations
#trace= TRUE will give output for step by step what its doing
#default is 2 dimensions, can put k=4 for 4 dimensions
gb.mds.05 #solution did not converge after 100 tries
summary(gb.mds.05)

gb.mds.mod05 #solution did not converge after 100 tries
summary(gb.mds.mod05)

#quick plot of results
stressplot(gb.mds.mod05, gb.bcd.mod05) #stressplot to show fit
ordiplot(gb.mds.mod05)

#overlay environmental variables on full data
envvec.nms<-envfit(gb.mds.05,all.env, na.rm=TRUE)
envvec.nms
plot(gb.mds.05)
plot(envvec.nms) #add vectors to previous ordination

#store scores in new dataframe
spscores1<-scores(gb.mds.05,display="sites",choices=1)
spscores2<-scores(gb.mds.05,display="sites",choices=2)
tplots<-dat.05[,4]
tplot_levels<-levels(tplots)
spscoresall<-data.frame(tplots,spscores1,spscores2)

#make nicer plot colored based burn and graze treatment in 2006
#help(ordiplot)
#first, set colors and shapes
cols1<- dat.05 %>% 
  dplyr::select(burn) %>% 
  mutate(color = "black", color = ifelse(burn == "burned", "red", 
                                         ifelse(burn =="unburned", "black", color))) #colors based on thermal group
Lcols <- rep(c("Red", "Black")) #colors for the legend
shapes <- dat.05 %>% 
  dplyr::select(graze) %>%
  mutate(shape = 17, shape = ifelse(graze == "grazed", 8, 
                                    ifelse(graze == "ungrazed", 17,shape))) #shapes based on year 
Lshapes <-rep(c(8,17))#shapes for legend
#shapes<-shapes %>% 
#  mutate(time="Pre-Fire", 
#         time= ifelse(year==2004, "Fire",
#               ifelse(year>=2005&year<2014, "Post-Fire",
#              ifelse(year>=2014, "2014 and after",time)))) 

#moderate
cols1<- mod.dat.05 %>% 
  dplyr::select(burn) %>% 
  mutate(color = "black", color = ifelse(burn == "burned", "red", 
                                         ifelse(burn =="unburned", "black", color))) #colors based on thermal group
Lcols <- rep(c("Red", "Black")) #colors for the legend
shapes <- mod.dat.05 %>% 
  dplyr::select(graze) %>%
  mutate(shape = 17, shape = ifelse(graze == "grazed", 8, 
                                    ifelse(graze == "ungrazed", 17,shape))) #shapes based on year 
Lshapes <-rep(c(8,17))#shapes for legend
#make the plot
bio.plot <- ordiplot(gb.mds.05, choices=c(1,2), type = "none", xlim = c(-3.0,2.0), ylim = c(-1.0,1.5))   #Set up the plot
points(spscoresall$NMDS1,spscoresall$NMDS2,col=cols1$color,pch=shapes$shape) 
#text(gb.mds.06, display = "species", cex=0.5, col="grey30") #label species
legend("topleft",legend=levels(as.factor(cols1$burn)), col=Lcols, pch=15, cex=0.9,inset=0.1,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)
legend("bottomleft",legend=levels(as.factor(shapes$graze)), col="black", pch=Lshapes, cex=0.9,inset=0.1,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)

#make nicer plot colored based on thermal
#help(ordiplot)
#first, set colors and shapes
cols1<- dat.05 %>% 
  dplyr::select(thermal) %>% 
  mutate(color = "black", 
         color = ifelse(thermal == "cool", "purple", 
                        ifelse(thermal=="very cool", "blue",
                               ifelse(thermal=="moderate", "orange", 
                                      ifelse(thermal=="warm", "yellow",     
                                             ifelse(thermal=="very warm", "red", color))))))
#moderate
mycol1<- rgb(0, 0, 255, max = 255, alpha = 0, names = "blue50")
cols1<- dat.05 %>% 
  dplyr::select(thermal) %>% 
  mutate(color = mycol1, 
         color = ifelse(thermal=="moderate", "orange", color))
                                    
Lcols <- rep(c("Purple", "Orange", "Blue", "Red", "Yellow"))
Lcols <- rep(c("Orange"))

shapes <- dat.05 %>% 
  dplyr::select(graze) %>%
  mutate(shape = 17, shape = ifelse(graze == "grazed", 8, 
                                    ifelse(graze == "ungrazed", 17,shape))) #shapes based on year 
#shapes<-shapes %>% 
#  mutate(time="Pre-Fire", 
#         time= ifelse(year==2004, "Fire",
#               ifelse(year>=2005&year<2014, "Post-Fire",
#              ifelse(year>=2014, "2014 and after",time)))) 

Lshapes <-rep(c(8,17))#shapes for legend
#make the plot
bio.plot <- ordiplot(gb.mds.05, choices=c(1,2), type = "none")   #Set up the plot
points(spscoresall$NMDS1,spscoresall$NMDS2,col=cols1$color,pch=shapes$shape) 
#plot(envvec.nms, col="green")
#text(gb.mds.06, display = "species", cex=0.5, col="grey30") #label species
#legend("topleft",legend=levels(as.factor(cols1$thermal)), col=Lcols, pch=15, cex=0.9,inset=0.1,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)
legend("bottomleft",legend=levels(as.factor(shapes$graze)), col="black", pch=Lshapes, cex=0.9,inset=0.1,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)

############################
#Indicator species for 2005#
############################
trt.dat.05 <- dat.05 %>%
  unite(treatment, c(burn, graze), remove=FALSE, sep = " ")
library(indicspecies)
trt_isa_05 = multipatt(cover.gb.05, trt.dat.05$treatment, control=how(nperm=999))
summary(trt_isa_05)



