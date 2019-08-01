## Set your datpath!! (file in Data_cleaning)

library(tidyverse)
library(readr)
library(vegan)
library(MASS)
library(dplyr)
library(RVAideMemoire) #for posthoc tests on permanova

# Import csv file, transform to wide data, call it data
dat <- read_csv(paste(datpath_clean, "/bugrdat.csv", sep=""))
env <- read_csv(paste(datpath_clean, "/NTN-CA66-deposition.csv", sep=""))
env<-env %>% rename(year=yr)
temp <- read_csv(paste(datpath_clean, "/san_jose_clim.csv", sep=""),skip=9)
temp <- temp %>% rename(year=DATE, temp=TAVG) %>% dplyr::select(year, temp)
env<-merge(env,temp)
dat<-as.data.frame(dat)
#dat2<-dat %>% mutate(cover=cover+0.00001) #a work around for data with lots of zeros
data<- dat %>% dplyr::select(-X1, -spcode) %>% spread(spname, cover)
data[is.na(data)] <- 0 #replace NAs with 0 (species not counted in plots have NAs when wide dataset created)
levels(as.factor(dat$quadratNew))
levels(as.factor(dat$thermal))
levels(as.factor(dat$spname)) #any to remove?

#merge env to match full data
all.env<-merge(env,data) %>% dplyr::select(NH4, NO3, totalN, ppt, temp)

data %>% 
  group_by(thermal) %>%
  summarise(no_rows = length(thermal)) #note very uneven

str(data)

data$ID <- seq.int(nrow(data))
plotnames<-data[,1]
cover.Bio<- data %>% dplyr::select(-c(1:4)) #wide data with ID columns removed, only species/cover for NMDS
rownames(cover.Bio)<-plotnames
#check for empty rows

cover.Biodrop<-cover.Bio[rowSums(cover.Bio[, (1:157)]) ==0, ] #no empty rows, next step not needed
#cover.Biodrop<-cover.Bio[rowSums(cover.Bio[, (1:157)])  >0 ]#remove empty rows


#if needed, relativize by row or column or calculate presence/absence
#cover.rowsums <- rowSums(cover.Bio [1:157])
#cover.relrow <- data.frame(cover.Bio /cover.rowsums)
#cover.colmax<-sapply(cover.Bio ,max)
#cover.relcolmax <- data.frame(sweep(cover.Bio ,2,cover.colmax,'/'))
#cover.pa <- cover.Bio %>% mutate_each(funs(ifelse(.>0,1,0)), 1:57)


######################
#2. NMS
#see also section 2.1 of vegan tutorial: 
#http://cc.oulu.fi/~jarioksa/opetus/metodi/vegantutor.pdf
######################

#make bray-curtis dissimilarity matrix
spp.bcd <- vegdist(cover.Bio)

#quick run to check out NMS ordination
spp.mds0 <-isoMDS(spp.bcd) #runs nms only once
spp.mds0  #by default 2 dimensions returned, stress is 6.4, converged
ordiplot(spp.mds0) #ugly

#prefer to run multiple NMS ordinations
#with different starting configurations, and choose the best
#this function does that, and in addition does several other steps too including: 
#standardizing the data (though fuction call below turns this off with autotransform=F)
#calculating distance matrix (default bray-curtis)
#running NMDS with random starts
#rotation of axes to maximize variance of site scores on axis 1
#calculate species scores based on weighted averaging
#help(metaMDS)

spp.mds<-metaMDS(cover.Bio, trace = TRUE, autotransform=T, trymax=100, k=6) #runs several with different starting configurations
#trace= TRUE will give output for step by step what its doing
#default is 2 dimensions, can put k=4 for 4 dimensions
spp.mds #solution did not converge after 100 tries
summary(spp.mds)

#quick plot of results
stressplot(spp.mds, spp.bcd) #stressplot to show fit
ordiplot(spp.mds)
spscores1<-scores(spp.mds,display="sites",choices=1)
spscores2<-scores(spp.mds,display="sites",choices=2)
tplots<-data[,4]
tplot_levels<-levels(tplots)
spscoresall<-data.frame(tplots,spscores1,spscores2)

#overlay environmental variables on full data
envvec.nms<-envfit(spp.mds,all.env, na.rm=TRUE)
envvec.nms
plot(spp.mds)
plot(envvec.nms) #add vectors to previous ordination

#make nicer plot colored based on thermal, shapes on pre/post fire
#help(ordiplot)
#first, set colors and shapes
cols1<- data %>% dplyr::select(thermal) %>% mutate(color = "black", 
       color = ifelse(thermal == "cool", "purple", 
                      ifelse(thermal=="moderate", "orange", 
                             ifelse(thermal=="very cool", "blue",
                                    ifelse(thermal=="very warm", "red", color))))) #colors based on thermal group
Lcols <- rep(c("Black", "Purple", "Orange", "Blue", "Red")) #colors for the legend
shapes <- data %>% dplyr::select(year) %>%
  mutate(shape = 1, shape = ifelse(year == "2004", 8, 
      ifelse(year>="2005", 16, shape))) #shapes based on year 
shapes<-shapes %>% mutate(time="Pre-Fire", 
                        time= ifelse(year==2004, "Fire",
                          ifelse(year>=2005, "Post-Fire", time))) 
Lshapes <-rep(c(8,16,1))#shapes for legend
#make the plot
bio.plot <- ordiplot(spp.mds, choices=c(1,2), type = "none")   #Set up the plot
points(spscoresall$NMDS1,spscoresall$NMDS2,col=cols1$color,pch=shapes$shape) 
plot(envvec.nms)
text(spp.mds, display = "species", cex=0.5, col="grey30") #label species
legend("bottomright",legend=levels(as.factor(cols1$thermal)), col=Lcols, pch=15, cex=0.9,inset=0.1,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)
legend("topright",legend=levels(as.factor(shapes$time)), col="black", pch=Lshapes, cex=0.9,inset=0.1,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)

#separate NMDS plots by thermal designation
#subset data by thermal
cool<-data %>% subset (thermal=="cool") 
mod<-data %>% subset (thermal=="moderate")
vcool<-data %>% subset(thermal=="very cool") 
vwarm<-data %>% subset(thermal=="very warm") 

cover.cool<-data %>% subset (thermal=="cool") %>% dplyr::select(-c(1:4))
cover.mod<-data %>% subset (thermal=="moderate") %>% dplyr::select(-c(1:4))
cover.vcool<-data %>% subset(thermal=="very cool") %>% dplyr::select(-c(1:4))
cover.vwarm<-data %>% subset(thermal=="very warm") %>% dplyr::select(-c(1:4))


##start with COOL
#make bray-curtis dissimilarity matrix
spp.bcd.cool <- vegdist(cover.cool)

spp.mds.cool<-metaMDS(cover.cool, trace = TRUE, autotransform=T, trymax=100, k=2) #runs several with different starting configurations
#trace= TRUE will give output for step by step what its doing
#default is 2 dimensions, can put k=4 for 4 dimensions
spp.mds.cool #solution did not converge after 100 tries
summary(spp.mds.cool)

#quick plot of results
stressplot(spp.mds.cool, spp.bcd.cool) #stressplot to show fit
ordiplot(spp.mds.cool)
spscores1.cool<-scores(spp.mds.cool,display="sites",choices=1)
spscores2.cool<-scores(spp.mds.cool,display="sites",choices=2)
tplots.cool<-cool[,4]
tplot_levels_cool<-levels(tplots.cool)
spscoresall.cool<-data.frame(tplots.cool,spscores1.cool,spscores2.cool)

#overlay environmental variables on plot
cool.env<-merge(cool, env) %>% dplyr::select(NH4, NO3, totalN, ppt)
envvec.nms.cool<-envfit(spp.mds.cool,cool.env, na.rm=TRUE)
envvec.nms.cool
plot(spp.mds.cool)
plot(envvec.nms.cool) #add vectors to previous ordination

#COOL thermal only, shapes on pre/post fire
#help(ordiplot)
bio.plot.cool <- ordiplot(spp.mds.cool, choices=c(1,2), type = "none")   #Set up the plot
cols1<- rep(c("purple"))
shapes.cool <- cool %>% dplyr::select(year) %>%
  mutate(shape = 1, shape = ifelse(year == "2004", 8, 
                                   ifelse(year>="2005", 16, shape))) #shapes based on year 
shapes.cool<-shapes %>% mutate(time="Pre-Fire", 
                               time= ifelse(year==2004, "Fire",
                                            ifelse(year>=2005, "Post-Fire", time))) 
Lshapes <-rep(c(8,16,1))#shapes for legend
points(spscoresall.cool$NMDS1,spscoresall.cool$NMDS2,col=cols1,pch=shapes.cool$shape) 
plot(envvec.nms.cool)
text(spp.mds.cool, display = "species", cex=0.5, col="grey30") #label species
legend("topright",legend=levels(as.factor(shapes$time)), col="black", pch=Lshapes, cex=0.9,inset=0.1,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)

##start with VERY COOL
#make bray-curtis dissimilarity matrix
spp.bcd.vcool <- vegdist(cover.vcool)

spp.mds.vcool<-metaMDS(cover.vcool, trace = TRUE, autotransform=T, trymax=100, k=2) #runs several with different starting configurations
#trace= TRUE will give output for step by step what its doing
#default is 2 dimensions, can put k=4 for 4 dimensions
spp.mds.vcool #solution converged after 20 tries, stress is 22.4%
summary(spp.mds.vcool)

#quick plot of results
stressplot(spp.mds.vcool, spp.bcd.vcool) #stressplot to show fit
ordiplot(spp.mds.vcool)
spscores1.vcool<-scores(spp.mds.vcool,display="sites",choices=1)
spscores2.vcool<-scores(spp.mds.vcool,display="sites",choices=2)
tplots.vcool<-vcool[,4]
tplot_levels_vcool<-levels(tplots.vcool)
spscoresall.vcool<-data.frame(tplots.vcool,spscores1.vcool,spscores2.vcool)

#overlay environmental variables on plot
vcool.env<-merge(vcool, env) %>% dplyr::select(NH4, NO3, totalN, ppt,temp)
envvec.nms.vcool<-envfit(spp.mds.vcool,vcool.env, na.rm=TRUE)
envvec.nms.vcool
plot(spp.mds.vcool)
plot(envvec.nms.vcool) #add vectors to previous ordination

#COOL thermal only, shapes on pre/post fire
#help(ordiplot)
#set colors and shapes
colsvc<- rep(c("blue"))
shapes.vcool <- vcool %>% dplyr::select(year) %>%
  mutate(shape = 1, shape = ifelse(year == "2004", 8, 
                                   ifelse(year>="2005"& year<="2013", 16, 
                                          ifelse(year>="2014", 15, shape)))) #shapes based on year 
shapes.vcool<-shapes.vcool %>% mutate(time="Pre-Fire", 
                               time= ifelse(year==2004, "Fire",
                                            ifelse(year>=2005 & year<="2013", "Post-Fire",
                                                   ifelse(year>=2014, "2014 and later", time)))) 
Lshapesvc <-rep(c(15,8,16,1))#shapes for legend
bio.plot.vcool <- ordiplot(spp.mds.vcool, choices=c(1,2), type = "none")   #Set up the plot
points(spscoresall.vcool$NMDS1,spscoresall.vcool$NMDS2,col=colsvc,pch=shapes.vcool$shape) 
plot(envvec.nms.vcool, col="green")
text(spp.mds.vcool, display = "species", cex=0.5, col="grey30") #label species
legend("topright",legend=levels(as.factor(shapes.vcool$time)), col="black", pch=Lshapesvc, cex=0.9,inset=0.05,bty="n",y.intersp=0.25,x.intersp=0.4,pt.cex=1.1)

