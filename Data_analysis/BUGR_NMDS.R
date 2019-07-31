## Set your datpath!! (file in Data_cleaning)

library(tidyverse)
library(readr)
library(vegan)
library(MASS)
library(dplyr)
library(RVAideMemoire) #for posthoc tests on permanova

# Import csv file, transform to wide data, call it data
dat <- read_csv(paste(datpath_clean, "/bugrdat.csv", sep=""))
dat<-as.data.frame(dat)
#dat2<-dat %>% mutate(cover=cover+0.00001) #a work around for data with lots of zeros
data<- dat2 %>% dplyr::select(-X1, -spcode) %>% spread(spname, cover)
data[is.na(data)] <- 0
levels(as.factor(dat$quadratNew))

levels(as.factor(dat$thermal))

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


#if needed, relativize by row or colum or calculate presence/absence
cover.rowsums <- rowSums(cover.Bio [1:157])
cover.relrow <- data.frame(cover.Bio /cover.rowsums)
#cover.colmax<-sapply(cover.Bio ,max)
#cover.relcolmax <- data.frame(sweep(cover.Bio ,2,cover.colmax,'/'))
#cover.pa <- cover.Bio %>% mutate_each(funs(ifelse(.>0,1,0)), 1:57)


######################
#2. NMS
#see also section 2.1 of vegan tutorial: 
#http://cc.oulu.fi/~jarioksa/opetus/metodi/vegantutor.pdf
######################

#make bray-curtis dissimilarity matrix
spp.bcd <- vegdist(cover.relrow)

#run NMS ordination
spp.mds0 <-isoMDS(spp.bcd) #runs nms only once
spp.mds0  #by default 2 dimensions returned, stress is 7.5, converged
ordiplot(spp.mds0)

#prefer to run multiple NMS ordinations
#with different starting configurations, and choose the best
#this function does that, and in addition does several other steps too
#including: 
#standardizing the data (though fuction call below turns this off with autotransform=F)
#calculating distance matrix (default bray-curtis)
#running NMDS with random starts
#rotation of axes to maximize variance of site scores on axis 1
#calculate species scores based on weighted averaging

help(metaMDS)
spp.mds<-metaMDS(cover.relrow, trace = TRUE, autotransform=F, trymax=100, k=6) #runs several with different starting configurations
#trace= TRUE will give output for step by step what its doing
#default is 2 dimensions, can put k=4 for 4 dimensions
spp.mds #solution converged after 49 tries, stress = 3.0
summary(spp.mds)

#plot results
stressplot(spp.mds, spp.bcd) #stressplot to show fit
ordiplot(spp.mds)
spscores1<-scores(spp.mds,display="sites",choices=1)
spscores2<-scores(spp.mds,display="sites",choices=2)
tplots<-data[,4]
tplot_levels<-levels(tplots)
spscoresall<-data.frame(tplots,spscores1,spscores2)

#plots colored based on treatment #help(ordiplot)
bio.plot <- ordiplot(spp.mds, choices=c(1,2), type = "none")   #Set up the plot
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
                          ifelse(year>=2006, "Post-Fire", time)))
Lshapes <-rep(c(8,16,1))#shapes for legend
points(spscoresall$NMDS1,spscoresall$NMDS2,col=cols1$color,pch=shapes$shape) 
#text(spp.mds, display = "species", cex=0.5, col="grey30") #label species
legend("bottomright",legend=levels(as.factor(cols1$thermal)), col=Lcols, pch=15, cex=0.9,inset=0.1,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)
legend("topright",legend=levels(as.factor(shapes$time)), col="black", pch=Lshapes, cex=0.9,inset=0.1,bty="n",y.intersp=0.5,x.intersp=0.8,pt.cex=1.1)

