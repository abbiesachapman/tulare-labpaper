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

# Import csv file, transform to wide data
all.dat <- read_csv(paste(datpath_clean, "/alldatsptrt.csv", sep=""))
all.dat<-as.data.frame(all.dat)
all.dat<-all.dat %>% dplyr::select(-status, -type)

#import env data (litter, cowpies, rocks, etc)
env <- read_csv(paste(datpath_clean, "/envdat.csv", sep=""))
all.dat<- right_join(all.dat,env, by=c("quadratNew", "year"))
all.dat<- all.dat %>% filter(rock<75)

#create wide data, first filter so year is 2005 to 2012
all.data<- all.dat %>% dplyr::select(-X1.x, -spcode, -bare,-rock,-litter,-cowpie,-gopher) %>% filter(year>2004 & year < 2013) %>% spread(spname, cover)
all.data[is.na(all.data)] <- 0 #replace NAs with 0 (species not counted in plots have NAs when wide dataset created)
levels(as.factor(all.dat$quadratNew))
levels(as.factor(all.dat$thermal))
levels(as.factor(all.dat$spname)) #any to remove?
levels(as.factor(all.dat$burn))
levels(as.factor(all.dat$graze))

#############################
#Test effects of burning#####
#########MODERATE############
#######2005-2008#############

#create wide data, first filter so year is 2005 to 2008
mod.data.early<- all.dat %>% dplyr::select(-X1.x, -spcode, -bare, -rock, -litter, -cowpie, -gopher) %>% filter(year>2004 & year < 2009) %>% filter(thermal=="moderate") %>% spread(spname, cover)
mod.data.early[is.na(mod.data.early)] <- 0 #replace NAs with 0 (species not counted in plots have NAs when wide dataset created)

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
cover.mod.early<- mod.data.early.nozero %>% dplyr::select(-c(1:7), -Unknown) #wide data with ID columns removed, only species/cover for NMDS
rownames(cover.mod.early)<-plotnames

#check for empty rows
cover.Biodrop.mod.early<-cover.mod.early[rowSums(cover.mod.early[, (1:68)]) ==0, ] #if no empty rows, next step not needed
#cover.Biodrop.gb<-cover.gb[rowSums(cover.gb[, (1:157)])  >0 ]#remove empty rows

#if needed, relativize by row or column or calculate presence/absence
cover.rowsums.me <- rowSums(cover.mod.early [1:68])
cover.relrow.me <- data.frame(cover.mod.early /cover.rowsums.me)
#cover.colmax<-sapply(cover.mod.early ,max)
#cover.relcolmax <- data.frame(sweep(cover.Bio ,2,cover.colmax,'/'))
#cover.pa <- cover.Bio %>% mutate_each(funs(ifelse(.>0,1,0)), 1:57)
mod.data.early$covercheck<-cover.rowsums.me#remove plots with very low cover?

######################
#NMDS FOR 
######################
#make bray-curtis dissimilarity matrix
mod.bcd.early <- vegdist(cover.relrow.me)
dimcheckMDS(cover.relrow.me, distance = "bray", k = 8, trymax = 20, autotransform = F) #check for optimal dimensions - choose when starts to flatten out
mod.mds.early<-metaMDS(cover.relrow.me, distance="bray", trace = TRUE, noshare=0.02, autotransform=F, trymax=100, k=6) #runs several with different starting configurations
mod.mds.early #solution did not converge after 100 tries, try 1000 more runs
mod.mds.early<-metaMDS(cover.relrow.me, distance="bray", previous.best = mod.mds.early, noshare=0.02, trace = TRUE, autotransform=T, trymax=1000, k=8)
#mod.mds.early<-metaMDS(cover.relrow.me, distance="bray", previous.best = mod.mds.early, noshare=0.02, trace = TRUE, autotransform=T, trymax=5000, k=8)
#mod.mds.early<-metaMDS(cover.relrow.me, distance="bray", previous.best = mod.mds.early, noshare=0.02, trace = TRUE, autotransform=T, trymax=10000, k=8)
summary(mod.mds.early)

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
  xlim(-1,1)+
  ylim(-1,1)+
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
  xlim(-1,1)+
  ylim(-1,1)+
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
  xlim(-1,1)+
  ylim(-1,1)+
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
  xlim(-1,1)+
  ylim(-1,1)+
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
  xlim(-1,1)+
  ylim(-1,1)+
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
#cover.rowsums.05 <- rowSums(cover.2005 [1:156])
#cover.relrow.05 <- data.frame(cover.2005/cover.rowsums.05)
mod.bcd.05 <- vegdist(cover.2005)
permanova05<-adonis(cover.2005~mod.2005$treatment,perm=1000, method="bray")
permanova05
pairwise.perm.manova(mod.bcd.05,mod.2005$treatment, nperm=1000) #all three treatments differ in 2005
mod_isa_05 = multipatt(cover.2005, mod.2005$treatment, control=how(nperm=999))
summary(mod_isa_05)

mod.2006<-subset(mod.data.early, year==2006)
cover.2006<-mod.2006 %>% dplyr::select(-quadratNew,-treatment,-thermal,-burn,-graze, -year, -transect)
#cover.rowsums.06 <- rowSums(cover.2006 [1:156])
#cover.relrow.06 <- data.frame(cover.2006/cover.rowsums.06)
mod.bcd.06 <- vegdist(cover.2006)
permanova06<-adonis(cover.2006~mod.2006$treatment, perm=1000, method="bray")
permanova06
pairwise.perm.manova(mod.bcd.06,mod.2006$treatment, nperm=1000) #all three treatments differ in 2006
mod_isa_06 = multipatt(cover.2006, mod.2006$treatment, control=how(nperm=999))
summary(mod_isa_06)

mod.2007<-subset(mod.data.early, year==2007)
cover.2007<-mod.2007 %>% dplyr::select(-quadratNew,-treatment,-thermal,-burn,-graze, -year, -transect)
#cover.rowsums.07 <- rowSums(cover.2007 [1:156])
#cover.relrow.07 <- data.frame(cover.2007/cover.rowsums.07)
mod.bcd.07 <- vegdist(cover.2007)
permanova07<-adonis(cover.2007~mod.2007$treatment, perm=1000, method="bray")
permanova07
pairwise.perm.manova(mod.bcd.07,mod.2007$treatment, nperm=1000) #ungrazed communities are same, both differ from grazed
mod_isa_07 = multipatt(cover.2007, mod.2007$treatment, control=how(nperm=999))
summary(mod_isa_07)

mod.2008<-subset(mod.data.early, year==2008)
cover.2008<-mod.2008 %>% dplyr::select(-quadratNew,-treatment,-thermal,-burn,-graze, -year, -transect)
#cover.rowsums.08 <- rowSums(cover.2008 [1:156])
#cover.relrow.08 <- data.frame(cover.2008/cover.rowsums.08)
mod.bcd.08 <- vegdist(cover.2008)
permanova08<-adonis(cover.2008~mod.2008$treatment, perm=1000, method="bray")
permanova08
pairwise.perm.manova(mod.bcd.08,mod.2008$treatment, nperm=1000) #ungrazed communities are same, both differ from grazed
mod_isa_08 = multipatt(cover.2008, mod.2008$treatment, control=how(nperm=999))
summary(mod_isa_08)

#####################
#successional vectors on summarized MODERATE data for burn (2005-2008)
####################
mod_yr_burn<-all.dat %>% dplyr::group_by(thermal, burn, graze, year, spname) %>% filter(thermal == "moderate", year>2004 & year<2009) %>% 
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
vec.mds<-metaMDS(cover.relrow, distance="bray", trace = TRUE, autotransform = T, noshare=0.02, trymax=100, k=4) #runs several with different starting configurations
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
mod_vec_isa_early = multipatt(cover.yr, mod_yr_burn$treatment, control=how(nperm=999))
summary(mod_vec_isa_early)

############################
###Make vector plot again in ggplot
###########################

#to plot indicator species (from above) on plot
species.e<-as.data.frame(vec.mds$species)
species.e$name<-row.names(species.e)
spc.e<- species.e %>% filter(name == "Trifolium.depauperatum"| name=="Crassula.connata"| name=="Calandrinia.ciliata"| name == "Rigiopappus.leptoclodus" |name=="Poa.secunda.ssp..secunda"| name=="Festuca.bromoides"|name=="Koeleria.macrantha"| name=="Galium.aparine"| name == "Rigiopappus.leptoclodus" |name=="Festuca.myuros"| name=="Layia.gaillardiodes"| name=="Silene.gallica"|name=="Athysanus.pusilus"|name=="Sisyrinchium.bellum"|name=="Epilobium.sp."| name=="Chlorogalum.pomeridianum"|name=="Sanicula.bipinnatifida"|name=="Lessingia.micradenia.glabratai"|name=="Triteleia.laxa"| name=="Allium.serra"|name=="Plantago.erecta"|name=="Lasthenia.californica"|name=="Aphanes.occidentalis"|name=="Erodium.cicutarium"|name=="Gilia.tricolor"|name=="Lepidium.nitidum"| name=="Hemizonia.congesta"| name=="Castilleja.densiflora"| name=="Microseris.douglasii"|name=="Agoseris.heterophylla"|name=="Brodiaea.spp."|name=="Hordeum.murinum ssp..leporinum"|name=="Muilla.maritima"|name=="Festuca.perennis")
#spc.e<- species.e %>% filter(name == "Trifolium depauperatum"| name=="Crassula connata"| name=="Calandrinia ciliata"| name == "Rigiopappus leptoclodus" |name=="Poa secunda ssp. secunda"| name=="Festuca bromoides"|name=="Koeleria macrantha"| name=="Galium aparine"| name == "Rigiopappus leptoclodus" |name=="Festuca myuros"| name=="Layia gaillardiodes"| name=="Silene gallica"|name=="Athysanus pusilus"|name=="Sisyrinchium bellum"|name=="Epilobium sp."| name=="Chlorogalum pomeridianum"|name=="Sanicula bipinnatifida"|name=="Lessingia micradenia glabratai"|name=="Triteleia laxa"| name=="Allium serra"|name=="Plantago erecta"|name=="Lasthenia californica"|name=="Aphanes occidentalis"|name=="Erodium cicutarium"|name=="Gilia tricolor"|name=="Lepidium nitidum"| name=="Hemizonia congesta"| name=="Castilleja densiflora"| name=="Microseris douglasii"|name=="Agoseris heterophylla"|name=="Brodiaea spp."|name=="Hordeum murinum ssp. leporinum"|name=="Muilla maritima"|name=="Festuca perennis")

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

################
##Test effects of grazing reintroduction
###MODERATE#####
###late years####
#create wide data, first filter so year is 2008 to 2012
mod.data.late<- all.dat %>% dplyr::select(-X1.x, -X1.y, -gopher, -bare, -rock, -litter, -cowpie, -spcode) %>% filter(year>2007 & year < 2013) %>% filter(thermal=="moderate") %>% spread(spname, cover)
mod.data.late[is.na(mod.data.late)] <- 0 #replace NAs with 0 (species not counted in plots have NAs when wide dataset created)

#check count of graze and burn factor
mod.data.late %>% 
  group_by(graze, burn, transect) %>%
  summarise(no_rows = length(graze))

plotnames<-mod.data.late[,1]
cover.mod.late<- mod.data.late %>% dplyr::select(-c(1:6)) #wide data with ID columns removed, only species/cover for NMDS
rownames(cover.mod.late)<-plotnames

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

mod.mds.late<-metaMDS(cover.relrow.ml, trace = TRUE, autotransform=F, noshare=0.01, trymax=100, k=5) #runs several with different starting configurations
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

#create plot in ggplot 
fig1b<-ggplot(subset(spscoresall.mod.l, year==2008), aes(x=NMDS1, y=NMDS2, fill=treatment, shape=as.factor(year)))+
  geom_point(cex=4.5)+
  ggtitle("b) 2008: grazing reintroduced")+
  xlim(-1,1)+
  ylim(-1,1)+
  scale_shape_manual(values=c(21),guide = guide_legend(title = "Year"))+
  scale_fill_manual(values=c("grey0", "grey85", "grey100"), guide = guide_legend(title = "Treatment"), #change legend title
                     labels=c("Burned & Grazed", "Burned & Ungrazed", "Unburned & Ungrazed"))+ #change labels in the legend)+
  theme(plot.title = element_text(color="black", face="bold.italic"))+
  theme(legend.position="none")
fig1b

fig1c<-ggplot(subset(spscoresall.mod.l, year==2009), aes(x=NMDS1, y=NMDS2, fill=treatment, shape=as.factor(year)))+
  geom_point(cex=4.5)+
  ggtitle("c) 2009: 1 year post-grazing")+
  xlim(-1,1)+
  ylim(-1,1)+
  scale_shape_manual(values=c(21),guide = guide_legend(title = "Year"))+
  scale_fill_manual(values=c("grey0", "grey85", "grey100"), guide = guide_legend(title = "Treatment"), #change legend title
                     labels=c("Burned & Grazed", "Burned & Ungrazed", "Unburned & Ungrazed"))+ #change labels in the legend)+
  theme(plot.title = element_text(color="black", face="bold.italic"))+
  theme(legend.position="none")
fig1c
#fig1c+geom_text(subset(spscoresall.mod.e, year==2009), mapping=aes(x=NMDS1, y=NMDS2, label=mod.2009$transect), cex=4)

fig1d<-ggplot(subset(spscoresall.mod.l, year==2010), aes(x=NMDS1, y=NMDS2, fill=treatment, shape=as.factor(year)))+
  geom_point(cex=4.5)+
  ggtitle("d) 2010: 2 years post-grazing")+
  xlim(-1,1)+
  ylim(-1,1)+
  scale_shape_manual(values=c(21),guide = guide_legend(title = "Year"))+
  scale_fill_manual(values=c("grey0", "grey85", "grey100"), guide = guide_legend(title = "Treatment"), #change legend title
                     labels=c("Burned & Grazed", "Burned & Ungrazed", "Unburned & Ungrazed"))+ #change labels in the legend)+
  theme(legend.position="none")+
  theme(plot.title = element_text(color="black", face="bold.italic"))
fig1d

fig1e<-ggplot(subset(spscoresall.mod.l, year==2011), aes(x=NMDS1, y=NMDS2, fill=treatment, shape=as.factor(year)))+
  geom_point(cex=4.5)+
  ggtitle("e) 2011: 3 years post-grazing")+
  xlim(-1,1)+
  ylim(-1,1)+
  scale_shape_manual(values=c(21),guide = guide_legend(title = "Year"))+
  scale_fill_manual(values=c("grey0", "grey85", "grey100"), guide = guide_legend(title = "Treatment"), #change legend title
                     labels=c("Burned & Grazed", "Burned & Ungrazed", "Unburned & Ungrazed"))+ #change labels in the legend)+
  theme(legend.position="none")+
  theme(plot.title = element_text(color="black", face="bold.italic"))
#theme(legend.position="bottom", legend.title=element_text(size=11), legend.text=element_text(size=10), axis.text=element_text(size=8), axis.title=element_text(size=11))+
#theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))
fig1e

fig1f<-ggplot(subset(spscoresall.mod.l, year==2012), aes(x=NMDS1, y=NMDS2, fill=treatment, shape=as.factor(year)))+
  geom_point(cex=4.5)+
  ggtitle("f) 2012: 4 years post-grazing")+
  xlim(-1,1)+
  ylim(-1,1)+
  scale_shape_manual(values=c(21),guide = guide_legend(title = "Year"))+
  scale_fill_manual(values=c("grey0", "grey85", "grey100"), guide = guide_legend(title = "Treatment"), #change legend title
                     labels=c("Burned & Grazed", "Burned & Ungrazed", "Unburned & Ungrazed"))+ #change labels in the legend)+
  theme(legend.position="none")+
  theme(plot.title = element_text(color="black", face="bold.italic"))
#theme(legend.position="bottom", legend.title=element_text(size=11), legend.text=element_text(size=10), axis.text=element_text(size=8), axis.title=element_text(size=11))+
#theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))
fig1f

fig1g<-ggplot(spscoresall.mod.l, aes(x=NMDS1, y=NMDS2, fill=spscoresall.mod.l$treatment, shape=as.factor(mod.data.late$year)))+
  geom_point(cex=4.5,pch = 21)+
  ggtitle("g)")+
  xlim(-1,1)+
  ylim(-1,1)+
  scale_fill_manual(values=c("grey0", "grey85", "grey100"), guide = guide_legend(title = "Treatment:"), #change legend title
                     labels=c("Burned & Grazed", "Burned & Ungrazed", "Unburned & Ungrazed"))+ #change labels in the legend)+
  #theme(legend.position="none")+
  theme(plot.title = element_text(color="black", face="bold.italic"))+
  theme(legend.position="left", legend.title=element_text(size=20), legend.text=element_text(size=17), axis.text=element_text(size=6), axis.title=element_text(size=11))+
  theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))
fig1g

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(fig1g)

##test for differences in treatment by year#######
mod.2008<-subset(mod.data.late, year==2008)
cover.2008<-mod.2008 %>% dplyr::select(-quadratNew,-treatment,-thermal,-burn,-graze, -year, -transect)
#cover.rowsums.05 <- rowSums(cover.2005 [1:156])
#cover.relrow.05 <- data.frame(cover.2005/cover.rowsums.05)
mod.bcd.08 <- vegdist(cover.2008)
permanova08<-adonis(cover.2008~mod.2008$treatment,perm=1000, method="bray")
permanova08
pairwise.perm.manova(mod.bcd.08,mod.2008$treatment, nperm=1000) #ungrazed communities are the same
mod_isa_08 = multipatt(cover.2008, mod.2008$treatment, control=how(nperm=999))
summary(mod_isa_08)

mod.2009<-subset(mod.data.late, year==2009)
cover.2009<-mod.2009 %>% dplyr::select(-quadratNew,-treatment,-thermal,-burn,-graze, -year, -transect)
#cover.rowsums.06 <- rowSums(cover.2006 [1:156])
#cover.relrow.06 <- data.frame(cover.2006/cover.rowsums.06)
mod.bcd.09 <- vegdist(cover.2009)
permanova09<-adonis(cover.2009~mod.2009$treatment, perm=1000, method="bray")
permanova09
pairwise.perm.manova(mod.bcd.09,mod.2009$treatment, nperm=1000) #all communities differ
mod_isa_09 = multipatt(cover.2009, mod.2009$treatment, control=how(nperm=999))
summary(mod_isa_09)

mod.2010<-subset(mod.data.late, year==2010)
cover.2010<-mod.2010 %>% dplyr::select(-quadratNew,-treatment,-thermal,-burn,-graze, -year, -transect)
#cover.rowsums.07 <- rowSums(cover.2007 [1:156])
#cover.relrow.07 <- data.frame(cover.2007/cover.rowsums.07)
mod.bcd.10 <- vegdist(cover.2010)
permanova10<-adonis(cover.2010~mod.2010$treatment, perm=1000, method="bray")
permanova10
pairwise.perm.manova(mod.bcd.10,mod.2010$treatment, nperm=1000) #all communities differ
mod_isa_10 = multipatt(cover.2010, mod.2010$treatment, control=how(nperm=999))
summary(mod_isa_10)

mod.2011<-subset(mod.data.late, year==2011)
cover.2011<-mod.2011 %>% dplyr::select(-quadratNew,-treatment,-thermal,-burn,-graze, -year, -transect)
#cover.rowsums.08 <- rowSums(cover.2008 [1:156])
#cover.relrow.08 <- data.frame(cover.2008/cover.rowsums.08)
mod.bcd.11 <- vegdist(cover.2011)
permanova11<-adonis(cover.2011~mod.2011$treatment, perm=1000, method="bray")
permanova11
pairwise.perm.manova(mod.bcd.11,mod.2011$treatment, nperm=1000) #all communities differ
mod_isa_11 = multipatt(cover.2011, mod.2011$treatment, control=how(nperm=999))
summary(mod_isa_11)

mod.2012<-subset(mod.data.late, year==2012)
cover.2012<-mod.2012 %>% dplyr::select(-quadratNew,-treatment,-thermal,-burn,-graze, -year, -transect)
cover.rowsums.12 <- rowSums(cover.2012 [1:156])
cover.relrow.12 <- data.frame(cover.2012/cover.rowsums.12)
mod.bcd.12 <- vegdist(cover.2012) #replace cover.12 with cover.relrow.12
permanova12<-adonis(cover.2012~mod.2012$treatment, perm=1000, method="bray")
permanova12
pairwise.perm.manova(mod.bcd.12,mod.2012$treatment, nperm=1000) #all communities differ
mod_isa_12 = multipatt(cover.2012, mod.2012$treatment, control=how(nperm=999))
summary(mod_isa_12)

#####################
#successional vectors on summarized MODERATE data for graze (2008-2012)
####################
mod_yr_graze<-all.dat %>% dplyr::group_by(thermal, burn, graze, year, spname) %>% filter(thermal == "moderate", year>2007 & year<2013) %>% dplyr::select(-X1.x, -X1.y, -gopher, -bare, -rock, -litter, -cowpie)%>%
  summarize(mean=mean(cover))%>% arrange(burn)%>%  arrange(graze)%>% arrange(year)%>%
  spread(spname, mean) 
mod_yr_graze[is.na(mod_yr_graze)] <- 0 

cover.yr <- mod_yr_graze %>% ungroup %>% dplyr::select(-thermal,-burn,-graze, -year)
cover.rowsums <- rowSums(cover.yr [1:156])
cover.relrow <- data.frame(cover.yr /cover.rowsums)

#make bray-curtis dissimilarity matrix
vec.bcd <- vegdist(cover.relrow)
dimcheckMDS(cover.relrow)#check for optimal dimensions
#starts to flatten out around 5 or 6

#NMDS 
vec.mds<-metaMDS(cover.relrow, distance="bray", trace = TRUE, autotransform = T, noshare=0.02, trymax=100, k=5) #runs several with different starting configurations
#trace= TRUE will give output for step by step what its doing
#default is 2 dimensions, can put k=4 for 4 dimensions
vec.mds #solution converged after 20 tries
summary(vec.mds)

#quick plot of results
stressplot(vec.mds, vec.bcd) #stressplot to show fit
ordiplot(vec.mds)

#store scores in new dataframe
spscores1.vec<-scores(vec.mds,display="sites",choices=1)
spscores2.vec<-scores(vec.mds,display="sites",choices=2)
year<-mod_yr_graze$year
burn<-mod_yr_graze$burn
spscoresall.vec<-data.frame(burn,year,spscores1.vec,spscores2.vec)

#make plot to show successional vectors
#first, set colors and shapes
mod_yr_graze <- mod_yr_graze %>%
  unite(treatment, c(burn, graze), remove=FALSE, sep = " ")

############################
###Make vector plot again in ggplot
###########################

#to plot indicator species (from above) on plot
species.e<-as.data.frame(vec.mds$species)
species.e$name<-row.names(species.e)
spc.e<- species.e %>% filter(name == "Agoseris.heterophylla"| name=="Microseris.douglasii"| name=="Rigiopappus.leptoclodus"| name == "Deschampsia.danthoniodes" |name=="Stellaria.nitens"| name=="Cryptantha.flaccida"|name=="Sisyrinchium.bellum"| name=="Athysanus.pusilus"| name == "Chlorogalum.pomeridianum" |name=="Epilobium.sp."| name=="Muilla.maritima"| name=="Amsinckia.intermedia"|name=="Sanicula.bipinnatifida"|name=="Allium.serra"|name=="Lasthenia.californica"| name=="Calandrinia.ciliata" | name=="Trifolium.depauperatum" | name=="Gilia.tricolor" | name=="Bromus.madritensis" | name=="Festuca.myuros"| name =="Hordeum.murinum.ssp..leporinum" | name== "Festuca.microstachys"| name=="Calystegia.subacaulis" | name == "Festuca.perennis" | name == "Eriogonum.nudum")

vec1<-ggplot(spscoresall.vec, aes(x=NMDS1, y=NMDS2))+
  geom_path(arrow=arrow(), aes(col=mod_yr_graze$treatment))+
  scale_color_manual(values = c("grey0", "grey60", "grey85")) +
  geom_point(cex=6, pch = 21, aes(fill=mod_yr_graze$treatment))+
  ggtitle("a) Change over time in response to grazing from 2008-2012")+
  scale_fill_manual(values=c("grey0", "grey85", "grey100"), guide = guide_legend(title = "Treatment"), #change legend title
                     labels=c("Burned & Grazed", "Burned & Ungrazed", "Unburned & Ungrazed"))+ #change labels in the legend)+
  theme(legend.position="none")+
  theme(plot.title = element_text(color="black", face="bold.italic"))
#theme(legend.position=c(0.8,0.8), legend.title=element_text(size=14), legend.text=element_text(size=12), axis.text=element_text(size=16), axis.title=element_text(size=16))+
#theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))
vec1

vec1<-vec1+geom_text(data=spc.e, mapping=aes(x=MDS1, y=MDS2, label=name), cex=3)
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
             c(6,6,7,7),
             c(6,6,7,7),
             c(6,6,7,7))
grid.arrange(vec1, fig1b, fig1c, fig1d, fig1e, fig1f, mylegend, layout_matrix = lay) #put panel together
#save as 800w x 1200l