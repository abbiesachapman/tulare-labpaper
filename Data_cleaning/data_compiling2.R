# manually set the datpath using the data_pathway file 

# then source the raw data
source("Data_cleaning/dataformatting_2001.R")
source("Data_cleaning/dataformatting_2002.R")
source("Data_cleaning/dataformatting_2003.R")
source("Data_cleaning/dataformatting_2004.R")
source("Data_cleaning/dataformatting_2005.R")
source("Data_cleaning/dataformatting_2006.R")
source("Data_cleaning/dataformatting_2007.R")
source("Data_cleaning/dataformatting_2008.R")
source("Data_cleaning/dataformatting_2009.R")
source("Data_cleaning/dataformatting_2010.R")
source("Data_cleaning/dataformatting_2011.R")
source("Data_cleaning/dataformatting_2012.R")
source("Data_cleaning/dataformatting_2013.R")
source("Data_cleaning/dataformatting_2014.R")
source("Data_cleaning/dataformatting_2015.R")
source("Data_cleaning/dataformatting_2016.R")
source("Data_cleaning/dataformatting_2017.R")
source("Data_cleaning/dataformatting_2018.R")

# rename keys
a <- key2001
b <- key2002
c <- key2003
d <- key2004
e <- key2005
f <- key2006
g <- key2007
h <- key2008
i <- key2009
j <- key2010
k <- key2011
l <- key2012
m <- key2013
n <- key2014
o <- key2015
p <- key2016
q <- key2017
r <- key2018

# create outersect and intersect function  
outersect <- function(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r) {
  big.vec <- c(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r)
  duplicates <- big.vec[duplicated(big.vec)]
  setdiff(big.vec, unique(duplicates))
}

intersect2 <- function(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r) {
  big.vec <- c(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r)
  duplicates <- big.vec[duplicated(big.vec)]
  intersect(big.vec, unique(duplicates))
}

# number of different quadrats across all years
length(outersect(unique(a$quadrat),unique(b$quadrat),unique(c$quadrat),
                 unique(d$quadrat),unique(e$quadrat),unique(f$quadrat),
                 unique(g$quadrat),unique(h$quadrat),unique(i$quadrat),
                 unique(j$quadrat),unique(k$quadrat),unique(l$quadrat),
                 unique(m$quadrat),unique(n$quadrat),unique(o$quadrat),
                 unique(p$quadrat),unique(q$quadrat),unique(r$quadrat)))

# names of different quadrats across all years
outersect(unique(a$quadrat),unique(b$quadrat),unique(c$quadrat),
          unique(d$quadrat),unique(e$quadrat),unique(f$quadrat),
          unique(g$quadrat),unique(h$quadrat),unique(i$quadrat),
          unique(j$quadrat),unique(k$quadrat),unique(l$quadrat),
          unique(m$quadrat),unique(n$quadrat),unique(o$quadrat),
          unique(p$quadrat),unique(q$quadrat),unique(r$quadrat))

# number of same quadrats across all years
length(intersect2(unique(key2001$quadrat),unique(key2002$quadrat),
                 unique(key2003$quadrat),unique(key2004$quadrat),
                 unique(key2005$quadrat),unique(key2006$quadrat),
                 unique(key2007$quadrat),unique(key2008$quadrat),
                 unique(key2009$quadrat),unique(key2010$quadrat),
                 unique(key2011$quadrat),unique(key2012$quadrat),
                 unique(key2013$quadrat),unique(key2014$quadrat),
                 unique(key2015$quadrat),unique(key2016$quadrat),
                 unique(key2017$quadrat),unique(key2018$quadrat)))

# names of same quadrats across all years
intersect2(unique(a$quadrat),unique(b$quadrat),unique(c$quadrat),
          unique(d$quadrat),unique(e$quadrat),unique(f$quadrat),
          unique(g$quadrat),unique(h$quadrat),unique(i$quadrat),
          unique(j$quadrat),unique(k$quadrat),unique(l$quadrat),
          unique(m$quadrat),unique(n$quadrat),unique(o$quadrat),
          unique(p$quadrat),unique(q$quadrat),unique(r$quadrat))
