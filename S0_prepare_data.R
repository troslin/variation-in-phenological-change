library(geosphere)

setwd("C:/LocalData/OVASKAIN/all stuff/manuscripts/Submitted/Russian phenology shifts")
datadir = "chronicle-of-nature-calendar.v1.0.5"
data=read.csv(file=paste0(datadir,"/Phenology.csv"))

print(dim(data))
data=data[data$quality=="OK",]
print(dim(data))

sum(data$taxon=="Cuculus canorus" & data$eventtype=="1st song")
sel = (data$taxon=="Cuculus canorus" & data$eventtype=="1st occurrence")
sum(sel)
data[sel,]$eventtype = "1st song"

sum(data$taxon=="Taraxacum officinale" & data$eventtype=="onset of blooming")
newtaxon = unique(data[data$taxon == "Taraxacum officinale",]$taxonidentifier)
sel = (data$taxon=="Taraxacum" & data$eventtype=="onset of blooming")
sum(sel)
data[sel,]$taxon = "Taraxacum officinale"
data[sel,]$taxonidentifier = newtaxon

sum(data$taxon=="Rosa acicularis" & data$eventtype=="onset of blooming")
newtaxon = unique(data[data$taxon == "Rosa acicularis",]$taxonidentifier)
sel = (data$taxon=="Rosa" & data$eventtype=="onset of blooming")
sum(sel)
data[sel,]$taxon = "Rosa acicularis"
data[sel,]$taxonidentifier = newtaxon

sum(data$taxon=="Betula pendula" & data$eventtype=="leaf fall end")
newtaxon = unique(data[data$taxon == "Betula pendula",]$taxonidentifier)
sel = (data$taxon=="Betula" & data$eventtype=="leaf fall end")
sum(sel)
data[sel,]$taxon ="Betula pendula"
data[sel,]$taxonidentifier = newtaxon

sum(data$taxon=="Betula pendula" & data$eventtype=="onset of sap bleeding")
newtaxon = unique(data[data$taxon == "Betula pendula",]$taxonidentifier)
sel = (data$taxon=="Betula" & data$eventtype=="onset of sap bleeding")
sum(sel)
data[sel,]$taxon ="Betula pendula"
data[sel,]$taxonidentifier = newtaxon

sum(data$taxon=="Betula pendula" & data$eventtype=="onset of leaf unfolding")
newtaxon = unique(data[data$taxon == "Betula pendula",]$taxonidentifier)
sel = (data$taxon=="Betula" & data$eventtype=="onset of leaf unfolding")
sum(sel)
data[sel,]$taxon ="Betula pendula"
data[sel,]$taxonidentifier = newtaxon

sum(data$taxon=="Betula pendula" & data$eventtype=="onset of autumn colouring")
newtaxon = unique(data[data$taxon == "Betula pendula",]$taxonidentifier)
sel = (data$taxon=="Betula" & data$eventtype=="onset of autumn colouring")
sum(sel)
data[sel,]$taxon ="Betula pendula"
data[sel,]$taxonidentifier = newtaxon

#ONLY 4% OF DATA IS BEFORE 1960 SO DROP THAT
mean(data$year>1959)
data = data[data$year>1959,]

n  = dim(data)[1]
speciesevent = rep(NA,n)
for(i in 1:n){
   speciesevent[i] = paste(data$taxon[i],as.character(data$eventtype[i]))
}

nold = 10^10
n  = dim(data)[1]
#DROP SAMPLING UNITS WITH LESS THAN 10 EVENTS
#DROP EVENTS WITH LESS THAN 100 SAMPLING UNITS
#CONTINUE UNTIL NOTHING DROPPED
while(n<nold){
   print(n)
   nold = n
   siteyear = rep(NA,n)
   speciesevent = rep(NA,n)
   for(i in 1:n){
      siteyear[i] = paste(data$studysite[i],as.character(data$year[i]))
      speciesevent[i] = paste(data$taxon[i],as.character(data$eventtype[i]))
   }
   ta.siteyear = table(siteyear)
   ta.speciesevent = table(speciesevent)
   sel.sy = ta.siteyear>=10
   sel.se = ta.speciesevent>=100
   sel = (siteyear %in% rownames(ta.siteyear[sel.sy])) & (speciesevent %in% rownames(ta.speciesevent[sel.se]))
   data = data[sel,]
   n  = dim(data)[1]
}

# FIXING EVENTS THAT GO OVER DECEMBER-JANUARY
speciesevents = names(ta.speciesevent)
nse = length(speciesevents)
miso = rep(0,nse)
maso = rep(0,nse)
mis = rep(0,nse)
mas = rep(0,nse)
for(i in 1:nse){
   sel = which(speciesevent==speciesevents[i])
   doy = data[sel,]$dayofyear
   mi = min(doy)
   ma = max(doy)
   miso[i] = mi
   maso[i] = ma
   mis[i] = mi
   mas[i] = ma
   if(mi<90 & ma>275){
      print(c(i,mi,ma))
      if(sum(doy>275)<sum(doy<90)){
         sel2 = which(doy>180)
         data[sel[sel2],]$dayofyear = data[sel[sel2],]$dayofyear-365
         data[sel[sel2],]$year = data[sel[sel2],]$year+1
         doy = data[sel,]$dayofyear
      } else {
         sel2 = which(doy<180)
         data[sel[sel2],]$dayofyear = data[sel[sel2],]$dayofyear+365
         data[sel[sel2],]$year = data[sel[sel2],]$year-1
         doy = data[sel,]$dayofyear
      }
      mi = min(doy)
      ma = max(doy)
      mis[i] = mi
      mas[i] = ma      
   }
}
par(mfrow=c(1,2))
plot(miso,maso)
abline(0,1)
plot(mis,mas)
abline(0,1)

species = data$taxon
event = data$eventtype
for(i in 1:n){
   speciesevent[i] = paste(species[i],event[i])
}
siteyears = unique(siteyear)
nsu = length(siteyears)
speciesevents = unique(speciesevent)
ns = length(speciesevents)
species = species[match(speciesevents,speciesevent)]
event = event[match(speciesevents,speciesevent)]


# CONSTRUCT THE Y-MATRIX
Y = matrix(NA,nrow=nsu,ncol=ns)
colnames(Y) = speciesevents
rownames(Y) = siteyears
site = rep(NA,nsu)
year = rep(NA,nsu)
for (i in 1:nsu){
   sel = which(siteyear==siteyears[i])
   site[i] = as.character(unique(data$studysite[sel]))
   year[i] = unique(data$year[sel])
   Y[i,match(speciesevent[sel],speciesevents)] = data[sel,]$dayofyear
}


# REDUCE Y SO THAT ALL EVENTS ARE RECORDED AT LEAST 10 TIMES
# IN AT LEAST 10 SITES,
# AND THAT ALL SAMPLING UNITS HAVE AT LEAST 10 EVENTS
# REMOVE SITE-EVENT PAIRS WITH LESS THAN 10 DATA POINTS
dim(Y)
sum(!is.na(Y))
zold = c(10^10,10^10)
z = dim(Y)
while(sum(z)<sum(zold)){
   print(z)
   zold = z
   #SELECT EVENTS THAT ARE RECORDED AT LEAST 10 TIMES IN AT LEAST 10 SITES
   sites = unique(site)
   
   co = matrix(0,length(sites),ns)
   for(i in 1:length(sites)){
      tmp = !is.na(Y[site == sites[i],])
      if(is.null(dim(tmp))){
         co[i,] = co[i,] + tmp
      } else {
         co[i,] = co[i,] + (colSums(tmp))
      }
      # REMOVE SITE-EVENT PAIRS WITH LESS THAN 10 DATA POINTS
      xx = co[i,]>0 & co[i,]<10
      if(length(xx)>0){
         Y[site == sites[i],xx] = NA
      }
   }
   sel = (colSums(co>=10)>=10)
   Y = Y[,sel]
   speciesevents = speciesevents[sel]
   species = species[sel]
   event = event[sel]
   ns = length(speciesevents)
   
   #SELECT SAMPLING UNITS WITH AT LEAST 10 EVENTS
   sel = rowSums(!is.na(Y))>=10
   Y = Y[sel,]
   site = site[sel]
   year = year[sel]

# DROP SAMPLING UNIT IF AT MOST 25 KM FROM ANOTHER ONE
   siteinfo = read.csv(file="Study sites with temp.csv")
   ma = match(sites,siteinfo$Study.site)
   siteinfo = siteinfo[ma,]
   longlat = matrix(NA,nrow = dim(siteinfo)[1],ncol = 2)
   longlat[,1] = siteinfo$Longitude
   longlat[,2] = siteinfo$Latitude
   di = distm(x=longlat)/1000
   diag(di) = max(di,na.rm = T)
   mi = min(di,na.rm = T)
   if(mi<25){
      ca = which(di==mi,arr.ind = T)
      ca = ca[1,]
      dr1 = sum(!is.na(Y[site==sites[ca[1]],])) # how many data points are lost if missing the site
      dr2 = sum(!is.na(Y[site==sites[ca[2]],]))
      if(dr1>dr2){
         dr = ca[2]
      } else {
         dr = ca[1]
      }
      sel = which(site!=sites[dr])
      Y = Y[sel,]
      site = site[sel]
      year = year[sel]
   }
   z = dim(Y)
}

species = as.character(species)
event = as.character(event)

da = read.csv(file = "data/EventsWithClassifiers.csv")
ind = rep(NA,length(speciesevents))
for (i in 1:length(speciesevents)){
   z = which(da$speciesevents==speciesevents[i])
   if(length(z)>0) {ind[i] = z}
}
mi =  which(is.na(ind))
if(length(mi)>0) {otso}

da = da[ind,]
if(mean(speciesevents==da$speciesevents)==1){
   print("events in same order")
} else{
   print("problem with event order")
}
sel = (!da$TypeOfEvent=="Anthropogenic")
da = da[sel,]
Y = Y[,sel]
species = species[sel]
event = event[sel]

group = da$TrophicLevel

# REMOVING OLD VERSIONS OF TEMPETURATURE DATA
siteinfo$temp_1_temp_2_temp_3=NULL
siteinfo$temp_4_temp_5_temp_6=NULL
siteinfo$temp_7_temp_8_temp_9=NULL
siteinfo$temp_10_temp_11_temp_12=NULL

save(Y,site,year,species,event,siteinfo,group,file="data/data")
