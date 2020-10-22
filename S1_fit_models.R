# This script fits HMSC models to the phenology data
# The script assumes that the current folder includes:
# -the R-file "data" under the subfolder "data"
# -the subfolder "models" where the fitted hmsc-models will be placed

localDir = "."
ModelDir = file.path(localDir, "models")
DataDir = file.path(localDir, "data")
library(Hmsc)

load(file=file.path(DataDir,"data"))
# the datafile includes the following objects:
# Y: the matrix of phenology data, with events as columns, and sampling units (site-year combinations) as rows
# site: the site of each sampling unit
# year: the year of each sampling unit
# species: the species of each column
# event: the event type of each column
# group: the grouping of the species
# siteinfo: site-specific predictors

# Ordering the columns according to the mean day of the event
mdoy=colMeans(Y,na.rm = TRUE)
ord = order(mdoy)
mdoy = mdoy[ord]
Y = Y[,ord]
species = species[ord]
event = event[ord]
group = group[ord]

# Computing the photoperiod from latitude and longitude
ny = dim(Y)[1]
ns = dim(Y)[2]
sites = unique(site)
nsites = length(sites)
nsites2 = dim(siteinfo)[1]
photoperiod = rep(NA,nsites2)
for (i in 1:nsites2){
  tmp = daylength(lat = siteinfo$Latitude[i],doy = 1:365)
  photoperiod[i] = max(tmp)-min(tmp)
}
siteinfo$photoperiod = photoperiod

# Preparing the variables needed to define the HMSC model

studyDesign = data.frame(site=as.factor(site),year=as.factor(year),sample=as.factor(1:ny))
temp = rep(NA,ny)
lat = rep(NA,ny)
phper = rep(NA,ny)
chil = rep(NA,ny)
for(i in 1:nsites){
  temp[site==sites[i]]=siteinfo$tempAvg[which(siteinfo$Study.site==sites[i])]
  lat[site==sites[i]]=siteinfo$Latitude[which(siteinfo$Study.site==sites[i])]
  phper[site==sites[i]]=siteinfo$photoperiod[which(siteinfo$Study.site==sites[i])]
  chil[site==sites[i]]=siteinfo$tempNegAvg[which(siteinfo$Study.site==sites[i])]
}

metemp = mean(temp)
melat = mean(lat)
mephper = mean(phper)
mechil = mean(chil)
meyear = mean(year)
temp = temp-metemp
lat = lat-melat
phper = phper - mephper
chil = chil-mechil  
year = year-meyear

XData = data.frame(year,temp,lat,phper,chil)

mY=mean(Y,na.rm = T)
sdY=sd(Y,na.rm = T)
Y=Y-mY
Y=Y/sdY
save(meyear,metemp,melat,mephper,mechil,mY,sdY,file = "models/scalings.RData")

TrData = data.frame(mdoy,group)
rL.sample = HmscRandomLevel(units = levels(studyDesign$sample))
rL.site = HmscRandomLevel(units = levels(studyDesign$site))
rL.year = HmscRandomLevel(units = levels(studyDesign$year))

XFormulas = c(~year*temp,~year*lat,~year*phper,~year*chil)
TrFormula = ~group * (cos(2*pi*mdoy/365) + sin(2*pi*mdoy/365))

# Fitting the four alternative HMSC models with increasing thinning
samples = 25
nChains = 4
for (thin in c(1,10,100,1000)){
  for (model in 1:4){
    m = Hmsc(Y=Y, XData = XData,  XFormula = XFormulas[[model]],
             TrData = TrData,TrFormula = TrFormula,
             distr="normal",
             studyDesign=studyDesign, ranLevels={list(site=rL.site,year=rL.year,sample=rL.sample)})
    m = sampleMcmc(m, samples = samples, thin=thin,
                   adaptNf=rep(ceiling(0.4*samples*thin),m$nr), 
                   transient = ceiling(0.5*samples*thin),
                   nChains = nChains)
    filename = file.path(ModelDir, paste("model_", as.character(model),
                                         "thin_", as.character(thin),
                                         "_samples_", as.character(samples),
                                         "_chains_",as.character(nChains),
                                         ".Rdata",sep = ""))
    save(m,file=filename)
  }
}
