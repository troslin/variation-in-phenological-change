localDir = "."
ModelDir = file.path(localDir, "models")
DataDir = file.path(localDir, "data")
library(Hmsc)

load(file=file.path(DataDir,"data"))
#Y,site,year,species,event,siteinfo,group

sum(!is.na(Y))
mdoy=colMeans(Y,na.rm = TRUE)
ord = order(mdoy)
mdoy = mdoy[ord]
Y = Y[,ord]
species = species[ord]
event = event[ord]
group = group[ord]

ny = dim(Y)[1]
ns = dim(Y)[2]
sites = unique(site)
nsites = length(sites)
plot(siteinfo$Longitude,siteinfo$Latitude)
siteinfo$temp=(siteinfo$temp_1_temp_2_temp_3+
                 siteinfo$temp_4_temp_5_temp_6+
                 siteinfo$temp_7_temp_8_temp_9+
                 siteinfo$temp_10_temp_11_temp_12)/4
studyDesign = data.frame(site=as.factor(site),year=as.factor(year))
temp = rep(NA,ny)
for(i in 1:nsites){
  temp[site==sites[i]]=siteinfo$temp[which(siteinfo$Study.site==sites[i])]
}
metemp = mean(temp)
temp = temp-metemp
meyear = mean(year)
year = year-meyear

XData = data.frame(year,temp)

mY=mean(Y,na.rm = T)
sdY=sd(Y,na.rm = T)
Y=Y-mY
Y=Y/sdY
save(meyear,metemp,mY,sdY,file = "models/scalings.RData")

hist(Y)
hist(year)
hist(mdoy)
TrData = data.frame(mdoy,group)
rL.site = HmscRandomLevel(units = levels(studyDesign$site))
rL.year = HmscRandomLevel(units = levels(studyDesign$year))

XFormula = ~year*temp
TrFormula = ~group * (cos(2*pi*mdoy/365) + sin(2*pi*mdoy/365))
samples = 1000
nChains = 2
for (thin in c(1,10,100,1000)){
  m = Hmsc(Y=Y, XData = XData,  XFormula = XFormula,
           TrData = TrData,TrFormula = TrFormula,
           distr="normal",
           studyDesign=studyDesign, ranLevels={list(site=rL.site,year=rL.year)})
  m = sampleMcmc(m, samples = samples, thin=thin,
                 adaptNf=rep(ceiling(0.4*samples*thin),2), 
                 transient = ceiling(0.5*samples*thin),
                 nChains = nChains)
  filename = file.path(ModelDir, paste("model_thin_", as.character(thin),
                                       "_samples_", as.character(samples),
                                       "_chains_",as.character(nChains),
                                       ".Rdata",sep = ""))
  save(m,file=filename)
}
