localDir = "."
ModelDir = file.path(localDir, "models")
DataDir = file.path(localDir, "data")
library(Hmsc)
library(corrplot)
# SET DIRECTORIES AND LOAD LIBRARIES (END) ####################################################################

# COMPUTE MODEL FIT (START) ####################################################################
load(file=file.path(DataDir,"data")) #Y,site,year,species,event,group,TrophLevel
load(file = "models/scalings.RData") # metemp,meyear,mY,sdY
thin = 100
samples = 1000
nChains = 2
filename = file.path(ModelDir, paste("model_thin_", as.character(thin),"_samples_", as.character(samples),
                                     "_chains_",as.character(nChains),
                                     ".Rdata",sep = ""))
load(filename)

predY= computePredictedValues(m,expected = TRUE)
MF = evaluateModelFit(m,predY=predY)
mylm=lm(MF$R2~poly(m$TrData$mdoy,degree = 2))
plot(m$TrData$mdoy,MF$R2)
lines(m$TrData$mdoy,fitted.values(mylm))


VP = computeVariancePartitioning(m,group = c(1,2,1,2),groupnames = c("temp","year*temp"))
tiff(paste0("panels//VP.tiff"))
plotVariancePartitioning(m, VP=VP)
dev.off()
# PLOT VARIANCE PARTITIONING (END) #########################################################


# PLOT EFFECTS OF COVARIATES (START) #########################################################
postBeta = getPostEstimate(m, parName="Beta")
postBeta$mean = postBeta$mean*sdY
postBeta$mean[1,] = postBeta$mean[1,]+mY
tiff(paste0("panels//Beta.tiff"))
postBeta2=postBeta
postBeta2$support[3,] = postBeta$support[2,]
postBeta2$mean[3,] = postBeta$mean[2,]
postBeta2$support[2,] = postBeta$support[3,]
postBeta2$mean[2,] = postBeta$mean[3,]
plotBeta(m, post=postBeta2, supportLevel = 0.95,param="Sign",spNamesNumbers=c(FALSE,FALSE))
dev.off()
c(mean(postBeta$support[2,]>0.95),mean(postBeta$support[2,]<0.05))
c(mean(postBeta$support[4,]>0.95),mean(postBeta$support[4,]<0.05))
#EC, EW, LC, LW
c(mean(postBeta$support[4,]>0.95 & postBeta$mean[2,]<0),
  mean(postBeta$support[4,]<0.05 & postBeta$mean[2,]<0),
  mean(postBeta$support[4,]<0.05 & postBeta$mean[2,]>0),
  mean(postBeta$support[4,]>0.95 & postBeta$mean[2,]>0))
res$intercept= postBeta$mean[1,]
res$year= postBeta$mean[2,]
res$temp= postBeta$mean[3,]
res$yeartemp = postBeta$mean[4,]
write.csv(res,file="panels/results.csv")

postGamma = getPostEstimate(m, parName="Gamma")
postGamma$mean = postGamma$mean*sdY
postGamma$mean[1,1] = postGamma$mean[1,1]+mY
tiff(paste0("panels//Gamma.tiff"))
plotGamma(m, post=postGamma, supportLevel = 0.95, param="Sign", trNamesNumbers=c(TRUE,TRUE))
dev.off()
# PLOT EFFECTS OF COVARIATES (END) #########################################################

# PLOT ASSOCIATION NETWORKS (START) #########################################################
OmegaCor = computeAssociations(m)
supportLevel = 0.95
for (r in 1:m$nr){
   plotOrder = corrMatOrder(OmegaCor[[r]]$mean,order="AOE")
   plotOrder = 1:m$ns
   toPlot = ((OmegaCor[[r]]$support>supportLevel) + (OmegaCor[[r]]$support<(1-supportLevel))>0)*sign(OmegaCor[[r]]$mean)
   colnames(toPlot)=rep("",m$ns)
   rownames(toPlot)=rep("",m$ns)
   tiff(paste0("panels//Omega_",m$rLNames[[r]],".tiff"))
   corrplot(toPlot[plotOrder,plotOrder], method = "color",
            col=colorRampPalette(c("blue","white","red"))(3),
            mar=c(0,0,1,0))
   dev.off()
}
# PLOT ASSOCIATION NETWORKS (END) #########################################################

#> m$trNames
#[1] "(Intercept)"                                   
#[2] "groupPrimary_consumer"                         
#[3] "groupProducer"                                 
#[4] "groupSaprotroph"                               
#[5] "groupSecondary_consumer"                       
#[6] "cos(2 * pi * mdoy/365)"                        
#[7] "sin(2 * pi * mdoy/365)"                        
#[8] "groupPrimary_consumer:cos(2 * pi * mdoy/365)"  
#[9] "groupProducer:cos(2 * pi * mdoy/365)"          
#[10] "groupSaprotroph:cos(2 * pi * mdoy/365)"        
#[11] "groupSecondary_consumer:cos(2 * pi * mdoy/365)"
#[12] "groupPrimary_consumer:sin(2 * pi * mdoy/365)"  
#[13] "groupProducer:sin(2 * pi * mdoy/365)"          
#[14] "groupSaprotroph:sin(2 * pi * mdoy/365)"        
#[15] "groupSecondary_consumer:sin(2 * pi * mdoy/365)"

g = postGamma$mean
doy=1:365
tiff(paste0("panels//Beta_vs_date.tiff"),width = 15,height = 15,units = "cm",res = 300)
par(mfrow=c(2,2),mar=c(3,3,3,3))
ls = levels(m$TrData$group)
cols = c("black","orange","green","grey","red")
for(ii in 1:3){
   i=c(2,1,3)[ii]
   sy = rep(1,m$ns)
   sy[(postBeta$support[i+1,]>0.95) |  (postBeta$support[i+1,]<0.05)] = 16
   sel = m$TrData$group==ls[1]
   plot(m$TrData$mdoy[sel],postBeta$mean[i+1,sel],xlab="",ylab = "",
        xlim = c(0,365), ylim = c(min(postBeta$mean[i+1,]),max(postBeta$mean[i+1,])),pch=sy[sel])
   for(t in 2:length(ls)){
      sel = m$TrData$group==ls[t]
      points(m$TrData$mdoy[sel],postBeta$mean[i+1,sel],col=cols[t],pch=sy[sel])
   }
   abline(0,0)
   sel = m$TrData$group==ls[1]
   tmp = m$TrData$mdoy[sel]
   sdoy = floor(min(tmp)):ceiling(max(tmp))
   lines(sdoy,g[i+1,1] +
            g[i+1,6]*cos(sdoy*2*pi/365) +
            g[i+1,7]*sin(sdoy*2*pi/365),
         col=cols[1]) #baseline = abiotic
   sel = m$TrData$group==ls[3]
   tmp = m$TrData$mdoy[sel]
   sdoy = floor(min(tmp)):ceiling(max(tmp))
   lines(sdoy,g[i+1,1] + g[i+1,3] +
            (g[i+1,6]+g[i+1,9])*cos(sdoy*2*pi/365) +
            (g[i+1,7]+g[i+1,13])*sin(sdoy*2*pi/365),
         col=cols[3]) #producer
   sel = m$TrData$group==ls[5]
   tmp = m$TrData$mdoy[sel]
   sdoy = floor(min(tmp)):ceiling(max(tmp))
   lines(sdoy,g[i+1,1] + g[i+1,3] +
            (g[i+1,6]+g[i+1,11])*cos(sdoy*2*pi/365) +
            (g[i+1,7]+g[i+1,15])*sin(sdoy*2*pi/365),
         col=cols[5]) #secondary consumer
}

sy = rep(1,m$ns)
sy[(postBeta$support[4,]>0.95) |  (postBeta$support[4,]<0.05)] = 16
xmax = max(abs(postBeta$mean[2,]))
ymax = max(abs(postBeta$mean[4,]))
sel = m$TrData$group==ls[1]
plot(postBeta$mean[2,sel],postBeta$mean[4,sel],xlab="",ylab = "",
     xlim = c(-xmax,xmax),
     ylim = c(-ymax,ymax),pch=sy[sel])
for(t in 2:length(ls)){
   sel = m$TrData$group==ls[t]
   points(postBeta$mean[2,sel],postBeta$mean[4,sel],col=cols[t],pch=sy[sel])
}
abline(h=0)
abline(v=0)
dev.off()

tiff(paste0("panels//R2.tiff"))
mdoy = m$TrData$mdoy
gr = m$TrData$group
mylm=lm(MF$R2~gr*(cos(2*pi*mdoy/365)+sin(2*pi*mdoy/365)))
plot(m$TrData$mdoy,MF$R2,pch=16,ylim=c(0,1))
for(t in 2:length(ls)){
   sel = m$TrData$group==ls[t]
   points(m$TrData$mdoy[sel],MF$R2[sel],col=cols[t],pch=16)
}
for(t in c(1,3,5)){
   sel = m$TrData$group==ls[t]
   tmp = m$TrData$mdoy[sel]
   sdoy = floor(min(tmp)):ceiling(max(tmp))
   lines(sdoy,predict(mylm,newdata = data.frame(mdoy=sdoy,gr=ls[t])),col=cols[t])
}
dev.off()


for(cc in 1:2){
   if(cc==1){xx = which(postBeta$support[4,]<0.05)}
   if(cc==2){xx = which(postBeta$support[4,]>0.95)}
   for(i in 1:length(xx)){
      x = xx[i]
      y = mY + sdY*m$Y[,x]
      sel = !is.na(y)
      y = y[sel]
      year = meyear+m$XData$year
      year = year[sel]
      site = m$studyDesign$site
      site = site[sel]
      sites = unique(site)
      b = postBeta$mean[,x]
      xd = m$XData[sel,]
      temps = xd$temp
      mitemp = min(xd$temp)
      matemp = max(xd$temp)
      hmscslopes = round(c(b[2]+mitemp*b[4],b[2]+matemp*b[4]),2)
      mylm = lm(y~year+site+year:temps)
      co = mylm$coefficients
      mlslopes = round(c(co[2]+mitemp*co[length(co)],co[2]+matemp*co[length(co)]),2)
      tiff(paste0("panels//TDS//",as.character(cc),"_",as.character(i),".tiff"),width = 15,height = 10,units = "cm",res = 300)
      
      plot(year,y,col = site,
           main=paste0(m$spNames[x],
                       "\n mdoy: ",as.character(round(m$TrData$mdoy[x],1)),
                       "\nHMSC: (",as.character(hmscslopes[1]),",",
                       as.character(hmscslopes[2]),")",
                       "; ML: (",as.character(mlslopes[1]),",",
                       as.character(mlslopes[2]),")")
      )
      for(ssi in 1:length(sites)){
         ss = sites[ssi]
         yy = year[site==ss]
         miy = min(yy)
         may = max(yy)
         tt = mean(xd[site==ss,]$temp)
         pp = predict(mylm,newdata = data.frame(year=c(miy,may),site=c(as.character(ss),as.character(ss)),temps=c(tt,tt)))
         lines(c(miy,may),pp,col=ssi)
      }
      dev.off()
   }
}
