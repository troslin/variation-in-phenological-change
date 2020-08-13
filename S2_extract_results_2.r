# This script extracts the results from fitted HMSC models
# The script assumes that the current folder includes:
# -the R-file "data" under the subfolder "data"
# -the subfolder "models" where the fitted hmsc-models have been placed
# -the subfolder "panels" where the results will be placed

localDir = "."
ModelDir = file.path(localDir, "models")
DataDir = file.path(localDir, "data")
library(Hmsc)
library(corrplot)

# Read and set metadata and
load(file=file.path(DataDir,"data")) #Y,site,year,species,event,group,TrophLevel
load(file = "models/scalings.RData") # metemp,melat,mephper,mechil,meyear,mY,sdY
thin = 1000
samples = 250
nChains = 4

# Loop over the four models and extract results from them
ge.beta = list()
es.beta = list()
for(model in 1:4){
   filename = file.path(ModelDir, paste("model_",as.character(model),
                                        "thin_", as.character(thin),"_samples_", as.character(samples),
                                        "_chains_",as.character(nChains),
                                        ".Rdata",sep = ""))
   load(filename)
   
   # MCMC convergence
   mpost = convertToCodaObject(m)
   ge.beta[[model]] = gelman.diag(mpost$Beta,multivariate=FALSE)$psrf
   es.beta[[model]] = effectiveSize(mpost$Beta)

   # Model fit   
   tiff(paste0("panels//R2_model_",as.character(model),".tiff"))
   predY= computePredictedValues(m,expected = TRUE)
   MF = evaluateModelFit(m,predY=predY)
   tiff(paste0("panels//R2_model_",as.character(model),".tiff"))
   mdoy = m$TrData$mdoy
   gr = m$TrData$group
   mylm=lm(MF$R2~gr*(cos(2*pi*mdoy/365)+sin(2*pi*mdoy/365)))
   plot(m$TrData$mdoy,MF$R2,pch=16,ylim=c(0,1))
   ls = levels(m$TrData$group)
   cols = c("black","orange","green","grey","red")
   for(t in 2:length(ls)){
      sel = m$TrData$group==ls[t]
      points(m$TrData$mdoy[sel],MF$R2[sel],col=cols[t],pch=16)
   }
   for(t in c(1,3,5)){
      sel = m$TrData$group==ls[t]
      tmp = m$TrData$mdoy[sel]
      sdoy = floor(min(tmp)):ceiling(max(tmp))
      preds=predict(mylm,newdata = data.frame(mdoy=sdoy,gr=ls[t]),se.fit = TRUE)
      lines(sdoy,preds$fit,col=cols[t])
      lines(sdoy,preds$fit+preds$se.fit,col=cols[t],lty=2)
      lines(sdoy,preds$fit-preds$se.fit,col=cols[t],lty=2)
   }
   dev.off()
   
   # Variance partitioning    
   if(model==1){groupnames = c("temp","year*temp")}
   if(model==2){groupnames = c("lat","year*lat")}
   if(model==3){groupnames = c("phper","year*phper")}
   if(model==4){groupnames = c("chil","year*chil")}
   VP = computeVariancePartitioning(m,group = c(1,2,1,2),groupnames = groupnames)
   tiff(paste0("panels//VP_model_",as.character(model),".tiff"))
   plotVariancePartitioning(m, VP=VP)
   dev.off()
   
   # Effects of covariates
   
   postBeta = getPostEstimate(m, parName="Beta")
   postBeta$mean = postBeta$mean*sdY
   postBeta$mean[1,] = postBeta$mean[1,]+mY
   tiff(paste0("panels//Beta_model_",as.character(model),".tiff"))
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
   res=NULL
   res$intercept= postBeta$mean[1,]
   res$year= postBeta$mean[2,]
   res$temp= postBeta$mean[3,]
   res$yeartemp = postBeta$mean[4,]
   write.csv(res,file=paste0("panels/results_model_",as.character(model),".csv"))
   
   # Association networks
   OmegaCor = computeAssociations(m)
   supportLevel = 0.95
   for (r in 1:m$nr){
      plotOrder = corrMatOrder(OmegaCor[[r]]$mean,order="AOE")
      plotOrder = 1:m$ns
      toPlot = ((OmegaCor[[r]]$support>supportLevel) + (OmegaCor[[r]]$support<(1-supportLevel))>0)*sign(OmegaCor[[r]]$mean)
      colnames(toPlot)=rep("",m$ns)
      rownames(toPlot)=rep("",m$ns)
      tiff(paste0("panels//Omega_",m$rLNames[[r]],"_model_",as.character(model),".tiff"))
      corrplot(toPlot[plotOrder,plotOrder], method = "color",
               col=colorRampPalette(c("blue","white","red"))(3),
               mar=c(0,0,1,0))
      dev.off()
   }
   
   # Plotting group-specific predictions
   postGamma = getPostEstimate(m, parName="Gamma")
   postGamma$mean = postGamma$mean*sdY
   postGamma$mean[1,1] = postGamma$mean[1,1]+mY
   # To understand how the predictions are generated, it is necessary
   # to know in which order the traits of the phenological events 
   # are coded in the fitted models:
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
   tiff(paste0("panels//Beta_vs_date_model_",as.character(model),".tiff"),
        width = 15,height = 15,units = "cm",res = 300)
   par(mfrow=c(2,2),mar=c(3,3,3,3))
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
   
   tiff(paste0("panels//group_shifts_",as.character(model),".tiff"),
        width = 8.13,height = 9.6,units = "cm",res = 300)
   me = c(metemp,melat,mephper,mechil)[model]
   for(t in c(1,3,5)){
      fgroup = ls[t]
      for(fmdoy in c(100,250)){
         xx = c(min(m$X[,3]),max(m$X[,3]))
         resp =  g[,1] +
            g[,6]*cos(2*pi*fmdoy/365) +
            g[,7]*sin(2*pi*fmdoy/365)
         
         #   if(fgroup=="Abiotic"){} nothing, as reference level
         
         if(fgroup=="Producer"){resp =  resp + g[,3] +
            g[,9]*cos(2*pi*fmdoy/365) +
            g[,13]*sin(2*pi*fmdoy/365)}
         
         if(fgroup=="Secondary_consumer"){resp =  resp + g[,5] +
            g[,11]*cos(2*pi*fmdoy/365) +
            g[,15]*sin(2*pi*fmdoy/365)}
         
         if(t==1 && fmdoy==100){plot(xx+me,resp[2] + resp[4]*xx,type="l",
                                     col=cols[t],
                                     xlab=m$covNames[3],
                                     ylab="phenological shift",
                                     ylim = c(-0.4,0.4))
            abline(0,0,col="gray")
         } else {
            # solid line for fmdoy=100, dashed line for fmoy=250
            lines(xx+me,resp[2] + resp[4]*xx,type="l",col=cols[t],
                  lty=if(fmdoy==250){2}else{1})
         }
      }
   }
   dev.off()
}

# Plot MCMC convergence
tiff(paste0("panels//mixing.tiff"),
     width = 15,height = 10,units = "cm",res = 300)
par(mfrow=c(1,2))
boxplot(ge.beta,names=c("T","L","P","C"),ylim=c(0,6.5),
        ylab = "Potential scale reduction factor")
abline(h=1,col="red")
boxplot(es.beta,names=c("T","L","P","C"),ylim=c(0,2000),
        ylab = "Effective sample size")
abline(h=1000,col="red")
dev.off()
