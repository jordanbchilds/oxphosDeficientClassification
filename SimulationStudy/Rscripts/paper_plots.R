# install.packages(c("data.table", "dplyr", "readr", "tidyr", "plyr", "devtools"))
library("data.table")
library("dplyr")
library("readr")
library("tidyr")
library("devtools")
library("parallel")
library("stringr")
library("bayestestR")
library("HDInterval")

# install.packages("../../../analysis2Dmito", type="source", repos=NULL)
library("analysis2Dmito")

parameterFile = file.path("..", "Data", "trueParameters_synData04.txt")
outputFolder = file.path("..", "OutputMaxESS", "stan_sampler_synData04")
freqOutputFolder = file.path("..", "Output", "frequentistLinearModel_synData04")
dataFolder = file.path("..", "Data", "synData04")

dir.create(file.path("..", "PDF"), showWarnings=FALSE)
dir.create(file.path("..", "PDF", basename(outputFolder)), showWarnings=FALSE)

alphaCol = function (col, alpha) {
  if (is.numeric(col)) {
    rgbVals = col2rgb(palette()[i])
  } 
  if (is.character(col)) {
    rgbVals = col2rgb(col)
  }
  return( rgb(red=rgbVals[1] / 255, green=rgbVals[2] / 255, blue=rgbVals[3] / 255, alpha=alpha) )
}

mitochan = "VDAC"
channels = c("NDUFB8", "MTCO1", "CYB")
nChan = length(channels)

ctrlIDs = paste0("C0", 1:4)
ptsIDs = c(paste0("P0", 1:7), "P09", paste0("P", 10:12))

roots = list.files(outputFolder, pattern="__POST.txt")
trueParameters = as.data.frame( fread( parameterFile) ) 

### --- ---
### SINGLE MYOFIBRE MISCLASSIFICATIONS
### --- ---

mitoplotCount = 0
for (patID in ptsIDs) {
  for (chan in channels) {
    data = fread( file.path(dataFolder, paste0(chan, "__", patID, "__postSummaryDATA.txt")))
    data = as.data.frame( data )
    dataPat = data[data$sbjType=="patient", ]
    
    trueDefCount = sum( dataPat$deficientStatus[dataPat$Channel==chan] )
    bayesDefCount = sum( dataPat$meanPosteriorDeficientProb[dataPat$Channel==chan] > 0.5 )
    
    if (trueDefCount != bayesDefCount) {
      mitoplotCount = mitoplotCount + 1
    }
  }
}

tiff (file.path("..", "PDF", basename(outputFolder), "misclassifications.tiff"), 
     width=9, height = 2.5 * floor(mitoplotCount/3) + 1, units="in", res=300)
{
  opp = par(mfrow=c(floor(mitoplotCount/3)+1,3), mar=c(4,4,2,2))
  
  for (patID in ptsIDs) {
    for (chan in channels) {
      
      data = fread( file.path(dataFolder, paste0(chan, "__", patID, "__postSummaryDATA.txt")))
      data = as.data.frame( data )
      dataPat = data[data$sbjType=="patient", ]
      
      trueDefCount = sum( dataPat$deficientStatus[dataPat$Channel==chan] )
      bayesDefCount = sum( dataPat$meanPosteriorDeficientProb[dataPat$Channel==chan] > 0.5 )
      freqDefCount = sum(dataPat$frequentistClassification[dataPat$Channel==chan])
      
      if (trueDefCount != bayesDefCount) {
        postpred = as.data.frame( fread( file.path(outputFolder, paste0(chan, "__", patID, "__POSTPRED.txt"))))
        
        plot(NA, 
             xlab=paste0("log(", mitochan, ")"), xlim=range(data$Value[data$Channel==mitochan]) + c(-1,1),
             ylab=paste0("log(", chan, ")"), ylim=range(data$Value[data$Channel==chan]) + c(-1,1),
             main=patID)
        
        points(data$Value[data$sbjType=="control" & data$Channel==mitochan],
               data$Value[data$sbjType=="control" & data$Channel==chan], 
               pch=20, col=alphaBlack(0.01))
        
        points(data$Value[data$sbjType=="patient" & data$Channel==mitochan],
               data$Value[data$sbjType=="patient" & data$Channel==chan], 
               pch=20, col=classcols(dataPat$meanPosteriorDeficientProb[dataPat$Channel==chan], 0.1))
        
        lines(postpred$mitochan, postpred$lwrNorm, lty=2, lwd=2, col=alphaGreen(1.0))
        lines(postpred$mitochan, postpred$medNorm, lty=1, lwd=2, col=alphaGreen(1.0))
        lines(postpred$mitochan, postpred$uprNorm, lty=2, lwd=2, col=alphaGreen(1.0))
        
        trueDefFibres = data$obsID[data$sbjType=="patient" & data$deficientStatus==1]
        bayesDefFibres = data$obsID[which(data$meanPosteriorDeficientProb>0.5)]
        
        falseDef = bayesDefFibres[ !(bayesDefFibres %in% trueDefFibres)]
        falseHealthy = trueDefFibres[ !(trueDefFibres %in% bayesDefFibres) ]
        
        points(data$Value[data$obsID%in%falseHealthy & data$Channel==mitochan],
               data$Value[data$obsID%in%falseHealthy & data$Channel==chan],
               pch=18, cex=2.0, col="gold" )
        
        points(data$Value[data$obsID%in%falseDef & data$Channel==mitochan],
               data$Value[data$obsID%in%falseDef & data$Channel==chan],
               pch=15, cex=2.0, col="gold" )
      }
    }
  }
  plot(NA, 
       xlim=c(0,1), ylim=c(0,1),
       xlab="", ylab="",
       axes=FALSE)
  legend("topleft", legend=c("False negative", "False positive"), 
         pch=c(18, 15), col="gold", bty="n", bg=NA, 
         cex=2.0)
  
  par(opp)
}
dev.off()

png (file.path("..", "PDF", basename(outputFolder), "misclassifications.png"), 
     width=2250, height= 700 * floor(mitoplotCount/3) + 1, units="px", res=300, pointsize=11)
{
  opp = par(mfrow=c( floor(mitoplotCount/3) + 1,3), mar=c(4,4,2,2))
  for (patID in ptsIDs) {
    for (chan in channels) {
      
      data = fread( file.path(dataFolder, paste0(chan, "__", patID, "__postSummaryDATA.txt")))
      data = as.data.frame( data )
      dataPat = data[data$sbjType=="patient", ]
      
      trueDefCount = sum( dataPat$deficientStatus[dataPat$Channel==chan] )
      bayesDefCount = sum( dataPat$meanPosteriorDeficientProb[dataPat$Channel==chan] > 0.5 )
      
      if (trueDefCount != bayesDefCount) {
        postpred = as.data.frame( fread( file.path(outputFolder, paste0(chan, "__", patID, "__POSTPRED.txt"))))
        
        plot(NA, 
             xlab=paste0("log(", mitochan, ")"), xlim=range(data$Value[data$Channel==mitochan]),
             ylab=paste0("log(", chan, ")"), ylim=range(data$Value[data$Channel==chan & data$sbjType=="control"]) + c(-2,2),
             main=patID)
        
        points(data$Value[data$sbjType=="control" & data$Channel==mitochan],
               data$Value[data$sbjType=="control" & data$Channel==chan], 
               pch=20, col=alphaBlack(0.01))
        
        points(data$Value[data$sbjType=="patient" & data$Channel==mitochan],
               data$Value[data$sbjType=="patient" & data$Channel==chan], 
               pch=20, col=classcols(dataPat$meanPosteriorDeficientProb[dataPat$Channel==chan], 0.1))
        
        lines(postpred$mitochan, postpred$lwrNorm, lty=2, lwd=2, col=alphaGreen(1.0))
        lines(postpred$mitochan, postpred$medNorm, lty=1, lwd=2, col=alphaGreen(1.0))
        lines(postpred$mitochan, postpred$uprNorm, lty=2, lwd=2, col=alphaGreen(1.0))
        
        trueDefFibres = data$obsID[data$sbjType=="patient" & data$deficientStatus==1]
        bayesDefFibres = data$obsID[which(data$meanPosteriorDeficientProb>0.5)]
        
        falseDef = bayesDefFibres[ !(bayesDefFibres %in% trueDefFibres)]
        falseHealthy = trueDefFibres[ !(trueDefFibres %in% bayesDefFibres) ]
        
        points(data$Value[data$obsID%in%falseHealthy & data$Channel==mitochan],
               data$Value[data$obsID%in%falseHealthy & data$Channel==chan],
               pch=18, cex=2.0, col="gold" )
        
        points(data$Value[data$obsID%in%falseDef & data$Channel==mitochan],
               data$Value[data$obsID%in%falseDef & data$Channel==chan],
               pch=15, cex=2.0, col="gold" )
      }
    }
  }
  plot(NA, 
       xlim=c(0,1), ylim=c(0,1),
       xlab="", ylab="",
       axes=FALSE)
  legend("topleft", legend=c("False negative", "False positive"), 
         pch=c(18, 15), col="gold", bty="n", bg=NA, 
         cex=1.25)
  par(opp)
}
dev.off()

### ----------------------------------------------------------------------------
### COMPARE BAYESIAN POSTERIOR WITH TRUE PoD AND FREQUENTIST ESTIMATE
### ----------------------------------------------------------------------------

png(file.path("..", "PDF", basename(outputFolder), "piComparison_synData.png"), 
    width=2250, height=1900, units="px", res=300, pointsize=11)
{
  opp = par(mfrow=c(1,1), mar=c(4,4,1,1))
  
  trueParameters = as.data.frame(fread(parameterFile))
  outsideInterval = list()
  
  piPosterior_diff = list()
  freqEst = list()
  freqEst_diff = list()
  freqProb = list()
  
  for (patID in ptsIDs) {
    for (chan in channels) {
      root = paste0(chan, "__", patID)
      fn = file.path(outputFolder, paste0(root, "__POST.txt"))
      fn_freq = file.path(freqOutputFolder, paste0(root, "__POST.txt"))
      post = fread( fn )
      freqPost = fread( fn_freq )
      
      piEst = sum( freqPost[freqPost$Channel==chan & freqPost$sampleID==patID, "frequentistClass"] ) / sum( freqPost$Channel==chan & freqPost$sampleID==patID )
      
      freqEst[[root]] = piEst
      freqEst_diff[[root]] = piEst - trueParameters[trueParameters$Channel==chan & trueParameters$sampleID==patID, "probdiff"]
      
      piPosterior_diff[[root]] = post$probdiff - trueParameters[trueParameters$Channel==chan & trueParameters$sampleID==patID, "probdiff"]
      hdInterval = hdi(piPosterior_diff[[root]], ci=0.99)
      outsideInterval[[root]] = as.logical( (0.0 < hdInterval["lower"]) | (0.0 > hdInterval["upper"]) )
      
      if (freqEst_diff[[root]] < 0.0) {
        freqProb[[root]] = sum(piPosterior_diff[[root]] < freqEst_diff[[root]]) / length(piPosterior_diff[[root]])
      } else {
        freqProb[[root]] = sum(piPosterior_diff[[root]] > freqEst_diff[[root]]) / length(piPosterior_diff[[root]])
      }
    }
  }
  outsideInterval = unlist(outsideInterval)
  
  tickLocation = (1:length(piPosterior_diff))[c(F,T,F)]
  plot(NA, 
       xlim = c(0, length(piPosterior_diff)), xlab="", 
       ylim = c(-0.2, 1.0), ylab="",
       axes=FALSE)
  for(i in seq(1, length(ptsIDs), 2)) {
    rect( xleft=(i-1)*3 + 0.5, xright=i*3 + 0.5, 
          ybottom=par("usr")[3], ytop=par("usr")[4], 
          border=NA, 
          col=alphaBlack(0.1))
  }
  
  stripchart (piPosterior_diff, add=TRUE,
              pch=20, 
              col=rep(c(alphaCol("lightseagreen", 0.01), alphaCol("deeppink3", 0.01), alphaCol("darkorange", 0.01)), length(ptsIDs)),
              vertical=TRUE, method="jitter", jitter=0.2, 
              ylab="", 
              axes=FALSE)
  
  points (1:length(freqEst_diff), 
          unlist(freqEst_diff), 
          pch=24, 
          col="black",
          bg=rep(c(alphaCol("lightseagreen", 1.0), alphaCol("deeppink3", 1.0), alphaCol("darkorange", 1.0)), length(ptsIDs))
  )
  points (x=c(1:length(piPosterior_diff))[unlist(freqProb)>0.01],
          y=rep(-0.2, sum(unlist(freqProb)>0.01)),
          pch=8, lwd=2)
  points (x=c(1:length(outsideInterval))[outsideInterval],
          y=rep(-0.2, sum(outsideInterval)),
          pch=4, lwd=2, col="red")
  
  
  axis (side=2)
  axis (side=1, 
        at=tickLocation, 
        labels=FALSE)
  text(x = tickLocation, y = par("usr")[3] - 0.05, 
       labels = ptsIDs, srt = 45, adj = 1, xpd = TRUE)
  abline(h=0, col="black", lwd=2, lty=2)
  
  legend(x=0, y=1, 
         bg="transparent",
         box.col=alphaBlack(0.0),
         legend=channels, 
         text.width = 10,
         pch=20, 
         col=c(alphaCol("lightseagreen", 0.7), 
               alphaCol("deeppink3", 0.7), 
               alphaCol("darkorange", 0.7)))
  text(x=0, y=1.0, labels=c("Bayesian distribution"), pos=4)
  
  legend(x=12, y=1.0, 
         bg="transparent", 
         box.col=alphaBlack(0.0),
         legend=channels, 
         text.width = 10,
         pch=24, 
         col="black",
         pt.bg=c(alphaCol("lightseagreen", 0.7), 
                 alphaCol("deeppink3", 0.7), 
                 alphaCol("darkorange", 0.7)) )
  text(x=12, y=1.0, labels=c("Frequentist point estimate"), pos=4)
  par(opp)
}
dev.off()

freqProb[ which(unlist(freqProb)>0.01) ]

if (sum(unlist(outsideInterval))==0) {
  message("There is no evidence of a significant difference between the ground-truth proportion of deficiency and Bayesian posterior.")
} else {
  for (rt in names(outsideInterval)) {
    if (outsideInterval[[rt]]) {
      message("The ground-truth proportion of deficiency for ", rt, " lie outside the 99% posterior HDI.")
    }
  }
}

message("Range in frequentist linear model's difference to ground-truth proportion of deficiency: ", 
        round(min(unlist(freqEst_diff)), 3), ", ", round(max(unlist(freqEst_diff)), 3) )

postMeans = unlist( lapply(piPosterior_diff, mean) )
message("Range of expected posterior difference to gound-truth: ", round(min(postMeans), 3), ", ", round(max(postMeans),3))

### ----------------------------------------------------------------------------
### EXAMPLE DATASET
### ----------------------------------------------------------------------------

alphaCol = function (col, alpha) {
  if (is.numeric(col)) {
    rgbVals = col2rgb(palette()[i])
  } 
  if (is.character(col)) {
    rgbVals = col2rgb(col)
  }
  return( rgb(red=rgbVals[1] / 255, green=rgbVals[2] / 255, blue=rgbVals[3] / 255, alpha=alpha) )
}



png (file.path("..", "PDF", basename(outputFolder), "exampleData.png"),
    width=2250, height=800, pointsize=11, res=300)
{
  op = par(mfrow=c(1,3), mar=c(4,4,1,1))
  patID = "P02"
  
  for (chan in channels) {
    root = paste0(chan, "__", patID)
    
    fn_data = file.path(dataFolder, paste0(root, "__synDATA.txt"))
    dat = fread(fn_data)
    
    thetaTrue = trueParameters[trueParameters$Channel==chan & trueParameters$sampleID==patID, ]
    
    posterior_fn = file.path(outputFolder,  paste0(root, "__POST.txt"))
    post = as.data.frame(fread(posterior_fn))
  
    {
      xLim_syn = range( c(dat$Value[dat$Channel==mitochan]) )  + c(-1, 1)
      yLim_syn = range( c(dat$Value[dat$Channel==chan]) ) + c(-1, 1)
      
      plot(NA,
           xlim = xLim_syn, ylim=yLim_syn, 
           xlab="log(VDAC)", 
           ylab=paste0("log(", chan, ")"), 
           main="")
      points(dat$Value[dat$sbjType=="control" & dat$Channel==mitochan], 
             dat$Value[dat$sbjType=="control" & dat$Channel==chan], 
             pch=20, col=alphaBlack(0.1))
      
      defObsID = dat$obsID[dat$deficientStatus==1]
      
      points(dat$Value[dat$sampleID==patID & dat$Channel==mitochan & !(dat$obsID%in%defObsID)], 
             dat$Value[dat$sampleID==patID & dat$Channel==chan & !(dat$obsID%in%defObsID)], 
             pch=20, col=alphaBlue(0.5))
      
      points(dat$Value[dat$sampleID==patID & dat$Channel==mitochan & dat$obsID%in%defObsID], 
             dat$Value[dat$sampleID==patID & dat$Channel==chan & dat$obsID%in%defObsID], 
             pch=20, col=alphaRed(0.5))
      
      points(dat$Value[dat$sampleID==patID & dat$Channel==mitochan & dat$obsID%in%defObsID], 
             dat$Value[dat$sampleID==patID & dat$Channel==chan & dat$obsID%in%defObsID], 
             pch=20, cex=0.2, col="yellow")
    }
  }
  par(op)
}
dev.off()

### ----------------------------------------------------------------------------
### EXAMPLE FREQUENTIST CLASSIFICATION
### ----------------------------------------------------------------------------


png (file.path("..", "PDF", basename(outputFolder), "exampleFrequentistClass.png"),
     width=2250, height=800, pointsize=11, res=300)
{
  op = par(mfrow=c(1,3), mar=c(4,4,1,1))
  patID = "P02"
  
  for (chan in channels) {
    root = paste0(chan, "__", patID)
    
    fn_data = file.path(dataFolder, paste0(root, "__postSummaryDATA.txt"))
    fn_pred = file.path(freqOutputFolder, paste0(root, "__PRED.txt"))
    
    dat = fread(fn_data)
    pred = fread(fn_pred)
    
    thetaTrue = trueParameters[trueParameters$Channel==chan & trueParameters$sampleID==patID, ]

    {
      xLim_syn = range( c(dat$Value[dat$Channel==mitochan]) )  + c(-1, 1)
      yLim_syn = range( c(dat$Value[dat$Channel==chan]) ) + c(-1, 1)
      
      plot(NA,
           xlim = xLim_syn, ylim=yLim_syn, 
           xlab=paste0("log(", mitochan, ")"), 
           ylab=paste0("log(", chan, ")"), 
           main="")
      points(dat$Value[dat$sbjType=="control" & dat$Channel==mitochan], 
             dat$Value[dat$sbjType=="control" & dat$Channel==chan], 
             pch=20, col=alphaBlack(0.1))
      
      
      true_defObsID = dat$obsID[which(dat$deficientStatus==1)]
      defObsID = dat$obsID[which(dat$frequentistClassification==1)]
      
      points(dat$Value[dat$sampleID==patID & dat$Channel==mitochan & !(dat$obsID%in%defObsID)], 
             dat$Value[dat$sampleID==patID & dat$Channel==chan & !(dat$obsID%in%defObsID)], 
             pch=20, col=alphaBlue(0.5))
      points(dat$Value[dat$sampleID==patID & dat$Channel==mitochan & dat$obsID%in%defObsID], 
             dat$Value[dat$sampleID==patID & dat$Channel==chan & dat$obsID%in%defObsID], 
             pch=20, col=alphaRed(0.5))
      points(dat$Value[dat$sampleID==patID & dat$Channel==mitochan & dat$obsID%in%true_defObsID], 
             dat$Value[dat$sampleID==patID & dat$Channel==chan & dat$obsID%in%true_defObsID], 
             pch=20, cex=0.33, col="gold")
      
      lines(pred$mitochan, pred$lwr, lty=2, lwd=2, col=alphaGreen(1.0))
      lines(pred$mitochan, pred$fit, lty=2, lwd=2, col=alphaGreen(1.0))
      lines(pred$mitochan, pred$upr, lty=2, lwd=2, col=alphaGreen(1.0))
    }
  }
  par(op)
}
dev.off()

### ----------------------------------------------------------------------------
### EXAMPLE POSTERIOR
### ----------------------------------------------------------------------------

png (file.path("..", "PDF", basename(outputFolder), "examplePosterior.png"), 
    width=2250, height=2000, pointsize=11, res=300)
{
  patID = "P02"
  chan = "NDUFB8"
  root = paste0(chan, "__", patID)
  
    trueParameters = as.data.frame( fread( parameterFile ) )
    thetaTrue = trueParameters[ trueParameters$sampleID==patID & trueParameters$Channel==chan, ]
    
    dataFp = file.path(dataFolder, paste0(root, "__synDATA.txt"))
    data = as.data.frame( fread( dataFp ) )
    
    fn_post = file.path(outputFolder, paste0(root, "__POST.txt"))
    post = as.data.frame( fread(fn_post) )
    
    fn_prior = file.path(outputFolder, paste0(root, "__PRIOR.txt"))
    prior = as.data.frame( fread(fn_prior) )
    
    fn_postpred = file.path(outputFolder, paste0(root, "__POSTPRED.txt"))
    postpred = as.data.frame( fread(fn_postpred) )
    
    fn_class = file.path(outputFolder, paste0(root, "__CLASSIF.txt"))
    class = apply( as.matrix(fread(fn_class) ), 2, mean)
    
    dataMats = getData_mats(data, channels=c(mitochan, chan), pts=patID, ctrlID=ctrlIDs)
    
    op = par(mfrow = c(3, 3), mar = c(4, 4, 1, 1), 
             cex.main = 1.0, cex.lab = 1.0, cex.axis = 1.0)
    {
      # ----------
      # --- slope
      # ----------
      
      {
        slopeNames = paste0("m[", 1:5, "]")
        slopeDensities = list()
        yMax = 0
        xMax = -1e99
        xMin = 1e99
        for (i in 1:5) {
          slopeDensities[[i]] = density(post[[slopeNames[i]]])
          yMax = ifelse(yMax < max(slopeDensities[[i]]$y), max(slopeDensities[[i]]$y), yMax)
          xMax = ifelse(xMax < max(slopeDensities[[i]]$x), max(slopeDensities[[i]]$x), xMax)
          xMin = ifelse(xMin > min(slopeDensities[[i]]$x), min(slopeDensities[[i]]$x), xMin)
        }
        
        xLim = range( unlist(c(xMin, xMax, thetaTrue[slopeNames]) ) )
        
        plot (NA, 
              xlab="Slope", xlim=xLim, 
              ylab="Density", ylim=c(0, yMax))
        lines( density(prior[["m_pred"]]), 
               lty=2, lwd=2, col=alphaPink(0.7))
        lines( density(post[["m_pred"]]), 
               lty=2, lwd=2, col=alphaGreen(0.7))
        for (i in 1:4) {
          lines(slopeDensities[[i]], col=alphaCol(i+1, 0.5))
          abline(v = thetaTrue[paste0("m[", i, "]")], col=alphaCol(i+1, 1.0), lwd=2, lty=3)
        }
        lines(slopeDensities[[5]], col=alphaGreen(1.0), lwd=2)
        abline(v = thetaTrue[paste0("m[5]")], col=alphaGreen(1.0), lwd=2, lty=3)
      }
      
      # --------------
      # --- intercept
      # --------------
      
      {
        interceptNames = paste0("c[", 1:5, "]")
        interceptDensities = list()
        yMax = 0
        xMax = -1e99
        xMin = 1e99
        for (i in 1:5) {
          interceptDensities[[i]] = density(post[[interceptNames[i]]])
          yMax = ifelse(yMax < max(interceptDensities[[i]]$y), max(interceptDensities[[i]]$y), yMax)
          xMax = ifelse(xMax < max(interceptDensities[[i]]$x), max(interceptDensities[[i]]$x), xMax)
          xMin = ifelse(xMin > min(interceptDensities[[i]]$x), min(interceptDensities[[i]]$x), xMin)
        }
        
        plot (NA, 
              xlab="Intercept", xlim=c(xMin, xMax), 
              ylab="Density", ylim=c(0, yMax))
        lines( density(prior[["c_pred"]]), 
               lty=2, lwd=2, col=alphaPink(0.7))
        lines( density(post[["c_pred"]]), 
               lty=2, lwd=2, col=alphaGreen(0.7))
        for (i in 1:4) {
          lines(interceptDensities[[i]], col=alphaCol(i+1, 0.5))
          abline(v = thetaTrue[interceptNames[i]], col=alphaCol(i+1, 1.0), lwd=2, lty=3)
        }
        lines(interceptDensities[[5]], col=alphaGreen(1.0), lwd=2)
        abline(v = thetaTrue[interceptNames[5]], col=alphaGreen(1.0), lwd=2, lty=3)}
      
      # -------------------------
      # --- proportion deficient
      # -------------------------
      
      {
        piDens = density(post$probdiff)
        
        curve( dunif(x, 0, 0.5), 
               col=alphaPink(0.5), lwd=2,
               xlab=bquote(pi),  xlim=c(-0.1, 0.6),
               ylab="Density", ylim=c(0, max(piDens$y)),
               main="")
        lines(piDens, col=alphaGreen(1.0), lwd=2)
        abline(v=thetaTrue["probdiff"], lty=3, lwd=2, col=alphaGreen(1.0))
      }
      
      # -----------------------
      # --- rest of parameters
      # -----------------------
      
      {
        paramNames = c("mu_m", "tau_m", "tau_norm", "mu_c", "tau_c")
        paramLabels = list(mu_m = bquote(mu[m]), 
                           tau_m = bquote(tau[m]), 
                           tau_norm = bquote(tau), 
                           mu_c = bquote(mu[c]), 
                           tau_c = bquote(tau[c]))
        
        for (param in paramNames) {
          postDens = density(post[, param])
          priorDens = density(prior[,param])
          xLim = range(c(postDens$x, priorDens$x))
          
          plot( density(prior[,param]), 
                col=alphaPink(0.5), lwd=2,
                xlab=paramLabels[[param]], xlim=xLim, 
                ylab="Density", ylim=c(0, max(postDens$y, priorDens$y)),
                main="")
          lines(postDens, col=alphaGreen(1.0), lwd=2)
          abline(v=thetaTrue[param], lty=3, lwd=2, col=alphaGreen(1.0))
        }
      }
      
      # ------------------------
      # --- classification plot
      # ------------------------
      
      {
        xLim = range( data$Value[data$Channel==mitochan]) + c(-1,1)
        yLim = range( data$Value[data$Channel==chan] ) + c(-1,1)
        
        plot(NA, 
             xlim=xLim, xlab=paste0("log(", mitochan, ")"),
             ylim=yLim, ylab=paste0("log(", chan, ")"))
        
        points(data[data$sbjType=="control" & data$Channel==mitochan, "Value"], 
               data[data$sbjType=="control" & data$Channel==chan, "Value"], 
               pch=20, col=alphaBlack(0.02))
        
        points(data[data$sampleID==patID & data$Channel==mitochan, "Value"], 
               data[data$sampleID==patID & data$Channel==chan, "Value"], 
               pch=20, col=classcols(class, alphaLevel=0.2))
        
        defObsID = data$obsID[data$deficientStatus==1]
        
        points(data$Value[data$sampleID==patID & data$Channel==mitochan & data$obsID%in%defObsID], 
               data$Value[data$sampleID==patID & data$Channel==chan & data$obsID%in%defObsID], 
               pch=20, col="yellow", cex=0.2)
        
        lines(postpred$mitochan, 
              postpred$medNorm, 
              lty=1, 
              col=alphaGreen(1.0), lwd=2)
        
        lines(postpred$mitochan, 
              postpred$lwrNorm, 
              lty=2, 
              col=alphaGreen(1.0), lwd=2)
        
        lines(postpred$mitochan, 
              postpred$uprNorm, 
              lty=2, 
              col=alphaGreen(1.0), lwd=2)
        
      }
    }
    par(op)
}
dev.off()

### ----------------------------------------------------------------------------
### TRAINING vs VALIDATION 
### ----------------------------------------------------------------------------


BayesProbabilityCalculator = function (dataPoint, post, gamma=0.0001) {
  
  expected = dataPoint[1] * post[,"m[5]"] + post[,"c[5]"]
  densHealthy = dnorm(dataPoint[2], mean=expected, sd=1/sqrt(post[,"tau_norm"]))
  densDeficient = dnorm(dataPoint[2], mean=expected, sd=1/sqrt(gamma))
  
  top = post[,"probdiff"] * densDeficient
  bottom = (1 - post[,"probdiff"]) * densHealthy + top
  return ( top / bottom )
}


dataFolder_cv = file.path("..", "Data", "synCVData03")

outputFolder_cv = file.path("..", "OutputMaxESS", "stan_sampler_synCVData03")
freqOutputFolder_cv = file.path("..", "Output", "frequentistLinearModel_synCVData03")

classifDiff = list()
classifSame = list()
for (chan in channels) {
  for (patID in ptsIDs) {
    root = paste0(chan, "__", patID)
    
    classifDiff[[root]] = c()
    classifSame[[root]] = c()
    
    fnData_og = file.path(dataFolder, paste0(root, "__postSummaryDATA.txt"))
    fnData_cv = file.path(dataFolder_cv, paste0(root, "__postSummaryDATA.txt"))
    fnPost_og = file.path(outputFolder, paste0(root, "__POST.txt"))
    fnPost_cv = file.path(outputFolder_cv, paste0(root, "__POST.txt"))
    

    data_og = fread( fnData_og )
    data_cv = fread( fnData_cv )
    post_og = as.data.frame( fread(fnPost_og) )
    post_cv = as.data.frame( fread(fnPost_cv) )
    
    obsIDs = unique(data_og$obsID[data_og$sbjType=="patient"])
    
    for (id in obsIDs) {
      dataFibre = data_og[data_og$obsID==id, ]
      dataPoint = c( dataFibre$Value[dataFibre$Channel==mitochan], 
                     dataFibre$Value[dataFibre$Channel==chan] )
      
      postProbDef_og = BayesProbabilityCalculator(dataPoint, post_og)
      postProbDef_cv = BayesProbabilityCalculator(dataPoint, post_cv)
      
      diff = postProbDef_og - postProbDef_cv
      
      interval = hdi(diff, 0.99)
      if (0<interval["lower"] | 0>interval["upper"]) {
        classifDiff[[root]] = c(classifDiff[[root]], id)
      } else {
        classifSame[[root]] = c(classifSame[[root]], id)
      }
    }
  }
}


### ----------------------------------------------------------------------------
### HOW MUCH DATA
### ----------------------------------------------------------------------------

totalPatientFibres = 0
for (chan in channels) {
  for (patID in ptsIDs) {
    root = paste0(chan, "__", patID)
    fnData = file.path(dataFolder, paste0(root, "__postSummaryDATA.txt"))
    data = fread( fnData )
    totalPatientFibres = totalPatientFibres + sum(data$sbjType=="patient" & data$Channel==chan)
  }
}

