
library(analysis2Dmito)
library("parallel")
library("dplyr")
library("readr")
library("tidyr")
library("data.table")
library("rstan")
library("stringr")

outputFolder = file.path("..", "Output", "frequentistLinearModel_synData03")
trueParameters_fp = file.path("..", "Data", "trueParameters_synData03.txt")
dataFolderPath = file.path("..", "Data", "synData03")

trueParameters = as.data.frame( fread( trueParameters_fp) ) 

dir.create(file.path(outputFolder), showWarnings = FALSE)

mitochan = "VDAC"
channels = c("NDUFB8", "CYB", "MTCO1")
nChan = length(channels)

ctrlIDs = paste0("C0", 1:4)
ptsIDs = c(paste0("P0", 1:7), "P09", paste0("P", 10:12))

for (chan in channels) {
  for (patID in ptsIDs) {
    root = paste0(chan, "__", patID)
    fp = file.path(dataFolderPath, paste0(root, "__synDATA.txt"))
    data = as.data.frame(fread(fp))
    
    dataWide = as.data.frame( pivot_wider(data, id_cols=c("fibreID", "sampleID", "obsID", "sbjType"), 
                                          names_from="Channel", 
                                          values_from=c("Value")) )
    
    
    ctrlDF = data.frame(mitochan = dataWide[[mitochan]][dataWide$sbjType=="control"], 
                        chan = dataWide[[chan]][dataWide$sbjType=="control"])
    
    mod = lm(chan ~ mitochan, data=ctrlDF)
    
    patDF = data.frame( obsID = dataWide$obsID[dataWide$sbjType=="patient"],
                        mitochan = dataWide[[mitochan]][dataWide$sbjType=="patient"],
                        chan = dataWide[[chan]][dataWide$sbjType=="patient"])
    
    synDF = data.frame(mitochan = seq(min(dataWide[[mitochan]])-2, max(dataWide[[mitochan]])+2, length.out=1000))
    
    synPred = predict.lm(mod, newdata=synDF, interval="prediction")
    
    pred = predict.lm(mod, newdata=patDF, interval="prediction")
    
    synPredDF = as.data.frame(synPred)
    synPredDF$mitochan = synDF$mitochan
    
    class = (patDF$chan < pred[,"lwr"]) | (patDF$chan > pred[,"upr"])
    
    names(class) = patDF$obsID
    
    dataWide$frequentistClass = 0
    
    dataWide[!is.na(match(dataWide$obsID, names(class))), "frequentistClass"] = class
    
    dataLong = pivot_longer(dataWide, cols=c(mitochan, chan), names_to="Channel", values_to="Value")
    
    write.table(dataLong, file=file.path(outputFolder, paste0(root, "__POST.txt")),
                col.names=TRUE, row.names=FALSE, sep="\t")
    
    write.table(synPredDF, file=file.path(outputFolder, paste0(root, "__PRED.txt")),
                col.names=TRUE, row.names=FALSE, sep="\t")
    
    pihat = sum(class) / length(class)
    {
      opp = par(mfrow=c(1,1))
      
      xLim = range( dataWide[[mitochan]] )
      yLim = range( dataWide[[chan]][dataWide$sbjType=="control"] ) + c(-3,3)
      
      xCtrl = data$Value[data$sbjType=="control" & data$Channel==mitochan]
      yCtrl = data$Value[data$sbjType=="control" & data$Channel==chan]
      
      plot(NA, 
           xlim=xLim, xlab=paste0("log(", mitochan, ")"), 
           ylim=yLim, ylab=paste0("log(", chan, ")"), 
           main=bquote(pi*":"~.(round(trueParameters[trueParameters$Channel==chan & trueParameters$sampleID==patID, "probdiff"], 3))*",  "*hat(pi)*":"~.(round(pihat,3))))
      title(main=patID, outer=TRUE, line=-1)
      points(dataWide[[mitochan]][data$sbjType=="control"],
             dataWide[[chan]][data$sbjType=="control"], 
             pch=20, col=alphaBlack(0.05))
      lines(synDF$mitochan, 
            synPred[,"lwr"], lty=2, col=alphaPink(1.0))
      lines(synDF$mitochan, 
            synPred[,"upr"], lty=2, col=alphaPink(1.0))
      
      points(dataWide[[mitochan]][dataWide$sbjType=="patient" & dataWide$frequentistClass==0], 
             dataWide[[chan]][dataWide$sbjType=="patient" & dataWide$frequentistClass==0],
             pch=20, col=alphaBlue(0.1))
      points(dataWide[[mitochan]][dataWide$frequentistClass==1], dataWide[[chan]][dataWide$frequentistClass==1],
             pch=20, col=alphaRed(0.1))
      par(opp)
    }
  }
}

# ### --- --- ---
# ### PLOTS
# ### --- --- ---

# --- --- --- ---
# --- CLASSIFICATION PLOT
# --- --- --- ---

pdf (file=file.path("..", "PDF", basename(outputFolder), "classifications.pdf"), 
     width=9, height=4)
{
  opp = par(mfrow=c(1,3), mar=c(4,4,4,1))
  
  for (patID in ptsIDs) {
    for (chan in channels) {
      root = paste0(chan, "__", patID)
      fp = file.path(outputFolder, paste0(root, "__POST.txt"))
      fp_pred = gsub("__POST.txt", "__PRED.txt", fp)
      
      data = as.data.frame(fread(fp))
      pred = as.data.frame(fread(fp_pred))
      
        xLim = range( data[data$Channel==mitochan, "Value"] )
        yLim = range( data[data$Channel==chan & data$sbjType=="control", "Value"] ) + c(-3,3)
        
        plot(NA, 
             xlim=xLim, xlab=paste0("log(", mitochan, ")"), 
             ylim=yLim, ylab=paste0("log(", chan, ")"), 
             main="")
        
        points(data[data$Channel==mitochan & data$sbjType=="control", "Value"],
               data[data$Channel==chan & data$sbjType=="control", "Value"], 
               pch=20, col=alphaBlack(0.05))
        lines(pred$mitochan, pred$lwr, lty=2, col=alphaGreen(1.0), lwd=2)
        lines(pred$mitochan, pred$upr, lty=2, col=alphaGreen(1.0), lwd=2)
        
        dataWide = as.data.frame( pivot_wider(data, id_cols=c("fibreID", "sampleID", "obsID", "sbjType", "frequentistClass"), 
                                              names_from="Channel", 
                                              values_from=c("Value")) )
        
        nPat = sum(dataWide$sbjType=="patient")
        propDef = sum(dataWide$sbjType=="patient" & dataWide$frequentistClass==1) / nPat

        points(dataWide[dataWide$sbjType=="patient" & dataWide$frequentistClass==0, mitochan], 
               dataWide[dataWide$sbjType=="patient" & dataWide$frequentistClass==0, chan],
               pch=20, col=alphaBlue(0.1))
        points(dataWide[dataWide$sbjType=="patient" & dataWide$frequentistClass==1, mitochan], 
               dataWide[dataWide$sbjType=="patient" & dataWide$frequentistClass==1, chan],
               pch=20, col=alphaRed(0.1))
        
        title(main=bquote(pi*":"~.(round(trueParameters[trueParameters$Channel==chan & trueParameters$sampleID==patID, "probdiff"], 3))*",  "*hat(pi)*":"~.(round(propDef,3))), outer=FALSE)
    }
    title(main=bquote(.(patID)*", N: "*.(nPat)), outer=TRUE, line=-1)
  }
  par(opp)
}
dev.off()









