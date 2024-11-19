library("tidyr")
# ------------------------------------------------------------------------------
# FREQUENTIST LINEAR REGRESSION 
# ------------------------------------------------------------------------------

dir.create(file.path("..", "Output"), showWarnings = FALSE)
dir.create(file.path("..", "Output", "frequentistAnalysis"), showWarnings = FALSE)

mitochan = "VDAC"
channels = c("NDUFB8", "CYB", "MTCO1")
nChan = length(channels)

raw_data = read.csv(file.path("..", "Data", "Data_prepped.csv"), header=TRUE)

raw_data = raw_data[,c("ID", "patient_id", mitochan, channels)]
colnames(raw_data) = c("fibreID", "sampleID", mitochan, channels)

data_lng = pivot_longer(raw_data, cols=c(mitochan, channels), names_to="Channel")
data_lng = as.data.frame(data_lng)

colnames(data_lng) = c("fibreID", "sampleID", "Channel", "Value")

data = data_lng
data$Value = log(data$Value)

sbjIDs = unique(data$sampleID)
ctrlIDs = grep("C", sbjIDs, value=TRUE)
ptsIDs = sbjIDs[!(sbjIDs %in% ctrlIDs)]
data$sbjType = "patient"
data$sbjType[data$sampleID%in%ctrlIDs] = "control"

data$obsID = paste(data$sampleID, data$fibreID, sep="_")

data$Classification = NA

dataCtrl = data[data$sbjType=="control", ]

dataWide = pivot_wider(data, id_cols=c("obsID", "fibreID", "sampleID", "sbjType"), 
                       names_from="Channel", values_from="Value")
dataWide$NDUFB8_class = NA
dataWide$CYB_class = NA
dataWide$MTCO1_class = NA

dataWide = as.data.frame(dataWide)

for (chan in channels) {
  obsID = dataWide[dataWide$sbjType=="control", "obsID"]
  xCtrl = dataWide[dataWide$sbjType=="control", mitochan]
  yCtrl = dataWide[dataWide$sbjType=="control", chan]
  
  ctrl_df = data.frame(obsID=obsID, mitochan=xCtrl, chan=yCtrl)
  mod = lm (chan~mitochan, data=ctrl_df)
  
  for (patID in ptsIDs) {
    xPat = dataWide[dataWide$sampleID==patID, c("obsID", mitochan, chan)]
    colnames(xPat) = c("obsID", "mitochan", "chan")
    patPred = predict.lm(mod, newdata=xPat, interval="prediction")
    
    dataWide[!is.na(match(dataWide$obsID, xPat$obsID)), paste0(chan, "_class")] = (xPat$chan > patPred[,"upr"]) | (xPat$chan < patPred[,"lwr"])
  }
  
  xSyn = seq(from=min(data$Value) - 1.0, to=max(data$Value) + 1.0, length.out=1000)
  synDF = data.frame(mitochan=xSyn)
  predInterval = as.data.frame( predict.lm(mod, newdata=synDF, interval="prediction") )
  predInterval$mitochan = xSyn
  write.table(predInterval, file=file.path("../Output", "frequentistAnalysis", paste0(chan, "__POSTPRED.txt")),
              col.names=TRUE, row.names=FALSE, sep="\t")
}


write.table(dataWide, file=file.path("..", "Data", "data_prepped_freqClassif.txt"), 
            col.names=TRUE, row.names=FALSE, sep="\t")






