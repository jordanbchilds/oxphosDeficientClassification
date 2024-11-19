
library(analysis2Dmito)
library("parallel")
library("dplyr")
library("readr")
library("tidyr")
library("data.table")
library("rstan")
library("stringr")
library("coda")
library("mcmcse")

dir.create(file.path("..", "Output"), showWarnings = FALSE)

dataFolder = file.path("..", "Data", "synData01")
bayesFolder = file.path("..", "Output", "bayesInference_synData01")
maxESSFolder = file.path("..", "Output", "bayesMaxESS_synData01")
freqFolder = file.path("..", "Output", "frequentistAnalysis_synData01")

dir.create(bayesFolder, showWarnings = FALSE)
dir.create(maxESSFolder, showWarnings = FALSE)
dir.create(freqFolder, showWarnings = FALSE)

validate = FALSE

mitochan = "VDAC"
channels = c("NDUFB8", "CYB", "MTCO1")
nChan = length(channels)

ctrlIDs = paste0("C0", 1:4)
ptsIDs = c(paste0("P0", 1:7), "P09", paste0("P", 10:12))

# -----------------------
# --- BAYESIAN ANALYSIS
# -----------------------
dir.create(file.path(bayesFolder), showWarnings = FALSE)

nChains = 3

data_list = list()
for (chan in channels) {
  for (patID in ptsIDs) {
    fp = file.path(dataFolder, paste0(chan, "__", patID, "__synDATA.txt"))
    data = as.data.frame(fread(fp))
    
    if (validate) {
      validationObsIDs = data$obsID[which(data$validationSet==1)]
      data = data[!(data$obsID%in%validationObsIDs), ]
    }
    
    dataMats = getData_mats(
      data,
      pts = patID,
      channels = c(mitochan, chan),
      ctrlID = ctrlIDs,
      getIndex = TRUE
    )
    
    for (i in 1:nChains) {
      root = paste0(chan, "__", patID, "__chain", str_pad(i, width=2, side="left", pad="0"))
      data_list[[root]] = dataMats
    }
  }
}

ncores = 6
cl  = makeCluster(ncores)
{
  output = parLapply(
    cl,
    data_list,
    stan_inference, 
    warmup=10, #20000, 
    iter=15) #22000)
}
stopCluster(cl)

for (rt in names(output)) {
  list_saver(output[[rt]], file.path(bayesFolder, rt), rootSep="__")
}

# ---------------------------
# --- SAVE MAX ESS OUTPUT
# ---------------------------

for (chan in channels) {
  for (pat in ptsIDs) {
    root = paste0(chan, "__", pat)
    fpsChains = list.files(file.path(bayesFolder), pattern=root, full.names=TRUE)
    
    postFiles = grep("__POST.txt", fpsChains, value=TRUE)
    
    minESS = vector("numeric", length=length(postFiles))
    names(minESS) = postFiles
    multiESS = vector("numeric", length=length(postFiles))
    names(multiESS) = postFiles
    for (fp in postFiles) {
      post = fread(fp)
      minESS[fp] = min( effectiveSize(post) )
      # multiESS[fp] = multiESS( post )
    }
    fnKeep = names(multiESS)[which.max(multiESS)]
    rootKeep = gsub("__POST.txt", "", fnKeep)
    
    fpsKeep = grep(rootKeep, fpsChains, value=TRUE)
    
    file.copy(from=fpsKeep,
              to=file.path(maxESSFolder, gsub("__chain..", "", basename(fpsKeep))),
              overwrite=TRUE)
  }
}

# ---------------------------
# --- FREQUENTIST ANALYSIS
# ---------------------------

for (chan in channels) {
  for (patID in ptsIDs) {
    root = paste0(chan, "__", patID)
    fp = file.path(dataFolder, paste0(root, "__synDATA.txt"))
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
    
    write.table(dataLong, file=file.path(freqFolder, paste0(root, "__POST.txt")),
                col.names=TRUE, row.names=FALSE, sep="\t")
    
    write.table(synPredDF, file=file.path(freqFolder, paste0(root, "__PRED.txt")),
                col.names=TRUE, row.names=FALSE, sep="\t")
  }
}

# -----------------------
# --- COMBINE ANALYSIS
# -----------------------


BayesProbabilityCalculator = function (dataPoint, post, gamma=0.0001) {
  
  expected = dataPoint[1] * post[,"m[5]"] + post[,"c[5]"]
  densHealthy = dnorm(dataPoint[2], mean=expected, sd=1/sqrt(post[,"tau_norm"]))
  densDeficient = dnorm(dataPoint[2], mean=expected, sd=1/sqrt(gamma))
  
  top = post[,"probdiff"] * densDeficient
  bottom = (1 - post[,"probdiff"]) * densHealthy + top
  return ( top / bottom )
}

for (chan in channels) {
  for (pat in ptsIDs) {
    root = paste0(chan, "__", pat)
    
    dataFn = file.path(dataFolder, paste0(root, "__synDATA.txt"))
    postFn = file.path(maxESSFolder, paste0(root, "__POST.txt"))
    freqFn = file.path(freqFolder, paste0(root, "__POST.txt"))
    
    post = as.data.frame( fread ( postFn ) )
    freq = as.data.frame( fread( freqFn ) )
    data = as.data.frame( fread ( dataFn ) )
    
    data$meanPosteriorDeficientProb = NA
    data$varPosteriorDeficientProb = NA
    data$frequentistClassification = NA
    
    obsIDs = unique(data$obsID[data$sbjType=="patient"])
    
    for (id in obsIDs) {
      data$frequentistClassification[data$obsID==id & data$Channel==chan] = freq$frequentistClass[freq$obsID==id & freq$Channel==chan]
      
      dataFibre = data[data$obsID==id, ]
      dataPoint = c( dataFibre[dataFibre$Channel==mitochan, "Value"], 
                     dataFibre[dataFibre$Channel==chan, "Value"] )
      tt = BayesProbabilityCalculator(dataPoint, post)
      data[data$obsID==id & data$Channel==chan, "meanPosteriorDeficientProb"] = mean(tt)
      data[data$obsID==id & data$Channel==chan, "varPosteriorDeficientProb"] = var(tt)
    }
    
    write.table(data, file=gsub("synDATA", "postSummaryDATA", dataFn), 
                col.names=TRUE, row.names=FALSE, sep="\t")
  }
}



























