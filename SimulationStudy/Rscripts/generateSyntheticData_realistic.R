

### --- create synthetic data for each patient in the dataset.


library(analysis2Dmito)
library("parallel")
library("dplyr")
library("readr")
library("tidyr")
library("data.table")
library("rstan")
library("stringr")
library("truncnorm")

outputFolder = file.path("..", "Data", "synCVData03")
parameterFolder = file.path("..", "Data", "trueParameters_synCVData03.txt")

dir.create(outputFolder, showWarnings = FALSE)

dataIncreaseFactor = 1 # 1.5
dataShiftFactor = 0 # 4
sdIncreaseFactor = 1 # 2


shiftDefFibresDown = 11
defFibresPrecisionShift = 10

mitochanTransform = function (x) {
  x*dataIncreaseFactor + dataShiftFactor
}

mitochan = "VDAC"
channels = c("NDUFB8", "CYB", "MTCO1")
nChan = length(channels)

raw_data = read.csv(file.path("..", "..", "..", "Data_prepped.csv"), header=TRUE)

raw_data = raw_data[,c("ID", "patient_id", mitochan, channels)]
colnames(raw_data) = c("fibreID", "sampleID", mitochan, channels)

data_lng = pivot_longer(raw_data, cols=c(mitochan, channels), 
                        names_to="Channel", values_to="Value")
data_lng = as.data.frame(data_lng)

data = data_lng
data$Value = log(data$Value)
data$sbjType = "patient"
data$sbjType[data$sampleID%in%c("C01", "C02", "C03", "C04")] = "control"

sbjIDs = unique(data$sampleID)
ctrlIDs = grep("C", sbjIDs, value=TRUE)
ptsIDs = sbjIDs[!(sbjIDs %in% ctrlIDs)]

nSamples = length(ctrlIDs) + 1
sampleSize = table(data$sampleID) / 4

set.seed(10101)

inferProportion = 1.0 # 0.8
validationProportion = 1 - inferProportion

trueParameters = data.frame(sampleID = rep(ptsIDs, each=length(channels)), 
                            Channel = rep(channels, length(ptsIDs)),
                            probdiff = NA, 
                            mu_m = NA, 
                            mu_c = NA, 
                            tau_m = NA, 
                            tau_c = NA, 
                            tau_norm = NA, 
                            tau_def = 0.0001)
for (i in 1:nSamples) {
  trueParameters[, paste0("m[", i, "]")] = NA
  trueParameters[, paste0("c[", i, "]")] = NA
}

for (patID in ptsIDs) {
  for (chan in channels) {
    # --- create a synthetic data frame for all sample-channel combinations  
    syntheticData = data[data$sampleID %in% c(ctrlIDs, patID) & data$Channel%in%c(mitochan, chan), ]
    syntheticData$obsID = paste0(syntheticData$sampleID, "_", syntheticData$fibreID)
    syntheticData$trueValue = syntheticData$Value
    syntheticData[syntheticData$Channel==chan, "Value"] = NA
    syntheticData$deficientStatus = 0
    
    # --- load posterior
    fn = file.path("..","..", "OutputMaxESS", "stan_sampler", paste0(chan, "__", patID, "__POST.txt"))
    post = fread(fn)
    # --- posterior expectation
    theta0 = colMeans(post)
    # --- save 
    parameterIndex = which(trueParameters$Channel==chan & trueParameters$sampleID==patID)
    
    # --- low level parameters
    # --- --- set to random draws from posterior
    trueParameters[parameterIndex, "mu_m"] = sample(post$mu_m, 1)
    trueParameters[parameterIndex, "mu_c"] = sample(post$mu_c, 1)
    trueParameters[parameterIndex, "tau_m"] = sample(post$tau_m, 1)
    trueParameters[parameterIndex, "tau_c"] = sample(post$tau_c, 1)
    trueParameters[parameterIndex, "tau_norm"] = sample(post$tau_norm, 1)
    trueParameters[parameterIndex, "probdiff"] = sample(post$probdiff, 1)
    
    # staticParameters = c("tau_norm", "mu_m", "mu_c", "tau_m", "tau_c", "probdiff")
    # trueParameters[parameterIndex, staticParameters] = theta0[staticParameters]
    
    trueParameters[parameterIndex, "tau_m"] = trueParameters[parameterIndex, "tau_m"] # / (sdIncreaseFactor^2)
    trueParameters[parameterIndex, "tau_c"] = trueParameters[parameterIndex, "tau_c"] # / (sdIncreaseFactor^2)
    trueParameters[parameterIndex, "tau_norm"] = trueParameters[parameterIndex, "tau_norm"] / (sdIncreaseFactor^2)
    
    # --- randomly draw higher level parameters depending on low level parameters  
    trueParameters[parameterIndex, paste0("m[", 1:nSamples, "]")] = rtruncnorm(nSamples, 
                                                                               a = 0.1,
                                                                               mean=trueParameters[parameterIndex, "mu_m"], 
                                                                               sd=1/sqrt(trueParameters[parameterIndex, "tau_m"]) )
    trueParameters[parameterIndex, paste0("c[", 1:nSamples, "]")] = rnorm(nSamples, 
                                                                          mean=trueParameters[parameterIndex, "mu_c"], 
                                                                          sd=1/sqrt(trueParameters[parameterIndex, "tau_c"]) )
    
    # simulate latent deficient state
    deficientStatus = rbinom(sampleSize[patID], 1, trueParameters[parameterIndex, "probdiff"] )
    names(deficientStatus) = syntheticData$obsID[syntheticData$Channel==chan & syntheticData$sampleID==patID]
    
    deficientFibresID = names(deficientStatus)[as.logical(deficientStatus)]
    healthyFibresID = names(deficientStatus)[!as.logical(deficientStatus)]
    
    syntheticData[syntheticData$obsID%in%deficientFibresID & syntheticData$Channel==chan, "deficientStatus"] = 1
    
    syntheticData$Value = mitochanTransform(syntheticData$trueValue )
    
    for (i in 1:length(ctrlIDs)) {
      xData = syntheticData[syntheticData$sampleID==ctrlIDs[i] & syntheticData$Channel==mitochan, "Value"]
      yMean = trueParameters[parameterIndex, paste0("m[", i, "]")] * xData + trueParameters[parameterIndex, paste0("c[", i, "]")]
      
      syntheticData[syntheticData$sampleID==ctrlIDs[i] & syntheticData$Channel==chan, "Value"] = rnorm(sampleSize[ctrlIDs[i]], 
                                                                                                       mean = yMean, 
                                                                                                       sd = 1 / sqrt(trueParameters[parameterIndex, "tau_norm"]))
    }
    
    xData = syntheticData[syntheticData$sampleID==patID & syntheticData$Channel==mitochan, "Value"]
    names(xData) = syntheticData[syntheticData$sampleID==patID & syntheticData$Channel==mitochan, "obsID"]
    
    yMean = trueParameters[parameterIndex, paste0("m[", nSamples, "]")] * xData + trueParameters[parameterIndex, paste0("c[", nSamples, "]")]
    
    yMean = ifelse(syntheticData[syntheticData$sampleID==patID & syntheticData$Channel==chan, "deficientStatus"], 
                   yMean - shiftDefFibresDown/sqrt(trueParameters[parameterIndex, "tau_norm"]), yMean)
    
    tauPat = ifelse(syntheticData[syntheticData$sampleID==patID & syntheticData$Channel==chan, "deficientStatus"], 
                    trueParameters[parameterIndex, "tau_norm"] / defFibresPrecisionShift, trueParameters[parameterIndex, "tau_norm"] )
    
    yNoise = rnorm(sampleSize[patID], 
                   mean = yMean,
                   sd = 1 / sqrt(tauPat))
    names(yNoise) = names(xData)
    
    syntheticData[!is.na(match(syntheticData$obsID, names(yNoise))) & syntheticData$Channel==chan, "Value"] = yNoise 
    
    syntheticData$validationSet = 0
        
    sampleIDs = unique(syntheticData$sampleID)
    
    for (sbjID in sampleIDs) {
      nSbj = sum(syntheticData$sampleID[syntheticData$Channel == mitochan] == sbjID)
      validation_obsID = sample(unique(syntheticData$obsID[syntheticData$sampleID == sbjID]), floor(validationProportion * nSbj))
      syntheticData$validationSet[syntheticData$obsID %in% validation_obsID] = 1
    }

    write.table(syntheticData, file=file.path(outputFolder, paste0(chan, "__", patID, "__", "synDATA.txt")),
                col.names=TRUE, row.names=FALSE, sep="\t")
  }
}

write.table(trueParameters, file=parameterFolder, 
            col.names=TRUE, row.names=FALSE, sep="\t")








