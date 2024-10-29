
mitochan = "VDAC"
channels = c("NDUFB8", "CYB", "MTCO1")
nChan = length(channels)

raw_data = read.csv(file.path("..", "..", "Data_prepped.csv"), header=TRUE)

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

BayesProbabilityCalculator = function (dataPoint, post, gamma=0.0001) {
  
  expected = dataPoint[1] * post[,"m[5]"] + post[,"c[5]"]
  densHealthy = dnorm(dataPoint[2], mean=expected, sd=sqrt(1/post[,"tau_norm"]))
  densDeficient = dnorm(dataPoint[2], mean=expected, sd=sqrt(1/gamma))
  
  top = post[,"probdiff"] * densDeficient
  bottom = (1 - post[,"probdiff"]) * densHealthy + top
  return ( top / bottom )
}

data$Classification = NA


freqClass_fn = file.path("..", "Output", "frequentist_linReg", "allData__CLASS.txt")
freqClass = as.data.frame( fread( freqClass_fn ) )

manualClass_fn = file.path("..", "..", "dat_with_class_prepped.txt")
manualClass = as.data.frame( fread( manualClass_fn ) )

### --- --- ---
### confusion matrix - FREQUENTIST
### --- --- ---
for (patID in ptsIDs) {
  for (chan in channels) {
    freqClass_tmp = freqClass[freqClass$sampleID==patID, c("sampleID", "fibreID", "obsID", "sbjType", chan, paste0(chan, "_class"))] # WIDE FROM
    manualClass_tmp = manualClass[manualClass$sampleID==patID & manualClass$Channel==chan, ]
    
    freqClass_patChan = freqClass_tmp[order(freqClass_tmp$fibreID), ]
    manualClass_patChan = manualClass_tmp[order(manualClass_tmp$fibreID), ]
    
    tt = table(freqClass_patChan[,paste0(chan, "_class")], manualClass_patChan$classJOINT)
    print (chan)
    print (tt)
    print ("-----------------")
    print ("-----------------")
    
  }
}

### --- --- ---
### confusion matrix - FREQUENTIST
### --- --- ---
bayesClass = as.data.frame( fread( file.path( postFolder, "data_prepped_withClassif.txt" ) ) )

for (patID in ptsIDs) {
  for (chan in channels) {
    bayesClass_tmp = bayesClass[bayesClass$sampleID==patID & bayesClass$Channel==chan, ] # WIDE FROM
    manualClass_tmp = manualClass[manualClass$sampleID==patID & manualClass$Channel==chan, ]
    
    bayesClass_patChan = bayesClass_tmp[order(bayesClass_tmp$fibreID), ]
    bayesClass_patChan$classif = bayesClass_patChan[,"meanDefProb"]>0.5
    manualClass_patChan = manualClass_tmp[order(manualClass_tmp$fibreID), ]
    
    tt = table(bayesClass_patChan$classif, manualClass_patChan$classJOINT)
    print (chan)
    print (tt)
    print ("-----------------")
    print ("-----------------")
    
  }
}

head(bayesClass_patChan)
dim(bayesClass_patChan)


