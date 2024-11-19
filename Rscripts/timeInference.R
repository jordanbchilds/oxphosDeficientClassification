
###
# TIME INFERENCE
###

library(analysis2Dmito)
library("parallel")
library("dplyr")
library("readr")
library("tidyr")
library("data.table")
library("rstan")
library("stringr")


dir.create(file.path("..", "Output"), showWarnings = FALSE)

mitochan = "VDAC"
channels = c("NDUFB8", "CYB", "MTCO1")
nChan = length(channels)

raw_data = read.csv(file.path("..", "Data", "Data_prepped.csv"), header=TRUE)

raw_data = raw_data[,c("ID", "patient_id", mitochan, channels)]
colnames(raw_data) = c("fibreID", "sampleID", mitochan, channels)

data_lng = pivot_longer(raw_data, cols=c(mitochan, channels), 
                        names_to="Channel", values_to="Value")
data_lng = as.data.frame(data_lng)

data = data_lng
data$Value = log(data$Value)

sbjIDs = unique( data$sampleID )
ctrlIDs = grep("C0", sbjIDs, value=TRUE)
ptsIDs = sbjIDs[!(sbjIDs %in% ctrlIDs) ]

data_list = list()
for (chan in channels) {
  for (patID in ptsIDs) {
    root = paste0(chan, "__", patID)
    data_list[[root]] = getData_mats(
      data,
      pts = patID,
      channels = c(mitochan, chan),
      ctrlID = ctrlIDs,
      getIndex = TRUE
    )
  }
}

timeDF = data.frame(Channel = rep(channels, each=length(ptsIDs)), 
                    sampleID = rep(ptsIDs, 3))
timeDF$dataID = paste0(timeDF$Channel, "__", timeDF$sampleID)
timeDF$inferenceTime = NA

for (id in timeDF$dataID) {
  time = system.time({ stan_inference(data_list[[id]],
                                      stan_inference, 
                                      warmup=20000, 
                                      iter=22000) })
  timeDF[timeDF$dataID==id, "inferenceTime"] = time["elapsed"]
}

write.table(timeDF, file=file.path("..", "Output", "bayesInferenceTime.txt"), 
            col.names=TRUE, row.names=FALSE, sep="\t")
