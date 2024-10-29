
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

folderName = file.path("..", "Output")
dataFolderPath = file.path("..", "Data", "synData")

dir.create(folderName, showWarnings = FALSE)

mitochan = "VDAC"
channels = c("NDUFB8", "CYB", "MTCO1")
nChan = length(channels)

ctrlIDs = paste0("C0", 1:4)
ptsIDs = c(paste0("P0", 1:7), "P09", paste0("P", 10:12))


data_list = list()
for (chan in channels) {
  for (patID in ptsIDs) {
    fp = file.path(dataFolderPath, paste0(chan, "__", patID, "__synData.txt"))
    
    data = as.data.frame(fread(fp))
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
timeDF$dataID = paste0(timeDF$Channel, "_", timeDF$sampleID)
timeDF$inferenceTime = NA

for (id in timeDF$dataID) {
  startTime = Sys.time()
  parLapply(
    cl,
    data_list[[id]],
    stan_inference, 
    warmup=20000, 
    iter=22000)
  endTime = Sys.time()
  timeDF[timeDF$dataID==id, "inferenceTime"] = endTime - startTime
}

write.table(timeDF, file=file.path("..", "inferenceTime.txt"), 
            col.names=TRUE, row.names=FALSE, sep="\t")
