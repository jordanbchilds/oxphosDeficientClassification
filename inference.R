
library(analysis2Dmito)
library("parallel")
library("dplyr")
library("readr")
library("tidyr")
library("data.table")
library("rstan")
library("stringr")

folderName = "Output_gamma01"

dir.create(file.path(folderName), showWarnings = FALSE)

mitochan = "VDAC"
channels = c("NDUFB8", "CYB", "MTCO1")
nChan = length(channels)

raw_data = read.csv("../Data_prepped.csv", header=TRUE)

raw_data = raw_data[,c("ID", "patient_id", mitochan, channels)]
colnames(raw_data) = c("fibreID", "sampleID", mitochan, channels)

data_lng = pivot_longer(raw_data, cols=c(mitochan, channels), 
                        names_to="Channel", values_to="Value")
data_lng = as.data.frame(data_lng)

data = data_lng
data$Value = log(data$Value)

sbjIDs = unique(data$sampleID)
ctrlIDs = grep("C", sbjIDs, value=TRUE)
ptsIDs = sbjIDs[!(sbjIDs %in% ctrlIDs)]

nChains = 3

data_list = list()
for (chan in channels) {
  for (pat in ptsIDs) {
    for (i in 1:nChains) {
      root = paste0(chan, "__", pat, "__chain", str_pad(i, width=2, side="left", pad="0"))
      data_list[[root]] = getData_mats(
        data,
        pts = pat,
        channels = c(mitochan, chan),
        ctrlID = ctrlIDs,
        getIndex = TRUE
      )
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
    warmup=40000, 
    iter=45000,
    parameterVals = list(tau_def=0.1))
}
stopCluster(cl)

for( rt in names(output) ){
  list_saver(output[[rt]], file.path(folderName, rt), rootSep="__")
}




