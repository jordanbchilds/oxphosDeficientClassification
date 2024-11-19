
library(analysis2Dmito)
library("parallel")
library("dplyr")
library("readr")
library("tidyr")
library("data.table")
library("rstan")
library("stringr")

dFlation =  5 # 5
tauFlation = 5 # 5

outputFolder = file.path("..", "Output", "bayesInference_widePrior")
dir.create(outputFolder, showWarnings = FALSE)

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

sbjIDs = unique(data$sampleID)
ctrlIDs = grep("C", sbjIDs, value=TRUE)
ptsIDs = sbjIDs[!(sbjIDs %in% ctrlIDs)]

nChains = 1 # 10

sampleList = list("NDUFB8" = ptsIDs,
                  "MTCO1" = ptsIDs, 
                  "CYB" = ptsIDs)

for (chan in channels) {
  ptsIDs = sampleList[[chan]]
  
  data_list = list()
  for (pat in ptsIDs) {
    dataMats = getData_mats(
      data,
      pts = pat,
      channels = c(mitochan, chan),
      ctrlID = ctrlIDs,
      getIndex = TRUE
    )
    
    indexID = unique(dataMats$indexCtrl)
    precEsts = vector("numeric", length=length(ctrlIDs))
    
    for (i in seq_along(indexID)) {
      mod = lm(dataMats$ctrl[dataMats$indexCtrl==indexID[i],2]~dataMats$ctrl[dataMats$indexCtrl==indexID[i],1])
      precEsts[i] = 1 / summary(mod)$sigma^2
    }
    
    prec_mu_m = 1 / (0.25^2 * dFlation)
    prec_mu_c = 1 / (0.25^2 * dFlation)
    
    tau_mode = mean(precEsts)
    tau_var = 10*tauFlation
    shape_tau = tau_mode^2 / tau_var
    rate_tau = tau_mode / tau_var
    
    tau_m_mode = 50
    tau_m_var = 50*tauFlation
    shape_tau_m = tau_m_mode^2 / tau_m_var
    rate_tau_m = tau_m_mode / tau_m_var
    
    tau_c_mode = 50
    tau_c_var = 50*tauFlation
    shape_tau_c = tau_c_mode^2 / tau_c_var
    rate_tau_c = tau_c_mode / tau_c_var
    
    for (i in 1:nChains) {
      root = paste0(chan, "__", pat, "__chain", str_pad(i, width=2, side="left", pad="0"))
      data_list[[root]] = dataMats
    }
  }
  
  ncores = detectCores() - 1
  cl  = makeCluster(ncores)
  {
    output = parLapply(
      cl,
      data_list,
      stan_inference, 
      warmup=10, # 20000, 
      iter=15, # 21000,
      parameterVals = list(prec_mu_m=prec_mu_m,
                           prec_mu_c=prec_mu_c, 
                           shape_tau=shape_tau,
                           rate_tau=rate_tau,
                           shape_tau_c=shape_tau_c,
                           rate_tau_c=rate_tau_c,
                           shape_tau_m=shape_tau_m,
                           rate_tau_m=rate_tau_m))
  }
  stopCluster(cl)
  
  for( rt in names(output) ){
    list_saver(output[[rt]], file.path(outputFolder, rt), rootSep="__")
  }
}


