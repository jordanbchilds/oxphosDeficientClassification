# install.packages(c("data.table", "dplyr", "readr", "tidyr", "plyr", "devtools"))
library("data.table")
library("dplyr")
library("readr")
library("tidyr")
library("devtools")
library("parallel")
library("stringr")

# install.packages("../../../analysis2Dmito", type="source", repos=NULL)
library("analysis2Dmito")

dir.create("PDF", showWarnings=FALSE)

mitochan = "VDAC"
channels = c("NDUFB8", "MTCO1", "CYB")
nChan = length(channels)
raw_data = read.csv(file.path("..", "..", "Data_prepped.csv"), header=TRUE)

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

postFolder = file.path("..", "OutputMaxESS", "stan_sampler")
outputFolder = file.path("..", "PDF", "OutputMaxESS", "stan_sampler")

dir.create(outputFolder, showWarnings=FALSE)

roots = list.files(file.path(postFolder), pattern="_POST.txt")

# ------ MCMC plot 
pdf(file.path("..", "PDF", outputFolder,  "MCMCplot.pdf"), width=13, height=8)
{
  for (rt in roots) {
    chanpat_str = substr(rt, 1, nchar(rt)-10 )
    chanpat_str = gsub("__chain..", "", chanpat_str)
    
    chanpat_lst = strsplit(chanpat_str, split="__")
    chan = chanpat_lst[[1]][1]
    pat = chanpat_lst[[1]][2]
    
    post = as.data.frame( fread( file.path(postFolder, rt) ) )
    post$lp__ = NULL
    
    fn_prior = file.path(postFolder, gsub("__POST.txt", "__PRIOR.txt", rt))
    prior = as.data.frame(fread(fn_prior))
    
    analysis2Dmito::MCMCplot(post, prior, nRow=3, lag=100, main_title=gsub("__POST.txt", "", rt))
  }
} 
dev.off()

# ------ postPlot
pdf(file.path("..", "PDF", outputFolder,  "postPlot.pdf"), width=13, height=8)
{
  for (rt in roots) {
    chanpat_str = substr(rt, 1, nchar(rt)-10 )
    chanpat_str = gsub("__chain..", "", chanpat_str)
    
    chanpat_lst = strsplit(chanpat_str, split="__")
    chan = chanpat_lst[[1]][1]
    pat = chanpat_lst[[1]][2]
    
    fn_post = file.path(postFolder, rt)
    post = as.data.frame( fread(fn_post) )
    
    fn_prior = file.path(postFolder, gsub("__POST.txt", "__PRIOR.txt", rt))
    prior = as.data.frame( fread(fn_prior) )
    
    fn_postpred = file.path(postFolder, gsub("__POST.txt", "__POSTPRED.txt", rt))
    postpred = as.data.frame( fread(fn_postpred) )
    
    fn_class = file.path(postFolder, gsub("__POST.txt", "__CLASSIF.txt", rt))
    class = apply( as.matrix(fread(fn_class) ), 2, mean)
    
    dataMats = getData_mats(data, channels=c(mitochan, chan), pts=pat, ctrlID=ctrlIDs)
    
    op = par(mfrow = c(3, 3), mar = c(4, 4, 1, 1), 
             cex.main = 2,cex.lab = 1.5,cex.axis = 1.5)
    analysis2Dmito::postPlot(
      post = post,
      prior = prior,
      postpred = postpred,
      classifs = class,
      dataMats = dataMats,
      var.names = c(
        "mu_m",
        "tau_m",
        "tau_norm",
        "mu_c",
        "tau_c",
        "probdiff",
        "m",
        "c"
      ),
      mitochan = paste0("log(", mitochan, ")"),
      chan = paste0("log(", chan, ")")
    )
    title(main=gsub("__POST.txt", "", rt), outer=TRUE, line=-1)
    par(op)
  }
}
dev.off()

### --- ---
### SAVE EXPECTED POSTERIOR PROBABILITIES OF DEFICIENCY PER MYOFIBRE
### --- ---

BayesProbabilityCalculator = function (dataPoint, post, gamma=0.0001) {
  
  expected = dataPoint[1] * post[,"m[5]"] + post[,"c[5]"]
  densHealthy = dnorm(dataPoint[2], mean=expected, sd=1/sqrt(post[,"tau_norm"]))
  densDeficient = dnorm(dataPoint[2], mean=expected, sd=1/sqrt(gamma))
  
  top = post[,"probdiff"] * densDeficient
  bottom = (1 - post[,"probdiff"]) * densHealthy + top
  return ( top / bottom )
}

dataClass = data
dataClass$obsID = paste0(dataClass$sampleID, "_", dataClass$fibreID)
dataClass$meanDefProb = NA

for (patID in ptsIDs) {
  for (chan in channels) {
    dataPat = dataClass[ (dataClass$Channel==chan | dataClass$Channel==mitochan) & dataClass$sampleID==patID, ]
    
    root = paste0(chan, "__", patID)
    
    post = as.data.frame( fread( file.path( postFolder, paste0(root, "__POST.txt")) ) )
    
    obsIDs = unique( dataPat$obsIDs )
    
    for (id in obsIDs) {
      dataFibre = dataPat[dataPat$obsID==id, ]
      dataPoint = c( dataFibre[dataFibre$Channel==mitochan, "Value"], 
                     dataFibre[dataFibre$Channel==chan, "Value"] )
      tt = BayesProbabilityCalculator(dataPoint, post)
      dataClass[dataClass$obsID==id & data$Channel==chan, "meanDefProb"] = mean(tt)
    }
  }
}

write.table(dataClass, file=file.path(postFolder, "data_prepped_withClassif.txt"), 
            sep="\t", col.names=TRUE, row.names=FALSE)






