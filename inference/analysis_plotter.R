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

######################
### LET'S START HERE
######################

outputFolder = "OutputMaxESS"
folder = "Output_gamma0000001"

dir.create(file.path("PDF", outputFolder, folder), showWarnings=FALSE)

roots = list.files(file.path(outputFolder, folder), pattern="_POST.txt")

# ------ MCMC plot 
pdf(file.path("PDF", outputFolder, folder,  "MCMCplot.pdf"), width=13, height=8)
{
  for (rt in roots) {
    chanpat_str = substr(rt, 1, nchar(rt)-9 )
    chanpat_str = gsub("_chain..", "", chanpat_str)
    
    chanpat_lst = strsplit(chanpat_str, split="__")
    chan = chanpat_lst[[1]][1]
    pat = chanpat_lst[[1]][2]
    
    post = as.data.frame( fread( file.path(outputFolder, folder, rt) ) )
    post$lp__ = NULL
    
    fn_prior = file.path(outputFolder, folder, gsub("_POST.txt", "_PRIOR.txt", rt))
    prior = as.data.frame(fread(fn_prior))
    
    analysis2Dmito::MCMCplot(post, prior, nRow=3, lag=100, main_title=gsub("_POST.txt", "", rt))
  }
} 
dev.off()

# ------ postPlot
pdf(file.path("PDF", outputFolder, folder,  "postPlot.pdf"), width=13, height=8)
{
  for (rt in roots) {
    chanpat_str = substr(rt, 1, nchar(rt)-10 )
    chanpat_str = gsub("__chain..", "", chanpat_str)
    
    chanpat_lst = strsplit(chanpat_str, split="__")
    chan = chanpat_lst[[1]][1]
    pat = chanpat_lst[[1]][2]
    
    fn_post = file.path(outputFolder, folder, rt)
    post = as.data.frame( fread(fn_post) )
    
    fn_prior = file.path(outputFolder, folder, gsub("__POST.txt", "__PRIOR.txt", rt))
    prior = as.data.frame( fread(fn_prior) )
    
    fn_postpred = file.path(outputFolder, folder, gsub("__POST.txt", "__POSTPRED.txt", rt))
    postpred = as.data.frame( fread(fn_postpred) )
    
    fn_class = file.path(outputFolder, folder, gsub("__POST.txt", "__CLASSIF.txt", rt))
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
