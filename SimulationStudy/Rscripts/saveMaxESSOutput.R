# --- ESS 
library("tidyr")
library("data.table")
library("coda")
library("mcmcse")
library("stringr")

dataFolder = file.path("..", "Data", "synData03")
postFolder = file.path("..", "Output", "stan_sampler_CV_synData03")
maxESSFolder = file.path("..", "OutputMaxESS", "stan_sampler_CV_synData03")

dir.create(maxESSFolder, showWarnings = FALSE)

fpsPost = list.files(postFolder, full.names=TRUE)

channels = c("CYB", "MTCO1", "NDUFB8")
ptsIDs = paste0("P", str_pad(width=2, side="left", pad="0", c(1:7, 9:12)))

# ----------------------------
# --- save highest ESS chain
# ----------------------------
for (chan in channels) {
  for (pat in ptsIDs) {
    root = paste0(chan, "__", pat)
    fpsChains = list.files(file.path(postFolder), pattern=root, full.names=TRUE)
    postFiles = grep("__POST.txt", fpsChains, value=TRUE)
    
    minESS = vector("numeric", length=length(postFiles))
    names(minESS) = postFiles
    multiESS = vector("numeric", length=length(postFiles))
    names(multiESS) = postFiles
    for (fp in postFiles) {
      post = fread(fp)
      minESS[fp] = min( effectiveSize(post) )
      multiESS[fp] = multiESS( post )
    }
    fnKeep = names(multiESS)[which.max(multiESS)]
    rootKeep = gsub("__POST.txt", "", fnKeep)
    
    fpsKeep = grep(rootKeep, fpsChains, value=TRUE)
    
    file.copy(from=fpsKeep,
              to=gsub(paste0(postFolder, .Platform$file.sep), paste0(maxESSFolder, .Platform$file.sep), gsub("__chain..", "", fpsKeep)),
              overwrite=TRUE)
  }
}

# ---------------------------------------------------
# --- calc posterior prob of deficiency per myofibre
# ---------------------------------------------------

mitochan = "VDAC"




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
    postFn = file.path(maxESSFolder, paste0(root, "__POST.txt"))
    dataFn = file.path(dataFolder, paste0(root, "__synDATA.txt"))
    
    post = as.data.frame( fread ( postFn ) )
    data = as.data.frame( fread ( dataFn ) )
    
    data$meanPosteriorDeficientProb = NA
    data$varPosteriorDeficientProb = NA
    
    obsIDs = unique(data$obsID[data$sbjType=="patient"])
    
    for (id in obsIDs) {
      dataFibre = data[data$obsID==id, ]
      dataPoint = c( dataFibre[dataFibre$Channel==mitochan, "Value"], 
                     dataFibre[dataFibre$Channel==chan, "Value"] )
      tt = BayesProbabilityCalculator(dataPoint, post)
      data[data$obsID==id & data$Channel==chan, "meanPosteriorDeficientProb"] = mean(tt)
      data[data$obsID==id & data$Channel==chan, "varPosteriorDeficientProb"] = var(tt)
    }
    
    write.table(data, file=gsub("synDATA", "postDATA", dataFn), 
                col.names=TRUE, row.names=FALSE, sep="\t")
  }
}

