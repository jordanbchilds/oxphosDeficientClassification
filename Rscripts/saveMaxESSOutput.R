# --- ESS 
library("tidyr")
library("data.table")
library("coda")
library("mcmcse")
library("stringr")

pathToOutput = file.path("..", "Output", "bayesInference")
pathToNewFile = file.path("..", "Output", "bayesMaxESS")

dir.create(pathToNewFile, showWarnings = FALSE)

fpsPost = list.files(pathToOutput, full.names=TRUE)

channels = c("CYB", "MTCO1", "NDUFB8")
ptsIDs = paste0("P", str_pad(width=2, side="left", pad="0", c(1:7, 9:12)))

for (chan in channels) {
  for (pat in ptsIDs) {
    root = paste0(chan, "__", pat)
    fpsChains = list.files(pathToOutput, pattern=root, full.names=TRUE)
    
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
    
    fpsMaxESS = grep(rootKeep, fpsChains, value=TRUE)
    fpsMove= gsub("__chain..", "", gsub(pathToOutput, pathToNewFile, fpsMaxESS))
    
    file.copy(from=fpsMaxESS, 
              to=fpsMove,
              overwrite=TRUE)
  }
}

