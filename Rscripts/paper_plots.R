

# install.packages(c("data.table", "dplyr", "readr", "tidyr", "plyr", "devtools"))
library("data.table")
library("dplyr")
library("readr")
library("tidyr")
library("devtools")

install_github("jordanbchilds/analysis2Dmito")
library("analysis2Dmito")

# --- getting data

outputFolder = file.path("..", "OutputMaxESS", "stan_sampler")
outputFolderFreq = file.path("..", "Output", "frequentist_linReg")

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

# --- patient data example
pat = "P09"

tiff(file.path("..", "Figures",  "TIFF", "pat_example.tiff"), 
     width=2250, height=900, units="px", res=300, pointsize=13,
     type="cairo")
{
  pat = "P09"
  op = par(cex.lab=1.3, cex.axis=1.3, mfrow=c(1,3), mar=c(4,4,1,0.5))
  for( chan in channels ){
    xlm = range( data[data$sampleID%in% c(ctrlIDs, pat) & data$Channel==mitochan, "Value"] ) + c(0,0.5)
    ylm = range( data[data$sampleID%in% c(ctrlIDs, pat) & data$Channel!=mitochan, "Value"] ) + c(0,0.5)
    
    plot(NA, xlim=xlm, ylim=ylm, 
         xlab=paste0("log(", mitochan,")"), ylab=paste0("log(", chan, ")") )
    points( data[data$sampleID%in% ctrlIDs & data$Channel==mitochan, "Value"],
            data[data$sampleID%in% ctrlIDs & data$Channel==chan, "Value"],
            col=alphaBlack(0.05), bg=alphaBlack(0.05), 
            pch=21, cex=0.3)
    points( data[data$sampleID==pat & data$Channel==mitochan, "Value"],
            data[data$sampleID==pat & data$Channel==chan, "Value"],
            cex=0.5,
            pch=21, col=alphaCol(0.1, "hotpink"), bg=alphaCol(0.1, "hotpink"))
    
    if(chan == channels[1])
      legend("topleft", 
             col=c(alphaBlack(0.7), alphaCol(0.7, "hotpink")), 
             pch=20, 
             legend=c("Control", "patient"))
  }
  par(op)
}
dev.off()

png(file.path("..", "Figures",  "PNG", "pat_example.png"), 
     width=2250, height=900, units="px", res=300, pointsize=13,
    type="cairo")
{
  pat = "P09"
  op = par(cex.lab=1.3, cex.axis=1.3, mfrow=c(1,3), mar=c(4,4,1,0.5))
  for( chan in channels ){
    xlm = range( data[data$sampleID%in% c(ctrlIDs, pat) & data$Channel==mitochan, "Value"] ) + c(0,0.5)
    ylm = range( data[data$sampleID%in% c(ctrlIDs, pat) & data$Channel!=mitochan, "Value"] ) + c(0,0.5)
    
    plot(NA, xlim=xlm, ylim=ylm, 
         xlab=paste0("log(", mitochan,")"), ylab=paste0("log(", chan, ")") )
    points( data[data$sampleID%in% ctrlIDs & data$Channel==mitochan, "Value"],
            data[data$sampleID%in% ctrlIDs & data$Channel==chan, "Value"],
            col=alphaBlack(0.05), bg=alphaBlack(0.05), 
            pch=21, cex=0.3)
    points( data[data$sampleID==pat & data$Channel==mitochan, "Value"],
            data[data$sampleID==pat & data$Channel==chan, "Value"],
            cex=0.5,
            pch=21, col=alphaCol(0.1, "hotpink"), bg=alphaCol(0.1, "hotpink"))
    
    if(chan == channels[1])
      legend("topleft", 
             col=c(alphaBlack(0.7), alphaCol(0.7, "hotpink")), 
             pch=20, 
             legend=c("Control", "patient"))
  }
  par(op)
}
dev.off()

# ------------------------------
# --- frequentist classification
# ------------------------------
patID = "P09"

mdat = as.data.frame( fread( file.path("..", "..", "dat_with_class_prepped.txt") ) )

mdat_pat = mdat[ mdat$sampleID==patID, ]

tiff(file.path("..", "Figures",  "TIFF", "freq_example.tiff"),
     width=2250, height=900, units="px", res=300, pointsize=13,
     type="quartz")
{
  pat = "P09"
  colVec = c(alphaBlue(0.5), alphaRed(0.5))
  
  op = par(cex.lab=1.3, cex.axis=1.3, mfrow=c(1,3), mar=c(4,4,1,0.5))
  for( chan in channels ){
    freq_class = as.data.frame( fread( file.path(outputFolderFreq, paste0(chan, "__", pat, "__CLASSIF.txt") ) ) )
    freq_postpred = as.data.frame( fread(file.path(outputFolderFreq, paste0(chan, "__POSTPRED.txt" ))))
    
    xlm = range( data[data$sampleID%in% c(ctrlIDs, pat) & data$Channel==mitochan, "Value"] ) + c(0,0.5)
    ylm = range( data[data$sampleID%in% c(ctrlIDs, pat) & data$Channel!=mitochan, "Value"] ) + c(0,0.5)
    
    pi_est = round(colMeans(freq_class), 3)
    
    plot(NA, xlim=xlm, ylim=ylm, 
         xlab=paste0("log(", mitochan,")"), ylab=paste0("log(", chan, ")"),
         main=bquote(pi*"="*.(pi_est)))
    points( data[data$sampleID%in% ctrlIDs & data$Channel==mitochan, "Value"],
            data[data$sampleID%in% ctrlIDs & data$Channel==chan, "Value"],
            pch=21, col=alphaBlack(0.05), bg=alphaBlack(0.05), cex=0.5)
    
    points( data[data$sampleID==pat & data$Channel==mitochan, "Value"],
            data[data$sampleID==pat & data$Channel==chan, "Value"],
            pch=21, col=colVec[freq_class[[1]]+1], bg=colVec[freq_class[[1]]+1],
            cex=0.5)
    
    manDef_cellID = mdat_pat[ mdat_pat$Channel==chan & mdat_pat$classJOINT==1, "fibreID"]
    
    manDef_x = mdat_pat[ mdat_pat$cell_id%in%manDef_cellID & mdat_pat$Channel==mitochan, "Value"]
    manDef_y = mdat_pat[ mdat_pat$cell_id%in%manDef_cellID & mdat_pat$Channel==chan, "Value"]
    
    points( log(manDef_x), log(manDef_y), pch=20, cex=0.3, col="yellow")
    
    lines(freq_postpred$mitochan, freq_postpred$fit, 
          col=alphaGreen(1.0), lty=1, lwd=2)
    lines(freq_postpred$mitochan, freq_postpred$lwr, 
          col=alphaGreen(1.0), lty=2, lwd=2)
    lines(freq_postpred$mitochan, freq_postpred$upr, 
          col=alphaGreen(1.0), lty=2, lwd=2)
    
    if( chan == channels[1] ) {
      legend("topleft", 
             col=c(alphaBlack(0.7), alphaBlue(0.7), alphaRed(0.7), alphaCol(0.7, name="yellow")),
             pch=20, 
             legend=c("Control", "Healthy", "Deficient", "Manual"),
             bg="transparent", 
             box.col=alphaBlack(0.0))
    }
    
    
  }
  par(op)
}
dev.off()

png(file.path("..", "Figures",  "PNG", "freq_example.png"), 
     width=2250, height=900, units="px", res=300, pointsize=13,
     type="quartz")
{{
  pat = "P09"
  colVec = c(alphaBlue(0.5), alphaRed(0.5))
  
  op = par(cex.lab=1.3, cex.axis=1.3, mfrow=c(1,3), mar=c(4,4,1,0.5))
  for( chan in channels ){
    freq_class = as.data.frame( fread( file.path(outputFolderFreq, paste0(chan, "__", pat, "__CLASSIF.txt") ) ) )
    freq_postpred = as.data.frame( fread(file.path(outputFolderFreq, paste0(chan, "__POSTPRED.txt" ))))
    
    xlm = range( data[data$sampleID%in% c(ctrlIDs, pat) & data$Channel==mitochan, "Value"] ) + c(0,0.5)
    ylm = range( data[data$sampleID%in% c(ctrlIDs, pat) & data$Channel!=mitochan, "Value"] ) + c(0,0.5)
    
    pi_est = round(colMeans(freq_class), 3)
    
    plot(NA, xlim=xlm, ylim=ylm, 
         xlab=paste0("log(", mitochan,")"), ylab=paste0("log(", chan, ")"),
         main=bquote(pi*"="*.(pi_est)))
    points( data[data$sampleID%in% ctrlIDs & data$Channel==mitochan, "Value"],
            data[data$sampleID%in% ctrlIDs & data$Channel==chan, "Value"],
            pch=21, col=alphaBlack(0.05), bg=alphaBlack(0.05), cex=0.5)
    
    points( data[data$sampleID==pat & data$Channel==mitochan, "Value"],
            data[data$sampleID==pat & data$Channel==chan, "Value"],
            pch=21, col=colVec[freq_class[[1]]+1], bg=colVec[freq_class[[1]]+1],
            cex=0.5)
    
    manDef_cellID = mdat_pat[ mdat_pat$Channel==chan & mdat_pat$classJOINT==1, "fibreID"]
    
    manDef_x = mdat_pat[ mdat_pat$cell_id%in%manDef_cellID & mdat_pat$Channel==mitochan, "Value"]
    manDef_y = mdat_pat[ mdat_pat$cell_id%in%manDef_cellID & mdat_pat$Channel==chan, "Value"]
    
    points( log(manDef_x), log(manDef_y), pch=20, cex=0.3, col="yellow")
    
    lines(freq_postpred$mitochan, freq_postpred$fit, 
          col=alphaGreen(1.0), lty=1, lwd=2)
    lines(freq_postpred$mitochan, freq_postpred$lwr, 
          col=alphaGreen(1.0), lty=2, lwd=2)
    lines(freq_postpred$mitochan, freq_postpred$upr, 
          col=alphaGreen(1.0), lty=2, lwd=2)
    
    if( chan == channels[1] ) {
      legend("topleft", 
             col=c(alphaBlack(0.7), alphaBlue(0.7), alphaRed(0.7), alphaCol(0.7, name="yellow")),
             pch=20, 
             legend=c("Control", "Healthy", "Deficient", "Manual"),
             bg="transparent", 
             box.col=alphaBlack(0.0))
    }
    
    
  }
  par(op)
}}
dev.off()

# --------------------------------
# --- bayes classification example
# --------------------------------
pat = "P09"

mdat = as.data.frame( fread( file.path("..", "..", "dat_with_class_prepped.txt") ) )
mdat_pat = mdat[ mdat$patient_id==pat, ]

tiff(file.path("..", "Figures",  "TIFF", "bayes_example.tiff"), 
     width=2250, height=900, units="px", res=300, pointsize=13, 
     type="cairo")
{
  pat = "P09"
  op = par(cex.lab=1.3, cex.axis=1.3, mfrow=c(1,3), mar=c(4,4,1,0.5))
  for (chan in channels) {
    bayes_class_mat = as.matrix( fread( file.path(outputFolder , paste0(chan, "__", pat, "__CLASSIF.txt") )) )
    bayes_class = colMeans( bayes_class_mat )
    
    bayes_postpred = as.data.frame(fread(file.path(outputFolder, paste0(chan, "__", pat, "__POSTPRED.txt") )))
    
    xlm = range( data[data$sampleID%in% c(ctrlIDs, pat) & data$Channel==mitochan, "Value"] ) + c(0,0.5)
    ylm = range( data[data$sampleID%in% c(ctrlIDs, pat) & data$Channel!=mitochan, "Value"] ) + c(0,0.5)
    
    pi_est = round( colMeans( fread( file.path(outputFolder, paste0(chan, "__", pat, "__POST.txt" )))[, "probdiff"] ), 3)
    
    plot(NA, xlim=xlm, ylim=ylm, 
         xlab=paste0("log(", mitochan,")"), ylab=paste0("log(", chan, ")"),
         main=bquote("E("*pi*"|Y)="*.(pi_est)))
    points( data[data$sampleID%in% ctrlIDs & data$Channel==mitochan, "Value"],
            data[data$sampleID%in% ctrlIDs & data$Channel==chan, "Value"],
            pch=20, col=alphaBlack(0.05), cex=0.5)
    
    points( data[data$sampleID==pat & data$Channel==mitochan, "Value"],
            data[data$sampleID==pat & data$Channel==chan, "Value"],
            pch=20, col=classcols(bayes_class, alphaLevel=0.5), cex=0.5)
    
    manDef_cellID = mdat_pat[ mdat_pat$channel==chan & mdat_pat$classJOINT==1, "fibreID"]
    
    manDef_x = mdat_pat[ mdat_pat$fibreID%in%manDef_cellID & mdat_pat$Channel==mitochan, "Value"]
    manDef_y = mdat_pat[ mdat_pat$fibreID%in%manDef_cellID & mdat_pat$Channel==chan, "Value"]
    
    points( log(manDef_x), log(manDef_y), pch=20, cex=0.3, col="yellow")
    
    lines(bayes_postpred$mitochan, bayes_postpred$medNorm, 
          col=alphaGreen(1.0), lty=1, lwd=2)
    lines(bayes_postpred$mitochan, bayes_postpred$lwrNorm, 
          col=alphaGreen(1.0), lty=2, lwd=2)
    lines(bayes_postpred$mitochan, bayes_postpred$uprNorm, 
          col=alphaGreen(1.0), lty=2, lwd=2)
    
    if( chan == channels[1] )
      legend("topleft", 
             col=c(alphaBlack(0.7), alphaBlue(0.7), alphaRed(0.7), alphaCol(0.7, name="yellow")),
             pch=20, 
             legend=c("Control", "Healthy", "Deficient", "Manual"),
             bg="transparent", 
             box.col=alphaBlack(0.0))
  }
  par(op)
}
dev.off()

png(file.path("..", "Figures",  "PNG", "bayes_example.png"),
     width=2250, height=900, units="px", res=300, pointsize=13, 
     type="cairo")
{{
  pat = "P09"
  op = par(cex.lab=1.3, cex.axis=1.3, mfrow=c(1,3), mar=c(4,4,1,0.5))
  for (chan in channels) {
    bayes_class_mat = as.matrix( fread( file.path(outputFolder , paste0(chan, "__", pat, "__CLASSIF.txt") )) )
    bayes_class = colMeans( bayes_class_mat )
    
    bayes_postpred = as.data.frame(fread(file.path(outputFolder, paste0(chan, "__", pat, "__POSTPRED.txt") )))
    
    xlm = range( data[data$sampleID%in% c(ctrlIDs, pat) & data$Channel==mitochan, "Value"] ) + c(0,0.5)
    ylm = range( data[data$sampleID%in% c(ctrlIDs, pat) & data$Channel!=mitochan, "Value"] ) + c(0,0.5)
    
    pi_est = round( colMeans( fread( file.path(outputFolder, paste0(chan, "__", pat, "__POST.txt" )))[, "probdiff"] ), 3)
    
    plot(NA, xlim=xlm, ylim=ylm, 
         xlab=paste0("log(", mitochan,")"), ylab=paste0("log(", chan, ")"),
         main=bquote("E("*pi*"|Y)="*.(pi_est)))
    points( data[data$sampleID%in% ctrlIDs & data$Channel==mitochan, "Value"],
            data[data$sampleID%in% ctrlIDs & data$Channel==chan, "Value"],
            pch=20, col=alphaBlack(0.05), cex=0.5)
    
    points( data[data$sampleID==pat & data$Channel==mitochan, "Value"],
            data[data$sampleID==pat & data$Channel==chan, "Value"],
            pch=20, col=classcols(bayes_class, alphaLevel=0.5), cex=0.5)
    
    manDef_cellID = mdat_pat[ mdat_pat$channel==chan & mdat_pat$classJOINT==1, "fibreID"]
    
    manDef_x = mdat_pat[ mdat_pat$fibreID%in%manDef_cellID & mdat_pat$Channel==mitochan, "Value"]
    manDef_y = mdat_pat[ mdat_pat$fibreID%in%manDef_cellID & mdat_pat$Channel==chan, "Value"]
    
    points( log(manDef_x), log(manDef_y), pch=20, cex=0.3, col="yellow")
    
    lines(bayes_postpred$mitochan, bayes_postpred$medNorm, 
          col=alphaGreen(1.0), lty=1, lwd=2)
    lines(bayes_postpred$mitochan, bayes_postpred$lwrNorm, 
          col=alphaGreen(1.0), lty=2, lwd=2)
    lines(bayes_postpred$mitochan, bayes_postpred$uprNorm, 
          col=alphaGreen(1.0), lty=2, lwd=2)
    
    if( chan == channels[1] )
      legend("topleft", 
             col=c(alphaBlack(0.7), alphaBlue(0.7), alphaRed(0.7), alphaCol(0.7, name="yellow")),
             pch=20, 
             legend=c("Control", "Healthy", "Deficient", "Manual"),
             bg="transparent", 
             box.col=alphaBlack(0.0))
  }
  par(op)
}}
dev.off()

# -------------------------------
# --- prior-posteriors all in one
# -------------------------------
tiff(file.path("..", "Figures",  "TIFF", "prior_post.tiff"), 
     width=2250, height=600*4, units="px", res=300, pointsize=13, 
     type="cairo")
{
  pat = "P09"
  chan = "NDUFB8"
  
  op = par(cex.lab=1.3, cex.axis=1.3, mfrow=c(4,2), mar=c(4,4,1,0.5))
  
  post = as.data.frame( fread( file.path(outputFolder, paste0(chan, "__", pat, "__POST.txt") ) ) )
  prior = as.data.frame( fread( file.path(outputFolder, paste0(chan, "__", pat, "__PRIOR.txt") ) ) )
  
  xlb = list(probdiff=bquote(pi), tau_norm=bquote(tau))
  names(xlb) = c("probdiff", "tau_norm")
  for(para in c("probdiff", "tau_norm")){
    dpost = density(post[,para])
    dprior = density(prior[,para])
    plot(NA, 
         xlim=range(c(dpost$x, dprior$x)), ylim=range(c(dpost$y, dprior$y)),
         ylab="Density", xlab=xlb[[para]]
    )
    if(para == "probdiff") curve(dunif(x, 0, 0.5), col=alphaPink(0.9), add=TRUE, lwd=2)
    else lines(dprior, col=alphaPink(0.9), lwd=2)
    lines(dpost, col=alphaGreen(0.9), lwd=2 )
  }
  
  xlb = list(mu_m=bquote(mu[m]), mu_c=bquote(mu[c]), tau_m=bquote(tau[m]), tau_c=bquote(tau[c]))
  for(para in c("mu_m", "mu_c", "tau_m", "tau_c")){
    dpost = density(post[,para])
    dprior = density(prior[,para])
    plot(NA, 
         xlim=range(c(dpost$x, dprior$x)), ylim=range(c(dpost$y, dprior$y)),
         ylab="Density", xlab=xlb[[para]]
    )
    lines(dprior, col=alphaPink(0.9), lwd=2)
    lines(dpost, col=alphaGreen(0.9), lwd=2)
  }
    xlb = list(m="Slope", c="Intercept")
    for(para in c("m", "c")){
      xlm = NULL
      ylm = NULL
      
      dpost_list = list()
      for(i in 1:(length(ctrlIDs)+1) ){
        dpost_list[[i]] = density( post[, paste0(para, "[", i, "]")] )
        xlm = range( c(xlm, range(dpost_list[[i]]$x) ))
        ylm = range( c(ylm, range(dpost_list[[i]]$y) ))
      }
      
      plot(NA, 
           xlim=xlm, ylim=ylm,
           ylab="Density", xlab=xlb[[para]]
      )
      lines(density(prior[,paste0(para, "_pred")]), lty=3, col=alphaPink(1.0), lwd=2)
      lines(density(post[,paste0(para, "_pred")]), lty=3, col=alphaGreen(1.0), lwd=2)
      
      for(i in 1:length(ctrlIDs) ){
        lines( dpost_list[[i]],  lty=1, col=alphaGreen(0.4), lwd=1)
      }
      lines( dpost_list[[length(ctrlIDs)+1]], lty=1, col=alphaGreen(1.0), lwd=2)
    }
  par(op)
}
dev.off()

png(file.path("..", "Figures",  "PNG", "prior_post.png"),
     width=2250, height=600*4, units="px", res=300, pointsize=13, 
     type="cairo")
{
  pat = "P09"
  chan = "NDUFB8"
  
  op = par(cex.lab=1.3, cex.axis=1.3, mfrow=c(4,2), mar=c(4,4,1,0.5))
  
  post = as.data.frame( fread( file.path(outputFolder, paste0(chan, "__", pat, "__POST.txt") ) ) )
  prior = as.data.frame( fread( file.path(outputFolder, paste0(chan, "__", pat, "__PRIOR.txt") ) ) )
  
  xlb = list(probdiff=bquote(pi), tau_norm=bquote(tau))
  names(xlb) = c("probdiff", "tau_norm")
  for(para in c("probdiff", "tau_norm")){
    dpost = density(post[,para])
    dprior = density(prior[,para])
    plot(NA, 
         xlim=range(c(dpost$x, dprior$x)), ylim=range(c(dpost$y, dprior$y)),
         ylab="Density", xlab=xlb[[para]]
    )
    if(para == "probdiff") curve(dunif(x, 0.0, 0.5), col=alphaPink(0.9), add=TRUE, lwd=2)
    else lines(dprior, col=alphaPink(0.9), lwd=2)
    lines(dpost, col=alphaGreen(0.9), lwd=2 )
  }
  
  xlb = list(mu_m=bquote(mu[m]), mu_c=bquote(mu[c]), tau_m=bquote(tau[m]), tau_c=bquote(tau[c]))
  for(para in c("mu_m", "mu_c", "tau_m", "tau_c")){
    dpost = density(post[,para])
    dprior = density(prior[,para])
    plot(NA, 
         xlim=range(c(dpost$x, dprior$x)), ylim=range(c(dpost$y, dprior$y)),
         ylab="Density", xlab=xlb[[para]]
    )
    lines(dprior, col=alphaPink(0.9), lwd=2)
    lines(dpost, col=alphaGreen(0.9), lwd=2)
  }
  xlb = list(m="Slope", c="Intercept")
  for(para in c("m", "c")){
    xlm = NULL
    ylm = NULL
    
    dpost_list = list()
    for(i in 1:(length(ctrlIDs)+1) ){
      dpost_list[[i]] = density( post[, paste0(para, "[", i, "]")] )
      xlm = range( c(xlm, range(dpost_list[[i]]$x) ))
      ylm = range( c(ylm, range(dpost_list[[i]]$y) ))
    }
    
    plot(NA, 
         xlim=xlm, ylim=ylm,
         ylab="Density", xlab=xlb[[para]]
    )
    lines(density(prior[,paste0(para, "_pred")]), lty=3, col=alphaPink(1.0), lwd=2)
    lines(density(post[,paste0(para, "_pred")]), lty=3, col=alphaGreen(1.0), lwd=2)
    
    for(i in 1:length(ctrlIDs) ){
      lines( dpost_list[[i]],  lty=1, col=alphaGreen(0.4), lwd=1)
    }
    lines( dpost_list[[length(ctrlIDs)+1]], lty=1, col=alphaGreen(1.0), lwd=2)
  }
  par(op)
}
dev.off()


# --- better comparison of classifs
mdat = as.data.frame( fread( file.path("..", "..", "dat_with_class_prepped.txt") ) )

bayes_diff = list()
freq_diff = list()

postprob = vector("numeric", length=length(channels)*length(ptsIDs))
postprob_names = apply(expand.grid(channels, ptsIDs), 1, paste, collapse = "__")
names(postprob) = postprob_names

for (pat in ptsIDs) {
  for( chan in channels ){
    root = paste0(chan, "__", pat)
    
    fn_post = file.path("..", folder, outputFolder, paste0(root, "__POST.txt"))
    probdiff_post = fread( fn_post )$probdiff
    
    man_est = mean( mdat[ mdat$sampleID==pat & mdat$Channel==chan, ]$classJOINT  )
    
    diff_vec = probdiff_post - man_est
    
    fn_freq = file.path("..", "Output", "frequentist_linReg", paste0(root, "__CLASSIF.txt"))
    freq_class = fread( fn_freq )[[1]]
    
    bayes_diff[[root]] = diff_vec
    freq_diff[[root]] = mean(freq_class) - man_est
    postprob[root] = sum( diff_vec > freq_diff[[root]] ) / length(diff_vec)
  }
}

mean( unlist(freq_diff) )
mean( unlist(lapply(bayes_diff, median) ))

tiff(file.path("..", "Figures",  "TIFF", "proportion_comparison.tiff"), 
     width=2250, height=1900, units="px", res=300, pointsize=13,
     type="cairo")
{
  op = par(cex.lab=1.3, cex.axis=1.3, mar=c(4,4,1,0.5), xpd=FALSE)
  
  colVec = c(alphaCol(0.01, name="lightseagreen"), 
             alphaCol(0.01, name="deeppink3"), 
             alphaCol(0.01, name="darkorange"))
  
  myAt = c(1:44)[ 1:44 %% 4 !=0 ]
  myTicks = 4* 1:11 - 2
  plot(NA, xlim=c(1,43), ylim=c(-0.1,1), xaxt="n", 
       ylab="Difference from manual classification", xlab="Patient")
  
  for(i in seq(0,10, by=2)) rect( xleft=i*4, xright=(i+1)*4, ybottom=par("usr")[3], ytop=par("usr")[4], border=NA, col=alphaBlack(0.1))
  lines(x=c(0,44), c(0,0), lty=3, col="black", lwd=3)
  
  stripchart(bayes_diff, vertical=TRUE, pch=20, cex=0.5, 
             method="jitter", main="", ylim=c(0,1), xlim=c(1,43),
             xlab="Patient", ylab="Deficient Proportion", 
             col=rep(colVec, length(bayes_diff)), 
             at=myAt, xaxt="n", add=TRUE)
  axis(side=1, at=myTicks, labels=ptsIDs)
  points(myAt, unlist(freq_diff), pch=24, cex=1.2,
         col="black",
         bg=rep(c(alphaCol(0.9, name="lightseagreen"), 
                   alphaCol(0.9, name="deeppink3"), 
                   alphaCol(0.9, name="darkorange")), length(freq_diff)))
  #points(x=myAt[unlist(nonZero)], y=rep(-0.1, sum(unlist(nonZero))), pch=8)
  legend(x=0, y=1, 
         bg="transparent",
         box.col=alphaBlack(0.0),
         legend=channels, 
         text.width = 10,
         pch=20, 
         col=c(alphaCol(0.7, name="lightseagreen"), 
               alphaCol(0.7, name="deeppink3"), 
               alphaCol(0.7, name="darkorange")))
  text(x=0, y=1.0, labels=c("Bayesian distribution"), pos=4)
  
  legend(x=12, y=1.0, 
         bg="transparent", 
         box.col=alphaBlack(0.0),
         legend=channels, 
         text.width = 10,
         pch=24, 
         col="black",
         pt.bg=c(alphaCol(0.7, name="lightseagreen"), 
               alphaCol(0.7, name="deeppink3"), 
               alphaCol(0.7, name="darkorange")) )
  text(x=12, y=1.0, labels=c("Frequentist point estimate"), pos=4)

  par(op)
}
dev.off()

png(file.path("..", "Figures",  "PNG", "proportion_comparison.png"),
     width=2250, height=1900, units="px", res=300, pointsize=13,
     type="cairo")
{
  op = par(cex.lab=1.3, cex.axis=1.3, mar=c(4,4,1,0.5), xpd=FALSE)
  
  colVec = c(alphaCol(0.01, name="lightseagreen"), 
             alphaCol(0.01, name="deeppink3"), 
             alphaCol(0.01, name="darkorange"))
  
  myAt = c(1:44)[ 1:44 %% 4 !=0 ]
  myTicks = 4* 1:11 - 2
  plot(NA, xlim=c(1,43), ylim=c(-0.1,1), xaxt="n", 
       ylab="Difference from manual classification", xlab="Patient")
  
  for(i in seq(0,10, by=2)) rect( xleft=i*4, xright=(i+1)*4, ybottom=par("usr")[3], ytop=par("usr")[4], border=NA, col=alphaBlack(0.1))
  lines(x=c(0,44), c(0,0), lty=3, col="black", lwd=3)
  
  stripchart(bayes_diff, vertical=TRUE, pch=20, cex=0.5, 
             method="jitter", main="", ylim=c(0,1), xlim=c(1,43),
             xlab="Patient", ylab="Deficient Proportion", 
             col=rep(colVec, length(bayes_diff)), 
             at=myAt, xaxt="n", add=TRUE)
  axis(side=1, at=myTicks, labels=ptsIDs, cex.axis=1.0)
  points(myAt, unlist(freq_diff), pch=24, cex=1.2,
         col="black",
         bg=rep(c(alphaCol(0.9, name="lightseagreen"), 
                  alphaCol(0.9, name="deeppink3"), 
                  alphaCol(0.9, name="darkorange")), length(freq_diff)))
  #points(x=myAt[unlist(nonZero)], y=rep(-0.1, sum(unlist(nonZero))), pch=8)
  legend(x=0, y=1, 
         bg="transparent",
         box.col=alphaBlack(0.0),
         legend=channels, 
         text.width = 10,
         pch=20, 
         col=c(alphaCol(0.7, name="lightseagreen"), 
               alphaCol(0.7, name="deeppink3"), 
               alphaCol(0.7, name="darkorange")))
  text(x=0, y=1.0, labels=c("Bayesian distribution"), pos=4)
  
  legend(x=12, y=1.0, 
         bg="transparent", 
         box.col=alphaBlack(0.0),
         legend=channels, 
         text.width = 10,
         pch=24, 
         col="black",
         pt.bg=c(alphaCol(0.7, name="lightseagreen"), 
                 alphaCol(0.7, name="deeppink3"), 
                 alphaCol(0.7, name="darkorange")) )
  text(x=12, y=1.0, labels=c("Frequentist point estimate"), pos=4)
  
  par(op)
}
dev.off()

# --- gamma comparison - MEAN ABSOLUTE DIFFERENCE
mdat = as.data.frame( fread( file.path("../dat_with_class_prepped.txt") ) )

mad = list()
gamma_labs = c("0000001", "000001", "00001", "0001", "001", "01", "10")
gamma_vals = c("0.000001", "0.00001", "0.0001", "0.001", "0.01", "0.1", "1.0")

for (gam in gamma_labs) {
  fldr = paste0("Output_gamma", gam)
  diffs = NULL
  for( pat in ptsIDs ){
    for( chan in channels ){
      mclass = unlist( mdat[ mdat$channel==chan & mdat$patient_id==pat, "classJBC"] )
      root = paste0(chan, "__", pat)
      fn_post = file.path("OutputMaxESS", fldr, paste0(root, "__POST.txt"))
      
      post = as.matrix( fread( fn_post) )
      diffs = c(diffs, abs(post[,"probdiff"] - mean(mclass)))
    }
  }
  mad[[gam]] = mean(diffs)
}

tiff(file.path("..", "Figures",  "TIFF", "gamma_comparison.tiff"),
    width=2250, height=1000, units="px", res=300, pointsize=13, 
    compression="lzw")
{
  op = par(cex.lab=1.3, cex.axis=1.3, mfrow=c(1,1), mar=c(4,4,1,0.5))
  plot(1:length(gamma_vals), unlist(mad), ylim=c(0,0.1),
       ylab="Mean absolute difference", xlab=bquote(gamma), 
       xaxt='n', type='b', 
       pch=20)
  axis(side=1, at=1:length(gamma_vals), labels=gamma_vals)
  
  par(op)
}
dev.off()

png(file.path("..", "Figures",  "PNG", "gamma_comparison.png"), 
     width=2250, height=1000, units="px", res=300, pointsize=13, 
     type="cairo")
{
  op = par(cex.lab=1.0, cex.axis=1.0, mfrow=c(1,1), mar=c(4,4,1,0.5))
  plot(1:length(gamma_vals), unlist(mad), ylim=c(0,0.1),
       ylab="Mean absolute difference", xlab=bquote(gamma), 
       xaxt='n', type='b', 
       pch=20)
  axis(side=1, at=1:length(gamma_vals), labels=gamma_vals)
  
  par(op)
}
dev.off()

# ----------------------------------------
# --- prior comparison against OG posterior
# ----------------------------------------

folders = c("Output_widePrior_new", "Output_narrowPrior")
names(folders) = c("wide", "narrow")

pis_diff = list()
pis_diff[["wide"]] = list()
pis_diff[["narrow"]] = list()
withinZero = list()
for (pat in ptsIDs) {
  for (chan in channels) {
    root = paste0(chan, "__", pat)
    
    fn_post_og = file.path(outputFolder, paste0(root, "__POST.txt"))
    pi_post_og = as.matrix(fread(fn_post_og))[, "probdiff"]
    
    for (pp in names(folders)) {
      fn_post_two = gsub("stan_sampler", folders[pp], fn_post_og)
      pi_post_two = as.matrix(fread(fn_post_two))[, "probdiff"]
      pis_diff[[pp]][[root]] = pi_post_og - pi_post_two
      HDI = hdi(pis_diff[[pp]][[root]], ci=0.99)
      withinZero[[root]] = (0.0>HDI$CI_high) | (0.0<HDI$CI_low)
    }
  }
}

message ("Number of distributions where 0.0 does not fall within the 95% HDI: ", sum ( unlist(withinZero) ))

ks_pValues = vector("numeric", length=length(pis_diff[["wide"]]))
names(ks_pValues) = names(pis_diff[["wide"]])

withinZeroDirect = list()
for (root in names(pis_diff[["wide"]])) {
  diff = pis_diff[["wide"]][[root]] - pis_diff[["narrow"]][[root]]
  HDI = hdi(diff, ci=0.99)
  withinZeroDirect[[root]] = (0.0<HDI$CI_low) | (0.0>HDI$CI_high)
    
  ks_pValues[root] = ks.test(pis_diff[["wide"]][[root]], pis_diff[["narrow"]][[root]])$p.value
  
  densWide = density(pis_diff[["wide"]][[root]])
  densNarrow = density(pis_diff[["narrow"]][[root]])
  plot(NA, 
       xlab=bquote(pi), xlim = range(c(densWide$x, densNarrow$x)),
       ylab="Density", ylim = range(c(densWide$y, densNarrow$y)), 
       main=root)
  lines(densWide, col="blue", lwd=2)
  lines(densNarrow, col="red", lwd=2)
  text(x=par("usr")[2] * 0.9, y=par("usr")[4] * 0.9, 
        adj = c(1.0,0),
       labels=paste("p value:", round(ks_pValues[root], 3)))
  legend("topleft", legend=c("Wide", "Narrow"), col=c("blue", "red"), lty=1)

}


tiff(file.path("..", "Figures",  "TIFF", "prior_comparison.tff"),
     width=2250, height=1200, units="px", res=300, pointsize=13, 
     type="cairo")
{
  op = par( mfrow=c(2,1) )
  
  colVec = c(alphaCol(0.01, name="lightseagreen"), 
             alphaCol(0.01, name="deeppink3"), 
             alphaCol(0.01, name="darkorange"))
  
  myAt = c(1:44)[ 1:44 %% 4 !=0 ]
  myTicks = 4* 1:11 - 2
  
  opp = par(mar=c(1,4,3,2), xpd=TRUE)
  {
    plot(NA, xlim=c(1,43), ylim=c(-0.2,0.2), 
         xaxt="n", 
         ylab="", xlab="", 
         main="")
    for(i in seq(0,10, by=2)) rect( xleft=i*4, xright=(i+1)*4, ybottom=par("usr")[3], ytop=par("usr")[4], border=NA, col=alphaBlack(0.1))
    lines(c(0,44), c(0,0), lty=3, col="black", lwd=2)
    
    stripchart(pis_diff[["wide"]], vertical=TRUE, pch=20, cex=0.5, 
               method="jitter", 
               col=rep(colVec, length(pis_diff)), 
               at=myAt, xaxt='n',
               add=TRUE)
    text(46, -0.01, labels="Wide prior", srt=-90, cex=1.0)
    # axis(side=1, at=myTicks, labels=pts, las=1, cex=1.0)
    exact_xlim = par("usr")[1:2]
    exact_ylim = par("usr")[3:4]
    
    text(-1, 0.19, pos=4, 
         labels="Channels")
    text(c(0,8,13), rep(0.15,3), 
         labels=channels, pos=4)
    points(c(0,8,13), rep(0.15,3), 
           pch=20, col=c(alphaCol(0.7, name="lightseagreen"), 
                         alphaCol(0.7, name="deeppink3"), 
                         alphaCol(0.7, name="darkorange")))
    # legend(exact_xlim[2] + 2.0, exact_ylim[2], 
    #        legend=channels, 
    #        bg="transparent", 
    #        box.col=alphaBlack(0.0),
    #        title = "Channel",
    #        pch=20,  cex=1.0, 
    #        xpd=TRUE,
    #        col=c(alphaCol(0.7, name="lightseagreen"), 
    #              alphaCol(0.7, name="deeppink3"), 
    #              alphaCol(0.7, name="darkorange")) )
  }
  par(opp)
  
  opp = par(mar=c(4,4,0,2), xpd=TRUE)
  {
    plot(NA, xlim=c(1,43), ylim=c(-0.2,0.2), xaxt="n", 
         ylab="", xlab="Patient", cex.lab=1.0,
         main="")
    for(i in seq(0,10, by=2)) rect( xleft=i*4, xright=(i+1)*4, ybottom=par("usr")[3], ytop=par("usr")[4], border=NA, col=alphaBlack(0.1))
    lines(c(0,44), c(0,0), lty=3, col="black", lwd=2)
    
    stripchart(pis_diff[["narrow"]], vertical=TRUE, pch=20, cex=0.5, 
               method="jitter",
               col=rep(colVec, length(pis_diff)), 
               at=myAt, xaxt="n", add=TRUE)
    axis(side=1, at=myTicks, labels=ptsIDs)
    text(46, -0.01, labels="Narrow prior", srt=-90, cex=1.0)
  }
  par(opp)
  mtext("Difference distribution", 
        side=2, line=-1, outer=TRUE, adj=0.6, padj=0.5, cex=1.0)
  par(op)
}
dev.off()

png(file.path("..", "Figures",  "PNG", "prior_comparison.png"),
     width=2250, height=1900, units="px", res=300, pointsize=13, 
     type="cairo")
{
  op = par( mfrow=c(2,1) )
  
  colVec = c(alphaCol(0.01, name="lightseagreen"), 
             alphaCol(0.01, name="deeppink3"), 
             alphaCol(0.01, name="darkorange"))
  
  myAt = c(1:44)[ 1:44 %% 4 !=0 ]
  myTicks = 4* 1:11 - 2
  
  opp = par(mar=c(1,4,3,2), xpd=TRUE)
  {
    plot(NA, xlim=c(1,43), ylim=c(-0.2,0.2), 
         xaxt="n", 
         ylab="", xlab="", 
         main="")
    for(i in seq(0,10, by=2)) rect( xleft=i*4, xright=(i+1)*4, ybottom=par("usr")[3], ytop=par("usr")[4], border=NA, col=alphaBlack(0.1))
    lines(c(0,44), c(0,0), lty=3, col="black", lwd=2)
    
    stripchart(pis_diff[["wide"]], vertical=TRUE, pch=20, cex=0.5, 
               method="jitter", 
               col=rep(colVec, length(pis_diff)), 
               at=myAt, xaxt='n',
               add=TRUE)
    text(46, -0.01, labels="Wide prior", srt=-90, cex=1.0)
    # axis(side=1, at=myTicks, labels=pts, las=1, cex=1.0)
    exact_xlim = par("usr")[1:2]
    exact_ylim = par("usr")[3:4]
    
    text(-1, 0.19, pos=4, 
         labels="Channels")
    text(c(0,8,13), rep(0.15,3), 
         labels=channels, pos=4)
    points(c(0,8,13), rep(0.15,3), 
           pch=20, col=c(alphaCol(0.7, name="lightseagreen"), 
                         alphaCol(0.7, name="deeppink3"), 
                         alphaCol(0.7, name="darkorange")))
    # legend(exact_xlim[2] + 2.0, exact_ylim[2], 
    #        legend=channels, 
    #        bg="transparent", 
    #        box.col=alphaBlack(0.0),
    #        title = "Channel",
    #        pch=20,  cex=1.0, 
    #        xpd=TRUE,
    #        col=c(alphaCol(0.7, name="lightseagreen"), 
    #              alphaCol(0.7, name="deeppink3"), 
    #              alphaCol(0.7, name="darkorange")) )
  }
  par(opp)
  
  opp = par(mar=c(4,4,0,2), xpd=TRUE)
  {
    plot(NA, xlim=c(1,43), ylim=c(-0.2,0.2), xaxt="n", 
         ylab="", xlab="Patient", cex.lab=1.0,
         main="")
    for(i in seq(0,10, by=2)) rect( xleft=i*4, xright=(i+1)*4, ybottom=par("usr")[3], ytop=par("usr")[4], border=NA, col=alphaBlack(0.1))
    lines(c(0,44), c(0,0), lty=3, col="black", lwd=2)
    
    stripchart(pis_diff[["narrow"]], vertical=TRUE, pch=20, cex=0.5, 
               method="jitter",
               col=rep(colVec, length(pis_diff)), 
               at=myAt, xaxt="n", add=TRUE)
    axis(side=1, at=myTicks, labels=ptsIDs)
    text(46, -0.01, labels="Narrow prior", srt=-90, cex=1.0)
  }
  par(opp)
  mtext("Difference distribution", 
        side=2, line=-1, outer=TRUE, adj=0.6, padj=0.5, cex=1.0)
  par(op)
}
dev.off()

# --- ALL POSTERIOR PREDICTIVIONS - supplementary

tiff(file.path("..", "Figures",  "TIFF", "allPost_classifs_1.tiff"), 
     width=2250, height=2625, units="px", res=300, pointsize=13,
     type="cairo")
{
  op = par(cex.lab=1.3, cex.axis=1.3, mfrow=c(4,3), mar=c(4,4,2,0.5), xpd=FALSE)
  
  for (pat in ptsIDs[1:4]) {
    for (chan in channels) {
      root = paste0(chan, "__", pat)
      fn_post = file.path(outputFolder, paste0(root, "__POST.txt"))
      fn_postpred = file.path(outputFolder, paste0(root, "__POSTPRED.txt"))
      fn_classif = file.path(outputFolder, paste0(root, "__CLASSIF.txt"))
      
      class_mats = as.matrix( fread(fn_classif) )
      class = colMeans( class_mats ) 
      
      postpred = as.matrix( fread(fn_postpred) )
      post = as.matrix( fread(fn_post) )
      pi_mean = round(mean( post[,"probdiff"] ), 3)
      
      dataMats = getData_mats(data, channels=c(mitochan, chan), 
                              pts=pat, ctrlID=ctrlIDs)
      
      classif_plot(dataMats, 
                   classifs = class,
                   postpred = postpred, 
                   colAlpha=0.2, 
                   cex=0.5, 
                   xlab=paste0("log(", mitochan, ")"), 
                   ylab=paste0("log(", chan, ")"),
                   main=bquote("E("*pi*"|Y)="*.(pi_mean)))
      if (chan == "NDUFB8") {
        text(par("usr")[1], par("usr")[4], labels=pat, adj=c(-0.2,1.3), cex=1.5)
      }
    }
  }
  par(op)
}
dev.off()

tiff(file.path("..", "Figures",  "TIFF", "allPost_classifs_2.tiff"), 
     width=2250, height=2625, units="px", res=300, pointsize=13, 
     type="cairo")
{
  op = par(cex.lab=1.3, cex.axis=1.3, mfrow=c(4,3), mar=c(4,4,2,0.5), xpd=FALSE)
  
  for( pat in ptsIDs[5:8] ){
    for( chan in channels ){
      root = paste0(chan, "__", pat)
      fn_post = file.path(outputFolder, paste0(root, "__POST.txt"))
      fn_postpred = file.path(outputFolder, paste0(root, "__POSTPRED.txt"))
      fn_classif = file.path(outputFolder, paste0(root, "__CLASSIF.txt"))
      
      class_mats = as.matrix( fread(fn_classif) )
      class = colMeans( class_mats ) 
      
      postpred = as.matrix( fread(fn_postpred) )
      post = as.matrix( fread(fn_post) )
      pi_mean = round(mean( post[,"probdiff"] ), 3)
      
      dataMats = getData_mats(data, channels=c(mitochan, chan), 
                              pts=pat, ctrlID=ctrlIDs)
      
      classif_plot(dataMats, 
                   classifs = class,
                   postpred = postpred, 
                   colAlpha=0.2, 
                   cex=0.5, 
                   xlab=paste0("log(", mitochan, ")"), 
                   ylab=paste0("log(", chan, ")"),
                   main=bquote("E("*pi*"|Y)="*.(pi_mean)))
      if( chan == "NDUFB8" ){
        text(par("usr")[1], par("usr")[4], labels=pat, adj=c(-0.2,1.3), cex=1.5)
      }
    }
  }
  par(op)
}
dev.off()

tiff(file.path("..", "Figures",  "TIFF", "allPost_classifs_3.tiff"),
     width=2250, height=2625 *2/3, units="px", res=300, pointsize=13, 
     type="cairo")
{
  op = par(cex.lab=1.3, cex.axis=1.3, mfrow=c(3,3), mar=c(4,4,2,0.5), xpd=FALSE)
  
  for( pat in ptsIDs[9:11] ){
    for( chan in channels ){
      root = paste0(chan, "__", pat)
      fn_post = file.path(outputFolder, paste0(root, "__POST.txt"))
      fn_postpred = file.path(outputFolder, paste0(root, "__POSTPRED.txt"))
      fn_classif = file.path(outputFolder, paste0(root, "__CLASSIF.txt"))
      
      class_mats = as.matrix( fread(fn_classif) )
      class = colMeans( class_mats ) 
      
      postpred = as.matrix( fread(fn_postpred) )
      post = as.matrix( fread(fn_post) )
      pi_mean = round(mean( post[,"probdiff"] ), 3)
      
      dataMats = getData_mats(data, channels=c(mitochan, chan), 
                              pts=pat, ctrlID=ctrlIDs)
      
      classif_plot(dataMats, 
                   classifs = class,
                   postpred = postpred, 
                   colAlpha=0.2, 
                   cex=0.5, 
                   xlab=paste0("log(", mitochan, ")"), 
                   ylab=paste0("log(", chan, ")"),
                   main=bquote("E("*pi*"|Y)="*.(pi_mean)))
      if( chan == "NDUFB8" ){
        text(par("usr")[1], par("usr")[4], labels=pat, adj=c(-0.2,1.3), cex=1.5)
      }
    }
  }
  par(op)
}
dev.off()

png(file.path("..", "Figures",  "PNG", "allPost_classifs_1.png"),
     width=2250, height=2625, units="px", res=300, pointsize=13,
     type="cairo")
{
  op = par(cex.lab=1.3, cex.axis=1.3, mfrow=c(4,3), mar=c(4,4,2,0.5), xpd=FALSE)
  
  for( pat in ptsIDs[1:4] ){
    for( chan in channels ){
      root = paste0(chan, "__", pat)
      fn_post = file.path(outputFolder, paste0(root, "__POST.txt"))
      fn_postpred = file.path(outputFolder, paste0(root, "__POSTPRED.txt"))
      fn_classif = file.path(outputFolder, paste0(root, "__CLASSIF.txt"))
      
      class_mats = as.matrix( fread(fn_classif) )
      class = colMeans( class_mats ) 
      
      postpred = as.matrix( fread(fn_postpred) )
      post = as.matrix( fread(fn_post) )
      pi_mean = round(mean( post[,"probdiff"] ), 3)
      
      dataMats = getData_mats(data, channels=c(mitochan, chan), 
                              pts=pat, ctrlID=ctrlIDs)
      
      classif_plot(dataMats, 
                   classifs = class,
                   postpred = postpred, 
                   colAlpha=0.2, 
                   cex=0.5, 
                   xlab=paste0("log(", mitochan, ")"), 
                   ylab=paste0("log(", chan, ")"),
                   main=bquote("E("*pi*"|Y)="*.(pi_mean)))
      if( chan == "NDUFB8" ){
        text(par("usr")[1], par("usr")[4], labels=pat, adj=c(-0.2,1.3), cex=1.5)
      }
    }
  }
  par(op)
}
dev.off()

png(file.path("..", "Figures",  "PNG", "allPost_classifs_2.png"), 
     width=2250, height=2625, units="px", res=300, pointsize=13, 
     type="cairo")
{
  op = par(cex.lab=1.3, cex.axis=1.3, mfrow=c(4,3), mar=c(4,4,2,0.5), xpd=FALSE)
  
  for( pat in ptsIDs[5:8] ){
    for( chan in channels ){
      root = paste0(chan, "__", pat)
      fn_post = file.path(outputFolder, paste0(root, "__POST.txt"))
      fn_postpred = file.path(outputFolder, paste0(root, "__POSTPRED.txt"))
      fn_classif = file.path(outputFolder, paste0(root, "__CLASSIF.txt"))
      
      class_mats = as.matrix( fread(fn_classif) )
      class = colMeans( class_mats ) 
      
      postpred = as.matrix( fread(fn_postpred) )
      post = as.matrix( fread(fn_post) )
      pi_mean = round(mean( post[,"probdiff"] ), 3)
      
      dataMats = getData_mats(data, channels=c(mitochan, chan), 
                              pts=pat, ctrlID=ctrlIDs)
      
      classif_plot(dataMats, 
                   classifs = class,
                   postpred = postpred, 
                   colAlpha=0.2, 
                   cex=0.5, 
                   xlab=paste0("log(", mitochan, ")"), 
                   ylab=paste0("log(", chan, ")"),
                   main=bquote("E("*pi*"|Y)="*.(pi_mean)))
      if( chan == "NDUFB8" ){
        text(par("usr")[1], par("usr")[4], labels=pat, adj=c(-0.2,1.3), cex=1.5)
      }
    }
  }
  par(op)
}
dev.off()

png(file.path("..", "Figures",  "PNG", "allPost_classifs_3.png"),
     width=2250, height=2625 *2/3, units="px", res=300, pointsize=13, 
     type="cairo")
{
  op = par(cex.lab=1.3, cex.axis=1.3, mfrow=c(3,3), mar=c(4,4,2,0.5), xpd=FALSE)
  
  for( pat in ptsIDs[9:11] ){
    for( chan in channels ){
      root = paste0(chan, "__", pat)
      fn_post = file.path(outputFolder, paste0(root, "__POST.txt"))
      fn_postpred = file.path(outputFolder, paste0(root, "__POSTPRED.txt"))
      fn_classif = file.path(outputFolder, paste0(root, "__CLASSIF.txt"))
      
      class_mats = as.matrix( fread(fn_classif) )
      class = colMeans( class_mats ) 
      
      postpred = as.matrix( fread(fn_postpred) )
      post = as.matrix( fread(fn_post) )
      pi_mean = round(mean( post[,"probdiff"] ), 3)
      
      dataMats = getData_mats(data, channels=c(mitochan, chan), 
                              pts=pat, ctrlID=ctrlIDs)
      
      classif_plot(dataMats, 
                   classifs = class,
                   postpred = postpred, 
                   colAlpha=0.2, 
                   cex=0.5, 
                   xlab=paste0("log(", mitochan, ")"), 
                   ylab=paste0("log(", chan, ")"),
                   main=bquote("E("*pi*"|Y)="*.(pi_mean)))
      if( chan == "NDUFB8" ){
        text(par("usr")[1], par("usr")[4], labels=pat, adj=c(-0.2,1.3), cex=1.5)
      }
    }
  }
  par(op)
}
dev.off()

# --------------------------------
# --- pi posterior - supplementary
# --------------------------------
bayes_post = list()
for (pat in ptsIDs) {
  for (chan in channels) {
    root = paste0(chan, "__", pat)
    fn_post = file.path(outputFolder, paste0(root, "__POST.txt"))
    probdiff_post = as.data.frame( fread( fn_post ))[,"probdiff"]
    bayes_post[[root]] = probdiff_post
  }
}

tiff(file.path("..", "Figures",  "TIFF", "piPost.tiff"),
     width=2250, height=1200, units="px", res=300, pointsize=13, 
     type="cairo")
{
  op = par(cex.lab=1.2, cex.axis=1.2, mar=c(4,4,1,0.5), xpd=FALSE)
  
  colVec = c(alphaCol(0.01, name="lightseagreen"), 
             alphaCol(0.01, name="deeppink3"), 
             alphaCol(0.01, name="darkorange"))
  
  myAt = c(1:44)[ 1:44 %% 4 !=0 ]
  myTicks = 4* 1:11 - 2
  plot(NA, xlim=c(1,43), ylim=c(0.0,1), xaxt="n", 
       ylab="Deficiency proportion", xlab="Patient")
  
  for(i in seq(0,10, by=2)) rect( xleft=i*4, xright=(i+1)*4, ybottom=par("usr")[3], ytop=par("usr")[4], border=NA, col=alphaBlack(0.1))
  
  stripchart(bayes_post, vertical=TRUE, pch=20, cex=0.5, 
             method="jitter", main="", ylim=c(0,1), xlim=c(1,43),
             xlab="Patient", ylab="Deficient Proportion", 
             col=rep(colVec, length(bayes_post)), 
             at=myAt, xaxt="n", add=TRUE)
  axis(side=1, at=myTicks, labels=ptsIDs)
  legend("topleft", 
         bg = "transparent",
         legend=channels, 
         box.col=alphaBlack(0.0),
         text.width = 10,
         pch=20, 
         col=c(alphaCol(0.7, name="lightseagreen"), 
               alphaCol(0.7, name="deeppink3"), 
               alphaCol(0.7, name="darkorange")))
  text(x=45, y=1.0, labels=c("Bayesian distribution"), pos=4)
  rect(xleft=45, xright=55.7, 
       ybottom=0.82, ytop=1.025, 
       lwd=1, col=NA, border="black")
  par(op)
}
dev.off()

png(file.path("..", "Figures",  "PNG", "piPost.png"), 
     width=2250, height=1200, units="px", res=300, pointsize=13, 
     type="cairo")
{
  op = par(cex.lab=1.2, cex.axis=1.2, mar=c(4,4,1,0.5), xpd=FALSE)
  
  colVec = c(alphaCol(0.01, name="lightseagreen"), 
             alphaCol(0.01, name="deeppink3"), 
             alphaCol(0.01, name="darkorange"))
  
  myAt = c(1:44)[ 1:44 %% 4 !=0 ]
  myTicks = 4* 1:11 - 2
  plot(NA, xlim=c(1,43), ylim=c(0.0,1), xaxt="n", 
       ylab="Deficiency proportion", xlab="Patient")
  
  for(i in seq(0,10, by=2)) rect( xleft=i*4, xright=(i+1)*4, ybottom=par("usr")[3], ytop=par("usr")[4], border=NA, col=alphaBlack(0.1))
  
  stripchart(bayes_post, vertical=TRUE, pch=20, cex=0.5, 
             method="jitter", main="", ylim=c(0,1), xlim=c(1,43),
             xlab="Patient", ylab="Deficient Proportion", 
             col=rep(colVec, length(bayes_post)), 
             at=myAt, xaxt="n", add=TRUE)
  axis(side=1, at=myTicks, labels=ptsIDs)
  legend("topleft", 
         bg = "transparent",
         legend=channels, 
         box.col=alphaBlack(0.0),
         text.width = 10,
         pch=20, 
         col=c(alphaCol(0.7, name="lightseagreen"), 
               alphaCol(0.7, name="deeppink3"), 
               alphaCol(0.7, name="darkorange")))
  text(x=45, y=1.0, labels=c("Bayesian distribution"), pos=4)
  rect(xleft=45, xright=55.7, 
       ybottom=0.82, ytop=1.025, 
       lwd=1, col=NA, border="black")
  par(op)
}
dev.off()

# ----------------------------------------------------
# --- wider priors and narrower priors - supplementary
# ----------------------------------------------------

tiff(file.path("..", "Figures",  "TIFF", "narrow_wide_priors.tiff"), 
    width=2250/2, height=600*3/2, units="px", res=150, pointsize=13, 
    type="cairo")
{
  op = par(mfrow=c(3,2), mar=c(4,4,1,0.5))
  
  pat = "P09"
  chan = "NDUFB8"
  filename = paste0(chan, "__", pat, "__PRIOR.txt")
  fn_prior_n = file.path(gsub(basename(outputFolder), "Output_narrowPrior_new", outputFolder), filename)
  fn_prior_w = file.path(gsub(basename(outputFolder), "Output_widePrior_new", outputFolder), filename)
  
  prior_n = as.matrix( fread(fn_prior_n) )
  prior_w = as.matrix( fread(fn_prior_w) )
  
  xlabels = list(bquote(mu[m]), bquote(tau[m]), bquote(mu[c]), bquote(tau[c]), bquote(tau))
  names(xlabels) = c("mu_m", "tau_m", "mu_c", "tau_c", "tau_norm")
  
  for (pp in names(xlabels)) {
    dens_n = density( prior_n[,pp])
    dens_w = density( prior_w[,pp])
    
    plot(NA, ylab="Density", xlab=xlabels[[pp]],
         xlim=range(c(dens_n$x, dens_w$x)), 
         ylim=c(0, max(c(dens_n$y, dens_w$y))) )
    lines( dens_n, col="brown4" )
    lines( dens_w, col="darkolivegreen" )
  }
  
  plot(NA, xlim=c(0,1), ylim=c(0,1), axes=FALSE, xlab="", ylab="")
  legend("topleft", 
         legend=c("Narrow", "Wide"), lty=1, col=c("brown4", "darkolivegreen"), 
         box.col=alphaBlack(0.0))
  par(op)
  
}
dev.off()

png(file.path("..", "Figures",  "PNG", "narrow_wide_priors.png"), 
    width=2250, height=600*3, units="px", res=300, pointsize=13, 
    type="cairo")
{
  op = par(mfrow=c(3,2), mar=c(4,4,1,0.5))
  
  pat = "P09"
  chan = "NDUFB8"
  filename = paste0(chan, "__", pat, "__PRIOR.txt")
  fn_prior_n = file.path(gsub(basename(outputFolder), "Output_narrowPrior_new", outputFolder), filename)
  fn_prior_w = file.path(gsub(basename(outputFolder), "Output_widePrior_new", outputFolder), filename)
  
  prior_n = as.matrix( fread(fn_prior_n) )
  prior_w = as.matrix( fread(fn_prior_w) )
  
  xlabels = list(bquote(mu[m]), bquote(tau[m]), bquote(mu[c]), bquote(tau[c]), bquote(tau))
  names(xlabels) = c("mu_m", "tau_m", "mu_c", "tau_c", "tau_norm")
  
  for (pp in names(xlabels)) {
    dens_n = density( prior_n[,pp])
    dens_w = density( prior_w[,pp])
    
    plot(NA, ylab="Density", xlab=xlabels[[pp]],
         xlim=range(c(dens_n$x, dens_w$x)), 
         ylim=c(0, max(c(dens_n$y, dens_w$y))) )
    lines( dens_n, col="brown4" )
    lines( dens_w, col="darkolivegreen" )
  }
  
  plot(NA, xlim=c(0,1), ylim=c(0,1), axes=FALSE, xlab="", ylab="")
  legend("topleft", 
         legend=c("Narrow", "Wide"), lty=1, col=c("brown4", "darkolivegreen"), 
         box.col=alphaBlack(0.0))
  par(op)
  
}
dev.off()

pdf(file.path("..", "Figures",  "PDF", "pi_post_comparison.pdf"), 
    width=12, height=6)
{
  op = par(mfrow=c(1,3))
  for(pat in ptsIDs){
    for(chan in channels){
      file = paste0(chan, "__", pat, "__POST.txt")
      
      fn_post_n = file.path(folder, "Output_narrowPrior_new", file)
      fn_post_w = file.path(folder, "Output_widePrior_new", file)
      
      post_n = as.matrix( fread(fn_post_n) )
      post_w = as.matrix( fread(fn_post_w) )
      
      dens_n = density( post_n[,"probdiff"] )
      dens_w = density( post_w[,"probdiff"] )
      
      plot(NA, 
           xlim=range(c(dens_n$x, dens_w$x)), ylim=c(0, max(c(dens_n$y, dens_w$y))),
           xlab=bquote(pi), ylab="Density")
      lines(dens_n, col="brown3", lwd=2)
      lines(dens_w, col="darkolivegreen", lwd=2)
    }
  }
  par(op)
}
dev.off()

# ------------------------------
# --- bayesian classif vs manual
# ------------------------------
mdat = as.data.frame( fread( file.path("..", "..", "dat_with_class_prepped.txt") ) )

png(file.path("..", "Figures",  "PNG", "diff_in_class.png"), 
    width=2250, height=1200, units="px", res=300, pointsize=13, 
    type="cairo")
{
  pat = "P09"
  op = par(cex.lab=1.3, cex.axis=1.3, mfrow=c(1,3), mar=c(4,4,1,0.5))
  
  for( chan in channels ){
    man_class = mdat[mdat$Channel==chan & mdat$sampleID==pat, "classCONOR"]
    
    bayes_class_mat = as.matrix( fread( file.path(outputFolder, paste0(chan, "__", pat, "__CLASSIF.txt"))))
    bayes_class = colMeans( bayes_class_mat )
    
    class_diff = man_class != as.numeric(bayes_class>0.5)
    
    bayes_postpred = as.data.frame( fread(  file.path(outputFolder, paste0(chan, "__", pat, "__POSTPRED.txt") ) ) )
    
    xlm = range( data[data$sampleID%in% c(ctrlIDs, pat) & data$Channel==mitochan, "Value"] ) + c(0,0.5)
    ylm = range( data[data$sampleID%in% c(ctrlIDs, pat) & data$Channel!=mitochan, "Value"] ) + c(0,0.5)
    
    pi_est = round( colMeans(  fread(  file.path(outputFolder, paste0(chan, "__", pat, "__POST.txt")))[,"probdiff"]),3)
    
    plot(NA, xlim=xlm, ylim=ylm, 
         xlab=paste0("log(", mitochan,")"), ylab=paste0("log(", chan, ")"),
         main=bquote("E("*pi*"|Y)="*.(pi_est)))
    
    points( data[data$sampleID%in% ctrlIDs & data$Channel==mitochan, "Value"],
            data[data$sampleID%in% ctrlIDs & data$Channel==chan, "Value"],
            pch=20, col=alphaBlack(0.05), cex=0.3)
    
    xpat = data[data$sampleID==pat & data$Channel==mitochan, "Value"]
    ypat = data[data$sampleID==pat & data$Channel==chan, "Value"]
    points( xpat, ypat,
            pch=20, col=classcols(bayes_class, alphaLevel=0.5), cex=0.5)
    
    # lines(bayes_postpred$mitochan, bayes_postpred$medNorm, 
    #       col=alphaGreen(1.0), lty=1, lwd=2)
    # lines(bayes_postpred$mitochan, bayes_postpred$lwrNorm, 
    #       col=alphaGreen(1.0), lty=2, lwd=2)
    # lines(bayes_postpred$mitochan, bayes_postpred$uprNorm, 
    #       col=alphaGreen(1.0), lty=2, lwd=2)
    
    points( xpat[class_diff], ypat[class_diff],
            pch=18, col="orange", cex=1.0)
    
    if( chan == channels[1] )
      legend("topleft", 
             col=c(alphaBlack(0.7), alphaBlue(0.7), alphaRed(0.7)),
             pch=20, 
             legend=c("Control", "Healthy", "Deficient"),
             box.col=alphaBlack(0.0),
             bg="transparent")
  }
  par(op)
}
dev.off()







