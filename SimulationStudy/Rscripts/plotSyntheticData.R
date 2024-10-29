
mitochan = "VDAC"
channels = c("MTCO1", "NDUFB8", "CYB")

trueParameters = as.data.frame(fread(file.path("..", "Data", "trueParameters_synData03.txt")))
dataFolder = file.path("..", "Data", "synData03")
outputFolder = file.path("..", "OutputMaxESS", "stan_sampler_synData03")

### ---------------------
### plot synthetic data
### ---------------------

alphaCol = function (col, alpha) {
  if (is.numeric(col)) {
    rgbVals = col2rgb(palette()[i])
  } 
  if (is.character(col)) {
    rgbVals = col2rgb(col)
  }
  print ( rgbVals )
  return( rgb(red=rgbVals[1] / 255, green=rgbVals[2] / 255, blue=rgbVals[3] / 255, alpha=alpha) )
}

fps = list.files(dataFolder, full.names = TRUE, pattern="__synDATA")

pdf(file.path("..", "PDF", "data", paste0(basename(dataFolder), ".pdf")),
    width=9, height=5)
for (fp in fps) {
  dat = fread(fp)
  fn = strsplit(basename(fp), split="__")
  chan = fn[[1]][1]
  patID = fn[[1]][2]
  
  thetaTrue = trueParameters[trueParameters$Channel==chan & trueParameters$sampleID==patID, ]
  
  posterior_fn = file.path(outputFolder,  paste0(chan, "__", patID, "__POST.txt"))
  post = as.data.frame(fread(posterior_fn))
  
  op = par(mfrow=c(2,3), mar=c(4,4,4,1))
  {
    xLim_true = range( dat$trueValue[(dat$sbjType=="control" | dat$sampleID==patID) & dat$Channel==mitochan] ) + c(-1,1)
    yLim_true = range( dat$trueValue[(dat$sbjType=="control" | dat$sampleID==patID) & dat$Channel==chan] ) + c(-1,1)
    
    plot(NA,
         xlim=xLim_true, ylim=yLim_true, 
         xlab="log(VDAC)", 
         ylab=paste0("log(", chan, ")"), 
         main="Observed Data")
    
    points(dat$trueValue[dat$sbjType=="control" & dat$Channel==mitochan], 
           dat$trueValue[dat$sbjType=="control" & dat$Channel==chan], 
           pch=20, col=alphaBlack(0.1))
    points(dat$trueValue[dat$sampleID==patID & dat$Channel==mitochan], 
           dat$trueValue[dat$sampleID==patID & dat$Channel==chan],
           pch=20, col=alphaBlue(0.2))
    
    xLim_syn = range( dat$Value[dat$Channel==mitochan] ) + c(-1,1)
    yLim_syn= range( dat$Value[dat$Channel==chan] ) + c(-1,1)
    
    
    plot(NA,
         xlim = xLim_syn, ylim=yLim_syn, 
         xlab="log(VDAC)", 
         ylab=paste0("log(", chan, ")"), 
         main="Synthetic Data")
    points(dat$Value[dat$sbjType=="control" & dat$Channel==mitochan], 
           dat$Value[dat$sbjType=="control" & dat$Channel==chan], 
           pch=20, col=alphaBlack(0.05))
    points(dat$Value[dat$sampleID==patID & dat$Channel==mitochan], 
           dat$Value[dat$sampleID==patID & dat$Channel==chan],
          pch=20, col=alphaBlue(0.5))
    
    defObsID = dat$obsID[dat$deficientStatus==1]
    points(dat$Value[dat$sampleID==patID & dat$Channel==mitochan & dat$obsID%in%defObsID], 
           dat$Value[dat$sampleID==patID & dat$Channel==chan & dat$obsID%in%defObsID], 
           pch=20, cex=0.33, col="yellow")
    
    xLim_syn = range( dat$Value[(dat$sbjType=="control" | dat$sampleID==patID) & dat$Channel==mitochan] )
    yLim_syn = range( dat$Value[(dat$sbjType=="control" | dat$sampleID==patID) & dat$Channel==chan] )
    
    plot(NA,
         xlim = c(0,1), ylim=c(0,1), 
         xlab="", ylab="", main="", 
         axes=FALSE)

    
    defObsID = dat$obsID[dat$deficientStatus==1]
    
    points(dat$Value[dat$sampleID==patID & dat$Channel==mitochan & dat$obsID%in%defObsID], 
           dat$Value[dat$sampleID==patID & dat$Channel==chan & dat$obsID%in%defObsID], 
           pch=20, cex=0.33, col="yellow")
    
    title(main=patID, outer=TRUE, line=-1)
    
    ##################
    # parameter values
    ##################
    
    mPredDens = density(post$m_pred)
    plot (mPredDens, 
          xlab="Slope", xlim=range(c(thetaTrue[c("m[1]", "m[2]", "m[3]", "m[4]", "m[5]")], mPredDens$x)),
          ylab="Density", ylim=c(0, max(mPredDens$y)),
          main="",
          lty=2, col=alphaGreen(0.7), 
          )
    abline(v=thetaTrue[c("m[1]", "m[2]", "m[3]", "m[4]")], col=alphaGreen(0.7))
    abline(v=thetaTrue["m[5]"], col=alphaGreen(1.0), lwd=3)
    
    cPredDens = density(post$c_pred)
    plot (cPredDens, 
          xlab="Intercept", xlim=range(c(thetaTrue[c("c[1]", "c[2]", "c[3]", "c[4]", "c[5]")], cPredDens$x)),
          ylab="Density", ylim=c(0, max(cPredDens$y)),
          main="",
          lty=2, col=alphaGreen(0.7), 
    )
    abline(v=thetaTrue[c("c[1]", "c[2]", "c[3]", "c[4]")], col=alphaGreen(0.7))
    abline(v=thetaTrue["c[5]"], col=alphaGreen(1.0), lwd=3)
    
    plot(NA, 
         xlim=c(0,1), xlab="", 
         ylim=c(0,1), ylab="", 
         axes=FALSE)
    legend("topleft", 
           legend=c(bquote(tau == .(as.numeric(thetaTrue["tau_norm"]))), 
                    bquote(pi == .(as.numeric(thetaTrue["probdiff"]))),
                    bquote(mu[m] == .(as.numeric(thetaTrue["mu_m"]))),
                    bquote(tau[m] == .(as.numeric(thetaTrue["tau_m"]))),
                    bquote(mu[c] == .(as.numeric(thetaTrue["mu_c"]))),
                    bquote(tau[c] == .(as.numeric(thetaTrue["tau_c"]))) 
                    ),
           pch=NA, 
           bty='n')
    
  }
  par(op)
}
dev.off()

### ---------------------
### plot synthetic data - with validation
### ---------------------
pdf(file.path("..", "PDF", "data", paste0(basename(dataFolder), ".pdf")),
    width=9, height=5)
for (fp in fps) {
  dat = fread(fp)
  fn = strsplit(basename(fp), split="__")
  chan = fn[[1]][1]
  patID = fn[[1]][2]
  
  thetaTrue = trueParameters[trueParameters$Channel==chan & trueParameters$sampleID==patID, ]
  
  posterior_fn = file.path(outputFolder,  paste0(chan, "__", patID, "__POST.txt"))
  post = as.data.frame(fread(posterior_fn))
  
  op = par(mfrow=c(2,3), mar=c(4,4,4,1))
  {
    xLim_true = range( dat$trueValue[(dat$sbjType=="control" | dat$sampleID==patID) & dat$Channel==mitochan] ) + c(-1,1)
    yLim_true = range( dat$trueValue[(dat$sbjType=="control" | dat$sampleID==patID) & dat$Channel==chan] ) + c(-1,1)
    
    plot(NA,
         xlim=xLim_true, ylim=yLim_true, 
         xlab="log(VDAC)", 
         ylab=paste0("log(", chan, ")"), 
         main="Observed Data")
    
    points(dat$trueValue[dat$sbjType=="control" & dat$Channel==mitochan], 
           dat$trueValue[dat$sbjType=="control" & dat$Channel==chan], 
           pch=20, col=alphaBlack(0.1))
    points(dat$trueValue[dat$sampleID==patID & dat$Channel==mitochan], 
           dat$trueValue[dat$sampleID==patID & dat$Channel==chan],
           pch=20, col=alphaBlue(0.2))
    
    xLim_syn = range( dat$Value[dat$Channel==mitochan] ) + c(-1,1)
    yLim_syn= range( dat$Value[dat$Channel==chan] ) + c(-1,1)
    
    
    plot(NA,
         xlim = xLim_syn, ylim=yLim_syn, 
         xlab="log(VDAC)", 
         ylab=paste0("log(", chan, ")"), 
         main="Synthetic Data")
    points(dat$Value[dat$sbjType=="control" & dat$Channel==mitochan], 
           dat$Value[dat$sbjType=="control" & dat$Channel==chan], 
           pch=20, col=alphaBlack(0.05))
    
    defObsID = dat$obsID[dat$deficientStatus==1]
    
    points(dat$Value[dat$sampleID==patID & dat$Channel==mitochan & !(dat$obsID%in%defObsID)], 
           dat$Value[dat$sampleID==patID & dat$Channel==chan & !(dat$obsID%in%defObsID)], 
           pch=20, cex=1.0, col=alphaBlue(0.2))
    points(dat$Value[dat$sampleID==patID & dat$Channel==mitochan & dat$obsID%in%defObsID], 
           dat$Value[dat$sampleID==patID & dat$Channel==chan & dat$obsID%in%defObsID], 
           pch=20, cex=1.0, col=alphaRed(0.2))
    
    validationObsID = unique( dat$obsID[dat$sbjType=="patient" & dat$validationSet==1] )
    
    points(dat$Value[dat$sampleID==patID & dat$Channel==mitochan & (dat$obsID%in%validationObsID)],
           dat$Value[dat$sampleID==patID & dat$Channel==chan & (dat$obsID%in%validationObsID)],
           pch=20, col="green", cex=0.33
           )
    
    plot (NA, 
          xlim=c(0,1), xlab="",
          ylim=c(0,1), ylab="",
          axes=FALSE)
    legend("topleft", 
           legend = c(paste0("nPat: ",  nrow(dat[dat$sbjType=="patient" & dat$Channel==mitochan])), 
                      paste0("nFit: ", nrow(dat[dat$sbjType=="patient" & dat$Channel==mitochan & !(dat$obsID%in%validationObsID)])),
                      paste0("nVal: ", nrow(dat[dat$sbjType=="patient" & dat$Channel==mitochan & (dat$obsID%in%validationObsID)]) ) ),
           pch=NA, 
           bty='n')
    
    ##################
    # parameter values
    ##################
    
    mPredDens = density(post$m_pred)
    plot (mPredDens, 
          xlab="Slope", xlim=range(c(thetaTrue[c("m[1]", "m[2]", "m[3]", "m[4]", "m[5]")], mPredDens$x)),
          ylab="Density", ylim=c(0, max(mPredDens$y)),
          main="",
          lty=2, col=alphaGreen(0.7), 
    )
    abline(v=thetaTrue[c("m[1]", "m[2]", "m[3]", "m[4]")], col=alphaGreen(0.7))
    abline(v=thetaTrue["m[5]"], col=alphaGreen(1.0), lwd=3)
    
    cPredDens = density(post$c_pred)
    plot (cPredDens, 
          xlab="Intercept", xlim=range(c(thetaTrue[c("c[1]", "c[2]", "c[3]", "c[4]", "c[5]")], cPredDens$x)),
          ylab="Density", ylim=c(0, max(cPredDens$y)),
          main="",
          lty=2, col=alphaGreen(0.7), 
    )
    abline(v=thetaTrue[c("c[1]", "c[2]", "c[3]", "c[4]")], col=alphaGreen(0.7))
    abline(v=thetaTrue["c[5]"], col=alphaGreen(1.0), lwd=3)
    
    plot(NA, 
         xlim=c(0,1), xlab="", 
         ylim=c(0,1), ylab="", 
         axes=FALSE)
    legend("topleft", 
           legend=c(bquote(tau == .(as.numeric(thetaTrue["tau_norm"]))), 
                    bquote(pi == .(as.numeric(thetaTrue["probdiff"]))),
                    bquote(mu[m] == .(as.numeric(thetaTrue["mu_m"]))),
                    bquote(tau[m] == .(as.numeric(thetaTrue["tau_m"]))),
                    bquote(mu[c] == .(as.numeric(thetaTrue["mu_c"]))),
                    bquote(tau[c] == .(as.numeric(thetaTrue["tau_c"]))) 
           ),
           pch=NA, 
           bty='n')
    
  }
  par(op)
}
dev.off()

