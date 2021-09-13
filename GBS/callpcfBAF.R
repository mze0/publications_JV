callpcfBAF <- function(BAFph,gammaSC,gammaMC,plateau){
  source("fastPCF.R")
  print("PCF segmentation is applying...")
  SegGenomes <- NULL
  for(i in 4:ncol(BAFph)){

    if(sum(grep("E*_Bl*",colnames(BAFph)[i]))==1 | sum(grep("_MC",colnames(BAFph)[i]))==1){gamma=gammaSC}else{gamma=gammaMC}
    SegGenome <- NULL
    for(chr in unique(BAFph[,"Chr"])){
      logRChr <- BAFph[as.character(BAFph$Chr)==chr,i]

      #if(sum(is.na(logRChr))>=1){warning(paste("There are",sum(is.na(logRChr),"missing values in logR-values of Chr.",chr))}

      while(sum(is.na(logRChr))>=1){logRChr[which(is.na(logRChr))]<-logRChr[which(is.na(logRChr))-1]}

      sdev <- getMad(logRChr,k=plateau)
      res <- selectFastPcf(logRChr,3,gamma*sdev,T)
      SegChr <- res$yhat

      SegGenome <- c(SegGenome,SegChr)
    }#end chr loop
    

    SegGenomes<-cbind(SegGenomes,SegGenome)
    print(paste(colnames(BAFph)[i],"==> gamma",gamma, "is applied"))

  }#end file loop

  SegLogRs <- cbind(BAFph[,c("Name","Chr","Position")],SegGenomes)
  colnames(SegLogRs)<- colnames(BAFph)

  SegLogRs
}#end function
