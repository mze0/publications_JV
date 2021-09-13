meanwindow <- function(logRsRaw,GC,window,Func,Family,ParScore,outPath){
  
Rs<- logRsRaw[,-c(1:3)]

MedRs <- NULL
for(ind in colnames(Rs)){  
  MedR <- NULL
  for(chr in unique(logRsRaw$Chr)){
    
    MedRChr <-matrix(data=NA,sum(logRsRaw$Chr==chr))
    RsChr <- Rs[logRsRaw$Chr==chr,ind]
    for(i in (1+window):(length(RsChr)-window)){	
      MedRChr[i] <- eval(parse(text=paste(Func,"(na.omit(RsChr[(i-window):(i+ window)]))",sep="")))	
    }#end i loop
    
    MedRChr[1:window] <- eval(parse(text=paste(Func,"(na.omit(RsChr[1:window]))",sep="")))
    MedRChr[(length(RsChr)-window):length(RsChr)] <- eval(parse(text=paste(Func,"(na.omit(RsChr[(length(RsChr)-window):length(RsChr)]))",sep="")))
    MedR <- rbind(MedR,MedRChr)
    
    #print(paste("Chromosome",chr))
  }#end chr loop
  
  MedRs<- cbind(MedRs,MedR)
  print(paste(ind,"is normalized"))
  
}#end ind loop


logRsMed <- cbind(logRsRaw[,c(1:3)],MedRs)
colnames(logRsMed)<- colnames(logRsRaw)
#ParScore<-load(paste(outPath,"ParScore_",Family,".rda",sep=""))
#ParScore<-eval(parse(text=ParScore))


#--------------------------------------
#---------- window on GC

for(chr in unique(logRsRaw$Chr)){
  
  MedGcChr <-matrix(data=NA,sum(logRsRaw$Chr==chr))
  GcChr <- GC[logRsRaw$Chr==chr,6]
  for(i in (1+window):(length(GcChr)-window)){
    
    MedGcChr[i] <- eval(parse(text=paste(Func,"(na.omit(GcChr[(i-window):(i+ window)]))",sep="")))	
    
  }#end i loop
  
  MedGcChr[1:window] <- eval(parse(text=paste(Func,"(na.omit(GcChr[1:window]))",sep="")))
  MedGcChr[(length(GcChr)-window):length(GcChr)] <- eval(parse(text=paste(Func,"(na.omit(GcChr[(length(GcChr)-window):length(GcChr)]))",sep="")))
  
  if(chr==unique(logRsRaw$Chr)[1]){MedGc <- MedGcChr}else{MedGc <- rbind(MedGc,MedGcChr)}
  
  print(paste("Chromosome",chr))
}#end chr loop

GC2 <- cbind(logRsRaw[,c(1:3)],MedGc)
colnames(GC2)[4]<-"GC"

source("wavecorrtmean2.R")
logRs <- wavecorrtmean2(logRsMed,ParScore,GC2,Family,outPath)
logRs

}#end function

