#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%   	      					   Function: GC correction Median             	   																%%%%%%%
#
#Author: MZE (RGLab)
#cDate: 
#mDtae: 
#
#(<-): dataFile outPath
# 
#
#(->): GC-corrected file
#
#Description: 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wavecorrtmean2 <- function(logRs,SamplesPo2,GC,Family,outPath){

print("Applying wave correction with trimmed mean adjustment ...")


#logRs <- logRsRaw
#rownames(logRs) <- as.character(logRs[,"Name"])
#rownames(GC) <- as.character(GC[,"Name"])
rowsTot <- intersect(rownames(logRs),rownames(GC))

data <- logRs[rowsTot,]


GC <-data.frame(GC[rowsTot,])
GC[,"GC"] <- as.numeric(as.character(GC[,"GC"]))

for(ind in colnames(logRs)[-c(1:3)]){
	
	logrBe <- data[,ind]	
	while(sum(is.na(logrBe))>=1){logrBe[which(is.na(logrBe))]<-logrBe[which(is.na(logrBe))-1]}
	FitLogR <- loessFit(logrBe,GC[,"GC"])
	#logrAfNorm <- ((logrBe - FitLogR$fitted)-mean(na.omit(logrBe - FitLogR$fitted)))/sd(na.omit(logrBe - FitLogR$fitted))
	DipChrs <- na.omit(rownames(SamplesPo2[[ind]])[SamplesPo2[[ind]][,"Par"]==12])
	if(length(DipChrs)>=1){
		
		for(chr in DipChrs){
			DipPo <- which(data[,"Chr"]==chr)
			if(chr==DipChrs[1]){DipPos<- DipPo}else{DipPos<-c(DipPos,DipPo)}
			}
	logrAfNorm <- (logrBe - FitLogR$fitted)-mean(na.omit((logrBe - FitLogR$fitted)[DipPos]),trim=0.10)
	print(paste("timmed mean value is",mean(na.omit((logrBe - FitLogR$fitted)[DipPos]),trim=0.10)))
        print(paste("timmed mean value after meand and GC-correction is",mean(logrAfNorm[DipPos],trim=0.10)))
	print(paste("mean of GC-corrected logRs is",mean(logrBe - FitLogR$fitted)))

	}else if (sum(grep("E*_Bl*",ind))==0){
		
	logrAfNorm <- (logrBe - FitLogR$fitted)-mean(na.omit((logrBe - FitLogR$fitted)[-c(which(data[,2]=="X" | data[,2]=="Y" | data[,2]=="XY"))]),trim=0.10)
		print(paste(ind, "is a multi-cell sample"))
	}else{	
		logrAfNorm <- (logrBe - FitLogR$fitted)
		print(paste("No diploid chromosome in",ind))
	}
	
	data[,ind]<- logrAfNorm
 		
	print(ind)
}#end ind loop

#write.table(data,paste(outPath,Family,"_logRs_WaveCorr_Mean.txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")

print("GC-corrected file was written")

data

}#end function

