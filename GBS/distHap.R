#outPath <- "/uz/data/Admin/cme_genome_calculated/"
#Family = "PGD047_C1"

#Int <- rbind(c(0,0,0,"Pat"),c(7,117105838,117356025,"Mat"))

disthap<-function(Family,Int,outPath){

Haps <- read.table(paste(outPath,"/",Family,"_Itp2.hap",sep=""),header=T,sep="\t",stringsAsFactors =F)



Chroms <-Int[Int[,1]!="0",1]
Scs <- colnames(Haps)[grep("E*_Bl*",colnames(Haps))]


for(ind in Scs){
	
	for(chr in Chroms){
		
		IntChr <- Int[Int[,1]==chr,]	
		#ChrLength <- as.numeric(ChrsLength[ChrsLength[,1]==paste("chr",chr,sep=""),2])
		HapChr <- Haps[Haps[,"Chr"]==chr,c("Position",ind)]	
		UP <- HapChr[HapChr$Position<=as.numeric(IntChr[2]),]
		Down <- HapChr[HapChr$Position>=as.numeric(IntChr[3]),]

		UPRle <- rle(UP[,ind])
		DownRle <- rle(Down[,ind])

		#if(UPRle$values[length(UPRle$values)]==DownRle$values[1]){
			Dist <- cbind(ind,chr,UPRle$lengths[length(UPRle$values)],DownRle$lengths[1],UPRle$values[length(UPRle$values)], DownRle$values[1])#}

		
		if(chr == Chroms[1]){Dists <- Dist}else{Dists<- rbind(Dists,Dist)}
	
	}#end chr
	
	if(ind==Scs[1]){DistFam <- Dists }else{DistFam <-rbind(DistFam,Dists)}
	print(ind)
		
}#end ind loop

colnames(DistFam) <- c("Haplotype","Chr","LengthUp","LengthDown","ValueUp","ValueDown")

write.table(DistFam,paste(outPath,"DistanceToHR_",Family,".txt",sep=""),quote=F,sep="\t",col.names=T,row.names=F)
}
