# REQUIRED:
# /uz/data/Admin/cme_genome_raw/PGD_Exported/PGD043_C1.txt
# /uz/data/Admin/cme_genome_raw/PGD_Exported/PGD043_C1_Parameters.txt
# /uz/data/Admin/cme_genome_raw/PGD_Exported/PGD043_C1_Intervals.txt
# 
# sed -e 's/,/\./g' /uz/data/Admin/cme_genome_raw/PGD_Exported/PGD043_C1.txt > /uz/data/Admin/cme_genome_raw/PGD043_C1.adj



#--------------------------------------------------------
#-------      Global variables         ------------
library(limma)
library(signal)
library(plotrix)
library(MASS)
library(signal)


args <- commandArgs(TRUE)
Family <- args[1]
dataFile <- args[2]
parametersFile <- args[3]
IntervalFile <- args[4]
outPath <- args[5]
siCHILD_DIR <- args[6]
PGD_EXPORTED_DIR <- paste0(outPath,Family)

if (!grepl('/$',outPath)) {outPath <- paste(outPath,"/",sep="")}

print("Loading sources...")
load(paste(siCHILD_DIR,"Ideogram_hg19.rda",sep=""))
load(paste(siCHILD_DIR,"REF_24h_QC_illuminaCytoSNP12.rda",sep=""))
source(paste(siCHILD_DIR,"PlotCyto.R",sep=""))
source(paste(siCHILD_DIR,"fastPCF.R",sep=""))
source(paste(siCHILD_DIR,"po2.R",sep=""))
source(paste(siCHILD_DIR,"qcgtype.R",sep=""))
source(paste(siCHILD_DIR,"qcbyparents.R",sep=""))
source(paste(siCHILD_DIR,"avgint.R",sep=""))
source(paste(siCHILD_DIR,"callpcf.R",sep=""))
source(paste(siCHILD_DIR,"chrxhtyping1.R",sep=""))
source(paste(siCHILD_DIR,"testmedfilt.R",sep=""))
source(paste(siCHILD_DIR,"inthapnew1.R",sep=""))
source(paste(siCHILD_DIR,"patscore2.R",sep=""))
source(paste(siCHILD_DIR,"artscsnpspar.R",sep=""))
source(paste(siCHILD_DIR,"htypingAutOpt2.R",sep=""))
source(paste(siCHILD_DIR,"meanwindow.R",sep=""))
source(paste(siCHILD_DIR,"htypingOpt1.R",sep=""))
source(paste(siCHILD_DIR,"htypingOpt2.R",sep=""))
source(paste(siCHILD_DIR,"intphappropser3.R",sep=""))
source(paste(siCHILD_DIR,"phasebaf2.R",sep=""))
source(paste(siCHILD_DIR,"phasebaf4.R",sep=""))
source(paste(siCHILD_DIR,"chrspecbafplot.R",sep=""))
source(paste(siCHILD_DIR,"genomebafplotsc.R",sep=""))
source(paste(siCHILD_DIR,"distHap.R",sep=""))


#Default parameter settings
options(scipen=999)
Func = "mean"
Chroms <- c(1:22,"X")
SibPattern = "E*_Bl*"

print("Loading parameters and intervals...")
Params <- read.table(parametersFile,sep="\t",header=T,stringsAsFactors=F)
Parent = Params[Params$Param=="Parent","Value"]
Parent1 = paste(Parent,"_",Family,sep="")
GC_File = Params[Params$Param=="GC_File","Value"]

for(p in c("GC_File")){Params <- Params[Params$Param != p,]}

for(i in 1:nrow(Params)) { if(Params[i,"Param"]=="Parent" | Params[i,"Param"]=="Seed"  | Params[i,"Param"]=="GC_File") {eval(parse(text=paste(Params[i,"Param"],"='",Params[i,"Value"],"'",sep="")))} 
			   else {eval(parse(text=paste(Params[i,"Param"],"=",Params[i,"Value"])))}
			 }


#Intervalfile <-  paste(PGD_EXPORTED_DIR,Family,"_Intervals.txt",sep="")
if(ExcInt==1){Int <- read.table(IntervalFile,sep="\t",header=T,stringsAsFactors=F)}
Int <- Int[complete.cases(Int),]


# E02_Bl001.GType : Embryo2 blastomere1 (var nr of cols, x per blastomere)
print(paste("Reading data file:",dataFile))
data <- read.table(dataFile,sep="\t",header=T)

for(c in c(1,2,grep("GType",colnames(data)))){data[,c]<-as.character(data[,c])}
for(c in c(grep("Log.R.Ratio",colnames(data)))){data[,c]<-as.numeric(as.character(data[,c]))}
for(c in c(grep("B.Allele.Freq",colnames(data)))){data[,c]<-as.numeric(as.character(data[,c]))}
dataraw <- na.omit(data)


dataraw <- dataraw[grep("cnvi",as.character(dataraw$Name),invert=TRUE),]
dataraw <- dataraw[dataraw[,"Position"]!=0,]
dataraw <- dataraw[dataraw[,"Chr"]!="Y",]
dataraw <- dataraw[dataraw[,"Chr"]!="XY",]

GC <- read.table(GC_File,header=F,sep="\t")
rownames(dataraw) <- as.character(dataraw[,1])
rownames(GC) <- as.character(GC[,4])
rowsTot <- intersect(rownames(dataraw),rownames(GC))
GC <- GC[rowsTot,]
dataraw <- dataraw[rowsTot,]


logRsRaw<- cbind(dataraw[,c("Name","Chr","Position")],dataraw[,grep(".Log.R.Ratio",colnames(dataraw))])
colnames(logRsRaw)[-c(1:3)] <- paste(gsub(".Log.R.Ratio","",colnames(logRsRaw)[-c(1:3)]),"_",Family,sep="")	
write.table(logRsRaw,paste(outPath,Family,"_logRsRaw.txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")

BAFs <- cbind(dataraw[,c("Name","Chr","Position")],dataraw[,grep(".B.Allele.Freq",colnames(dataraw))])
colnames(BAFs)[-c(1:3)] <- paste(gsub(".B.Allele.Freq","",colnames(BAFs)[-c(1:3)]),"_",Family,sep="")
write.table(BAFs,paste(outPath,Family,"_BAFs.txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")


Gtypes<- cbind(dataraw[,c("Name","Chr","Position")],dataraw[,grep(".GType",colnames(dataraw))])
colnames(Gtypes) <- colnames(logRsRaw)
write.table(Gtypes,paste(outPath,Family,"_0.75.gtp",sep=""),col.names=T,row.names=F,quote=F,sep="\t")


ChrPos <- Gtypes[,c("Chr","Position")]
QC <- qcgtype(Gtypes,ChrPos,Family,SibPattern,outPath) 
QCbyParents <- qcbyparents(Gtypes,SibPattern)
SnpSpecArts <- artscsnpspar(Gtypes,SibPattern,outPath)
SnpSpecArts[,"Position"] <- as.numeric(SnpSpecArts[,"Position"])

dataPo <- po2(Gtypes,Family,SibPattern,outPath)

ParScore <- patscore2(dataPo,QC,Chroms,Gtypes,Family) 

names(ParScore)<-gsub("E00_Bl000_","",names(ParScore))

save(ParScore,file=paste(outPath,"ParScore_",Family,".rda",sep=""))


library(limma)
window=Win/2
logRs <- meanwindow(logRsRaw,GC,window,Func,Family,ParScore,outPath)
write.table(logRs,paste(outPath,Family,"_logRsAvgWindow.txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")

data <- dataraw
gamma=gammaBAF

library(signal)

if(Seed=="Grandparents"){
	Haps <- htypingOpt1(Gtypes,dataPo,ParScore,Parent1,SibPattern,Family,outPath)
	print("-------------------------------------")
	print("Option one haplotyping was applied...")
	print("-------------------------------------")
	PhBAF <- phasebaf2(Gtypes,dataraw,Family,Parent1,gamma,outPath)
}else{
	Haps <- htypingOpt2(Gtypes,Window,Int,dataPo,ParScore,SibPattern,Family,outPath)
	print("-------------------------------------")
	print("Option two haplotyping was applied...")
	print("-------------------------------------")
	PhBAF <- phasebaf4(Gtypes,dataraw,Family,gamma,ParScore,outPath)
}

dataHap <- Haps[["dataHap"]]
dataHapRaw <- Haps[["dataHapRaw"]] 

write.table(dataHap,paste(outPath,Family,"_Itp.hap",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
write.table(dataHapRaw,paste(outPath,Family,"_Raw.hap",sep=""),col.names=T,row.names=F,quote=F,sep="\t")

AvgLogRs <- avgint(logRs,Int,Family,outPath)
SegLogRs <- callpcf(logRs,gammaSC,gammaMC,plateau,Family,outPath)


library(MASS)
Intp <- intphappropser3(outPath,Family,ParScore,SibPattern)

disthap(Family,Int,outPath)

# PpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPp
# Plotting

chrspecbafplot(dataHap,dataHapRaw,dataPo,BAFs,logRs,logRsSeg,P1,P1Seg,P2,P2Seg,M1,M1Seg,M2,M2Seg,ChrsLengths,ideogram,Family,outPath,Chroms,Int) # Chr specific plotting
genomebafplot(BAFs,logRs,logRsSeg,P1,P1Seg,P2,P2Seg,M1,M1Seg,M2,M2Seg,ChrsLengths,ideogram,Family,outPath)

#PpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPp

