library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
library(VariantAnnotation)
library(vcfR)
library(intrval)


pct <- 100

meb <- c("GM")
#,"PGD230","PGD111","PGD058","PGD051","PGD088","PGD146")
	# "GBS_PGD189_C3","GBS_PGD090_C3","GBS_PGD440_C1","GBS_PGD068_C3","GBS_PGD330_C1")     
seed <- c("Mother")


for (i in 1:length(meb)) {
	id <- meb[i]
	print(id)

	wd <- "/uz/data/hydra/genomicscore/raw/HiSeqComputed/new/diagnostics/uz/cme_pgd/gcpu/samples/"

	# wd <- "/uz/data/avalok/symbiosys/gcpi_d_uz/cme_pgd/gcpu/runs/"
	fl <- paste0(wd,"/GBS_PGD_output/merge.vcf.pct",pct,"/",id,".CHR.merge.vcf")	# for PGD230 ect
	# fl <- paste0(wd,"/GBS_PGD_output/",id,"/merge/",id,".CHR.merge.freebayers.vcf") # 	paste0(wd,"/GBS_PGD_output/",id,"/merge/",id,".CHR.merge.vcf") # gatk

	configure_fn <- paste0(wd,"GBS_conf/",id,"_conf.txt")
	CN <- read.table(configure_fn, header=F, stringsAsFactors=F)
	cnv1 <- NULL
	cnv2 <- cnv1
	for (f in 1:nrow(CN)){
		fn1 <- paste0(wd,"GBS_",CN[f,1],"/align_bwa/current/build_bam/current/QDNAseq/current/result/", "GBS_",CN[f,1],"_100_copyNumbersSmooth.bed")
		fn2 <- paste0(wd,"GBS_",CN[f,1],"/align_bwa/current/build_bam/current/QDNAseq/current/result/", "GBS_",CN[f,1],"_100_seg.bed")
		b1 <- read.table(fn1, stringsAsFactors=F, skip=1)
		b2 <- read.table(fn2, stringsAsFactors=F, skip=1)
		cnv1 <- cbind(cnv1,b1[,5])
		cnv2 <- cbind(cnv2,b2[,5])
	}
	colnames(cnv1) <- paste0(CN[,4],"_",id)
	colnames(cnv2) <- paste0(CN[,4],"_",id)
	cnv1 <- cbind(b1[,1:3],cnv1) 	
	colnames(cnv1)[1:3] <- c("Name","Chr","Position")
	cnv2 <- cbind(b1[,1:3],cnv2) 
	colnames(cnv2)[1:3] <- c("Name","Chr","Position")
	
	out_dir <- file.path(wd,"/GBS_PGD_output/", id)
	if (!file.exists(out_dir)){dir.create(file.path(wd, id),showWarnings = FALSE)}	
	
	id <- gsub("GBS_","",id)
	out_fn1 <- paste0(out_dir,"/output/",id,"_logRrawGBS.txt")
	out_fn2 <- paste0(out_dir,"/output/",id,"_logRsegGBS.txt")
	print(out_fn2)
	write.table(cnv1, out_fn1, row.names=F, col.names=T, quote=F, sep="\t")
	write.table(cnv2, out_fn2, row.names=F, col.names=T, quote=F, sep="\t")
	
	chooseCRANmirror(ind=8)

	if("VariantAnnotation" %in% rownames(installed.packages()) == FALSE) {install.packages("VariantAnnotation")}
	if("vcfR" %in% rownames(installed.packages()) == FALSE) {install.packages("vcfR")}
	if("SNPlocs.Hsapiens.dbSNP144.GRCh37" %in% rownames(installed.packages()) == FALSE) {
	source("https://bioconductor.org/biocLite.R")
	biocLite("SNPlocs.Hsapiens.dbSNP144.GRCh37")
}




dps <- function(X){
	Y <- do.call("rbind",strsplit(X,","))[,1:2]
	Y2 <- as.numeric(Y[,2])/(as.numeric(Y[,1])+as.numeric(Y[,2]))
	Y2
} #end function


SNParray <- NULL
hap <- NULL


for (chr in c(1:22,"X","Y")) { 
	# print(chr)
	vcf_fl <- gsub("CHR",chr,fl)
	print(vcf_fl)
	vcf <- readVcf(vcf_fl, "hg19")

	Name <- names(vcf)
	Position <- start(vcf)
	Chr <- rep(chr,length(Position))
	res <- genotypeToSnpMatrix(vcf,uncertain=FALSE)
	snparray <- t(as(res$genotype, "character"))
	
# generate GT value
if (chr %in% c("X","Y") ) {
        GT <- (snparray)
	male_ncol <- NULL
	female_ncol <- NULL
        male_ncol <-  grep("12877_MC",colnames(geno(vcf)$GT))
        male_ncol <-  c(male_ncol,grep("12882_MC",colnames(geno(vcf)$GT)))
        male_ncol <-  c(male_ncol,grep("_F",colnames(geno(vcf)$GT)))
        male_ncol <-  c(male_ncol,grep("_MGF",colnames(geno(vcf)$GT)))
        male_ncol <-  c(male_ncol,grep("_PGF",colnames(geno(vcf)$GT)))
        GT_sex <- (geno(vcf)$GT)
        GT_sex[,male_ncol][GT_sex[,male_ncol]==0] <- "A/A"
        GT_sex[,male_ncol][GT_sex[,male_ncol]==1] <- "B/B"
        GT_sex <- as.data.frame(GT_sex)
        GT_sex[,-c(male_ncol)] <- as.data.frame(GT[,-c(male_ncol)])
        GT <- GT_sex
} else {
        GT <- as.data.frame(snparray)
}

	for(c in 1:ncol(GT)){GT[,c] <- gsub("/","",GT[,c]);GT[,c] <- gsub("NA","NC",GT[,c])}


	# define column names of GT
	CN <- read.table(configure_fn, header=F, stringsAsFactors=F)

	# in case of onePGT samples, the sample folder names are GC071417_AAACATCG, need to replace it to the real name
	mycolumns <- colnames(GT)
	myobject <- CN[,2]
	for(i in myobject){
		mycolumns[grepl(i, mycolumns)] <- CN[CN[,2] %in% i,1]
	}
	colnames(GT) <- paste0("GBS_",mycolumns)


	for (i in 1:length(colnames(GT))){
		colnames(GT)[i] <- gsub("GBS_","",colnames(GT)[i])
		colnames(GT)[i] <- gsub("$",".GType",CN[CN[,1] %in% colnames(GT)[i],4])	# the 4th column in the conf.txt is the colnames of each feature.
	}


	# generate BAF
	vcf <- read.vcfR (vcf_fl, "hg19")
	A <- extract.gt(vcf,element='AD')
	BAF <- apply(A,2,dps)
	colnames(BAF)<- gsub(".GType",".B Allele Freq",colnames(GT))


	# generate logR
	logR <- BAF
	logR <- matrix(runif((ncol(BAF)*nrow(BAF)), 0, 1), ncol=ncol(BAF))
	colnames(logR) <- gsub(".B Allele Freq",".Log R Ratio",colnames(BAF))

	SNParray <- rbind(SNParray,cbind(Name,Chr,Position,GT,BAF,logR))

}

SNParray[SNParray == "."] <- "NC"
SNParray_org <- SNParray


sub <- SNParray[(SNParray$Father.GType == "AB" & (SNParray$Mother.GType == "AA" | SNParray$Mother.GType == "BB")) | (SNParray$Mother.GType == "AB" & (SNParray$Father.GType == "AA" | SNParray$Father.GType == "BB")),]
#b_seed <- gsub("$",".GType",seed[i])
#o_seed <- gsub("$",".GType",ifelse(seed[i] == "Mother","Father","Mother"))
sample_name_F <- sample(as.character(sub[(sub["Father.GType"] == "AA" ),1]),size=(nrow(sub[(sub["Mother.GType"] == "AA" ),])-nrow(sub[(sub["Mother.GType"] == "BB" ),])))
sample_name_M <- sample(as.character(sub[(sub["Mother.GType"] == "AA" ),1]),size=(nrow(sub[(sub["Father.GType"] == "AA" ),])-nrow(sub[(sub["Father.GType"] == "BB" ),])))
SNParray <- (SNParray[SNParray[,1] %ni% c(sample_name_F,sample_name_M),])


id<-gsub("GBS_","",id)

out_fn <- paste0(out_dir,"/",id,".adj")
print(out_fn)


write.table(SNParray,out_fn,row.names=F, col.names=T,quote=F, sep="\t")



}


