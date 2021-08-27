args <- commandArgs(TRUE)


bamDir <- args[1]	
OutDir <- args[2]	
bin_size <- as.numeric(args[3])
disease_loci <- args[4]


	library("BSgenome.Hsapiens.UCSC.hg19")
	library(Matrix)
	library(sqldf)
	library(reshape2)
	library(reshape)
	library(matrixStats)
	library(Biostrings)
	library(Rsamtools)
	library(Matrix)
	library(dplyr)
	library(BiocGenerics)
	library(Biobase)
	library(DNAcopy)
	library(CGHregions)
	library(BSgenome.Hsapiens.UCSC.hg19)
	library(QDNAseq.hg19)
	library(QDNAseq)



dir.create(file.path(OutDir))


chrs <- c(1:22,"X","Y")
bins <- getBinAnnotations(binSize=bin_size)
chrLength <- seqlengths(Hsapiens)

print(" loading done ! < -----------")
bamFile <- list.files(bamDir,"\\.bam$")
bamID <- head(unlist(strsplit(tail(unlist(strsplit(bamFile, "/")),1),"[.]")),1)


bamcoverage <- function (bamfile) {
	  # read in the bam file
	  bam <- scanBam(bamfile)[[1]] # the result comes in nested lists
	  # filter reads without match position
	  ind <- ! is.na(bam$pos)
	  ## remove non-matches, they are not relevant to us
	  bam <- lapply(bam, function(x) x[ind])
	  ranges <- IRanges(start=bam$pos, width=bam$qwidth, names=make.names(bam$qname, unique=TRUE))
	  ## names of the bam data frame:
	  ## "qname"  "flag"   "rname"  "strand" "pos"    "qwidth"
	  ## "mapq"   "cigar"  "mrnm"   "mpos"   "isize"  "seq"    "qual"
	  ## construc: genomic ranges object containing all reads
	  ranges <- GRanges(seqnames=Rle(bam$rname), ranges=ranges, strand=Rle(bam$strand), flag=bam$flag, readid=bam$rname )
	  ## returns a coverage for each reference sequence (aka. chromosome) in the bam file
	  return (mean(coverage(ranges)))      
}  




loop_samples <- function(bins,InputFile,BIN_S,bamID) {

print(paste0("============= InputFile:", InputFile))
	OutFile <- paste0(OutDir,"/",bamID,"_","avgCov",".txt")
	write.table(bamcoverage(InputFile),OutFile,sep="\t",col.names=F)

print(paste0("============= applyFilters"))
	readCounts <- applyFilters(readCounts,residual=TRUE, blacklist=TRUE)

print(paste0("============= estimateCorrection"))
	readCounts <- estimateCorrection(readCounts)

print(paste0("============= applyFilters correctBins normaizeBins smoothOutlierBins"))
	readCounts <- applyFilters(readCounts,residual=TRUE, blacklist=TRUE, chromosomes=NA)	# includes both X & Y
	copyNumbers <- correctBins(readCounts)
	copyNumbersNormalized <- normalizeBins(copyNumbers)
	copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)

print(paste0("============= segmentBins & normalizeSegmentedBins"))
	copyNumbersSegmented <- segmentBins(copyNumbersSmooth, transformFun="sqrt")
	copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)

print(paste0("============= exportBins"))
     	OutFile <- paste0(OutDir,"/",bamID,"_",BIN_S,"_seg.bed")
     	exportBins(copyNumbersSegmented,file = OutFile, format="bed", filter=T,type=c("segments"))

print(paste0("============= plot genome-wide ============="))
	OutFile <- paste0(OutDir,"/",bamID,"_",BIN_S,".pdf")
	print(OutFile)
	pdf(OutFile,paper='a4r')
        plot(copyNumbersSegmented, plot.type = "w", main = paste0(unlist(strsplit(bamID, "-"))[1]," ",unlist(strsplit(bamID, "-"))[2],",bin=",BIN_S)) 
	dev.off()

	OutFile <- paste0(OutDir,"/",bamID,"_",BIN_S,"header.csv")
	print(OutFile)
	write.table(copyNumbersSegmented@phenoData@data, file = OutFile, row.names = T, sep = "\t", quote = FALSE, append = FALSE)

print(paste0("============= plot per_chr ============="))
	f.data <- as(featureData(copyNumbersSegmented), "data.frame")
	OutFile <- paste0(OutDir,"/",bamID,"_","byChr",BIN_S,".pdf")
	pdf(OutFile)
	
	if (!is.na((grep("^([0-9]+)",disease_loci)) || (grep("X",disease_loci) ==1 ) || (grep("Y",disease_loci) ==1))) {
		dl <- read.table(text=disease_loci,col.names=c('chr','pos'))
	} else if (!is.na(((grep("^N",disease_loci))==1) || ((grep("^GC",disease_loci))==1))) {
		dl <- NA
	}

		
	for (chromosome in unique(f.data$chromosome)) {
		select <- which(f.data$chromosome == chromosome)
		print(paste("Plotting chromosome:", chromosome, sep = " "))
		## Looping over samples, make a plot per sample, per chromosome
		for (c in 1:ncol(copyNumbersSegmented)){
			sample.name <- sampleNames(phenoData(copyNumbersSegmented))[c] # You can use this, for plot title or plot name
			sample.name <- bamID
			plot(copyNumbersSegmented[select, c], ylim=c(-2,2), main=paste0(unlist(strsplit(bamID, "-"))[1]," ",unlist(strsplit(bamID, "-"))[2]," Chr",chromosome,",bin=",BIN_S))         

		}
			
		if (is.data.frame(dl)) {
			if(chromosome %in% dl[,1]) {
				abline(v=dl[dl$chr == chromosome,2],lty=2,col="orange")	
			}
		}
	}
	dev.off()

}


print(paste0("============= run QDNAseq main ============="))
loop_samples(bins,paste(bamDir,bamFile,sep="/"),bin_size,bamID)



