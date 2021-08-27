args <- commandArgs(TRUE)


bamDir <- args[1]	
OutDir <- args[2]	
bin_size <- as.numeric(args[3])
disease_loci <- args[4]
opt <- as.character(args[6])	# singlecell = yes


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


if (opt != "yes") {
buildCNVtable <- function(copyNumbersOutFile, cytobandfile, bamID, BIN_S) {
	cytobands <- read.csv(cytobandfile, header=F,sep="\t",na.strings=c("","NA"), comment.char = "#", as.is=TRUE, stringsAsFactors=FALSE)
	cytobands <- cbind(as.data.frame(apply(cytobands[,c(1:5)], 2, as.character)), stringsAsFactors=FALSE)
	cytobandsColNames<-c("chr","start","end","cytob","value")
	colnames(cytobands)<-(cytobandsColNames)
	cytobands$chr <- as.character(gsub("chr", "", cytobands$chr))

	copyNumbersSmooth <- read.csv(copyNumbersOutFile, header=F,sep="\t",na.strings=c("","NA"), comment.char = "#", skip = 1, as.is=TRUE, stringsAsFactors=FALSE)
	head(copyNumbersSmooth)
	ListColNames<-c("chr","start","end","Nameshort","value","strand")
	colnames(copyNumbersSmooth)<-(ListColNames)
	copyNumbersSmooth$variable <- paste0(bamID,"_",BIN_S)
	dft3 <- f_add_cluster_v2(copyNumbersSmooth)
	dft3 <- dft3[complete.cases(dft3),]
	dft3 <- cbind(as.data.frame(apply(dft3[,c(1:6)], 2, as.character)), as.data.frame(apply(dft3[,c(7:8)], 2, as.numeric), stringsAsFactors=FALSE))

	####################################add cytoband
	#build SQL
	cytobands <- cbind(as.data.frame(apply(cytobands[,c(1,4:5)], 2, as.character)), as.data.frame(apply(cytobands[,c(2:3)], 2, as.numeric), stringsAsFactors=FALSE))
	sql2 <- sprintf("select A.*, GROUP_CONCAT(B.cytob) as cytobands
    	            from dft3 A, cytobands B
        	        where ( A.chr like B.chr and (A.start+1) >= B.start*1 and (A.start+1) < B.end*1 and A.end*1 > B.start*1 and A.end*1 < B.end*1 ) or
            	    ( A.chr like B.chr and (A.start+1) >= B.start*1 and A.end*1 > B.end*1 and (A.start+1) < B.end*1 ) or
                	( A.chr like B.chr and (A.start+1) < B.start*1 and A.end*1 < B.end*1 and A.end*1 >= B.start*1 )
                	group by A.variable, A.chr, A.start, A.end;")

	dft3cyto <- sqldf(sql2)

	####################################group clusters together by cluster id and calculate stats
	#build SQL
	sql <- sprintf("select A.chr, min(A.start*1) as start, max(A.end*1) as end, max(A.end*1)-min(A.start*1) as Size, A.variable as sample, A.cluster, (SUM(A.value)/count(A.value)) as AVG_score,
				   GROUP_CONCAT(printf('%%.2f', A.value)) as scores,
				   GROUP_CONCAT(A.cytobands) as cytobands,
				   GROUP_CONCAT(A.Nameshort) as regions, count(A.Nameshort) as regcount,
				   (CASE WHEN (SUM(A.value)/count(A.value)) >= 0.3 THEN 'GAIN'
				   WHEN (SUM(A.value)/count(A.value)) <= -0.3 THEN 'LOSS'
				   ELSE 'NotDef' END) as DupDel
				   from dft3cyto A
				   group by A.variable, A.chr, A.cluster;")

	sqlres <- sqldf(sql)
	sqlres$Size <- paste(format(round(sqlres$Size / 1e6, 1), trim = TRUE), "Mb")
	sqlres$startV <-  format(sqlres$start,big.mark=",", trim=TRUE)
	sqlres$endV <- format(sqlres$end,big.mark=",", trim=TRUE)

	sqlres <- cbind(as.data.frame(apply(sqlres[,c(1:14)], 2, as.character)), stringsAsFactors=FALSE)
    dft3cyto <- cbind(as.data.frame(apply(dft3cyto[,c(1:9)], 2, as.character)), stringsAsFactors=FALSE)

	####################################get first and last cytoband
	#build SQL
	sql <- sprintf("select A.chr, A.start, A.end, A.Size, A.Sample, A.AVG_score, A.regcount, A.DupDel,'arr[hg19]'||A.chr||B.cytobands||C.cytobands||'('||A.startV||'-'||A.endV||')' as ISCN_notation,
                       B.cytobands as minCyto, C.cytobands as maxCyto, A.scores, A.regions
              from sqlres A, dft3cyto B, dft3cyto C
              where A.sample like B.variable and A.chr like B.chr and A.cluster*1 == B.cluster*1 and A.start*1 == B.start*1 and
                A.sample like C.variable and A.chr like C.chr and A.cluster*1 == C.cluster*1 and A.end*1 == C.end*1 ;")

	sqlres_out <- sqldf(sql)
	return(sqlres_out)
}



#####
#function to get cluster numbering for each sample
#multisample data requires sorting upfront on dft2
#columns chr,start,end,Nameshort,rcovmean,covsd,covsem,covIQR,rmean,akkdist,sd,sem,IQR,variable,value

f_add_cluster_v2 <- function(clccdf) {
  library(dplyr)
  library(data.table)
  clccdf$chr<-as.character(clccdf$chr)
  clccdf$sample<-as.character(clccdf$variable)
  clccdf$samplec<-as.numeric(as.factor(clccdf$variable))
  clccdf$chrc<-as.numeric(as.factor(clccdf$chr))

  #add preceeding and following value
  dft <- clccdf %>%
    group_by(samplec + chrc ) %>%
    mutate(nscore = lead(value, 1),
           bscore = lag(value, 1))
  dft$samechr<-1
  dft$newchr<-2
  dft$switchchr <- ifelse(!is.na(dft$bscore), dft$samechr, dft$newchr)

  ######################testing
  #l low nl new low ll last low
  #h high nh new high lh last high
  dft$hln <-       ifelse(dft$value>=0.375 & dft$bscore>=0.375,dft$hln<-"h",
                          ifelse(dft$value>=0.375 & dft$bscore<0.375,dft$hln<-"nh",
                                 ifelse(dft$value>=0.375 & dft$nscore>=0.375,dft$hln<-"h",
                                        ifelse(dft$value>=0.375 & dft$nscore<0.375,dft$hln<-"h",
                                               ifelse(dft$value<=-0.375 & dft$bscore<=-0.375,dft$hln<-"l",
                                                      ifelse(dft$value<=-0.375 & dft$bscore>-0.375,dft$hln<-"nl",
                                                             ifelse(dft$value<=-0.375 & dft$nscore<=-0.375,dft$hln<-"l",
                                                                    ifelse(dft$value<=-0.375 & dft$nscore>-0.375,dft$hln<-"l","n"))))))))

  dft$hlnn <-       ifelse(dft$value>=0.375 & dft$switchchr == 2, dft$hlnn<-"nh",
                            ifelse(dft$value<=-0.375 & dft$switchchr == 2, dft$hlnn<-"nl",dft$hln))

  dft$hln <- dft$hlnn
  #sfCat(head(dft$hln, "\n"))
  #store row index for all
  dft$rowNa <- rownames(dft)
  dft <- as.data.frame(dft)

  #subset the seeds of cluster and add clusternumber
  dfts <- dft[which(dft$hln == "nh" | dft$hln == "nl"),]
  dfts$rowNa <- rownames(dfts)
  dfts$seqnh <- seq.int(nrow(dfts))

  #subset the cluster extensions
  dfthl <- dft[which(dft$hln == "h" | dft$hln == "l"),]
  dfthl$rowNa <- rownames(dfthl)

  #add cluster number of closest seed
  dfthl <- cbind(dfthl, as.integer(sapply(dfthl$rowNa, function(x) GetSmallerClusterID(x, dfts), USE.NAMES=F)))
  colnames(dfthl)<-c(colnames(dfthl[1:ncol(dfthl)-1]),"seqnh")

  #merge data
  dfttom <- rbind(dfts,dfthl)
  dftm <- merge(dft, dfttom, by="rowNa", all.x=TRUE)
  #drop columns
  dftm <- as.data.frame(dftm)
  dftm <- dftm[,c("chr.x","start.x","end.x","Nameshort.x","strand.x","variable.x","value.x","seqnh")]
  colnames(dftm) <- c("chr","start","end","Nameshort","strand","variable","value","cluster")
  #return
  return(dftm)
}


#helper function to get the smaller cluster ID used in f_add_cluster_v2

GetSmallerClusterID <- function(x, dfts) {
  clind <- sum(as.numeric(dfts$rowNa) <= as.numeric(x))
  clnr <- as.integer(dfts$seqnh[clind])
  return(clnr)
}
}


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
print(paste0("============= binReadCounts"))
	readCounts <- binReadCounts(bins, bamfiles=InputFile)
	OutFile <- paste0(OutDir,"/",bamID,"_",BIN_S,"_","highlightFilters",BIN_S,".pdf")
	pdf(OutFile)
	plot(readCounts, logTransform=FALSE, ylim=c(-200, 1000))
	highlightFilters(readCounts, logTransform=FALSE,residual=TRUE, blacklist=TRUE)
	dev.off()

print(paste0("============= applyFilters"))
	readCounts <- applyFilters(readCounts,residual=TRUE, blacklist=TRUE)

print(paste0("============= estimateCorrection"))
	readCounts <- estimateCorrection(readCounts)

print(paste0("============= plot noise"))
	OutFile <- paste0(OutDir,"/",bamID,"_","noisePlot",BIN_S,".pdf")
	pdf(OutFile)
	noisePlot(readCounts)
	dev.off()
print(paste0("============= applyFilters correctBins normaizeBins smoothOutlierBins"))
	readCounts <- applyFilters(readCounts,residual=TRUE, blacklist=TRUE, chromosomes=NA)	# includes both X & Y
	copyNumbers <- correctBins(readCounts)

	copyNumbersNormalized <- normalizeBins(copyNumbers)
	copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)
print(paste0("============= copyNumbersSmooth"))
     	OutFile <- paste0(OutDir,"/",bamID,"_",BIN_S,"_copyNumbersSmooth.bed")
     	exportBins(copyNumbersSmooth,file = OutFile)
     	OutFile <- paste0(OutDir,"/",bamID,"_",BIN_S,"_copyNumbersSmooth.igv")
     	exportBins(copyNumbersSmooth,file = OutFile, format="igv")

print(paste0("~~~~~~~~~~~~~11111111111111~~~~~~~~~~~~~~~"))
print(str(copyNumbersSmooth))
	if (length(InputFile) > 1) {	### compare to reference, but the reference has to be one file.
		copyNumbersSmooth <- compareToReference(copyNumbersSmooth,c(2,FALSE))
	}
print(paste0("~~~~~~~~~~~~222222222222222~~~~~~~~~~~~~~~~"))


print(paste0("============= segmentBins & normalizeSegmentedBins"))
	copyNumbersSegmented <- segmentBins(copyNumbersSmooth, transformFun="sqrt")
	copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)

print(paste0("============= exportBins"))
     	OutFile <- paste0(OutDir,"/",bamID,"_",BIN_S,"_seg.bed")
     	exportBins(copyNumbersSegmented,file = OutFile, format="bed", filter=T,type=c("segments"))
    	OutFile <- paste0(OutDir,"/",bamID,"_",BIN_S,"_seg.igv")
     	exportBins(copyNumbersSegmented,file = OutFile, format="igv", filter=T,type=c("segments"))


#print(paste0("============= output copyNumbersSegmented into rds file"))
    	OutFile <- paste0(OutDir,"/",bamID,"_",BIN_S,"before.rds")
	saveRDS(copyNumbersSegmented, file=OutFile)

print(paste0("============= plot genome-wide ============="))
	OutFile <- paste0(OutDir,"/",bamID,"_",BIN_S,".pdf")
	print(OutFile)
	pdf(OutFile,paper='a4r')
        plot(copyNumbersSegmented, plot.type = "w", main = paste0(unlist(strsplit(bamID, "-"))[1]," ",unlist(strsplit(bamID, "-"))[2],",bin=",BIN_S)) 
	dev.off()

	OutFile <- paste0(OutDir,"/",bamID,"_",BIN_S,"header.csv")
	print(OutFile)
	write.table(copyNumbersSegmented@phenoData@data, file = OutFile, row.names = T, sep = "\t", quote = FALSE, append = FALSE)


	cmd <- paste0("/uz/data/hydra/shared_app/apps/imagemagick/7.0.7-27/bin/convert -flatten -density 250 -quality 100 ", paste0(OutDir,"/",bamID,"_",BIN_S,".pdf"), " ",   paste0(OutDir,"/",bamID,"_",BIN_S,".jpeg") )
	print(cmd)
	system(cmd)

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

#save(copyNumbersSegmented, file = paste(outPath, "", Family, ".rda", sep = ""))
	#print(paste0("============= set cutoff & plot ============="))
	#	copyNumbersCalled <- try(callBins(copyNumbersSegmented),silent=T)
	#	if (inherits(copyNumbersCalled,'try-error')) {
	#		copyNumbersCalled <- try(callBins(copyNumbersSegmented,method="cutoff",cutoffs=c(-0.1, 0.1)),silent=T)
	#	}

	#	cgh <- makeCgh(copyNumbersCalled)
	#	regions <- CGHregions(cgh)
	#	OutFile <- paste0(OutDir,"/",bamID,"_",BIN_S,"_after_regions.bed")
	#	exportBins(regions,file = OutFile, format="bed", filter=T,type=c("calls"))

    #    OutFile <- paste0(OutDir,"/",bamID,"_",BIN_S,"copyNumbersSmooth.Rdata")
	#	save(copyNumbersSmooth, file = OutFile)

    #    OutFile <- paste0(OutDir,"/",bamID,"_",BIN_S,".Rdata")
	#	save.image( file = OutFile)

		#OutFile <- paste0(OutDir,"/",bamID,"_",BIN_S,".vcf")
		#exportBins(copyNumbersSmooth,file = OutFile, format="vcf", filter=T)


	#	OutFile <- paste0(OutDir,"/",bamID,"_","Calls",BIN_S,".pdf")
	#	pdf(OutFile)
	#	plot(copyNumbersCalled)
	#	dev.off()



	if (opt != "yes") {
	print(paste0("============= exportBins: Smooth"))
		copyNumbersSmoothOutFile<-paste0(OutDir,"/",bamID,"_",BIN_S,"_smooth.bed")
	 	exportBins(copyNumbersSmooth,file = copyNumbersSmoothOutFile, format="bed", filter=T)
	    #get table in CNV table format
		outCNVtableSmooth <- buildCNVtable(copyNumbersSmoothOutFile, cytobandfile, bamID, BIN_S)
		outn <- paste0(OutDir,"/",bamID,"_",BIN_S,"_Smooth_CNVs.txt")
		write.table(outCNVtableSmooth, file = outn, row.names = FALSE, sep = "\t", quote = FALSE, append = FALSE)

		outCNVtableSegmented <- buildCNVtable(copyNumbersSegmentedOutFile, cytobandfile, bamID, BIN_S)
		outn <- paste0(OutDir,"/",bamID,"_",BIN_S,"_Segmented_CNVs.txt")
		write.table(outCNVtableSegmented, file = outn, row.names = FALSE, sep = "\t", quote = FALSE, append = FALSE)
	}

}


print(paste0("============= run QDNAseq main ============="))
loop_samples(bins,paste(bamDir,bamFile,sep="/"),bin_size,bamID)



