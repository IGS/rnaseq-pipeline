#! /usr/local/packages/R-3.3.1/bin/Rscript --slave --vanilla

#
# Usage: deseq2.R $sample_info $out_dir 
#
# tcreasy@som.umaryland.edu 
# modified by thodges@som.umaryland.edu to use DESeq2v1.18 instead of DESeq 
#modified by Apaala Chatterjee to use in pipeline
###############################################################################

##
## sample_info: absolute path of sample info file. Sample info file format: sample_name\tphenotype\tabsolute path of read count file, no header.
## out_dir: absolute path of output directory
## annotation_file: optional. Absolute path of annotation file, first column of annotation file must be gene ID used for read counting, tab-delimited.
## 
##

rm(list=ls(all=T))
cat("\n***** Starting DESeq2 (v) analysis ******\n\n")

# input arguments
args <- commandArgs(TRUE)
print(length(args))
stopifnot(length(args) > 1 )

# sample comparison input file
sample.info <- read.delim(args[1], header=FALSE, sep="\t")

# output directory
output.dir <- args[2]

# set the working directory
setwd(output.dir)

## combine single file into a data matrix and output counting statistics
data.tab <- read.delim( as.character(sample.info[1,3]), header=F, sep="\t" )
#idx <- which(data.tab=="__no_feature")
idx <- grep("__no_feature", data.tab$V1)[1]
count.stat <- matrix(c("sum counted", "total", sum(data.tab[,2][1:(idx-1)]), sum(data.tab[,2])), nrow=2, ncol=2)
count.stat <- rbind(data.tab[idx:nrow(data.tab),], count.stat)
colnames(count.stat) <- c("Stat", as.character(sample.info[1,1]))
data.tab <- data.tab[1:(idx-1),]

for (i in 2:nrow(sample.info)) {
	d <- read.delim(as.character(sample.info[i,3]), header=F, sep="\t")
	tmp <- matrix(c(sum(d[,2][1:(idx-1)]), sum(d[,2])), nrow=2, ncol=1)
	tmp <- data.frame(x=c(d[idx:nrow(d),2], tmp))
	n <- colnames(count.stat)
	count.stat <- cbind(count.stat, tmp)
	colnames(count.stat) <-c (n, as.character(sample.info[i,1]))
	d <- d[1:(idx-1),]
	data.tab <- merge(data.tab, d, by=1)
}


colnames(data.tab) <- c("ID", as.character(sample.info[,1]))
write.table(data.tab, file.path("all_counts"), na="", quote=F, row.names=F, sep="\t")
write.table(count.stat, file.path("count_stat"), na="", quote=F, row.names=F, sep="\t")


# remove genes without reads for all samples
data.tab <- data.tab[which(rowSums(data.tab[,2:ncol(data.tab)])>0),]
write.table(data.tab, file.path("all_counts_noZero"), na="", quote=F, row.names=F, sep="\t")
cat(paste("\n* There are ", ncol(data.tab)-1, " samples and ", nrow(data.tab), " genes used for DE analysis\n", sep=""))


# import DESeq2 and gplots
suppressMessages( library("DESeq2") )
suppressMessages( library("RColorBrewer") )
suppressMessages( library("gplots") )

# phenotypes to compare
pheno <- unique(as.character(sample.info[,2]))
cat(paste("\n* Phenotypes found: ", toString(pheno), "\n", sep=""))
d <- data.tab[,2:ncol(data.tab)]
rownames(d) <- data.tab$ID

# condition information from sample info file
condition <- factor(as.character(sample.info[,2]))
#condition=sort(condition, decreasing=T)
print("Condition printed")
print(condition)
# setup colData for creating DESeqData object
colData <- data.frame(row.names = colnames(d), condition)

# Create DESeqDataSetFromMatrix object 
dds = DESeqDataSetFromMatrix(
	countData = d,
	colData = colData, 
	design = ~ condition)
	
## DE analysis - DESeq() runs DESeq2 algorithm:

		## estimating size factors
		## estimating dispersions
		## gene-wise dispersion estimates
		## mean-dispersion relationship
		## final dispersion estimates
		## fitting model and testing

		dds = DESeq(dds)
		
		# To see Size factors for normalization
		colData(dds)		
	
		#info about the matrix and design
		head(dds)
		resultsNames(dds)


# create an output file name for the output PDF
pdf.name <- paste(pheno[1], "-", pheno[2], ".pdf", sep="")
pdf(pdf.name)


# variance testing
cat("\n* Estimating variance...\n")
# rlog transformation for clustering
rld <- rlog(dds)

# color palette for plots
hmcol = colorRampPalette(brewer.pal(9, "RdBu"))(100)

# output to tab file
out <- cbind(rownames(assay(rld)), assay(rld))
colnames(out) <- c("ID", colnames(assay(rld)))
write.table(out, file.path("all_counts_noZero_normalized"), na="", sep="\t", quote=F, row.names=F)

# Heatmap showing clustering of samples
dists = dist( t( assay(rld) ) )
mat = as.matrix( dists )
rownames(mat) = colnames(mat)
dist.title <- paste("Sample Clustering", sep="")
heatmap.2(mat, col=rev(hmcol), trace="none", main=dist.title, margin=c(13,13), cexRow=0.8, cexCol=0.8, keysize=1.0)


for (k in 1:(length(pheno)-1)) {
	
	for (m in (k+1):length(pheno)) {
		
	cat(paste("\n* Running DESeq2 algorithm for: ", pheno[k], " vs ", pheno[m], "\n", sep=""))
    #print(type(k))				
	# To get results
	res <- results(dds, contrast=c("condition",pheno[k] ,pheno[m]),cooksCutoff=FALSE)

	#Results object
	#Note: I need to discuss further with Amol about the independentFiltering parameter
    res <- results(dds, independentFiltering=FALSE, contrast=c("condition",pheno[k] , pheno[m]),cooksCutoff=FALSE)
    	    
    #Information on results dataframe
    mcols(res, use.names =T)
    
    # order output by FDR
    res <- res[order(res$padj),]
    	
    cat("\n* Results Snippet: res\n")
    print(head(res))
    	
	#summarize some basic tallies using the summary function
    summary(res)	
       
    # plot the results using FDR=0.05 as the cutoff
    ma.title <- paste("DEG MA Plot", " (FDR < 0.05)", sep="")
    plotMA(res, main=ma.title, xlab="Mean of Normalized Counts", ylab=paste("LFC: ", pheno[k], " VS ", pheno[m], sep="")) 
    
    #Results no longer includes the Mean per gene for each condition so it needs to be added to the results
    #Extracted the counts (normalized by size factors) for two conditions of interest
    #From the counts slot of the dds object
    
    Read.Count.k <- rowMeans(counts(dds,normalized=TRUE)[,dds$condition == pheno[k]])
    Read.Count.m <- rowMeans(counts(dds,normalized=TRUE)[,dds$condition == pheno[m]])
    read.counts.both <- merge(Read.Count.m, Read.Count.k, by.x=0, by.y=0)
    rownames(read.counts.both) <- read.counts.both[,1]
    read.counts.both <- read.counts.both[,-1]
    colnames(read.counts.both) <- c("Read.Count.m", "Read.Count.k")
    
    #create dataframe of res oject
    res.df <- as.data.frame(res)
    #Merge with read.counts.both
    resdata <- merge(read.counts.both, res.df, by.x=0, by.y=0)
    resdata <- resdata[,c(1,4,2,3,5:9)]
    colnames(resdata) <- c("Feature.ID", "Read.Count.All", "Read.Count.m", "Read.Count.k",  "log2FoldChange","lfcSE", "stat", "p.Value", "padj")	
    #Note:  Check with Amol to see if you need to: 1. add a FC column and 2. remove or keep the lfcSE, and stat columns
    
    # get read counts for each group for the top 30 most significant DEGs
    # order output by absolute LFC
    resdata <- resdata[order(-abs(resdata$log2FoldChange)),]
    sig.genes = resdata[!is.na(resdata$padj<=0.05),]
    print(dim(sig.genes))
    print("*********************")
    print(head(sig.genes))
    if(nrow(sig.genes) < 2) {
      sig.genes <- resdata[resdata$pval<=0.05,]
      print(dim(sig.genes))
    }
    
    if(nrow(sig.genes) > 30) {
      sig.genes <- sig.genes[1:30,]
      print(dim(sig.genes))
    }
    
    cat("\n* Results Snippet: sig.genes\n")
    print(head(sig.genes))
    
		
    #Prepare a matrix of sig.genes for heatmap	
    read.counts.sig <- cbind(as.numeric(sig.genes[,3]), as.numeric(sig.genes[,4]))
    colnames(read.counts.sig) <- c(pheno[m], pheno[k])
    rownames(read.counts.sig) <- c(sig.genes[,1])
    
    cat("\n* Results Snippet: read.counts.1\n")
    print(head(read.counts.sig))
    
    # draw heatmap of normalized read counts for the significant genes of each sample
    sig.title <- paste("Top Significant DEGs", " (per condition)", sep="")		
    heatmap.2(read.counts.sig, col=hmcol, trace="none", main=sig.title, margin=c(13,13), cexRow=0.8, cexCol=0.8, keysize=1.0)
    
    read.counts.sig <- sig.genes[,c(1,5)]
    colnames(read.counts.sig) <- c("ID", "LFC")
    read.counts.sig <- merge(read.counts.sig, out, by.x=1, by.y=1)
    read.counts.sig <- read.counts.sig[order(-abs(read.counts.sig$LFC)),]
    
    cat("\n* Results Snippet: read.counts.2\n")
    print(head(read.counts.sig))
    #####Changed order of k and m to indicate deseq2 calculations of condition 2 over 1
    write.table(read.counts.sig, file.path(paste(pheno[k], "_vs_", pheno[m], ".top30.counts.txt", sep="")), na="", quote=F, row.names=F, sep="\t")
    
    hmap <- read.delim(file.path(paste(pheno[k], "_vs_", pheno[m], ".top30.counts.txt", sep="")), header=T, sep="\t" )
    
    
    cat("\n* Results Snippet: hmap.1\n")
    print(head(hmap))
    
    hmap <- hmap[,c(3:ncol(hmap))]
    colnames(hmap) <- c(colnames(read.counts.sig)[3:ncol(read.counts.sig)])
    rownames(hmap) <- c(read.counts.sig[,1])
    hmap <- data.matrix(hmap)
    
    cat("\n* Results Snippet: hmap.2\n")
    print(head(hmap))
    
    # draw heatmap of normalized read counts for the significant genes of each sample
    sig.title <- paste("Top Significant DEGs",  " (per sample)", sep="\n")		
    heatmap.2(hmap, col=hmcol, trace="none", main=sig.title, margin=c(13,13), cexRow=0.8, cexCol=0.8, keysize=1.0)
    
    # Change column names for clarity and brevity
    # Deseq2 calculates LFC condition 2 vs 1. Renaming changed NOV 30
    colnames(resdata) <- c("Feature.ID", "Read.Count.All", paste("Read.Count.", pheno[m], sep=""), paste("Read.Count.", pheno[k], sep=""), paste("LFC(", pheno[k], "/", pheno[m], ")", sep=""),"lfcSE", "stat", "p.Value", "FDR")	
    #colnames(resdata) <- c("Feature.ID", "Read.Count.All", paste("Read.Count.", pheno[m], sep=""), paste("Read.Count.", pheno[k], sep=""), paste("LFC(", pheno[m], "/", pheno[k], ")", sep=""),"lfcSE", "stat", "p.Value", "FDR")
    # write data to tsv file. Replaced m with k to change name.
    write.table(resdata, file.path(paste(pheno[k], "_vs_", pheno[m], ".de_genes.txt", sep="")), na="", quote=F, row.names=F, sep="\t")

	}
}

data.tab <- NULL
d <- NULL
dev.off()

cat("\n\n* Garbage Collection Information\n\n")
gc()

cat("\n\n***** DEG Analysis Complete *****\n\n")






















