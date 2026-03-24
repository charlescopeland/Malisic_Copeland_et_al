#!/usr/bin/Rscript
#author: Tak Lee
#modified by Charles Copeland
#usage: Rscript DESeq2_wsva_PCA_diffs.R
#will produce PCA plots and do the DE analysis after the sva (batch correction). However you might want to do the DE analysis once the number of the surrogate variable from sva is fixed
# remove the if(FALSE){} statement to run the DE analysis
library(sva)
library(DESeq2)
library(ggplot2)
library(tidyverse)
library(limma)
library(foreach)
library(stringr)
library(readxl)
#setwd('\\\\fs-bio.mpipz.mpg.de\\netscratch\\dep_psl\\grp_psl\\CharlesNet\\VAPE15-4')
setwd('N:/dep_psl/grp_psl/CharlesNet/MIRO02-01')

suffix = "proc"
sampleGrouping = "condition"
interaction = length(str_split(sampleGrouping, " \\* ")[[1]]) > 1
figureType = ".pdf"
plotAllSVs = F
overwrite = F


#countFile <- paste("N:/dep_psl/grp_psl/CharlesNet/VAPE15-14/featurecounts/featureCounts_paper.tsv", sep="")
countFile <- "N:/dep_psl/grp_psl/CharlesNet/MIRO02-01/featurecounts/featureCounts_nofrax.tsv"

full <- read.delim(countFile, header = T,  comment.char = '#')

cts <- full[3:nrow(full), 7:(ncol(full))]
rownames(cts) <- full[3:nrow(full),1]

cts <- mutate_all(cts, function(x) as.integer(as.character(x)))

sampleList <- read_excel("MIRO02-01_samples.xlsx", sheet = "samples")
mySamples <- subset(sampleList, Sample %in% colnames(cts))

condition <- factor(mySamples[["Condition"]], levels = c("IT466_mock", "IT466_sid", "IT466_frax", 
                                                "tnrega_mock", "tnrega_sid", "drega_mock",
                                                "drega_sid"))
genotype <- factor(mySamples[["Genotype"]],levels = c("IT466", "tnrega", "drega"))
treatment <- factor(mySamples[["Treatment"]], levels=c("mock", "sid", "frax"))
sampleNames = mySamples[["Sample"]]
sampleInfo <- data.frame(sampleNames = sampleNames, condition = condition, genotype = genotype, treatment = treatment)
rownames(sampleInfo) <- sampleInfo$sampleNames
sampleInfo["type"] <- "paired-end"

if (setequal(rownames(sampleInfo), colnames(cts))) {
  sampleInfo_reordered <- sampleInfo[match(colnames(cts), rownames(sampleInfo)), , drop = FALSE]
  sampleInfo <- sampleInfo_reordered
}
print ("Row and Col positions of NA values")
which(is.na(cts), arr.ind=TRUE)
cts <- mutate_all(cts, ~coalesce(.,0))
initialDesign <- as.formula(paste("~", sampleGrouping))

ddsMat <- DESeqDataSetFromMatrix(cts, sampleInfo, design = initialDesign)
#keep <- rowSums(counts(ddsMat) >= 8) >= 3
#dds <- dds[keep, ]

dds <- estimateSizeFactors(ddsMat)
dat  <- counts(dds, normalized = TRUE)

if (overwrite) {
  write.table(dat, paste("deseq/", suffix, "_normalized_counts.tsv", sep=""))
}

###############filtering the lowly expressed genes################

geneidx  <- rowMeans(dat) > 1
dat <- dat[geneidx,]

###############finding the surrogate variables################
mod  <- model.matrix(initialDesign, colData(dds))
mod0 <- model.matrix(~ 1,colData(dds))
#manual detection of surrogate variable
#svseq <- svaseq(dat, mod, mod0, n.sv = 6)

#initializing with dummy data
dmicnt <- rep(0,50)
#auto detection of surrogate variables, this is repeated 5 times since it is a heuristic process. The best surrogate variable is taken from 5 iterations.
for (i in 0:5){
	svobj <- sva(dat, mod, mod0)
	dmicnt[svobj$n.sv] <- dmicnt[svobj$n.sv] +1
}
svnum <- which.max(dmicnt)
print(svnum)

plotTheme = theme(panel.grid = element_blank(),
                  panel.background = element_blank(),
                  plot.background = element_blank(), 
                  legend.background = element_blank())

if (plotAllSVs) {
  lowestSV = 0
} else {
  lowestSV = svnum
}
#plotting PCA plots for all the number of surrogate variables. Observe the PCA plot and find the variable number that groups your replicates together but does not group other samples that were orignially separated. This will usually lead to the best surrogate number from above.
for (k in unique(c(0, lowestSV:svnum))){
	print(k)
	svseq <- svaseq(dat, mod, mod0, n.sv = k)
	ddssva <- dds 
	fmla <- "~" 
	if (k > 0) {
	  for (j in 1:k){
		  newcol <- svseq$sv[,j]
		  colData(ddssva) <- cbind(colData(ddssva), newcol)
		  names(colData(ddssva))[j+length(sampleInfo)] <- paste0("SV",j)
		  if (j == 1){fmla=paste0(fmla,paste0("SV",j))} else {
			  fmla=paste(fmla,paste0("SV",j),sep="+")}}
		  fmla=paste(fmla, sampleGrouping, sep="+")
	  
	
	} else {
	  fmla = initialDesign
	}
	print(fmla)
	design(ddssva) <- as.formula(fmla)
	samples <- paste(colData(ddssva)[ ,"condition"], sep = "-")#[ , sampleGrouping]
	design <- model.matrix(~samples)

	################plotting PCA#################
	#vst normalization of data
	#vsd <- vst(ddssva)[geneidx,]
	vsd <- varianceStabilizingTransformation(dds, fitType = "local")[geneidx,]
	if (k > 0) {
	  rldData <- assay(vsd) %>%
		    removeBatchEffect(covariates = svseq$sv, design = design)
	} else {
	  rldData <- assay(vsd)
	}
	sampleIdx <- 1:length(samples)
	#two different treatments, group1 and group2, Use this wehn you have two (or more) different groups of data and would like to plot the PCAs separately
	#group1idx <- sampleIdx[colnames(rldData) %in% c("R179wt-HKR179wthiOD_rep1","R179wt-HKR179wthiOD_rep2","R179wt-HKR179wthiOD_rep3","R179wt-HKR179wthiOD_rep4","R179dssAB-HKR179wthiOD_rep1","R179dssAB-HKR179wthiOD_rep2","R179dssAB-HKR179wthiOD_rep3","axenic-HKR179wthiOD_rep1","axenic-HKR179wthiOD_rep2","axenic-HKR179wthiOD_rep3","axenic-HKR179wthiOD_rep4","axenic-MgSO4_rep1","axenic-MgSO4_rep2","axenic-MgSO4_rep3","axenic-MgSO4_rep4","axenic-R179wthiOD_rep1","axenic-R179wthiOD_rep2","axenic-R179wthiOD_rep3","axenic-R179wthiOD_rep4","axenic-R179dssABhiOD_rep1","axenic-R179dssABhiOD_rep2","axenic-R179dssABhiOD_rep3","axenic-R179dssABhiOD_rep4")]

	#group2idx <- sampleIdx[!colnames(rldData) %in% c("R179wt-HKR179wthiOD_rep1","R179wt-HKR179wthiOD_rep2","R179wt-HKR179wthiOD_rep3","R179wt-HKR179wthiOD_rep4","R179dssAB-HKR179wthiOD_rep1","R179dssAB-HKR179wthiOD_rep2","R179dssAB-HKR179wthiOD_rep3","axenic-HKR179wthiOD_rep1","axenic-HKR179wthiOD_rep2","axenic-HKR179wthiOD_rep3","axenic-HKR179wthiOD_rep4","axenic-MgSO4_rep1","axenic-MgSO4_rep2","axenic-MgSO4_rep3","axenic-MgSO4_rep4","axenic-R179wthiOD_rep1","axenic-R179wthiOD_rep2","axenic-R179wthiOD_rep3","axenic-R179wthiOD_rep4","axenic-R179dssABhiOD_rep1","axenic-R179dssABhiOD_rep2","axenic-R179dssABhiOD_rep3","axenic-R179dssABhiOD_rep4")]
	#only when you have subsets
	#group <- c("all","subset1","subset2")
	group <- c("")
	#PCA just for the top 500 genes
	top <- 500
	allgn <- length(rownames(rldData))
	allidx <- list(sampleIdx)#,group1idx,group2idx)
	for(i in 1:length(allidx)){
		idx <- allidx[[i]]
		pca <- prcomp(t(rldData[, idx]))
		length(rownames(rldData))
		for (n in c(top,allgn)){
		#for top 500 genes that contribute the most to the variance of PC1 and PC2.
			p1 <- order(abs(pca$rotation[,1]),decreasing=TRUE)[1:n]
			p2 <- order(abs(pca$rotation[,2]),decreasing=TRUE)[1:n]
			p3 <- order(abs(pca$rotation[,3]),decreasing=TRUE)[1:n]
			g1 <- names(pca$rotation[,1][p1])
			g2 <- names(pca$rotation[,2][p2])
			g3 <- names(pca$rotation[,3][p3])
			p1p2 <- union(g1,g2)
			p2p3 <- union(g1,g3)
			pca2 <- prcomp(t(rldData[p1p2, idx]))
			pca3 <- prcomp(t(rldData[p2p3, idx]))
			pc1 <- pca2$x[,1]
			pc2 <- pca2$x[,2]
			pc23 <- pca3$x[,1]
			pc33 <- pca3$x[,2]
			percentVar <- pca2$sdev^2/sum(pca2$sdev^2)
			percentVar <- round(100 * percentVar)

			pcaData <- data.frame(PC1 = pc1, PC2 = pc2, Genotype = colData(dds)[idx, 2], Bacteria = colData(dds)[idx, 4], ID = rownames(colData(dds))[idx])
			ggplot(pcaData, aes(x = PC1, y = PC2, colour = Genotype, shape = Bacteria)) + #, label = ID)) +
				geom_point(size = 4) +
				xlab(paste0("PC1: ",percentVar[1],"% variance")) +
				ylab(paste0("PC2: ",percentVar[2],"% variance")) +
				stat_ellipse(aes(x = PC1, y = PC2, group = interaction(Genotype, Bacteria)), type = 't', linetype = 2, level = 0.8) +
				coord_fixed(1) +
				theme_classic() +
				geom_text(aes(label=ID),vjust=2, size = 3) +
				theme(plot.title = element_text(hjust = 0.5, size = 12, face = 'bold'),
				      legend.text.align = 0,
				      axis.text = element_text(size = 18),
				      axis.title = element_text(size = 18),
				      legend.text=element_text(size= 18),
				      legend.title = element_text(size = 18),
				      panel.grid = element_blank(),
				      panel.background = element_blank(),
				      plot.background = element_blank(), 
				      legend.background = element_blank())
			if (overwrite) {
  			if(n == allgn){
  				outpdf <- paste("deseq/", suffix, group[i],"_PCA_sva",k, figureType, sep="")
  				ggsave(outpdf, width = 30,height=30, units = "cm")
  			} else {
  				outpdf <- paste("deseq/", suffix, group[i],"_PCA_top",n,"_sva",k, figureType, sep="")
  				ggsave(outpdf, width = 30,height=30, units = "cm")
  			}
  			#outjpg <- paste(group[i],"_PCA_top500.jpg",sep="")
  			#ggsave(outjpg, width = 30,height=30)
  			}

		}
	}
}

##################Differential Expression of genes###################
if( !interaction){
ddssva <- DESeq(ddssva[geneidx, ])
#modify your contrast list. The first element is the control (denominator)
useSV = F #the surrogate variable did not seem to make a difference so it was not used
if (!useSV) {
    design(ddssva) <- as.formula(paste("~", sampleGrouping))
}

simpleContrasts <- list(c("IT466_mock", "tnrega_mock"),
                        c("IT466_mock", "drega_mock"),
                        c("tnrega_mock", "drega_mock"),
                        c("IT466_mock", "IT466_sid"),
                        c("IT466_mock", "IT466_frax"),
                        c("tnrega_mock", "tnrega_sid"),
                        c("drega_mock", "drega_sid"),
                        c("IT466_sid", "tnrega_sid"),
                        c("IT466_sid", "drega_sid"),
                        c("tnrega_sid", "drega_sid"))
simpleContrasts <- list(c("IT466_sid", "tnrega_mock"),
                        c("IT466_sid", "drega_mock"))
combinedContrasts <- list(c("Col_mock_12", "Col_R569_12"), 
                         c("Col_mock_24", "Col_R569_24"),
                         c("Van_mock_12", "Van_R569_12"),
                         c("Van_mock_24", "Van_R569_24"),
                         c("Col_mock_120", "Col_R569_120"),
                         c("Van_mock_120", "Van_R569_120"),
                         c("Van_RW_24", "Van_R569_24"),
                         c("Col_RW_24", "Col_R569_24"),
                         c("Van_R569_12","Van_R569_48"),
                         c("Col_R569_12", "Col_R569_48"))


contrasts <- simpleContrasts

for(i in 1:length(contrasts)){
  
  
  ctrl <- contrasts[[i]][1]
	trt <- contrasts[[i]][2]
	res <- results(ddssva, contrast=c(sampleGrouping ,trt,ctrl))
	outf <- paste("deseq/", ctrl,".VS.",trt,"_diff_", suffix, ".csv",sep="")
	write.csv(res,outf)
	#print(outf)
	}
} else if (interaction) {
  ddssva <- DESeq(ddssva[geneidx, ])
  design(ddssva) <- ~ genotype * treatment
  ddssva <- DESeq(ddssva)
  
  #contraster function from https://www.r-bloggers.com/2024/05/a-guide-to-designs-and-contrasts-in-deseq2/
  contraster <- function(dds,    # should contain colData and design
                       group1, # list of character vectors each with 2 or more items 
                       group2, # list of character vectors each with 2 or more items
                       weighted = F
  ){
  
  
  mod_mat <- model.matrix(design(dds), colData(dds))
  
  grp1_rows <- list()
  grp2_rows <- list()
  
  
  for(i in 1:length(group1)){
    
    grp1_rows[[i]] <- colData(dds)[[group1[[i]][1]]] %in% group1[[i]][2:length(group1[[i]])]
    
  }
  
  
  for(i in 1:length(group2)){
    
    grp2_rows[[i]] <- colData(dds)[[group2[[i]][1]]] %in% group2[[i]][2:length(group2[[i]])]
    
  }
  
  grp1_rows <- Reduce(function(x, y) x & y, grp1_rows)
  grp2_rows <- Reduce(function(x, y) x & y, grp2_rows)
  
  mod_mat1 <- mod_mat[grp1_rows, ,drop=F]
  mod_mat2 <- mod_mat[grp2_rows, ,drop=F]
  
  if(!weighted){
    
    mod_mat1 <- mod_mat1[!duplicated(mod_mat1),,drop=F]
    mod_mat2 <- mod_mat2[!duplicated(mod_mat2),,drop=F]
    
  }
  
  return(colMeans(mod_mat1)-colMeans(mod_mat2))
  
  
}
  
  outf <- paste("deseq/", suffix, "_", paste(str_split(sampleGrouping, " \\* ")[[1]], collapse = "_"), "_.VS.sid_drega.csv",sep="")
  
  diff1 <-  contraster(ddssva,
                       group1 = list(c("genotype", "IT466"), c("treatment", "sid")),
                       group2 = list(c("genotype", "drega"), c("treatment", "sid")))

  diff2 <-  contraster(ddssva,
                       group1 = list(c("genotype", "IT466"), c("treatment", "mock")),
                       group2 = list(c("genotype", "drega"), c("treatment", "mock")))
  
  res1 <- results(ddssva, 
                  contrast = diff1-diff2)
  write.csv(res1,outf)
  
  
}



