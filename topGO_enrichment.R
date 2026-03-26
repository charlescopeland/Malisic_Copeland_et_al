#### Begin ####
#!/usr/bin/Rscript
#author: Tak Lee
#using topGO for structured GO term enrichment analysis
#usage: Rscript topGO.R filelist_with_list_of_genes
library(topGO)
library(ggplot2)
library(dplyr)
  #input arabidopsis go term


if (Sys.info()["nodename"] == "DAHLIA") {
  setwd("\\\\fs-bio.mpipz.mpg.de\\biodata\\dep_psl\\grp_psl\\Charles")
} else {
  setwd("N:\\dep_psl\\grp_psl\\CharlesNet\\MIRO02-01\\functions")
}
g2go <- readMappings(file = "N:\\dep_psl\\grp_psl\\CharlesNet\\MIRO02-01\\functions\\IT466_gene_to_go_shared.map")
geneNames <- names(g2go)
go2g <- inverseList(g2go)
#### Data ####

#fileName <- "N:\\dep_psl\\grp_psl\\CharlesNet\\VAPE16\\deseq\\48_interaction_paper_genotype_bacteria_interaction.csv"
#geneColumn = "Gene"
#prefix <- "48_interaction_paper_down"



letterSize <- 10
combinedTheme = theme(axis.line = element_line(colour = "Black", linewidth = 0.5),
                      panel.grid = element_blank(),
                      panel.background = element_blank(), #element_rect(fill = "transparent"),
                      plot.background = element_blank(), #element_rect(fill = "transparent", color = NA),
                      text = element_text(size = letterSize, colour = "black"),
                      axis.text = element_text(size = letterSize, colour = "black"),
                      axis.text.y = element_text(vjust = 0.3),
                      axis.text.x = element_text(vjust = 0.5, hjust = 0.5),
                      legend.background = element_blank(), #element_rect(fill = "transparent", color = NA),
                      #legend.title = element_blank(),
                      legend.text = element_text(size = letterSize ),
                      strip.background = element_blank(),
                      strip.text = element_text(angle = 90))

#input of list of files with genes as a list

#for each genelist, do GO term enrichment analysis
#genes <- read.csv(fileName)
#genes <- subset(genes, padj < 0.05 & log2FoldChange < 0)
#genes <- genes[[geneColumn]]
#geneList <- factor(as.integer(geneNames %in% genes))
#names(geneList) <- geneNames
fileList = read.csv("N:\\dep_psl\\grp_psl\\CharlesNet\\MIRO02-01\\kmeans\\MIRO02-01_all_sig_genes_hierarch_deseq_normalized_k7_cluster_genes.fof", header = F)
for (file in fileList$V1) {
  geneList = read.csv(file, header = F)
  geneList <- factor(as.integer(geneNames %in% geneList[["V1"]]))
  names(geneList) <- geneNames
  if (length(levels(geneList)) > 1){
    mergedat = data.frame()
  	for (categ in c("BP", "MF", "CC")){
  	  figName = paste(file, "_shared_GO_enr_",categ,".pdf",sep="")
  	  outfilen <- paste(file,"_shared_GO_enr_",categ,".txt",sep="")
  	  
  	  GOdata <- new("topGOdata",ontology=categ,allGenes=geneList,annot = annFUN.gene2GO, gene2GO = g2go)
  		testout <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
  		
  		allRes <- GenTable(GOdata, Fis = testout,topNodes=100)
  		#newAllRes <- allRes %>% mutate(Fis = ifelse(Fis < 1e-30, 1e-30, as.numeric(Fis)))
  		
  		allRes$Fis <- as.numeric(allRes$Fis)
  		allRes$Fis[is.na(allRes$Fis)] <- 1e-30
  		pfilt <- allRes[allRes$Fis < 0.01,]
  		names(pfilt) <- c("GOid","GOterm","Annotated","NgenesInTerm","Expected","SignificanceP")
  		pfilt$GenesInTermRatio <- pfilt$NgenesInTerm/pfilt$Annotated
  		pfilt$Enrichment <- pfilt$NgenesInTerm/pfilt$Expected
  		
  		allGO <- genesInTerm(GOdata)
  		genesinTerm <- lapply(allGO,function(x) x[geneList[x] == 1] )
  		annotdegs <- c()
  		for (gg in pfilt$GOid){
  			 annotg <- unlist(genesinTerm[gg])
  			 ll <- length(annotg)
  			 annotgenes <- paste(annotg,collapse=",")
  			 annotdegs <- append(annotdegs,annotgenes)
  		}
  		pfilt$GOannotdegs <- annotdegs
  	  #write out information table
  		write.table(pfilt,file=outfilen,sep="\t",quote=F,row.names=F)
  		title = ''
  		if (categ == "CC") { title = "Cellular_Component"} else if (categ == "BP") { title = "Biological_Process"} else { title = "Molecular_Function" }
  		m = (min(-log10(pfilt$SignificanceP)) + max(-log10(pfilt$SignificanceP)))/2
  		pfilt$GOs <- paste(pfilt$GOid,pfilt$GOterm)
  		pfilt$GOs <- factor(pfilt$GOs, levels = pfilt$GOs[order(pfilt$GenesInTermRatio)])
  		#plot enrichment dot plots
  		
  		gofig <- ggplot(pfilt, aes(x = log2(Enrichment), y = GOs, color = -log10(SignificanceP), size = NgenesInTerm)) +
  		      geom_point() +
  		      scale_color_gradient(low = "#fdaeff", high = "#8800aa")+ #scale_color_gradient2(low="#eeaa00",mid="#ee5500",high="#aa0000",midpoint=m) +
  		      xlab("log2 enrichment in number of genes vs expected") +
  		      ylab("GOterms") +
  		      ggtitle(paste(title,"enrichments")) +
  		      scale_size_continuous(range = c(1, 8)) +
  		      combinedTheme +
  		      theme(panel.grid.major.y = element_line(colour = "lightgrey", linetype = 2))
  		gofig
  		
  		ggsave(figName, 
  		       plot = gofig , 
  		       width = (27),
  		       units = "cm")
  		    
  		    #pfilt$categoryy <- rep(categ,dim(pfilt)[1])
  		    #mergedat <- rbind(mergedat,pfilt)
  
  	}
  }
}
	#plot enrichment of all "BP","CC","MF" terms in a single plot
	    #mergedpdf = gsub(".txt","_GOenrichment_all.pdf",f)
	    #mergedat$GOs <- factor(mergedat$GOs, levels = mergedat$GOs[order(mergedat$categoryy)])
	    #pdf(file=mergedpdf,width=16,height=40)
	    #plot <- ggplot(mergedat, aes(x = categoryy, y = GOs, color = -log10(SignificanceP), size = log10(NgenesInTerm))) +
	    #  geom_point() +
	    #  scale_color_gradient2(low="#ffcc00",mid="#ff6600",high="#e60000",midpoint=m) +
	    #  xlab("GOterm category") +
	    #  ylab("GOterms") +
	    #  ggtitle("All GO term category enrichments") +
	    #  scale_size_continuous(range = c(1, 10)) +
	    #  theme_classic()
	    #print(plot)
	    #dev.off()
	    #mfilen = gsub(".txt","_GOenrichment_all.txt",f)
	    #mergedat$GOs <- NULL
	    #write.table(mergedat,file=mfilen,sep="\t",quote=F,row.names=F)
	

