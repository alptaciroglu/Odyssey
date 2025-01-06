# for(i in 1:10){
#   rnGen <- runif(n=159,min=0,max=1)
#   kngpSave[[2]] <- kngpSave[[2]] + (kngpSave[[2]] * rnGen)
#   KNGP_results <- KNGPAlgorithm(kngpSave[[1]], kngpSave[[2]], kngpSave[[3]], kngpSave[[4]])
#   if(i == 1){
#     KNGPDT <- KNGP_results$KNGPResult[order(KNGP_results$KNGPResult$postProb, decreasing = TRUE),]
#   }
#   else{
#     tmpRes <- KNGP_results$KNGPResult[order(KNGP_results$KNGPResult$postProb, decreasing = TRUE),]
#     KNGPDT <- cbind(KNGPDT, tmpRes)
#   }
# }
# 
# 
# 
# gotableList <- list()
# GOTable.tmp <- unique(unlist(genetoGO[unique(saveNodeNames)]))
# 
# 
# for(i in 1:length(GOTable.tmp)){
#   checkRedun <- unique(unlist(gotoGene[GOTable.tmp[i]]))
#   gotableList[[i]] <- intersect(checkRedun, saveNodeNames)
#   names(gotableList)[i] <- GOTable.tmp[i]
#   if(length(intersect(checkRedun, saveNodeNames)) == length(saveNodeNames)){
#     GOTable.tmp <- GOTable.tmp[-i]
#     print('Redundant GO Term deleted...')
#   }
# }
# 
# mart = useMart("ensembl")
# mart=useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
# nodeNames <- 'PMAIP1'
# biomaRtGO <- biomartr::getGO(organism = "Homo sapiens", 
#                                           genes    = nodeNames,
#                                           filters  = "hgnc_symbol")
# biomaRtGO <- as.data.table(biomaRtGO)
# 
# 
# 
# 
# #######################################################################
# library("biomaRt")
# 
# martObj <- useMart("ensembl")
# martObj <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
# listAttributes(mart = martObj) # Choose data types you want to download
# 
# go <- getBM(attributes=c("hgnc_id", "go_id", "namespace_1003"), mart = martObj)
# 
# go <- go[go[,3]!="",]; go[,3] <- as.character(go[,3])
# write.table(go, "GOannotationsBiomart_mod.txt", quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
# 
# ### Incomplete, downloaded from Biomart instead...
# 
# catdb <- makeCATdb(myfile="GOannotationsBiomart_mod.txt", lib=NULL, org="", colno=c(1,2,3), idconv=NULL)
# catdb
# 
# ## Create catDB from Bioconductor annotation package
# # catdb <- makeCATdb(myfile=NULL, lib="ath1121501.db", org="", colno=c(1,2,3), idconv=NULL)
# 
# ## AffyID-to-GeneID mappings when working with AffyIDs 
# # affy2locusDF <- systemPipeR:::.AffyID2GeneID(map = "ftp://ftp.arabidopsis.org/home/tair/Microarrays/Affymetrix/affy_ATH1_array_elements-2010-12-20.txt", download=TRUE)
# # catdb_conv <- makeCATdb(myfile="GOannotationsBiomart_mod.txt", lib=NULL, org="", colno=c(1,2,3), idconv=list(affy=affy2locusDF))
# # systemPipeR:::.AffyID2GeneID(catdb=catdb_conv, affyIDs=c("244901_at", "244902_at"))
# 
# ## Next time catDB can be loaded from file
# save(catdb, file="catdb.RData") 
# load("catdb.RData")
# 
# testOdyssey <- c('PMAIP1', 'TP53', 'ULK3')
# 
# goResultsMF <- GOHyperGAll(catdb=catdb, gocat="MF", sample=testOdyssey, Nannot=2)
# goResultsMF <- goResultsMF[goResultsMF$Padj < 0.05,]
# 
# goResultsBP <- GOHyperGAll(catdb=catdb, gocat="BP", sample=testOdyssey, Nannot=2)
# goResultsBP <- goResultsBP[goResultsBP$Padj < 0.05,]
# 
# goResultsCC <- GOHyperGAll(catdb=catdb, gocat="CC", sample=testOdyssey, Nannot=2)
# goResultsCC <- goResultsCC[goResultsCC$Padj < 0.05,]
# 
# goResultsComplete <- rbind(goResultsMF, goResultsBP, goResultsCC)
# goResultsComplete <- goResultsComplete[order(goResultsComplete$Padj),]

###########################
# aList <- list()
# 
# for(node in NodeNamesVectorServer){
#   
#   if(length(intersect(names(aList), node)) == 0){
#     
#     part1 <- unique(saveNetworkFile[saveNetworkFile$Gene == node,]$miRNA)
#     part2 <- unique(saveNetworkFile[saveNetworkFile$miRNA == node,]$Gene)
#     aList[[node]] <- c(part1, part2)
#     
#   }
#   
# }

############################

exampleDataList <- c("GSE35389_1", "GSE35389_2", "GSE88721", "GSE39061_1", "GSE39061_2", "GSE39061_3", "GSE49697_1", "GSE49697_2", "GSE25402", "GSE32539_1", "GSE32539_2", "GSE32539_3", "GSE34681_1", "GSE34681_2", "GSE38617", "GSE40321", "GSE59702_1", "GSE59702_2", "GSE81867", "GSE104268_1", "GSE104268_2", "GSE90604")
limmaListmiRNA <- list()
limmaListmRNA <- list()
load('miRNAexprsDB.RData')
load('mRNAexprsDB.RData')
load('miRNAphenoDB.RData')
load('mRNAphenoDB.RData')

for(i in 1:length(exampleDataList)){
  
  tmpDataCode <- unlist(strsplit(exampleDataList[i], "_"))[1]
  mrnaFile <- mRNAexprsDB[[tmpDataCode]]
  mirnaFile <- miRNAexprsDB[[tmpDataCode]]
  
  if(exampleDataList[i] == 'GSE25402'){
    
    ControlsmRNA <- mrnaFile[,grepl('Nonobese', mRNAphenoDB[['GSE25402']]$title, ignore.case = FALSE)]
    TreatmentsmRNA <- mrnaFile[,grepl('Obese', mRNAphenoDB[['GSE25402']]$title, ignore.case = FALSE)]
    ControlsmiRNA <- mirnaFile[,grepl('Nonobese', miRNAphenoDB[['GSE25402']]$title, ignore.case = FALSE)]
    TreatmentsmiRNA <- mirnaFile[,grepl('Obese', miRNAphenoDB[['GSE25402']]$title, ignore.case = FALSE)]
  }
  
  if(exampleDataList[i] == 'GSE32539_1'){
    
    ControlsmRNA <- mrnaFile[,grepl('control', mRNAphenoDB[['GSE32539']]$`final diagnosis:ch1`, ignore.case = FALSE)]
    TreatmentsmRNA <- mrnaFile[,grepl('IPF/UIP', mRNAphenoDB[['GSE32539']]$`final diagnosis:ch1`, ignore.case = FALSE)]
    ControlsmiRNA <- mirnaFile[,grepl('Control', miRNAphenoDB[['GSE32539']]$`final diagnosis:ch1`, ignore.case = FALSE)]
    TreatmentsmiRNA <- mirnaFile[,grepl('IPF/UIP', miRNAphenoDB[['GSE32539']]$`final diagnosis:ch1`, ignore.case = FALSE)]
  }
  
  if(exampleDataList[i] == 'GSE32539_2'){
    
    ControlsmRNA <- mrnaFile[,grepl('control', mRNAphenoDB[['GSE32539']]$`final diagnosis:ch1`, ignore.case = FALSE)]
    TreatmentsmRNA <- mrnaFile[,grepl('NSIP', mRNAphenoDB[['GSE32539']]$`final diagnosis:ch1`, ignore.case = FALSE)]
    ControlsmiRNA <- mirnaFile[,grepl('Control', miRNAphenoDB[['GSE32539']]$`final diagnosis:ch1`, ignore.case = FALSE)]
    TreatmentsmiRNA <- mirnaFile[,grepl('NSIP', miRNAphenoDB[['GSE32539']]$`final diagnosis:ch1`, ignore.case = FALSE)]
  }
  
  if(exampleDataList[i] == 'GSE32539_3'){
    
    ControlsmRNA <- mrnaFile[,grepl('control', mRNAphenoDB[['GSE32539']]$`final diagnosis:ch1`, ignore.case = FALSE)]
    TreatmentsmRNA <- mrnaFile[,grepl('RB-ILD', mRNAphenoDB[['GSE32539']]$`final diagnosis:ch1`, ignore.case = FALSE)]
    ControlsmiRNA <- mirnaFile[,grepl('Control', miRNAphenoDB[['GSE32539']]$`final diagnosis:ch1`, ignore.case = FALSE)]
    TreatmentsmiRNA <- mirnaFile[,grepl('RB-ILD', miRNAphenoDB[['GSE32539']]$`final diagnosis:ch1`, ignore.case = FALSE)]
  }
  
  if(exampleDataList[i] == 'GSE34681_1'){
    
    ControlsmRNA <- mrnaFile[,grepl('control siRNA', mRNAphenoDB[['GSE34681']]$`treatment:ch1`, ignore.case = FALSE)]
    TreatmentsmRNA <- mrnaFile[,grepl('Ars2 siRNA', mRNAphenoDB[['GSE34681']]$`treatment:ch1`, ignore.case = FALSE)]
    ControlsmiRNA <- mirnaFile[,grepl('control siRNA', miRNAphenoDB[['GSE34681']]$`treatment:ch1`, ignore.case = FALSE)]
    TreatmentsmiRNA <- mirnaFile[,grepl('Ars2 siRNA', miRNAphenoDB[['GSE34681']]$`treatment:ch1`, ignore.case = FALSE)]
  }
  
  if(exampleDataList[i] == 'GSE34681_2'){
    
    ControlsmRNA <- mrnaFile[,grepl('control siRNA', mRNAphenoDB[['GSE34681']]$`treatment:ch1`, ignore.case = FALSE)]
    TreatmentsmRNA <- mrnaFile[,grepl('DGCR8 siRNA', mRNAphenoDB[['GSE34681']]$`treatment:ch1`, ignore.case = FALSE)]
    ControlsmiRNA <- mirnaFile[,grepl('control siRNA', miRNAphenoDB[['GSE34681']]$`treatment:ch1`, ignore.case = FALSE)]
    TreatmentsmiRNA <- mirnaFile[,grepl('DGCR8 siRNA', miRNAphenoDB[['GSE34681']]$`treatment:ch1`, ignore.case = FALSE)]
  }
  
  if(exampleDataList[i] == 'GSE38617'){
    
    ControlsmRNA <- mrnaFile[,grepl('healthy', mRNAphenoDB[['GSE38617']]$`disease state:ch1`, ignore.case = FALSE)]
    TreatmentsmRNA <- mrnaFile[,grepl('oral lichen planus', mRNAphenoDB[['GSE38617']]$`disease state:ch1`, ignore.case = FALSE)]
    ControlsmiRNA <- mirnaFile[,grepl('healthy', miRNAphenoDB[['GSE38617']]$`disease state:ch1`, ignore.case = FALSE)]
    TreatmentsmiRNA <- mirnaFile[,grepl('oral lichen planus', miRNAphenoDB[['GSE38617']]$`disease state:ch1`, ignore.case = FALSE)]
  }
  
  if(exampleDataList[i] == 'GSE39061_1'){
    
    ControlsmRNA <- mrnaFile[,grepl('Confluent', mRNAphenoDB[['GSE39061']]$`culture stage:ch1`, ignore.case = FALSE)]
    TreatmentsmRNA <- mrnaFile[,grepl('Day 28', mRNAphenoDB[['GSE39061']]$`culture stage:ch1`, ignore.case = FALSE)]
    ControlsmiRNA <- mirnaFile[,grepl('Confluent', miRNAphenoDB[['GSE39061']]$`culture stage:ch1`, ignore.case = FALSE)]
    TreatmentsmiRNA <- mirnaFile[,grepl('Day 28', miRNAphenoDB[['GSE39061']]$`culture stage:ch1`, ignore.case = FALSE)]
  }
  
  if(exampleDataList[i] == 'GSE39061_2'){
    
    ControlsmRNA <- mrnaFile[,grepl('Subconfluent', mRNAphenoDB[['GSE39061']]$`culture stage:ch1`, ignore.case = FALSE)]
    TreatmentsmRNA <- mrnaFile[,grepl('Confluent', mRNAphenoDB[['GSE39061']]$`culture stage:ch1`, ignore.case = FALSE)]
    ControlsmiRNA <- mirnaFile[,grepl('Subconfluent', miRNAphenoDB[['GSE39061']]$`culture stage:ch1`, ignore.case = FALSE)]
    TreatmentsmiRNA <- mirnaFile[,grepl('Confluent', miRNAphenoDB[['GSE39061']]$`culture stage:ch1`, ignore.case = FALSE)]
  }
  
  if(exampleDataList[i] == 'GSE39061_3'){
    
    ControlsmRNA <- mrnaFile[,grepl('Subconfluent', mRNAphenoDB[['GSE39061']]$`culture stage:ch1`, ignore.case = FALSE)]
    TreatmentsmRNA <- mrnaFile[,grepl('Day 28', mRNAphenoDB[['GSE39061']]$`culture stage:ch1`, ignore.case = FALSE)]
    ControlsmiRNA <- mirnaFile[,grepl('Subconfluent', miRNAphenoDB[['GSE39061']]$`culture stage:ch1`, ignore.case = FALSE)]
    TreatmentsmiRNA <- mirnaFile[,grepl('Day 28', miRNAphenoDB[['GSE39061']]$`culture stage:ch1`, ignore.case = FALSE)]
  }
  
  if(exampleDataList[i] == 'GSE40321'){
    
    ControlsmRNA <- mrnaFile[,grepl('46,XY', mRNAphenoDB[['GSE40321']]$`genotype:ch1`, ignore.case = FALSE, fixed = TRUE)]
    TreatmentsmRNA <- mrnaFile[,grepl('47,XY,+8', mRNAphenoDB[['GSE40321']]$`genotype:ch1`, ignore.case = FALSE, fixed = TRUE)]
    ControlsmiRNA <- mirnaFile[,grepl('46,XY', miRNAphenoDB[['GSE40321']]$`genotype:ch1`, ignore.case = FALSE, fixed = TRUE)]
    TreatmentsmiRNA <- mirnaFile[,grepl('47,XY,+8', miRNAphenoDB[['GSE40321']]$`genotype:ch1`, ignore.case = FALSE, fixed = TRUE)]
  }
  
  if(exampleDataList[i] == 'GSE49697_1'){
    
    mirnaphenoTmp <- miRNAphenoDB[['GSE49697']][!grepl('poolBC', miRNAphenoDB[['GSE49697']]$title),]
    mirnaFile <- mirnaFile[,colnames(mirnaFile) %in% rownames(mirnaphenoTmp)]
    ControlsmRNA <- mrnaFile[,grepl('_US-48h', mRNAphenoDB[['GSE49697']]$title, ignore.case = FALSE, fixed = TRUE)]
    TreatmentsmRNA <- mrnaFile[,grepl('_S-48h', mRNAphenoDB[['GSE49697']]$title, ignore.case = FALSE, fixed = TRUE)]
    ControlsmiRNA <- mirnaFile[,grepl('_US-48h', mirnaphenoTmp$title, ignore.case = FALSE, fixed = TRUE)]
    TreatmentsmiRNA <- mirnaFile[,grepl('_S-48h', mirnaphenoTmp$title, ignore.case = FALSE, fixed = TRUE)]
  }
  
  if(exampleDataList[i] == 'GSE49697_2'){
    
    mirnaphenoTmp <- miRNAphenoDB[['GSE49697']][!grepl('poolBC', miRNAphenoDB[['GSE49697']]$title),]
    mirnaFile <- mirnaFile[,colnames(mirnaFile) %in% rownames(mirnaphenoTmp)]
    ControlsmRNA <- mrnaFile[,grepl('_US-24h', mRNAphenoDB[['GSE49697']]$title, ignore.case = FALSE, fixed = TRUE)]
    TreatmentsmRNA <- mrnaFile[,grepl('_S-24h', mRNAphenoDB[['GSE49697']]$title, ignore.case = FALSE, fixed = TRUE)]
    ControlsmiRNA <- mirnaFile[,grepl('_US-24h', mirnaphenoTmp$title, ignore.case = FALSE, fixed = TRUE)]
    TreatmentsmiRNA <- mirnaFile[,grepl('_S-24h', mirnaphenoTmp$title, ignore.case = FALSE, fixed = TRUE)]
  }
  
  if(exampleDataList[i] == 'GSE59702_1'){
    
    ControlsmRNA <- mrnaFile[,grepl('of fusion negative tumor', mRNAphenoDB[['GSE59702']]$title, ignore.case = FALSE, fixed = TRUE)]
    TreatmentsmRNA <- mrnaFile[,grepl('Fusion negative tumor', mRNAphenoDB[['GSE59702']]$title, ignore.case = FALSE, fixed = TRUE)]
    ControlsmiRNA <- mirnaFile[,grepl('of fusion negative tumor', miRNAphenoDB[['GSE59702']]$title, ignore.case = FALSE, fixed = TRUE)]
    TreatmentsmiRNA <- mirnaFile[,grepl('Fusion negative tumor', miRNAphenoDB[['GSE59702']]$title, ignore.case = FALSE, fixed = TRUE)]
  }
  
  if(exampleDataList[i] == 'GSE59702_2'){
    
    ControlsmRNA <- mrnaFile[,grepl('of fusion positive tumor', mRNAphenoDB[['GSE59702']]$title, ignore.case = FALSE, fixed = TRUE)]
    TreatmentsmRNA <- mrnaFile[,grepl('Fusion positive tumor', mRNAphenoDB[['GSE59702']]$title, ignore.case = FALSE, fixed = TRUE)]
    ControlsmiRNA <- mirnaFile[,grepl('of fusion positive tumor', miRNAphenoDB[['GSE59702']]$title, ignore.case = FALSE, fixed = TRUE)]
    TreatmentsmiRNA <- mirnaFile[,grepl('Fusion positive tumor', miRNAphenoDB[['GSE59702']]$title, ignore.case = FALSE, fixed = TRUE)]
  }
  
  if(exampleDataList[i] == 'GSE104268_1'){
    
    ControlsmRNA <- mrnaFile[,grepl('Control Triplicate', mRNAphenoDB[['GSE104268']]$title, ignore.case = FALSE, fixed = TRUE)]
    TreatmentsmRNA <- mrnaFile[,grepl('GSE Triplicate', mRNAphenoDB[['GSE104268']]$title, ignore.case = FALSE, fixed = TRUE)]
    ControlsmiRNA <- mirnaFile[,grepl('Control Triplicate', miRNAphenoDB[['GSE104268']]$title, ignore.case = FALSE, fixed = TRUE)]
    TreatmentsmiRNA <- mirnaFile[,grepl('GSE Triplicate', miRNAphenoDB[['GSE104268']]$title, ignore.case = FALSE, fixed = TRUE)]
  }
  
  if(exampleDataList[i] == 'GSE104268_2'){
    
    ControlsmRNA <- mrnaFile[,grepl('Control Triplicate', mRNAphenoDB[['GSE104268']]$title, ignore.case = FALSE, fixed = TRUE)]
    TreatmentsmRNA <- mrnaFile[,grepl('TSA Replicate', mRNAphenoDB[['GSE104268']]$title, ignore.case = FALSE, fixed = TRUE)]
    ControlsmiRNA <- mirnaFile[,grepl('Control Triplicate', miRNAphenoDB[['GSE104268']]$title, ignore.case = FALSE, fixed = TRUE)]
    TreatmentsmiRNA <- mirnaFile[,grepl('TSA Replicate', miRNAphenoDB[['GSE104268']]$title, ignore.case = FALSE, fixed = TRUE)]
  }
  
  if(exampleDataList[i] == 'GSE81867'){
    
    ControlsmRNA <- mrnaFile[,c(4,5,6)]
    TreatmentsmRNA <- mrnaFile[,c(1,2,3)]
    ControlsmiRNA <- mirnaFile[,c(4,5,6)]
    TreatmentsmiRNA <- mirnaFile[,c(1,2,3)]
  }
  
  if(exampleDataList[i] == 'GSE90604'){
    
    ControlsmRNA <- mrnaFile[,!grepl('Glioblastoma', mRNAphenoDB[['GSE90604']]$title, ignore.case = FALSE, fixed = TRUE)]
    TreatmentsmRNA <- mrnaFile[,grepl('Glioblastoma', mRNAphenoDB[['GSE90604']]$title, ignore.case = FALSE, fixed = TRUE)]
    ControlsmiRNA <- mirnaFile[,!grepl('Glioblastoma', miRNAphenoDB[['GSE90604']]$title, ignore.case = FALSE, fixed = TRUE)]
    TreatmentsmiRNA <- mirnaFile[,grepl('Glioblastoma', miRNAphenoDB[['GSE90604']]$title, ignore.case = FALSE, fixed = TRUE)]
  }
  
  if(exampleDataList[i] == 'GSE35389_1'){
    
    ControlsmRNA <- mrnaFile[,grepl('normal melanocyte', mRNAphenoDB[['GSE35389']]$title, ignore.case = FALSE, fixed = TRUE)]
    TreatmentsmRNA <- mrnaFile[,grepl('melanoma cell', mRNAphenoDB[['GSE35389']]$title, ignore.case = FALSE, fixed = TRUE)]
    ControlsmiRNA <- mirnaFile[,grepl('normal melanocyte', miRNAphenoDB[['GSE35389']]$title, ignore.case = FALSE, fixed = TRUE)]
    TreatmentsmiRNA <- mirnaFile[,grepl('melanoma cell', miRNAphenoDB[['GSE35389']]$title, ignore.case = FALSE, fixed = TRUE)]
  }
  
  if(exampleDataList[i] == 'GSE35389_2'){
    
    ControlsmRNA <- mrnaFile[,grepl('normal melanocyte', mRNAphenoDB[['GSE35389']]$title, ignore.case = FALSE, fixed = TRUE)]
    TreatmentsmRNA <- mrnaFile[,grepl('exosome', mRNAphenoDB[['GSE35389']]$title, ignore.case = FALSE, fixed = TRUE)]
    ControlsmiRNA <- mirnaFile[,grepl('normal melanocyte', miRNAphenoDB[['GSE35389']]$title, ignore.case = FALSE, fixed = TRUE)]
    TreatmentsmiRNA <- mirnaFile[,grepl('exosome', miRNAphenoDB[['GSE35389']]$title, ignore.case = FALSE, fixed = TRUE)]
  }
  
  if(exampleDataList[i] == 'GSE88721'){
    
    ControlsmRNA <- as.data.frame(mrnaFile[,grepl('Meningial Cells', mRNAphenoDB[['GSE88721']]$title, ignore.case = FALSE, fixed = TRUE)])
    TreatmentsmRNA <- mrnaFile[,grepl('Meningioma', mRNAphenoDB[['GSE88721']]$title, ignore.case = FALSE, fixed = TRUE)]
    colnames(ControlsmRNA) <- 'GSM2344707'
    rownames(ControlsmRNA) <- rownames(TreatmentsmRNA)
    ControlsmiRNA <- as.data.frame(mirnaFile[,grepl('Meningial Cells', miRNAphenoDB[['GSE88721']]$GSE88721.title, ignore.case = FALSE, fixed = TRUE)])
    TreatmentsmiRNA <- mirnaFile[,grepl('Meningioma', miRNAphenoDB[['GSE88721']]$GSE88721.title, ignore.case = FALSE, fixed = TRUE)]
    colnames(ControlsmiRNA) <- 'GSM2344692'
    rownames(ControlsmiRNA) <- rownames(TreatmentsmiRNA)
  }
  
  mirnaFile <- cbind(ControlsmiRNA, TreatmentsmiRNA)
  design.miRNA <- model.matrix(~ 0+factor(c(rep(1, ncol(ControlsmiRNA)), rep(2, ncol(TreatmentsmiRNA)))))
  colnames(design.miRNA) <- c("Control", "Treatment")
  fit <- lmFit(mirnaFile, design.miRNA)
  cont.matrix <- makeContrasts(TreatvsCtrl=Treatment-Control, levels=design.miRNA)
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2)
  mirna.limma <- topTable(fit2, adjust="BH", number = nrow(mirnaFile))
  
  load('targetscan-mir-gene-hsa.mirna.RData')
  mirna.limma$TargetScan <- 'No'
  commonmiRNA <- intersect(names(PredictedTargetsmiRNA), rownames(mirna.limma))
  mirna.limma[rownames(mirna.limma) %in% commonmiRNA,]$TargetScan <- 'Yes'
  load('miRNet-mir-gene-hsa.mirna.RData')
  mirna.limma$miRNet <- 'No'
  commonmiRNA <- intersect(names(PredictedTargetsmiRNA), rownames(mirna.limma))
  mirna.limma[rownames(mirna.limma) %in% commonmiRNA,]$miRNet <- 'Yes'
  
  mrnaFile <- cbind(ControlsmRNA, TreatmentsmRNA)
  design.mRNA <- model.matrix(~ 0+factor(c(rep(1, ncol(ControlsmRNA)), rep(2, ncol(TreatmentsmRNA)))))
  colnames(design.mRNA) <- c("Control", "Treatment")
  fit <- lmFit(mrnaFile, design.mRNA)
  cont.matrix <- makeContrasts(TreatvsCtrl=Treatment-Control, levels=design.mRNA)
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2)
  mrna.limma <- topTable(fit2, adjust="BH", number = nrow(mrnaFile))
  
  load('targetscan-mir-gene-hsa.gene.RData')
  mrna.limma$TargetScan <- 'No'
  commonmRNA <- intersect(names(PredictedTargetsGene), rownames(mrna.limma))
  mrna.limma[rownames(mrna.limma) %in% commonmRNA,]$TargetScan <- 'Yes'
  load('miRNet-mir-gene-hsa.gene.RData')
  mrna.limma$miRNet <- 'No'
  commonmRNA <- intersect(names(PredictedTargetsGene), rownames(mrna.limma))
  mrna.limma[rownames(mrna.limma) %in% commonmRNA,]$miRNet <- 'Yes'
  
  limmaListmiRNA[[exampleDataList[i]]] <- mirna.limma
  limmaListmRNA[[exampleDataList[i]]] <- mrna.limma
}