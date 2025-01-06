load(file='GenePlatformAnnotations.RData')
load(file='miRNA.id_hs.RData')

load('miRNAexprsDB.RData')
load('mRNAexprsDB.RData')
load('miRNAphenoDB.RData')
load('mRNAphenoDB.RData')
load('examplePhenoData.RData')
load('limmaListmiRNA.RData')
load('limmaListmRNA.RData')
load('mirnet_interactiondb.RData')

################################################################

parseGSEmRNA<- function(GSEmRNA, biomartVar){
    ## This function is created to convert probe ids to gene names in the mRNA expression file.
    ## Also expression difference of the mRNA data is calculated here
    
    GSEmRNA <- cbind(GSEmRNA, rownames(GSEmRNA))
    colnames(GSEmRNA)[ncol(GSEmRNA)] <- "probe.ids"
    GSEmRNA <- merge(biomartVar, GSEmRNA, by="probe.ids")
    ## Annotate the input with hgnc symbols taken from biomart
    GSEmRNA <- GSEmRNA[,-1]
    ## GSEmRNA[,2] <- as.numeric(as.character(GSEmRNA[,2]))
    ## GSEmRNA[,3] <- as.numeric(as.character(GSEmRNA[,3]))
    GSEmRNA <- aggregate(.~hgnc_symbol, data=GSEmRNA, mean)
    ## Some genes are represented with multiple probes.
    ## In order to merge them we use aggregate function and take mean of the probes that represent a gene
    rownames(GSEmRNA) <- GSEmRNA$hgnc_symbol
    GSEmRNA <- GSEmRNA[-1,-1]
    ## Parse the dataset for some more
    ## hgnc_symbols <- GSEmRNA[,1]
    ## GSEmRNA <- GSEmRNA$Treatment - GSEmRNA$Control
    ## Calculate the expression difference
    ## GSEmRNA <- cbind(as.character(hgnc_symbols), as.numeric(as.character(GSEmRNA)))
    ## colnames(GSEmRNA) <- c("hgnc_symbols","DiffmRNA")
    return(GSEmRNA)
    
}

#########################################################

ProcessExpressionData <- function(miRNADatabase, mRNAExpressionFile, miRNAExpressionFile, ExampleData, ExampleDataSelection, GenePlatformAnnotations, miRNA.identifiers, miRNAexprsDB, mRNAexprsDB, limmaListmiRNA, limmaListmRNA){
    ## This function is written to execute when side button in shiny application is clicked on.
    ## Requires few parameters i.e. QueryName holds the user query name entered, miRNADatabase holds the type of MTI database
    ## i.e. mirnet-targetscan, unionorintersect holds either union or intersect format of mirnet-targetscan combination is used, NULL if no combination database is being processed,
    ## mRNAExpressionFile and miRNAExpressionFile holds the expression data selected within Odyssey or uploaded by user,
    ## sep and sep2 holds the delimination format of the expression files,
    ## ExampleData stores a boolean value to check if Expression Data is selected from example data within Odyssey,
    ## ExampleDataSelection holds the name of the example data selected by the user, returns NULL if user uploads her own data instead.
    
    datacheck.mrna <- 'testedData'
    datacheck.mirna <- 'testedData'
    ControlsmRNA <- NULL
    TreatmentsmRNA <- NULL
    ControlsmiRNA <- NULL
    TreatmentsmiRNA <- NULL
    
    if(ExampleData == 'Upload Data'){
        ## Read mRNA expression data
        try(mrnaFile <- read.table(mRNAExpressionFile$datapath, sep='\t', header=TRUE,row.names=1, comment.char='!'), silent=TRUE)
        ## Read miRNA expression data
        try(mirnaFile <- read.table(miRNAExpressionFile$datapath, sep='\t', header=TRUE,row.names=1, comment.char='!'), silent=TRUE)
        try(range.mrnaFile <- range(mrnaFile), silent=TRUE)
        try(range.mirnaFile <- range(mirnaFile), silent=TRUE)
        
        if(!exists('mrnaFile')){
            if(!is.null(mRNAExpressionFile)){
                datacheck.mrna <- 'Data is not formatted correctly!'
            }
            else{
                datacheck.mrna <- 'Proceed the analysis without expression data!'
            }
        }
        if(!exists('mirnaFile')){
            if(!is.null(miRNAExpressionFile)){
                datacheck.mirna <- 'Data is not formatted correctly!'
            }
            else{
                datacheck.mirna <- 'Proceed the analysis without expression data!'
            }
        }
        
        if(exists('range.mrnaFile')){
            if(range.mrnaFile[2] > 20){
                mrnaFile <- log2(mrnaFile)
                try(range.mrnaFile <- range(mrnaFile), silent=TRUE)
            }
            if(range.mrnaFile[2] > 20){
                datacheck.mrna  <- 'Data is not formatted correctly!'
            }
        }
        else if(!exists('range.mrnaFile') & is.null(mRNAExpressionFile)){
            datacheck.mrna <- 'Proceed the analysis without expression data!'
        }
        else if(!exists('range.mrnaFile')){
            datacheck.mrna  <- 'Data is not formatted correctly!'
        }
        
        if(exists('range.mirnaFile')){
            if(range.mirnaFile[2] > 20){
                mirnaFile <- log2(mirnaFile)
                try(range.mirnaFile <- range(mirnaFile), silent=TRUE)
            }
            if(range.mirnaFile[2] > 20){
                datacheck.mirna  <- 'Data is not formatted correctly!'
            }
        }
        else if(!exists('range.mirnaFile') & is.null(miRNAExpressionFile)){
            datacheck.mirna <- 'Proceed the analysis without expression data!'
        }
        else if(!exists('range.mirnaFile')){
            datacheck.mirna  <- 'Data is not formatted correctly!'
        }
        
    }
    
    else if(ExampleData == 'Example Data') {
        
        ExampleDataSelectionStripped <- unlist(strsplit(ExampleDataSelection, '_'))[1]
        if(ExampleDataSelectionStripped == 'GSE35389'){
            mrnaFile <- mRNAexprsDB$GSE35389
            ## If example data is chosen, read the example mRNA data
            mirnaFile <- miRNAexprsDB$GSE35389
            ## If example data is chosen, read the example miRNA data
        }
        
        else if(ExampleDataSelectionStripped == 'GSE88721'){
            mrnaFile <- mRNAexprsDB$GSE88721
            mirnaFile <- miRNAexprsDB$GSE88721
        }
        
        else if(ExampleDataSelectionStripped == 'GSE39061'){
            mrnaFile <- mRNAexprsDB$GSE39061
            mirnaFile <- miRNAexprsDB$GSE39061
        }
        
        else if(ExampleDataSelectionStripped == 'GSE49697'){
            mrnaFile <- mRNAexprsDB$GSE49697
            mirnaFile <- miRNAexprsDB$GSE49697
        }
        
        else if(ExampleDataSelectionStripped == 'GSE25402'){
            mrnaFile <- mRNAexprsDB$GSE25402
            mirnaFile <- miRNAexprsDB$GSE25402
        }
        
        else if(ExampleDataSelectionStripped == 'GSE32539'){
            mrnaFile <- mRNAexprsDB$GSE32539
            mirnaFile <- miRNAexprsDB$GSE32539
        }
        
        else if(ExampleDataSelectionStripped == 'GSE34681'){
            mrnaFile <- mRNAexprsDB$GSE34681
            mirnaFile <- miRNAexprsDB$GSE34681
        }
        
        else if(ExampleDataSelectionStripped == 'GSE38617'){
            mrnaFile <- mRNAexprsDB$GSE38617
            mirnaFile <- miRNAexprsDB$GSE38617
        }
        
        else if(ExampleDataSelectionStripped == 'GSE40321'){
            mrnaFile <- mRNAexprsDB$GSE40321
            mirnaFile <- miRNAexprsDB$GSE40321
        }
        
        else if(ExampleDataSelectionStripped == 'GSE59702'){
            mrnaFile <- mRNAexprsDB$GSE59702
            mirnaFile <- miRNAexprsDB$GSE59702
        }
        
        else if(ExampleDataSelectionStripped == 'GSE81867'){
            mrnaFile <- mRNAexprsDB$GSE81867
            mirnaFile <- miRNAexprsDB$GSE81867
        }
        
        else if(ExampleDataSelectionStripped == 'GSE104268'){
            mrnaFile <- mRNAexprsDB$GSE104268
            mirnaFile <- miRNAexprsDB$GSE104268
        }
        
        else if(ExampleDataSelectionStripped == 'GSE90604'){
            mrnaFile <- mRNAexprsDB$GSE90604
            mirnaFile <- miRNAexprsDB$GSE90604
        }
        
    }
    else if(ExampleData == 'No Data'){
        datacheck.mrna <- 'Proceed the analysis without expression data!'
        datacheck.mirna <- 'Proceed the analysis without expression data!'
    }
    
    if(datacheck.mrna == 'Data is not formatted correctly!' | datacheck.mrna == 'Proceed the analysis without expression data!'){
        mrnaFile <- 'No Data'
    }
    else if(datacheck.mrna == 'testedData'){
        if(ExampleData == 'Upload Data'){
            mrnaFile <- parseGSEmRNA(mrnaFile, GenePlatformAnnotations)
            mrnaFile <- round_df(mrnaFile)
        }
    }
    
    if(datacheck.mirna == 'Data is not formatted correctly!' | datacheck.mirna == 'Proceed the analysis without expression data!'){
        mirnaFile <- 'No Data'
    }
    else if(datacheck.mirna == 'testedData'){
        mirnaFile <- mirnaFile[!grepl('-mir', rownames(mirnaFile)),]
        rownames(mirnaFile) <- make.unique(toupper(rownames(mirnaFile)))
        mirnaFile <- round_df(mirnaFile)
        if(ExampleData != 'Example Data'){
            dimensionsExpmiRNA <- nrow(mirnaFile)
            mirnaFile.annotated <- merge(mirnaFile, miRNA.identifiers, by=0)
            if(nrow(mirnaFile.annotated) / dimensionsExpmiRNA * 100 >= 5){
                mirnaFile <- mirnaFile.annotated
                mirnaFile <- mirnaFile[,-1]
                nSamples <- ncol(mirnaFile)
                mirnaFile <- mirnaFile[,-(nSamples-1)]
                
                if(sum(duplicated(mirnaFile$v22)) > 0){
                    mirnaFile <- aggregate(.~v22, data = mirnaFile, mean)
                    rownames(mirnaFile) <- mirnaFile$v22
                    mirnaFile <- mirnaFile[,-1]
                }
                else {
                    rownames(mirnaFile) <- mirnaFile$v22
                    mirnaFile <- mirnaFile[,-(nSamples-1)]
                }
            }
            else {
                mirnaFile <- cbind(rownames(mirnaFile), mirnaFile)
                colnames(mirnaFile)[1] <- 'input_id'
                ## mirnaFile$input_id <- gsub(pattern = '_ST', replacement = '', mirnaFile$input_id)
                ## mirnaFile$input_id <- gsub(pattern = '-STAR', replacement = '*', mirnaFile$input_id)
                mirnaFile.annotated <- merge(mirnaFile, miRNA.identifiers, by='input_id')
                if(nrow(mirnaFile.annotated) / dimensionsExpmiRNA * 100 >= 5){
                    mirnaFile <- mirnaFile.annotated
                    mirnaFile <- mirnaFile[,-1]
                    if(sum(duplicated(mirnaFile$v22)) > 0){
                        mirnaFile <- aggregate(.~v22, data = mirnaFile, mean)
                        rownames(mirnaFile) <- mirnaFile$v22
                        mirnaFile <- mirnaFile[,-1]
                    }
                    else {
                        nSamples <- ncol(mirnaFile)
                        rownames(mirnaFile) <- mirnaFile$v22
                        mirnaFile <- mirnaFile[,-nSamples]
                    }
                }
                else{
                    datacheck.mrna == 'Data is not formatted correctly!'
                    mirnaFile <- 'No Data'
                }
            }
        }
    }
    
    if(ExampleData == 'Example Data'){
        geneIDs <- rownames(mrnaFile)
        ## Set the example data columns as treatments and controls
        
        if(ExampleDataSelection == 'GSE25402'){
            
            ControlsmRNA <- mrnaFile[,grepl('Nonobese', mRNAphenoDB[['GSE25402']]$title, ignore.case = FALSE)]
            TreatmentsmRNA <- mrnaFile[,grepl('Obese', mRNAphenoDB[['GSE25402']]$title, ignore.case = FALSE)]
            ControlsmiRNA <- mirnaFile[,grepl('Nonobese', miRNAphenoDB[['GSE25402']]$title, ignore.case = FALSE)]
            TreatmentsmiRNA <- mirnaFile[,grepl('Obese', miRNAphenoDB[['GSE25402']]$title, ignore.case = FALSE)]
        }
        
        if(ExampleDataSelection == 'GSE32539_1'){
            
            ControlsmRNA <- mrnaFile[,grepl('control', mRNAphenoDB[['GSE32539']]$`final diagnosis:ch1`, ignore.case = FALSE)]
            TreatmentsmRNA <- mrnaFile[,grepl('IPF/UIP', mRNAphenoDB[['GSE32539']]$`final diagnosis:ch1`, ignore.case = FALSE)]
            ControlsmiRNA <- mirnaFile[,grepl('Control', miRNAphenoDB[['GSE32539']]$`final diagnosis:ch1`, ignore.case = FALSE)]
            TreatmentsmiRNA <- mirnaFile[,grepl('IPF/UIP', miRNAphenoDB[['GSE32539']]$`final diagnosis:ch1`, ignore.case = FALSE)]
        }
        
        if(ExampleDataSelection == 'GSE32539_2'){
            
            ControlsmRNA <- mrnaFile[,grepl('control', mRNAphenoDB[['GSE32539']]$`final diagnosis:ch1`, ignore.case = FALSE)]
            TreatmentsmRNA <- mrnaFile[,grepl('NSIP', mRNAphenoDB[['GSE32539']]$`final diagnosis:ch1`, ignore.case = FALSE)]
            ControlsmiRNA <- mirnaFile[,grepl('Control', miRNAphenoDB[['GSE32539']]$`final diagnosis:ch1`, ignore.case = FALSE)]
            TreatmentsmiRNA <- mirnaFile[,grepl('NSIP', miRNAphenoDB[['GSE32539']]$`final diagnosis:ch1`, ignore.case = FALSE)]
        }
        
        if(ExampleDataSelection == 'GSE32539_3'){
            
            ControlsmRNA <- mrnaFile[,grepl('control', mRNAphenoDB[['GSE32539']]$`final diagnosis:ch1`, ignore.case = FALSE)]
            TreatmentsmRNA <- mrnaFile[,grepl('RB-ILD', mRNAphenoDB[['GSE32539']]$`final diagnosis:ch1`, ignore.case = FALSE)]
            ControlsmiRNA <- mirnaFile[,grepl('Control', miRNAphenoDB[['GSE32539']]$`final diagnosis:ch1`, ignore.case = FALSE)]
            TreatmentsmiRNA <- mirnaFile[,grepl('RB-ILD', miRNAphenoDB[['GSE32539']]$`final diagnosis:ch1`, ignore.case = FALSE)]
        }
        
        if(ExampleDataSelection == 'GSE34681_1'){
            
            ControlsmRNA <- mrnaFile[,grepl('control siRNA', mRNAphenoDB[['GSE34681']]$`treatment:ch1`, ignore.case = FALSE)]
            TreatmentsmRNA <- mrnaFile[,grepl('Ars2 siRNA', mRNAphenoDB[['GSE34681']]$`treatment:ch1`, ignore.case = FALSE)]
            ControlsmiRNA <- mirnaFile[,grepl('control siRNA', miRNAphenoDB[['GSE34681']]$`treatment:ch1`, ignore.case = FALSE)]
            TreatmentsmiRNA <- mirnaFile[,grepl('Ars2 siRNA', miRNAphenoDB[['GSE34681']]$`treatment:ch1`, ignore.case = FALSE)]
        }
        
        if(ExampleDataSelection == 'GSE34681_2'){
            
            ControlsmRNA <- mrnaFile[,grepl('control siRNA', mRNAphenoDB[['GSE34681']]$`treatment:ch1`, ignore.case = FALSE)]
            TreatmentsmRNA <- mrnaFile[,grepl('DGCR8 siRNA', mRNAphenoDB[['GSE34681']]$`treatment:ch1`, ignore.case = FALSE)]
            ControlsmiRNA <- mirnaFile[,grepl('control siRNA', miRNAphenoDB[['GSE34681']]$`treatment:ch1`, ignore.case = FALSE)]
            TreatmentsmiRNA <- mirnaFile[,grepl('DGCR8 siRNA', miRNAphenoDB[['GSE34681']]$`treatment:ch1`, ignore.case = FALSE)]
        }
        
        if(ExampleDataSelection == 'GSE38617'){
            
            ControlsmRNA <- mrnaFile[,grepl('healthy', mRNAphenoDB[['GSE38617']]$`disease state:ch1`, ignore.case = FALSE)]
            TreatmentsmRNA <- mrnaFile[,grepl('oral lichen planus', mRNAphenoDB[['GSE38617']]$`disease state:ch1`, ignore.case = FALSE)]
            ControlsmiRNA <- mirnaFile[,grepl('healthy', miRNAphenoDB[['GSE38617']]$`disease state:ch1`, ignore.case = FALSE)]
            TreatmentsmiRNA <- mirnaFile[,grepl('oral lichen planus', miRNAphenoDB[['GSE38617']]$`disease state:ch1`, ignore.case = FALSE)]
        }
        
        if(ExampleDataSelection == 'GSE39061_1'){
            
            ControlsmRNA <- mrnaFile[,grepl('Confluent', mRNAphenoDB[['GSE39061']]$`culture stage:ch1`, ignore.case = FALSE)]
            TreatmentsmRNA <- mrnaFile[,grepl('Day 28', mRNAphenoDB[['GSE39061']]$`culture stage:ch1`, ignore.case = FALSE)]
            ControlsmiRNA <- mirnaFile[,grepl('Confluent', miRNAphenoDB[['GSE39061']]$`culture stage:ch1`, ignore.case = FALSE)]
            TreatmentsmiRNA <- mirnaFile[,grepl('Day 28', miRNAphenoDB[['GSE39061']]$`culture stage:ch1`, ignore.case = FALSE)]
        }
        
        if(ExampleDataSelection == 'GSE39061_2'){
            
            ControlsmRNA <- mrnaFile[,grepl('Subconfluent', mRNAphenoDB[['GSE39061']]$`culture stage:ch1`, ignore.case = FALSE)]
            TreatmentsmRNA <- mrnaFile[,grepl('Confluent', mRNAphenoDB[['GSE39061']]$`culture stage:ch1`, ignore.case = FALSE)]
            ControlsmiRNA <- mirnaFile[,grepl('Subconfluent', miRNAphenoDB[['GSE39061']]$`culture stage:ch1`, ignore.case = FALSE)]
            TreatmentsmiRNA <- mirnaFile[,grepl('Confluent', miRNAphenoDB[['GSE39061']]$`culture stage:ch1`, ignore.case = FALSE)]
        }
        
        if(ExampleDataSelection == 'GSE39061_3'){
            
            ControlsmRNA <- mrnaFile[,grepl('Subconfluent', mRNAphenoDB[['GSE39061']]$`culture stage:ch1`, ignore.case = FALSE)]
            TreatmentsmRNA <- mrnaFile[,grepl('Day 28', mRNAphenoDB[['GSE39061']]$`culture stage:ch1`, ignore.case = FALSE)]
            ControlsmiRNA <- mirnaFile[,grepl('Subconfluent', miRNAphenoDB[['GSE39061']]$`culture stage:ch1`, ignore.case = FALSE)]
            TreatmentsmiRNA <- mirnaFile[,grepl('Day 28', miRNAphenoDB[['GSE39061']]$`culture stage:ch1`, ignore.case = FALSE)]
        }
        
        if(ExampleDataSelection == 'GSE40321'){
            
            ControlsmRNA <- mrnaFile[,grepl('46,XY', mRNAphenoDB[['GSE40321']]$`genotype:ch1`, ignore.case = FALSE, fixed = TRUE)]
            TreatmentsmRNA <- mrnaFile[,grepl('47,XY,+8', mRNAphenoDB[['GSE40321']]$`genotype:ch1`, ignore.case = FALSE, fixed = TRUE)]
            ControlsmiRNA <- mirnaFile[,grepl('46,XY', miRNAphenoDB[['GSE40321']]$`genotype:ch1`, ignore.case = FALSE, fixed = TRUE)]
            TreatmentsmiRNA <- mirnaFile[,grepl('47,XY,+8', miRNAphenoDB[['GSE40321']]$`genotype:ch1`, ignore.case = FALSE, fixed = TRUE)]
        }
        
        if(ExampleDataSelection == 'GSE49697_1'){
            
            mirnaphenoTmp <- miRNAphenoDB[['GSE49697']][!grepl('poolBC', miRNAphenoDB[['GSE49697']]$title),]
            mirnaFile <- mirnaFile[,colnames(mirnaFile) %in% rownames(mirnaphenoTmp)]
            ControlsmRNA <- mrnaFile[,grepl('_US-48h', mRNAphenoDB[['GSE49697']]$title, ignore.case = FALSE, fixed = TRUE)]
            TreatmentsmRNA <- mrnaFile[,grepl('_S-48h', mRNAphenoDB[['GSE49697']]$title, ignore.case = FALSE, fixed = TRUE)]
            ControlsmiRNA <- mirnaFile[,grepl('_US-48h', mirnaphenoTmp$title, ignore.case = FALSE, fixed = TRUE)]
            TreatmentsmiRNA <- mirnaFile[,grepl('_S-48h', mirnaphenoTmp$title, ignore.case = FALSE, fixed = TRUE)]
        }
        
        if(ExampleDataSelection == 'GSE49697_2'){
            
            mirnaphenoTmp <- miRNAphenoDB[['GSE49697']][!grepl('poolBC', miRNAphenoDB[['GSE49697']]$title),]
            mirnaFile <- mirnaFile[,colnames(mirnaFile) %in% rownames(mirnaphenoTmp)]
            ControlsmRNA <- mrnaFile[,grepl('_US-24h', mRNAphenoDB[['GSE49697']]$title, ignore.case = FALSE, fixed = TRUE)]
            TreatmentsmRNA <- mrnaFile[,grepl('_S-24h', mRNAphenoDB[['GSE49697']]$title, ignore.case = FALSE, fixed = TRUE)]
            ControlsmiRNA <- mirnaFile[,grepl('_US-24h', mirnaphenoTmp$title, ignore.case = FALSE, fixed = TRUE)]
            TreatmentsmiRNA <- mirnaFile[,grepl('_S-24h', mirnaphenoTmp$title, ignore.case = FALSE, fixed = TRUE)]
        }
        
        if(ExampleDataSelection == 'GSE59702_1'){
            
            ControlsmRNA <- mrnaFile[,grepl('of fusion negative tumor', mRNAphenoDB[['GSE59702']]$title, ignore.case = FALSE, fixed = TRUE)]
            TreatmentsmRNA <- mrnaFile[,grepl('Fusion negative tumor', mRNAphenoDB[['GSE59702']]$title, ignore.case = FALSE, fixed = TRUE)]
            ControlsmiRNA <- mirnaFile[,grepl('of fusion negative tumor', miRNAphenoDB[['GSE59702']]$title, ignore.case = FALSE, fixed = TRUE)]
            TreatmentsmiRNA <- mirnaFile[,grepl('Fusion negative tumor', miRNAphenoDB[['GSE59702']]$title, ignore.case = FALSE, fixed = TRUE)]
        }
        
        if(ExampleDataSelection == 'GSE59702_2'){
            
            ControlsmRNA <- mrnaFile[,grepl('of fusion positive tumor', mRNAphenoDB[['GSE59702']]$title, ignore.case = FALSE, fixed = TRUE)]
            TreatmentsmRNA <- mrnaFile[,grepl('Fusion positive tumor', mRNAphenoDB[['GSE59702']]$title, ignore.case = FALSE, fixed = TRUE)]
            ControlsmiRNA <- mirnaFile[,grepl('of fusion positive tumor', miRNAphenoDB[['GSE59702']]$title, ignore.case = FALSE, fixed = TRUE)]
            TreatmentsmiRNA <- mirnaFile[,grepl('Fusion positive tumor', miRNAphenoDB[['GSE59702']]$title, ignore.case = FALSE, fixed = TRUE)]
        }
        
        if(ExampleDataSelection == 'GSE104268_1'){
            
            ControlsmRNA <- mrnaFile[,grepl('Control Triplicate', mRNAphenoDB[['GSE104268']]$title, ignore.case = FALSE, fixed = TRUE)]
            TreatmentsmRNA <- mrnaFile[,grepl('GSE Triplicate', mRNAphenoDB[['GSE104268']]$title, ignore.case = FALSE, fixed = TRUE)]
            ControlsmiRNA <- mirnaFile[,grepl('Control Triplicate', miRNAphenoDB[['GSE104268']]$title, ignore.case = FALSE, fixed = TRUE)]
            TreatmentsmiRNA <- mirnaFile[,grepl('GSE Triplicate', miRNAphenoDB[['GSE104268']]$title, ignore.case = FALSE, fixed = TRUE)]
        }
        
        if(ExampleDataSelection == 'GSE104268_2'){
            
            ControlsmRNA <- mrnaFile[,grepl('Control Triplicate', mRNAphenoDB[['GSE104268']]$title, ignore.case = FALSE, fixed = TRUE)]
            TreatmentsmRNA <- mrnaFile[,grepl('TSA Replicate', mRNAphenoDB[['GSE104268']]$title, ignore.case = FALSE, fixed = TRUE)]
            ControlsmiRNA <- mirnaFile[,grepl('Control Triplicate', miRNAphenoDB[['GSE104268']]$title, ignore.case = FALSE, fixed = TRUE)]
            TreatmentsmiRNA <- mirnaFile[,grepl('TSA Replicate', miRNAphenoDB[['GSE104268']]$title, ignore.case = FALSE, fixed = TRUE)]
        }
        
        if(ExampleDataSelection == 'GSE81867'){
            
            ControlsmRNA <- mrnaFile[,c(4,5,6)]
            TreatmentsmRNA <- mrnaFile[,c(1,2,3)]
            ControlsmiRNA <- mirnaFile[,c(4,5,6)]
            TreatmentsmiRNA <- mirnaFile[,c(1,2,3)]
        }
        
        if(ExampleDataSelection == 'GSE90604'){
            
            ControlsmRNA <- mrnaFile[,!grepl('Glioblastoma', mRNAphenoDB[['GSE90604']]$title, ignore.case = FALSE, fixed = TRUE)]
            TreatmentsmRNA <- mrnaFile[,grepl('Glioblastoma', mRNAphenoDB[['GSE90604']]$title, ignore.case = FALSE, fixed = TRUE)]
            ControlsmiRNA <- mirnaFile[,!grepl('Glioblastoma', miRNAphenoDB[['GSE90604']]$title, ignore.case = FALSE, fixed = TRUE)]
            TreatmentsmiRNA <- mirnaFile[,grepl('Glioblastoma', miRNAphenoDB[['GSE90604']]$title, ignore.case = FALSE, fixed = TRUE)]
        }
        
        if(ExampleDataSelection == 'GSE35389_1'){
            
            ControlsmRNA <- mrnaFile[,grepl('normal melanocyte', mRNAphenoDB[['GSE35389']]$title, ignore.case = FALSE, fixed = TRUE)]
            TreatmentsmRNA <- mrnaFile[,grepl('melanoma cell', mRNAphenoDB[['GSE35389']]$title, ignore.case = FALSE, fixed = TRUE)]
            ControlsmiRNA <- mirnaFile[,grepl('normal melanocyte', miRNAphenoDB[['GSE35389']]$title, ignore.case = FALSE, fixed = TRUE)]
            TreatmentsmiRNA <- mirnaFile[,grepl('melanoma cell', miRNAphenoDB[['GSE35389']]$title, ignore.case = FALSE, fixed = TRUE)]
        }
        
        if(ExampleDataSelection == 'GSE35389_2'){
            
            ControlsmRNA <- mrnaFile[,grepl('normal melanocyte', mRNAphenoDB[['GSE35389']]$title, ignore.case = FALSE, fixed = TRUE)]
            TreatmentsmRNA <- mrnaFile[,grepl('exosome', mRNAphenoDB[['GSE35389']]$title, ignore.case = FALSE, fixed = TRUE)]
            ControlsmiRNA <- mirnaFile[,grepl('normal melanocyte', miRNAphenoDB[['GSE35389']]$title, ignore.case = FALSE, fixed = TRUE)]
            TreatmentsmiRNA <- mirnaFile[,grepl('exosome', miRNAphenoDB[['GSE35389']]$title, ignore.case = FALSE, fixed = TRUE)]
        }
        
        if(ExampleDataSelection == 'GSE88721'){
            
            ControlsmRNA <- as.data.frame(mrnaFile[,grepl('Meningial Cells', mRNAphenoDB[['GSE88721']]$title, ignore.case = FALSE, fixed = TRUE)])
            TreatmentsmRNA <- mrnaFile[,grepl('Meningioma', mRNAphenoDB[['GSE88721']]$title, ignore.case = FALSE, fixed = TRUE)]
            colnames(ControlsmRNA) <- 'GSM2344707'
            rownames(ControlsmRNA) <- rownames(TreatmentsmRNA)
            ControlsmiRNA <- as.data.frame(mirnaFile[,grepl('Meningial Cells', miRNAphenoDB[['GSE88721']]$GSE88721.title, ignore.case = FALSE, fixed = TRUE)])
            TreatmentsmiRNA <- mirnaFile[,grepl('Meningioma', miRNAphenoDB[['GSE88721']]$GSE88721.title, ignore.case = FALSE, fixed = TRUE)]
            colnames(ControlsmiRNA) <- 'GSM2344692'
            rownames(ControlsmiRNA) <- rownames(TreatmentsmiRNA)
        }
    }
    
    return(list(mirnaFile, mrnaFile, datacheck.mirna, datacheck.mrna, miRNAphenoDB, mRNAphenoDB, TreatmentsmiRNA, ControlsmiRNA, TreatmentsmRNA, ControlsmRNA))
}

ProcessDGExMTIdb <- function(mirnaFile, mrnaFile, datacheck.mirna, datacheck.mrna, miRNADatabase, unionorintersect, ExampleData, TreatmentsmiRNA, ControlsmiRNA, TreatmentsmRNA, ControlsmRNA){
    
    ## Query is converted to all upper cases for convenient search in the database.
    ## Errors caused by mixing lower case characters with upper case characters are prevented.
    downloadHandle <- ''
    
    if(ExampleData == 'Upload Data'){
        if(mrnaFile != 'No Data'){
            geneIDs <- rownames(mrnaFile)
        }
    }
    
    if(mrnaFile != 'No Data'){
        if(is.null(ControlsmRNA) | is.null(TreatmentsmRNA)){
            return('NoSampleSelection')
        }
        
        if(ncol(as.matrix(ControlsmRNA)) == 0 | ncol(as.matrix(TreatmentsmRNA)) == 0){
            return('NoSampleSelection')
        }
        
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
        mrna.limma <- round_df(mrna.limma)
        mrna.limma <- as.data.frame.matrix(mrna.limma)
        mrnaFile2 <- cbind(rownames(mrna.limma), mrna.limma$logFC)
        
        colnames(mrnaFile2) <- c('Gene','DiffmRNA')
        mrnaFile2 <- as.data.frame.matrix(mrnaFile2)
        ## Set column names to be used later with merge function.
        ## merge function requires common column names for objects to be merged...
    }
    else{
        mrnaFile2 <- 0
        mrna.limma <- 'No Data'
    }
    
    if(mirnaFile != 'No Data'){
        ## If multiple columns are selected as treatment or control samples, process these for differential expression analysis...
        ## Mean of rows are taken here for that purpose. Median is another option that could be used..
        if(is.null(ControlsmiRNA) | is.null(TreatmentsmiRNA)){
            return('NoSampleSelection')
        }
        
        miRNAIDs <- rownames(mirnaFile)
        if(ncol(as.matrix(ControlsmiRNA)) == 0 | ncol(as.matrix(TreatmentsmiRNA)) == 0){
            return('NoSampleSelection')
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
        mirna.limma <<- mirna.limma
        commonmiRNA <- intersect(names(PredictedTargetsmiRNA), rownames(mirna.limma))
        mirna.limma[rownames(mirna.limma) %in% commonmiRNA,]$TargetScan <- 'Yes'
        load('miRNet-mir-gene-hsa.mirna.RData')
        mirna.limma$miRNet <- 'No'
        commonmiRNA <- intersect(names(PredictedTargetsmiRNA), rownames(mirna.limma))
        mirna.limma[rownames(mirna.limma) %in% commonmiRNA,]$miRNet <- 'Yes'
        mirna.limma <- round_df(mirna.limma)
        mirna.limma <- as.data.frame.matrix(mirna.limma)
        mirnaFile2 <- cbind(rownames(mirna.limma), mirna.limma$logFC)
        
        colnames(mirnaFile2) <- c("miRNA", "DiffmiRNA")
        mirnaFile2 <- as.data.frame.matrix(mirnaFile2)
        
    }
    else {
        mirnaFile2 <- 0
        mirna.limma <- 'No Data'
    }
    
    return(list(mirnaFile, mrnaFile, mirnaFile2, mrnaFile2, mirna.limma, mrna.limma))
    ## Return processed variables ...
}

#########################################################

ProcessSIF <- function(QueryName, querySelection, mirnaFile, mrnaFile, mirnaFile2, mrnaFile2, SecondaryInteractions, miRNADatabase, unionorintersect, PredictedTargetsmiRNA, PredictedTargetsGene) {
    ## Expression file is merged with Cytoscape file prepared in the previous function via this script...
    
    completeDataUnmatch <- ''
    QueryName <- toupper(QueryName)
    ## Convert query to all upper case to handle errors caused by lower case and upper case mixing by user query text box ...
    if(querySelection == "NotFound"){return()}
    
    if(miRNADatabase == "TargetScan"){
        ## Load a different MTIdb and GO annotation file based on user selection...
        
        load(file='targetscan-mir-gene-hsa.gene.RData')
        load(file='targetscan-mir-gene-hsa.mirna.RData')
    }
    ## Read the targetscan database
    
    else if(miRNADatabase == "miRNet"){
        
        load(file='miRNet-mir-gene-hsa.gene.RData')
        load(file='miRNet-mir-gene-hsa.mirna.RData')
    }
    
    else if(miRNADatabase == "TargetScan&miRNet"){
        
        if(unionorintersect == "Intersect"){
            
            load(file='mirnet.targetscan.intersect.gene.RData')
            load(file='mirnet.targetscan.intersect.mirna.RData')
            
            downloadHandle <- 'Intersect'
        }
        
        else if(unionorintersect == "Union"){
            
            load(file='mirnet.targetscan.union.gene.RData')
            load(file='mirnet.targetscan.union.mirna.RData')
            
            downloadHandle <- 'Union'
        }
    }
    
    queryMasterChecker <- NULL
    tmpIDchecker <- TRUE
    multipleMatches <- FALSE
    queryModal <- NULL
    
    if(querySelection == 'Gene'){
        if(mrnaFile == 'No Data'){
            if(is.null(PredictedTargetsGene[[QueryName]])){
                queryMasterChecker <- names(PredictedTargetsGene)[grep(QueryName, names(PredictedTargetsGene))]
                if(length(queryMasterChecker) > 1){
                    multipleMatches <- TRUE
                    queryModal <- 'Query is not present in selected Interaction Database! Select from similar identifiers or try a different query. Here are a few similar matches ...'
                    return(list(queryMasterChecker, multipleMatches, queryModal))
                }
                else if(length(queryMasterChecker) == 0){
                    queryModal <- 'Query is not present in selected Interaction Database! Please try a different query...'
                    multipleMatches <- 'NoMatch'
                    return(list(queryMasterChecker, multipleMatches, queryModal))
                }
                else if(length(queryMasterChecker) == 1){
                    ## Not for one expression match but instead for MTIdb one match.... 
                    queryModal <- paste('Query is not present in selected Interaction Database! Proceede with close match ', queryMasterChecker, '?', sep='') 
                    multipleMatches <- 'OneExpressionMatch'
                    QueryName <- queryMasterChecker
                }
            }
            else {
                queryMasterChecker <- QueryName
            }
        }
        else {
            if(!is.null(PredictedTargetsGene[[QueryName]])){
                queryMasterChecker <- intersect(names(PredictedTargetsGene)[grep(QueryName, names(PredictedTargetsGene))], rownames(mrnaFile))
                if(length(queryMasterChecker) == 0){
                    queryModal <- 'Query is not present in expression data identifiers! Proceeding without full expression data benefits...'
                    multipleMatches <- 'NoExpressionMatch'
                }
            }
            else{
                if(length(grep(QueryName, names(PredictedTargetsGene))) > 1) {
                    
                    queryMasterChecker <- intersect(names(PredictedTargetsGene)[grep(QueryName, names(PredictedTargetsGene))], rownames(mrnaFile))
                    if(length(queryMasterChecker) > 1){
                        queryModal <- 'Query is not present in selected Interaction Database! Select from similar identifiers or try a different query. Here are a few similar matches ...'
                        multipleMatches <- TRUE
                        return(list(queryMasterChecker, multipleMatches, queryModal))
                    }
                    else if(length(queryMasterChecker) == 0){
                        queryModal <- 'Query is not present in selected Interaction Database! Please try a different query...'
                        multipleMatches <- 'NoMatch'
                        return(list(queryMasterChecker, multipleMatches, queryModal))
                    }
                    else if(length(queryMasterChecker) == 1){
                        queryModal <- paste('Query is not present in selected Interaction Database! Proceede with close match: ', queryMasterChecker, '?', sep='') 
                        multipleMatches <- 'OneExpressionMatch'
                        QueryName <- queryMasterChecker
                    }
                }
                if(length(grep(QueryName, names(PredictedTargetsGene))) == 1) {
                    queryMasterChecker <- intersect(names(PredictedTargetsGene)[grep(QueryName, names(PredictedTargetsGene))], rownames(mrnaFile))
                    
                    if(length(queryMasterChecker) == 1){
                        queryModal <- paste('Query is not present in selected Interaction Database! Proceede with close match: ', queryMasterChecker, '?', sep='')
                        QueryName <- queryMasterChecker
                        multipleMatches <- 'OneExpressionMatch'
                    }
                    else if(length(queryMasterChecker) == 0){
                        queryModal <- 'Query is not present in selected Interaction Database! Please try a different query...'
                        multipleMatches <- 'NoMatch'
                        return(list(queryMasterChecker, multipleMatches, queryModal))
                    }
                }
                else if(length(grep(QueryName, names(PredictedTargetsGene))) == 0) {
                    tmpIDchecker <- FALSE
                }
            }
        }
    }
    else if(querySelection == 'miRNA'){
        if(mirnaFile == 'No Data'){
            if(is.null(PredictedTargetsmiRNA[[QueryName]])){
                queryMasterChecker <- names(PredictedTargetsmiRNA)[grep(QueryName, names(PredictedTargetsmiRNA))]
                if(length(queryMasterChecker) > 1){
                    queryModal <- 'Query is not present in selected Interaction Database! Select from similar identifiers or try a different query. Here are a few similar matches ...'
                    multipleMatches <- TRUE
                    return(list(queryMasterChecker, multipleMatches, queryModal))
                }
                else if(length(queryMasterChecker) == 0){
                    queryModal <- 'Query is not present in selected Interaction Database! Please try a different query...'
                    multipleMatches <- 'NoMatch'
                    return(list(queryMasterChecker, multipleMatches, queryModal))
                }
                else if(length(queryMasterChecker) == 1){
                    ## Not for one expression match but instead for MTIdb one match.... 
                    queryModal <- paste('Query is not present in selected Interaction Database! Proceede with close match ', queryMasterChecker, '?', sep='') 
                    multipleMatches <- 'OneExpressionMatch'
                    QueryName <- queryMasterChecker
                }
            }
            else{
                queryMasterChecker <- QueryName
            }
        }
        else {
            if(!is.null(PredictedTargetsmiRNA[[QueryName]])){
                queryMasterChecker <- intersect(names(PredictedTargetsmiRNA)[grep(QueryName, names(PredictedTargetsmiRNA))], rownames(mirnaFile))
                if(length(queryMasterChecker) == 0){
                    queryModal <- 'Query is not present in expression data identifiers! Proceeding without full expression data benefits...'
                    multipleMatches <- 'NoExpressionMatch'
                }
            }
            else{
                if(length(grep(QueryName, names(PredictedTargetsmiRNA))) > 1) {
                    
                    queryMasterChecker <- intersect(names(PredictedTargetsmiRNA)[grep(QueryName, names(PredictedTargetsmiRNA))], rownames(mirnaFile))
                    if(length(queryMasterChecker) > 1){
                        queryModal <- 'Query is not present in selected Interaction Database! Select from similar identifiers or try a different query. Here are a few similar matches ...'
                        multipleMatches <- TRUE
                        return(list(queryMasterChecker, multipleMatches, queryModal))
                    }
                    
                    else if(length(queryMasterChecker) == 0){
                        queryModal <- 'Query is not present in selected Interaction Database! Please try a different query...'
                        multipleMatches <- 'NoMatch'
                        return(list(queryMasterChecker, multipleMatches, queryModal))
                    }
                    else if(length(queryMasterChecker) == 1){
                        queryModal <- paste('Query is not present in selected Interaction Database! Proceede with close match: ', queryMasterChecker, '?', sep='') 
                        multipleMatches <- 'OneExpressionMatch'
                        QueryName <- queryMasterChecker
                    }
                }
                if(length(grep(QueryName, names(PredictedTargetsmiRNA))) == 1) {
                    queryMasterChecker <- intersect(names(PredictedTargetsmiRNA)[grep(QueryName, names(PredictedTargetsmiRNA))], rownames(mirnaFile))
                    
                    if(length(queryMasterChecker) == 1){
                        queryModal <- paste('Query is not present in selected Interaction Database! Proceede with close match: ', queryMasterChecker, '?', sep='')
                        QueryName <- queryMasterChecker
                        multipleMatches <- 'OneExpressionMatch'
                    }
                    else if(length(queryMasterChecker) == 0){
                        queryModal <- 'Query is not present in selected Interaction Database! Please try a different query...'
                        multipleMatches <- 'NoMatch'
                        return(list(queryMasterChecker, multipleMatches, queryModal))
                    }
                }
                else if(length(grep(QueryName, names(PredictedTargetsmiRNA))) == 0) {
                    tmpIDchecker <- FALSE
                }
            }
        }
    }
    
    if(tmpIDchecker == FALSE) {
        if(querySelection == 'miRNA'){
            tmpMismatch <- names(PredictedTargetsmiRNA)[agrep(QueryName, names(PredictedTargetsmiRNA), ignore.case=TRUE)]
            tmpMismatch <- intersect(tmpMismatch, rownames(mirnaFile))
        }
        else if(querySelection == 'Gene'){
            tmpMismatch <- names(PredictedTargetsGene)[agrep(QueryName, names(PredictedTargetsGene), ignore.case=TRUE)]
            tmpMismatch <- intersect(tmpMismatch, rownames(mrnaFile))
        }
        tmpMismatch <- tmpMismatch[1:5]
        tmpMismatch <- tmpMismatch[!is.na(tmpMismatch)]
        
        if(length(tmpMismatch) > 0){
            queryModal <- 'Query is not present in selected Interaction Database! Select from similar identifiers or try a different query. Here are a few similar matches ...'
            multipleMatches <- TRUE
            return(list(tmpMismatch, multipleMatches, queryModal))
        }
        else {
            queryModal <- 'Query is not present in selected Interaction Database! Please try a different query...'
            multipleMatches <- 'NoMatch'
            return(list(tmpMismatch, multipleMatches, queryModal))
        }
    }
    
    if(!is.null(PredictedTargetsGene[[QueryName]])) {
        ## Check whether query has a correspondence in the database.
        ## Query can be matched to transposed or normal table of targetscan database depending on whether mRNA or miRNA query is entered
        ExtractQueryFromMTIdb <- PredictedTargetsGene[[QueryName]]
        ExtractQueryFromMTIdb <- as.data.frame.matrix(cbind(QueryName, ExtractQueryFromMTIdb))
    }
    else if(!is.null(PredictedTargetsmiRNA[[QueryName]])) {
        ExtractQueryFromMTIdb <- PredictedTargetsmiRNA[[QueryName]]
        ExtractQueryFromMTIdb <- as.data.frame.matrix(cbind(QueryName, ExtractQueryFromMTIdb))
    }
    else{
        ExtractQueryFromMTIdb <- NULL
        if(querySelection == 'miRNA'){
            tmpMismatch <- intersect(names(PredictedTargetsmiRNA)[agrep(QueryName, names(PredictedTargetsmiRNA), ignore.case=TRUE)], rownames(mirnaFile))
        }
        else if(querySelection == 'Gene'){
            tmpMismatch <- intersect(names(PredictedTargetsGene)[agrep(QueryName, names(PredictedTargetsGene), ignore.case=TRUE)], rownames(mrnaFile))
        }
        
        tmpMismatch <- tmpMismatch[1:5]
        tmpMismatch <- tmpMismatch[!is.na(tmpMismatch)]
        querySelection <- "NotFound"
        
        if(length(tmpMismatch) > 0){
            queryModal <- 'No exact matches found! Select from similar identifiers or try a different query. Here are a few similar matches ...'
            multipleMatches <- TRUE
            return(list(tmpMismatch, multipleMatches, queryModal))
        }
        else {
            queryModal <- 'No matches found! Please try a different query...'
            multipleMatches <- 'NoMatch'
            return(list(tmpMismatch, multipleMatches, queryModal))
        }
    }
    
    if(miRNADatabase == "TargetScan&miRNet"){
        if(unionorintersect == "Union"){
            wholeNetworkIDs <- sapply(strsplit(as.matrix(ExtractQueryFromMTIdb[,2]),'/'),c)[1,]
        }
        else{
            wholeNetworkIDs <- as.matrix(ExtractQueryFromMTIdb[,2])[,1]
        }
    }
    else{
        wholeNetworkIDs <- as.matrix(ExtractQueryFromMTIdb[,2])[,1]
    }
    
    
    if(querySelection == 'miRNA'){
        wholeNetwork <- lapply(wholeNetworkIDs,function(x)PredictedTargetsGene[[x]])
    }
    else if(querySelection == 'Gene'){
        wholeNetwork <- lapply(wholeNetworkIDs,function(x)PredictedTargetsmiRNA[[x]])
    }
    
    names(wholeNetwork) <- wholeNetworkIDs
    wholeNetwork <- wholeNetwork[!sapply(wholeNetwork,is.null)]
    wholeNetwork <- lapply(names(wholeNetwork), function(x)cbind(x,wholeNetwork[[x]]))
    CytoscapeFile <- as.data.frame.matrix(do.call(rbind, wholeNetwork))
    
    if(querySelection == 'miRNA'){
        CytoscapeFile <- as.data.frame.matrix(data.frame(miRNA=CytoscapeFile[,2],Gene=CytoscapeFile[,1]))
        colnames(ExtractQueryFromMTIdb) <- c('miRNA', 'Gene')
    }
    else if(querySelection == 'Gene'){
        ExtractQueryFromMTIdb <- as.data.frame.matrix(data.frame(miRNA=ExtractQueryFromMTIdb[,2], Gene=ExtractQueryFromMTIdb[,1]))
        CytoscapeFile <- as.data.frame.matrix(data.frame(miRNA=CytoscapeFile[,1],Gene=CytoscapeFile[,2]))
    }
    
    if(miRNADatabase == "TargetScan&miRNet"){
        ## After this part of the script convertion of MTIdb into Cytoscape ready file format completed without the addition of Expression File into the prospective network...
        if(unionorintersect == "Union"){
            ## Cover if combination of databases are being used...
            if(querySelection == "miRNA"){
                dataOriginCyto.mid <- strsplit(as.matrix(CytoscapeFile[1]),"/")
                dataOriginCyto <- as.matrix(unlist(lapply(dataOriginCyto.mid,function(x)x[2])))
                CytoscapeFile[1] <-  as.matrix(unlist(lapply(dataOriginCyto.mid,function(x)x[1])))
                dataOriginTable.mid <- strsplit(as.matrix(ExtractQueryFromMTIdb[2]),"/")
                dataOriginTable <- as.matrix(unlist(lapply(dataOriginTable.mid,function(x)x[2])))
                ExtractQueryFromMTIdb[2] <- as.matrix(unlist(lapply(dataOriginTable.mid,function(x)x[1])))
                CytoscapeFile <- as.data.frame.matrix(cbind(CytoscapeFile, dataOriginCyto))
                colnames(CytoscapeFile) <- c('miRNA', 'Gene', 'DataOrigin')
                
                ExtractQueryFromMTIdb <- as.data.frame.matrix(cbind(ExtractQueryFromMTIdb, dataOriginTable))
                colnames(ExtractQueryFromMTIdb) <- c('miRNA', 'Gene', 'DataOrigin')
            }
            else if(querySelection == "Gene"){
                dataOriginCyto.mid <- strsplit(as.matrix(CytoscapeFile[2]),"/")
                dataOriginCyto <- as.matrix(unlist(lapply(dataOriginCyto.mid,function(x)x[2])))
                CytoscapeFile[2] <-  as.matrix(unlist(lapply(dataOriginCyto.mid,function(x)x[1])))
                dataOriginTable.mid <- strsplit(as.matrix(ExtractQueryFromMTIdb[1]),"/")
                dataOriginTable <- as.matrix(unlist(lapply(dataOriginTable.mid,function(x)x[2])))
                ExtractQueryFromMTIdb[1] <-  as.matrix(unlist(lapply(dataOriginTable.mid,function(x)x[1])))
                CytoscapeFile <- as.data.frame.matrix(cbind(CytoscapeFile, dataOriginCyto))
                
                colnames(CytoscapeFile) <- c('miRNA', 'Gene', 'DataOrigin')
                
                ExtractQueryFromMTIdb <- as.data.frame.matrix(cbind(ExtractQueryFromMTIdb, dataOriginTable))
                colnames(ExtractQueryFromMTIdb) <- c('miRNA', 'Gene', 'DataOrigin')
                
                colSwap <- as.matrix(ExtractQueryFromMTIdb[1])
                ExtractQueryFromMTIdb[1] <- as.matrix(ExtractQueryFromMTIdb[2])
                ExtractQueryFromMTIdb[2] <- colSwap
            }
        }
    }
    
    CytoscapeFile <- rbind(do.call(data.frame, ExtractQueryFromMTIdb), do.call(data.frame, CytoscapeFile))
    CytoscapeFile <- unique(as.data.frame.matrix(CytoscapeFile))
    
    
    if(mrnaFile != 'No Data'){
        if(nrow(merge(CytoscapeFile, mrnaFile2, by='Gene')) == 0){
            CytoscapeFile <- merge(CytoscapeFile, mrnaFile2, by='Gene', all.x = TRUE)
            CytoscapeFile$DiffmRNA <- 0
            mrnaFile2 <- 0
            completeDataUnmatch <- 'Gene'
        }
        else {
            CytoscapeFile <- merge(CytoscapeFile, mrnaFile2, by='Gene', all.x = TRUE)
        }
        nullGenes <- unique(CytoscapeFile[is.na(CytoscapeFile$DiffmRNA),]$Gene)
        CytoscapeFile <- as.data.frame.matrix(CytoscapeFile)
        CytoscapeFile[is.na(CytoscapeFile$DiffmRNA),]$DiffmRNA <- 0
        CytoscapeFile$DiffmRNA <- round(as.numeric(as.character(CytoscapeFile$DiffmRNA)), 3)
    }
    else{
        CytoscapeFile$DiffmRNA <- 0
        nullGenes <- NULL
    }
    
    if(mirnaFile != 'No Data'){
        if(nrow(merge(CytoscapeFile, mirnaFile2, by="miRNA")) == 0){
            CytoscapeFile <- merge(CytoscapeFile, mirnaFile2, by="miRNA", all.x = TRUE)
            CytoscapeFile$DiffmiRNA <- 0
            mirnaFile2 <- 0
            completeDataUnmatch <- 'miRNA'
        }
        else{
            CytoscapeFile <- merge(CytoscapeFile, mirnaFile2, by="miRNA", all.x = TRUE)
        }
        nullmiRNAs <- unique(CytoscapeFile[is.na(CytoscapeFile$DiffmiRNA),]$miRNA)
        CytoscapeFile <- as.data.frame.matrix(CytoscapeFile)
        CytoscapeFile[is.na(CytoscapeFile$DiffmiRNA),]$DiffmiRNA <- 0
        CytoscapeFile$DiffmiRNA <- round(as.numeric(as.character(CytoscapeFile$DiffmiRNA)), 3)
    }
    else{
        CytoscapeFile$DiffmiRNA <- 0
        nullmiRNAs <- NULL
    }
    
    CytoscapeFile <- as.data.frame.matrix(CytoscapeFile)
    
    if(querySelection == 'miRNA'){
        ## Store a primer query object for reloading if user decides to choose different filtering parameters later on.
        ## Used for increased performance, in other case reprocessing the files would have been required...
        PrimerQuery <- CytoscapeFile[CytoscapeFile$miRNA == QueryName,]
    }
    else if(querySelection == 'Gene'){
        PrimerQuery <- CytoscapeFile[CytoscapeFile$Gene == QueryName,]
    }
    
    if(SecondaryInteractions == FALSE){
        ## If comprehensive network option is not used, whole network is comprised of initial query user entered
        ## Network is not expanded to second degree interactions..
        CytoscapeFile <- PrimerQuery
    }
    
    ## Gene related GO terms are extracted here based on the query.
    ## Gene queries or miRNA queries both return a GO term table for the corresponding genes in the network...
    
    DiffmiRNAStatic <- range(as.numeric(as.character(CytoscapeFile$DiffmiRNA)))
    ## Differential expression value distribution is calculated here for max and min values for both miRNA expression file and mRNA expression file...
    DiffmRNAStatic <- range(as.numeric(as.character(CytoscapeFile$DiffmRNA)))
    quantilemRNA <- CytoscapeFile$DiffmRNA[CytoscapeFile$DiffmRNA != 0]
    quantilemiRNA <- CytoscapeFile$DiffmiRNA[CytoscapeFile$DiffmiRNA != 0]
    
    quantile00mRNA <- quantile(as.numeric(as.character(quantilemRNA)), 0)
    quantile05mRNA <- quantile(as.numeric(as.character(quantilemRNA)), 0.05)
    ## 5 percent and 95 percent quantiles are calculated for automated filtering selections suggested by Odyssey...
    quantile95mRNA <- quantile(as.numeric(as.character(quantilemRNA)), 0.95)
    quantile100mRNA <- quantile(as.numeric(as.character(quantilemRNA)), 1.0)
    quantile00miRNA <- quantile(as.numeric(as.character(quantilemiRNA)), 0)
    quantile05miRNA <- quantile(as.numeric(as.character(quantilemiRNA)), 0.05)
    quantile95miRNA <- quantile(as.numeric(as.character(quantilemiRNA)), 0.95)
    quantile100miRNA <- quantile(as.numeric(as.character(quantilemiRNA)), 1.0)
    quantileValues <- list(quantile00miRNA, quantile05miRNA, quantile95miRNA, quantile100miRNA, quantile00mRNA, quantile05mRNA, quantile95mRNA, quantile100mRNA)
    names(quantileValues) <- c('miRNA00', 'miRNA05', 'miRNA95', 'miRNA100', 'mRNA00', 'mRNA05', 'mRNA95', 'mRNA100')
    
    return(list(queryMasterChecker, multipleMatches, queryModal, PrimerQuery, CytoscapeFile, DiffmRNAStatic, DiffmiRNAStatic, mirnaFile2, mrnaFile2, quantileValues, completeDataUnmatch, nullmiRNAs, nullGenes))
    
}

#########################################################

## Filtering of the nodes based on expression difference is executed here.
## User can either proceed with suggested filtering intervals based on quantile calculations or can select to use a more or less strict filtering interval.

ProcessFilteringNetworkDGex <- function(CytoscapeFile, slidermiRNA, slidermRNA, InitQueryFilter, PrimerQuery){
    CytoscapeFileInitial <- CytoscapeFile
    
    if(sum(as.numeric(unique(CytoscapeFile$DiffmiRNA))) != 0){
        CytoscapeFile <- CytoscapeFile[which(as.numeric(as.character(CytoscapeFile$DiffmiRNA)) <= slidermiRNA[1] | as.numeric(as.character(CytoscapeFile$DiffmiRNA)) >= slidermiRNA[2]),]
    }
    if(sum(as.numeric(unique(CytoscapeFile$DiffmRNA))) != 0){
        CytoscapeFile <- CytoscapeFile[which(as.numeric(as.character(CytoscapeFile$DiffmRNA)) <= slidermRNA[1] | as.numeric(as.character(CytoscapeFile$DiffmRNA)) >= slidermRNA[2]),]
    }
    
    ## Firstly resulting network can be selected to stay intact or can be added to filtering process.
    ## This is useful because sometimes filtering process filters out most of the query related nodes and leaves the user with a network not interesting for conducted particular analysis...
    if(sum(as.numeric(unique(CytoscapeFile$DiffmiRNA))) != 0){
        PrimerQueryFiltered <- PrimerQuery[which(as.numeric(as.character(PrimerQuery$DiffmiRNA)) <= slidermiRNA[1] | as.numeric(as.character(PrimerQuery$DiffmiRNA)) >= slidermiRNA[2]),]
    }
    else{
        PrimerQueryFiltered <- PrimerQuery
    }
    if(sum(as.numeric(unique(CytoscapeFile$DiffmRNA))) != 0){
        PrimerQueryFiltered <- PrimerQueryFiltered[which(as.numeric(as.character(PrimerQueryFiltered$DiffmRNA)) <= slidermRNA[1] | as.numeric(as.character(PrimerQueryFiltered$DiffmRNA)) >= slidermRNA[2]),]
    }
    
    if(InitQueryFilter == TRUE) {
        CytoscapeFile <- rbind(CytoscapeFile, PrimerQueryFiltered)
    }
    else {
        CytoscapeFile <- rbind(CytoscapeFile, PrimerQuery)
    }
    
    CytoscapeFile <- as.data.frame.matrix(unique(CytoscapeFile))
    return(list(CytoscapeFile, CytoscapeFileInitial))
}
#########################################################

ProcessFilteringNetworkDegree <- function(CytoscapeFile, degreeFilterUI, degreeFilterValue){
    CytoDegree <- CytoscapeFile
    
    miRNACounts <- as.data.frame(table(CytoscapeFile$miRNA))
    colnames(miRNACounts) <- c('miRNA', 'miRNACount')
    geneCounts <- as.data.frame(table(CytoscapeFile$Gene))
    colnames(geneCounts) <- c('Gene', 'geneCount')
    CytoscapeFile <- merge(CytoscapeFile, miRNACounts, by='miRNA')
    CytoscapeFile <- merge(CytoscapeFile, geneCounts, by='Gene')
    sliderLenGene <- length(CytoscapeFile$geneCount)
    sliderLenmiRNA <- length(CytoscapeFile$miRNACount)
    sliderParameters <- c(CytoscapeFile$geneCount, CytoscapeFile$miRNACount)
    names(sliderParameters) <- c(rep('Gene', sliderLenGene), rep('miRNA', sliderLenmiRNA))
    sliderParameters <- sliderParameters[!grepl(names(which.max(sliderParameters)), x = names(sliderParameters))]
    sliderParameters <- unique(sliderParameters)
    sliderParameters <- sliderParameters[sliderParameters < max(sliderParameters)]
    sliderParametersMax <- max(sliderParameters)
    sliderParametersMin <- degreeFilterValue
    
    if(degreeFilterUI == TRUE){
        CytoscapeFile <- CytoscapeFile[CytoscapeFile$miRNACount > degreeFilterValue,]
        CytoscapeFile <- CytoscapeFile[CytoscapeFile$geneCount > degreeFilterValue,]
    }
    else{
        CytoscapeFile <- CytoDegree
    }
    
    if(is.null(CytoscapeFile$DataOrigin)){
        CytoscapeFile <- CytoscapeFile[,-c(5,6)]
    }
    else {
        CytoscapeFile <- CytoscapeFile[,-c(6,7)]
    }
    
    miRNACounts <- as.data.frame(table(CytoscapeFile$miRNA))
    if(nrow(miRNACounts) > 0){
        colnames(miRNACounts) <- c('miRNA', 'miRNACount')
    }
    geneCounts <- as.data.frame(table(CytoscapeFile$Gene))
    if(nrow(geneCounts) > 0){
        colnames(geneCounts) <- c('Gene', 'geneCount')
    }
    if(nrow(miRNACounts) > 0 & nrow(geneCounts) > 0){
        CytoscapeFile <- merge(CytoscapeFile, miRNACounts, by='miRNA')
        CytoscapeFile <- merge(CytoscapeFile, geneCounts, by='Gene')
        return(list(CytoscapeFile, sliderParametersMax, sliderParametersMin))
    }
}

########################################################################################################


#Title: KNGP Algorithm
#Author: Chad Kimmel

#Summary: The following is the code for the KNGP Algorithm.  The user must input 4 sources of input.  
#Go to bottom of this script for an explanation of the inputs needed by the user and the main function call.  

#These packages should be installed
require(foreach) #install.packages("foreach")
require(pROC) #install.packages("pROC")

#The main KGNPAlgorithm Call
KNGPAlgorithm <- function(LinkKnowledgeMatrix,NodeKnowledgeVector,NodeNamesVector,rootNodes)
{
    #Obtain the transition probability matrix. 
    TransitionProbMatrix <- getTransitionProbMatrix(LinkKnowledgeMatrix,apply(LinkKnowledgeMatrix,1,sum,na.rm=TRUE))
    rm(LinkKnowledgeMatrix); gc()
    
    #Obtain the root node vector index
    rootNodeVectorIndex <- foreach(i=1:length(rootNodes),.combine=c) %do% {which(NodeNamesVector==rootNodes[i])}
    
    #List of f values to test.  
    f_vector <- c(0,0.5,1,2,5,10,50,100,1000,10000)
    
    best_f <- find_best_f(TransitionProbMatrix,NodeKnowledgeVector,rootNodeVectorIndex,f_vector)
    
    #Compute the prior prob vector using best_f
    priorProbVector <- getPriorProbVector(rootNodeVectorIndex,NodeKnowledgeVector,best_f)
    
    #Obtain the posterior prob vector.  
    posteriorProbVector <- inference(TransitionProbMatrix,matrix(priorProbVector,nrow=length(priorProbVector),ncol=1))
    
    #Pass back data frame
    result <- data.frame(cbind(name=NodeNamesVector,postProb=posteriorProbVector))
    names(result) <- c("name","postProb")
    result$name <- as.character(result$name)
    result$postProb <- as.numeric(as.character(result$postProb))
    result$name <- ifelse(result$name %in% rootNodes,paste0(result$name,"*"),result$name)
    return(list(KNGPResult=result,best_f=best_f))
    
}

#The find_best_f procedure.  
find_best_f <- function(TransitionProbMatrix,NodeKnowledgeVector,rootNodeVectorIndex,f_vector)
{
    best_f <- NULL
    best_AUC <- -100000
    
    nonRootNodes <- setdiff(1:dim(TransitionProbMatrix)[1],rootNodeVectorIndex) #Get the vector of all non-root nodes
    
    for (f in f_vector)
    {
        f_result <- foreach(i=1:length(rootNodeVectorIndex), .combine=rbind) %do%
            {
                rootIndex <- rootNodeVectorIndex[i] #The root index for the given loop
                proteinsToRank <- c(rootIndex, sample(nonRootNodes,getSampleSize(length(nonRootNodes))) ) #Create the list of proteins to rank.
                
                newRootNodes <- setdiff(rootNodeVectorIndex,rootIndex) #The new root nodes minus the rootIndex
                
                priorProbVector <- getPriorProbVector(newRootNodes,NodeKnowledgeVector,f) #Compute the prior prob vector
                posteriorProbVector <- inference(TransitionProbMatrix,matrix(priorProbVector,nrow=length(priorProbVector),ncol=1)) #Compute the posterior prob vector.
                data.frame(cbind(score=posteriorProbVector[,1][proteinsToRank],actual=c(1,rep(0,length(proteinsToRank)-1))))
            }
        
        auc <- as.numeric(auc(f_result$actual,f_result$score, quiet=TRUE))
        if (auc>best_AUC) {best_AUC=auc;best_f=f}
        
    }
    
    return(best_f)
    
}

#The PageRank inference procedure.  
inference <- function(TransitionProbMatrix,priorProbMatrix)
{
    B = 0.5
    threshold = 0.00001
    delta = 1
    Po_prev = matrix(rep(0,dim(priorProbMatrix)[1]),nrow=dim(priorProbMatrix)[1],ncol=1)
    i = 1
    
    while (delta>threshold)
    {
        Po_curr <- ((1-B)*(TransitionProbMatrix %*% Po_prev)) + (B*priorProbMatrix)
        delta <- abs(sum(Po_curr) - sum(Po_prev))
        Po_prev <- Po_curr
        
    }
    
    return(Po_prev)
}



#Get the transition probability matrix
getTransitionProbMatrix <- function(TransitionProbMatrix,sumMatrix)
{
    
    for (c in 1:dim(TransitionProbMatrix)[1])
    {
        if (sumMatrix[c]!=0) TransitionProbMatrix[,c] <- TransitionProbMatrix[,c]/sumMatrix[c] else
            TransitionProbMatrix[,c] <- rep(0,dim(TransitionProbMatrix)[1])  
    }
    
    return(TransitionProbMatrix)
}

#Get the prior probability matrix.  
getPriorProbVector <- function(rootNodes, NodeKnowledgeVector, f)
{
    newPriorValueVector <- foreach(s=1:length(NodeKnowledgeVector), .combine=c) %do% 
        {
            if (s %in% rootNodes) NodeKnowledgeVector[s]*f else
                NodeKnowledgeVector[s]
        }
    return(newPriorValueVector/sum(newPriorValueVector))
}

#Get the sample size to pull from.
getSampleSize <- function(size)
{
    if (size>=99) return(99) else
        return(size)
}

#########################################################

ProcessvisNetworkShiny <- function(CytoscapeFile, miRNADatabase, unionorintersect, sliderColorSpectrum, nullmiRNAs, nullGenes) {
    ## Network is created here in visual based on Cytoscape file prepared in the previous functions.
    
    ids.mtx <- c(CytoscapeFile$miRNA ,CytoscapeFile$Gene)
    type.mtx <- c(rep('miRNA', length(CytoscapeFile$miRNA)), rep('Gene', length(CytoscapeFile$Gene)))
    tbl.nodes <- unique(data.frame(id=ids.mtx, group=type.mtx))
    tbl.nodes$color <- rep('NA', nrow(tbl.nodes))
    ##  tbl.nodes$shape <- rep('NA', nrow(tbl.nodes))
    
    source.mtx <- CytoscapeFile$miRNA
    target.mtx <- CytoscapeFile$Gene
    
    if(length(source.mtx) == 0){
        return(list(NA, NA))    
    }
    
    tbl.edges <- unique(data.frame(from=source.mtx, to=target.mtx, arrows='to'))
    
    if(miRNADatabase == "TargetScan&miRNet"){
        if(unionorintersect == "Union"){
            ## Different edge colors are given based on Data Origin.
            tblSave <<- tbl.edges
            cytoSave <<- CytoscapeFile
            tbl.edges$color <- rep('#00FF00', nrow(tbl.edges))
            tbl.edges$color[which(grepl('TARGETSCAN', CytoscapeFile$DataOrigin))] <- 'rgb(0, 0, 170)'
            tbl.edges$color[which(grepl('MIRNET',  CytoscapeFile$DataOrigin))] <- 'rgb(221, 0, 0)'
            tbl.edges$color[which(grepl('BOTH',  CytoscapeFile$DataOrigin))] <- 'rgb(0, 200, 0)'
        }
    }
    
    if(!is.null(CytoscapeFile$DiffmiRNA)){
        ## Node color column is set based on expression value difference in the expression files for miRNA and gene ..
        miRNA_factor <- factor(as.character(as.matrix(CytoscapeFile$miRNA)))
        tbl.nodes$color[grep('HSA-', tbl.nodes$id)] <- as.numeric(CytoscapeFile$DiffmiRNA[match(tbl.nodes$id[grep('HSA-', tbl.nodes$id)], miRNA_factor)])
    }
    if(!is.null(CytoscapeFile$DiffmRNA)){
        mRNAs <- tbl.nodes$id[!grepl('HSA-', tbl.nodes$id)]
        mRNA_factor <- factor(as.character(as.matrix(CytoscapeFile$Gene)))
        tbl.nodes$color[!grepl('HSA-', tbl.nodes$id)] <- as.numeric(CytoscapeFile$DiffmRNA[match(mRNAs, mRNA_factor)])
    }
    if(!is.null(CytoscapeFile$DiffmiRNA) || !is.null(CytoscapeFile$DiffmRNA)){
        ## Sometimes nodes can have below minimum or above maximum differential expression values for intervals set for node colors.
        ## These nodes are colored by border colors within the color spectrum chosen.
        ## Color spectrum spreads from shades of green to shades of red.
        
        tmp.nodes1 <- tbl.nodes[tbl.nodes$id %in% nullmiRNAs,]
        tmp.nodes2 <- tbl.nodes[tbl.nodes$id %in% nullGenes,]
        tmp.nodes <- rbind(tmp.nodes1, tmp.nodes2)
        rm(tmp.nodes1)
        rm(tmp.nodes2)
        
        tbl.nodes <- tbl.nodes[!tbl.nodes$id %in% tmp.nodes$id,]
        if(nrow(tmp.nodes) > 0){
            tmp.nodes$color <- '#BFEFFF'
            tmp.nodes$group <- paste("Nodes", "\n", "without", "\n", "Expression", "\n", "Data", sep="", collapse="")
        }
        
        lowExprs <- which(as.numeric(tbl.nodes$color) < sliderColorSpectrum[1])
        highExprs <- which(as.numeric(tbl.nodes$color) > sliderColorSpectrum[2])
        
        temp <- cut(as.numeric(tbl.nodes$color), breaks = seq(sliderColorSpectrum[1], sliderColorSpectrum[2], length.out = 100), include.lowest = TRUE, ordered_result = FALSE)
        colorsForGraph <- colorRampPalette(c('green', 'grey', 'red'))(99)[temp]
        ## Nodes that are not differentially expressed in control vs treatment samples are colored grey
        tbl.nodes$color <- colorsForGraph
        tbl.nodes$group <- NA
        tbl.nodes[grep('HSA-', tbl.nodes$id),]$group <- 'miRNA'
        tbl.nodes[!grepl('HSA-', tbl.nodes$id),]$group <- 'Gene'
        tbl.nodes <- rbind(tbl.nodes, tmp.nodes)
        rm(tmp.nodes)
        
        if(length(lowExprs) > 0){
            tbl.nodes[lowExprs,]$color <- '#008000'
        }
        if(length(highExprs) > 0){
            tbl.nodes[highExprs,]$color <- '#FF0000'
        }
        
        if(!is.null(CytoscapeFile$DiffmiRNA) & !is.null(CytoscapeFile$DiffmRNA)){
            CytoscapeFile$Correlate <- as.numeric(as.character(CytoscapeFile$DiffmiRNA)) * as.numeric(as.character(CytoscapeFile$DiffmRNA))
            
            negCorrList = NA
            posCorrList = NA
            for(i in 1:nrow(CytoscapeFile)){
                if(CytoscapeFile[i,]$Correlate > 0){
                    posCorrList <- c(posCorrList, CytoscapeFile[i,]$miRNA, CytoscapeFile[i,]$Gene)
                }
                else{
                    negCorrList <- c(negCorrList, CytoscapeFile[i,]$miRNA, CytoscapeFile[i,]$Gene)
                }
            }
            
            posCorrList <- posCorrList[complete.cases(posCorrList)]
            negCorrList <- negCorrList[complete.cases(negCorrList)]
            totalCorrList <- c(posCorrList, negCorrList)
            corrList <- setdiff(posCorrList, negCorrList)
            tbl.nodes$colorSave <- tbl.nodes$color
            tbl.nodes$colorCorr <- tbl.nodes$color
            try(tbl.nodes[tbl.nodes$id %in% corrList,]$colorCorr <- '#EAE9E9', silent = TRUE)
        }
        
    }
    else {
        ## If there are no expression file chosen or uploaded, all nodes are colored with designated color.
        tbl.nodes$color[grep('HSA-', tbl.nodes$id)] <- 'rgb(255,102,102)'
        tbl.nodes$color[as.numeric(setdiff(tbl.nodes$id, grep('HSA-', tbl.nodes$id)))] <- 'rgb(118,33,192)'
    }
    
    graph <- graph.data.frame(tbl.edges, directed = T)
    degree_value <- igraph::degree(graph, mode = "all")
    tbl.nodes$value <- log2(degree_value[match(tbl.nodes$id, names(degree_value))]) + 100
    tbl.nodes$label <- tbl.nodes$id
    tbl.nodes <- unique(tbl.nodes)
    
    totalCorr.df <- as.data.frame(plyr::count(totalCorrList))
    colnames(totalCorr.df) <- c('id', 'totalCounts')
    negCorr.df <- as.data.frame(plyr::count(negCorrList))
    colnames(negCorr.df) <- c('id', 'negCounts')
    posCorr.df <- as.data.frame(plyr::count(posCorrList))
    colnames(posCorr.df) <- c('id', 'posCounts')
    
    tbl.nodes <- merge(tbl.nodes, totalCorr.df, by='id', all.x = TRUE)
    tbl.nodes <- merge(tbl.nodes, negCorr.df, by='id', all.x = TRUE)
    
    if(sum(is.na(tbl.nodes$negCounts)) > 0){
        tbl.nodes[is.na(tbl.nodes$negCounts),]$negCounts <- 0
    }
    tbl.nodes$totalCounts <- (tbl.nodes$negCounts / tbl.nodes$totalCounts) * 100
    tbl.nodes <- tbl.nodes[,-9]
    
    tbl.nodes$label <- paste(tbl.nodes$id, ' (', round(tbl.nodes$totalCounts), '%)', sep='')
    tbl.nodes$title <- tbl.nodes$label
    
    return(list(tbl.nodes, tbl.edges, negCorr.df, posCorr.df, totalCorr.df))
    
}


###############################################################################################

# load(file='GenePlatformAnnotations.RData')
# load(file='miRNA.id_hs.RData')
# 
# addExampleData <- function(coords = getwd(), miRNAPlatforms = c('GPL8786', 'GPL21572'), mRNAPlatforms = c('GPL570', 'GPL96', 'GPL97', 'GPL6244', 'GPL17692')){
#     
#     exampleDataList <- list.files()
#     miRNAFiles <- exampleDataList[unlist(sapply(miRNAPlatforms, function(x){grep(exampleDataList, pattern = x, fixed = TRUE)}))]
#     mRNAFiles <- exampleDataList[unlist(sapply(mRNAPlatforms, function(x){grep(exampleDataList, pattern = x, fixed = TRUE)}))]
#     exampleDB.mirnas <- list()
#     exampleDB.mrnas <- list()
# 
#     for(i in 1:length(miRNAFiles)){
# 
#         singlemirnaFile <- read.table(miRNAFiles[i], sep='\t', header=TRUE,row.names=1, comment.char='!')
#         singlemirnaFile <- singlemirnaFile[complete.cases(singlemirnaFile),]
# 
#         range.singlemirnaFile <- range(singlemirnaFile)
#         if(range.singlemirnaFile[2] > 20){
#             singlemirnaFile <- log2(singlemirnaFile)
#         }
# 
#         dimensionsExpmiRNA <- nrow(singlemirnaFile)
#         singlemirnaFile <- cbind(rownames(singlemirnaFile), singlemirnaFile)
#         colnames(singlemirnaFile)[1] <- 'input_id'
#         mirnaFile.annotated <- merge(singlemirnaFile, miRNA.identifiers, by='input_id')
#         if(nrow(mirnaFile.annotated) / dimensionsExpmiRNA * 100 >= 5){
#             singlemirnaFile <- mirnaFile.annotated
#             singlemirnaFile <- singlemirnaFile[,-1]
#             nSamples <- ncol(singlemirnaFile)-1
#             ## singlemirnaFile <- singlemirnaFile[,-(nSamples-1)]
# 
#             if(sum(duplicated(singlemirnaFile$v22)) > 0){
#                 singlemirnaFile <- aggregate(.~v22, data = singlemirnaFile, mean)
#                 rownames(singlemirnaFile) <- singlemirnaFile$v22
#                 singlemirnaFile <- singlemirnaFile[,-1]
#             }
#             else {
#                 rownames(singlemirnaFile) <- singlemirnaFile$v22
#                 singlemirnaFile <- singlemirnaFile[,-(nSamples+1)]
#             }
#             tempObject <- unlist(strsplit(miRNAFiles[i], '-'))[[1]]
#             exampleDB.mirnas[[i]] <- singlemirnaFile
#             names(exampleDB.mirnas)[i] <- tempObject
#         }
#         else {
#             singlemirnaFile <- cbind(rownames(singlemirnaFile), singlemirnaFile)
#             colnames(singlemirnaFile)[1] <- 'input_id'
#             singlemirnaFile$input_id <- gsub(pattern = '_st', replacement = '', singlemirnaFile$input_id)
#             singlemirnaFile$input_id <- gsub(pattern = '-star', replacement = '*', singlemirnaFile$input_id)
#             mirnaFile.annotated <- merge(singlemirnaFile, miRNA.identifiers, by='input_id')
#             if(nrow(mirnaFile.annotated) / dimensionsExpmiRNA * 100 >= 5){
#                 singlemirnaFile <- mirnaFile.annotated
#                 singlemirnaFile <- singlemirnaFile[,-1]
#                 if(sum(duplicated(singlemirnaFile$v22)) > 0){
#                     singlemirnaFile <- aggregate(.~v22, data = singlemirnaFile, mean)
#                     rownames(singlemirnaFile) <- singlemirnaFile$v22
#                     singlemirnaFile <- singlemirnaFile[,-1]
#                 }
#                 else {
#                     nSamples <- ncol(singlemirnaFile)
#                     rownames(singlemirnaFile) <- singlemirnaFile$v22
#                     singlemirnaFile <- singlemirnaFile[,-nSamples]
#                 }
#                 tempObject <- unlist(strsplit(miRNAFiles[i], '-'))[[1]]
#                 exampleDB.mirnas[[i]] <- singlemirnaFile
#                 names(exampleDB.mirnas)[i] <- tempObject
#             }
#             else{
#                 datacheck.mrna == 'Data is not formatted correctly!'
#                 singlemirnaFile <- 'No Data'
#                 print('Passed on miRNA Data')
#             }
#         }
#     }
# 
#     for(i in 1:length(mRNAFiles)){
# 
#         singlemrnaFile <- read.table(mRNAFiles[i], sep='\t', header=TRUE,row.names=1, comment.char='!')
#         singlemrnaFile <- singlemrnaFile[complete.cases(singlemrnaFile),]
#         range.singlemrnaFile <- range(singlemrnaFile)
#         if(range.singlemrnaFile[2] > 20){
#             singlemrnaFile <- log2(singlemrnaFile)
#         }
#         dimensionsExpmRNA <- nrow(singlemrnaFile)
#         mrnaFile.annotated <- parseGSEmRNA(singlemrnaFile, GenePlatformAnnotations)
#         if(nrow(mrnaFile.annotated) / dimensionsExpmRNA * 100 >= 5){
#             tempObject2 <- unlist(strsplit(mRNAFiles[i], '-'))[[1]]
#             exampleDB.mrnas[[i]] <- mrnaFile.annotated
#             names(exampleDB.mrnas)[i] <- tempObject2
#         }
#         else{
#             print('Passed on mRNA Data')
#         }
# 
#     }
#     return(list(exampleDB.mirnas, exampleDB.mrnas))
# }


##########################################################################

applyPriorityOnNetwork <- function(nodeData, KNGPData) {
    
    colnames(nodeData)[1] <- 'id'
    colnames(KNGPData)[1] <- 'id'
    KGNPData <- as.data.frame.matrix(KNGPData)
    KNGPData <- KNGPData[,-3]
    KNGPData$id <- gsub(KNGPData$id, pattern = '*', replacement = '', fixed=TRUE)
    tmpNodeDT <- merge(nodeData, KNGPData, by='id')
    tmpNodeDT$value <- tmpNodeDT$Involvement_Score * nrow(tmpNodeDT) * 20
    tmpNodeDT <- tmpNodeDT[c(-6,-7)]
    tmpNodeDT$group <- 'Low'
    quantile90Nodes <- quantile(as.numeric(as.character(tmpNodeDT$value)), 0.90)
    quantile75Nodes <- quantile(as.numeric(as.character(tmpNodeDT$value)), 0.75)
    tmpNodeDT[tmpNodeDT$value >= quantile75Nodes,]$group <- 'Medium - 25%' 
    tmpNodeDT[tmpNodeDT$value >= quantile90Nodes,]$group <- 'High - 10%'
    tmpNodeDT <- tmpNodeDT[match(nodeData$id, tmpNodeDT$id),]
    
    return(tmpNodeDT)
}


##########################################################################

# library(GEOquery)
# miRNApdDB <- list()
# mRNApdDB <- list()
# exampleDataVector <- names(miRNAexprsDB)
# 
# phenoDataOdyssey <- function(exampleDataVector){
#     
#     for(i in 1:length(exampleDataVector)){ 
#         tryCatch(tmpGSE <- getGEO(exampleDataVector[i], GSEMatrix=TRUE), error=function(e)NA)
#         if(!is.na(tmpGSE)){
#             if(exampleDataVector[i] == 'GSE40321'){
#                 tmppdmRNA <- pData(phenoData(tmpGSE[[2]]))
#                 tmppdmiRNA <- pData(phenoData(tmpGSE[[3]]))
#                 
#                 mRNApdDB[[exampleDataVector[i]]] <- tmppdmRNA
#                 miRNApdDB[[exampleDataVector[i]]] <- tmppdmiRNA
#             }
#             else{
#                 
#                 tmppdmRNA <- pData(phenoData(tmpGSE[[1]]))
#                 tmppdmiRNA <- pData(phenoData(tmpGSE[[2]]))
#                 
#                 mRNApdDB[[exampleDataVector[i]]] <- tmppdmRNA
#                 miRNApdDB[[exampleDataVector[i]]] <- tmppdmiRNA
#                 
#             }
#         }
#         else{
#             print(paste('Encountered an error with ', exampleDataVector[i], ' Skipping to the next dataset', sep=''))    
#         }
#         
#         print(paste('Completed ', i, ' of ', length(exampleDataVector), sep=''))
#     }
#     
#     return(list(mRNApdDB, miRNApdDB))
# }



#########################################################################################

layoutNameConverter <- function(layoutInput){
    
    if(layoutInput == "Star-Shape"){
        return("layout_as_star")
    }
    else if(layoutInput == 'Reingold-Tilford'){
        return("layout_as_tree")
    }
    else if(layoutInput == 'Vertices On A Circle'){
        return("layout_in_circle")
    }
    else if(layoutInput == 'Nicely'){
        return("layout_nicely")
    }
    else if(layoutInput == 'On Grid'){
        return("layout_on_grid")
    }
    else if(layoutInput == 'On Sphere'){
        return("layout_on_sphere")
    }
    else if(layoutInput == 'Randomly'){
        return("layout_randomly")
    }
    else if(layoutInput == 'Davidson-Harel'){
        return("layout_with_dh")
    }
    else if(layoutInput == 'Fruchterman-Reingold'){
        return("layout_with_fr")
    }
    else if(layoutInput == 'GEM force-directed'){
        return("layout_with_gem")    
    }
    else if(layoutInput == 'Graphopt'){
        return("layout_with_graphopt")
    }
    else if(layoutInput == 'Kamada-Kawai'){
        return("layout_with_kk")
    }
    else if(layoutInput == 'Large Graph Layout'){
        return("layout_with_lgl")
    }
    else if(layoutInput == 'Multidimensional Scaling'){
        return("layout_with_mds")
    }
    else if(layoutInput == 'Sugiyama'){
        return("layout_with_sugiyama")
    }
    
}

#########################################################################################

round_df <- function(x, digits = 3) {
    # round all numeric variables
    # x: data frame 
    # digits: number of digits to round
    numeric_columns <- sapply(x, mode) == 'numeric'
    x[numeric_columns] <-  round(x[numeric_columns], digits)
    return(x)
}

#########################################################################################

sliderInput2 <- function(inputId, label, min, max, value, step=NULL, from_min, from_max){
    x <- sliderInput(inputId, label, min, max, value, step)
    x$children[[2]]$attribs <- c(x$children[[2]]$attribs, 
                                 "data-from-min" = from_min, 
                                 "data-from-max" = from_max, 
                                 "data-from-shadow" = TRUE)
    x
}

#### toupper fn ####

# miRNADBs <- list.files()
# 
# for(i in 1:length(miRNADBs)){
# 
#     load(miRNADBs[i])
#     if(i %% 2 == 0){
#         names(PredictedTargetsmiRNA) <- toupper(names(PredictedTargetsmiRNA))
#         PredictedTargetsmiRNA <- lapply(PredictedTargetsmiRNA, function(x){unique(toupper(x))})
#         save(PredictedTargetsmiRNA, file = miRNADBs[i])
#     }
#     else{
#         PredictedTargetsGene <- lapply(PredictedTargetsGene, function(x){unique(toupper(x))})
#         names(PredictedTargetsGene) <- toupper(names(PredictedTargetsGene))
#         save(PredictedTargetsGene, file = miRNADBs[i])
#     }
# }
