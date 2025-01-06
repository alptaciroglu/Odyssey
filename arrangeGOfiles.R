setwd('GO_annotations/')

gofiles <- list.files(pattern='RData')

for(i in 1:length(gofiles)){
    
    load(gofiles[i])
    if(length(grep('GOtomiRNA',gofiles[i])) == 1){
      
        GOtomiRNA <- sampleJSON
        save(GOtomiRNA, file=gofiles[i])
      
    }
    
    else{
      
        miRNAtoGO <- sampleJSON
        save(miRNAtoGO, file=gofiles[i])
      
    }
  
}

goList = list()
geneList = list()

for(i in 1:nrow(GenetoGO.ulti)){
  
  if(as.character(unlist(GenetoGO.ulti[i,1])) %in% names(goList)){
    goList[[as.character(unlist(GenetoGO.ulti[i,1]))]] <- unique(c(goList[[as.character(unlist(GenetoGO.ulti[i,1]))]], as.character(unlist(GenetoGO.ulti[i,2]))))
  }
  else {
    goList[[as.character(unlist(GenetoGO.ulti[i,1]))]] <- as.character(unlist(GenetoGO.ulti[i,2]))
  }
  
  
  if(as.character(unlist(GenetoGO.ulti[i,2])) %in% names(geneList)){
    geneList[[as.character(unlist(GenetoGO.ulti[i,2]))]] <- unique(c(geneList[[as.character(unlist(GenetoGO.ulti[i,2]))]], as.character(unlist(GenetoGO.ulti[i,1]))))
  }
  else {
    geneList[[as.character(unlist(GenetoGO.ulti[i,2]))]] <- as.character(unlist(GenetoGO.ulti[i,1]))
  }
  
  
  if(100%%(i/nrow(GenetoGO.ulti)*100) == 0){
    print(paste(i/nrow(GenetoGO.ulti)*100, ' is completed!', sep=''))
  }
}
