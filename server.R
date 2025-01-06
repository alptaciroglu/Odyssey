
source('functionsArsenal.R', local=TRUE)
## Source external scripts in Shiny server so they can be used

shinyServer(function(input, output, session) {
  
  hideTab(inputId = "sidetabs", target = 'NetworkOptions')
  hideTab(inputId = "tabs", target = 'panelExprs')
  hideTab(inputId = "tabs", target = 'panelCyto')
  hideTab(inputId = "tabs", target = 'panelNetwork')
  hideTab(inputId = "tabs", target = 'panelKNGP')
  hideTab(inputId = "sidetabs", target = 'InfoBox')
  hideTab(inputId = "tabs", target = "panelDGEx")
  hideTab(inputId = "sidetabs", target = 'SampleSelection')
  hideTab(inputId = "sidetabs", target = 'QuerySelection')
  shinyjs::hide('networkArrangement')
  
  floor_dec <- function(x, level=1) round(x - 5*10^(-level-1), level)
  ## In-house functions to be used in filtering
  ceiling_dec <- function(x, level=1) round(x + 5*10^(-level-1), level)
  
  network.data <- reactiveValues(multipleMatch = FALSE, removeUIcounter=0, removeNoMatchUI=0, removeExpressionNoMatch=0, ExpressionFile = NULL ,CytoscapeFile = NULL, NetworkFile = NULL, index = NULL, mrnaFile = NULL, mirnaFile = NULL,  noExprs = data.table('No Data were found'),  nodeList = NULL, networkNodes = NULL, networkEdges = NULL, TreatmentsmiRNAColumns = 'None', ControlsmiRNAColumns = 'None', TreatmentsmRNAColumns = 'None', ControlsmRNAColumns = 'None', sliderParameters = 1, miRNAUploadChecker = NULL, mRNAUploadChecker = NULL, KNGPModalBool = FALSE, prioritizationModalBool = FALSE, iterator = 0, mRNAphenoDB = mRNAphenoDB, miRNAphenoDB = miRNAphenoDB, corrSliderPre = '')
  ## Major collection of reactive values..
  ## It should be noted that more reactive application gets it will consume more resources.
  
  
  ################# Simple Action Button Collection ################
  ##################################################################
  
  observeEvent(input$BacktoSidePanel,{
    
    ## Change user visible interface panels according to analysis process
    showTab(inputId = "sidetabs", target = 'DataSelection')
    hideTab(inputId = "sidetabs", target = 'SampleSelection')
    hideTab(inputId = "tabs", target = 'panelExprs')
    updateTabsetPanel(session, "sidetabs", selected = "DataSelection")
    network.data$TreatmentsmiRNAColumns <- 'None' 
    network.data$ControlsmiRNAColumns <- 'None' 
    network.data$TreatmentsmRNAColumns <- 'None' 
    network.data$ControlsmRNAColumns <- 'None'
    
  })
  
  observeEvent(input$BacktoExprsPanel,{
    
    ## Change visible user interface panels according to analysis process
    showTab(inputId = "tabs", target = 'panelExprs')
    showTab(inputId = "sidetabs", target = 'SampleSelection')
    hideTab(inputId = "sidetabs", target = 'QuerySelection')
    hideTab(inputId = "tabs", target = 'panelDGEx')
    updateTabsetPanel(session, "sidetabs", selected = "SampleSelection")
    updateTabsetPanel(session, "tabs", selected = "panelExprs")
    
    if(network.data$mirnaFile == 'No Data' && network.data$mrnaFile == 'No Data'){
      click('BacktoSidePanel')
    }
  })
  
  observeEvent(input$BacktoDGExPanel,{
    
    ## Change visible user interface panels according to analysis process
    hideTab(inputId = "tabs", target = 'panelCyto')
    hideTab(inputId = "sidetabs", target = 'NetworkOptions')
    showTab(inputId = "tabs", target = "panelDGEx")
    showTab(inputId = "sidetabs", target = "QuerySelection")
    updateTabsetPanel(session, "sidetabs", selected = "QuerySelection")
    updateTabsetPanel(session, "tabs", selected = "panelDGEx")
    
  })
  
  observeEvent(input$BacktoCyto,{
    
    ## Change visible user interface panels according to analysis process
    showTab(inputId = "tabs", target = 'panelCyto')
    showTab(inputId = "sidetabs", target = 'NetworkOptions')
    hideTab(inputId = "tabs", target = 'panelNetwork')
    hideTab(inputId = "tabs", target = 'panelKNGP')
    hideTab(inputId = "sidetabs", target = 'InfoBox')
    updateTabsetPanel(session, "sidetabs", selected = "NetworkOptions")
    updateTabsetPanel(session, "tabs", selected = "panelCyto")
    shinyjs::hideElement(id= "networkLayoutPanel")
  })
  
  observeEvent(input$BacktoNetwork,{
    ## Change user visible interface panels according to analysis process
    updateTabsetPanel(session, "tabs", selected = "panelNetwork")
  })
  
  observeEvent(input$okServerProtection, {
    removeModal()
  })
  
  observeEvent(input$okquerySelection, {
    removeModal()
  })
  
  observeEvent(input$okdataUpload, {
    removeModal()
  })
  
  observeEvent(input$okSampleSelection, {
    removeModal()
  })
  
  observeEvent(input$okPrioritization, {
    network.data$prioritizationModalBool <- TRUE
    removeModal()
  })
  
  observeEvent(input$okcompleteDataUnmatchModal, {
    removeModal()
  })
  
  observeEvent(input$refresh,{
    ## At the end of the one session, user can click on refresh button for a new analysis without changing expression data selection
    hideTab(inputId = "tabs", target = 'panelCyto')
    hideTab(inputId = "tabs", target = 'panelNetwork')
    hideTab(inputId = "tabs", target = 'panelKNGP')
    hideTab(inputId = "sidetabs", target = 'InfoBox')
    hideTab(inputId = "tabs", target = 'panelExprs')
    showTab(inputId = "tabs", target = 'panelDGEx')
    showTab(inputId = "sidetabs", target = 'QuerySelection')
    updateTabsetPanel(session, "sidetabs", selected = "QuerySelection")
  })
  
  observeEvent(input$okKNGPModal, {
    network.data$KNGPModalBool <- TRUE
    removeModal()
  })
  
  observeEvent(input$okNullNetwork, {
    removeModal()
  })
  
  observeEvent(input$filterTypeSelection,{
    if(input$ShowSIF == 0){return()}
    
    if(input$filterTypeSelection == 'Sliders'){
      updateSliderInput(
        session = session,
        inputId = 'slidermiRNA',
        value = c(as.numeric(input$textSlidermiRNALeft), as.numeric(input$textSlidermiRNARight))
      )
      updateSliderInput(
        session = session,
        inputId = 'slidermRNA',
        value = c(as.numeric(input$textSlidermRNALeft), as.numeric(input$textSlidermRNARight))
      )
    }
    
    if(input$filterTypeSelection == 'Numbers'){
      updateNumericInput(
        session = session,
        inputId = 'textSlidermiRNALeft',
        value = as.numeric(input$slidermiRNA[1])
      )
      updateNumericInput(
        session = session,
        inputId = 'textSlidermiRNARight',
        value = as.numeric(input$slidermiRNA[2])
      )
      updateNumericInput(
        session = session,
        inputId = 'textSlidermRNALeft',
        value = as.numeric(input$slidermRNA[1])
      )
      updateNumericInput(
        session = session,
        inputId = 'textSlidermRNARight',
        value = as.numeric(input$slidermRNA[2])
      )
    }
    
  })
  
  ############## Simple Action Button Collection ENDS ##############
  ##################################################################
  
  #################### Reactivity Collection #######################
  ##################################################################
  
  observeEvent(input$miRNAExpressionFile$datapath,{
    network.data$miRNAExpressionUploaded <- input$miRNAExpressionFile
    shinyjs::enable('ShowExpressionData')
  })
  
  observeEvent(input$mRNAExpressionFile$datapath,{
    network.data$mRNAExpressionUploaded <- input$mRNAExpressionFile
    shinyjs::enable('ShowExpressionData')
  })
  output$mirnaUI <- reactive({
    network.data$mirnaFile
  })
  
  outputOptions(output, 'mirnaUI', suspendWhenHidden = FALSE)
  
  output$mrnaUI <- reactive({
    network.data$mrnaFile
  })
  
  outputOptions(output, 'mrnaUI', suspendWhenHidden = FALSE)
  
  observeEvent(input$miRNA_Controls,{
    ## Expression data upload requires assigning control and treatment samples.
    ## This assignment is done via column wise selection in data.table of expression data.
    ## Next four functions are used for control and treatment sample selection for microRNA and mRNA data...
    shiny::validate(shiny::need(network.data$mirnaFile, message = "Please upload or download a dataset first"))
    if(is.null(network.data$mirnaFile) && input$ExampleData != 'Example Data'){return()}
    ControlsmiRNA <- as.data.frame(network.data$mirnaFile[,input$ExpressionFilemiRNA_columns_selected])
    rownames(ControlsmiRNA) <- rownames(network.data$mirnaFile)
    ControlsmiRNA <- as.data.frame.matrix(ControlsmiRNA)
    network.data$ControlsmiRNAColumns <- colnames(network.data$mirnaFile)[input$ExpressionFilemiRNA_columns_selected]
    colnames(ControlsmiRNA) <- network.data$ControlsmiRNAColumns
    network.data$ControlsmiRNA <- ControlsmiRNA
  })
  
  observeEvent(input$miRNA_Treatments,{
    shiny::validate(shiny::need(network.data$mirnaFile, message = "Please upload or download a dataset first"))
    if(is.null(network.data$mirnaFile) && input$ExampleData != 'Example Data'){return()}
    TreatmentsmiRNA <- as.data.frame(network.data$mirnaFile[,input$ExpressionFilemiRNA_columns_selected])
    rownames(TreatmentsmiRNA) <- rownames(network.data$mirnaFile)
    TreatmentsmiRNA <- as.data.frame.matrix(TreatmentsmiRNA)
    network.data$TreatmentsmiRNAColumns <- colnames(network.data$mirnaFile)[input$ExpressionFilemiRNA_columns_selected]
    colnames(TreatmentsmiRNA) <- network.data$TreatmentsmiRNAColumns
    network.data$TreatmentsmiRNA <- TreatmentsmiRNA
  })
  
  observeEvent(input$Gene_Controls,{
    shiny::validate(shiny::need(network.data$mrnaFile, message = "Please upload or download a dataset first"))
    if(is.null(network.data$mrnaFile) && input$ExampleData != 'Example Data'){return()}
    ControlsmRNA <- as.data.frame(network.data$mrnaFile[,input$ExpressionFilemRNA_columns_selected])
    rownames(ControlsmRNA) <- rownames(network.data$mrnaFile)
    ControlsmRNA <- as.data.frame.matrix(ControlsmRNA)
    network.data$ControlsmRNAColumns <- colnames(network.data$mrnaFile)[input$ExpressionFilemRNA_columns_selected]
    colnames(ControlsmRNA) <- network.data$ControlsmRNAColumns
    network.data$ControlsmRNA <- ControlsmRNA
  })
  
  observeEvent(input$Gene_Treatments,{
    shiny::validate(shiny::need(network.data$mrnaFile, message = "Please upload or download a dataset first"))
    if(is.null(network.data$mrnaFile) && input$ExampleData != 'Example Data'){return()}
    TreatmentsmRNA <- as.data.frame(network.data$mrnaFile[,input$ExpressionFilemRNA_columns_selected])
    rownames(TreatmentsmRNA) <- rownames(network.data$mrnaFile)
    TreatmentsmRNA <- as.data.frame.matrix(TreatmentsmRNA)
    network.data$TreatmentsmRNAColumns <- colnames(network.data$mrnaFile)[input$ExpressionFilemRNA_columns_selected]
    colnames(TreatmentsmRNA) <- network.data$TreatmentsmRNAColumns
    network.data$TreatmentsmRNA <- TreatmentsmRNA
  })
  
  observe({
    input$ExampleData
    input$ExampleDataSelection
    hideTab(inputId = 'tabs', target = 'panelExprs')
  })
  
  observeEvent(input$tabs,{
    
    if(input$tabs == 'panelNetwork'){
      shinyjs::showElement(id= "networkLayoutPanel")
      shinyjs::showElement(id= "networkCorrelationPhysicsPanel")
    }
    if(input$tabs == 'panelKNGP'){
      shinyjs::hideElement(id= "networkLayoutPanel")
      shinyjs::hideElement(id= "networkCorrelationPhysicsPanel")
    }
    
  })
  
  observeEvent(input$filterTypeSelection,{
    
    if(input$filterTypeSelection == 'Sliders'){
      shinyjs::hideElement('textSlidermiRNALeftUI')
      shinyjs::hideElement('textSlidermiRNARightUI')
      shinyjs::hideElement('textSlidermRNALeftUI')
      shinyjs::hideElement('textSlidermRNARightUI')
      shinyjs::showElement('slidermiRNAUI')
      shinyjs::showElement('slidermRNAUI')
    }
    else if(input$filterTypeSelection == 'Numbers'){
      shinyjs::hideElement('slidermiRNAUI')
      shinyjs::hideElement('slidermRNAUI')
      shinyjs::showElement('textSlidermiRNALeftUI')
      shinyjs::showElement('textSlidermiRNARightUI')
      shinyjs::showElement('textSlidermRNALeftUI')
      shinyjs::showElement('textSlidermRNARightUI')
    }
    
  })
  
  observeEvent(input$correlationSlider,{
    input$networkLayout
    input$correlationCheckbox
    
    networkNodes <- network.data$networkMaster[[1]]
    networkNodes$color <- networkNodes$colorSave
    networkNodes$totalCountsSave <- networkNodes$totalCounts
    
    if(input$correlationCheckbox == TRUE){
      if(input$correlationSlider == '0%'){
        correlationSliderValue <- 0
        correlationSliderFormatted <- ''
      }
      else{
        correlationSliderFormatted <- gsub(pattern = '%', replacement = '', input$correlationSlider)
        correlationSliderFormatted <- unlist(strsplit(correlationSliderFormatted, ' '))
      }
      
      if(correlationSliderFormatted[1] == 'Positive'){
        correlationSliderValue <- as.numeric(correlationSliderFormatted[2])
        networkNodes$totalCounts <- networkNodes$totalCountsSave
        networkNodes$totalCounts <- 100 - networkNodes$totalCounts
      }
      else if(correlationSliderFormatted[1] == 'Negative'){
        correlationSliderValue <- as.numeric(correlationSliderFormatted[2])
        networkNodes$totalCounts <- networkNodes$totalCountsSave
      }
      
      if(correlationSliderValue > 0) {
        if(nrow(networkNodes[networkNodes$totalCounts < correlationSliderValue,]) == 0){
          networkNodes$color <- networkNodes$colorSave
        }
        else {
          networkNodes[networkNodes$totalCounts < correlationSliderValue,]$color <- "#EAE9E9"
          visNetworkProxy("visNetworkUI") %>% visUpdateNodes(networkNodes, updateOptions = FALSE) %>% visFit()
        }
      }
      else{
        networkNodes$color <- networkNodes$colorSave
        visNetworkProxy("visNetworkUI") %>% visUpdateNodes(networkNodes, updateOptions = FALSE) %>% visFit()
      }
    }
    else{
      if(input$correlationCheckbox == FALSE){
        networkNodes$color <- networkNodes$colorSave
        visNetworkProxy("visNetworkUI") %>% visUpdateNodes(networkNodes, updateOptions = FALSE) %>% visFit()
      }
    }
    
    network.data$colorLayoutSaved <- networkNodes$color
    network.data$sliderSelectionSaved <- input$correlationSlider
    network.data$iterator <- network.data$iterator + 1
  })
  
  observeEvent(input$sliderColorSpectrum,{
    
    networkNodes <- network.data$networkMaster[[1]]
    networkNodes$color <- networkNodes$colorSave
    
    lowExprs <- which(as.numeric(networkNodes$color) < input$sliderColorSpectrum[1])
    highExprs <- which(as.numeric(networkNodes$color) > input$sliderColorSpectrum[2])
    
    temp <- cut(as.numeric(networkNodes$color), breaks = seq(input$sliderColorSpectrum[1], input$sliderColorSpectrum[2], length.out = 100), include.lowest = TRUE, ordered_result = FALSE)
    colorsForGraph <- colorRampPalette(c('green', 'grey', 'red'))(99)[temp]
    ## Nodes that are not differentially expressed in control vs treatment samples are colored grey
    networkNodes$color <- colorsForGraph
    
    if(length(lowExprs) > 0){
      networkNodes[lowExprs,]$color <- '#008000'
    }
    if(length(highExprs) > 0){
      networkNodes[highExprs,]$color <- '#FF0000'
    }
    visNetworkProxy("visNetworkUI") %>% visUpdateNodes(networkNodes, updateOptions = FALSE) %>% visFit()
    
  })
  
  
  observeEvent(input$correlationCheckbox,{
    updateSliderTextInput(session, inputId = "correlationSlider", selected = '0%')
  })
  
  observeEvent(input$physicsEngineNetwork,{
    visNetworkProxy("visNetworkUI") %>%
      visFit()
  })
  
  observeEvent(input$removeDataUpload,{
    
    shinyjs::reset('miRNAExpressionFile')
    shinyjs::reset('mRNAExpressionFile')
    network.data$mirnaFile <- 'No Data'
    network.data$mrnaFile <- 'No Data'
    network.data$mirnaFile2 <- 0
    network.data$mrnaFile2 <- 0
    network.data$mRNAExpressionUploaded <- NULL
    network.data$miRNAExpressionUploaded <- NULL
    
  })
  
  plotInput <- reactive({
    
    ggscatter(
      network.data$networkPlotFile, x = "logFC", y = "Correlation",
      color = "Classification", palette = "Dark2",
      add = "reg.line"
    ) +
      facet_wrap(~Classification) +
      stat_cor(label.y = -10) +
      stat_regline_equation(label.y = -15) +
      ggtitle(paste(network.data$QueryMasterID, ' ', input$ExampleDataSelection, sep=''))
    
  })
  
  ################### Reactivity Collection END ####################
  ##################################################################
  
  #################### UI Modules Collection #######################
  ##################################################################
  
  output$ExprsInfomiRNAText <- renderUI({
    ExprsInfomiRNABox <- processExprsInfomiRNABoxTextServer(network.data$mirnaFile, network.data$TreatmentsmiRNA, network.data$ControlsmiRNA)
    lapply(ExprsInfomiRNABox, function(x)
      tags$div(style = "font-family:'Open Sans'; font-weight: 200; font-size:125%; line-height: 1.5; text-decoration: underline; margin:2px; color: #464646;", x)
    )
  })
  
  output$ExprsInfomiRNANumbers <- renderUI({
    ExprsInfomiRNABox <- processExprsInfomiRNABoxNumbersServer(network.data$mirnaFile, network.data$TreatmentsmiRNA, network.data$ControlsmiRNA)
    lapply(ExprsInfomiRNABox, function(x)
      tags$div(style = "font-family:'Open Sans'; font-weight: 200; font-size:125%; line-height: 1.5; margin:2px; color: #464646;", x)
    )
  })
  
  output$ExprsInfomRNAText <- renderUI({
    ExprsInfomRNATextBox <- processExprsInfomRNABoxTextServer(network.data$mrnaFile, network.data$TreatmentsmRNA, network.data$ControlsmRNA)
    lapply(ExprsInfomRNATextBox, function(x)
      tags$div(style = "font-family:'Open Sans'; font-weight: 200; font-size:125%; line-height: 1.5; text-decoration: underline; margin:2px; color: #464646;", x)
    )
  })
  
  output$ExprsInfomRNANumbers <- renderUI({
    ExprsInfomRNANumbersBox <- processExprsInfomRNABoxNumbersServer(network.data$mrnaFile, network.data$TreatmentsmRNA, network.data$ControlsmRNA)
    lapply(ExprsInfomRNANumbersBox, function(x)
      tags$div(style = "font-family:'Open Sans'; font-weight: 200; font-size:125%; line-height: 1.5; margin:2px; color: #464646;", x)
    )
  })
  
  output$SelectionmiRNA <- renderUI({
    ## miRNA Sample Selections are rendered here
    input$miRNA_Controls
    input$miRNA_Treatments
    input$ShowExpressionData
    if(input$ExampleData == 'No Data'){return()}
    
    if(input$ExampleData == 'Example Data'){
      mirnaSampleSelections <- processmiRNASampleSelectionServer(network.data$mirnaFile, input$ExampleData, colnames(network.data$ControlsmiRNA), colnames(network.data$TreatmentsmiRNA))
    }
    else if(input$ExampleData == 'Upload Data'){
      mirnaSampleSelections <- processmiRNASampleSelectionServer(network.data$mirnaFile, input$ExampleData, network.data$ControlsmiRNAColumns, network.data$TreatmentsmiRNAColumns)
    }
    
    lapply(mirnaSampleSelections, function(x)
      if (x[[2]] == 'headerMini')  { tags$div(style = "font-family: 'Open Sans'; font-weight: 200; font-size:130%; line-height: 1.25; text-decoration: underline; margin:2px; color: #464646;", tags$ul(tags$b(x[[1]]))) }
      else { tags$div(style = "font-family:'Open Sans'; font-weight: 200; font-size:130%; line-height: 1.25; margin:2px; color: #464646;", tags$ul(tags$li(x[[1]]))) }
    )
  })
  
  output$SelectionmRNA <- renderUI({
    ## mRNA Sample Selections are rendered here
    input$Gene_Controls
    input$Gene_Treatments
    input$ShowExpressionData
    
    if(input$ExampleData == 'No Data'){return()}
    
    if(input$ExampleData == 'Example Data'){
      mrnaSampleSelections <- processmRNASampleSelectionServer(network.data$mrnaFile, input$ExampleData, colnames(network.data$ControlsmRNA), colnames(network.data$TreatmentsmRNA))
    }
    else if(input$ExampleData == 'Upload Data'){
      mrnaSampleSelections <- processmRNASampleSelectionServer(network.data$mrnaFile, input$ExampleData, network.data$ControlsmRNAColumns, network.data$TreatmentsmRNAColumns)
    }
    
    
    lapply(mrnaSampleSelections, function(x)
      if (x[[2]] == 'headerMini')  { tags$div(style = "font-family: 'Open Sans'; font-weight: 200; font-size:130%; line-height: 1.25; text-decoration: underline; margin:2px; color: #464646;", tags$ul(tags$b(x[[1]]))) }
      else { tags$div(style = "font-family:'Open Sans'; font-weight: 200; font-size:130%; line-height: 1.25; margin:2px; color: #464646;", tags$ul(tags$li(x[[1]]))) }
    )
  })
  
  output$KNGPInfo <- renderUI({
    input$Prioritize
    KNGPInfo <- processKNGPInfo(input$nodeKnowledge, network.data$KNGP_root, input$edgeKnowledge, network.data$KNGP_bestf, network.data$KNGPDT)
    
    lapply(KNGPInfo, function(x)
      if (x[[2]] == 'headerMini')  { tags$div(style = "font-family: 'Open Sans'; font-weight: 200; font-size:130%; line-height: 1.25; text-decoration: underline; margin:2px; color: #464646;", tags$ul(tags$b(x[[1]]))) }
      else { tags$div(style = "font-family:'Open Sans'; font-weight: 200; font-size:130%; line-height: 1.25; margin:2px; color: #464646;", tags$ul(tags$li(x[[1]]))) }
    )
  })
  
  output$Information <- renderUI({
    ## Information box is rendered here ...
    lapply(network.data$infoBox, function(x)
      if(x[[2]] == 'headerMaxi')  { tags$div(style = "font-family: 'Open Sans'; font-weight: 400; font-size:130%; line-height: 1.25; text-decoration:underline; margin:2px; color: #550000;", tags$ul(tags$strong(x[[1]]))) }
      else if (x[[2]] == 'headerMini')  { tags$div(style = "font-family: 'Open Sans'; font-weight: 200; font-size:130%; line-height: 1.25; text-decoration: underline; margin:2px; color: #464646;", tags$ul(tags$b(x[[1]]))) }
      else { tags$div(style = "font-family:'Open Sans'; font-weight: 200; font-size:130%; line-height: 1.25; margin:2px; color: #464646;", tags$ul(tags$li(x[[1]]))) }
    )
  })
  
  output$InfoNumbers <- renderUI({
    ## Information box is rendered here ...
    lapply(network.data$infoBox, function(x)
      if (grepl('header', x[[2]]) == FALSE) { tags$div(style = "font-family: 'Open Sans'; font-weight: 200; font-size:130%; line-height: 1.25; margin:2px; color: #464646;", tags$ul(x[[2]])) }
      else { tags$div(style = "font-family: 'Open Sans'; font-weight: 200; font-size:130%; line-height: 1.25; margin:2px; color: #464646;", tags$ul(tags$br()))}
    )
  })
  
  output$queryBoxText <- renderUI({
    lapply(network.data$queryBox, function(x)
      if(x[[2]] == 'headerMaxi')  { tags$div(style = "font-family: 'Open Sans'; font-weight: 400; font-size:130%; line-height: 1.25; text-decoration:underline; margin:2px; color: #550000;", tags$ul(tags$strong(x[[1]]))) }
      else if (x[[2]] == 'headerMini')  { tags$div(style = "font-family: 'Open Sans'; font-weight: 200; font-size:130%; line-height: 1.25; text-decoration: underline; margin:2px; color: #464646;", tags$ul(tags$b(x[[1]]))) }
      else { tags$div(style = "font-family:'Open Sans'; font-weight: 200; font-size:130%; line-height: 1.25; margin:2px; color: #464646;;", tags$ul(tags$li(x[[1]]))) }
    )
  })
  
  output$queryBoxValues <- renderUI({
    lapply(network.data$queryBox, function(x)
      if (grepl('header', x[[2]]) == FALSE) { tags$div(style = "font-family: 'Open Sans'; font-weight: 200; font-size:130%; line-height: 1.25; margin:2px; color: #464646;", tags$ul(x[[2]])) }
      else { tags$div(style = "font-family: 'Open Sans'; font-weight: 200; font-size:130%; line-height: 1.25; margin:2px; color: #464646;", tags$ul(tags$br()))}
    )
  })  
  
  output$prospectiveText <- renderUI({
    
    input$GenerateNetwork
    ## Information box is rendered here ...
    tmpProspectives <- processProspectivesServerTextShiny()
    lapply(tmpProspectives, function(x)
      if(grepl('BREAKLINE', x) == FALSE){
        tags$div(style = "font-family:'Open Sans'; font-weight: 200; font-size:17px; color: #464646; padding-top:10px;", x)
      }
      else {
        tags$div(tags$ul(tags$br()))
      }
    )
  })
  
  output$prospectiveNumbers <- renderUI({
    
    input$GenerateNetwork
    ## Information box is rendered here ...
    tmpProspectives <- processProspectivesServerNumbersShiny()
    lapply(tmpProspectives, function(x)
      if(grepl('BREAKLINE', x) == FALSE){
        tags$div(style = "font-family:'Open Sans'; font-weight: 400; font-size:17px; color: #550000; padding-top:10px; text-decoration: underline;", x)
      }
      else {
        tags$div(tags$ul(tags$br()))
      }
    )
  })
  
  output$histmiRNA <- renderPlot({
    ## Render histograms here to help out with filtering.
    ## Expression value distributions are displayed in histograms
    input$ShowSIF
    if(is.null(input$ShowSIF)){
      return()
    }
    shiny::validate(
      shiny::need(unique(network.data$CytoscapeFile$DiffmiRNA) != 0, message = 'Data not loaded!')
    )
    miRNA.hist <- as.numeric(as.character(unique(cbind(network.data$CytoscapeFile$miRNA, network.data$CytoscapeFile$DiffmiRNA))[,2]))
    miRNA.hist <- miRNA.hist[miRNA.hist != 0]
    hist(miRNA.hist , main = 'logFC values (miRNA)', breaks = 30, col=rgb(0.3, 0.5, 1, 0.4), ylab='', xlab='')
  })
  
  output$histmRNA <- renderPlot({
    ## Histogram for mRNA data, for expression value distribution..
    input$ShowSIF
    if(is.null(input$ShowSIF)){
      return()
    }
    shiny::validate(
      shiny::need(unique(network.data$CytoscapeFile$DiffmRNA) != 0, message = 'Data not loaded!')
    )
    Gene.hist <- as.numeric(as.character(unique(cbind(network.data$CytoscapeFile$Gene, network.data$CytoscapeFile$DiffmRNA))[,2]))
    Gene.hist <- Gene.hist[Gene.hist != 0]
    hist(Gene.hist , main = 'logFC values (mRNA)', breaks = 30, col=rgb(0.3, 0.5, 1, 0.4), ylab='', xlab='')
  })
  
  output$slidermiRNAUI <- renderUI({
    ## Slider bar rendering for filtering based on expression value, miRNA data
    sliderInput("slidermiRNA", label =  HTML("<p style= 'font-family:Open Sans; color:#464646; font-size:20px; font-weight:500; padding-left:10px;'> Filter miRNA nodes by logFC values </p>"), width= "100%", min = ceiling_dec(network.data$quantileValues$miRNA00, 1), max = floor_dec(network.data$quantileValues$miRNA100, 1), value=c(ceiling_dec(network.data$quantileValues$miRNA05, 1), floor_dec(network.data$quantileValues$miRNA95, 1)), step=0.1)
  })
  
  output$slidermRNAUI <- renderUI({
    ## Slider bar rendering for filtering based on expression value, mRNA data
    sliderInput("slidermRNA", label =  HTML("<p style= 'font-family:Open Sans; color:#464646; font-size:20px; font-weight:500; padding-left:10px;'> Filter mRNA nodes by logFC values </p>"), width= "100%", min = ceiling_dec(network.data$quantileValues$mRNA00, 1), max = floor_dec(network.data$quantileValues$mRNA100, 1), value=c(ceiling_dec(network.data$quantileValues$mRNA05, 1), floor_dec(network.data$quantileValues$mRNA95, 1)), step=0.1)
  })
  
  output$textSlidermiRNALeftUI <- renderUI({
    numericInput("textSlidermiRNALeft", label = HTML("<p style='font-family:Open Sans; color:#464646; font-size:18px; font-weight:500; padding-left:20px; padding-top:15px;'>Lower logFC - miRNA</p>"), value= ceiling_dec(network.data$quantileValues$miRNA05, 1))
  })
  
  output$textSlidermiRNARightUI <- renderUI({
    numericInput("textSlidermiRNARight", label = HTML("<p style='font-family:Open Sans; color:#464646; font-size:18px; font-weight:500; padding-left:20px; padding-top:15px;'>Upper logFC - miRNA</p>"), value= floor_dec(network.data$quantileValues$miRNA95, 1))
  })
  
  output$textSlidermRNALeftUI <- renderUI({
    numericInput("textSlidermRNALeft", label = HTML("<p style='font-family:Open Sans; color:#464646; font-size:18px; font-weight:500; padding-left:20px; padding-top:15px;'>Lower logFC - mRNA</p>"), value= ceiling_dec(network.data$quantileValues$mRNA05, 1))
  })
  
  output$textSlidermRNARightUI <- renderUI({
    numericInput("textSlidermRNARight", label = HTML("<p style='font-family:Open Sans; color:#464646; font-size:18px; font-weight:500; padding-left:20px; padding-top:15px;'>Upper logFC - mRNA</p>"), value= floor_dec(network.data$quantileValues$mRNA95, 1))
  })
  
  output$sliderDegreeUI <- renderUI({
    ## Slider bar rendering for filtering based on expression value, mRNA data
    sliderInput("sliderDegree", label =  HTML("<p style= 'font-family:Open Sans; color:#464646; font-size:20px; font-weight:500; padding-left:10px;'> Filter Nodes by degree </p>"), width= "100%", min = 1, max = network.data$sliderParameterMax, value=1, step=1)
  })
  
  output$correlationSliderUI <- renderUI({
    shinyWidgets::sliderTextInput("correlationSlider", label = '',  width= "100%", force_edges = FALSE, choices = c(paste('Positive ', seq(100, 5, -5), '%', sep=''), '0%' , paste('Negative ', seq(5, 100, 5), '%', sep='')), selected = '0%')
  })
  
  output$prioritizationSliderDegreeUI <- renderUI({
    ## Slider bar rendering for prioritization based on node degree
    sliderInput("prioritizationSliderDegree", label =  HTML("<p style= 'font-family:Open Sans; color:#464646; font-size:16px; font-weight:500; padding-left:10px;'> Determine Root Nodes by degree </p>"), width= "100%", min = 1, max = network.data$sliderParameterMax, value = round(quantile(x = c(1, network.data$sliderParameterMax), 0.75)), step=1)
  })
  
  output$chooseRootNode <- renderUI({
    
    if(toupper(network.data$QueryMasterID) %in% toupper(network.data$nodeListForPrioritization)){
      selectedNode <- toupper(network.data$QueryMasterID)
    }
    else{
      selectedNode = ''
    }
    ## Network nodes are rendered here to select for prioritization  ...
    selectInput('rootNodeInput',  HTML("<p style= 'font-family:Open Sans; color:#464646; font-size:16px; font-weight:500; padding-left:10px;'> Select node(s) for prioritization! </p>"), network.data$nodeListForPrioritization, selected = selectedNode, multiple=TRUE, width = '600px')
  })
  
  output$choose_node <- renderUI({
    if(is.null(network.data$nodeList)){
      selectInput("selectName", HTML("<p style= 'font-family:Open Sans; color:#464646; font-size:20px; font-weight:500; padding-left:10px;'> Select Node by ID:  </p>"), choices = c(""), multiple=TRUE, width = '600px')
    }
    else{
      selectInput("selectName", HTML("<p style= 'font-family:Open Sans; color:#464646; font-size:20px; font-weight:500; padding-left:10px;'> Select Node by ID:  </p>"), choices = c("", network.data$nodeList), multiple=TRUE, width = '600px')
    }
  })
  
  ################## UI Modules Collection ENDS ####################
  ##################################################################
  
  ##################### renderDT Collection  #######################
  ##################################################################
  
  output$ExampleDataInfoDTUI <- DT::renderDataTable({
    network.data$exampleDataInfoBox <- processExampleDataInfoBoxServer(examplePhenoData)
    DT::datatable(network.data$exampleDataInfoBox, rownames = FALSE, colnames = '', escape = FALSE, options= list(pageLength=5, dom = 't', columnDefs = list(list(autoWidth  = T))))
  })
  
  output$ExpressionFilemiRNA <- DT::renderDataTable({
    ## Visualize the expression data of miRNAs with DT package as a table
    input$ShowExpressionData
    if (input$ShowExpressionData == 0) {return()}
    
    if(input$ExampleData == 'Example Data'){
      DT::datatable(network.data$mirnaFile, selection = list(target='column'),rownames = TRUE, options= list(pageLength=2, lengthMenu = c(2, 5, 10, 20, 100), columnDefs = list(list(autoWidth  = T))), colnames = c('ID' = 1))
    }
    else if (input$ExampleData != 'Example Data' & network.data$mirnaFile != 'No Data'){
      DT::datatable(network.data$mirnaFile, selection = list(target='column'),rownames = TRUE, options= list(pageLength=2, lengthMenu = c(2, 5, 10, 20, 100), columnDefs = list(list(autoWidth  = T))), colnames = c('ID' = 1))
    }
    else {
      DT::datatable(data.table (' Expression data was not found .. '), rownames = FALSE, colnames = c('ID' = 1))
    }
  })
  
  
  output$ExpressionFilemRNA <- DT::renderDataTable ({
    ## Visualize the expression data of mRNAs with DT package as a table
    input$ShowExpressionData
    if (input$ShowExpressionData== 0) {return()}
    if(input$ExampleData == 'Example Data'){
      DT::datatable(network.data$mrnaFile, selection = list(target='column'),rownames = TRUE, options= list(pageLength=2, lengthMenu = c(2, 5, 10, 20, 100), columnDefs = list(list(autoWidth  = T))), colnames = c('ID' = 1))
    }
    else if (input$ExampleData != 'Example Data' & network.data$mrnaFile != 'No Data'){
      DT::datatable(network.data$mrnaFile, selection = list(target='column'), rownames = TRUE, options= list(pageLength=2, lengthMenu = c(2, 5, 10, 20, 100), columnDefs = list(list(autoWidth  = T))), colnames = c('ID' = 1))
    }
    else {
      DT::datatable(data.table(' Expression data was not found .. '), rownames = FALSE, colnames = c('ID' = 1))
    }
  })
  
  
  output$CytoscapeFileDT <- DT::renderDataTable({
    ## Visualize the Cytoscape file input with DT package as a table.
    ## Depending on the absence or presence of the expression datas table is generated
    if (input$ShowSIF== 0) {return()}
    network.data$CytoscapeFile <- as.data.frame.matrix(network.data$CytoscapeFile)
    if(is.null(network.data$mrnaFile) && is.null(network.data$mirnaFile)){
      return(DT::datatable(network.data$CytoscapeFile[,-3:-4],selection = list(target='column'), rownames = TRUE, escape = FALSE, options= list(pageLength=10), colnames = c('ID' = 1)))
    }
    else if(is.null(network.data$mrnaFile)){
      network.data$CytoscapeFile$DiffmiRNA <- round(as.numeric(as.character(network.data$CytoscapeFile$DiffmiRNA)), 3)
      return(DT::datatable(network.data$CytoscapeFile[,-4],selection = list(target='column'), rownames = TRUE, escape = FALSE, options= list(pageLength=10), colnames = c('ID' = 1)))
    }
    else if(is.null(network.data$mirnaFile)){
      network.data$CytoscapeFile$DiffmRNA <- round(as.numeric(as.character(network.data$CytoscapeFile$DiffmRNA)), 3)
      return(DT::datatable(network.data$CytoscapeFile[,-3],selection = list(target='column'), rownames = TRUE, escape = FALSE, options= list(pageLength=10), colnames = c('ID' = 1)))
    }
    else{
      network.data$CytoscapeFile$DiffmiRNA <- round(as.numeric(as.character(network.data$CytoscapeFile$DiffmiRNA)), 3)
      network.data$CytoscapeFile$DiffmRNA <- round(as.numeric(as.character(network.data$CytoscapeFile$DiffmRNA)), 3)
      return(DT::datatable(network.data$CytoscapeFile,selection = list(target='column'), rownames = TRUE, escape = FALSE, options= list(pageLength=10), colnames = c('ID' = 1)))}
  })
  
  output$KNGPUI <- DT::renderDataTable ({
    input$Prioritize
    ## Visualize the expression data of mRNAs with DT package as a table
    if(nrow(network.data$KNGPDT) > 0){
      return(DT::datatable(network.data$KNGPDT, rownames = TRUE, options= list(pageLength=10), colnames = c('ID' = 1)))
    }
    else{
      return(DT::datatable('Knowledge Network Gene Prioritization algorithm did not bring any results ... '), rownames = FALSE, colnames = c('ID' = 1))
    }
  })
  
  output$phenoDTmiRNA <- DT::renderDataTable ({
    input$ShowExpressionData
    if(input$ShowExpressionData == 0){
      return()
    }
    Sys.sleep(time = 0.25)
    ## Visualize the expression data of mRNAs with DT package as a table
    if(nrow(network.data$phenoDTmiRNA) > 0){
      return(DT::datatable(network.data$phenoDTmiRNA, rownames = TRUE, escape = FALSE, options= list(processing=FALSE, pageLength=2, lengthMenu = c(2, 5, 10), columnDefs = list(list(autoWidth  = T,
                                                                                                                                                                                      targets = network.data$targetsColmiRNA,
                                                                                                                                                                                      render = JS(
                                                                                                                                                                                        "function(data, type, row, meta) {",
                                                                                                                                                                                        "return type === 'display' && data.length > 20 ?",
                                                                                                                                                                                        "'<span title=\"' + data + '\">' + data.substr(0, 20) + '...</span>' : data;",
                                                                                                                                                                                        "}")))), callback = JS('table.page(3).draw(false);'), colnames = c('ID' = 1)))
    }
    else{
      DT::datatable(data.table (' Phenodata was not found .. '), rownames = FALSE, colnames = c('ID' = 1))
    }
  })
  
  output$phenoDTmRNA <- DT::renderDataTable ({
    input$ShowExpressionData
    if(input$ShowExpressionData == 0){
      return()
    }
    Sys.sleep(time = 0.25)
    ## Visualize the expression data of mRNAs with DT package as a table
    if(nrow(network.data$phenoDTmRNA) > 0){
      return(DT::datatable(network.data$phenoDTmRNA, rownames = TRUE, escape = FALSE, options= list(processing=FALSE, pageLength=2, lengthMenu = c(2, 5, 10), columnDefs = list(list(autoWidth  = T,
                                                                                                                                                                                     targets = network.data$targetsColmRNA,
                                                                                                                                                                                     render = JS(
                                                                                                                                                                                       "function(data, type, row, meta) {",
                                                                                                                                                                                       "return type === 'display' && data.length > 20 ?",
                                                                                                                                                                                       "'<span title=\"' + data + '\">' + data.substr(0, 20) + '...</span>' : data;",
                                                                                                                                                                                       "}")))), callback = JS('table.page(3).draw(false);'), colnames = c('ID' = 1)))
    }
    else{
      DT::datatable(data.table (' Phenodata was not found .. '), rownames = FALSE, colnames = c('ID' = 1))
    }
  })
  
  output$DGExDTmiRNAUI <- DT::renderDataTable ({
    input$ShowQuery
    if(input$ShowQuery == 0){return()}
    if(input$ExampleData != 'No Data') {
      if(network.data$mirnaFile == 'No Data'){return(DT::datatable(data.table (' Expression data was not found .. '), rownames = FALSE, colnames = c('ID' = 1)))}
      if(nrow(network.data$DGExDTmiRNA) > 0){
        return(DT::datatable(network.data$DGExDTmiRNA, rownames = TRUE, options= list(pageLength=5), colnames = c('ID' = 1)))
      }
    }
    else{
      return(DT::datatable(data.table (' Expression data was not found .. '), rownames = FALSE, colnames = c('ID' = 1)))
    }
  })
  
  output$DGExDTmRNAUI <- DT::renderDataTable ({
    input$ShowQuery
    if(input$ShowQuery == 0){return()}
    if(input$ExampleData != 'No Data') {
      if(network.data$mrnaFile == 'No Data'){return(DT::datatable(data.table (' Expression data was not found .. '), rownames = FALSE, colnames = c('ID' = 1)))}
      if(nrow(network.data$DGExDTmRNA) > 0){
        return(DT::datatable(network.data$DGExDTmRNA, rownames = TRUE, options= list(pageLength=5), colnames = c('ID' = 1)))
      }
    }
    else{
      return(DT::datatable(data.table (' Expression data was not found .. '), rownames = FALSE, colnames = c('ID' = 1)))
    }
  })
  
  ################### renderDT Collection ENDS #####################
  ##################################################################
  
  ##################### Download Collection  #######################
  ##################################################################
  
  output$downloadLimmamiRNA <- downloadHandler(
    ## File name is randomly generated
    filename = function() { paste(input$ExampleDataSelection, '_limmaResults', '.txt', sep='') },
    content = function(filename) {
      write.table(network.data$DGExDTmiRNA, filename, sep = '\t')
    }
  )
  
  output$downloadLimmamRNA <- downloadHandler(
    filename = function() { paste(input$ExampleDataSelection, '_limmaResults','.txt', sep='') },
    content = function(filename) {
      write.table(network.data$DGExDTmRNA, filename, sep = '\t')
    }
  )
  
  output$downloadNetworkProperties <- downloadHandler(
    filename = function() { paste(network.data$QueryMasterID, "_", input$miRNADatabase, "_miRNA(", network.data$filterTypemiRNA[1], "-", network.data$filterTypemiRNA[2], ")", "_mRNA(", network.data$filterTypemRNA[1], "-", network.data$filterTypemRNA[2], ")_", input$ExampleDataSelection, '_networkInfo','.txt', sep='') },
    content = function(filename) {
      write.table(network.data$networkInfoFile, filename, sep = '\t', row.names = FALSE)
    }
  )
  
  output$downloadNetworkInteractions <- downloadHandler(
    filename = function() { paste(network.data$QueryMasterID, "_", input$miRNADatabase, "_miRNA(", network.data$filterTypemiRNA[1], "-", network.data$filterTypemiRNA[2], ")", "_mRNA(", network.data$filterTypemRNA[1], "-", network.data$filterTypemRNA[2], ")_", input$ExampleDataSelection, '_networkInteractions','.txt', sep='') },
    content = function(filename) {
      write.table(network.data$visNetworkFile, filename, sep = '\t', row.names = FALSE)
    }
  )
  
  output$downloadNetworkPlot <- downloadHandler(
    filename = function() { paste(network.data$QueryMasterID, "_", input$miRNADatabase, "_miRNA(", network.data$filterTypemiRNA[1], "-", network.data$filterTypemiRNA[2], ")", "_mRNA(", network.data$filterTypemRNA[1], "-", network.data$filterTypemRNA[2], ")_", input$ExampleDataSelection, '_correlationPlot','.png', sep='') },
    content = function(filename) {
      plotShiny <- plotInput()
      plotShiny <- annotate_figure(plotShiny, top = text_grob(paste("Fisher's test p value for Genes: ", round(as.numeric(network.data$GenesFisher$p.value), digits = 3), "     Fisher's test p value for miRNAs: ", round(as.numeric(network.data$miRNAsFisher$p.value), digits = 3), sep=''), color = "red"))
      ggsave(filename, plot = plotShiny, device = "png")
    }
  )
  
  output$downloadNetwork <- downloadHandler(
    
    filename = function() {
      paste(network.data$networkDownloadName, '.html', sep='')
    },
    content = function(con) {
      layoutUgly <- layoutNameConverter(input$networkLayout)
      if(input$networkArrangement == 'Initial State'){
        
        visNetwork(network.data$Nodes, network.data$Edges, main = "Odyssey Network", submain = network.data$QueryMasterID, height = '800px', width = '1600px') %>%
          visIgraphLayout(layout = layoutUgly, physics = input$physicsEngineNetwork, type ='full') %>%
          visNodes(scaling = list(label = list(enabled = T)), borderWidth = 2, borderWidthSelected = 4, color = list(border = network.data$Nodes$color)) %>%
          visGroups(groupname = "miRNA", shape = "triangle", color = 'grey') %>%
          visGroups(groupname = "Gene", color = 'grey') %>%
          visGroups(groupname = paste("Nodes", "\n", "without", "\n", "Expression", "\n", "Data", sep = "", collapse = ""), shape = 'square', color = '#BFEFFF') %>%
          visLegend(width = 0.10, position = 'right', zoom = FALSE) %>%
          visExport(label = "Export as .png", name = network.data$networkDownloadName, style=paste0("font-family: Open Sans;
                  width: 8%;
                  background: linear-gradient(#6d7070, #474949 50%, #3d3f3f);
                  border: 1px solid #2e2f2f;
                  border-radius: 4px;
                  color: White;
                  font-size: 16px;
                  background: #3d3f3f;")) %>%
          visOptions(autoResize = TRUE) %>% visSave(con, background = '#F5F5F5')
        
      }
      
      else if(input$networkArrangement == 'Prioritized'){
        
        visNetwork(network.data$Nodes, network.data$Edges, main = "Odyssey Network", submain = network.data$QueryMasterID, height = '800px', width = '1600px') %>%
          visIgraphLayout(layout = layoutUgly, physics = input$physicsEngineNetwork, type ='full') %>%
          visNodes(borderWidth = 2, borderWidthSelected = 4, color=list(border=network.data$Nodes$color)) %>%
          visGroups(groupname = "High - 10%", shape = "diamond") %>%
          visGroups(groupname = "Medium - 25%", shape = "hexagon") %>%
          visGroups(groupname = "Low", shape = 'rectangle') %>%
          visLegend(width = 0.10, position = 'right', zoom = FALSE) %>%
          visClusteringByGroup(groups = c("Low", "Medium - 25%")) %>%
          visExport(label = "Export as .png", name = network.data$networkDownloadName,
                    style=paste0("font-family: Open Sans;
                  width: 8%;
                  background: linear-gradient(#6d7070, #474949 50%, #3d3f3f);
                  border: 1px solid #2e2f2f;
                  border-radius: 4px;
                  color: White;
                  font-size: 16px;
                  background: #3d3f3f;")) %>%
          visOptions(autoResize = TRUE) %>% visSave(con, background = '#F5F5F5')
      }
    }
  )
  
  output$downloadCyto <- downloadHandler(
    ## File name is randomly generated
    filename = function() { paste(network.data$QueryMasterID, " ", input$miRNADatabase, " ", network.data$downloadHandle, '.txt', sep='') },
    content = function(filename) {
      write.table(network.data$CytoscapeFile, filename, sep = '\t')
    }
  )
  
  ################### Download Collection ENDS #####################
  ##################################################################
  
  ################### small function Collection ####################
  ##################################################################
  
  processExprsInfomiRNABoxTextServer <- function(mirnaFile, TreatmentsmiRNA, ControlsmiRNA){
    
    if(mirnaFile == 'No Data'){return()}
    infoBox1 <- list(paste("Number of samples in data:\n\n ", sep=" "))
    infoBox2 <- list(paste("Number of selected treatment samples:\n\n ", sep=" "))
    infoBox3 <- list(paste("Number of selected control samples:\n\n ", sep=" "))
    return(list(infoBox1, infoBox2, infoBox3))
    
  }
  
  processExprsInfomiRNABoxNumbersServer <- function(mirnaFile, TreatmentsmiRNA, ControlsmiRNA){
    
    if(mirnaFile == 'No Data'){return()}
    infoBox1 <- list(ncol(network.data$mirnaFile))
    infoBox2 <- list(ncol(network.data$TreatmentsmiRNA))
    infoBox3 <- list(ncol(network.data$ControlsmiRNA))
    return(list(infoBox1, infoBox2, infoBox3))
    
  }
  
  processExprsInfomRNABoxTextServer <- function(mrnaFile, TreatmentsmRNA, ControlsmRNA){
    
    if(mrnaFile == 'No Data'){return()}
    infoBox1 <- list(paste("Number of samples in data:\n\n ", sep=" "))
    infoBox2 <- list(paste("Number of selected treatment samples:\n\n ", sep=" "))
    infoBox3 <- list(paste("Number of selected control samples:\n\n ", sep=" "))
    return(list(infoBox1, infoBox2, infoBox3))
    
  }
  
  processExprsInfomRNABoxNumbersServer <- function(mrnaFile, TreatmentsmRNA, ControlsmRNA){
    
    if(mrnaFile == 'No Data'){return()}
    infoBox1 <- list(ncol(network.data$mrnaFile))
    infoBox2 <- list(ncol(network.data$TreatmentsmRNA))
    infoBox3 <- list(ncol(network.data$ControlsmRNA))
    return(list(infoBox1, infoBox2, infoBox3))
    
  }
  
  processmiRNASampleSelectionServer <- function(mirnaFile, ExampleData, ControlsmiRNA, TreatmentsmiRNA){
    
    if(mirnaFile == 'No Data'){return()}
    else if(ExampleData != 'No Data'){
      infoBox0.5 <- list(paste('These control samples are selected:\n\n', sep = " "), 'headerMini')
      infoBox1 <- list(paste(ControlsmiRNA, sep = "", collapse = ", "), '')
      infoBox1.5 <- list(paste('These treatment samples are selected:\n\n', sep = " "), 'headerMini')
      infoBox2 <- list(paste(TreatmentsmiRNA, sep = "", collapse = ", "), '')
      return(list(infoBox0.5, infoBox1, infoBox1.5, infoBox2))
    }
  }
  
  processmRNASampleSelectionServer <- function(mrnaFile, ExampleData, ControlsmRNA, TreatmentsmRNA){
    
    if(mrnaFile == 'No Data'){return()}
    else if(ExampleData != 'No Data'){
      infoBox0.5 <- list(paste('These control samples are selected:\n\n', sep = " "), 'headerMini')
      infoBox1 <- list(paste(ControlsmRNA, sep = "", collapse = ", "), '')
      infoBox1.5 <- list(paste('These treatment samples are selected:\n\n', sep = " "), 'headerMini')
      infoBox2 <- list(paste(TreatmentsmRNA, sep = "", collapse = ", "), '')
      return(list(infoBox0.5, infoBox1, infoBox1.5, infoBox2))
    }
  }
  
  processQueryBoxServerOdyssey <- function(DiffmRNAStatic, DiffmiRNAStatic) {
    ## Function for information box processing.
    ## Information box helps user recall the options chosen during analysis.
    ## Information box is displayed as a summary at the end of the analysis.
    
    if(!is.null(network.data$QueryMasterID) ){
      infoBox0.75 <- list(paste("User Query & Expression Data","\n\n", sep=" "), 'headerMaxi')
      infoBox1 <- list(paste("Your query is: ","\n\n", sep=" "), network.data$QueryMasterID)
    }
    else {
      infoBox0.75 <- NA
      infoBox1 <- NA
    }
    if(input$ShowSIF != 0){
      if(!is.null(network.data$mirnaFile)){
        infoBox2 <- list(paste("Range of values of miRNA data is:\n\n ",sep=" "), paste(round(DiffmiRNAStatic, digits=2), collapse=' , '))
      }
      else {
        infoBox2 <- NA
      }
      if(!is.null(network.data$mrnaFile)){
        infoBox3 <- list(paste("Range of values of mRNA data is:\n\n ", sep=" "), paste(round(DiffmRNAStatic, digits=2), collapse=' , '))
      }
      else {
        infoBox3 <- NA
      }
      infoBox <- list(infoBox0.75, infoBox1, infoBox2, infoBox3)
      queryBox <- infoBox[!is.na(infoBox)]
    }
    return(queryBox)
  }
  
  processExampleDataInfoBoxServer <- function(examplePhenoData) {
    
    exampleDataInfoDT <- data.table()
    if(input$ExampleData == 'Example Data'){
      
      examplePhenoSelection <- examplePhenoData[grep(input$ExampleDataSelection, rownames(examplePhenoData)),]
      Headers <- c('Selection:', 'Description:', 'Sample Comparison:', 'GEO URL', 'PubMed URL')
      GEO_URL <- paste0("<a href='",as.character(examplePhenoSelection$`GEO URL`),"'>",unlist(strsplit(input$ExampleDataSelection, '_'))[1],"</a>")
      PubMed_URL <- paste0("<a href='",as.character(examplePhenoSelection$`PubMed URL`),"'>",unlist(strsplit(input$ExampleDataSelection, '_'))[1],"</a>")
      Values <- c(input$ExampleDataSelection, as.character(examplePhenoSelection$Description), paste(as.character(examplePhenoSelection$Control), ' vs. ', as.character(examplePhenoSelection$Treatment), sep=''), GEO_URL, PubMed_URL)
      
      exampleDataInfoDT$Headers <- Headers
      exampleDataInfoDT$Values <- Values
      
    }
    
    return(exampleDataInfoDT)
  }
  
  processInfoBoxServerOdyssey <- function() {
    ## Similar to function above. Some HTML formatting is done on the information box
    
    if(input$ExampleData == 'Example Data'){
      infoBox4.5 <- list(paste("Network Options","\n\n", sep=" "), 'headerMaxi')
      infoBox5 <- list(paste("Dataset selected:", "\n", sep=" "), input$ExampleDataSelection)
    }
    else {
      infoBox4.5 <- NA
      infoBox5 <- NA
    }
    infoBox6 <- list(paste("Include Secondary Interactions: ", "\n\n",sep=" "), input$SecondaryInteractions)
    infoBox8 <- list(paste("Filtering nodes from initial query: ", "\n\n", sep=" "), input$InitQueryFilter)
    
    if(input$miRNADatabase == "TargetScan&miRNet"){
      infoBox10 <- list(paste("Database selected: ", "\n\n", sep=" "),paste(input$miRNADatabase, "\n\n", sep=" "))
      infoBox10.5 <- list(paste("Combination method: ", "\n\n", sep=" "),paste(input$unionorintersect, "\n\n", sep=" "))
      infoBox11 <- list(paste("User Notice!", sep= " "), 'headerMini')
      infoBox12 <- list(paste("Edge colors define interaction information:", sep= " "), 'headerMini')
      infoBox13 <- list(paste("TargetScan database ONLY are represented by: ", sep= " "), 'Blue')
      infoBox14 <- list(paste("miRNet database ONLY are represented by: ", sep= " "), 'Red')
      infoBox15 <- list(paste("TargetScan&miRNet BOTH are represented by: ", sep= " "), 'Green')
    }
    else{
      infoBox10 <- list(paste("Database selected: ", "\n\n", sep=""),input$miRNADatabase)
      infoBox10.5 <- NA
      infoBox11 <- NA
      infoBox12 <- NA
      infoBox13 <- NA
      infoBox14 <- NA
      infoBox15 <- NA
    }
    
    infoBox <- list(infoBox4.5,infoBox5,infoBox6,infoBox8,infoBox10,infoBox10.5,infoBox11,infoBox12,infoBox13,infoBox14,infoBox15)  ## infoBox9 removed because layouts became deprecated
    infoBox <- infoBox[!is.na(infoBox)]
    return(infoBox)
  }
  
  processProspectivesServerTextShiny <- function() {
    
    prospectiveBox1 <- list(paste('Number of nodes: ', sep=' '))
    prospectiveBox1.5 <- list(paste('\n\n', 'BREAKLINE', sep = ' '))
    prospectiveBox2 <- list(paste('Number of edges: ', sep = ' '))
    return(list(prospectiveBox1, prospectiveBox1.5, prospectiveBox2)) 
  }
  
  processProspectivesServerNumbersShiny <- function() {
    
    if(input$filterTypeSelection == 'Sliders'){
      filterTypemiRNA <- input$slidermiRNA
      filterTypemRNA <- input$slidermRNA
    }
    else if(input$filterTypeSelection == 'Numbers'){
      filterTypemiRNA <- c(input$textSlidermiRNALeft, input$textSlidermiRNARight)
      filterTypemRNA <- c(input$textSlidermRNALeft, input$textSlidermRNARight)
    }
    
    prospectiveNetwork <- ProcessFilteringNetworkDGex(network.data$CytoscapeFile, filterTypemiRNA, filterTypemRNA, input$InitQueryFilter, network.data$PrimerQuery)
    prospectiveNetwork <- ProcessFilteringNetworkDegree(prospectiveNetwork[[1]], input$DegreeFilter, input$sliderDegree)
    network.data$sliderParameterMax <- prospectiveNetwork[[2]]
    prospectiveNetwork <- prospectiveNetwork[[1]]
    
    prospectiveNodeNumber <- length(unique(c(prospectiveNetwork[,1],prospectiveNetwork[,2])))
    prospectiveEdgeNumber <- nrow(prospectiveNetwork)
    
    prospectiveBox1 <- list(prospectiveNodeNumber)
    prospectiveBox1.5 <- list('BREAKLINE')
    prospectiveBox2 <- list(prospectiveEdgeNumber)
    return(list(prospectiveBox1, prospectiveBox1.5, prospectiveBox2)) 
  }
  
  processKNGPInfo <- function(nodeKnowledge, rootNodes, edgeKnowledge, bestf, KNGPDT){
    
    infoBox1 <- list(paste('Root node selection:\n\n', sep = " "), 'headerMini')
    infoBox2 <- list(paste(nodeKnowledge, sep = "", collapse = ", "), '')
    infoBox3 <- list(paste('These are the root nodes:\n\n', sep = " ", collapse = ",  "), 'headerMini')
    infoBox4 <- list(paste(rootNodes, sep = "", collapse = ", "), '')
    infoBox5 <- list(paste('Edge Knowledge Parameter:\n\n', sep = " "), 'headerMini')
    infoBox6 <- list(paste(edgeKnowledge, sep = "", collapse = ", "), '')
    infoBox7 <- list(paste('Results were obtained for bestf value of:\n\n', sep = " "), 'headerMini')
    infoBox8 <- list(paste(bestf, sep = "", collapse = ", "), '')
    
    if(nrow(KNGPDT) < 1){
      infoBox1 <- list(paste('Prioritization did not bring any results! Consider choosing a different setting.. \n\n'), 'headerMini')
      return(list(infoBox1))
    }
    else{
      return(list(infoBox1, infoBox2, infoBox3, infoBox4, infoBox5, infoBox6, infoBox7, infoBox8))
    }
  }
  
  serverProtectionModal <- function(failed = FALSE) {
    modalDialog(
      if(failed)
        div(tags$b('Please try another query...', style = "color: red;")),
      
      footer = tagList(
        actionButton("okServerProtection", "OK")
      )
    )
  }
  
  nullNetworkModal <- function(failed = FALSE) {
    modalDialog(
      if(failed)
        div(tags$b('No nodes remained in the network... Try a less stringent filtering!', style = "color: red;")),
      
      footer = tagList(
        actionButton("okNullNetwork", "OK")
      )
    )
  }
  
  querySelectionModal <- function(failed = FALSE) {
    modalDialog(
      if(failed)
        div(tags$b(network.data$queryModal, style = "color: red;")),
      
      footer = tagList(
        actionButton("okquerySelection", "OK")
      )
    )
  }
  
  dataUploadModal <- function(failed = FALSE) {
    modalDialog(
      if(failed)
        div(tags$b('Your data is being processed... This process usually takes < 1 min, depending on the size of the data!', style = "color: red;")),
      
      footer = tagList(
        actionButton("okdataUpload", "OK")
      )
    )
  }
  
  sampleSelectionModal <- function(failed = FALSE) {
    modalDialog(
      if(failed)
        div(tags$b("Data upload requires Control and Treatment sample selection!", style = "color: red;")),
      
      footer = tagList(
        actionButton("okSampleSelection", "OK")
      )
    )
  }
  
  completeDataUnmatchModal <- function(failed = FALSE) {
    modalDialog(
      if(failed)
        div(tags$b(network.data$completeDataUnmatchMessage, style = "color: red;")),
      
      footer = tagList(
        actionButton("okcompleteDataUnmatchModal", "OK")
      )
    )
  }
  
  prioritizationModal <- function(failed = FALSE) {
    modalDialog(
      if(failed)
        div(tags$b("Prioritization takes < 1 mins depending on the size of the network...", style = "font-family:courier; color: black; font-size: 20px;")),
      
      footer = tagList(
        actionButton("okPrioritization", "OK")
      )
    )
  }
  
  KNGPModal <- function(failed = FALSE) {
    modalDialog(
      if(failed)
        div(tags$b("Check Prioritization tab for further information!", style = "font-family:courier; color: black; font-size: 20px;")),
      
      footer = tagList(
        actionButton("okKNGPModal", "OK")
      )
    )
  }
  
  ################### small function Collection ENDS ###############
  ##################################################################
  
  ##################### Modules Collection #########################
  ##################################################################
  
  ###################### ShowExpressionData  #######################
  ##################################################################
  
  observeEvent(input$ShowExpressionData,{
    
    if(input$ExampleData == 'Example Data'){
      showTab(inputId = "subtabs", target = 'subTabPheno')
      if(input$ExampleDataSelection == "GSE81867"){
        network.data$phenoDTmiRNA <- data.table()
        network.data$phenoDTmRNA <- data.table()
      }
      else{
        exampleDataSelectionSplit <- unlist(strsplit(input$ExampleDataSelection, '_'))[1]
        network.data$phenoDTmRNA <- network.data$mRNAphenoDB[[exampleDataSelectionSplit]]
        network.data$targetsColmRNA <- 1:ncol(network.data$phenoDTmRNA)
        network.data$phenoDTmiRNA <- network.data$miRNAphenoDB[[exampleDataSelectionSplit]]
        network.data$targetsColmiRNA <- 1:ncol(network.data$phenoDTmiRNA)
      }
    }
    else{
      hideTab(inputId = "subtabs", target = 'subTabPheno')
    }
    
    network.data$TreatmentsmiRNAColumns <- 'None' 
    network.data$ControlsmiRNAColumns <- 'None' 
    network.data$TreatmentsmRNAColumns <- 'None' 
    network.data$ControlsmRNAColumns <- 'None'
    
    if(input$ExampleData == 'Upload Data'){
      if(!is.null(network.data$miRNAExpressionUploaded) | !is.null(network.data$mRNAExpressionUploaded)){
        showModal(dataUploadModal(failed = TRUE))
      }
    }
    
    ExpressionDataBundle <- ProcessExpressionData(input$miRNADatabase, network.data$mRNAExpressionUploaded, network.data$miRNAExpressionUploaded, input$ExampleData, input$ExampleDataSelection, GenePlatformAnnotations, miRNA.identifiers, miRNAexprsDB, mRNAexprsDB, limmaListmiRNA, limmaListmRNA)
    
    network.data$mirnaFile <- ExpressionDataBundle[[1]]
    network.data$mrnaFile <- ExpressionDataBundle[[2]]
    network.data$datacheck.mirna <- ExpressionDataBundle[[3]]
    network.data$datacheck.mrna <- ExpressionDataBundle[[4]]
    network.data$miRNAphenoDB <- ExpressionDataBundle[[5]]
    network.data$mRNAphenoDB <- ExpressionDataBundle[[6]]
    network.data$TreatmentsmiRNA <- ExpressionDataBundle[[7]]
    network.data$ControlsmiRNA <- ExpressionDataBundle[[8]]
    network.data$TreatmentsmRNA <- ExpressionDataBundle[[9]]
    network.data$ControlsmRNA <- ExpressionDataBundle[[10]]
    
    network.data$ExpressionDataBundlePreviousSession <- list(network.data$mRNAExpressionUploaded, network.data$miRNAExpressionUploaded, input$ExampleData, input$ExampleDataSelection)
    
    hideTab(inputId = "sidetabs", target = 'DataSelection')
    showTab(inputId = "sidetabs", target = 'SampleSelection')
    showTab(inputId = "tabs", target = 'panelExprs')
    updateTabsetPanel(session, "sidetabs", selected = "SampleSelection")
    updateTabsetPanel(session, "tabs", selected = "panelExprs")
    shinyjs::hideElement('miRNAExprsDataInfo')
    shinyjs::hideElement('mRNAExprsDataInfo')
    
    if(network.data$mirnaFile == 'No Data' && network.data$mrnaFile == 'No Data'){
      click("ShowQuery")
    }
    if(network.data$mirnaFile != 'No Data'){
      shinyjs::showElement('miRNAExprsDataInfo')
    }
    if(network.data$mrnaFile != 'No Data'){
      shinyjs::showElement('mRNAExprsDataInfo')
    }
    
  })
  
  #################### ShowExpressionData ENDS #####################
  ##################################################################
  
  ########################## ShowQuery #############################
  ##################################################################
  
  observeEvent(input$ShowQuery,{
    
    DGExBundle <- ProcessDGExMTIdb(network.data$mirnaFile, network.data$mrnaFile, network.data$datacheck.mirna, network.data$datacheck.mrna, input$miRNADatabase, input$unionorintersect, input$ExampleData, network.data$TreatmentsmiRNA, network.data$ControlsmiRNA, network.data$TreatmentsmRNA, network.data$ControlsmRNA)
    
    if(DGExBundle[[1]] == 'NoSampleSelection'){
      showModal(sampleSelectionModal(failed = TRUE))
      return()
    }
    
    network.data$mirnaFile <- DGExBundle[[1]]
    network.data$mrnaFile <- DGExBundle[[2]]
    network.data$mirnaFile2 <- DGExBundle[[3]]
    network.data$mrnaFile2 <- DGExBundle[[4]]
    network.data$DGExDTmiRNA <- DGExBundle[[5]]
    network.data$DGExDTmRNA <- DGExBundle[[6]]
    
    hideTab(inputId = "sidetabs", target = 'SampleSelection')
    hideTab(inputId = "tabs", target = 'panelExprs')
    showTab(inputId = "sidetabs", target = 'QuerySelection')
    showTab(inputId = "tabs", target = 'panelDGEx')
    updateTabsetPanel(session, "sidetabs", selected = "QuerySelection")
    updateTabsetPanel(session, "tabs", selected = "panelDGEx")
    shinyjs::showElement('downloadLimmamiRNA')
    shinyjs::showElement('downloadLimmamRNA')
    
    if(input$ExampleData == "Example Data"){
      network.data$DGExDTmiRNA[,c(1:6)] <- round(network.data$DGExDTmiRNA[,c(1:6)], digits=3)
      network.data$DGExDTmRNA[,c(1:6)] <- round(network.data$DGExDTmRNA[,c(1:6)], digits=3)
    }
    else if(input$ExampleData == "Upload Data"){
      if(network.data$mirnaFile != 'No Data'){
        network.data$DGExDTmiRNA[,c(1:6)] <- round(network.data$DGExDTmiRNA[,c(1:6)], digits=3)
      }
      else{
        shinyjs::hideElement('downloadLimmamiRNA')
      }
      if(network.data$mrnaFile != 'No Data'){
        network.data$DGExDTmRNA[,c(1:6)] <- round(network.data$DGExDTmRNA[,c(1:6)], digits=3)
      }
      else{
        shinyjs::hideElement('downloadLimmamRNA')
      }
    }
    else if(input$ExampleData == 'No Data'){
      shinyjs::hideElement('downloadLimmamiRNA')
      shinyjs::hideElement('downloadLimmamRNA')
    }
    
  })
  
  ######################## ShowQuery ENDS ##########################
  ##################################################################
  
  ######################## ShowSIF Module ##########################
  ##################################################################
  
  observeEvent(input$ShowSIF,{
    ## Generate cytoscape file to be used for network creation.
    
    if(nchar(input$id) > 20){
      showModal(serverProtectionModal(failed = TRUE))
      updateTextInput(session, "id", value = '')
      return()
    }
    
    if(input$id == ""){
      network.data$queryModal <- HTML('Please enter a valid Gene or microRNA ID!')
      showModal(querySelectionModal(failed = TRUE))
      return()
    }
    
    network.data$QueryMasterID <- input$id
    SIFBundle <- ProcessSIF(network.data$QueryMasterID, input$querySelection, network.data$mirnaFile, network.data$mrnaFile, network.data$mirnaFile2, network.data$mrnaFile2, input$SecondaryInteractions, input$miRNADatabase, input$unionorintersect)
    queryMasterChecker <- SIFBundle[[1]]
    network.data$multipleMatch <- SIFBundle[[2]]
    network.data$queryModal <- SIFBundle[[3]]
    
    if(network.data$multipleMatch == TRUE) {
      closeMatches <- SIFBundle[[1]]
      network.data$queryModal <- HTML(paste(network.data$queryModal,"<br>", paste(closeMatches, collapse="<br>")))
      showModal(querySelectionModal(failed = TRUE))
      return()
    }
    
    else if(network.data$multipleMatch == 'NoMatch'){
      network.data$queryModal <- SIFBundle[[3]]
      showModal(querySelectionModal(failed = TRUE))
      return()
    }
    
    network.data$PrimerQuery <- SIFBundle[[4]]
    network.data$CytoscapeFile <- SIFBundle[[5]]
    network.data$DiffmRNAStatic <- SIFBundle[[6]]
    network.data$DiffmiRNAStatic <- SIFBundle[[7]]
    network.data$mirnaFile2 <- SIFBundle[[8]]
    network.data$mrnaFile2 <- SIFBundle[[9]]
    network.data$quantileValues <- SIFBundle[[10]]
    completeDataUnmatch <- SIFBundle[[11]]
    network.data$nullmiRNAs <- SIFBundle[[12]]
    network.data$nullGenes <- SIFBundle[[13]]
    
    pubmed_join <- dplyr::inner_join(network.data$CytoscapeFile[,c(1,2)], mirnet_interactiondb[,c(1,2)])
    mirna_pubmed <- pubmed_join$miRNA
    gene_pubmed <- pubmed_join$Gene
    mirnet_mirna_pubmed <- mirnet_interactiondb[mirnet_interactiondb$miRNA %in% mirna_pubmed,]
    pubmed_join <- mirnet_mirna_pubmed[mirnet_mirna_pubmed$Gene %in% gene_pubmed,]
    network.data$CytoscapeFile <- merge(network.data$CytoscapeFile, pubmed_join, by = c('miRNA', 'Gene'))
    
    network.data$CytoscapeFile$PubmedID <- lapply(lapply(strsplit(as.character(network.data$CytoscapeFile$PubmedID), '|', fixed = T), function(x){if(!is.na(x)){paste0("<a href=https://www.ncbi.nlm.nih.gov/pubmed/", x, " target='_blank'>", x, "</a>")}}), function(y){paste(y, collapse = '<br>')})
    
    if(network.data$multipleMatch == 'OneExpressionMatch'){
      showModal(querySelectionModal(failed = TRUE))
      updateTextInput(session, "id", value = queryMasterChecker)
      network.data$multipleMatch <- FALSE
    }
    else if(network.data$multipleMatch == 'NoExpressionMatch'){
      showModal(querySelectionModal(failed = TRUE))
    }
    
    if(completeDataUnmatch != TRUE){
      if(completeDataUnmatch == 'miRNA'){
        network.data$completeDataUnmatchMessage <- 'None of the miRNA ID(s) is present in expression data identifiers! Proceeding without full expression data benefits...'
        showModal(completeDataUnmatchModal(failed = TRUE)) 
      }
      else if(completeDataUnmatch == 'Gene'){
        network.data$completeDataUnmatchMessage <- 'None of the Gene ID(s) is present in expression data identifiers! Proceeding without full expression data benefits...'
        showModal(completeDataUnmatchModal(failed = TRUE)) 
      }
    }
    
    network.data$queryBox <- processQueryBoxServerOdyssey(network.data$DiffmRNAStatic, network.data$DiffmiRNAStatic)
    network.data$infoBox <- processInfoBoxServerOdyssey()
    
    hideTab(inputId = "sidetabs", target = 'QuerySelection')
    hideTab(inputId = "tabs", target = 'panelDGEx')
    showTab(inputId = "sidetabs", target = 'NetworkOptions')
    showTab(inputId = "tabs", target = 'panelCyto')
    
    updateTabsetPanel(session, "tabs", selected = "panelCyto")
    ## Change the tab on the top to Cytoscape file visualization tab
    updateTabsetPanel(session, "sidetabs", selected = "NetworkOptions")
    ## Change the tab on the top to Cytoscape file visualization tab
    
  }, priority = 5)
  
  ###################### ShowSIF Module ENDS #######################
  ##################################################################
  
  #################### Generate Network Module #####################
  ##################################################################
  
  observeEvent(input$GenerateNetwork,{
    
    hideTab(inputId = "tabs", target = 'panelKNGP')
    updateSelectInput(session, "networkArrangement", selected = 'Initial State')
    updateSliderTextInput(session, inputId = "correlationSlider", selected = '0%')
    updateCheckboxInput(session, 'correlationCheckbox', value = FALSE)
    shinyjs::hide('networkArrangement')
    ## Network generation is executed here ...
    
    if(is.null(dim(network.data$mirnaFile2)) & is.null(dim(network.data$mrnaFile2))){
      if(network.data$mirnaFile2 == 0 & network.data$mrnaFile2 == 0){
        shinyjs::hideElement(id = "panelPrioritize")
      }
    }
    
    else if(input$ExampleData == 'No Data'){
      shinyjs::hideElement(id = "panelPrioritize")
    }
    else{
      shinyjs::showElement(id = "panelPrioritize")
    }
    hideTab(inputId = "sidetabs", target = 'NetworkOptions')
    hideTab(inputId = "tabs", target = 'panelCyto')
    showTab(inputId = "sidetabs", target = 'InfoBox')
    showTab(inputId = "tabs", target = 'panelNetwork')
    if(input$physicsEngineNetwork == FALSE){
      shinyjs::showElement(id= "networkLayoutPanel")
    }
    updateTabsetPanel(session, "tabs", selected = "panelNetwork")
    ## Change the tab on the top to Network visualization tab..
    updateTabsetPanel(session, "sidetabs", selected = "InfoBox")
    ## Change the tab on the top to Network visualization tab..
    
    network.data$queryBox <- processQueryBoxServerOdyssey(network.data$DiffmRNAStatic, network.data$DiffmiRNAStatic)
    network.data$infoBox <- processInfoBoxServerOdyssey()
    
    if(input$filterTypeSelection == 'Sliders'){
      network.data$filterTypemiRNA <- input$slidermiRNA
      network.data$filterTypemRNA <- input$slidermRNA
    }
    else if(input$filterTypeSelection == 'Numbers'){
      network.data$filterTypemiRNA <- c(input$textSlidermiRNALeft, input$textSlidermiRNARight)
      network.data$filterTypemRNA <- c(input$textSlidermRNALeft, input$textSlidermRNARight)
    }
    
    NetworkFile <- ProcessFilteringNetworkDGex(network.data$CytoscapeFile, network.data$filterTypemiRNA, network.data$filterTypemRNA, input$InitQueryFilter, network.data$PrimerQuery)
    NetworkFile <- ProcessFilteringNetworkDegree(NetworkFile[[1]], input$DegreeFilter, input$sliderDegree)
    NetworkFile <- NetworkFile[[1]]
    
    shiny::validate(
      need(dim(NetworkFile)[1] > 0, "No nodes remained in the network. Please consider using a less strict logFC filtering")
    )
    
    nodeNames <- unique(c(NetworkFile$miRNA, NetworkFile$Gene))
    GOTable.tmp <- unique(unlist(network.data$genetoGO[unique(nodeNames)]))
    
    tmpmiRNA <- data.frame(unique(cbind(NetworkFile$miRNA, NetworkFile$DiffmiRNA)))
    colnames(tmpmiRNA) <- c('Node', 'Expression')
    tmpmRNA <- data.frame(unique(cbind(NetworkFile$Gene, NetworkFile$DiffmRNA)))
    colnames(tmpmRNA) <- c('Node', 'Expression')
    tempNetworkObject <- rbind(tmpmiRNA, tmpmRNA)
    network.data$nodeListForPrioritization <- as.character(tempNetworkObject$Node)
    
    network.data$GO_Table <- GOTable.tmp
    rm(GOTable.tmp)
    
    network.data$NetworkFile <- NetworkFile
    network.data$visNetworkFile <- NetworkFile
    
  })
  
  ################# Generate Network Module ENDS ###################
  ##################################################################
  
  ################## Visualize Network Module  #####################
  ##################################################################
  
  output$visNetworkUI <- renderVisNetwork({
    
    input$GenerateNetwork
    input$Prioritize
    
    if(input$miRNADatabase == 'TargetScan&miRNet'){
      dbNameForDL <- paste(input$miRNADatabase, '_', input$unionorintersect, sep='')
    }
    else{
      dbNameForDL <- input$miRNADatabase
    }
    if(input$ExampleData == 'Example Data'){
      dataSelectionForDL <- input$ExampleDataSelection
    }
    else{
      dataSelectionForDL <- ''
    }
    
    if(input$networkArrangement == "Initial State"){
      network.data$networkDownloadName <- paste(network.data$QueryMasterID, "_", dbNameForDL, "_miRNA(", network.data$filterTypemiRNA[1], "-", network.data$filterTypemiRNA[2], ")", "_mRNA(", network.data$filterTypemRNA[1], "-", network.data$filterTypemRNA[2], ")_", dataSelectionForDL, '_Original', sep='')
    }
    
    
    if(input$networkArrangement == "Prioritized"){
      if(input$edgeKnowledge == 'Negative Correlation'){
        edgeKnowForDL <- 'NegCorr'
      }
      else if(input$edgeKnowledge == 'Positive Correlation'){
        edgeKnowForDL <- 'PosCorr'
      }
      if(input$nodeKnowledge == 'Nodes by ID'){
        nodeKnowForDL <- paste(input$rootNodeInput, sep = '', collapse = '_') 
      }
      else if(input$nodeKnowledge == 'Nodes by degree'){
        nodeKnowForDL <- paste('prioritizedbyNodeDegree_', input$prioritizationSliderDegree, sep = '')
      }
      nodeKnowForDL <- stringr::str_trunc(nodeKnowForDL, 30)
      network.data$networkDownloadName <- paste(network.data$QueryMasterID, "_", dbNameForDL, "_miRNA(", network.data$filterTypemiRNA[1], "-", network.data$filterTypemiRNA[2], ")", "_mRNA(", network.data$filterTypemRNA[1], "-", network.data$filterTypemRNA[2], ")_", nodeKnowForDL, '_',  edgeKnowForDL, '_', dataSelectionForDL, '_Prioritized', sep='')
    }
    
    network.data$networkMaster <- ProcessvisNetworkShiny(network.data$visNetworkFile, input$miRNADatabase, input$unionorintersect, input$sliderColorSpectrum, network.data$nullmiRNAs, network.data$nullGenes)
     
    networkNodes <- network.data$networkMaster[[1]]
    networkEdges <- network.data$networkMaster[[2]]
    if(is.null(dim(networkNodes)) & is.null(dim(networkEdges))){
      showModal(nullNetworkModal(failed = TRUE))
      hideTab(inputId = "tabs", target = 'panelNetwork')
      showTab(inputId = "tabs", target = 'panelCyto')
      updateTabsetPanel(session, "sidetabs", selected = "NetworkOptions")
      updateTabsetPanel(session, "tabs", selected = "panelCyto")
      return()
    }
    else{
      showTab(inputId = "tabs", target = 'panelNetwork')
    }
    
    if(input$ExampleData == 'Example Data'){
      tmpGene <- unique(as.data.frame.matrix(cbind(network.data$NetworkFile$Gene, network.data$NetworkFile$geneCount)))
      tmpmiRNA <- unique(as.data.frame.matrix(cbind(network.data$NetworkFile$miRNA, network.data$NetworkFile$miRNACount)))
      colnames(tmpGene) <- c('ID', 'Degree')
      colnames(tmpmiRNA) <- c('ID', 'Degree')
      tmpNetwork <- rbind(tmpGene, tmpmiRNA)
      tmpNetwork <- unique(tmpNetwork)
      tmpNetwork2 <- cbind(as.character(networkNodes[[1]]), networkNodes[[7]])
      colnames(tmpNetwork2) <- c('ID', 'Correlation')
      tmpNetwork <- merge(tmpNetwork, tmpNetwork2, by='ID')
      rownames(tmpNetwork) <- tmpNetwork$ID
      tmpNetwork <- tmpNetwork[,-1]
      tmpLimma <- rbind(network.data$DGExDTmiRNA, network.data$DGExDTmRNA)
      tmpNetwork$Correlation <- gsub("[\\(\\)]", "", regmatches(tmpNetwork$Correlation, gregexpr("\\(.*?\\)", tmpNetwork$Correlation)))
      network.data$networkInfoFile <- merge(tmpNetwork, tmpLimma, by=0)
      networkPlotFile <- network.data$networkInfoFile
      
      networkPlotFile$Correlation <- gsub(pattern = '%', replacement = '', networkPlotFile$Correlation)
      networkPlotFile <- networkPlotFile[,c(-8,-9,-10)]
      networkPlotFile$Classification <- NA
      networkPlotFile[grep('HSA-', networkPlotFile$Row.names),]$Classification <- 'miRNA'
      networkPlotFile[!grepl('HSA-', networkPlotFile$Row.names),]$Classification <- 'Gene'
      networkPlotFile$Correlation <- as.numeric(networkPlotFile$Correlation)
      networkPlotFile$logFC <- as.numeric(networkPlotFile$logFC)
      network.data$networkPlotFile <- networkPlotFile
      mirnaTable <- networkPlotFile[networkPlotFile$Classification == 'miRNA',]
      mirnaTable <- as.data.frame.matrix(cbind(mirnaTable$Correlation, mirnaTable$logFC))
      geneTable <- networkPlotFile[networkPlotFile$Classification == 'Gene',]
      geneTable <- as.data.frame.matrix(cbind(geneTable$Correlation, geneTable$logFC))
      colnames(mirnaTable) <- c('Correlation', 'logFC')
      colnames(geneTable) <- c('Correlation', 'logFC')
      
      indNeg <- which(mirnaTable$Correlation > 50)
      indPos <- which(mirnaTable$Correlation <= 50)
      indDown <- which(mirnaTable$logFC <= 0)
      indUp <- which(mirnaTable$logFC > 0)
      try(mirnaTable$Correlation[indNeg] <- 'NegativeCorr')
      try(mirnaTable$Correlation[indPos] <- 'PositiveCorr')
      try(mirnaTable$logFC[indUp] <- 'UpRegulated')
      try(mirnaTable$logFC[indDown] <- 'DownRegulated')
      
      indNeg <- which(geneTable$Correlation > 50)
      indPos <- which(geneTable$Correlation <= 50)
      indDown <- which(geneTable$logFC <= 0)
      indUp <- which(geneTable$logFC > 0)
      try(geneTable$Correlation[indNeg] <- 'NegativeCorr')
      try(geneTable$Correlation[indPos] <- 'PositiveCorr')
      try(geneTable$logFC[indUp] <- 'UpRegulated')
      try(geneTable$logFC[indDown] <- 'DownRegulated')
      
      tmpMatrix <-  matrix(c('NegativeCorr', 'NegativeCorr', 'PositiveCorr', 'PositiveCorr', 'DownRegulated', 'UpRegulated', 'DownRegulated', 'UpRegulated'), nrow=4, ncol=2)
      colnames(tmpMatrix) <- c('Correlation', 'logFC')
      geneTable <- rbind(geneTable, tmpMatrix)
      mirnaTable <- rbind(mirnaTable, tmpMatrix)
      mirnaTable <- table(mirnaTable)
      geneTable <- table(geneTable)
      Genes <- geneTable - 1
      miRNAs <- mirnaTable - 1
      network.data$GenesFisher <- fisher.test(Genes)
      network.data$miRNAsFisher <- fisher.test(miRNAs)
    }
    
    if(!is.null(input$correlationSlider) & input$correlationCheckbox != FALSE){
      if(input$correlationSlider == network.data$sliderSelectionSaved){
        if(length(networkNodes$color) == length(network.data$colorLayoutSaved)){
          networkNodes$color <- network.data$colorLayoutSaved
        }
      }
    }
    
    if(network.data$iterator == 1){
      networkNodes$color <- networkNodes$colorSave
    }
    
    if(input$networkArrangement == 'Initial State'){
      ## Network is being generated here ...
      shinyjs::showElement(id = 'downloadNetworkProperties')
      shinyjs::showElement(id = 'downloadNetworkInteractions')
      shinyjs::showElement(id = 'downloadNetworkPlot')
      layoutUgly <- layoutNameConverter(input$networkLayout)
      network.data$Nodes <- networkNodes
      network.data$Edges <- networkEdges
      
      network.data$visNetworkShiny <- visNetwork(networkNodes, networkEdges, main = "Odyssey Network", submain = network.data$QueryMasterID, height = '800px', width = '1600px') %>%
        visIgraphLayout(layout = layoutUgly, physics = input$physicsEngineNetwork, type ='full') %>%
        visNodes(scaling = list(label = list(enabled = T)), borderWidth = 2, borderWidthSelected = 4, color = list(border = networkNodes$color)) %>%
        visGroups(groupname = "miRNA", shape = "triangle", color = 'grey') %>%
        visGroups(groupname = "Gene", color = 'grey') %>%
        visGroups(groupname = paste("Nodes", "\n", "without", "\n", "Expression", "\n", "Data", sep = "", collapse = ""), shape = 'square', color = '#BFEFFF') %>%
        visLegend(width = 0.10, position = 'right', zoom = FALSE) %>%
        visInteraction(hover = TRUE, multiselect = TRUE, navigationButtons = TRUE, hideEdgesOnDrag = TRUE, keyboard = TRUE, tooltipDelay = 300, tooltipStyle = 'position: fixed; visibility:hidden; padding:5px; white-space:nowrap;
font-family: Open Sans; font-size:18px; font-color:#0F0F0F;') %>%
        visExport(label = "Export as .png", name = network.data$networkDownloadName, style=paste0("font-family: Open Sans;
                  width: 8%;
                  background: linear-gradient(#6d7070, #474949 50%, #3d3f3f);
                  border: 1px solid #2e2f2f;
                  border-radius: 4px;
                  color: White;
                  font-size: 16px;
                  background: #3d3f3f;")) %>%
        visOptions(highlightNearest = TRUE, manipulation = TRUE,
                   nodesIdSelection = list(enabled = TRUE), autoResize = TRUE)
      network.data$visNetworkShiny
    }
    else if(input$networkArrangement == 'Prioritized'){
      
      updateCheckboxInput(session, "correlationCheckbox", value = FALSE)
      shinyjs::hideElement(id = 'downloadNetworkProperties')
      shinyjs::hideElement(id = 'downloadNetworkInteractions')
      shinyjs::hideElement(id = 'downloadNetworkPlot')
      networkNodes <- network.data$networkMaster[[1]]
      networkNodes <- networkNodes[,-4]
      
      networkNodes <- applyPriorityOnNetwork(networkNodes, network.data$KNGPDT)
      networkEdges <- network.data$networkMaster[[2]]
      layoutUgly <- layoutNameConverter(input$networkLayout)
      network.data$Nodes <- networkNodes
      network.data$Edges <- networkEdges
      
      network.data$visNetworkShiny <- visNetwork(networkNodes, networkEdges, main = "Odyssey Network", submain = network.data$QueryMasterID, height = '800px', width = '1600px') %>%
        visIgraphLayout(layout = layoutUgly, physics = input$physicsEngineNetwork, type ='full') %>%
        visNodes(borderWidth = 2, borderWidthSelected = 4, color=list(border=networkNodes$color)) %>%
        visGroups(groupname = "High - 10%", shape = "diamond") %>%
        visGroups(groupname = "Medium - 25%", shape = "hexagon") %>%
        visGroups(groupname = "Low", shape = 'rectangle') %>%
        visLegend(width = 0.10, position = 'right', zoom = FALSE) %>%
        visClusteringByGroup(groups = c("Low", "Medium - 25%")) %>%
        visInteraction(hover = TRUE, multiselect = TRUE, navigationButtons = TRUE, hideEdgesOnDrag = TRUE, keyboard = TRUE, tooltipDelay = 300, tooltipStyle = 'position: fixed; visibility:hidden; padding:5px; white-space:nowrap;
font-family: Open Sans; font-size:18px; font-color:#0F0F0F;') %>%
        visExport(label = "Export as .png", name = network.data$networkDownloadName, style=paste0("font-family: Open Sans;
                  width: 8%;
                  background: linear-gradient(#6d7070, #474949 50%, #3d3f3f);
                  border: 1px solid #2e2f2f;
                  border-radius: 4px;
                  color: White;
                  font-size: 16px;
                  background: #3d3f3f;")) %>%
        visOptions(highlightNearest = list(enabled = TRUE), manipulation = TRUE, nodesIdSelection = list(enabled = TRUE), collapse = list(enabled = TRUE, fit = TRUE), autoResize = TRUE)
      network.data$visNetworkShiny
    }
  })
  
  ################# Visualize Network Module ENDS ##################
  ##################################################################
  
  ################## Prioritize Network Module  ####################
  ##################################################################
  
  observeEvent({input$Prioritize
    input$nodeKnowledge
    input$edgeKnowledge
    input$prioritizationSliderDegree
    1
  },{
    ## GO term selection based graph expansion is executed here ...
    
    if(input$Prioritize == 0){
      return()
    }
    
    edgeList <- list()
    NetworkFile <- network.data$NetworkFile
    tmpmiRNA <- data.frame(unique(cbind(NetworkFile$miRNA, NetworkFile$DiffmiRNA)))
    colnames(tmpmiRNA) <- c('Node', 'Expression')
    tmpmRNA <- data.frame(unique(cbind(NetworkFile$Gene, NetworkFile$DiffmRNA)))
    colnames(tmpmRNA) <- c('Node', 'Expression')
    tempNetworkObject <- rbind(tmpmiRNA, tmpmRNA)
    
    NodeNamesVectorServer <- as.character(tempNetworkObject$Node)
    network.data$nodeListForPrioritization <- NodeNamesVectorServer
    NodeKnowledgeVectorServer <- abs(as.numeric(as.character(tempNetworkObject$Expression)))
    rm(tempNetworkObject)
    
    NetworkFile$miRNACount
    
    for(node in NodeNamesVectorServer){
      
      if(length(intersect(names(edgeList), node)) == 0){
        
        part1 <- unique(NetworkFile[NetworkFile$Gene == node,]$miRNA)
        part2 <- unique(NetworkFile[NetworkFile$miRNA == node,]$Gene)
        edgeList[[node]] <- c(part1, part2)
        
      }
      
    }
    
    network.data$KNGPDT <- data.table()
    if(network.data$prioritizationModalBool == FALSE){
      network.data$prioritizationModalBool <- TRUE
      showModal(prioritizationModal(failed = TRUE))
    }
    
    if(input$nodeKnowledge == "Nodes by ID"){
      rootNodesServer <- unique(input$rootNodeInput)
    }
    else if(input$nodeKnowledge == 'Nodes by degree'){
      edgeNumbers <- lapply(edgeList, function(x){length(x)})
      rootNodesServer <- unique(names(edgeList[which(as.numeric(edgeNumbers) > as.numeric(input$prioritizationSliderDegree))]))
    }
    
    if(input$edgeKnowledge == 'Negative Correlation'){
      
      negCorr.df <- network.data$networkMaster[[3]]
      totalCorr.df <- network.data$networkMaster[[5]]
      
      LinkKnowledgeMatrixServer <- matrix(nrow = length(NodeNamesVectorServer), ncol= length(NodeNamesVectorServer))
      for(i in 1:length(NodeNamesVectorServer)){
        for(j in 1:length(NodeNamesVectorServer)){
          if(length(grep(NodeNamesVectorServer[j],as.character(unlist(edgeList[NodeNamesVectorServer[i]])))) > 0){
            
            totalDegree <- length(unlist(edgeList[NodeNamesVectorServer[i]])) + length(unlist(edgeList[NodeNamesVectorServer[j]]))
            totalDegree2 <- as.numeric(totalCorr.df[totalCorr.df$id %in% NodeNamesVectorServer[i],]$totalCounts) + as.numeric(totalCorr.df[totalCorr.df$id %in% NodeNamesVectorServer[j],]$totalCounts)
            negDegree <- as.numeric(negCorr.df[negCorr.df$id %in% NodeNamesVectorServer[i],]$negCounts) + as.numeric(negCorr.df[negCorr.df$id %in% NodeNamesVectorServer[j],]$negCounts)
            if(length(negDegree) == 0){
              negDegree <- 0
            }
            valTemp <- negDegree / totalDegree
            LinkKnowledgeMatrixServer[i,j] <- valTemp
          }
          else{
            LinkKnowledgeMatrixServer[i,j] <- 0
          }
        }
      }
    }
    
    else if(input$edgeKnowledge == 'Positive Correlation'){
      
      posCorr.df <- network.data$networkMaster[[5]]
      totalCorr.df <- network.data$networkMaster[[5]]
      
      LinkKnowledgeMatrixServer <- matrix(nrow = length(NodeNamesVectorServer), ncol= length(NodeNamesVectorServer))
      for(i in 1:length(NodeNamesVectorServer)){
        for(j in 1:length(NodeNamesVectorServer)){
          if(length(grep(NodeNamesVectorServer[j],as.character(unlist(edgeList[NodeNamesVectorServer[i]])))) > 0){
            
            totalDegree <- length(unlist(edgeList[NodeNamesVectorServer[i]])) + length(unlist(edgeList[NodeNamesVectorServer[j]]))
            totalDegree2 <- as.numeric(totalCorr.df[totalCorr.df$id %in% NodeNamesVectorServer[i],]$totalCounts) + as.numeric(totalCorr.df[totalCorr.df$id %in% NodeNamesVectorServer[j],]$totalCounts)
            posDegree <- as.numeric(posCorr.df[posCorr.df$id %in% NodeNamesVectorServer[i],]$posCounts) + as.numeric(posCorr.df[posCorr.df$id %in% NodeNamesVectorServer[j],]$posCounts)
            if(length(posDegree) == 0){
              posDegree <- 0
            }
            valTemp <- posDegree / totalDegree
            LinkKnowledgeMatrixServer[i,j] <- valTemp
          }
          else{
            LinkKnowledgeMatrixServer[i,j] <- 0
          }
        }
      }
      
    }
    
    KNGP_results <- KNGPAlgorithm(LinkKnowledgeMatrixServer, NodeKnowledgeVectorServer, NodeNamesVectorServer, rootNodesServer)
    
    network.data$KNGP_bestf <- KNGP_results$best_f
    network.data$KNGP_root <- rootNodesServer
    network.data$KNGPDT <- KNGP_results$KNGPResult[order(KNGP_results$KNGPResult$postProb, decreasing = TRUE),]
    network.data$KNGPDT <- cbind(network.data$KNGPDT, 1:nrow(network.data$KNGPDT))
    network.data$KNGPDT <- data.table(network.data$KNGPDT)
    colnames(network.data$KNGPDT) <- c('Node_ID', 'Involvement_Score', 'Node Rank #')
    
    
    network.data$KNGPDT$Involvement_Score <- round(network.data$KNGPDT$Involvement_Score, 3)
    showTab(inputId = "tabs", target = 'panelKNGP')
    updateSelectInput(session, "networkArrangement", selected = 'Prioritized')
    shinyjs::show('networkArrangement')
    
    if(network.data$KNGPModalBool == FALSE){
      network.data$KNGPModalBool <- TRUE
      showModal(KNGPModal(failed = TRUE))
    }
    
  })
  
  ############### Prioritize Network Module  ENDS ##################
  ##################################################################
  
})
