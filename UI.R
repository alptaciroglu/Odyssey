library(shinyBS)
library(shiny)
library(DT)
library(shinyjs)
library(visNetwork)
library(data.table)
library(plyr)
library(dplyr)
library(tidyr)
library(htmlwidgets)
library(devtools)
library(R.utils)
library(igraph)
library(shinythemes)
library(shinyWidgets)
library(stringr)
library(limma)
library(ggplot2)
library(ggpubr)
library(shinycssloaders)

source('functionsArsenal.R', local=TRUE)
source('global.R', local=TRUE)


shinyUI(fluidPage(theme = shinytheme("spacelab"),
                  
                  chooseSliderSkin(skin = "Flat", color = 'DimGray'),
                  
                  tags$script("function disableButton() { document.getElementById('ShowExpressionData').disabled = true; }"),
                  
                  tags$head(tags$style(HTML('.shiny-output-error-validation { font-family:courier; color:#B54B4B;'))),
                  ## Error messages are modified to have different color and font
                  
                  tags$style(type='text/css', ".selectize-input { padding-left: 20px; border: 1px solid #080808;} .selectize-dropdown { padding-left: 10px; }.selectize-input, .selectize-input input { color: #333333; font-family: Open Sans; font-size: 16px; line-height: 20px; -webkit-font-smoothing: inherit;}.selectize-control.single .selectize-input.input-active { background: #FFAAAA;}"),
                  ## Some user input based modules are decorated with CSS elements
                  
                  tags$head(
                    ## Open Sans font type is imported
                    tags$style(HTML("@import url('//fonts.googleapis.com/css?family=Open+Sans:400,700');
      h1 {
        font-family: 'Open Sans', cursive;
        font-weight: 500;
        line-height: 2.0;
        color: #464646;
        font-size:400%;
      }

    "))
                  ),
                  
                  sidebarPanel(width=3,
                               
                               tags$style(".topimg {
                            margin-bottom:25px;
                            margin-right:15px;
                            float:left;
                          }"),
                               
                               tagList(tags$div(tags$style(type ="text/css", "a{color:#550000; font-family:Open Sans; font-size:110%; font-weight:600; text-decoration: underline;}","body {background-color: #fff; }")), tags$div(class="topimg", HTML('<a href="https://diagnose.shinyapps.io/odyssey/">', '<img src="Logo.jpg", height="100px", margin-bottom:10px"/></a>')),  tags$div(style="font-family: 'Open Sans'; font-weight: 1000; font-size:300%; margin:10px; color: #464646;", " Odyssey "), tagList(tags$div(style = "float:right; font-color:aqua; font-size:20px; font-weight:500; text-decoration:underline;", icon("download"), a("Tutorial",target="_blank",href="Odyssey_Tutorial.pdf")))),
                               
                               ## Header to hold logo, title
                               
                               ## Sidebar panel have 4 tabs
                               tabsetPanel(id="sidetabs",
                                           
                                           tabPanel("Data Selection", value="DataSelection", icon=icon('book'),
                                                    ## Data selection or upload executed here
                                                    wellPanel(
                                                      ## Example data selection or upload
                                                      fluidRow(
                                                        column(12,
                                                               selectInput('ExampleData', HTML("<p style= 'font-family:Open Sans; color:#464646; font-size:20px; font-weight:500; padding-left:10px;'> Expression Data Choice </p>"), choices=c("Example Data", "Upload Data", "No Data"), selected="Example Data"),
                                                               tags$head(tags$style(HTML("#ExampleData {padding-left:20px;}"))),
                                                               
                                                               conditionalPanel(condition="input.ExampleData == 'Example Data'",
                                                                                selectInput("ExampleDataSelection", HTML("<p style = 'font-family:Open Sans; color:#464646; font-size:20px; font-weight:500; padding-left:20px; padding-top:15px;'>  Which example data you want to use ? </p>"),
                                                                                            choices=c("GSE35389_1", "GSE35389_2", "GSE88721", "GSE39061_1", "GSE39061_2", "GSE39061_3", "GSE49697_1", "GSE49697_2", "GSE25402", "GSE32539_1", "GSE32539_2", "GSE32539_3", "GSE34681_1", "GSE34681_2", "GSE38617", "GSE40321", "GSE59702_1", "GSE59702_2", "GSE81867",
                                                                                                      "GSE104268_1", "GSE104268_2", "GSE90604"), selected="GSE38617")
                                                               ),
                                                               conditionalPanel(condition= "input.ExampleData == 'Upload Data'",
                                                                                ## Data upload part is enabled only when no example data is being used.
                                                                                fluidRow(
                                                                                  column(12,
                                                                                         HTML('<p style="font-family:open sans; color:#443E7E"> Uploading an expression file is optional </p>'),
                                                                                         fileInput('miRNAExpressionFile', HTML("<p style = 'font-family:Open Sans; color:#464646; font-size:20px; font-weight:500; padding-left:20px; padding-top:15px;'> Upload your miRNA expression file </p>"),
                                                                                                   ## First file input for miRNA expression.
                                                                                                   accept=c('text/csv',
                                                                                                            'text/comma-separated-values,text/plain')) %>%
                                                                                           {temp = .
                                                                                           temp$children[[2]]$children[[1]]$children[[1]]$children[[2]]$attribs$onchange <- "disableButton()"
                                                                                           temp},
                                                                                         fileInput('mRNAExpressionFile',  HTML("<p style = 'font-family:Open Sans; color:#464646; font-size:20px; font-weight:500; padding-left:20px; padding-top:15px;'> Upload your mRNA expression file </p>"),
                                                                                                   ## Second file input for mRNA expression.
                                                                                                   accept=c('text/csv',
                                                                                                            'text/comma-separated-values,text/plain')) %>%
                                                                                           {temp = .
                                                                                           temp$children[[2]]$children[[1]]$children[[1]]$children[[2]]$attribs$onchange <- "disableButton()"
                                                                                           temp}
                                                                                  ),
                                                                                  column(3,
                                                                                         actionButton('removeDataUpload','Reset!'),
                                                                                         tags$head(tags$style(HTML('#removeDataUpload{width: 100%; float:right; margin-top: 25px; }'))),tags$label(tags$style(HTML('#removeDataUpload{color:White; font-family:Open Sans; font-size: 18px;}'))),
                                                                                         bsTooltip("removeDataUpload", HTML('<p style="font-family:courier; background-color: black; font-size: 15px; color: white;"> Removes uploaded data from server! </p>'), placement = "top", trigger = "hover", options = list(container = "body"))
                                                                                  )
                                                                                )
                                                               ),
                                                               conditionalPanel(condition ="input.ExampleData == 'Example Data'",
                                                                                wellPanel(
                                                                                  fluidRow(
                                                                                    column(12,
                                                                                           div(style = 'overflow-x: scroll; font-family:Montserrat; font-size:120%;', DT::dataTableOutput("ExampleDataInfoDTUI"))
                                                                                    )
                                                                                  )
                                                                                )
                                                               )
                                                        )
                                                      )
                                                    ),
                                                    actionButton("ShowExpressionData","Proceed"),
                                                    ## Enter initial parameters and files then press the action button to visualize input files
                                                    
                                                    tags$head(tags$style(HTML('#ShowExpressionData{width: 50%; }'))),tags$label(tags$style(HTML('#ShowExpressionData{color:White; font-family:Open Sans; font-size: 18px;}')))
                                           ),
                                           tabPanel("Data Information", value = "SampleSelection", icon = icon('asterisk'),
                                                    conditionalPanel(condition = "input.ExampleData != 'No Data'",
                                                                     wellPanel(id = 'miRNAExprsDataInfo', HTML('<p style="font-family:Open Sans; font-size:20px; text-decoration: underline; margin:5px; color:#550000; font-weight:bold"> miRNA Data Information </p>'),
                                                                               fluidRow(
                                                                                 column(8,
                                                                                        htmlOutput("ExprsInfomiRNAText")
                                                                                        ## Show selected miRNA control and treatment samples here.
                                                                                 ),
                                                                                 column(4,
                                                                                        htmlOutput("ExprsInfomiRNANumbers")
                                                                                 )
                                                                               )
                                                                     ),
                                                                     wellPanel(id = 'mRNAExprsDataInfo', HTML('<p style="font-family:Open Sans; font-size:20px; text-decoration: underline; margin:5px; color: #550000; font-weight:bold"> mRNA Data Information </p>'),
                                                                               fluidRow(
                                                                                 column(8,
                                                                                        htmlOutput("ExprsInfomRNAText")
                                                                                        ## Show selected mRNA control and treatment samples here.
                                                                                 ),
                                                                                 column(4,
                                                                                        htmlOutput("ExprsInfomRNANumbers")
                                                                                 )
                                                                               )
                                                                     )
                                                    ),
                                                    conditionalPanel(condition = "output.mirnaUI != 'No Data'",
                                                                     wellPanel(HTML('<p style="font-family:Open Sans; font-size:20px; text-decoration: underline; margin:5px; color: #550000; font-weight:bold"> miRNA Sample Selections </p>'),
                                                                               fluidRow(
                                                                                 column(12,
                                                                                        htmlOutput("SelectionmiRNA")
                                                                                        ## Show selected miRNA control and treatment samples here.
                                                                                 )
                                                                               )
                                                                     )
                                                    ),
                                                    conditionalPanel(condition = "output.mrnaUI != 'No Data'",
                                                                     wellPanel(HTML('<p style="font-family:Open Sans; font-size:20px; text-decoration: underline; margin:5px; color: #550000; font-weight:bold"> mRNA Sample Selections </p>'),
                                                                               fluidRow(
                                                                                 column(12,
                                                                                        htmlOutput("SelectionmRNA")
                                                                                        ## Show selected mRNA control and treatment samples here.
                                                                                 )
                                                                               )
                                                                     )
                                                    )
                                           ),
                                           tabPanel("Query Selection", value = "QuerySelection", icon = icon('question-circle'),
                                                    wellPanel(
                                                      fluidRow(
                                                        column(6,
                                                               textInput("id", label = HTML("<p style='font-family:Open Sans; color:#464646; font-size:14px; font-weight:500; padding-left:20px; padding-top:15px;'> Gene or Micro-rna ID: </p>"), value=""),
                                                               tags$style("#id {font-family:Open Sans; color:#464646; font-size:24px; font-weight:500; height:76%;}"),
                                                               bsTooltip("id", HTML('<p style="font-family:courier; background-color: black; font-size: 15px; color: white;"> e.g. hsa-miR-31-5P, PMAIP1... Odyssey uses Official Gene Names and miRBase identifiers v22 </p>'), placement = "top", trigger = "hover", options = list(container = "body"))
                                                        ),
                                                        column(6,
                                                               shinyWidgets::radioGroupButtons(inputId = "querySelection", justified = TRUE, individual = TRUE, size = 'normal', label = HTML("<p style='font-family:Open Sans; color:#464646; font-size:15px; font-weight:500; padding-top:15px; padding-left:15px;'> Querying for:  </p>"), choiceNames = list(HTML("<p style= 'font-family:Open Sans; color:white; font-size:12px; font-weight:500; margin-left:10px;'> miRNA </p>"), HTML("<p style= 'font-family:Open Sans; color:white; font-size:12px; font-weight:500; margin-left:10px;'> Gene </p>")), choiceValues = c("miRNA", "Gene"), status = "primary", checkIcon = list(yes = tags$i(class = "fa fa-check-square", style = "color: white"), no = tags$i(class = "fa fa-square-o", style = "color: white")))
                                                        )
                                                      )
                                                    ),
                                                    ## Code works both for microRNA name and for gene id.
                                                    wellPanel(
                                                      fluidRow(
                                                        column(12,
                                                               selectInput("miRNADatabase", HTML("<p style= 'font-family:Open Sans; color:#464646; font-size:20px; font-weight:500; padding-left:20px; padding-top:15px;'> Interaction database </p>"),
                                                                           choices= c("TargetScan", "miRNet","TargetScan&miRNet"), selected = "miRNet"),
                                                               bsTooltip("miRNADatabase", HTML('<p style="font-family:courier; background-color: black; font-size: 15px; color: white;"> TargetScan for predicted, miRNet for experimentally validated interactions </p>'), placement = "top", trigger = "hover", options = list(container = "body")),
                                                               conditionalPanel(condition = "input.miRNADatabase == 'TargetScan&miRNet'",
                                                                                selectInput("unionorintersect", HTML("<p style = 'font-family:Open Sans; color:#464646; font-size:20px; font-weight:500; padding-left:20px; padding-top:15px;'> Select 'intersect' for common ids present in both databases, or select 'union' for ids that is present in at least one of the databases </p>"), choices= c("Intersect", "Union"), selected="Intersect")
                                                               )
                                                               ),
                                                        column(12,
                                                               checkboxInput('SecondaryInteractions', HTML("<p style= 'font-family:Open Sans; color:#464646; font-size:20px; font-weight:500; padding-left:10px;'> Include Secondary Interactions </p>"), TRUE),
                                                               tags$head(tags$style(HTML("#SecondaryInteractions {padding-left:20px;}"))),
                                                               ## Check to include targets of target in the network
                                                               bsTooltip("SecondaryInteractions", HTML('<p style="font-family:courier; background-color: black; font-size: 15px; color: white;"> Check the box to extend the network to second degree interactions </p>'), placement = "top", trigger = "hover", options = list(container = "body"))
                                                        ),
                                                        column(12,
                                                               HTML(paste0("<p style= 'font-family:Open Sans; color:Black; font-size:14px; font-weight:500; text-decoration: underline; float:right'> Last Updated <br/> 30th Jan 2022 </p>"))
                                                               )
                                                      )
                                                    )
                                           ),
                                           tabPanel("Network Options",value="NetworkOptions", icon=icon('cogs'),
                                                    wellPanel(
                                                      fluidRow(
                                                        column(6,
                                                               checkboxInput('InitQueryFilter', HTML("<p style= 'font-family:Open Sans; color:#464646; font-size:17px; font-weight:500; padding-left:10px;'> Filter initial query nodes </p>"), TRUE),
                                                               bsTooltip("InitQueryFilter", HTML('<p style="font-family:courier; background-color: black; font-size: 15px; color: white;"> Your initial query interactions will always show independant of the filters you apply unless chosen otherwise </p>'), placement = "top", trigger = "hover", options = list(container = "body")),
                                                               checkboxInput('DegreeFilter', HTML("<p style= 'font-family:Open Sans; color:#464646; font-size:17px; font-weight:500; padding-left:10px;'> Filter network by degree </p>"), FALSE),
                                                               bsTooltip("DegreeFilter", HTML('<p style="font-family:courier; background-color: black; font-size: 15px; color: white;"> Degree filter applied after Expression based filter  </p>'), placement = "top", trigger = "hover", options = list(container = "body")),
                                                               conditionalPanel(condition="input.DegreeFilter == true",
                                                                                uiOutput('sliderDegreeUI')
                                                                                
                                                               )
                                                        ),
                                                        column(4,
                                                               htmlOutput('prospectiveText')
                                                        ),
                                                        column(2,
                                                               htmlOutput('prospectiveNumbers')
                                                               )
                                                      )
                                                    ),
                                                    conditionalPanel(condition = 'input.ShowSIF',
                                                                     ## User interface is gradually made visible as the analysis proceeds.
                                                                     wellPanel(
                                                                       conditionalPanel(condition= "output.mirnaUI != 'No Data' || output.mrnaUI != 'No Data'",
                                                                                        fluidRow(
                                                                                          column(12,
                                                                                                 shinyWidgets::radioGroupButtons(inputId = "filterTypeSelection", justified = TRUE, individual = TRUE, size = 'normal', label = HTML("<p style='font-family:Open Sans; color:#464646; font-size:20px; font-weight:500; padding-top:15px; padding-left:15px;'> Filter using:  </p>"), choiceNames = list(HTML("<p style= 'font-family:Open Sans; color:white; font-size:18px; font-weight:500; margin-left:10px;'> Sliders </p>"), HTML("<p style= 'font-family:Open Sans; color:white; font-size:18px; font-weight:500; margin-left:10px;'> Numbers </p>")), choiceValues = c("Sliders", "Numbers"), status = "primary", checkIcon = list(yes = tags$i(class = "fa fa-check-square", style = "color: white"), no = tags$i(class = "fa fa-square-o", style = "color: white")))
                                                                                          )
                                                                                        )
                                                                       ),
                                                                       fluidRow(
                                                                         column(12,
                                                                                conditionalPanel(condition= "output.mirnaUI != 'No Data'",
                                                                                                 column(6,
                                                                                                        uiOutput('slidermiRNAUI')
                                                                                                 )
                                                                                                 ## miRNAs with expression value that is inside this range of slider will be filtered out
                                                                                ),
                                                                                conditionalPanel(condition= "output.mrnaUI != 'No Data'",
                                                                                                 column(6,
                                                                                                        uiOutput('slidermRNAUI')
                                                                                                 )
                                                                                                 ## mRNAs with expression value that is inside this range of slider will be filtered out
                                                                                )
                                                                         )
                                                                       ),
                                                                       fluidRow(
                                                                         column(12,
                                                                                conditionalPanel(condition= "output.mirnaUI != 'No Data'",
                                                                                                 column(3,
                                                                                                        uiOutput('textSlidermiRNALeftUI')
                                                                                                 ),
                                                                                                 column(3,
                                                                                                        uiOutput('textSlidermiRNARightUI')
                                                                                                 )
                                                                                ),
                                                                                conditionalPanel(condition= "output.mrnaUI != 'No Data'",
                                                                                                 column(3,
                                                                                                        uiOutput('textSlidermRNALeftUI')
                                                                                                 ),
                                                                                                 column(3,
                                                                                                        uiOutput('textSlidermRNARightUI')
                                                                                                 )
                                                                                )
                                                                         )
                                                                       )
                                                                     ),
                                                                     wellPanel(
                                                                       fluidRow(
                                                                         column(12,
                                                                                conditionalPanel(condition= "output.mirnaUI != 'No Data'",
                                                                                                 column(6,
                                                                                                        plotOutput("histmiRNA") ## Histograms are displayed
                                                                                                 )
                                                                                ),
                                                                                conditionalPanel(condition= "output.mrnaUI != 'No Data'",
                                                                                                 column(6,
                                                                                                        plotOutput("histmRNA")  ## Histograms are displayed
                                                                                                 )
                                                                                )
                                                                         )
                                                                       ),
                                                                       HTML('<p style="font-family:Open Sans; color:#443E7E"> Filtering intervals are automatically adjusted according to 10% & 90% quantile values </p>')
                                                                     )
                                                    )
                                           ),
                                           tabPanel('Network Statistics', value = 'InfoBox', icon = icon('info'),
                                                    wellPanel(
                                                      fluidRow(
                                                        column(8, htmlOutput('queryBoxText')),
                                                        column(4, htmlOutput('queryBoxValues'))
                                                      )
                                                    ),
                                                    tags$div(id='divInfo',
                                                             wellPanel(
                                                               fluidRow(
                                                                 column(8,htmlOutput("Information")),
                                                                 ## Show selected miRNA control and treatment samples here along with more user dependant selections
                                                                 column(4,htmlOutput("InfoNumbers"))
                                                                 ## Show selected miRNA control and treatment samples here along with more user dependant selections
                                                               )
                                                             )
                                                    ),
                                                    conditionalPanel(condition= "input.ShowSIF",
                                                                     ## GO term selection is executed here. HTML and CSS formatting of the user interface is applied.
                                                                     wellPanel(id = 'panelPrioritize', HTML("<p style='font-family: Open Sans; font-weight: 400; font-size:130%; line-height: 1.25; text-decoration:underline; margin-left:20px; margin-top:2px; color: #550000;'> <b>Prioritization Settings</b> </p>"),
                                                                               fluidRow(
                                                                                 column(6,
                                                                                        selectInput("nodeKnowledge", HTML("<p style= 'font-family:Open Sans; color:#464646; font-size:16px; font-weight:500; padding-left:20px; padding-top:15px; '> Root Node Parameters </p>"),
                                                                                                    choices= c("Nodes by ID", "Nodes by degree"), selected = "Nodes by ID"),
                                                                                        bsTooltip("nodeKnowledge", HTML('<p style="font-family:courier; background-color: black; font-size: 15px; color: white;"> Select root node(s) manually or by degree for prioritization. </p>'), placement = "top", trigger = "hover", options = list(container = "body"))
                                                                                 ),
                                                                                 column(6,
                                                                                        selectInput("edgeKnowledge", HTML("<p style= 'font-family:Open Sans; color:#464646; font-size:16px; font-weight:500; padding-left:20px; padding-top:15px; '> Edge Knowledge Parameters </p>"),
                                                                                                    choices= c("Negative Correlation", "Positive Correlation"), selected = "Negative Correlation"),
                                                                                        bsTooltip("edgeKnowledge", HTML('<p style="font-family:courier; background-color: black; font-size: 15px; color: white;"> Select edge knowledge parameter. Negative Correlation and Positive Correlation prioritizes by the negatively correlated and positively correlated node pairs respectively. </p>'), placement = "top", trigger = "hover", options = list(container = "body"))
                                                                                 ),
                                                                                 conditionalPanel(condition = "input.nodeKnowledge == 'Nodes by ID'",
                                                                                                  column(6,
                                                                                                         uiOutput("chooseRootNode"))
                                                                                 ),
                                                                                 conditionalPanel(condition="input.nodeKnowledge == 'Nodes by degree'",
                                                                                                  column(6,
                                                                                                         uiOutput('prioritizationSliderDegreeUI')
                                                                                                  )
                                                                                 ),
                                                                                 column(6,
                                                                                        actionButton('Prioritize','Prioritize!'),
                                                                                        tags$head(tags$style(HTML('#Prioritize{width: 100%; float:left; margin-top: 25px; }'))),tags$label(tags$style(HTML('#Prioritize{color:White; font-family:Open Sans; font-size: 24px;}'))),
                                                                                        bsTooltip("Prioritize", HTML('<p style="font-family:courier; background-color: black; font-size: 15px; color: white;"> Prioritize network with selected settings! </p>'), placement = "top", trigger = "hover", options = list(container = "body")))
                                                                               )
                                                                     )
                                                    )
                                           )
                               )
                  ),
                  
                  mainPanel(width=9,
                            shinyjs::useShinyjs(),
                            tabsetPanel(id='tabs',
                                        tabPanel("Data Visualization",value="panelExprs", icon=icon('columns'),
                                                 bsTooltip("panelExprs", HTML('<p style="font-family:courier; background-color: black; font-size: 15px; color: white;"> Expression data is displayed here </p>'), placement = "top", trigger = "hover", options = list(container = "body")),
                                                 ## Tab for expression files
                                                 
                                                 conditionalPanel(condition=" $('html').hasClass('shiny-busy')",
                                                                  tags$div(style="float:right; padding-right:30px; padding-top:10px;",
                                                                           tags$img(src="achilles_loading.gif",height=75,width=75)),
                                                                  tags$div(style="float:right; padding-right:30px; padding-top:10px; font-family: 'Open Sans'; font-weight: 200; font-size:125%; margin:10px; color: #464646;", "Processing ...")
                                                 ),
                                                 tabsetPanel(id='subtabs',
                                                             tabPanel('Expression Data', value='subTabExprs', icon = icon('table'),
                                                                      conditionalPanel(condition= "output.mirnaUI != 'No Data'",
                                                                                       wellPanel(HTML('<p style="font-family:Open Sans; font-size:18px; text-decoration: underline; margin:5px; color:#550000; font-weight:bold"> Expression-miRNA </p>'),
                                                                                                 ## Visualize miRNA expression file here
                                                                                                 fluidRow(
                                                                                                   column(12,
                                                                                                          div(style = 'overflow-x: scroll; font-family:Montserrat; font-size:120%;', DT::dataTableOutput("ExpressionFilemiRNA"))
                                                                                                   )),
                                                                                                 conditionalPanel(condition = "input.ExampleData != 'Example Data' & output.mirnaUI != 'No Data'",
                                                                                                                  fluidRow(
                                                                                                                    column(6,
                                                                                                                           actionButton("miRNA_Controls","miRNA Controls"),
                                                                                                                           ## Action button to select control samples of the miRNA expression file
                                                                                                                           tags$head(tags$style(HTML('#miRNA_Controls{width: 45%; }'))),tags$label(tags$style(HTML('#miRNA_Controls {color:White; font-family:Open Sans; font-size: 18px;}')))
                                                                                                                    ),
                                                                                                                    column(6,
                                                                                                                           actionButton("miRNA_Treatments","miRNA Treatments"),
                                                                                                                           ##Action button to select treatment samples of the miRNA expression file
                                                                                                                           tags$head(tags$style(HTML('#miRNA_Treatments{width: 45%; }'))),tags$label(tags$style(HTML('#miRNA_Treatments {color:White; font-family:Open Sans; font-size: 18px;}')))
                                                                                                                    )
                                                                                                                  ))
                                                                                       )),
                                                                      conditionalPanel(condition= "output.mrnaUI != 'No Data'",
                                                                                       wellPanel(HTML('<p style="font-family:Open Sans; font-size:18px; text-decoration: underline; margin:5px; color:#550000; font-weight:bold"> Expression-mRNA </p>'),
                                                                                                 fluidRow(
                                                                                                   column(12,
                                                                                                          div(style = 'overflow-x: scroll; font-family:Montserrat; font-size:120%;', DT::dataTableOutput("ExpressionFilemRNA")))),
                                                                                                 ## Visualize mRNA expression file here
                                                                                                 conditionalPanel(condition = "input.ExampleData != 'Example Data' & output.mrnaUI != 'No Data'",
                                                                                                                  fluidRow(
                                                                                                                    column(6,
                                                                                                                           actionButton("Gene_Controls","Gene Controls"),
                                                                                                                           ##Action button to select control samples of the mRNA expression file
                                                                                                                           tags$head(tags$style(HTML('#Gene_Controls{width: 45%; }'))),tags$label(tags$style(HTML('#Gene_Controls {color:White; font-family:Open Sans; font-size: 18px;}')))
                                                                                                                    ),
                                                                                                                    column(6,
                                                                                                                           actionButton("Gene_Treatments","Gene Treatments"),
                                                                                                                           ##Action button to select treatment samples of the mRNA expression file
                                                                                                                           tags$head(tags$style(HTML('#Gene_Treatments{width:45%; }'))),tags$label(tags$style(HTML('#Gene_Treatments{color:White; font-family:Open Sans; font-size:18px;}'))))
                                                                                                                  )
                                                                                                 )
                                                                                       )
                                                                      )
                                                             ),
                                                             tabPanel("PhenoData", value="subTabPheno", icon = icon('file-medical-alt'), 
                                                                      conditionalPanel(condition= "input.ExampleData == 'Example Data'",
                                                                                       wellPanel(HTML('<p style="font-family:Open Sans; font-size:18px; text-decoration: underline; margin:5px; color:#550000; font-weight:bold"> phenoData-miRNA </p>'),
                                                                                                 fluidRow(
                                                                                                   column(12,
                                                                                                          div(style = 'overflow-x: scroll; font-family:Montserrat; font-size:120%;', DT::dataTableOutput("phenoDTmiRNA"))))
                                                                                       ),
                                                                                       wellPanel(HTML('<p style="font-family:Open Sans; font-size:18px; text-decoration: underline; margin:5px; color:#550000; font-weight:bold"> phenoData-mRNA </p>'),
                                                                                                 fluidRow(
                                                                                                   column(12,
                                                                                                          div(style = 'overflow-x: scroll; font-family:Montserrat; font-size:120%;', DT::dataTableOutput("phenoDTmRNA"))))
                                                                                       )
                                                                      )
                                                             )
                                                 ),
                                                 conditionalPanel(condition= "input.ShowExpressionData",
                                                                  ## Select the treatment and control samples and press the action button to visualize the the file that will be used to draw the network on Cytoscape.
                                                                  column(6,
                                                                         actionButton("BacktoSidePanel","Back"),
                                                                         tags$head(tags$style(HTML('#BacktoSidePanel{float:left; width:50%; border:15px;}'))),tags$label(tags$style(HTML('#BacktoSidePanel{color:White; font-family:Open Sans; font-size:22px;}')))
                                                                  )),
                                                 conditionalPanel(condition= "input.ShowExpressionData",
                                                                  ## Select the treatment and control samples and press the action button to visualize the the file that will be used to draw the network on Cytoscape.
                                                                  column(6,
                                                                         actionButton("ShowQuery","Forward"),
                                                                         tags$head(tags$style(HTML('#ShowQuery{float:right; width:65%; border:15px;}'))),tags$label(tags$style(HTML('#ShowQuery{color:White; font-family:Open Sans; font-size:22px;}')))
                                                                  )
                                                 )
                                        ),
                                        tabPanel('Differential Expression', value='panelDGEx', icon = icon('filter'),
                                                 wellPanel(id = 'miRNADGExWellPanel', HTML('<p style="font-family:Open Sans; font-size:18px; text-decoration: underline; margin:5px; color:#550000; font-weight:bold"> Differential Expression Results - miRNA </p>'),
                                                           fluidRow(
                                                             column(12,
                                                                    div(style = 'overflow-x: scroll; font-family:Montserrat; font-size:120%;', DT::dataTableOutput("DGExDTmiRNAUI")))),
                                                           fluidRow(
                                                             column(12,
                                                                    downloadButton('downloadLimmamiRNA', 'Download'),
                                                                    tags$head(tags$style(HTML('#downloadLimmamiRNA{width: 30%; }'))),tags$label(tags$style(HTML('#downloadLimmamiRNA{color:White; font-family:Open Sans; font-size: 18px;}')))))
                                                 ),
                                                 wellPanel(id = 'mRNADGExWellPanel', HTML('<p style="font-family:Open Sans; font-size:18px; text-decoration: underline; margin:5px; color:#550000; font-weight:bold"> Differential Expression Results - mRNA </p>'),
                                                           fluidRow(
                                                             column(12,
                                                                    div(style = 'overflow-x: scroll; font-family:Montserrat; font-size:120%;', DT::dataTableOutput("DGExDTmRNAUI")))),
                                                           fluidRow(
                                                             column(12,
                                                                    downloadButton('downloadLimmamRNA', 'Download'),
                                                                    tags$head(tags$style(HTML('#downloadLimmamRNA{width: 30%; }'))),tags$label(tags$style(HTML('#downloadLimmamRNA{color:White; font-family:Open Sans; font-size: 18px;}')))))
                                                 ),
                                                 conditionalPanel(condition= "input.ShowExpressionData",
                                                                  ## Select the treatment and control samples and press the action button to visualize the the file that will be used to draw the network on Cytoscape.
                                                                  column(6,
                                                                         actionButton("BacktoExprsPanel","Back"),
                                                                         tags$head(tags$style(HTML('#BacktoExprsPanel{float:left; width:50%; border:15px;}'))),tags$label(tags$style(HTML('#BacktoExprsPanel{color:White; font-family:Open Sans; font-size:22px;}')))
                                                                  )),
                                                 conditionalPanel(condition= "input.ShowExpressionData",
                                                                  ## Select the treatment and control samples and press the action button to visualize the the file that will be used to draw the network on Cytoscape.
                                                                  column(6,
                                                                         actionButton("ShowSIF","Forward"),
                                                                         tags$head(tags$style(HTML('#ShowSIF{float:right; width:65%; border:15px;}'))),tags$label(tags$style(HTML('#ShowSIF{color:White; font-family:Open Sans; font-size:22px;}')))
                                                                  )
                                                 )
                                        ),
                                        tabPanel("Network File",value="panelCyto", icon=icon('database'),
                                                 ## Tab for cytoscape input file
                                                 
                                                 conditionalPanel(condition=" $('html').hasClass('shiny-busy')",
                                                                  tags$div(style="float:right; padding-right:30px; padding-top:10px;",
                                                                           tags$img(src="achilles_loading.gif",height=75,width=75)),
                                                                  tags$div(style="float:right; padding-right:30px; padding-top:10px; font-family: 'Open Sans'; font-weight: 200; font-size:125%; margin:10px; color: #464646;", "Processing ...")
                                                 ),
                                                 
                                                 conditionalPanel(condition= "input.ShowSIF",
                                                                  wellPanel(HTML('<p style="font-family:Open Sans; font-size:18px; text-decoration: underline; margin:5px; color:#550000; font-weight:bold"> Network File </p>'),
                                                                            fluidRow(
                                                                              column(12,
                                                                                     div(style = 'overflow-x: scroll; font-family:Montserrat; font-size:120%;', DT::dataTableOutput("CytoscapeFileDT")))),
                                                                            fluidRow(
                                                                              column(12,
                                                                                     ## Visualize Cytoscape input file here
                                                                                     downloadButton('downloadCyto', 'Download'),
                                                                                     tags$head(tags$style(HTML('#downloadCyto{width: 30%; }'))),tags$label(tags$style(HTML('#downloadCyto{color:White; font-family:Open Sans; font-size: 18px;}')))
                                                                              )
                                                                            )
                                                                  ),
                                                                  conditionalPanel(condition= "input.ShowSIF",
                                                                                   fluidRow(
                                                                                     column(6,
                                                                                            actionButton("BacktoDGExPanel","Back"),
                                                                                            ## Enter initial parameters and files then press the action button to visualize input files
                                                                                            tags$head(tags$style(HTML('#BacktoDGExPanel{width:50%; border:15px;}'))),tags$label(tags$style(HTML('#BacktoDGExPanel{color:White; font-family:Open Sans; font-size: 22px;}')))),
                                                                                     column(6,
                                                                                            actionButton("GenerateNetwork","Generate Network"),
                                                                                            ## Press the action button to draw the network
                                                                                            tags$head(tags$style(HTML('#GenerateNetwork{width:65%; float:right; border:15px;}'))),tags$label(tags$style(HTML('#GenerateNetwork{color:White; font-family:Open Sans; font-size:22px;}')))
                                                                                     )
                                                                                   )
                                                                  )
                                                                  
                                                 )
                                        ),
                                        tabPanel('Network', value = 'panelNetwork', icon=icon('puzzle-piece'),
                                                 ## Network display is executed here
                                                 
                                                 conditionalPanel(condition=" $('html').hasClass('shiny-busy')",
                                                                  tags$div(style="float:right; padding-right:30px; padding-top:10px;",
                                                                           tags$img(src="achilles_loading.gif",height=75,width=75)),
                                                                  tags$div(style="float:right; padding-right:30px; padding-top:10px; font-family: 'Open Sans'; font-weight: 200; font-size:125%; margin:10px; color: #464646;", "Processing ...")
                                                 ),
                                                 
                                                 conditionalPanel(condition = 'input.GenerateNetwork',
                                                                  wellPanel( 
                                                                    fluidRow(
                                                                      column(12,
                                                                             dropdownButton(inputId = 'networkDropdown',
                                                                                            conditionalPanel(condition= "output.mirnaUI != 'No Data' & output.mrnaUI != 'No Data' & input.networkArrangement == 'Initial State'",
                                                                                                             wellPanel(id = 'networkCorrelationPhysicsPanel',
                                                                                                                       fluidRow(
                                                                                                                         column(12, 
                                                                                                                                checkboxInput('correlationCheckbox', HTML("<p style= 'font-family:Open Sans; color:#464646; font-size:20px; font-weight:500; padding-left:10px;'> Change node color by negatively/positively correlated edges </p>"), FALSE),
                                                                                                                                bsTooltip("correlationCheckbox", HTML('<p style="font-family:courier; background-color: black; font-size: 15px; color: white;"> Positively associate node pairs are shown in grey </p>'), placement = "top", trigger = "hover", options = list(container = "body")),
                                                                                                                                conditionalPanel(condition = "input.correlationCheckbox",
                                                                                                                                                 uiOutput('correlationSliderUI')
                                                                                                                                )
                                                                                                                         )
                                                                                                                       )
                                                                                                             )
                                                                                            ),
                                                                                            wellPanel(
                                                                                              fluidRow(
                                                                                                column(12,
                                                                                                       checkboxInput('physicsEngineNetwork', HTML("<p style= 'font-family:Open Sans; color:#464646; font-size:20px; font-weight:500; padding-left:10px;'> Physics engine active on network </p>"), FALSE),
                                                                                                       bsTooltip("physicsEngineNetwork", HTML('<p style="font-family:courier; background-color: black; font-size: 15px; color: white;"> Check to activate! Disabled physics engine increases the performance </p>'), placement = "top", trigger = "hover", options = list(container = "body"))
                                                                                                )
                                                                                              ),
                                                                                              fluidRow(
                                                                                                column(12,
                                                                                                       conditionalPanel(condition = "input.physicsEngineNetwork == false",
                                                                                                                        selectInput("networkLayout", HTML("<p style= 'font-family:Open Sans; color:#464646; font-size:20px; font-weight:500; padding-left:20px; padding-top:15px;'> Change network layout </p>"),
                                                                                                                                    choices= c("Star-Shape", "Reingold-Tilford", "Vertices On A Circle", "Nicely", "On Grid", "On Sphere", "Randomly", "Davidson-Harel", "Fruchterman-Reingold", "GEM force-directed", "Graphopt", "Kamada-Kawai", "Large Graph Layout", "Multidimensional Scaling", "Sugiyama"), selected= "Kamada-Kawai"),
                                                                                                                        bsTooltip("networkLayout", HTML('<p style="font-family:courier; background-color: black; font-size: 15px; color: white;"> Layout algorithms are implemented from igraph package of R </p>'), placement = "top", trigger = "hover", options = list(container = "body")),
                                                                                                       )
                                                                                                )
                                                                                              )
                                                                                            ),
                                                                                            conditionalPanel(condition = "input.Prioritize",
                                                                                                             wellPanel(id = "networkLayoutPanel",
                                                                                                                       fluidRow(
                                                                                                                         column(12, 
                                                                                                                                selectInput("networkArrangement", HTML("<p style= 'font-family:Open Sans; color:#464646; font-size:20px; font-weight:500; padding-left:20px; padding-top:15px;'> Change network arrangement </p>"),
                                                                                                                                            choices= c("Initial State", "Prioritized"), selected = "Initial State"))
                                                                                                                         
                                                                                                                       )
                                                                                                             )
                                                                                            ),
                                                                                            wellPanel(
                                                                                              fluidRow(
                                                                                                column(12,
                                                                                                       downloadButton('downloadNetwork', 'Export as .html'),
                                                                                                       tags$head(tags$style(HTML('#downloadNetwork{width: 100%;}'))),tags$label(tags$style(HTML('#downloadNetwork{color:White; font-family:Open Sans; font-size: 16px;}'))),
                                                                                                       bsTooltip("downloadNetwork", HTML('<p style="font-family:courier; background-color: black; font-size: 15px; color: white;"> Click to export network as an interactive .html file </p>'), placement = "top", trigger = "hover", options = list(container = "body")),
                                                                                                ),
                                                                                                column(12,
                                                                                                       downloadButton('downloadNetworkInteractions', 'Network Interactions'),
                                                                                                       tags$head(tags$style(HTML('#downloadNetworkInteractions{width: 100%;}'))),tags$label(tags$style(HTML('#downloadNetworkInteractions{color:White; font-family:Open Sans; font-size: 16px;}'))),
                                                                                                       bsTooltip("downloadNetworkInteractions", HTML('<p style="font-family:courier; background-color: black; font-size: 15px; color: white;"> Click to download network interactions as a .txt file </p>'), placement = "top", trigger = "hover", options = list(container = "body")),
                                                                                                )
                                                                                              ),
                                                                                              fluidRow(
                                                                                                conditionalPanel(condition = 'input.ExampleData == "Example Data"',
                                                                                                                 column(12,
                                                                                                                        downloadButton('downloadNetworkProperties', 'Network Properties'),
                                                                                                                        tags$head(tags$style(HTML('#downloadNetworkProperties{width: 100%;}'))),tags$label(tags$style(HTML('#downloadNetworkProperties{color:White; font-family:Open Sans; font-size: 16px;}'))),
                                                                                                                        bsTooltip("downloadNetworkProperties", HTML('<p style="font-family:courier; background-color: black; font-size: 15px; color: white;"> Click to download network properties as a .txt file </p>'), placement = "top", trigger = "hover", options = list(container = "body")),
                                                                                                                 ),
                                                                                                                 column(12,
                                                                                                                        downloadButton('downloadNetworkPlot', 'Correlation Plot'),
                                                                                                                        tags$head(tags$style(HTML('#downloadNetworkPlot{width: 100%;}'))),tags$label(tags$style(HTML('#downloadNetworkPlot{color:White; font-family:Open Sans; font-size: 16px;}'))),
                                                                                                                        bsTooltip("downloadNetworkPlot", HTML('<p style="font-family:courier; background-color: black; font-size: 15px; color: white;"> Click to download Correlation ~ logFC plot as a .png file </p>'), placement = "top", trigger = "hover", options = list(container = "body"))
                                                                                                                 )
                                                                                                )
                                                                                              )
                                                                                            ),
                                                                                            wellPanel(
                                                                                              fluidRow(
                                                                                                column(12,
                                                                                                       sliderTextInput("sliderColorSpectrum", label =  HTML("<p style= 'font-family:Open Sans; color:#464646; font-size:20px; font-weight:500; padding-left:20px;'> Color spectrum of nodes by logFC values </p>"), choices = seq(from = -8, to = 8, by = 0.5), from_min = -8, from_max = -3, to_min = 3, to_max = 8, selected = c(-5,5)),
                                                                                                       bsTooltip("sliderColorSpectrum", HTML('<p style="font-family:courier; background-color: black; font-size: 15px; color: white;"> Nodes will be colored from blue to red by the interval chosen </p>'), placement = "top", trigger = "hover", options = list(container = "body"))
                                                                                                       ##The end values of slider represents the expression value difference in which nodes should be colored with specific colors
                                                                                                       ## i.e. Red for upregulation and Blue for downregulation and spectrum in between for middle values.
                                                                                                )
                                                                                              )
                                                                                            ),
                                                                                            circle = TRUE, status = "danger",
                                                                                            icon = icon("gear"), width = "300px",
                                                                                            right = TRUE,
                                                                                            tooltip = tooltipOptions(title = "Click to see network options !", placement = 'left')
                                                                             ),
                                                                             tags$head(tags$style(HTML('#networkDropdown{float:right; border:15px;}')))
                                                                      )
                                                                    ),
                                                                    fluidRow(
                                                                      column(12,
                                                                             visNetworkOutput('visNetworkUI', height= '900px')
                                                                      )
                                                                    )
                                                                  ),
                                                                  fluidRow(
                                                                    column(12,
                                                                           actionButton("BacktoCyto","Back"),
                                                                           tags$head(tags$style(HTML('#BacktoCyto{width:30%; float:left; border:15px;}'))),tags$label(tags$style(HTML('#BacktoCyto{color:White; font-family:Open Sans; font-size:22px;}'))),
                                                                           
                                                                           actionButton("refresh","Restart", icon=icon('refresh')),
                                                                           tags$head(tags$style(HTML('#refresh{width:55%; float:right; border:15px;}'))),tags$label(tags$style(HTML('#refresh{color:White; font-family:Open Sans; font-size:22px;}')))
                                                                    )
                                                                  )
                                                 )
                                        ),
                                        tabPanel('Prioritization', value = 'panelKNGP', icon=icon('code-branch'),
                                                 wellPanel(HTML('<p style="font-family:Open Sans; font-size:20px; text-decoration: underline; margin:5px; color:#550000; font-weight:bold"> Prioritization Results </p>'),
                                                           fluidRow(
                                                             column(12,
                                                                    div(style = 'overflow-x: scroll; font-family:Montserrat; font-size:120%;', DT::dataTableOutput("KNGPUI"))
                                                             )
                                                           )
                                                 ),
                                                 wellPanel(HTML('<p style="font-family:Open Sans; font-size:20px; text-decoration: underline; margin:5px; color:#550000; font-weight:bold"> Prioritization Info </p>'),
                                                           fluidRow(
                                                             column(6,
                                                                    htmlOutput('KNGPInfo')
                                                             )
                                                           )
                                                 ),
                                                 fluidRow(
                                                   column(12,
                                                          actionButton("BacktoNetwork","Back"),
                                                          tags$head(tags$style(HTML('#BacktoNetwork{width:30%; float:left; border:15px;}'))),tags$label(tags$style(HTML('#BacktoNetwork{color:White; font-family:Open Sans; font-size:22px;}')))
                                                   )
                                                 )
                                        )
                            )
                  )
)
)
