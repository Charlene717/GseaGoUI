

##### Load Packages  #####
  #### Basic installation ####
  ## Check whether the installation of those packages is required from basic
  Package.set <- c("tidyverse","ggplot2","Seurat","SeuratData","patchwork","plyr","eoffice","DT","shiny","shinyFiles")
  for (i in 1:length(Package.set)) {
    if (!requireNamespace(Package.set[i], quietly = TRUE)){
      install.packages(Package.set[i])
    }
  }
  ## Load Packages
  lapply(Package.set, library, character.only = TRUE)
  rm(Package.set,i)


#### BiocManager installation ####
  ## Set the desired organism
  # organism = "org.Hs.eg.db" ## c("org.Hs.eg.db","org.Mm.eg.db")   ##  c("org.Dm.eg.db")

  ## Check whether the installation of those packages is required from BiocManager
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  Package.set <- c("clusterProfiler","enrichplot","pathview") # c(organism,"fgsea","clusterProfiler","enrichplot","pathview")
  for (i in 1:length(Package.set)) {
    if (!requireNamespace(Package.set[i], quietly = TRUE)){
      BiocManager::install(Package.set[i])
    }
  }
  ## Load Packages
  lapply(Package.set, library, character.only = TRUE)
  rm(Package.set,i)

  # options(stringsAsFactors = FALSE)


# Sys.setlocale(category = "LC_ALL", locale = "UTF-8")


##### Function setting #####
  ## Call function
  source("FUN_Beautify_ggplot.R")
  source("FUN_Find_Markers.R")
  source("FUN_VolcanoPlot.R")
  source("FUN_GSEA_LargeGeneSet.R")
  source("FUN_GSEA_ggplot.R")
  source("FUN_ggPlot_vline.R")

  source("FUN_DistrPlot.R")
  source("FUN_Group_GE.R")
  source("FUN_DEG_Analysis.R")
  source("FUN_GSEA_ANAL.R")
  source("FUN_GSEA_ForOFFL.R")

##### UI ########
  ui <- navbarPage("GseaGoUI",
###############################################################################################################################################
      navbarMenu("Enrichment Analysis",
#####*********** Basic Enrichment Analysis ***********#####
         tabPanel("Basic",
                  fluidPage(
                    # https://stackoverflow.com/questions/57037758/r-shiny-how-to-color-margin-of-title-panel
                    titlePanel(h1("GseaGoUI",
                                  style='background-color:#e6d5f2;
                                         color:#474973;
                                         font-weight: 500;
                                         font-family: Arial Black;
                                         line-height: 1.2;
                                         padding-left: 15px'
                                 )
                               ),

                    sidebarLayout(
                      sidebarPanel(
                        fileInput("File_GeneExp", "Choose GeneExp File", accept = ".tsv", multiple = T),
                        fileInput("File_Anno", "Choose Annotation File", accept = ".tsv", multiple = T),
                        fileInput("File_GeneSet", "Choose GeneSet Files", accept = ".txt", multiple = F),
                        hr(),
                        # https://stackoverflow.com/questions/39196743/interactive-directory-input-in-shiny-app-r
                        tags$div(h2("Choose the path to save the results",
                                    style= 'color:#474973;
                                            font-size: 1.7rem; line-height: 1.7rem;
                                            font-family: Arial Black;
                                            padding-left: 0px'),
                                   tags$label("Save Path", class="btn btn-primary",
                                   tags$input(id = "fileIn", webkitdirectory = TRUE, type = "file", style="display: none;", onchange="pressed()"))),
                      ),

                      # mainPanel(textOutput(outputId="BestSent")),
                      mainPanel(
                        tabsetPanel(
                          navbarMenu("Data preprocessing",
                                     tabPanel("Fliter by phenotype",
                                              textInput("FliterByPhenotype", label = "Pheno section","sample_type"),
                                     ),
                                     tabPanel("Fliter by GeneExp",
                                              textInput("FliterByGeneExp", label = "Gene name","TP53"),
                                     )
                          ),

                          navbarMenu("Group setting",
                                     tabPanel("Group by Pheno",
                                              column(6,
                                                     fluidPage(plotOutput("DistPlt2"))
                                              ),
                                              column(3,
                                                     h3("Group by Phenotype"),
                                                     hr(),
                                                     textInput("PhenoColSet", label = "Phenotype","sample_type"),
                                                     textInput("PhenoType1Set", label = "Group 1","Recurrent Tumor"),
                                                     textInput("PhenoType2Set", label = "Group 2","Primary Tumor"),
                                                     hr(),
                                                     actionButton(inputId="DistPlot2", label="Set group", icon=icon(name = "gears")) # https://fontawesomeicons.com/

                                              )
                                     ),
                                     tabPanel("Group by GeneExp",
                                                column(6,
                                                       fluidPage(plotOutput("DistPlt"))
                                                ),
                                                column(3,
                                                       # br(),
                                                       h3("Group by Gene Expression"),
                                                       hr(),
                                                       textInput("GeneNameSet", label = "Gene name","TP53"),
                                                       selectizeInput("GroupByGeneStats", label = "Cutoff of GeneExp",
                                                                      choices = list("Mean" = "Mean", "Mean+1SD" = "Mean1SD", "Mean+2SD" = "Mean2SD", "Mean+3SD" = "Mean3SD",
                                                                                     "Median" = "Median" , "Quartiles" = "Quartiles",
                                                                                     "Customize Cutoff" = "CustomizeCutoff"),
                                                                      selected = "Mean1SD"),
                                                       column(6,
                                                       textInput("UpBoundGeneExp", label = "Upper Cutoff","1")),
                                                       column(6,
                                                       textInput("LowBoundGeneExp", label = "Lower Cutoff","1")),
                                                       column(12,hr()),
                                                       actionButton(inputId="DistPlot", label="Set group", icon=icon(name = "gears")) # https://fontawesomeicons.com/
                                                )
                                              )
                          ),

                          navbarMenu("DEG setting",
                                     tabPanel("Basic setting",
                                              h3("Filter"),
                                              hr(),
                                              textInput("LogFCSet", label = "LogFC Cutoff","1"),
                                              textInput("PvalueCSet", label = "P Value Cutoff","0.05"),
                                              hr(),
                                              actionButton(inputId="RunDEG", label="Run DEG", icon=icon(name = "gears")) # https://fontawesomeicons.com/
                                     ),
                                     tabPanel("Advance setting",
                                              textInput("XXXX2", label = "Gene name","TP53"),
                                     ),
                                     tabPanel("Filter setting",
                                              textInput("XXXX3", label = "Gene name","TP53"),
                                     )

                          ),

                          tabPanel("GSEA setting",
                                   h3("GSEA Analysis"),
                                   hr(),
                                   textInput("GSEASet_NumPermu", label = "Number of Permutations","1000"),
                                   selectizeInput("GSEASet_PermuType", label = "Permutation type",
                                                  choices = list("PhenoType" = "GSEASet_PhenoType",
                                                                 "GeneSet" = "GSEASet_GeneSet"),
                                                  selected = "GSEASet_GeneSet"),
                                   textInput("GSEASet_MaxGSize", label = "Max geneset size","500"),
                                   textInput("GSEASet_MinGSize", label = "Min geneset size","15"),
                                   hr(),
                                   actionButton("RunOFL", "Official files", icon=icon(name = "fas fa-file-download")),
                                   actionButton("RunGSEA", "Run GSEA", icon=icon(name = "gears")) # https://fontawesomeicons.com/
                          ),

                          tabPanel("ORA setting",
                                   h3("ORA Enrichment Analysis"),
                                   hr(),
                                   textInput("ORASet_MinOverlap", label = "Min Overlap","3"),
                                   textInput("ORASet_PValue", label = "P Value Cutoff","0.05"),
                                   textInput("ORASet_MinEnrich", label = "Min Enrichment","1.5"),
                                   hr(),
                                   actionButton("RunOFL", "Official files", icon=icon(name = "fas fa-file-download")),
                                   actionButton("RunGO", "Run ORA", icon=icon(name = "gears")) # https://fontawesomeicons.com/
                          )
                        )

                      ),


                      position = c("left")
                    ),

                    ##### Summary Page ##################################################################
                    tabsetPanel(
                      # tabPanel("Summary",
                      #          fluidPage(
                      #           # plotOutput("HisFig"),
                      #           # tableOutput("SumTable"))
                      #            fluidPage(fluidRow(dataTableOutput("SumTable"))))
                      #          ),

                      ##### Text search Page #####
                      navbarMenu("View Input",
                                 tabPanel("GeneExp",
                                          fluidPage(fluidRow(dataTableOutput("GeneExpOut")))
                                 ),
                                 tabPanel("Annotation",
                                          fluidPage(fluidRow(dataTableOutput("AnnoOut")))
                                 ),
                                 tabPanel("Gene sets",
                                          fluidPage(fluidRow(dataTableOutput("GenesetsOut")))
                                 )
                      ),

                      ##### Data preprocessing Result Page #####
                      navbarMenu("Cleaned data",
                                 tabPanel("GeneExp",
                                          fluidPage(fluidRow(dataTableOutput("GeneExpCleanOut")))
                                 ),
                                 tabPanel("Annotation",
                                          fluidPage(fluidRow(dataTableOutput("AnnoCleanOut")))
                                 )
                      ),

                      ##### Analysis Result Page: DEG Analysis #####
                      navbarMenu("DEG Analysis",
                                 tabPanel("DEG MTX",
                                          fluidPage(fluidRow(dataTableOutput("DEGMTX")))
                                 ),
                                 tabPanel("Volcano plot",
                                          fluidPage(fluidRow(plotOutput("VolcanoPlot")))
                                 )
                      ),

                      ##### Analysis Result Page: GSEA Analysis #####
                      navbarMenu("GSEA Analysis",
                                 tabPanel("Bar plot",
                                          fluidPage(plotOutput("GSEABarPlot",height = "800px", width = "1500px"))
                                 ),
                                 tabPanel("Dot plot",
                                          fluidPage(plotOutput("GSEADotPlot",height = "800px", width = "1000px"))
                                 ),
                                 tabPanel("UpSet Plot",
                                          fluidPage(plotOutput("UpSetPlot",height = "800px", width = "1500px"))
                                 ),
                                 tabPanel("Enrichment plot",
                                          fluidPage(plotOutput("GseaPlot",height = "1000px", width = "1500px"))
                                 ),
                                 tabPanel("Overlay Enrichment plot",
                                          fluidPage(plotOutput("OverlayGseaPlot",height = "600px", width = "800px"))
                                 ),
                                 # tabPanel("SentFreq Dimension Reduction",
                                 #          fluidPage(img(src = "Monocle3_UMAP.PNG",
                                 #                        height = "450px", width = "1700px", align = "center"), br())
                                 # )
                      ),
                      ##### Analysis Result Page: GO Analysis #####
                      navbarMenu("ORA Analysis",
                                 tabPanel("Bar plot",
                                          fluidPage(plotOutput("GOBarPlot",height = "800px", width = "1500px"))
                                 ),
                                 tabPanel("Dot plot",
                                          fluidPage(plotOutput("GODotPlot",height = "800px", width = "1000px"))
                                 ),
                                 tabPanel("Network plot",
                                          fluidPage(plotOutput("GONetworkPlot",height = "800px", width = "1000px"))
                                 ),
                                 tabPanel("Network plot2",
                                          fluidPage(plotOutput("GONetworkPlot2",height = "800px", width = "1000px"))
                                 ),
                                 tabPanel("UpSet Plot",
                                          fluidPage(plotOutput("GOUpSetPlot",height = "800px", width = "1500px"))
                                 )
                      )
                    )
                  ),
         ),
#####*********** MultiGroup Enrichment Analysis ***********#####
         tabPanel("MultiGroup")
      ),
###############################################################################################################################################
         tabPanel("Custom Gene Sets"),
###############################################################################################################################################
         navbarMenu("Visualization",
                    tabPanel("From official GSEA results"),
                    tabPanel("From official Metascape results"),
                    tabPanel("From official GO results")
         ),
###############################################################################################################################################

         navbarMenu("About",
                    tabPanel("Summary"),
                    tabPanel("Tutorial"),
                    tabPanel("Citations"),
                    tabPanel("Release"),
                    tabPanel("Terms of Use"),
                    tabPanel("Contact")
                    )
  )



