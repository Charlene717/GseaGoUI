





##### Load Packages  #####
  #### Basic installation ####
  ## Check whether the installation of those packages is required from basic
  Package.set <- c("tidyverse","ggplot2","Seurat","SeuratData","patchwork","plyr","eoffice","DT")
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
         tabPanel("Enrichment Analysis",
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
#
#                         # textInput("word_select", label = "Gene name","TP53"),
#                         # actionButton("RUNDEG", "DEG"),
#                         actionButton("RunOFL", "OFL"),
#                         actionButton("RunGSEA", "GSEA"),
#                         # actionButton("RUNGSEA", "ORA")
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
                                                     fluidPage(plotOutput("DistPlt"))
                                              ),
                                              column(3,
                                                     hr(),
                                                     textInput("PhenoColSet", label = "Group by Phenotype","sample_type"),
                                                     textInput("PhenoType1Set", label = "Group by Phen Type1","Recurrent Tumor"),
                                                     textInput("PhenoType2Set", label = "Group by Phen Type2","Primary Tumor"),
                                                     hr(),
                                                     actionButton(inputId="DistPlot", label="Run",
                                                                  icon=icon(name = "gears")) # https://fontawesomeicons.com/palette

                                              )
                                     ),
                                     tabPanel("Group by GeneExp",
                                                column(6,
                                                       fluidPage(plotOutput("DistPlt2"))
                                                ),
                                                column(3,
                                                       # br(),
                                                       hr(),
                                                       textInput("GeneNameSet", label = "Group by GeneExp","TP53"),
                                                       selectizeInput("GroupByGeneStats", label = "Group by GeneExp Thr",
                                                                      choices = list("Mean1SD" = "Mean1SD", "Quartiles" = "Quartiles"),
                                                                      selected = "Mean1SD"),
                                                       hr(),
                                                       actionButton(inputId="DistPlot2", label="Run",
                                                                    icon=icon(name = "gears")) # https://fontawesomeicons.com/palette
                                                )
                                              )
                          ),

                          navbarMenu("DEG setting",
                                     tabPanel("DEG Analysis",
                                              textInput("LogFCSet", label = "LogFC","1"),
                                              textInput("PvalueCSet", label = "Pval","0.05"),
                                              hr(),
                                              actionButton(inputId="RunDEG", label="Run",
                                                           icon=icon(name = "gears")) # https://fontawesomeicons.com/palette
                                     ),
                                     tabPanel("Advance setting",
                                              textInput("XXXX2", label = "Gene name","TP53"),
                                     )

                          ),

                          tabPanel("GSEA setting",
                                   actionButton("RunOFL", "OFL"),
                                   actionButton("RunGSEA", "GSEA"),


                          ),

                          tabPanel("GO setting",
                                   actionButton("RunOFL", "OFL"),
                                   actionButton("RunGO", "GO"),


                          )
                        )

                      ),


                      position = c("left")
                    ),

                    ##### Summary Page #####
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

                      navbarMenu("DEG Analysis",
                                 tabPanel("DEG MTX",
                                          fluidPage(fluidRow(dataTableOutput("DEGMTX")))
                                 ),
                                 tabPanel("Volcano plot",
                                          fluidPage(fluidRow(plotOutput("VolcanoPlot")))
                                 )
                      ),

                      ##### Analysis search Page #####
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
                      ##### Analysis search Page #####
                      navbarMenu("GO Analysis",
                                 tabPanel("Bar plot",
                                          fluidPage(plotOutput("GOBarPlot",height = "800px", width = "1500px"))
                                 ),
                                 tabPanel("Dot plot",
                                          fluidPage(plotOutput("GODotPlot",height = "800px", width = "1000px"))
                                 ),
                                 tabPanel("UpSet Plot",
                                          fluidPage(plotOutput("GOUpSetPlot",height = "800px", width = "1500px"))
                                 )
                      )
                    )
                  )
         ),
###############################################################################################################################################
         tabPanel("Custom Gene Sets"),
###############################################################################################################################################
         tabPanel("MultiGroup Enrichment Analysis"),

###############################################################################################################################################
         tabPanel("Visualization"),
###############################################################################################################################################

         navbarMenu("About",
                    tabPanel("Summary"),
                    tabPanel("Tutorial"),
                    tabPanel("Citations"),
                    tabPanel("Release"),
                    tabPanel("Terms of Use"),
                    tabPanel("Contact"),
                    tabPanel("Sub-Component B")

                    )
  )



# ui =
#
#   fluidPage(
#     # https://stackoverflow.com/questions/57037758/r-shiny-how-to-color-margin-of-title-panel
#     titlePanel(h1("GseaGoUI",
#                   style='background-color:#e6d5f2;
#                          color:#474973;
#                          font-weight: 500;
#                          font-family: Arial Black;
#                          line-height: 1.2;
#                          padding-left: 15px'
#                   )
#                ),
#
#     sidebarLayout(
#       sidebarPanel(
#         fileInput("File_GeneExp", "Choose GeneExp File", accept = ".tsv", multiple = T),
#         fileInput("File_Anno", "Choose Annotation File", accept = ".tsv", multiple = T),
#         fileInput("File_GeneSet", "Choose GeneSet Files", accept = ".txt", multiple = F),
#
#         # textInput("word_select", label = "Gene name","TP53"),
#         # actionButton("RUNDEG", "DEG"),
#         actionButton("RunOFL", "OFL"),
#         actionButton("RunGSEA", "GSEA"),
#         # actionButton("RUNGSEA", "ORA")
#       ),
#
#       # mainPanel(textOutput(outputId="BestSent")),
#       mainPanel(
#         tabsetPanel(
#           navbarMenu("Data preprocessing",
#                      tabPanel("Fliter by phenotype",
#                               textInput("FliterByPhenotype", label = "Pheno section","sample_type"),
#                      ),
#                      tabPanel("Fliter by GeneExp",
#                               textInput("FliterByGeneExp", label = "Gene name","TP53"),
#                      )
#           ),
#
#           navbarMenu("Basic setting",
#                      tabPanel("Grouping",
#                           column(6,
#                             fluidPage(plotOutput("DistPlt"))
#                             ),
#                           column(3,
#                              selectizeInput("GroupMethod", label = "Group by",
#                                             choices = list("GeneExp" = "GroupByGene", "Phenotype" = "GroupByPheno"),
#                                             selected = "GeneExp"),
#                              # br(),
#                              # fluidPage(img(src = "Monocle3_UMAP.PNG",
#                              #               height = "450px", width = "1700px", align = "center"))
#                              hr(),
#                               textInput("GeneNameSet", label = "Group by GeneExp","TP53"),
#                               selectizeInput("GroupByGeneStats", label = "Group by GeneExp Thr",
#                                              choices = list("Mean1SD" = "Mean1SD", "Quartiles" = "Quartiles"),
#                                              selected = "Mean1SD"),
#                              hr(),
#                              actionButton(inputId="DistPlot", label="Run",
#                                           icon=icon(name = "gears")) # https://fontawesomeicons.com/palette
#                           ),
#                           column(3,
#                             hr(),
#                             textInput("PhenoColSet", label = "Group by Phenotype","sample_type"),
#                             textInput("PhenoType1Set", label = "Group by Phen Type1","Recurrent Tumor"),
#                             textInput("PhenoType2Set", label = "Group by Phen Type2","Primary Tumor"),
#                           ),
#                      ),
#                      tabPanel("DEG Analysis",
#                               textInput("LogFCSet", label = "LogFC","1"),
#                               textInput("PvalueCSet", label = "Pval","0.05"),
#                               hr(),
#                               actionButton(inputId="RunDEG", label="Run",
#                                            icon=icon(name = "gears")) # https://fontawesomeicons.com/palette
#                               )
#                     ),
#
#           navbarMenu("GSEA setting",
#                      tabPanel("Basic setting",
#                               textInput("XXXX1", label = "Pheno section","sample_type"),
#                               ),
#                      tabPanel("Advance setting",
#                               textInput("XXXX2", label = "Gene name","TP53"),
#                               ),
#                      tabPanel("GSEA",
#                               fluidPage(fluidRow(dataTableOutput("XXX3")))
#                      )
#
#                     ),
#
#           navbarMenu("Geneset setting",
#                      tabPanel("Com",
#                               textInput("word_selectXXX", label = "Gene name","TP53"),
#                               ),
#                      tabPanel("Fliter",
#                               textInput("word_selectXXX2", label = "Gene name","TP53"),
#                              )
#                     )
#         )
#
#       ),
#
#
#       position = c("left")
#     ),
#
#     ##### Summary Page #####
#     tabsetPanel(
#       # tabPanel("Summary",
#       #          fluidPage(
#       #           # plotOutput("HisFig"),
#       #           # tableOutput("SumTable"))
#       #            fluidPage(fluidRow(dataTableOutput("SumTable"))))
#       #          ),
#
#       ##### Text search Page #####
#       navbarMenu("View input",
#                  tabPanel("GeneExp",
#                           fluidPage(fluidRow(dataTableOutput("GeneExpOut")))
#                  ),
#                  tabPanel("Annotation",
#                           fluidPage(fluidRow(dataTableOutput("AnnoOut")))
#                  ),
#                  tabPanel("Gene sets",
#                           fluidPage(fluidRow(dataTableOutput("GenesetsOut")))
#                  )
#       ),
#
#       navbarMenu("DEG Analysis",
#                    tabPanel("DEG MTX",
#                             fluidPage(fluidRow(dataTableOutput("DEGMTX")))
#                           ),
#                    tabPanel("Volcano plot",
#                             fluidPage(fluidRow(plotOutput("VolcanoPlot")))
#                           )
#                 ),
#
#       ##### Analysis search Page #####
#       navbarMenu("GSEA Analysis",
#                  tabPanel("Bar plot",
#                           fluidPage(plotOutput("GSEABarPlot",height = "800px", width = "1500px"))
#                           ),
#                  tabPanel("Dot plot",
#                           fluidPage(plotOutput("GSEADotPlot",height = "800px", width = "1000px"))
#                          ),
#                  tabPanel("UpSet Plot",
#                           fluidPage(plotOutput("UpSetPlot",height = "800px", width = "1500px"))
#                  ),
#                  tabPanel("Enrichment plot",
#                           fluidPage(plotOutput("GseaPlot",height = "1000px", width = "1500px"))
#                  ),
#                  tabPanel("Overlay Enrichment plot",
#                           fluidPage(plotOutput("OverlayGseaPlot",height = "600px", width = "800px"))
#                  ),
#                  # tabPanel("SentFreq Dimension Reduction",
#                  #          fluidPage(img(src = "Monocle3_UMAP.PNG",
#                  #                        height = "450px", width = "1700px", align = "center"), br())
#                  # )
#               ),
#       ##### Analysis search Page #####
#       navbarMenu("GO Analysis",
#                  tabPanel("Bar plot",
#                           fluidPage(plotOutput("GOBarPlot",height = "800px", width = "1500px"))
#                  ),
#                  tabPanel("Dot plot",
#                           fluidPage(plotOutput("GODotPlot",height = "800px", width = "1000px"))
#                  ),
#                  tabPanel("UpSet Plot",
#                           fluidPage(plotOutput("GOUpSetPlot",height = "800px", width = "1500px"))
#                  )
#       )
#     )
#   )
