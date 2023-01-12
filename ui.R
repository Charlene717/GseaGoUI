##### Load Packages  #####
source("FUN_Package_InstLoad.R")
FUN_Basic.set <- c("tidyverse","ggplot2","Seurat","SeuratData","patchwork","plyr","eoffice","DT","shiny","shinyFiles")
FUN_BiocManager.set <- c("clusterProfiler","enrichplot","pathview")
## Set the desired organism
# organism = "org.Hs.eg.db" ## c("org.Hs.eg.db","org.Mm.eg.db","org.Dm.eg.db")
# c(organism,"fgsea")

FUN_Package_InstLoad(Basic.set = FUN_Basic.set, BiocManager.set = FUN_BiocManager.set)
rm(FUN_Basic.set, FUN_BiocManager.set)

##### Function setting #####
  ## Call function
  source("FUN_GSEA_ANAL.R"); source("FUN_GSEA_LargeGeneSet.R"); source("FUN_GSEA_ggplot.R"); source("FUN_GSEA_ForOFFL.R");
  source("FUN_Find_Markers.R"); source("FUN_DEG_Analysis.R"); source("FUN_VolcanoPlot.R");
  source("FUN_ggPlot_vline.R"); source("FUN_DistrPlot.R"); source("FUN_Group_GE.R"); source("FUN_Beautify_ggplot.R");

##### UI ########
ui <- navbarPage(
                 # p(strong(code("GseaGoUI", style='background-color:#a7c1d9; color:#0c3d6b; font-family: Arial Black'))),
                 p(strong(img(src = "GSEAGOUI.png", height = "90px", width = "85px",align = "center"))),
                 # h1("GseaGoUI", style='background-color:#a7c1d9; color:#0c3d6b; font-family: Arial Black' ),

###############################################################################################################################################
        #####********************** Enrichment Analysis **********************#####
        navbarMenu
        ("Enrichment Analysis",
          # img(src = "GSEAGOUI.png", height = "80px", width = "75px", align = "center"),
          #####********************** Basic Enrichment Analysis **********************#####
          tabPanel
          ( "Basic",

            fluidPage
            (
              # https://stackoverflow.com/questions/57037758/r-shiny-how-to-color-margin-of-title-panel
              titlePanel(
                          h3(br(),"GseaGoUI"),
                          # h1("GseaGoUI", style='background-color:#e6d5f2; color:#474973; font-weight: 500; font-family: Arial Black; line-height: 1.2; padding-left: 15px')
                        ),

              sidebarLayout
              (
                #####*********** Input Page ***********#####
                sidebarPanel
                (
                  fileInput("File_GeneExp", "Choose GeneExp File", accept = ".tsv", multiple = T),
                  fileInput("File_Anno", "Choose Annotation File", accept = ".tsv", multiple = T),
                  fileInput("File_GeneSet", "Choose GeneSet Files", accept = ".txt", multiple = F),
                  hr(),
                  # https://stackoverflow.com/questions/39196743/interactive-directory-input-in-shiny-app-r
                  tags$div
                  (
                    h2
                    ("Choose the path to save the results",
                       style= 'color:#474973; font-size: 1.7rem; line-height: 1.7rem;
                               font-family: Arial Black; padding-left: 0px'
                    ),
                    tags$label
                    ("Save Path", class="btn btn-primary",
                      tags$input(id = "fileIn", webkitdirectory = TRUE, type = "file", style="display: none;", onchange="pressed()")
                    )
                  ),
                ),

                #####*********** Setting Page ***********#####
                # mainPanel(textOutput(outputId="BestSent")),
                mainPanel
                (
                  tabsetPanel
                  (
                    #####*********** Data preprocessing ***********#####
                    tabPanel
                    ("Data preprocessing",
                      column
                      ( 4, h3("Fliter by Phenotype"),
                        textInput("FltByPhenoSet1", label = "Phenotype section","gender"),
                        textInput("FltByPhenoSet2", label = "Group section","FEMALE"), # c("Recurrent Tumor", "Primary Tumor")
                        h3("Fliter by Gene Expression"),
                        textInput("FltByGESetName", label = "Gene Name","TP53"),
                        column(5,textInput("FltByGESetMax", label = "Max GeneExp","10"),),
                        column(5,textInput("FltByGESetMin", label = "Min GeneExp","2"),),
                        br(),
                        actionButton(inputId="ActFltSet", label="Filter", icon=icon(name = "filter-circle-xmark")) # https://fontawesomeicons.com/
                      ),
                      column
                      ( 3,h3("UpdataGeneName"),
                        # radioButtons("UpdataGeneName", h3("Updata Gene Name"),
                        #               choices = list("Yes" = 1, "No" = 2), selected = 1),
                        actionButton(inputId="ActUPDGName", label="Update", icon=icon(name = "refresh")) # https://fontawesomeicons.com/
                      )

                    ),

                    #####*********** Group setting ***********#####
                    navbarMenu
                    ("Group setting",
                       tabPanel
                       ("Group by Pheno",
                         column(6, fluidPage(plotOutput("DistPlt2"))),
                         column
                         (3, h3("Group by Phenotype"),
                          hr(),
                          textInput("GrpSetPhenoCol", label = "Phenotype","sample_type"),
                          textInput("GrpSetPhenoType1", label = "Group 1","Recurrent Tumor"),
                          textInput("GrpSetPhenoType2", label = "Group 2","Primary Tumor"),
                          hr(),
                          actionButton(inputId="DistPlot2", label="See Dist", icon=icon(name = "photo")), # https://fontawesomeicons.com/
                          actionButton(inputId="CheckPheno", label="Check", icon=icon(name = "check-square-o")) # https://fontawesomeicons.com/
                         )
                       ),
                       tabPanel
                       ("Group by GeneExp",
                         column(6, fluidPage(plotOutput("DistPlt"))),
                         column
                         (4, # br(),
                          h3("Group by Gene Expression"),
                          hr(),
                          textInput("GeneNameSet", label = "Gene name","CLU"),
                          selectizeInput
                          ("GroupByGeneStats", label = "Cutoff of GeneExp",
                            choices = list("Mean" = "Mean", "Mean+1SD" = "Mean1SD", "Mean+2SD" = "Mean2SD", "Mean+3SD" = "Mean3SD",
                                           "Median" = "Median" , "Quartiles" = "Quartiles",
                                           "Customize Cutoff" = "CustomizeCutoff"),
                            selected = "Mean1SD"),
                          column(6, textInput("GrpSetGEUpCF", label = "Upper Cutoff","1")),
                          column(6, textInput("GrpSetGELwCF", label = "Lower Cutoff","1")),
                          column(12,hr()),
                          actionButton(inputId="DistPlot", label="See Dist", icon=icon(name = "photo")), # https://fontawesomeicons.com/
                          actionButton(inputId="CheckGE", label="Check", icon=icon(name = "check-square-o")) # https://fontawesomeicons.com/

                         )
                       )
                    ),

                    #####*********** DEG setting ***********#####
                    tabPanel
                    ("DEG setting",
                      h3("Filter"),
                      hr(),
                      textInput("LogFCSet", label = "LogFC Cutoff","1"),
                      textInput("PvalueSet", label = "P Value Cutoff","0.05"),
                      hr(),
                      actionButton(inputId="RunDEG", label="Run DEG", icon=icon(name = "gears")) # https://fontawesomeicons.com/
                    ),

                    #####*********** GSEA setting ***********#####
                    tabPanel
                    ("GSEA setting",
                      column
                      (6, h3("GSEA Analysis"),
                       hr(),
                       selectizeInput("GSEASet_Grp", label = "Group by",
                                       choices = list("Phenotype" = "GSEAGroupbyPheno",
                                                      "Gene Expression" = "GSEAGroupbyGeneExp"),
                                       selected = "GroupbyPheno"),

                       textInput("GSEASet_PermuNum", label = "Number of Permutations","1000"),
                       selectizeInput
                       ("pAdjustMethod", label = "pAdjust Method",
                         choices = list("BH Method" = "BH",
                                        "BY Method" = "BY",
                                        "bonferroni Method" = "bonferroni",
                                        "FDR Method" = "fdr",
                                        "hommel Method" = "hommel",
                                        "holm Method" = "holm",
                                        "hochberg Method" = "hochberg",
                                        "None" = "none"),
                         selected = "GSEASet_GeneSet"
                       ),
                       column(3, textInput("GSEASet_MinGSize", label = "Min geneset size","15")),
                       column(3, textInput("GSEASet_MaxGSize", label = "Max geneset size","500"))
                      ),
                      column
                      (6, h3("Visualization settings"),
                       hr(),
                       textInput("GSEASet_TopGS", label = "TOP Gene Sets","10"),
                       textInput("GSEASet_BottomGS", label = "Bottom Gene Sets","10"),
                       textInput("GSEASet_CustomizedGS", label = "Customized Gene Sets",""),
                       br(),
                       actionButton("RunOFLGSEA", "Official files", icon=icon(name = "file-download")),
                       actionButton("RunGSEA", "Run GSEA", icon=icon(name = "gears")) # https://fontawesomeicons.com/
                      )
                    ),

                    #####*********** ORA setting ***********#####
                    tabPanel
                    ("ORA setting",
                             column
                             (6, h3("ORA Enrichment Analysis"),
                              hr(),
                              selectizeInput("ORASet_Grp", label = "Group by",
                                              choices = list("Phenotype" = "ORAGroupbyPheno",
                                                             "Gene Expression" = "ORAGroupbyGeneExp"),
                                              selected = "GroupbyPheno"),
                              textInput("ORASet_MinOverlap", label = "Min Overlap","3"),
                              textInput("ORASet_PValue", label = "P Value Cutoff","0.05"),
                              textInput("ORASet_MinEnrich", label = "Min Enrichment","1.5"),
                             ),
                             column
                             (6, h3("Visualization settings"),
                              hr(),
                              textInput("GOSet_TopGS", label = "TOP Gene Sets","10"),
                              textInput("GOSet_BottomGS", label = "Bottom Gene Sets","10"),
                              textInput("GOSet_CustomizedGS", label = "Customized Gene Sets",""),
                              br(),
                              actionButton("RunOFLGO", "Official file", icon=icon(name = "file-download")),
                              actionButton("RunGO", "Run ORA", icon=icon(name = "gears")) # https://fontawesomeicons.com/
                             )
                    )
                    #####*********** *********** ***********#####
                  )
                ),
                #####*********** ************ ***********#####
                position = c("left")

              ),

              #####*********** Output Page ***********#####
              tabsetPanel
              (
                ##### View Input Page #####
                navbarMenu
                ("View Input",
                 tabPanel("GeneExp",fluidPage(fluidRow(dataTableOutput("GeneExpOut")))),
                 tabPanel("Annotation", fluidPage(fluidRow(dataTableOutput("AnnoOut")))),
                 tabPanel("Gene sets", fluidPage(fluidRow(dataTableOutput("GenesetsOut"))))
                ),

                ##### Data preprocessing Result Page #####
                navbarMenu
                ("Cleaned data",
                  tabPanel("GeneExp", br(),actionButton(inputId="SaveCDGE", label="SaveTSV", icon=icon(name = "file-text-o")), # https://fontawesomeicons.com/
                           hr(),
                           fluidPage(fluidRow(dataTableOutput("GeneExpCleanOut")))),
                  tabPanel("Annotation", br(),actionButton(inputId="SaveCDAnno", label="SaveTSV", icon=icon(name = "file-text-o")), # https://fontawesomeicons.com/
                           hr(),
                           fluidPage(fluidRow(dataTableOutput("AnnoCleanOut"))))
                ),

                ##### Analysis Result Page: DEG Analysis #####
                navbarMenu
                ("DEG Analysis",
                  tabPanel("DEG MTX",br(),actionButton(inputId="SaveDEG", label="SaveTSV", icon=icon(name = "file-text-o")), # https://fontawesomeicons.com/
                           hr(),
                           fluidPage(fluidRow(dataTableOutput("DEGMTX")))),
                  tabPanel("Volcano plot",br(),
                           actionButton(inputId="SaveVolcPDF", label="SavePDF", icon=icon(name = "file-pdf")), # https://fontawesomeicons.com/
                           actionButton(inputId="SaveVolcRData", label="SaveRData", icon=icon(name = "file-download")), # https://fontawesomeicons.com/
                           hr(),
                           fluidPage(fluidRow(plotOutput("VolcPlt",height = "500px", width = "500px"))))
                ),

                ##### Analysis Result Page: GSEA Analysis #####
                navbarMenu
                ("GSEA Analysis",
                  tabPanel("Bar plot", fluidPage(plotOutput("GSEABarPlot",height = "800px", width = "1500px"))),
                  tabPanel("Dot plot", fluidPage(plotOutput("GSEADotPlot",height = "800px", width = "1000px"))),
                  tabPanel("UpSet Plot", fluidPage(plotOutput("UpSetPlot",height = "800px", width = "1500px"))),
                  tabPanel("Enrichment plot", fluidPage(plotOutput("GseaPlot",height = "1000px", width = "1500px"))),
                  tabPanel("Overlay Enrichment plot", fluidPage(plotOutput("OverlayGseaPlot",height = "600px", width = "800px"))),
                  # # tabPanel("SentFreq Dimension Reduction", fluidPage(img(src = "Monocle3_UMAP.PNG", height = "450px", width = "1700px", align = "center"), br()))
                  # tabPanel("Bar plot", fluidPage(br(),img(src = "BarPlot.png", height = "700px", width = "600px", align = "center"),br())),
                  # tabPanel("Dot plot", fluidPage(br(),img(src = "Dotplot.png", height = "700px", width = "600px", align = "center"),br())),
                  # tabPanel("UpSet Plot", fluidPage(br(),img(src = "UpSetPlot.PNG", height = "700px", width = "900px", align = "center"),br())),
                  # tabPanel("Enrichment plot", fluidPage(br(),img(src = "Enrichment plot.PNG", height = "700px", width = "600px", align = "center"),br())),
                  # tabPanel("Overlay Enrichment plot", fluidPage(br(),img(src = "Overlay Enrichment plot.PNG", height = "700px", width = "600px", align = "center"),br())),
                  tabPanel("Gene Pathway Heatmap", fluidPage(br(),img(src = "GenePathwayHeatmap.PNG", height = "1400px", width = "1200px", align = "center"),br())),

                ),
                ##### Analysis Result Page: GO Analysis #####
                navbarMenu
                ("ORA Analysis",
                  tabPanel("Bar plot", fluidPage(plotOutput("GOBarPlot",height = "800px", width = "1500px"))),
                  tabPanel("Dot plot", fluidPage(plotOutput("GODotPlot",height = "800px", width = "1000px"))),
                  tabPanel("Network plot", fluidPage(plotOutput("GONetworkPlot",height = "800px", width = "1000px"))),
                  tabPanel("Network plot2", fluidPage(plotOutput("GONetworkPlot2",height = "800px", width = "1000px"))),
                  tabPanel("UpSet Plot", fluidPage(plotOutput("GOUpSetPlot",height = "800px", width = "1500px")))
                )
              )
              #####*********** *********** ***********#####
            ),

            # style = "background-color: #DEEBF7",skin = "blue",
            # tags$head(tags$style('.headerrow{height:8vh; background-color:#267dff}')),
            tags$style
            (HTML("
             .navbar-default .navbar-brand {color:white;}
             .navbar-default .navbar-brand:hover {color:white;}
             .navbar { background-color:#275682;}
             .navbar-default .navbar-nav > li > a {color:white;}
             .navbar-default .navbar-nav > .active > a,
             .navbar-default .navbar-nav > .active > a:focus,
             .navbar-default .navbar-nav > .active > a:hover {color:black;background-color:white;}
             .navbar-default .navbar-nav > li > a:hover {color:white;background-color:#275682;text-decoration}
            "))
          ),

          #####*********** MultiGroup Enrichment Analysis ***********#####
          tabPanel("MultiGroup")
        ),
###############################################################################################################################################
        #####********************** Custom Gene Sets **********************#####
        tabPanel("Custom Gene Sets"),
###############################################################################################################################################
        #####********************** Visualization **********************#####
        navbarMenu
        ("Visualization",
          tabPanel("From official GSEA results"),
          tabPanel("From official Metascape results"),
          tabPanel("From official GO results")
        ),
###############################################################################################################################################
        #####********************** About **********************#####
        navbarMenu
        ("About",
          tabPanel("Summary"),
          tabPanel("Tutorial"),
          tabPanel("Citations"),
          tabPanel("Release"),
          tabPanel("Terms of Use"),
          tabPanel("Contact")
        )

###############################################################################################################################################

  )



