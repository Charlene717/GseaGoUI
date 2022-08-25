rm(list = ls())
options(shiny.maxRequestSize=300*1024^2)

# ##### Current path and new folder setting* #####
# ProjectName = "TCGA"
# Sampletype = "LGG"
# TarGene_name = "TP53"
# Version = paste0(Sys.Date(),"_",ProjectName,"_",Sampletype,"_",TarGene_name)
# Save.Path = paste0(getwd(),"/",Version)
# ## Create new folder
# if (!dir.exists(Save.Path)){
#   dir.create(Save.Path)
# }

server <- function(input, output, session){
  ##### Reactive #####
    #### Main reactive ####
    GeneExp_ReA = reactive({
      # GeneExp.df <- input$File_GeneExp$datapath %>%
      #               read.table(header=T, row.names = 1, sep="\t")
      GeneExp.df <- input$File_GeneExp$datapath %>%
        read.csv(sep = "\t",
                 row.names = 1,
                 check.names = F) %>%
        as.data.frame()
    })

    AnnoOut_ReA = reactive({
      Anno.df <- input$File_Anno$datapath %>%
        read.csv(sep = "\t",
                 row.names = 1,
                 check.names = F) %>%
        as.data.frame()
    })

    GeneSet_ReA = reactive({
      GeneSet.df <- input$File_GeneSet$datapath %>%
        read.csv(sep = "\t",
                 header = F,
                 check.names = F) %>%
        as.data.frame()
    })

    #### Basic setting ####
    GeneName_ReA = reactive({
      Keyword = input$GeneNameSet
    })

  ##### Analysis #####
    #### Dist ####
    Plot.Dist <- eventReactive(c(input$DistPlot), {
      GeneExp.df <- input$File_GeneExp$datapath %>%
        read.csv(sep = "\t",
                 row.names = 1,
                 check.names = F) %>%
        as.data.frame()
      Plot.DistrPlot <- FUN_DistrPlot(GeneExp.df,
                                      TarGeneName = input$GeneNameSet, GroupMode = list(Mode="Mean",SD=1),
                                      Save.Path = Save.Path, SampleName = SampleName)
      # Plot.DistrPlot_SD_Q <- Plot.DistrPlot[["TGeneDen_SD_Q.p"]]
      # Plot.DistrPlot_SD_Q
      })

    #### DEG ####
    DE_Ext.df <- eventReactive(c(input$RunDEG), {
      GeneExp.df <- input$File_GeneExp$datapath %>%
        read.csv(sep = "\t",
                 row.names = 1,
                 check.names = F) %>%
        as.data.frame()
      Anno.df <- input$File_Anno$datapath %>%
        read.csv(sep = "\t",
                 row.names = 1,
                 check.names = F) %>%
        as.data.frame()

      # AnnoSet.lt <- list(GroupType = "sample_type", GroupCompare = c("Primary Tumor","Recurrent Tumor") )
      # Thr.lt <- list(LogFC = c("logFC",1), pVal = c("PValue",0.05) )

      DEG_ANAL.lt <- FUN_DEG_Analysis(GeneExp.df, Anno.df,
                                      GroupType = input$PhenoColSet ,GroupCompare = c(input$PhenoType1Set,input$PhenoType2Set),
                                      ThrSet = list(LogFC = c("logFC", input$LogFCSet),
                                                    pVal = c("PValue",input$PvalueCSet) ),
                                      TarGeneName = input$GeneNameSet, GroupMode = Mode_Group, SampleID = "X_INTEGRATION",
                                      Save.Path = Save.Path, SampleName = SampleName, AnnoName = "AvB")
      DE_Extract.df <- DEG_ANAL.lt[["DE_Extract.df"]]
      rownames(DE_Extract.df) <- seq(1:nrow(DE_Extract.df))
      DE_Extract.df
    })

    #### GSEA ####
    GSEAResult.lt <- eventReactive(c(input$RunGSEA), {
      GeneExp.df <- input$File_GeneExp$datapath %>%
        read.csv(sep = "\t",
                 row.names = 1,
                 check.names = F) %>%
        as.data.frame()
      Anno.df <- input$File_Anno$datapath %>%
        read.csv(sep = "\t",
                 row.names = 1,
                 check.names = F) %>%
        as.data.frame()

      # AnnoSet.lt <- list(GroupType = "sample_type", GroupCompare = c("Primary Tumor","Recurrent Tumor") )
      # Thr.lt <- list(LogFC = c("logFC",1), pVal = c("PValue",0.05) )

      DEG_ANAL.lt <- FUN_DEG_Analysis(GeneExp.df, Anno.df,
                                      GroupType = input$PhenoColSet ,GroupCompare = c(input$PhenoType1Set,input$PhenoType2Set),
                                      ThrSet = list(LogFC = c("logFC", input$LogFCSet),
                                                    pVal = c("PValue",input$PvalueCSet) ),
                                      TarGeneName = input$GeneNameSet, GroupMode = Mode_Group, SampleID = "X_INTEGRATION",
                                      Save.Path = Save.Path, SampleName = SampleName, AnnoName = "AvB")
      DE_Extract.df <- DEG_ANAL.lt[["DE_Extract.df"]]


      GeneSet.df <- input$File_GeneSet$datapath %>%
        read.csv(sep = "\t",
                 header = F,
                 check.names = F) %>%
        as.data.frame()

      GSEA_Result.lt <- FUN_GSEA_ANAL(DE_Extract.df, pathwayGeneSet = GeneSet.df,
                                      TarGeneName = input$GeneNameSet, GroupMode = Mode_Group,
                                      ThrSet = list(LogFC = c("logFC", input$LogFCSet),
                                                   pVal = c("PValue",input$PvalueCSet) ),
                                      Species = "Homo sapiens", # Speices type can check by msigdbr_species()
                                      Save.Path = Save.Path, SampleName = SampleName, AnnoName = "Path")
      GSEA_Result.lt

    })
  ##### Output #####
    #### View input data ####
    # output$GeneExpOut <- renderTable({
    output$GeneExpOut <- renderDataTable({
      if (length(input$File_GeneExp)>0){
        datatable(GeneExp_ReA() %>% as.data.frame(),
                  options = list(searchHighlight = TRUE, pageLength = 50))
      }else{
        print("Welcome to GseaGoUI!" %>% as.data.frame())
      }
      })

    output$AnnoOut <- renderDataTable({
      # GeneExp_ReA()[c(1:20),] %>% as.data.frame()
      datatable(AnnoOut_ReA() %>% as.data.frame(),
                options = list(searchHighlight = TRUE, pageLength = 50))
    })

    output$GenesetsOut <- renderDataTable({
      # GeneExp_ReA()[c(1:20),] %>% as.data.frame()
      datatable(GeneSet_ReA() %>% as.data.frame(),
                options = list(searchHighlight = TRUE, pageLength = 50))
    })

    #### Dist. ####
    output$DistPlt <- renderPlot({
      if (length(input$File_GeneExp)>0){
          Plot.Dist.lt <- Plot.Dist()
          Plot.Dist.Plot <- Plot.Dist.lt[["TGeneDen_SD_Q.p"]]
          Plot.Dist.Plot
      }else{
        print("Welcome to GseaGoUI!" %>% as.data.frame())
      }
    })
      # output$DistPlt <- renderPlot({
      #   Plot.Dist.lt <- Plot.Dist()
      #   Plot.Dist.Plot <- Plot.Dist.lt[["TGeneDen_SD_Q.p"]]
      #   Plot.Dist.Plot
      # })

    #### DEG ####
    output$DEGMTX <- renderDataTable({
      datatable(DE_Ext.df() %>% as.data.frame(),
                options = list(searchHighlight = TRUE, pageLength = 50))
    })

    #### GSEA ####
    output$GSEABarPlot = renderPlot({
      GSEAResult <- GSEAResult.lt()
      Plot.GSEABar.Plot <-  GSEAResult[["GSEABar_Plot"]]
      Plot.GSEABar.Plot
    })
    output$GSEADotPlot <- renderPlot({
      GSEAResult <- GSEAResult.lt()
      Plot.GSEADot.Plot <-  GSEAResult[["GSEADot_Plot"]]
      Plot.GSEADot.Plot
    })

    output$UpSetPlot <- renderPlot({
      GSEAResult <- GSEAResult.lt()
      Plot.GSEAUpSet.Plot <-  GSEAResult[["UpSet_Plot"]]
      Plot.GSEAUpSet.Plot
    })

    output$GseaPlot <- renderPlot({
      GSEAResult <- GSEAResult.lt()
      Plot.GSEAUpSet.Plot <-  GSEAResult[["Gsea_Plot"]]
      Plot.GSEAUpSet.Plot
    })

    output$OverlayGseaPlot <- renderPlot({
      GSEAResult <- GSEAResult.lt()
      Plot.GSEAUpSet.Plot <-  GSEAResult[["OverlayGsea_Plot"]]
      Plot.GSEAUpSet.Plot
    })

  }
