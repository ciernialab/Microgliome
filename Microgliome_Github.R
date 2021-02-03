# MICROGLIOME
# Updated - JK 1.27.2020

# resources consulted:
# https://rstats.wtf/debugging-r-code.html
# https://shiny.rstudio.com/articles/layout-guide.html#:~:text=To%20create%20rows%20within%20the,use%20the%20column()%20function.
# https://stackoverflow.com/questions/31573610/r-shiny-updating-checkboxgroupinput-based-on-another-checkboxgroupinput


## set working directory
setwd("/Users/jennkim/Desktop/Microgliome/Microgliome/")

## load libraries
library(shiny)
library(shinyWidgets)
library(trackViewer) # load package
library(GenomicFeatures) # load GenomicFeatures to create TxDb from UCSC
library(org.Mm.eg.db) # load annotation database
library(Mus.musculus)

## make mm10 dataset
mm10_genes <- as.data.frame(org.Mm.egSYMBOL)
mm10_GRanges <- genes(TxDb.Mmusculus.UCSC.mm10.knownGene, single.strand.genes.only=FALSE)
mm10_GRanges <- as.list(mm10_GRanges)

# UPDATE AS NEEDED: define datalists for checkboxes
Ayata2018 <- c("CB_H3K27me3_ChIP","STR_H3K27me3_ChIP",
               "CB_input","STR_input")
Gosselin2017 <- c("Rep1_ATAC","Rep1_H3K4me2_ChIP","Rep2_H3K4me2_ChIP",
                  "Rep1_H3K27ac_ChIP","Rep2_H3K27ac_ChIP","Rep2_PU1_ChIP",
                  "Pooled_input")
Gyoneva2019 <- c("Female_Cx3cr1KO_H3K27ac_ChIP", "Female_Cx3cr1KO_Rbp2_ChIP", "Female_Cx3cr1KO_input",
                 "Female_WT_H3K27ac_ChIP","Female_WT_Rbp2_ChIP","Female_WT_input",
                 "Male_Cx3cr1kO_H3K27ac_ChIP", "Male_Cx3cr1KO_Rbp2_ChIP", "Male_Cx3cr1KO_Input",
                 "Male_WT_H3K27ac_ChIP","Male_WT_Rbp2_ChIP","Male_WT_input")
TestRun <- c("H3K27ac", "HDAC12ko", "H3K9ac")

## UPDATE AS NEEDED: compile master datalist:
data_list <- c(Ayata2018, Gosselin2017, Gyoneva2019, TestRun)


###### MICROGLIOME SHINY CODE STARTS HERE ######

## ui ----
ui <- fluidPage(
  titlePanel("Microgliome"),
  
  fluidRow(
    column(6,
           browseTracksOutput('plot1')),
    column(6,
           browseTracksOutput('plot2'))
    
  ),
  
  fluidRow(
    column(6,
           browseTracksOutput('plot3')),
    column(6,
           browseTracksOutput('plot4'))
  ),
  
  hr(),
  
  fluidRow(
    
    #gene search
    column(3,
           textInput("geneSearch",
                     label = h3("Type in gene of interest"),
                     value = "Abca2"),
           actionButton("add_graph", "Add Graph"),
           actionButton("reset_graph", "Reset Graphs"),
           br(),
           br()),
    column(3,
           checkboxGroupInput("Ayata2018",
                              label = h3("Ayata et al., 2018"), 
                              choices = Ayata2018),
           checkboxGroupInput("Gosselin2017", 
                              label = h3("Gosselin et al., 2017"), 
                              choices = Gosselin2017),
           br(),
           br()),
    column(3,
           checkboxGroupInput("Gyoneva2019",
                              label = h3("Gyoneva et al., 2019"),
                              choices = Gyoneva2019),
           checkboxGroupInput("TestRun",
                              label = h3("Test Run"),
                              choices = TestRun),
           br(),
           br())
    
    #range search
    #column(3,
    #       textInput("range",
    #                 label = h3("Indicate genomic range of interest"),
    #                 value = "chr3:108,276,000-108,285,000")),
    #dataset(s) selection
  )
)

## server ----
server <- function(input,output){
  
  current_graph <- reactive({
    print("a")
    
    query_gene <- intersect(input$geneSearch, mm10_genes$symbol)
    query_geneID <- mm10_genes$gene_id[mm10_genes$symbol==query_gene]
    #mm10_GRanges is a list of GRanges objects
    gr <- mget(query_geneID,as.environment(mm10_GRanges))
    gr <- makeGRangesFromDataFrame(gr)
    print("b")
    
    trs <- geneModelFromTxdb(TxDb.Mmusculus.UCSC.mm10.knownGene, org.Mm.eg.db, gr=gr)
    
    chooseData_output <- c(input$Ayata2018, input$Gosselin2017, input$Gyoneva2019, input$TestRun)
    #chooseData_output <- as.character(chooseData_output)
    #data_list = character, 
    #TracksToShow <- mget(data_list,as.environment(chooseData_output))
    TracksToShow <- intersect(chooseData_output, data_list)
    print("ab")
    #TracksToHide <- setdiff(data_list, chooseData_output)

    
    ## UPDATE AS NEEDED: create trackViewer-compatible datalists
    
    # test run
    H3K27ac = importScore(file.path("/Users/jennkim/Desktop/Microgliome/Microgliome/", "data/H3K27ac_Repl1.bw"), format="BigWig", ranges=gr)
    HDAC12ko = importScore(file.path("/Users/jennkim/Desktop/Microgliome/Microgliome/", "data/HDAC12KO_H3K27ac_Repl3.bw"), format="BigWig", ranges=gr)
    H3K9ac = importScore(file.path("/Users/jennkim/Desktop/Microgliome/Microgliome/", "data/WT_H3K9ac_Repl2.bw"), format="BigWig", ranges=gr)
    TestRun_trackViewer <- c(H3K27ac, HDAC12ko, H3K9ac)
    
    print("ab1")
    
    #Ayata 2018 - # are not UCSC-compatible? double-check
    CB_H3K27me3_ChIP = importScore(file.path("/Volumes/Seagate Backup Plus Drive/Microgliome/Data/", "Ayata2018/cbMG_H3K27me3.bw"), format = "BigWig", ranges=gr)
    STR_H3K27me3_ChIP = importScore(file.path("/Volumes/Seagate Backup Plus Drive/Microgliome/Data/", "Ayata2018/StrMG_H3K27me3.bw"), format = "BigWig", ranges=gr)
    #CTX_H3K27me3_ChIP = importScore(file.path("/Volumes/Seagate Backup Plus Drive/Microgliome/Data/", "Ayata2018/cxMg_H3K27me3.bw"), format = "BigWig", ranges=gr)
    CB_input = importScore(file.path("/Volumes/Seagate Backup Plus Drive/Microgliome/Data/", "Ayata2018/cbMG_Input.bw"), format = "BigWig", ranges=gr)
    STR_input = importScore(file.path("/Volumes/Seagate Backup Plus Drive/Microgliome/Data/", "Ayata2018/StrMG_Input.bw"), format = "BigWig", ranges=gr)
    #CTX_input = importScore(file.path("/Volumes/Seagate Backup Plus Drive/Microgliome/Data/", "Ayata2018/cxMg_Input.bw"), format = "BigWig", ranges=gr)
    Ayata2018_trackViewer <- c(CB_H3K27me3_ChIP, STR_H3K27me3_ChIP, CB_input, STR_input)
    
    print("ab12")
    
    #Gosselin 2017 - # are not UCSC-compatible? double-check
    Rep1_ATAC = importScore(file.path("/Volumes/Seagate Backup Plus Drive/Microgliome/Data/", "Gosselin2017/ATAC_Repl1.bw"), format = "BigWig", ranges=gr)
    #Rep2_ATAC = importScore(file.path("/Volumes/Seagate Backup Plus Drive/Microgliome/Data/", "Gosselin2017/ATAC_Repl2.bw"), format = "BigWig", ranges=gr)
    Rep1_H3K4me2_ChIP = importScore(file.path("/Volumes/Seagate Backup Plus Drive/Microgliome/Data/", "Gosselin2017/H3K4me2_Repl1.bw"), format = "BigWig", ranges=gr)
    Rep2_H3K4me2_ChIP = importScore(file.path("/Volumes/Seagate Backup Plus Drive/Microgliome/Data/", "Gosselin2017/H3K4me2_Repl2.bw"), format = "BigWig", ranges=gr)
    Rep1_H3K27ac_ChIP = importScore(file.path("/Volumes/Seagate Backup Plus Drive/Microgliome/Data/", "Gosselin2017/H3K27ac_Repl1.bw"), format = "BigWig", ranges=gr)
    Rep2_H3K27ac_ChIP = importScore(file.path("/Volumes/Seagate Backup Plus Drive/Microgliome/Data/", "Gosselin2017/H3K27ac_Repl2.bw"), format = "BigWig", ranges=gr)
    #Rep1_PU1_ChIP = importScore(file.path("/Volumes/Seagate Backup Plus Drive/Microgliome/Data/", "Gosselin2017/PU1_Repl1.bw"), format = "BigWig", ranges=gr)
    Rep2_PU1_ChIP = importScore(file.path("/Volumes/Seagate Backup Plus Drive/Microgliome/Data/", "Gosselin2017/PU1_Repl2.bw"), format = "BigWig", ranges=gr)
    Pooled_input = importScore(file.path("/Volumes/Seagate Backup Plus Drive/Microgliome/Data/", "Gosselin2017/Input_Pool.bw"), format = "BigWig", ranges=gr)
    Gosselin2017_trackViewer <- c(Rep1_ATAC, Rep1_H3K4me2_ChIP, Rep2_H3K4me2_ChIP, Rep1_H3K27ac_ChIP, Rep2_H3K27ac_ChIP,
                                  Rep2_PU1_ChIP, Pooled_input)
    
    print("ab123")
    
    #Gyoneva 2019
    Female_Cx3cr1KO_H3K27ac_ChIP = importScore(file.path("/Volumes/Seagate Backup Plus Drive/Microgliome/Data/", "Gyoneva2019/Cx3cr1KO_H3k27ac_F.bw"), format = "BigWig", ranges=gr)
    Female_Cx3cr1KO_Rbp2_ChIP = importScore(file.path("/Volumes/Seagate Backup Plus Drive/Microgliome/Data/", "Gyoneva2019/Cx3cr1KO_Rpb2_F.bw"), format = "BigWig", ranges=gr)
    Female_Cx3cr1KO_input = importScore(file.path("/Volumes/Seagate Backup Plus Drive/Microgliome/Data/", "Gyoneva2019/Cx3cr1KO_Input_F.bw"), format = "BigWig", ranges=gr)
    Female_WT_H3K27ac_ChIP = importScore(file.path("/Volumes/Seagate Backup Plus Drive/Microgliome/Data/", "Gyoneva2019/WT_H3k27ac_F.bw"), format = "BigWig", ranges=gr)
    Female_WT_Rbp2_ChIP = importScore(file.path("/Volumes/Seagate Backup Plus Drive/Microgliome/Data/", "Gyoneva2019/WT_Rpb2_F.bw"), format = "BigWig", ranges=gr)
    Female_WT_input = importScore(file.path("/Volumes/Seagate Backup Plus Drive/Microgliome/Data/", "Gyoneva2019/WT_Input_F.bw"), format = "BigWig", ranges=gr)
    
    Male_Cx3cr1kO_H3K27ac_ChIP = importScore(file.path("/Volumes/Seagate Backup Plus Drive/Microgliome/Data/", "Gyoneva2019/Cx3cr1KO_H3k27ac_M.bw"), format = "BigWig", ranges=gr)
    Male_Cx3cr1KO_Rbp2_ChIP = importScore(file.path("/Volumes/Seagate Backup Plus Drive/Microgliome/Data/", "Gyoneva2019/Cx3cr1KO_Rpb2_M.bw"), format = "BigWig", ranges=gr)
    Male_Cx3cr1KO_Input = importScore(file.path("/Volumes/Seagate Backup Plus Drive/Microgliome/Data/", "Gyoneva2019/Cx3cr1KO_Input_M.bw"), format = "BigWig", ranges=gr)
    Male_WT_H3K27ac_ChIP = importScore(file.path("/Volumes/Seagate Backup Plus Drive/Microgliome/Data/", "Gyoneva2019/WT_H3k27ac_M.bw"), format = "BigWig", ranges=gr)
    Male_WT_Rbp2_ChIP = importScore(file.path("/Volumes/Seagate Backup Plus Drive/Microgliome/Data/", "Gyoneva2019/WT_Rpb2_M.bw"), format = "BigWig", ranges=gr)
    Male_WT_input = importScore(file.path("/Volumes/Seagate Backup Plus Drive/Microgliome/Data/", "Gyoneva2019/WT_Input_M.bw"), format = "BigWig", ranges=gr)
    
    Gyoneva2019_trackViewer <- c(Female_Cx3cr1KO_H3K27ac_ChIP, Female_Cx3cr1KO_Rbp2_ChIP, Female_Cx3cr1KO_input,
                                 Female_WT_H3K27ac_ChIP,Female_WT_Rbp2_ChIP,Female_WT_input,
                                 Male_Cx3cr1kO_H3K27ac_ChIP, Male_Cx3cr1KO_Rbp2_ChIP, Male_Cx3cr1KO_Input,
                                 Male_WT_H3K27ac_ChIP,Male_WT_Rbp2_ChIP,Male_WT_input)
    print("ab1234")
    
    ## UPDATE AS NEEDED: master list of trackViewer datasets if you need to:
    data_list_trackViewer <- c(TestRun_trackViewer, Ayata2018_trackViewer, Gosselin2017_trackViewer, Gyoneva2019_trackViewer)
    names(data_list_trackViewer) <- data_list
    
    data <- mget(TracksToShow,as.environment(data_list_trackViewer))
    print("c")
    
    ## create styles by preset theme
    optSty <- optimizeStyle(trackList(trs, data), theme="col")
    trackList <- optSty$tracks
    viewerStyle <- optSty$style
    ## adjust the styles for this track
    ### rename the trackList for each track
    names(trackList)[1:3] <- paste0("Sort1: ", names(trackList)[1:3])
    print(names(trackList))
    #names(trackList)[4] <- input$chooseData
    print(input$chooseData)
    ### change the lab positions for gene model track to bottomleft
    setTrackStyleParam(trackList[[1]], "ylabpos", "bottomleft")
    setTrackStyleParam(trackList[[2]], "ylabpos", "bottomleft")
    setTrackStyleParam(trackList[[3]], "ylabpos", "bottomleft")
    ### change the color of gene model track
    setTrackStyleParam(trackList[[1]], "ylabgp", list(cex=1, col="red"))
    setTrackStyleParam(trackList[[2]], "ylabgp", list(cex=1, col="green"))
    setTrackStyleParam(trackList[[2]], "ylabgp", list(cex=1, col="blue"))
    ### remove the xaxis
    setTrackViewerStyleParam(viewerStyle, "xaxis", FALSE)
    
    browseTracks(trackList, gr = gr, viewerStyle = viewerStyle)
    
  })
  
  
  vals <- reactiveValues(p = 1, vp1 = NULL, vp2 = NULL, vp3 = NULL, vp4 = NULL)
  
  observeEvent(
    input$add_graph,
    {
      req(vals$p <= 4)
      switch (as.character(vals$p),
              "1" = {vals$vp1 = current_graph()},
              "2" = {vals$vp2 = current_graph()},
              "3" = {vals$vp3 = current_graph()},
              "4" = {vals$vp4 = current_graph()}
      )
      vals$p <- vals$p + 1
    }
  )
  
  observeEvent(
    input$reset_graph,
    {
      vals$p <- 1
      vals$vp1 <- NULL
      vals$vp2 <- NULL
      vals$vp3 <- NULL
      vals$vp4 <- NULL
    }
  )
  
  output$plot1 <- renderbrowseTracks({
    print("d")
    vals$vp1
  })
  
  output$plot2 <- renderbrowseTracks({
    print("e")
    vals$vp2
    
  })
  output$plot3 <- renderbrowseTracks({
    print("f")
    vals$vp3
  })
  
  output$plot4 <- renderbrowseTracks({
    print("g")
    vals$vp4
  })
}

# Run app ----
shinyApp(ui, server)






