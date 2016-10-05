library(shiny)
library(readr)
library(dplyr)
library(ReactomePA)
library(DOSE)
library(reactome.db)
library(igraph)
library(graphite)
library(GOSemSim)
library(AnnotationDbi)
library(png)
library(grDevices)
probe_list <- read_csv("PROBELIST.csv")
kegg_list <- read_csv("kegglist.csv")
network_table <- read_csv("network_table.csv")
gene_names <- unique(probe_list$ENTREZ_GENE)
gene_names <- as.character(gene_names)
mygeneList <- rep(100, length(gene_names))
names(mygeneList) <- gene_names

shinyServer(function(input, output) {
  Data <- reactive({
    input$goButton
    isolate({ 
      search_term <- trimws(input$search_in)
      search_out <- probe_list[apply(probe_list, 1, FUN = 
                                       function(x, y = search_term) {any(grepl(y, x, ignore.case = TRUE))}), ]
      search_out <- search_out %>% dplyr::select(CHEMICAL_ID, TARGET_NAME, UNIPROTKB, DATA_SOURCE)
      if(search_term == "Search Input"){return()}
      return(search_out)
    })
  })
  
  Data_interact <- reactive({
    input$goButton_interact
    isolate({ 
      search_term_interact <- trimws(input$search_in_interact)
      search_out_interact <- network_table[apply(network_table[, 1:2], 1, FUN = 
                                       function(x, y = search_term_interact) {any(grepl(y, x, ignore.case = TRUE))}), ]
      if(search_term_interact == "Target UNIPROT ID"){return()}
      return(search_out_interact)
    })
  })
  
  Data_interact_probes <- reactive({
    input$goButton_interact
    isolate({ 
      search_term_interact <- trimws(input$search_in_interact)
      search_out_interact <- network_table[apply(network_table[, 1:2], 1, FUN = 
                                                   function(x, y = search_term_interact) {any(grepl(y, x, ignore.case = TRUE))}), ]
      interacting_proteins <- c(search_out_interact$A, search_out_interact$B)
      interacting_proteins_out <- probe_list %>% dplyr::filter(UNIPROTKB %in% interacting_proteins) %>% 
        dplyr::select(CHEMICAL_ID, TARGET_NAME, UNIPROTKB, DATA_SOURCE)
      if(search_term_interact == "Target UNIPROT ID"){return()}
      return(interacting_proteins_out)
    })
  })
  
  
  KEGG <- reactive({
    input$goButton_kegg
    isolate({ 
      filename <- normalizePath(file.path('./KeggPlots', kegg_list[kegg_list$EntryName == input$kegg_in, ]$filename))
      return(list(src = filename,
                  alt = paste("Image number", input$n)))
      
    })
  })
  
  reactome <- reactive({
    input$goButton_reactome
    isolate({ 
      
      return(tryCatch({
        
        viewPathway(trimws(input$reactome_in), organism = "human", readable=TRUE, fixed = TRUE,
                    foldChange=mygeneList, vertex.label.cex = 0.75, col.bin = 1)

      }, error=function(e){}))
      
    })
  })
  
  output$probePlot <- renderTable({
    Data()
  }, caption = "Probes found from basic search:",
  caption.placement = getOption("xtable.caption.placement", "top"))
  
  output$probePlot_interact <- renderTable({
     Data_interact()
   }, caption = "Protein interaction Search. Relevant protein-protein interactions:",
  caption.placement = getOption("xtable.caption.placement", "top"))
   
  output$probePlot_interact_probes <- renderTable({
    Data_interact_probes()
  }, caption = "Probes that may indirectly modulate the target:",
  caption.placement = getOption("xtable.caption.placement", "top"))
  
  output$reactomePlot <- renderPlot({
    reactome()
  })
  
  output$keggPlot <- renderImage({
    KEGG()
  }, deleteFile = FALSE)
  
})