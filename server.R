
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(knitr)
source("helpers.R")
rand_num <- round(runif(1, 0, 9), 0)
if(rand_num == 7){
  list_of_files <- list.files(path = ".", pattern = ".png", full.names = TRUE)
  do.call(file.remove, list(list_of_files))
  list_of_files <- list.files(path = ".", pattern = ".xml", full.names = TRUE)
  do.call(file.remove, list(list_of_files))
}
shinyServer(function(input, output) {
  Data <- reactive({
    input$goButton
    isolate({ 
      search_term <- trimws(input$search_in)
      search_out <- probe_list[apply(probe_list, 1, FUN = 
                                       function(x, y = search_term) {any(grepl(y, x, ignore.case = TRUE))}), ]
      search_out <- search_out %>% dplyr::select(CHEMICAL_ID, target_name, UNIPROTKB)
      if(search_term == "Search Input"){return()}
      return(search_out)
    })
  })
  KEGG <- reactive({
    input$goButton_kegg
    isolate({ 
      if(!(input$kegg_in %in% c("","KEGG pathway id"))){
      kegg_pathway <- trimws(input$kegg_in)
      pathview(gene.data = gene_names, pathway.id = kegg_pathway, 
               species = "hsa", col.key = FALSE)
      return(list(src = paste(kegg_pathway, ".pathview.png", sep = ""),
         contentType = 'image/png',
         width = 1200,
         height = 750,
         alt = "KEGG Pathway Visualization"))
      }
    ####need a way to delete files afterwards, perhaps by putting a random number suffix on 
      ### file output then remove.files
      list_of_files <- list.files(path = ".", pattern = ".png", full.names = TRUE)
      do.call(file.remove, list(list_of_files))
    })
  })

  
  output$probePlot <- renderTable({
    Data()
  })
  
  output$reactomePlot <- renderPlot({
    observeEvent(input$goButton_reactome,
    isolate(viewPathway(trimws(input$reactome_in), organism = "human", fixed = FALSE, readable = TRUE,
                foldChange=mygeneList, vertex.label.cex = 0.75)))
  })
  
  output$keggPlot <- renderImage({
    KEGG()
  }, deleteFile = TRUE)
  
})
