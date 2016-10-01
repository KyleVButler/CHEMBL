
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(knitr)
source("helpers.R")
shinyServer(function(input, output) {



  
  output$probePlot <- renderTable({
    input$goButton
    search_out <- probe_list[apply(probe_list, 1, FUN = 
                                     function(x, y = isolate(input$search_in)) {any(grepl(y, x, ignore.case = TRUE))}), ]
    
    search_out %>% dplyr::select(CHEMICAL_ID, target_name, UNIPROTKB)
    
  })
  
  output$reactomePlot <- renderPlot({
    observeEvent(input$goButton_reactome,
    isolate(viewPathway(trimws(input$reactome_in), organism = "human", fixed = FALSE,
                foldChange=mygeneList, vertex.label.cex = 0.75)))
  })
  output$keggPlot <- renderImage({
    observeEvent(input$goButton_kegg,
                 isolate(pathview(gene.data = gene_names, pathway.id = trimws(input$kegg_in), 
                                  species = "hsa", col.key = FALSE)))
    list(src = paste(input$kegg_in, ".pathview.png", sep = ""),
         contentType = 'image/png',
         width = 400,
         height = 300,
         alt = "This is alternate text")
    files <- list.files()
    files <- files[grep(input$kegg_in, files)]
    file.remove(files)

  })
})
