library(shiny)
library(readr)
library(shinyIncubator)
library(dplyr)
library(tibble)
#library(ReactomePA)
#library(DOSE)
library(igraph)
#library(png)
#library(grDevices)
probe_list <- read_csv("PROBELIST.csv")
kegg_list <- read_csv("kegglist.csv")
network_table <- read_csv("network_table.csv")
gene_names <- unique(probe_list$ENTREZ_GENE)
gene_names <- as.character(gene_names)
mygeneList <- rep(100, length(gene_names))
names(mygeneList) <- gene_names
probe_uniprot <- as.character(unique(probe_list$UNIPROTKB))

shinyServer(function(input, output) {
  Data <- reactive({
    input$goButton
    isolate({ 
      search_term <- trimws(input$search_in)
      search_out <- probe_list[apply(probe_list, 1, FUN = 
                                       function(x, y = search_term) {any(grepl(y, x, ignore.case = TRUE))}), ]
      search_out <- search_out %>% dplyr::select(CHEMICAL_ID, TARGET_NAME, UNIPROTKB, DATA_SOURCE, is_agonist, orthogonal)
      if(search_term == "Search Input"){return()}
      return(search_out)
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

  


  
  output$probePlot <- renderTable({
    Data()
  }, caption = "Probes found from basic search:",
  caption.placement = getOption("xtable.caption.placement", "top"))
#---------------------------------------------------------------------------  
  output$probePlot_interact <- renderPlot({
    input$goButton_interact
    isolate({
    target <- trimws(input$search_in_interact)
    network_table <- filter(network_table, confidenceScore >= input$slider1)
    probe_uniprot <- as.character(unique(probe_list$UNIPROTKB))
    target_connections <- network_table[network_table$A == target | network_table$B == target, ]
    target_primary <- unique(c(target_connections$A, target_connections$B))
    target_connections <- network_table[network_table$A %in% target_primary | network_table$B %in% target_primary, ]
    relevant_probes <- probe_uniprot[probe_uniprot %in% unique(c(target_connections$A, target_connections$B))]
    probe_primary <- network_table[network_table$A %in% relevant_probes | network_table$B %in% relevant_probes, ]
    probe_primary <- unique(c(probe_primary$A, probe_primary$B))
    middle_nodes <- probe_primary[probe_primary %in% target_primary]
    #build attribute list
    name_list <- tibble(node = c(target_connections$A, target_connections$B), alias = c(target_connections$aliasA, target_connections$aliasB))
    name_list$char <- nchar(name_list$alias)
    name_list <- name_list %>% group_by(node) %>% arrange(char) %>% dplyr::slice(1) %>% data.frame()
    name_list$char[is.na(name_list$char)] <- 0
    name_list$alias[name_list$char > 8] <- name_list$node[name_list$char > 8]
    name_list$alias[name_list$char == 0] <- name_list$node[name_list$char == 0]
    relevant_nodes <- tibble(node = unique(c(middle_nodes, target, relevant_probes)), node_type = "middle")
    relevant_nodes <- left_join(relevant_nodes, name_list)
    relevant_nodes[relevant_nodes$node == target, 2] <- "target"
    relevant_nodes[relevant_nodes$node %in% relevant_probes, 2] <- "probe"
    target_connections <- network_table[network_table$A %in% relevant_nodes$node & network_table$B %in% relevant_nodes$node, ]
    g <- graph.data.frame(as.matrix(target_connections[, c(1, 2)]), vertices = relevant_nodes, directed = FALSE)
    V(g)[V(g)$node_type == "probe"]$color <- "red"
    V(g)[V(g)$node_type == "middle"]$color <- "yellow"
    V(g)[V(g)$node_type == "target"]$color <- "lightblue"
    V(g)$name <- V(g)$alias
    E(g)$weight <- target_connections$confidenceScore * 7
    plot(g, layout = layout.fruchterman.reingold, vertex.size = 4, vertex.label.cex = 1, vertex.label.dist=0.4, edge.curved = FALSE, 
         edge.width=E(g)$weight)
    })
   }, height = 600, width = 900 )
   
  output$probePlot_interact_table <- renderTable({
    input$goButton_interact
    isolate({
    target <- trimws(input$search_in_interact)
    network_table <- filter(network_table, confidenceScore >= input$slider1)
    probe_uniprot <- as.character(unique(probe_list$UNIPROTKB))
    target_connections <- network_table[network_table$A == target | network_table$B == target, ]
    target_primary <- unique(c(target_connections$A, target_connections$B))
    target_connections <- network_table[network_table$A %in% target_primary | network_table$B %in% target_primary, ]
    relevant_probes <- probe_uniprot[probe_uniprot %in% unique(c(target_connections$A, target_connections$B))]
    interacting_proteins_out <- probe_list %>% dplyr::filter(UNIPROTKB %in% relevant_probes) %>% 
      dplyr::select(CHEMICAL_ID, TARGET_NAME, UNIPROTKB, DATA_SOURCE)
    print(interacting_proteins_out)
    })

  }, caption = "Probes in network:",
  caption.placement = getOption("xtable.caption.placement", "top"))

  output$probePlot_interact_table_2 <- renderTable({
    input$goButton_interact
    isolate({
      target <- trimws(input$search_in_interact)
      network_table <- dplyr::filter(network_table, confidenceScore >= input$slider1)
      probe_uniprot <- as.character(unique(probe_list$UNIPROTKB))
      target_connections <- network_table[network_table$A == target | network_table$B == target, ]
      target_primary <- base::unique(c(target_connections$A, target_connections$B))
      target_connections <- network_table[network_table$A %in% target_primary | network_table$B %in% target_primary, ]
      relevant_probes <- probe_uniprot[probe_uniprot %in% base::unique(c(target_connections$A, target_connections$B))]
      probe_primary <- network_table[network_table$A %in% relevant_probes | network_table$B %in% relevant_probes, ]
      probe_primary <- base::unique(c(probe_primary$A, probe_primary$B))
      middle_nodes <- probe_primary[probe_primary %in% target_primary]
      #build attribute list
      name_list <- tibble(node = c(target_connections$A, target_connections$B), alias = c(target_connections$aliasA, target_connections$aliasB))
      name_list$char <- nchar(name_list$alias)
      name_list <- name_list %>% group_by(node) %>% arrange(char) %>% dplyr::slice(1) %>% data.frame()
      name_list$char[is.na(name_list$char)] <- 0
      name_list$alias[name_list$char > 8] <- name_list$node[name_list$char > 8]
      name_list$alias[name_list$char == 0] <- name_list$node[name_list$char == 0]
      relevant_nodes <- tibble(node = base::unique(c(middle_nodes, target, relevant_probes)), node_type = "middle")
      relevant_nodes <- left_join(relevant_nodes, name_list)
      relevant_nodes[relevant_nodes$node == target, 2] <- "target"
      relevant_nodes[relevant_nodes$node %in% relevant_probes, 2] <- "probe"
      target_connections <- network_table[network_table$A %in% relevant_nodes$node & network_table$B %in% relevant_nodes$node, ]
      print(target_connections)
    })
  }, caption = "Relevant Interactions:",
  caption.placement = getOption("xtable.caption.placement", "top"))  
  
  
  output$reactomePlot <- renderPlot({
      input$goButton_reactome
      isolate({ 
        tryCatch({
          print(viewPathway(trimws(input$reactome_in), organism = "human", readable=TRUE, fixed = TRUE,
                            foldChange=mygeneList, vertex.label.cex = 0.75, col.bin = 1))
        }, error=function(e){})
      })
  }, height = 500, width = 500)
  
  output$keggPlot <- renderImage({
    KEGG()
  }, deleteFile = FALSE)
  
  output$chembl_out <- renderUI({
    input$goButtonChembl
    isolate({
      chembl_id <- trimws(input$chembl_in)
      html_in <- '<object data="https://glados-ebitest.rhcloud.com/compound_report_card/CHEMBL25/embed/name_and_classification/" width="800px" height="350px"></object>'
      html_in <- gsub("CHEMBL25", chembl_id, html_in, ignore.case = TRUE)
      return(HTML(html_in))
    #cat('<object data="https://glados-ebitest.rhcloud.com/compound_report_card/', html_in, '/embed/name_and_classification/" width="800px" height="350px"></object>', sep = "")
      })
  })
  
})