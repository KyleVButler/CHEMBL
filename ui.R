
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
source("helpers.R")
shinyUI(fluidPage(

  # Application title
  titlePanel("Chemical Probes", p("")),

  # Sidebar with a slider input for number of bins
  sidebarLayout(position = "left",
    sidebarPanel(strong("Visualization and search functions\n\n"), h1(""),
      textInput("search_in",
                "\nSearch for probes by UNIPROT (e.g. P42345), Gene ontology 
                (e.g. GO:0031931)
                , or target name (e.g. mTOR)", 
                value = "Search Input"),
      actionButton("goButton", "Search"),
      p(""), 
      strong("Visualize Reactome pathway"),
      p("Enter a reactome pathway name, e.g. 'Regulation of TP53 Activity 
                through Methylation', or 'AKT phosphorylates targets in the cytosol' Targets with available probes appear in red. 
                Search term must be punctuated exactly as it appears at:"),
      a("www.reactome.org/PathwayBrowser/", href = "http://www.reactome.org/PathwayBrowser/"),
      textInput("reactome_in", label = "Reactome pathway:", value = ""),
      actionButton("goButton_reactome", "Display"),
      p(""), 
      strong("Visualize KEGG pathway"),
      p("Enter a KEGG pathway id, e.g. 'hsa04014' or 'hsa04068'"),
      textInput("kegg_in", label = "KEGG pathway:", value = ""),
      actionButton("goButton_kegg", "Display")
    ),

    # Show a plot of the generated distribution
    mainPanel(
      tableOutput("probePlot"),
      plotOutput("reactomePlot"),
      imageOutput("keggPlot")
    )
  )
))
