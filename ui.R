
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
source("helpers.R")
shinyUI(navbarPage("Chemical Probes App",
  tabPanel("Chemical Probe Search",
    sidebarLayout(position = "left",
      sidebarPanel(strong("Search"), p(""),
        p("Search for probes by UNIPROT (e.g. P42345), Gene ontology 
          (e.g. GO:0042393), or target name (e.g. mTOR)"),
        textInput("search_in", label = NULL, value = "Search Input"),
        actionButton("goButton", "Search")
      ),
      mainPanel(
        tableOutput("probePlot")
    )
  )
  ),
  tabPanel("Network Visualization",
           sidebarLayout(
             sidebarPanel(position = "below",
             strong("Visualize Reactome pathway"),
             p("Enter a formal reactome pathway name, e.g. 'Regulation of TP53 Activity 
                through Methylation', or 'AKT phosphorylates targets in the cytosol'. Targets with 
                available probes appear in red. 
                Search terms must be punctuated exactly as they appear at:"),
             a("www.reactome.org/PathwayBrowser/", href = "http://www.reactome.org/PathwayBrowser/"),
             textInput("reactome_in", label = NULL, value = "Reactome pathway"),
             actionButton("goButton_reactome", "Display"),
             p(""), 
             strong("Visualize KEGG pathway"),
             p("Enter a KEGG pathway id, e.g. hsa04014 or hsa05218. Targets with available probes 
                appear in red. See:"),
             a("www.genome.jp/kegg/pathway.html", href = "http://www.genome.jp/kegg/pathway.html"),
             textInput("kegg_in", label = NULL, value = "KEGG pathway id"),
             actionButton("goButton_kegg", "Display (may take up to one minute)")
           ),
           mainPanel(
             plotOutput("reactomePlot"),
             imageOutput("keggPlot")
           )
           )),
  tabPanel("Info",
           verbatimTextOutput("Information here")
           )
))
