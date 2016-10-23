

library(shiny)
library(shinyIncubator)
library(readr)
kegg_list <- read_csv("kegglist.csv")
probe_list <- read_csv("PROBELIST.csv")
shinyUI(navbarPage(
  "Chemical Probes App",
  tabPanel("Info",
           mainPanel(
             h2("The chemical probes app"), br(),
             img(src='golf.jpeg', align = "right", height = '200'), br(),
             em("Now that you have biological data, how will you choose the best small molecule for pharmacological experiments?"), br(), br(),
             div("This app will help life science researchers choose the best chemicals to modulate pathways and proteins of interest."),br(),br(),
             div("This app contains a comprehensive database of chemical probes. 
                Compounds were classified as probes if they met the criteria defined 
               by the ", a("Structural Genomics Consortium", href = "http://www.thesgc.org"), ("
               The properties of a probe, as defined by the SGC, are as follows:")), br(),
             p("1. In vitro potency of <100 nM at a single protein or two isoforms"),
             p("2. >30-fold selectivity vs other subfamilies"),
             p("3. Demonstration of on-target effect in cells at <1 μM"), br(),
             div("The Chembl database was analyzed to find probes that meet the above criteria. These probes were combined with user-validated probes 
              from ", a("www.chemicalprobes.org", href = "http://www.chemicalprobes.org/"), (" to create the database. The database can be used in the 
                                                                                          following ways, accessible by the navigation bar:")), br(),
             strong("Chemical Probe Search"), br(),
             p("Search for probes by chemical name, protein target, gene, kegg pathway id, reactome pathway id, or by gene ontology id"), br(),
             strong("Search by protein interaction networks"), br(),
             p("Many proteins are not directly targeted by a chemical probe. Input the target of interest into 
               the search bar, and the database will identify probes for any proteins in the local protein-protein interaction network. 
               Interactions were extracted from the Intact database. "),
             a("www.ebi.ac.uk/intact/", href = "http://www.ebi.ac.uk/intact/"), br(),
             strong("KEGG network visualization"), br(),
             p("View which nodes in KEGG networks can be targeted by probes. Useful for generating ideas for combination therapies"), br(),
             p("This website was created with R shiny and uses hadleyverse and bioconductor packages. Kegg graphs were created 
               with the pathview package."), br(),
             p("Contact: kylevbutler@gmail.com")
             )),
  tabPanel(
    "Chemical Probe Search",
    sidebarLayout(
      position = "left",
      sidebarPanel(
        strong("Basic Search"),
        p(
          "Search for probes by target name, UNIPROT, or Gene ontology
          (e.g. GO:0070577 or 'lysine-acetylated histone binding')"
        ),
        textInput("search_in", label = NULL, value = "P42345"),
        #selectInput("search_by", label = "Search By:", choices = 
        #              names(probe_list)[c(1,2,3,8,11,12,13,14,15)]),
        actionButton("goButton", "Search"),br(),br(),
        strong("Display CHEMBL information"),
        textInput("chembl_in", label = "Chembl ID", value = "CHEMBL1801204"),
        actionButton("goButtonChembl", "Display", icon = icon("th"))
        ),
        
      mainPanel(tableOutput("probePlot"), htmlOutput("chembl_out"))
    )
  ),
  tabPanel(
    "Target Network Search",
    sidebarLayout(
      position = "left",
      sidebarPanel(
        strong("Target network search"),
        p(
          "Search for probes by protein target network, using Uniprot IDs. 
          Finds proteins covered by probes in the local protein-protein interaction
          network. Proteins associated with probes will 
          be shown if they are separated from the target by one or two edges, and are  highlighted in red."
        ), 
        a("http://www.uniprot.org/", href = "http://www.uniprot.org/"),
        textInput("search_in_interact", label = NULL, value = "P01116"),
        p("Omit interactions below a given Intact confidence threshold to sparsify the graph."),
        sliderInput("slider1", label = h3("Interaction confidence threshold"), min = 0.1, 
                    max = 0.9, value = 0.36),
        actionButton("goButton_interact", "Search")
        
        ),
      mainPanel(plotOutput("probePlot_interact", width = "100%"), br(), br(), br(), br(), br(), br(), br(), tableOutput("probePlot_interact_table"),
                tableOutput("probePlot_interact_table_2"))
    )
  ),
  tabPanel(
    "KEGG Pathway Visualization",
    sidebarLayout(
      sidebarPanel(
        position = "below",
        # strong("Visualize Reactome pathway"),
        # p("Enter Reactome pathway name exactly as it appears in Reactome database. e.g. 'Regulation of TP53 Activity through Methylation'"),
        # a("reactome.org", href = "http://www.reactome.org/PathwayBrowser/"),
        # textInput("reactome_in", label = NULL, value = "Regulation of TP53 Activity through Methylation"),
        # actionButton("goButton_reactome", "Display"), br(), br(),
        strong("Visualize KEGG pathway"),
        p(
          "Select a KEGG pathway. Targets with available probes
          appear in red. See:"
        ),
        a("www.genome.jp/kegg/pathway.html", href = "http://www.genome.jp/kegg/pathway.html"),
        selectInput("kegg_in", label = NULL, choices = kegg_list$EntryName),
        actionButton("goButton_kegg", "Display")
        ),
      mainPanel(#plotOutput("reactomePlot",  width = "100%"), br(), br(), 
                imageOutput("keggPlot")
                )
      )
    )
  ))