

library(shiny)
library(readr)
kegg_list <- read_csv("kegglist.csv")
shinyUI(navbarPage(
  "Chemical Probes App",
  tabPanel("Info",
           mainPanel(
             strong("Information"),
             div("Here you will find a comprehensive database of chemical probes. Compounds were classified as probes if they met the criteria defined 
               by the ", a("Structural Genomics Consortium", href = "http://www.thesgc.org"), ("
               The properties of a probe, as defined by the SGC, are as follows:")), br(),
             p("1. In vitro potency of <100 nM at a single protein or two isoforms"),
             p("2. >30-fold selectivity vs other subfamilies"),
             p("3. Demonstration of on-target effect in cells at <1 μM"), br(),
             div("The Chembl database was analyzed to find probes that meet the above criteria. These probes were combined with user-validated probes 
              from ", a("chemicalprobes.org", href = "http://www.chemicalprobes.org/"), (" to create a database. The database can be used in the 
                                                                                          following ways, accessible by the navigation bar:")), br(),
             strong("Chemical Probe Search"), br(),
             em("Search for probes by protein target or by gene ontology id"), br(), br(),
             em("Search by protein interaction networks"), br(),
             p("Many proteins are not directly targeted by a chemical probe. Input the target of interest into 
               the search bar, and the database will identify probes for any proteins that are known to directly interact with the target of interest. 
               For example, searching for P04637 (TP53 uniprot id) will return BRD7 inhibitors, among others, as BRD7 is known to bind TP53. 
               This function can also be used to find proteins that can be indirectly modulated with a chemical probe."),
             p(""),
             strong("Network Visualization"), br(),
             em("Reactome network visualization"), br(),
             p("View which nodes in reactome networks can be targeted by probes"),
             em("KEGG network visualization"), br(),
             p("View which nodes in KEGG networks can be targeted by probes. Useful for generating ideas for combination therapies"), br(),
             strong("Chembl compounds can be found at:"),
             a("https://www.ebi.ac.uk/chembl/", href = "https://www.ebi.ac.uk/chembl/"), br()
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
        textInput("search_in", label = NULL, value = "Search Input"),
        actionButton("goButton", "Search"), br(), br(),
        strong("Protein interaction search"),
        p(
          "Search for probes by protein interaction network, using Uniprot IDs. Finds proteins that interact with the target, along with probes 
          for those proteins."
        ),
        a("http://www.uniprot.org/", href = "http://www.uniprot.org/"),
        textInput("search_in_interact", label = NULL, value = "Target UNIPROT ID"),
        actionButton("goButton_interact", "Search")
        ),
      mainPanel(tableOutput("probePlot"), tableOutput("probePlot_interact_probes"), tableOutput("probePlot_interact"))
    )
  ),
  tabPanel(
    "Network Visualization",
    sidebarLayout(
      sidebarPanel(
        position = "below",
        strong("Visualize Reactome pathway"),
        p(
          "Enter a formal reactome pathway name, e.g. 'Regulation of TP53 Activity
          through Methylation', or 'AKT phosphorylates targets in the cytosol'. Targets with
          available probes appear in red.
          Search terms must be punctuated exactly as they appear at:"
        ),
        a("www.reactome.org/PathwayBrowser/", href = "http://www.reactome.org/PathwayBrowser/"),
        textInput("reactome_in", label = NULL, value = "Regulation of TP53 Activity through Methylation"),
        actionButton("goButton_reactome", "Display"),
        p(""),
        strong("Visualize KEGG pathway"),
        p(
          "Select a KEGG pathway. Targets with available probes
          appear in red. See:"
        ),
        a("www.genome.jp/kegg/pathway.html", href = "http://www.genome.jp/kegg/pathway.html"),
        selectInput("kegg_in", label = NULL, choices = kegg_list$EntryName),
        actionButton("goButton_kegg", "Display")
        ),
      mainPanel(plotOutput("reactomePlot"), 
                imageOutput("keggPlot")
                )
      )
    )
  ))