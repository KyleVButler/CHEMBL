require(ReactomePA)
require(org.Hs.eg.db)
library(stringr)
library(readr)
gene_names <- unique(probe_list$ENTREZ_GENE)
gene_names <- as.character(gene_names)
mygeneList <- rep(100, length(gene_names))
names(mygeneList) <- gene_names

reactome_pathways <- read_csv("reactomepathways.csv")
# reactome_pathways$ReactomeID <- str_split(reactome_pathways$ReactomeID, " ", n = 2, simplify = TRUE)[, 2]
# reactome_pathways$ReactomeID <- str_split(reactome_pathways$ReactomeID, "Homo", n = 2, simplify = TRUE)[, 1]
# reactome_pathways$ReactomeID <- trimws(reactome_pathways$ReactomeID)
# reactome_pathways_filenames <- read_csv("reactomepathways.csv")
# reactome_pathways_filenames$ReactomeID <- str_split(reactome_pathways_filenames$ReactomeID, " ", n = 2, simplify = TRUE)[, 1]
# reactome_pathways_filenames$ReactomeID <- trimws(reactome_pathways_filenames$ReactomeID)
# reactome_pathways$filename <- reactome_pathways_filenames$ReactomeID
reactome_pathways 

setwd("/Users/kylebutler/Desktop/chembl/ReactomePlots")

for (i in 200:nrow(reactome_pathways)) {
  tryCatch({
    png(filename = paste(reactome_pathways[i, 2], ".png", sep = ""))
    viewPathway(reactome_pathways[i, ]$ReactomeID, organism = "human", readable=TRUE, 
                foldChange=mygeneList, vertex.label.cex = 0.75, col.bin = 1)
    dev.off()
  }, error=function(e){})
}
setwd("/Users/kylebutler/Desktop/chembl")

viewPathway("AKT phosphorylates targets in the cytosol", organism = "human", readable=TRUE, fixed = FALSE, foldChange=mygeneList, vertex.label.cex = 0.75, col.bin = 1)
