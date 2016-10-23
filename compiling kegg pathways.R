keggs <- read_csv("kegg pathways.csv")
keggs$Name <- str_split(keggs$Name, pattern = " - ", simplify = TRUE)[, 1]
keggs$EntryName <- str_c(keggs$Entry, keggs$Name, sep = " ", collapse = NULL)

gene_names <- unique(probe_list$ENTREZ_GENE)
gene_names <- as.character(na.omit(gene_names))
mygeneList <- rep(100, length(gene_names))
names(mygeneList) <- gene_names


library("pathview")

setwd("/Users/kylebutler/Desktop/chembl/KeggPlots")
for(i in 1:nrow(keggs)){
  pathway <- keggs[i, 1]
  pathview(gene.data = gene_names, pathway.id = pathway, species = "hsa", col.key = FALSE, kegg.native = TRUE)
}
files_to_keep <- list.files(pattern = "pathview")
files_to_remove <- list.files()
files_to_remove <- files_to_remove[!(files_to_remove %in% files_to_keep)]
file.remove(files_to_remove)
keggs$filename <- apply(keggs[,1], 1, FUN = function(x) paste(x, ".pathview.png", sep = ""))
setwd("/Users/kylebutler/Desktop/chembl")
write_csv(keggs, "kegglist.csv")
