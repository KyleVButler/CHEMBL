#append uniprot information to protein accessions
library(UniProt.ws)
up <- UniProt.ws(taxId=9606)
columns(up)
res2 <- UniProt.ws::select(up,
                          keys = unique(probe_list$UNIPROTKB),
                          columns = c("ENTREZ_GENE", "ENTRY-NAME"),
                          keytype = "UNIPROTKB")

View(res2)
res2 <- res2 %>% dplyr::group_by(UNIPROTKB) %>% dplyr::arrange(ENTREZ_GENE) %>% dplyr::slice(1)
res2 <- res2[!duplicated(res2$UNIPROTKB), ]
protein_information <- res2


###########
library(KEGGgraph)


###################################################
### code chunk number 2: remoteRetrieval (eval = FALSE)
###################################################

tmp <- tempfile()
pName <- "p53 signaling pathway"
pId <- mget(pName, KEGGPATHNAME2ID)[[1]]
retrieveKGML("04014", organism="hsa", destfile=tmp, method="wget", quiet=TRUE)


#KEGG pathway visualizaiton
library("pathview")
pathway <- "hsa04014" 
pv2 <- download.kegg(pathway.id = "04014", species = "hsa",
              file.type= "xml")
tmp <- tempfile()
pv <- pathview(gene.data = sample[, 1], pathway.id = pathway, species = "hsa", col.key = FALSE, kegg.native = FALSE, out.suffix = "x")
keggview.graph(plot.data.gene = pv$plot.data.gene, plot.data.cpd = pv$plot.data.cpd, path.graph = )
pathway <- "hsa04068"  
pathview(gene.data = res2$ENTREZ_GENE, pathway.id = pathway, species = "hsa", col.key = FALSE)

library(mygene)
protein_information <- tibble(ENTREZ_GENE = res2$ENTREZ_GENE, UNIPROTKB = res2$UNIPROTKB, REACTOME = "", 
                              GO_CC = "", GO_BP = "", GO_MF = "", KEGG = "")
for(i in 1:nrow(protein_information)){
  protein_information$REACTOME[i] <- paste(unlist(getGene(protein_information$ENTREZ_GENE[i], fields="all")$pathway$reactome), sep=" ", collapse=" ")
}
for(i in 1:nrow(protein_information)){
  protein_information$GO_CC[i] <- paste(unlist(getGene(protein_information$ENTREZ_GENE[i], fields="all")$go$CC), sep=" ", collapse=" ")
}
for(i in 1:nrow(protein_information)){
  protein_information$GO_BP[i] <- paste(unlist(getGene(protein_information$ENTREZ_GENE[i], fields="all")$go$BP), sep=" ", collapse=" ")
}
for(i in 1:nrow(protein_information)){
  protein_information$GO_MF[i] <- paste(unlist(getGene(protein_information$ENTREZ_GENE[i], fields="all")$go$MF), sep=" ", collapse=" ")
}
for(i in 1:nrow(protein_information)){
  protein_information$KEGG[i] <- paste(unlist(getGene(protein_information$ENTREZ_GENE[i], fields="all")$pathway$kegg), sep=" ", collapse=" ")
}
write_csv(protein_information, "PROTEIN_INFORMATION.csv")



# best to do the panther analysis here:
#http://pantherdb.org/

#ReactomePA
#noticed that for entrez_gene = 5170, mygene returns 77 reactome entries and uniprot.ws returns 17
require(ReactomePA)
require(org.Hs.eg.db)
library(tibble)
gene_names <- unique(probe_list$ENTREZ_GENE)
gene_names <- as.character(na.omit(gene_names))
mygeneList <- rep(100, length(gene_names))
names(mygeneList) <- gene_names
dev.new()
viewPathway(trimws(input$reactome_in), organism = "human", fixed = FALSE, readable = TRUE,
            foldChange=mygeneList, vertex.label.cex = 0.75)
png(filename = "reactomeimg.png")
viewPathway("AKT phosphorylates targets in the cytosol", organism = "human", readable=TRUE, 
                                foldChange=mygeneList, vertex.label.cex = 0.75)
dev.off()
de <- names(mygeneList)
x <- enrichPathway(gene=de,pvalueCutoff=0.05, readable=T)
enrichMap(x, layout=igraph::layout.kamada.kawai, vertex.label.cex = 0.75, fixed = FALSE)
cnetplot(x, categorySize="pvalue", foldChange=geneList, vertex.label.cex = 0.5)

x <- enrichPathway(gene=sample,pvalueCutoff=0.05, qvalueCutoff=0.05, readable=T)
head(summary(x))
plot(x, showCategory=5)



get_interacting_probes <- function(gene_id, probe_list, database_source){
  library(RefNet)
  refnet <- RefNet()
  probe_gene_ids <- as.character(unique(probe_list$UNIPROTKB))
  network_table <- interactions(refnet, species="9606", id=gene_id, provider= database_source)
  network_table$A <- str_split(network_table$A, ":", simplify = TRUE)[, 2]
  network_table$B <- str_split(network_table$B, ":", simplify = TRUE)[, 2]
  writeLines("Interactions found:")
  print(network_table %>% dplyr::filter(A %in% probe_gene_ids | B %in% probe_gene_ids) %>% dplyr::select(A, B, detectionMethod, publicationID) %>%
    dplyr::rename(Protein1 = A, Protein2 = B))
  writeLines("\n Probes found:")
  print(probe_list %>% dplyr::select(CHEMICAL_ID, UNIPROTKB, target_name, orthogonal, is_agonist) %>% dplyr::filter(UNIPROTKB %in% c(network_table$A, network_table$B)))
}

get_interacting_probes("Q01196", probe_list, "mentha")
get_interacting_probes("O43524", probe_list, "InnateDB-All")  #not many hits from this one
get_interacting_probes("O43524", probe_list, "IntAct")
