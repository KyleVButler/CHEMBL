#append uniprot information to protein accessions
library(UniProt.ws)
up <- UniProt.ws(taxId=9606)
res2 <- UniProt.ws::select(up,
                          keys = unique(probe_list$accession),
                          columns = c("ENTREZ_GENE", "GO"),
                          keytype = "UNIPROTKB")

View(res2)
res2 <- res2 %>% dplyr::group_by(UNIPROTKB) %>% dplyr::arrange(ENTREZ_GENE) %>% dplyr::slice(1)
res2 <- res2[!duplicated(res2$UNIPROTKB), ]
protein_information <- res2

#KEGG pathway visualizaiton
library("pathview")
pathway <- "hsa04014" 
pathview(gene.data = res2$ENTREZ_GENE, pathway.id = pathway, species = "hsa", col.key = FALSE)
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

sample <- tibble(ENTREZ_GENE = protein_information$ENTREZ_GENE)

sample  <- mutate(sample, Value = 100)
mygeneList <- as.numeric(sample$Value)
names(mygeneList) <- sample$ENTREZ_GENE
dev.new()
viewPathway("Regulation of TP53 Activity through Methylation", organism = "human", readable=TRUE, foldChange=mygeneList, vertex.label.cex = 0.75)
dev.new()
viewPathway("AKT phosphorylates targets in the cytosol", organism = "human", readable=TRUE, foldChange=mygeneList, vertex.label.cex = 0.75)

de <- names(mygeneList)
x <- enrichPathway(gene=de,pvalueCutoff=0.05, readable=T)
enrichMap(x, layout=igraph::layout.kamada.kawai, vertex.label.cex = 0.75, fixed = FALSE)
cnetplot(x, categorySize="pvalue", foldChange=geneList, vertex.label.cex = 0.5)

x <- enrichPathway(gene=sample,pvalueCutoff=0.05, qvalueCutoff=0.05, readable=T)
head(summary(x))
plot(x, showCategory=5)