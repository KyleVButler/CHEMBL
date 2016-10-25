binding 136
receptor activity 136	
signal transducer activity (GO:0004871)	67	
channel regulator activity (GO:0016247)	1	
catalytic activity (GO:0003824)	196	
antioxidant activity (GO:0016209)	1
transporter activity (GO:0005215)	39

panther <- tibble(Function = c("Binding", "Receptor", "Signal Transducer", "Channel Regulator", "Catalytic", "Antioxidant", "Transporter"), 
                  n = c(136, 136, 67, 1, 196, 1, 39))
library(ggplot2)
dev.new()
gg <- ggplot(panther, aes(x = Function, y = n, fill = Function)) + geom_bar(stat="identity") +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), text = element_text(size=20))
gg


#-------------------------    GO ANALYSIS   --------------------------------------

#could find overlap between all gene go membership and probe gene go membership
#get all gene names----------
library(UniProt.ws)
library(clusterProfiler)
library(org.Hs.eg.db)
library(tibble)
library(dplyr)
up <- UniProt.ws(taxId=9606)
egs <- as.character(unique(keys(up, "ENTREZ_GENE")))
gene <- as.character(unique(probe_list$ENTREZ_GENE))
ggo_all <- groupGO(gene     = egs,
               OrgDb    = "org.Hs.eg.db",
               ont      = "MF",
               level = 5,
               readable = TRUE)
length(unique(ggo_all@result$Description))
ggo_all_df <- ggo_all@result %>% select(ID, Description, Count) %>% mutate(Ratio_all = Count / length(egs)) %>% rename(count_all = Count)
ggo_probes <- groupGO(gene     = gene,
                   OrgDb    = "org.Hs.eg.db",
                   ont      = "MF",
                   level = 5,
                   readable = TRUE)
length(unique(ggo_probes@result$Description))
ggo_probes_df <- ggo_probes@result %>% select(ID, Description, Count) %>% mutate(Ratio_probes = Count / length(gene)) %>% rename(count_probes = Count)
#####  may wish to loop through levels above to get many go terms look at MF and BP
ggo_df <- left_join(ggo_all_df, ggo_probes_df) %>% mutate(relative_enrichment = Ratio_probes/Ratio_all) %>% arrange(desc(relative_enrichment))
head(ggo_df %>% filter(count_probes > 3 & count_all > 20))
tail(ggo_df %>% filter(count_probes > 3 & count_all > 20))
head(ggo_df %>% filter(count_all > 50))
tail(ggo_df %>% filter(count_all > 50))
#any correlation between number of probes and number of genes
cor(ggo_df$count_all, ggo_df$count_probes)



barplot(ggo, drop=TRUE, showCategory=20)

kk <- enrichKEGG(gene         = gene,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)
barplot(kk, drop=TRUE, showCategory=15)
head(kk@result)
