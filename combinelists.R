probe_list <- bind_rows(list(top_probes_onetarget, top_probes_twotarget))

length(unique(probe_list$pref_name))
length(unique(probe_list$accession))
length(unique(probe_list$molregno))
table(probe_list$is_agonist)
table(probe_list$orthogonal)

#2765 total human protein single protein targets in initial library
#476 targets are covered (582 probes) if i need to do a manual check of on-target cellular activity - 464 targets after pains filter and 525 potential probes
#382 targets (469 probes) are covered if i use bao to check for on-target cellular activity

#could filter out MW>1000
#PAINS filter
# filter out pains structures - looks like about 70 pains
library(ChemmineR)
smiset <- as(probe_list$canonical_smiles, "SMIset") 
cid(smiset) <- probe_list$molregno
write.SMI(smiset, file="sub.smi", cid=TRUE) 
xx <- c(319012, 1732967, 1839558, 60417, 1334262, 808386, 64904, 1145390, 203388, 
        476314, 832142, 378065, 514511, 59721, 1425178, 40961, 1520447, 1520448, 318799, 1839556, 1520731, 493127, 	
        36624, 357202, 351035, 320166, 1839559, 364904  )
# i copied molregno from the website into the variable xx
detach("package:ChemmineR", unload=TRUE)

probe_list  <- probe_list %>% filter(!(molregno %in% xx)) 
probe_list_temp <- probe_list
#
probe_list <- probe_list %>%  select(CHEMICAL_ID, pref_name, accession, orthogonal, is_agonist, canonical_smiles, group)
# sgc probes
sgcprobes <- read_csv("accessions.csv")
sgcprobes <- sgcprobes %>% select(CHEMICAL_ID, accession, pref_name) %>% mutate(group = 0)
#make full list
probe_list <- bind_rows(list(sgcprobes, probe_list))
probe_list <- rename(probe_list, UNIPROTKB = accession)
#make list for protein annotation

up <- UniProt.ws::UniProt.ws(taxId=9606)
res2 <- UniProt.ws::select(up,
                           keys = unique(probe_list$UNIPROTKB),
                           columns = c("ENTREZ_GENE", "ENTRY-NAME"),
                           keytype = "UNIPROTKB")

res2 <- res2[!duplicated(res2$UNIPROTKB), ]
detach("package:RefNet", unload=TRUE)
detach("package:AnnotationHub", unload=TRUE)
detach("package:UniProt.ws", unload=TRUE)
detach("package:PSICQUIC", unload=TRUE)
detach("package:AnnotationDbi", unload=TRUE)
detach("package:IRanges", unload=TRUE)
detach("package:S4Vectors", unload=TRUE)
detach("package:Biobase", unload=TRUE)
detach("package:BiocGenerics", unload=TRUE)

res2$`ENTRY-NAME` <- stringr::str_split(res2$`ENTRY-NAME`, "_", simplify = TRUE)[,1]

library(mygene)
protein_information <- tibble(ENTREZ_GENE = res2$ENTREZ_GENE, ENTRY_NAME = res2$'ENTRY-NAME', UNIPROTKB = res2$UNIPROTKB, REACTOME = "", 
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
detach("package:mygene", unload=TRUE)
detach("package:S4Vectors", unload=TRUE)
detach("package:Biobase", unload=TRUE)
detach("package:IRanges", unload=TRUE)
detach("package:AnnotationDbi", unload=TRUE)
detach("package:BiocGenerics", unload=TRUE)
#add protein information to probe list

probe_list <- left_join(probe_list, protein_information, by = c("UNIPROTKB" = "UNIPROTKB"))
probe_list <- probe_list  %>% mutate(source = "CHEMBL") %>% arrange(group, target_name)
probe_list$source[probe_list$group == 0] <- "chemicalprobes.org"
#remove duplicates - where do they come from?
probe_list <- probe_list[!duplicated(probe_list), ]

# write file
write_csv(probe_list, "PROBELIST.csv", na = "NA")
