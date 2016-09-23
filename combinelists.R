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
library(ChemmineOB)
smiset <- as(probe_list$canonical_smiles, "SMIset") 
cid(smiset) <- probe_list$molregno
write.SMI(smiset, file="sub.smi", cid=TRUE) 

# i copied molregno from the website into the variable xx
probe_list  <- probe_list %>% filter(!(molregno %in% xx)) 
probe_list_temp <- probe_list
#
probe_list <- probe_list %>%  select(CHEMICAL_ID, pref_name, accession, orthogonal, is_agonist, canonical_smiles, group)
# sgc probes
sgcprobes <- read_csv("accessions.csv")
sgcprobes <- sgcprobes %>% select(CHEMICAL_ID, accession, pref_name) %>% mutate(group = 0)
#make full list
probe_list <- bind_rows(list(sgcprobes, probe_list))
length(unique(probe_list$pref_name))
#make list for protein annotation
protein_accessions <- tibble(accession = unique(probe_list$accession))

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

#add protein information to probe list
probe_list <- left_join(probe_list, protein_information, by = c("accession" = "UNIPROTKB"))
probe_list <- probe_list %>% rename(UNIPROTKB = accession) %>% mutate(source = "CHEMBL") %>% arrange(group, pref_name) %>% rename(target_name = pref_name)
probe_list$source[probe_list$group == 0] <- "chemicalprobes.org"
#remove duplicates - where do they come from?
probe_list <- probe_list[!duplicated(probe_list), ]

# write file
write_excel_csv(probe_list, "PROBELIST.csv", na = "NA")
#write_excel_csv(probe_list, "PROBELIST_autocellcheck.csv", na = "NA")