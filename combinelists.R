probe_list <- bind_rows(list(top_probes_onetarget, top_probes_twotarget))
length(unique(probe_list$pref_name))
length(unique(probe_list$accession))
length(unique(probe_list$molregno))
table(probe_list$is_agonist)
#2765 total human protein single protein targets in initial library
#476 targets are covered (582 probes) if i need to do a manual check of on-target cellular activity - 464 targets after pains filter and 525 potential probes
#382 targets (469 probes) are covered if i use bao to check for on-target cellular activity

#could filter out MW>1000
#PAINS filter
# filter out pains structures - looks like about 70 pains
#library(ChemmineR)
#library(ChemmineOB)
#smiset <- as(probe_list$canonical_smiles, "SMIset") 
#cid(smiset) <- probe_list$molregno
#write.SMI(smiset, file="sub.smi", cid=TRUE) 

# i copied molregno from the website into the variable xx
probe_list  <- probe_list %>% filter(!(molregno %in% xx)) 

# sgc probes
sgcprobes <- read_csv("accessions.csv")
sgcprobes <- sgcprobes %>% select(CHEMICAL_ID, accession, pref_name) %>% mutate(group = 0)
#make full list
probe_list <- bind_rows(list(sgcprobes, probe_list))
probe_list  <- probe_list %>% dplyr::select(CHEMICAL_ID, accession, pref_name, group, canonical_smiles, is_agonist)
length(unique(probe_list$pref_name))
#make list for protein annotation
protein_accessions <- tibble(accession = unique(probe_list$accession))

#add protein information to probe list
probe_list <- left_join(probe_list, protein_information, by = c("accession" = "UNIPROTKB"))

# write file
write_excel_csv(probe_list, "PROBELIST.csv", na = "NA")
#write_excel_csv(probe_list, "PROBELIST_autocellcheck.csv", na = "NA")