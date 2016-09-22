library(dplyr)
library(tibble)
library(stringr)
library(readr)
chembl_db <- src_sqlite("/Users/kylebutler/Desktop/chembl_21_sqlite/chembl_21.db", create = TRUE)
#type .tables in sqlite3 to see tables
activities <- tbl(chembl_db, "activities")

#collect only compounds with some ic50 below 100 nM
activities <- activities %>% filter(standard_units == "nM" & standard_value < 100) %>% select(molregno)
cmpdstokeep <- collect(activities, n = Inf)
cmpdstokeep <- unique(cmpdstokeep$molregno)
activities <- tbl(chembl_db, "activities")
activities_collected <- activities %>% filter(molregno %in% cmpdstokeep) %>% 
  select(activity_id, assay_id, molregno, standard_value, standard_units, standard_type, data_validity_comment) 
activities_collected <- collect(activities_collected, n = Inf)
length(unique(cmpdstokeep)); length(unique(activities_collected$molregno))
dim(activities_collected)

#remove bad data
activities_collected <- activities_collected[is.na(activities_collected$data_validity_comment), ]
activities_collected <- filter(activities_collected, standard_value > 0)
dim(activities_collected)
#append compound ids and assay ids for checking
molecule_dictionary <- tbl(chembl_db, "molecule_dictionary")
cmpdstokeep <- unique(activities_collected$molregno)
join_vector <- molecule_dictionary %>% select(molregno, chembl_id) %>% filter(molregno %in% cmpdstokeep)
join_vector <- collect(join_vector, n = Inf)
activities_collected <- full_join(activities_collected, join_vector, by = "molregno")
dim(activities_collected)
filter(activities_collected, molregno == 1837801) 
assays <- tbl(chembl_db, "assays")
cmpdstokeep <- unique(activities_collected$assay_id)
join_vector <- assays %>% select(assay_id, assay_type, tid, bao_format, confidence_score, description) %>% filter(assay_id %in% cmpdstokeep)
join_vector <- collect(join_vector, n = Inf)
activities_collected <- full_join(activities_collected, join_vector, by = "assay_id")
dim(activities_collected)
filter(activities_collected, molregno == 1837801) 
target_dictionary <- tbl(chembl_db, "target_dictionary")
cmpdstokeep <- unique(activities_collected$tid)
join_vector <- target_dictionary %>% select(tid, target_type, pref_name, organism) %>%
  filter(tid %in% cmpdstokeep)
join_vector <- collect(join_vector, n = Inf)
activities_collected <- full_join(activities_collected, join_vector, by = "tid")

activities_collected <- activities_collected %>% filter(organism == "Homo sapiens" | is.na(organism)) %>% filter(assay_type %in% c("B", "F"))

#get only compounds with sub 100nM assay for a single protein human targets only
cmpdstokeep <- activities_collected %>% filter(standard_value < 100 & standard_units == "nM" & confidence_score == 9 & 
                                                 organism == "Homo sapiens" & target_type == "SINGLE PROTEIN") %>% select(molregno)
cmpdstokeep <- unique(cmpdstokeep$molregno)
activities_collected <- filter(activities_collected, molregno %in% cmpdstokeep)
#make sure the compound has at least one activity below 1000 nM in a cell based assay
cmpdstokeep <- activities_collected %>% filter(standard_value < 1000 & standard_units == "nM" & bao_format == "BAO_0000219") %>% select(molregno)
cmpdstokeep <- unique(cmpdstokeep$molregno)
activities_collected <- filter(activities_collected, molregno %in% cmpdstokeep)
dim(activities_collected)
length(unique(activities_collected$molregno))
#add component_id and accession
target_components <- tbl(chembl_db, "target_components")
join_vector <- target_components %>% select(tid, component_id) %>% filter(tid %in% activities_collected$tid)
join_vector <- collect(join_vector, n = Inf)
activities_collected <- left_join(activities_collected, join_vector, by = "tid")
#join accession
component_sequences <- tbl(chembl_db, "component_sequences")
join_vector <- component_sequences %>% select(accession, component_id) %>% filter(component_id %in% activities_collected$component_id)
join_vector <- collect(join_vector, n = Inf)
activities_collected <- left_join(activities_collected, join_vector, by = "component_id")

activities_collected <- activities_collected %>% filter(substr(pref_name, 1, 7) != "Cytochr") %>% filter(substr(pref_name, 1, 4) != "HERG")
activities_collected_temp <- activities_collected
dim(activities_collected)
length(unique(activities_collected$molregno))