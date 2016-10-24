library(dplyr)
library(tibble)
library(stringr)
library(readr)

selectivity_types <- c("Selectivity ratio", "Ratio IC50", "Ratio", "Ratio Ki", "Ratio EC50", "Fold selectivity", 
                              "Selectivity Index", "Selectivity index",
                              "Relative potency", "Ratio pIC50", "Ratio pKi", "Selectivity")

binding_types <- c("Activity", "EC50", "IC50", "Kb", "KB", "Kd", "Ki", "Kinact", "Potency", "Log Ki", 
                          "Log EC50", "pKb")
efficacy_types <- c("Emax", "max activation", "efficacy", "Efficacy")



chembl_db <- src_sqlite("chembl_22.db", create = TRUE)
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
activities_collected <- left_join(activities_collected, join_vector, by = "molregno")
dim(activities_collected)

assays <- tbl(chembl_db, "assays")
cmpdstokeep <- unique(activities_collected$assay_id)
join_vector <- assays %>% select(assay_id, assay_type, tid, bao_format, confidence_score, description) %>% filter(assay_id %in% cmpdstokeep)
join_vector <- collect(join_vector, n = Inf)
activities_collected <- left_join(activities_collected, join_vector, by = "assay_id")
dim(activities_collected)
target_dictionary <- tbl(chembl_db, "target_dictionary")
cmpdstokeep <- unique(activities_collected$tid)
join_vector <- target_dictionary %>% select(tid, target_type, pref_name, organism) %>%
  filter(tid %in% cmpdstokeep)
join_vector <- collect(join_vector, n = Inf)
activities_collected <- left_join(activities_collected, join_vector, by = "tid")

activities_collected <- activities_collected %>% filter(organism == "Homo sapiens" | is.na(organism)) %>% filter(assay_type %in% c("B", "F"))

#get only compounds with sub 100nM assay for a single protein human targets only
cmpdstokeep <- activities_collected %>% filter(standard_value < 100 & standard_units == "nM" & 
                                                 confidence_score %in% c(8,9) & 
                                                 organism == "Homo sapiens" & 
                                                 target_type == "SINGLE PROTEIN") %>% select(molregno)
cmpdstokeep <- unique(cmpdstokeep$molregno)
activities_collected <- filter(activities_collected, molregno %in% cmpdstokeep)
#make sure the compound has at least one activity below 1000 nM in a cell based assay

cmpdstokeep <- activities_collected %>% 
  filter(!grepl('insect|sf9|E Coli|escherichia|bacteria|baculovirus', description, ignore.case = TRUE)) %>% 
    filter(bao_format == "BAO_0000219" & standard_units == "nM" & standard_value < 1000 & confidence_score
           %in% c(4,5,6,7,8,9) & organism == "Homo sapiens") %>% 
  select(molregno) %>% distinct(molregno)

activities_collected <- filter(activities_collected, molregno %in% cmpdstokeep$molregno)
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

#remove some known pains structures and compounds i have found to be wrongly classified like 218012
xx <- c(319012, 1732967, 1839558, 60417, 1334262, 808386, 64904, 1145390, 203388, 
        476314, 832142, 378065, 514511, 59721, 1425178, 40961, 1520447, 1520448, 318799, 1839556, 1520731, 493127, 	
        36624, 357202, 351035, 320166, 218012)

activities_collected <- activities_collected %>% filter(substr(pref_name, 1, 7) != "Cytochr") %>% filter(substr(pref_name, 1, 4) != "HERG") %>% 
  filter(!(molregno %in% xx))
activities_collected_temp <- activities_collected
dim(activities_collected)
length(unique(activities_collected$molregno))