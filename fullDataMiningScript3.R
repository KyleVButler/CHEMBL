library(dplyr)
library(stringr)
library(readr)

# group the different types of activity measurements, for use later in the script
selectivity_types <- c("Selectivity ratio", "Ratio IC50", "Ratio", "Ratio Ki", "Ratio EC50", "Fold selectivity", 
                              "Selectivity Index", "Selectivity index", 
                              "Relative potency", "Ratio pIC50", "Ratio pKi", "Selectivity")
thermal_melting_types <- c("Thermal melting change", "Delta Tm")
binding_types <- c("Activity", "EC50", "IC50", "Kb", "KB", "Kd", "Ki", "Kinact", "Potency", "Log Ki", 
                          "Log EC50", "pKb")
single_point_types <- c("Inhibition", "Potency", "Activity", "Residual activity", "Residual Activity", "Delta Tm", "Thermal melting change")
efficacy_types <- c("Emax", "max activation", "efficacy", "Efficacy")

#Load chembl into R
#load chembl_22.db
chembl_db <- src_sqlite("chembl_22.db", create = TRUE)
activities <- tbl(chembl_db, "activities")

#keep only compounds that have at least one ic50 <= 100 nM and collect data into R variable:
#activities_collected
activities <- activities %>% filter(standard_units == "nM" & standard_value <= 100) %>% select(molregno)
cmpdstokeep <- collect(activities, n = Inf)
cmpdstokeep <- unique(cmpdstokeep$molregno)
activities <- tbl(chembl_db, "activities")
activities_collected <- activities %>% filter(molregno %in% cmpdstokeep) %>% 
  select(activity_id, assay_id, molregno, standard_value, standard_units, standard_type, data_validity_comment) 
activities_collected <- collect(activities_collected, n = Inf)

#remove bad data flagged in ChEMBL and data where the value is zero
activities_collected <- activities_collected[is.na(activities_collected$data_validity_comment), ]
activities_collected <- filter(activities_collected, standard_value > 0)

#append compound information and assay information to dataset
molecule_dictionary <- tbl(chembl_db, "molecule_dictionary")
cmpdstokeep <- unique(activities_collected$molregno)
join_vector <- molecule_dictionary %>% select(molregno, chembl_id) %>% filter(molregno %in% cmpdstokeep)
join_vector <- collect(join_vector, n = Inf)
activities_collected <- left_join(activities_collected, join_vector, by = "molregno")
assays <- tbl(chembl_db, "assays")
cmpdstokeep <- unique(activities_collected$assay_id)
join_vector <- assays %>% select(assay_id, assay_type, tid, bao_format, confidence_score, description) %>% filter(assay_id %in% cmpdstokeep)
join_vector <- collect(join_vector, n = Inf)
activities_collected <- left_join(activities_collected, join_vector, by = "assay_id")

#Also add target information to dataset
target_dictionary <- tbl(chembl_db, "target_dictionary")
cmpdstokeep <- unique(activities_collected$tid)
join_vector <- target_dictionary %>% select(tid, target_type, pref_name, organism) %>%
  filter(tid %in% cmpdstokeep)
join_vector <- collect(join_vector, n = Inf)
activities_collected <- left_join(activities_collected, join_vector, by = "tid")

#exclude admet data and concentrate analysis on human (or unclassified) proteins
activities_collected <- activities_collected %>% filter(organism == "Homo sapiens" | is.na(organism)) %>% filter(assay_type %in% c("B", "F"))

#keep compounds that have ic50 <=100nM in a high confidence human protein assay
cmpdstokeep <- activities_collected %>% filter(standard_value <= 100 & standard_units == "nM" & 
                                                 confidence_score == 9 & 
                                                 organism == "Homo sapiens" & 
                                                 target_type == "SINGLE PROTEIN") %>% select(molregno)
activities_collected <- filter(activities_collected, molregno %in% cmpdstokeep$molregno)

#keep compounds that have ic50 <1000nM in a cell based assay associated with a human protein
#cell based assays have bao_format == "BAO_0000219" 
#to make sure the assays are in mammalian cells, we exclude assays that have "insect", etc. in the description
cmpdstokeep <- activities_collected %>% 
  filter(!grepl('insect|sf9|E Coli|escherichia|bacteria|baculovirus', description, ignore.case = TRUE)) %>% 
    filter(bao_format == "BAO_0000219" & standard_units == "nM" & standard_value < 1000 & confidence_score
           %in% c(4,5,6,7,8,9) & organism == "Homo sapiens") %>% 
  select(molregno) %>% distinct(molregno)
activities_collected <- filter(activities_collected, molregno %in% cmpdstokeep$molregno)

#add component_id and accession to the data set
target_components <- tbl(chembl_db, "target_components")
join_vector <- target_components %>% select(tid, component_id) %>% filter(tid %in% activities_collected$tid)
join_vector <- collect(join_vector, n = Inf)
activities_collected <- left_join(activities_collected, join_vector, by = "tid")
component_sequences <- tbl(chembl_db, "component_sequences")
join_vector <- component_sequences %>% select(accession, component_id) %>% filter(component_id %in% activities_collected$component_id)
join_vector <- collect(join_vector, n = Inf)
activities_collected <- left_join(activities_collected, join_vector, by = "component_id")

#remove some known pains structures (listed in text) and a few compounds 
#i have found to be wrongly classified like 218012
xx <- c(218012, 33470, 684443, 1363818, 1839559, 
        558150, 341114, 495196, 1839557, 1145256, 1557952, 1430343, 49110, 464846, 1839560, 1732967, 	
        1839558, 808386, 1145390, 514511, 1839556, 202347, 1425178, 703577  )

activities_collected <- activities_collected %>% filter(substr(pref_name, 1, 7) != "Cytochr") %>% filter(substr(pref_name, 1, 4) != "HERG") %>% 
  filter(!(molregno %in% xx)) 

#display the most common types of bioactivity observations
names(table(activities_collected$standard_type))[as.vector(table(activities_collected$standard_type) > 100)]
table(activities_collected$standard_type)[as.vector(table(activities_collected$standard_type) > 100)]

#What about Log Ki, Log EC50, pKb - need to change to nM for standard_units == NA -- 
#chembl should have done this, use code below if not
#activities_collected[activities_collected$standard_type == "Log Ki" & is.na(activities_collected$standard_units), ]$standard_value <-
#  10^(-(activities_collected[activities_collected$standard_type == "Log Ki" & is.na(activities_collected$standard_units), ]$standard_value)) * 10^9
#activities_collected[activities_collected$standard_type == "Log Ki" & is.na(activities_collected$standard_units), ]$standard_units <- "nM"

#collect potential probes in compound_summary and annotate compounds with main target, 
#number of targets assayed against, and total number of observations
compound_summary <- activities_collected %>% count(molregno) %>% dplyr::rename(n_total = n)
selectivity_measurements <- activities_collected %>% group_by(molregno) %>% 
  filter(standard_type %in% selectivity_types & assay_type == "B") %>% 
  count(molregno) %>% dplyr::rename(n_selectivity = n)
potency_measurements <- activities_collected %>% group_by(molregno) %>% 
  filter(target_type == "SINGLE PROTEIN" & organism == "Homo sapiens" & assay_type == "B" & confidence_score > 7 & 
           standard_units %in% c("nM", "%", "degrees C") & 
           standard_type %in% c(binding_types, single_point_types)) %>% 
  distinct(pref_name) %>% count(molregno) %>% dplyr::rename(n_potency = n)
compound_summary <- left_join(compound_summary, selectivity_measurements)
compound_summary <- left_join(compound_summary, potency_measurements)
compound_summary[is.na(compound_summary)] <- 0

#keep only compounds with n_total > 2 and n_select > 1
compound_summary <- compound_summary %>% mutate(n_select = n_selectivity + n_potency) %>% select(molregno, n_total, n_select) %>% 
  filter(n_total > 2 & n_select > 1)

# make activities list with selected compounds 
activities_collected2 <- activities_collected %>% filter(molregno %in% compound_summary$molregno) 

#match each compound with its most potent target and keep only those that have one target under 100 nM
min_target <- activities_collected2 %>% 
  filter(target_type == "SINGLE PROTEIN" & organism == "Homo sapiens" & assay_type == "B" & confidence_score == 9 & 
           standard_units == "nM" & standard_type %in% binding_types & standard_value <= 100) %>%
  group_by(molregno, pref_name, accession) %>% summarize(min_value = min(standard_value)) %>% group_by(molregno) %>% 
  filter(n() == 1) 
compound_summary <- compound_summary %>% filter(molregno %in% min_target$molregno)
compound_summary <- right_join(compound_summary, min_target)
compound_summary <- compound_summary %>% rename(target_pref_name = pref_name, target_accession = accession)
activities_collected2 <- filter(activities_collected2, molregno %in% compound_summary$molregno)
activities_collected2 <- left_join(activities_collected2, compound_summary)

#check that a cell assay exists for the main target
cmpdstokeep <- activities_collected2 %>% filter(!grepl('insect|sf9|E Coli|escherichia|bacteria|baculovirus', description, 
                                                                             ignore.case = TRUE)) %>% 
  filter(bao_format == "BAO_0000219" & standard_units == "nM" & standard_value < 1000 & organism == "Homo sapiens") %>% 
  filter(target_accession == accession)
compound_summary <- compound_summary %>% filter(molregno %in% unique(cmpdstokeep$molregno))

#find if probes are agonist or not by checking if efficacy at main target is >25%
agonist_check <- activities_collected2 %>% filter(accession == target_accession & standard_units == "%" & standard_type %in% 
                                                    efficacy_types & standard_value > 25) 
compound_summary$is_agonist <- "NO"
compound_summary$is_agonist[compound_summary$molregno %in% agonist_check$molregno] <- "YES"

#check that no selectivity measurements containing target_accession are under 30 
cmpdstokeep <- activities_collected2 %>% 
  filter(standard_type %in% selectivity_types) %>% 
  filter(assay_type == "B" & organism == "Homo sapiens" & standard_value > 1 &
           standard_value < 30 & accession == target_accession) 
compound_summary <- compound_summary %>% filter(!(molregno %in% unique(cmpdstokeep$molregno)))

#at this point we only focus on selectivity, so we will remove any assays that contain the target_accession from
#compound_summary from activities_collected2

activities_collected2 <- activities_collected2 %>% filter(molregno %in% compound_summary$molregno & confidence_score %in% c(4,5,6,7,8,9)) %>%
  filter(!(is.na(accession)))

for(i in 1:nrow(compound_summary)){
  exempt_assays <- activities_collected2 %>% filter(molregno == compound_summary$molregno[i]) %>% 
    filter(accession == target_accession)
  activities_collected2 <- activities_collected2 %>% filter(!(molregno == compound_summary$molregno[i] & 
                                                                pref_name %in% exempt_assays$pref_name))
}


#check that ic50 type measurements have at least 30 fold selectivity
cmpdstokeep <- activities_collected2 %>% mutate(fold_selectivity = (standard_value/min_value)) %>% filter(
  organism == "Homo sapiens" & standard_units == "nM" & standard_type %in% binding_types & fold_selectivity < 30)
compound_summary <- compound_summary %>% filter(!(molregno %in% unique(cmpdstokeep$molregno)))

#check no thermal shift are over 5 for any other proteins
cmpdstokeep <- activities_collected2 %>% filter(molregno %in% compound_summary$molregno) %>%
  filter(target_type == "SINGLE PROTEIN" & organism == "Homo sapiens" & confidence_score %in% c(8, 9) & standard_units == "degrees C" & 
           standard_type %in% thermal_melting_types & standard_value > 5)
compound_summary <- compound_summary %>% filter(!(molregno %in% unique(cmpdstokeep$molregno)))

#check that no % activity type measurements are above 50 for targets without known ic50 values
cmpdstokeep <- activities_collected2 %>% filter(standard_type %in% c("Activity", "Residual activity", "Residual Activity") & 
                                                   standard_units == "%" & assay_type == "B" & target_type == "SINGLE PROTEIN" &
                                                   organism == "Homo sapiens" & standard_value < 50)
compound_summary <- compound_summary %>% filter(!(molregno %in% unique(cmpdstokeep$molregno)))

#check that no % inhibition type measurements are above 50 for targets without known ic50 values
cmpdstokeep <- activities_collected2 %>%  filter(standard_type == "Inhibition" & 
                                                   standard_units == "%" & assay_type == "B" & target_type == "SINGLE PROTEIN" &
                                                   organism == "Homo sapiens" & standard_value > 50)
compound_summary <- compound_summary %>% filter(!(molregno %in% unique(cmpdstokeep$molregno)))

# add chembl ids and smiles and inchi
molecule_dictionary <- tbl(chembl_db, "molecule_dictionary")
join_vector <- molecule_dictionary %>% select(molregno, chembl_id) %>% filter(molregno %in% compound_summary$molregno)
join_vector <- collect(join_vector, n = Inf)
compound_summary <- left_join(compound_summary, join_vector, by = "molregno")
compound_structures <- tbl(chembl_db, "compound_structures")
join_vector <- compound_structures %>% select(molregno, canonical_smiles, standard_inchi, 
                                              standard_inchi_key) %>% filter(molregno %in% compound_summary$molregno)
join_vector <- collect(join_vector, n = Inf)
compound_summary <- left_join(compound_summary, join_vector, by = "molregno")

#remove duplicated compounds
compound_summary <- compound_summary %>% group_by(molregno) %>% slice(1)

#remove those with no smiles
compound_summary <- compound_summary[!(is.na(compound_summary$canonical_smiles)), ]

#rank probes based on number of observations, and for each target, find tanimoto distance from best probe 

library(rcdk)


fun_get_tanimoto <- function(smilesin){
  query.mol <- parse.smiles(smilesin)[[1]]
  target.mols <- parse.smiles(smilesin)
  query.fp <- get.fingerprint(query.mol, type='maccs')
  target.fps <- lapply(target.mols, get.fingerprint, type='maccs')
  sims <- unname(unlist(lapply(target.fps, distance, fp2=query.fp, method='tanimoto')))
  return(sims)
}
compound_summary <- rename(compound_summary, pref_name = target_pref_name)
compound_summary <- compound_summary %>% group_by(is_agonist, pref_name) %>% 
  arrange(is_agonist, pref_name, desc(n_select), desc(n_total), min_value) %>% 
  mutate(tanimoto = fun_get_tanimoto(canonical_smiles))
detach("package:rcdk", unload=TRUE)
detach("package:fingerprint", unload=TRUE)
compound_summary$orthogonal <- "PRIMARY"
compound_summary$orthogonal[compound_summary$tanimoto < 0.8] <- "ORTHOGONAL"


g <- compound_summary %>% ggplot(aes(tanimoto), colors = blue, fill = blue) + 
  geom_histogram(breaks=c(seq(0, 1, by=0.05)), color = "blue", fill = "light blue") +
  ggtitle("Tanimoto coefficients for potential probes")
ggsave("tanimotoplot.png", plot = g, dpi = 300, height = 4, width = 4)

#make main list, containing the top primary and orthogonal probes for each target, for both agonists and non-agonists
top_probes_onetarget <- compound_summary %>% group_by(pref_name, is_agonist, orthogonal) %>% 
  arrange(desc(n_total), desc(n_select)) %>% slice(1)
top_probes_onetarget$accession <- top_probes_onetarget$target_accession

#designate as group 1 for one target
top_probes_onetarget <- top_probes_onetarget %>% mutate(group = 1) %>% rename(CHEMICAL_ID = chembl_id)


##############


###############

#now select for probes that have two homologous targets under 100 nM potency
#want to annotate compounds with main target, number of targets assayed against (n_select), and total number of observations
compound_summary <- activities_collected %>% count(molregno) %>% dplyr::rename(n_total = n)
selectivity_measurements <- activities_collected %>% group_by(molregno) %>% 
  filter(standard_type %in% selectivity_types & assay_type == "B" & confidence_score >3) %>% 
  count(molregno) %>% dplyr::rename(n_selectivity = n)
potency_measurements <- activities_collected %>% group_by(molregno) %>% 
  filter(target_type == "SINGLE PROTEIN" & organism == "Homo sapiens" & assay_type == "B" & confidence_score > 7 & 
           standard_units %in% c("nM", "%", "degrees C") & 
           standard_type %in% c(binding_types, single_point_types)) %>% 
  distinct(pref_name) %>% count(molregno) %>% dplyr::rename(n_potency = n)
compound_summary <- left_join(compound_summary, selectivity_measurements)
compound_summary <- left_join(compound_summary, potency_measurements)
compound_summary[is.na(compound_summary)] <- 0

#keep only compounds with n_total > 4 and n_potency > 5 - need to have activities at two targets, and two off-targets
compound_summary <- compound_summary %>% mutate(n_select = n_selectivity + n_potency) %>% select(molregno, n_total, n_select, n_potency) %>% 
  filter(n_total > 4 & n_potency > 3)

# make activities list with selected compounds
activities_collected2 <- activities_collected %>% filter(molregno %in% compound_summary$molregno) 

#match each compound with its most potent target and keep only those that have two targets under 100 nM
min_target <- activities_collected2 %>% 
  filter(target_type == "SINGLE PROTEIN" & organism == "Homo sapiens" & confidence_score == 9 & assay_type == "B" & 
           standard_units == "nM" & standard_type %in% binding_types & standard_value <= 100) %>%
  group_by(molregno, pref_name, accession) %>% summarize(min_value = min(standard_value)) %>% group_by(molregno) %>% 
  filter(n() == 2) 
compound_summary <- compound_summary %>% filter(molregno %in% min_target$molregno)
compound_summary <- left_join(compound_summary, min_target)
compound_summary <- compound_summary %>% filter(!(is.na(pref_name))) %>% rename(target_pref_name = pref_name, target_accession = accession)
activities_collected2 <- filter(activities_collected2, molregno %in% compound_summary$molregno)
activities_collected2 <- left_join(activities_collected2, compound_summary)
activities_collected2 <- activities_collected2 %>% filter(!(is.na(target_pref_name)))

#check that a cell assay exists for one of the main targets
cmpdstokeep <- activities_collected2 %>% filter(!grepl('insect|sf9|E Coli|escherichia|bacteria|baculovirus', description, 
                                                       ignore.case = TRUE)) %>% 
  filter(bao_format == "BAO_0000219" & standard_units == "nM" & standard_value < 1000 & organism == "Homo sapiens") %>% 
  filter(target_accession == accession)
compound_summary <- compound_summary %>% filter(molregno %in% unique(cmpdstokeep$molregno))

#check to see if the names are similar for the targets - this is a good shortcut to confirm they are homologous
namecheck1 <- compound_summary %>% group_by(molregno) %>% arrange(molregno, target_pref_name) %>% select(molregno, target_pref_name) %>% 
  slice(1) %>% rename(pref_one = target_pref_name)
namecheck2 <- compound_summary %>% group_by(molregno) %>% arrange(molregno, target_pref_name) %>% select(molregno, target_pref_name) %>% 
  slice(2) %>% rename(pref_two = target_pref_name)
namecheck <- left_join(namecheck1, namecheck2)
namecheck <- namecheck %>% mutate(check1 = (substr(pref_one, 1, 6)) == substr(pref_two, 1, 6)) %>% 
  mutate(check2 = (str_sub(pref_one, -6, -1)) == str_sub(pref_two, -6, -1))
namecheck <- namecheck %>% filter(check1 == TRUE | check2 == TRUE)
compound_summary <- compound_summary %>% filter(molregno %in% unique(namecheck$molregno))

#check that no selectivity measurements with one main target and one off target are > 30
cmpdstokeep <- unique(compound_summary$molregno)
cmpdstoremove <- rep(TRUE, length(cmpdstokeep))
for(i in 1:length(cmpdstokeep)){
  assay_subset <- activities_collected2 %>% filter(molregno == cmpdstokeep[i] & standard_type %in% selectivity_types & standard_value < 30 & 
                                                     standard_value >= 1) %>%
    filter(!(is.na(accession))) 
  assays_to_keep <- assay_subset %>% filter(accession == target_accession) 
  assay_subset <- assay_subset %>% filter(description %in% assays_to_keep$description)
  assay_subset <- assay_subset %>% filter(!(accession %in% unique(assays_to_keep$target_accession)))
  if(nrow(assay_subset) > 0){
    cmpdstoremove[i] <- FALSE
  }
}

cmpdstokeep <- cmpdstokeep[cmpdstoremove]
compound_summary <- compound_summary %>% filter(molregno %in% cmpdstokeep)

#at this point we only focus on selectivity, so we will remove any assays that contain the target_accession from
#compound_summary from activities_collected2

activities_collected2 <- activities_collected2 %>% filter(molregno %in% compound_summary$molregno & confidence_score %in% c(4,5,6,7,8,9)) %>%
  filter(!(is.na(accession)))

for(i in 1:nrow(compound_summary)){
  exempt_assays <- activities_collected2 %>% filter(molregno == compound_summary$molregno[i]) %>% 
    filter(accession == target_accession)
  activities_collected2 <- activities_collected2 %>% filter(!(molregno == compound_summary$molregno[i] & 
                                                                pref_name %in% exempt_assays$pref_name))
}

#check that ic50 type measurements have at least 30 fold selectivity
cmpdstokeep <- activities_collected2 %>% mutate(fold_selectivity = (standard_value/min_value)) %>% filter(
  organism == "Homo sapiens" & standard_units == "nM" & standard_type %in% binding_types & fold_selectivity < 30)
compound_summary <- compound_summary %>% filter(!(molregno %in% unique(cmpdstokeep$molregno)))


#check no thermal shift are over 5 for any other proteins
cmpdstokeep <- activities_collected2 %>% filter(molregno %in% compound_summary$molregno) %>%
  filter(target_type == "SINGLE PROTEIN" & organism == "Homo sapiens" & confidence_score %in% c(8, 9) & standard_units == "degrees C" & 
           standard_type %in% thermal_melting_types & standard_value > 5)
compound_summary <- compound_summary %>% filter(!(molregno %in% unique(cmpdstokeep$molregno)))

#check that no % activity type measurements are above 50 for targets without known ic50 values
cmpdstokeep <- activities_collected2 %>% filter(standard_type %in% c("Activity", "Residual activity", "Residual Activity") & 
                                                  standard_units == "%" & assay_type == "B" & target_type == "SINGLE PROTEIN" &
                                                  organism == "Homo sapiens" & standard_value < 50)
compound_summary <- compound_summary %>% filter(!(molregno %in% unique(cmpdstokeep$molregno)))

#check that no % inhibition type measurements are above 50 for targets without known ic50 values
cmpdstokeep <- activities_collected2 %>%  filter(standard_type == "Inhibition" & 
                                                   standard_units == "%" & assay_type == "B" & target_type == "SINGLE PROTEIN" &
                                                   organism == "Homo sapiens" & standard_value > 50)
compound_summary <- compound_summary %>% filter(!(molregno %in% unique(cmpdstokeep$molregno)))

#change the names back
compound_summary$pref_name <- compound_summary$target_pref_name
compound_summary$accession <- compound_summary$target_accession

#remove compounds where both targets are in one target list
cmpdstokeep <- compound_summary %>% filter(!(pref_name %in% top_probes_onetarget$pref_name))
compound_summary <- compound_summary %>% filter(molregno %in% cmpdstokeep$molregno)

#select probes with highest number of observations
cmpdstokeep <- compound_summary %>% group_by(pref_name) %>% 
  arrange(desc(n_total), desc(n_select)) %>% slice(1:2)
top_probes_twotarget <- compound_summary %>% filter(molregno %in% cmpdstokeep$molregno) 

# add chembl ids and smiles and inchi
molecule_dictionary <- tbl(chembl_db, "molecule_dictionary")
join_vector <- molecule_dictionary %>% select(molregno, chembl_id) %>% filter(molregno %in% top_probes_twotarget$molregno)
join_vector <- collect(join_vector, n = Inf)
top_probes_twotarget <- left_join(top_probes_twotarget, join_vector, by = "molregno")
compound_structures <- tbl(chembl_db, "compound_structures")
join_vector <- compound_structures %>% select(molregno, canonical_smiles, standard_inchi, standard_inchi_key) %>% 
  filter(molregno %in% top_probes_twotarget$molregno)
join_vector <- collect(join_vector, n = Inf)
top_probes_twotarget <- left_join(top_probes_twotarget, join_vector, by = "molregno")


#the two target probes are group 2
top_probes_twotarget <-  top_probes_twotarget %>% mutate(group = 2) %>% rename(CHEMICAL_ID = chembl_id)

#join the two data sets
probe_list <- bind_rows(list(top_probes_onetarget, top_probes_twotarget))

#display summary numbers
length(unique(probe_list$pref_name))
length(unique(probe_list$accession))
length(unique(probe_list$molregno))
table(probe_list$is_agonist)
table(probe_list$orthogonal)
median(probe_list$n_total)
median(probe_list$n_select)

#PAINS filter
# filter out pains structures 
require("ChemmineR")

smiset <- as(probe_list$canonical_smiles, "SMIset") 
cid(smiset) <- probe_list$molregno
write.SMI(smiset, file="sub.smi", cid=TRUE) 
#the sub.smi file should be submitted to the website: cbligand.org/PAINS/search_struct.php
#the ids of compounds with substructures that are forbidden (catechols, etc) are manually copied
#into the variable xx below
xx <- c(1839564, 1334262 )
detach("package:ChemmineR", unload=TRUE)
probe_list  <- probe_list %>% filter(!(molregno %in% xx)) 
probe_list <- probe_list %>%  select(CHEMICAL_ID, pref_name, accession, orthogonal, is_agonist, canonical_smiles, group)

####---------------- histograms of selectivity
hist_plots <- bind_rows(top_probes_onetarget %>% dplyr::select(molregno, n_total, n_select), top_probes_twotarget %>% dplyr::select(molregno, n_total, n_select))
hist_plots <- hist_plots %>% ungroup() %>% dplyr::select(molregno, n_total, n_select) %>% 
  dplyr::rename(`Selectivity measurements` = n_select, `Total activities` = n_total)
hist_plots <- hist_plots[!(duplicated(hist_plots)),]
hist_plots$`Total activities`[hist_plots$`Total activities`>100] <- 100
hist_plots$`Selectivity measurements`[hist_plots$`Selectivity measurements`>50] <- 50
g <- ggplot(data = hist_plots, aes(`Total activities`), colors = "blue") + 
  #stat_count(color = "blue", fill = "light blue") +
  geom_histogram(breaks=c(seq(0, 100, by=5)), color = "blue", fill = "light blue") +
  scale_x_continuous(limits=c(0, 100), breaks=c(seq(0, 100, by=10)), labels=c(seq(0,90, by=10), "100+")) + 
  annotate("text", x = 65, y = 100, label = "Median = 11", size = 6) + ggtitle("Total measurements") + xlab("") +
  theme(text = element_text(size = 16))
ggsave("totalactivities.png", plot = g, dpi = 300, height = 4, width = 4)
g <- ggplot(data = hist_plots, aes(`Selectivity measurements`), colors = "blue") + 
  #stat_count(color = "blue", fill = "light blue") +
  geom_histogram(breaks=c(seq(0, 50, by=1)), color = "blue", fill = "light blue") +
  scale_x_continuous(limits=c(0, 50), breaks=c(seq(0, 50, by=5)), labels=c(seq(0,45, by=5), "50+")) + 
  annotate("text", x = 35, y = 65, label = "Median = 4", size = 6) + ggtitle("Proteins assayed against") + xlab("") +
  theme(text = element_text(size = 16))
ggsave("selectivity plot.png", plot = g, dpi = 300, height = 4, width = 4)
median(hist_plots$`Selectivity measurements`)
median(hist_plots$`Total activities`)


# add probes from the chemical probes portal
sgcprobes <- read_csv("accessions.csv")
sgcprobes <- sgcprobes %>% select(CHEMICAL_ID, accession, pref_name) %>% mutate(group = 0)
#make full list
probe_list <- bind_rows(list(sgcprobes, probe_list))
probe_list <- rename(probe_list, UNIPROTKB = accession)

#get entrez ids and entry names for each target

require("UniProt.ws")

up <- UniProt.ws::UniProt.ws(taxId=9606)
res2 <- UniProt.ws::select(up,
                           keys = unique(probe_list$UNIPROTKB),
                           columns = c("ENTREZ_GENE", "ENTRY-NAME", "KEGG"),
                           keytype = "UNIPROTKB")
res2 <- res2[!duplicated(res2$UNIPROTKB), ]
res2$`ENTRY-NAME` <- stringr::str_split(res2$`ENTRY-NAME`, "_", simplify = TRUE)[,1]
res2$`ENTRY-NAME` <- stringr::str_split(res2$`ENTRY-NAME`, "_", simplify = TRUE)[,1]

# add the kegg and gene ontology information to targets

require("mygene")

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

#add protein information to probe list
probe_list <- left_join(probe_list, protein_information, by = c("UNIPROTKB" = "UNIPROTKB"))
probe_list <- probe_list  %>% mutate(source = "CHEMBL") %>% arrange(group, pref_name)
probe_list$source[probe_list$group == 0] <- "chemicalprobes.org"
#remove duplicates
probe_list <- probe_list[!duplicated(probe_list), ]
probe_list <- probe_list %>% dplyr::rename(TARGET_NAME = pref_name, DATA_SOURCE = source)
View(probe_list)
# write file to disk
write_csv(probe_list, "PROBELIST.csv", na = "NA")



length(unique(probe_list %>% filter(group != 0) %>% .$TARGET_NAME))
length(unique(probe_list %>% filter(group != 0) %>% .$UNIPROTKB))
length(unique(probe_list %>% filter(group != 0) %>% .$CHEMICAL_ID))
table(probe_list %>% filter(group != 0) %>% .$is_agonist)
table(probe_list %>% filter(group != 0) %>% .$orthogonal)
median(top_probes_onetarget$n_total)
median(top_probes_onetarget$n_select)

library(ggplot2)

