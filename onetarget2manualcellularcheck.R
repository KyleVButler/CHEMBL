#big loss of specificity is in checking the functional assay, could not include any results that have " E coli, baculovirus " or " Sf9 "  (include spaces)
#in the description   -others?  insect, bacteria, E Coli, should we exclude CHO? maybe not, some results have good data for that
#check some common probes and see how the funcitonal assays can be related to the target


library(dplyr)
library(tibble)
library(stringr)
library(readr)
#chembl_db <- src_sqlite("/Users/kylebutler/Desktop/chembl_21_sqlite/chembl_21.db", create = TRUE)
chembl_db <- src_sqlite("chembl_21.db", create = TRUE)

#only want to focus on common observations
names(table(activities_collected$standard_type))[as.vector(table(activities_collected$standard_type) > 100)]

#What about Log Ki, Log EC50, pKb - need to change to nM for standard_units == NA -- chembl should have done this, use code below if not
#activities_collected[activities_collected$standard_type == "Log Ki" & is.na(activities_collected$standard_units), ]$standard_value <-
#  10^(-(activities_collected[activities_collected$standard_type == "Log Ki" & is.na(activities_collected$standard_units), ]$standard_value)) * 10^9
#activities_collected[activities_collected$standard_type == "Log Ki" & is.na(activities_collected$standard_units), ]$standard_units <- "nM"

#want to annotate compounds with main target, number of targets assayed against, and total number of observations

compound_summary <- activities_collected %>% count(molregno) %>% dplyr::rename(n_total = n)
selectivity_measurements <- activities_collected %>% group_by(molregno) %>% 
  filter(standard_type %in% c("Selectivity ratio", "Ratio IC50", "Ratio", "Ratio Ki", "Ratio EC50", "Fold selectivity", "Selectivity Index", "Selectivity index",
                              "Relative potency", "Ratio pIC50", "Ratio pKi", "Selectivity") & assay_type == "B" & bao_format != "BAO_0000219") %>% 
  count(molregno) %>% dplyr::rename(n_selectivity = n)
#What about Log Ki, Log EC50, pKb - need to change to nM for standard_units == NA
potency_measurements <- activities_collected %>% group_by(molregno) %>% 
  filter(target_type == "SINGLE PROTEIN" & organism == "Homo sapiens" & confidence_score %in% c(8, 9) & standard_units %in% c("nM", "%") & 
           standard_type %in% c("Activity", "EC50", "IC50", "Kb", "KB", "Kd", "Ki", "Kinact", "Potency", "Log Ki", "Log EC50", "pKb", "Activity", "Inhibition")) %>% 
  distinct(pref_name) %>% count(molregno) %>% dplyr::rename(n_potency = n)
compound_summary <- full_join(compound_summary, selectivity_measurements)


compound_summary <- full_join(compound_summary, potency_measurements)
compound_summary[is.na(compound_summary)] <- 0

#keep only compounds with n_total >= 3 and n_select > 1
compound_summary <- compound_summary %>% mutate(n_select = n_selectivity + n_potency) %>% select(molregno, n_total, n_select) %>% 
  filter(n_total > 2 & n_select > 1)


# make activities list with selected compounds and remove cytochromes
activities_collected2 <- activities_collected %>% filter(molregno %in% compound_summary$molregno) %>% filter(substr(pref_name, 1, 7) != "Cytochr") 

#match each compound with its most potent target and keep only those that have one target under 100 nM
min_target <- activities_collected2 %>% 
  filter(target_type == "SINGLE PROTEIN" & organism == "Homo sapiens" & assay_type == "B" & confidence_score %in% c(8, 9) & standard_units == "nM" & 
           standard_type %in% c("Activity", "EC50", "IC50", "Kb", "KB", "Kd", "Ki", "Kinact", "Potency", "Log Ki", "Log EC50", "pKb") & standard_value < 100) %>%
  group_by(molregno, pref_name) %>% summarize(min_value = min(standard_value)) %>% group_by(molregno) %>% filter(n() == 1) 

compound_summary <- right_join(compound_summary, min_target)
compound_summary

#check that a cell assay exists for the main target
#activities_collected_subset <- right_join(activities_collected2, (compound_summary %>% select(molregno, pref_name) %>% rename(target = pref_name)))
#activities_collected_subset <- activities_collected_subset %>% filter(!grepl('insect|sf9|E Coli|escherichia|bacteria', description, ignore.case = TRUE)) %>% 
#  filter(bao_format == "BAO_0000219" & standard_units == "nM" & standard_value < 1000)
#activities_collected_subset <- activities_collected_subset[activities_collected_subset$pref_name == activities_collected_subset$target, ]
#compound_summary <- compound_summary %>% filter(molregno %in% unique(activities_collected_subset$molregno))
#compound_summary 

#check that no selectivity measurements are under 30
cmpdstokeep <- activities_collected2 %>% group_by(molregno) %>% 
  filter(standard_type %in% c("Selectivity ratio", "Ratio IC50", "Ratio", "Ratio Ki", "Ratio EC50", "Fold selectivity", "Selectivity Index", "Selectivity index",
                              "Relative potency", "Ratio pIC50", "Ratio pKi", "Selectivity") & assay_type == "B" & bao_format != "BAO_0000219" & standard_value < 30) 
compound_summary <- compound_summary %>% filter(!(molregno %in% unique(cmpdstokeep$molregno)))
compound_summary 

#check that ic50 type measurements have at least 30 fold selectivity
min_target <- activities_collected2 %>% filter(molregno %in% compound_summary$molregno) %>%
  filter(target_type == "SINGLE PROTEIN" & organism == "Homo sapiens" & confidence_score %in% c(8, 9) & standard_units == "nM" & 
           standard_type %in% c("Activity", "EC50", "IC50", "Kb", "KB", "Kd", "Ki", "Kinact", "Potency", "Log Ki", "Log EC50", "pKb")) %>%
  group_by(molregno, pref_name) %>% summarize(min_value = min(standard_value)) %>% group_by(molregno) %>% filter(n() > 1) %>% arrange(molregno, min_value)
min_target_check <- min_target %>% slice(1) %>% select(molregno, min_value) 
min_target_check2 <- min_target %>% slice(2) %>% select(molregno, min_value) %>% rename(second_value = min_value)
min_target_check <- full_join(min_target_check, min_target_check2)
min_target_check <- min_target_check %>% mutate(selectivity_ratio = second_value/min_value) %>% filter(selectivity_ratio < 30)
compound_summary <- compound_summary %>% filter(!(molregno %in% unique(min_target_check$molregno)))
compound_summary 
activities_collected2 <- activities_collected %>% filter(molregno %in% compound_summary$molregno)

#check that no % activity type measurements are above 50 for targets without known ic50 values
x <- rep(TRUE, nrow(compound_summary))
activities_collected_subset <- filter(activities_collected2, standard_type == "Activity" & standard_units == "%" & assay_type == "B" & target_type == "SINGLE PROTEIN")
activities_collected_subset2 <- filter(activities_collected2, target_type == "SINGLE PROTEIN" & organism == "Homo sapiens" & confidence_score %in% c(8, 9) & 
                                         standard_units == "nM" & 
                                         standard_type %in% c("Activity", "EC50", "IC50", "Kb", "KB", "Kd", "Ki", "Kinact", "Potency", "Log Ki", "Log EC50", "pKb"))
for(i in 1:nrow(compound_summary)){
  kilist <- filter(activities_collected_subset2, molregno == compound_summary$molregno[i])
  if(nrow(activities_collected_subset %>% filter(molregno == compound_summary$molregno[i] & standard_value < 50 & 
                                                 !(pref_name %in% kilist$pref_name))) >= 1){
    x[i] <- FALSE
  }
}
table(x)
compound_summary <- compound_summary[x, ]

#check that no % inhibition type measurements are above 50 for targets without known ic50 values
x <- rep(TRUE, nrow(compound_summary))
activities_collected_subset <- filter(activities_collected2, standard_type == "Inhibition" & standard_units == "%" & assay_type == "B" & target_type == "SINGLE PROTEIN")
activities_collected_subset2 <- filter(activities_collected2, target_type == "SINGLE PROTEIN" & organism == "Homo sapiens" & confidence_score %in% c(8, 9) & 
                                         standard_units == "nM" & 
                                         standard_type %in% c("Activity", "EC50", "IC50", "Kb", "KB", "Kd", "Ki", "Kinact", "Potency", "Log Ki", "Log EC50", "pKb"))
for(i in 1:nrow(compound_summary)){
  kilist <- filter(activities_collected_subset2, molregno == compound_summary$molregno[i])
  if(nrow(activities_collected_subset %>% filter(molregno == compound_summary$molregno[i] & standard_value > 50 & 
                                                 !(pref_name %in% kilist$pref_name))) >= 1){
    x[i] <- FALSE
  }
}
table(x)
compound_summary <- compound_summary[x, ]
compound_summary
length(unique(compound_summary$pref_name))

#find if probes are agonist or not
agonist_check <- compound_summary %>% select(pref_name, molregno) %>% rename(target_check = pref_name)
agonist_check <- left_join(activities_collected2, agonist_check)
agonist_check <- agonist_check %>% mutate(targets_identical = (target_check == pref_name)) %>% 
  filter(targets_identical == 1 & standard_units == "%" & standard_type %in% c("Emax", "max activation", "efficacy", "Efficacy") &
           standard_value > 25 & confidence_score %in% c(8,9)) 

compound_summary$is_agonist <- "NO"
table(compound_summary$is_agonist)
compound_summary$is_agonist[compound_summary$molregno %in% agonist_check$molregno] <- "YES"
table(compound_summary$is_agonist)

# add chembl ids and smiles
molecule_dictionary <- tbl(chembl_db, "molecule_dictionary")
join_vector <- molecule_dictionary %>% select(molregno, chembl_id) %>% filter(molregno %in% compound_summary$molregno)
join_vector <- collect(join_vector, n = Inf)
compound_summary <- left_join(compound_summary, join_vector, by = "molregno")
compound_structures <- tbl(chembl_db, "compound_structures")
join_vector <- compound_structures %>% select(molregno, canonical_smiles) %>% filter(molregno %in% compound_summary$molregno)
join_vector <- collect(join_vector, n = Inf)
compound_summary <- left_join(compound_summary, join_vector, by = "molregno")
#remove those with no smiles
compound_summary <- compound_summary[!(is.na(compound_summary$canonical_smiles)), ]

################## need to remove duplicated compounds
compound_summary <- compound_summary %>% group_by(molregno) %>% slice(1)



#rank probes based on number of observations, for each target, find tanimoto distance from best probe 
library(rcdk)

fun_get_tanimoto <- function(smilesin){
  query.mol <- parse.smiles(smilesin)[[1]]
  target.mols <- parse.smiles(smilesin)
  query.fp <- get.fingerprint(query.mol, type='maccs')
  target.fps <- lapply(target.mols, get.fingerprint, type='maccs')
  sims <- unname(unlist(lapply(target.fps, distance, fp2=query.fp, method='tanimoto')))
  return(sims)
}
compound_summary <- compound_summary %>% group_by(is_agonist, pref_name) %>% 
  arrange(is_agonist, pref_name, desc(n_total), desc(n_select), min_value) %>% mutate(tanimoto = fun_get_tanimoto(canonical_smiles))
View(agonists)
detach("package:rcdk", unload=TRUE)
detach("package:fingerprint", unload=TRUE)

high_select_probes <- compound_summary %>% group_by(is_aognist, pref_name, orthogonal) %>% 
  arrange(desc(n_select), desc(n_total), min_value) %>% slice(1)


#select probes with highest total number of observations and highest selectivity observations, and tiebreak on potency
high_select_probes <- compound_summary %>% group_by(pref_name) %>% arrange(pref_name, desc(n_select), desc(n_total), min_value) %>% slice(1)
hist(compound_summary$tanimoto)
#it is clear from distribution is bimodal and intersects at 0.75
high_total_probes_agonist <- compound_summary %>% filter(is_agonist == 1) %>% group_by(pref_name) %>% arrange(pref_name, desc(n_total), desc(n_select), min_value) %>% 
  slice(1)
high_total_probes_agonist_orthogonal <- compound_summary %>% filter(is_agonist == 1 & tanimoto < 0.75) %>% group_by(pref_name) %>% 
  arrange(pref_name, desc(n_total), desc(n_select), min_value) %>% slice(1)
high_total_probes_no_agonist <- compound_summary %>% filter(is_agonist != 1) %>% group_by(pref_name) %>% 
  arrange(pref_name, desc(n_total), desc(n_select), min_value) %>% slice(1)
high_total_probes_no_agonist_orthogonal <- compound_summary %>% filter(is_agonist != 1 & tanimoto < 0.75) %>% group_by(pref_name) %>% 
  arrange(pref_name, desc(n_total), desc(n_select), min_value) %>% slice(1)
top_probes_onetarget <- unique(rbind(high_total_probes_agonist, high_total_probes_no_agonist, high_total_probes_agonist_orthogonal, 
                                     high_total_probes_no_agonist_orthogonal))
View(top_probes_onetarget)
length(unique(top_probes_onetarget$pref_name))

#add accession to potential_probes
#join protein information
target_dictionary <- tbl(chembl_db, "target_dictionary")
join_vector <- target_dictionary %>% select(tid, pref_name, organism) %>% filter(organism == "Homo sapiens") %>%
  filter(pref_name %in% top_probes_onetarget$pref_name)
join_vector <- collect(join_vector, n = Inf)
top_probes_onetarget <- left_join(top_probes_onetarget, join_vector, by = "pref_name")
#join component id
target_components <- tbl(chembl_db, "target_components")
join_vector <- target_components %>% select(tid, component_id) %>% filter(tid %in% top_probes_onetarget$tid)
join_vector <- collect(join_vector, n = Inf)
top_probes_onetarget <- left_join(top_probes_onetarget, join_vector, by = "tid")
#join accession
component_sequences <- tbl(chembl_db, "component_sequences")
join_vector <- component_sequences %>% select(accession, component_id) %>% filter(component_id %in% top_probes_onetarget$component_id)
join_vector <- collect(join_vector, n = Inf)
top_probes_onetarget <- left_join(top_probes_onetarget, join_vector, by = "component_id")

top_probes_onetarget <- top_probes_onetarget %>% mutate(group = 1) %>% rename(CHEMICAL_ID = chembl_id)
