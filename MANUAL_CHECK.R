potential_probes_onetarget <- mutate(top_probes_onetarget, isprobe = "M")

manual_check <- potential_probes_onetarget %>% group_by(pref_name) %>% select(molregno, pref_name, CHEMICAL_ID)

manual_check$check <- rep(NA, nrow(manual_check))

for(i in seq(from = 33, to = nrow(manual_check), by = 20)){
  x <- data.frame((activities_collected %>% filter(molregno == manual_check$molregno[i]) %>% filter(organism == "Homo sapiens" | is.na(organism)) %>%
         select(molregno, standard_value, standard_units, standard_type, assay_type, pref_name, description, target_type, bao_format, confidence_score))) 
  View(x)
  print(manual_check[manual_check$molregno == manual_check$molregno[i], ])
  manual_check$check[i] <- readline("Is it a probe?")
}

##########
potential_probes_twotarget <- mutate(potential_probes_twotarget, isprobe = "M")

manual_check <- potential_probes_twotarget %>% group_by(pref_name, n) %>% select(molregno, pref_name, n, group, chembl_id)

manual_check$check <- rep(NA, nrow(manual_check))

for(i in 1:nrow(manual_check)){
  x <- data.frame((activities_collected %>% filter(molregno == manual_check$molregno[i]) %>% filter(organism == "Homo sapiens" | is.na(organism)) %>%
                     select(molregno, standard_value, standard_units, standard_type, assay_type, pref_name, description, target_type, bao_format, confidence_score))) 
  View(x)
  print(manual_check[manual_check$molregno == manual_check$molregno[i], ])
  manual_check$check[i] <- readline("Is it a probe?")
}
