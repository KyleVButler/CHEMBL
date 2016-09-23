# CHEMBL
These are R files for unbiased classification of chemical probes from the Chembl dataset. 
These require that chembl_21.db is running locally with sqlite3

the script makeinitiallist.R makes the main data file containing all compounds with in vitro potency < 100 nM

onetarget2manualcellularcheck.R compiles a probe list for probes with one protein target
twotarget2manualcellularcheck.R compiles a probe list for probes with two protein targets

In the above files, the option to automatically check on-target cell activity is commented out

combinelists.R combines the above with the chemicalprobes.org data into one list

protein_annotation.R is needed to add the protein annotations to the final list and also contains some functions for network visualization - this file needs further work
