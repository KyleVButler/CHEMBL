# CHEMBL
These are R files for unbiased classification of chemical probes from the Chembl dataset. 
These require that chembl_21.db is in the working directory

the script makeinitiallist.R makes the main data file containing all compounds with in vitro potency < 100 nM

onetarget2manualcellularcheck.R compiles a probe list for probes with one protein target
twotarget2manualcellularcheck.R compiles a probe list for probes with two protein targets

combinelists.R combines the above with the chemicalprobes.org data into one list and adds GO, kegg, entrez, uniprot annotations
