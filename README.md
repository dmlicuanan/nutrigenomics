# Milkfish nutrigenomics
nutriseqs.R is an examination of available sequence data for the Slc6a family of transporters, whose functions in relation to milkfish (*Chanos chanos*) nutrigenomics will be investigated. The script performs the following tasks: 
- Automate downloading of metazoan DNA sequences for genes Slc6a15, Slc6a16, Slc6a17, Slc6a18, Slc6a19, and Slc6a20
- Link taxonomic information to the downloaded sequences
- Examine the distribution of sequences across chordate classes
- Check which Slc6a genes have milkfish-labelled sequences
- Create a reduced sequence dataset such that only Actinopteri species with sequences for at least the four milkfish Slc6a genes are included
- Create a reduced sequence data set such that only Actinopteri species are included
- Create a reduced sequence data set such that only Mammalia species with sequences for all six Slc6a genes are included (serves as outgroup)
  
nutrigenomics_blast.sh builds a BLAST database from the milkfish genome and queries cut Slc6a sequences against the database

functions.R contains functions used in nutriseqs.R for reading and writing FASTA files

***

Author's note: 

Scripts were written in November 2021. Latest code edits: 2023-10-18.
