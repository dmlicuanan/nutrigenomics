# install.packages("rentrez")

# load relevant packages 
library("rentrez")
library("data.table")
library("XML")

# set working directory
setwd("D:/Documents/Repositories/nutrigenomics/RDS")
setwd("C://Users//MPGL//Documents//Dea//nutrigenomics")

# load other functions needed
source("D:/Documents/Repositories/echinoderms/functions.R")

# references:
# https://peerj.com/preprints/3179v2.pdf
# https://cran.r-project.org/web/packages/rentrez/vignettes/rentrez_tutorial.html 

# entrez_dbs() gives you a list of database names
dataBase <- entrez_dbs()

# info about databases
resInfo <- lapply(dataBase, function(x) entrez_db_summary(x))
resInfo[[3]]
resInfo[[5]]
resInfo[[20]]

# relevant databases are nucleotide and nuccore
# get info about nucleotide database
entrez_db_searchable("nucleotide")
entrez_db_searchable("nuccore")
entrez_db_searchable("gene")
# these are still probably relevant databases for searching mRNA sequences

# relevant searchable fields for database 'nuccore'
# ALL 	 All terms from all searchable fields 
# ORGN 	 Scientific and common names of organism, and all higher levels of taxonomy 
# ACCN 	 Accession number of sequence 
# GENE 	 Name of gene associated with sequence 
# PROT 	 Name of protein associated with sequence 
# PORG 	 Scientific and common names of primary organism, and all higher levels of taxonomy 



################ attempt 1 (downloading seqs)
# download for Slc6a16  (has fewest seqs)

# character string of gene
gene <- "Slc6a16"

# find the unique identifier associated with Slc6a16 in Gene database
# the NCBI database dealing with genetic loci (rather than particular sequences) is called ‘Gene'
slc_gene <- entrez_search(db="gene", 
                          term=paste("Metazoa[ORGN] AND ", gene,"[GENE]", sep =""), 
                          retmax =99999,
                          use_history = TRUE)

# we now have a unique ids for Slc6a16 in the Gene database
# we need to find records from the NCBI Nucleotide database that are associated with the Gene record
# the function entrez_link can be used to find cross-referenced records
# a single call to entrez_link can identify Slc6a sequences in the nucleotide database in general and in a number of restrictive subsets of that database.
nuc_links <- entrez_link(dbfrom="gene", 
                         id = slc_gene$ids, 
                         db="nuccore", 
                         cmd = "neighbor_history")

# the RefSeq RNA subset on the Nucleotide database contains a curated set of mRNA transcripts for different genes
# the unique identifiers contained in the gene_nuccore_refseqrna element correspond to the sequences we wish to download
# the function entrez_fetch allows users to retrieve complete records in a variety of formats
# here the sequences are retrieved in the standard ‘fasta’ format, and returned as a character vector with a single element.
raw_recs <- entrez_fetch(db="nuccore",
                         web_history = nuc_links$web_histories$gene_nuccore_refseqrna,
                         rettype="fasta")

# inspect first few characters of long string
substr(raw_recs, start = 1, stop = 5000)



################ attempt 2 (downloading seqs; for large requests)
# character string of gene
gene <- "Slc6a15"

# find the unique identifier associated with Slc6a16 in Gene database
slc_gene <- entrez_search(db="gene", 
                          term=paste("Metazoa[ORGN] AND ", gene,"[GENE]", sep =""), 
                          retmax =99999,
                          use_history = TRUE)

# we now have a unique ids for Slc6a15 in the Gene database
# using id in the entrez_link function returns the ff error:
# Error in entrez_check(response) : 
#  HTTP failure 414, the request is too large. For large requests, try using web history as described in the rentrez tutorial

# we have to chunk the IDs into smaller sets (300 IDs will be given to entrez_link at a time)

# create empty list where ids, links and seqs will be entered
ls <- list()

# number of subset that will be created
nsub <- ceiling(length(slc_gene$ids)/300)

# fill in ids
for(i in 1:nsub) {
  ls$id[[i]] <- slc_gene$ids[300*(i-1)+1:300*i]
}

# remove NA IDs from last subset
ls$id[[nsub]] <- ls$id[[nsub]][!is.na(ls$id[[nsub]])]

# fill in links and seqs in list
for(i in 1:nsub) {
  # entrez_link
  ls$nuc_links[[i]] <- entrez_link(dbfrom="gene", 
                                   id = ls$id[[i]], 
                                   db="nuccore", 
                                   cmd = "neighbor_history")
  
  # entrez_fetch
  ls$raw_recs[[i]] <- entrez_fetch(db = "nuccore",
                                   rettype = "fasta",
                                   web_history = ls$nuc_links[[i]]$web_histories$gene_nuccore_refseqrna)
  
  # make data frame from sequences (seqs separated by "\n\n") 
  ls$df[[i]] <- data.frame(seq = unlist(data.table::tstrsplit(ls$raw_recs[[i]], "\n\n")))
}

# combine seqs of subsets into one vector
df <- data.frame(seq = unlist(ls$df))

# split id and sequence data 
for(i in 1:nrow(df)) {
  # splits fasta string 
  splitTemp <- unlist(strsplit(df$seq[i], "\n"))
  # takes id (first element in splitTemp)
  df$id[i] <- splitTemp[1]
  # pastes sequences together
  df$seq[i] <- paste(splitTemp[-1], collapse = "")
}

# save RDS for now (data frame with id and seq column)
# saveRDS(df, paste0(gene, ".RDS"))

# perform code above for all genes 

# list of RDS files to read
files <- paste("Slc6a", 15:20, ".RDS", sep = "")

# read all RDS files
rds <- rbindlist(lapply(files, function (x) data.table(readRDS(x))))

# reorder columns
rds <- rds[, c(2,1)]

# remove ">" in id column
rds$id <- sub(">", "", rds$id)

# Slc6a 15:20 combined in one rds
# saveRDS(rds, "Slc6a.RDS")

# write fasta file for combined rds
# writeFasta(id = rds$id, seq = rds$seq, file = "Slc6a.fas")

# read compiled seqs
slc <- readRDS("Slc6a.RDS")
# 3120 rows but only 3116 unique ids

# check if there are slc seqs for Chanos
slc$id[grep("chanos", tolower(slc$id))]



################ attempt 3 (automating downloading seqs for all genes)
# we will also get taxonomic info of the seqs

# list where entrez_search results will be stored
# transcript names from Figure 1 doi: 10.3389/fphys.2018.00212
searchRes <- list(gene = paste("Slc6a", 15:20, sep = "")) 

# empty dataframe to dump searchRes info
searchResdf <- data.frame()

# find the unique identifier associated with each gene in Gene database
# the NCBI database dealing with genetic loci (rather than particular sequences) is called 'Gene'
for(i in 1:length(searchRes$gene)) {
  # entrez_search()
  searchRes$geneRes[[i]] <- entrez_search(db="gene", 
                                          term=paste("Metazoa[ORGN] AND ", searchRes$gene[i],"[GENE]", sep =""), 
                                          retmax =99999)
  # fill in searchResdf with gene and geneid
  temp <- data.frame(gene = searchRes$gene[i], geneid = searchRes$geneRes[[i]]$ids)
  searchResdf <- rbind(searchResdf, temp)
  
  # print gene if search done
  print(searchRes$gene[i])
}

# unique IDs from gene database
geneidUse <- unique(searchResdf$geneid)

# we now have a unique ids from the Gene database
# we have to chunk the IDs into smaller sets (100 IDs will be given to entrez_link at a time)

# list where entrez_link and entrez_fetch results will be stored
fetchRes <- list()
# elements:
# geneid = unique ids from gene database
# nuccid = ids from nuccore database, result of entrez_link
# fasta = dataframe for fasta (columns: seq, id)

# desired number of gene ids per subset
nseq <- 100
# number of subsets that will be created
nsub <- ceiling(length(geneidUse)/nseq)

# fill in gene ids in fetchRes
for(i in 1:nsub) { fetchRes$geneid[[i]] <- geneidUse[(nseq*(i-1)+1):(nseq*i)] }

# remove NA IDs from last subset
fetchRes$geneid[[nsub]] <- fetchRes$geneid[[nsub]][!is.na(fetchRes$geneid[[nsub]])]

# add nuccore ids as elements in list
# by BATCH
for(i in 1:nsub) {
  # entrez_link
  temp <- entrez_link(dbfrom = "gene", 
                      id = fetchRes$geneid[[i]], 
                      db = "nuccore")
  
  # add nuccore ids to list
  fetchRes$nuccid[[i]] <- temp$links$gene_nuccore_refseqrna
  
  # print progress
  print(paste0("Progress: ", round((i/nsub)*100, 1), "%"))
}


# start the clock!
ptm <- proc.time()

# we download sequences based on BATCHES of nuccore ids in fetchRes
# time: 62.81 secs
for(i in 1:nsub) {  
  # entrez_fetch
  temp <- entrez_fetch(db = "nuccore",
                       rettype = "fasta",
                       id =  fetchRes$nuccid[[i]])
  
  # make data frame from sequences (seqs separated by "\n\n") 
  df <- data.frame(seq = unlist(data.table::tstrsplit(temp, "\n\n")))
  
  # split id and sequence data of df
  for(j in 1:nrow(df)) {
    # splits fasta string 
    splitTemp <- unlist(strsplit(df$seq[j], "\n"))
    # takes id (first element in splitTemp)
    df$id[j] <- splitTemp[1]
    # pastes sequences together
    df$seq[j] <- paste(splitTemp[-1], collapse = "")
  }
  
  # save df to fetchRes
  fetchRes$fasta[[i]] <- df
  
  # print progress
  print(paste0("Progress: ", round((i/nsub)*100, 1), "%"))
}

# stop the clock
proc.time() - ptm

# save RDS
# saveRDS(fetchRes, "fetchRes.RDS")

# combine seqs of subsets into one vector
allfasta <- do.call("rbind", fetchRes$fasta)
# remove duplicated seqs and reorder columns
allfasta <- allfasta[!duplicated(allfasta$id), c("id","seq")]

# remove ">" in id column
allfasta$id <- sub(">", "", allfasta$id)

# write fasta file 
# writeFasta(id = allfasta$id, seq = allfasta$seq, file = "Slc6a_v2.fas")

# save RDS
# saveRDS(allfasta, "allfasta.RDS")

# preview fasta
# readLines("Slc6a_v2.fas", 3)


# now that we have sequences, we need the taxonomic info associated with the seqs
# so we create a dataframe where the nuccore id is linked to taxid and accession version
# taxid = key for taxonomic metadata
# accv = key for linking to fasta
seqtotax <- data.frame(nuccid = unique(unlist(fetchRes$nuccid)), taxid = NA, accv = NA)

# start the clock!
ptm <- proc.time()

# find the taxid and accession version for each nuccid
# not per batch (per nuccore id)
# time to run:  user (45.09), system (1.82), elapsed (2846.92 secs?)
for(i in 1:nrow(seqtotax)) {
  # entrez_summary
  temp <- entrez_summary(db = "nuccore",
                         id = seqtotax$nuccid[i])
  # taxid
  seqtotax$taxid[i] <- temp$taxid
  # accession version
  seqtotax$accv[i] <- temp$accessionversion
  
  # percent progress
  perc <- round((i/nrow(seqtotax))*100, 1)
  # print progress
  if(perc %in% seq(5, 100, 5)) { print(paste0("Progress: ", perc, "%")) } 
}

# stop the clock
proc.time() - ptm

# save seqtotax.RDS
# saveRDS(seqtotax, "seqtotax.RDS")

# read seqtotax.RDS
seqtotax <- readRDS("seqtotax.RDS")

# get unique taxids and add to data.frame
taxinfo <- data.frame(taxid = unique(na.omit(seqtotax$taxid)))

# taxonomic levels of interest
taxlevs <- c("kingdom", "phylum", "class","order", "family", "genus", "species", "subspecies")

# add columns to df for each taxonomic level of interest
taxinfo <- cbind(taxinfo, setNames(lapply(taxlevs, function(x) x=NA), taxlevs))

# start the clock!
ptm <- proc.time()

# fill in columns of taxinfo df
# time elapsed: 461.67 secs (7.69 mins)
for(i in 1:nrow(taxinfo)) {
  # entrez_fetch will output XML
  tax_rec <- entrez_fetch(db = "taxonomy",
                          id = taxinfo$taxid[i],
                          rettype = "xml", 
                          parsed = TRUE) 
  
  # get taxa 
  taxon <- xpathSApply(tax_rec, "//ScientificName", XML::xmlValue)
  # get rank of taxa
  rank <- xpathSApply(tax_rec, "//Rank", XML::xmlValue)
  # df of taxa and tank
  ranklist <- data.frame(taxon, rank)
  
  # paste taxa to taxinfo dataframe
  taxinfo[i,2:9] <- ranklist$taxon[match(taxlevs, ranklist$rank)]
  
  # percent progress
  perc <- round((i/nrow(taxinfo))*100, 1)
  
  # print progress
  if(perc %in% seq(5, 100, 5)) { print(paste0("Progress: ", perc, "%")) } 
}

# stop the clock
proc.time() - ptm

# save RDS
# saveRDS(taxinfo, "taxinfo.RDS")



################ distribution of seq across taxa

# readin relevant RDS
allfasta <- readRDS("allfasta.RDS")
seqtotax <- readRDS("seqtotax.RDS")
taxinfo <- readRDS("taxinfo.RDS")

# extract accession version from id
allfasta$accv <- unlist(lapply(allfasta$id, function(x) strsplit(x, " ")[[1]][1]))

# update with taxid
allfasta$taxid <- seqtotax$taxid[match(allfasta$accv, seqtotax$accv)]

# fill in taxonomic details in allfasta
allfasta <- merge(allfasta, taxinfo, by = "taxid" )

# save RDS
# saveRDS(allfasta, "allfasta.RDS")

# distribution across classes
table(allfasta[allfasta$phylum == "Chordata", ]$class)
barplot(table(allfasta[allfasta$phylum == "Chordata", ]$class))

################ reduce the sequences

# create two groups of sequences
# 1: represents all species from class Actinopteri
# 2: contains only species which have sequences for at least the 4 genes present represented in Chanos chanos: slc6a15, 17, 19, 20

# readin relevant RDS
allfasta <- readRDS("allfasta.RDS")

# check genes represented in Chanos chanos
allfasta[allfasta$species == "Chanos chanos","id"]
# genes: slc6a15, 17, 18, 20
# check length of sequences represented in Chanos chanos
nchar(allfasta[allfasta$species == "Chanos chanos","seq"])

# only retain those with sequence length between 1800-5000
redseq <- allfasta[(nchar(allfasta$seq) >= 1800) & (nchar(allfasta$seq) <= 5000) & (!is.na(allfasta$kingdom)),]

# add column for gene 
redseq$gene <- ifelse(grepl("slc6a15|B\\(0\\)AT2|member 15", redseq$id, ignore.case = T), "slc6a15",
                      ifelse(grepl("slc6a16|NTT5", redseq$id, ignore.case = T), "slc6a16",
                             ifelse(grepl("slc6a17|member 17", redseq$id, ignore.case = T), "slc6a17",
                                    ifelse(grepl("slc6a18|B\\(0\\)AT3|member 18", redseq$id, ignore.case = T), "slc6a18",
                                           ifelse(grepl("slc6a19|B\\(0\\)AT1|member 19", redseq$id, ignore.case = T), "slc6a19",
                                                  ifelse(grepl("slc6a20|XTRP3|member 20", redseq$id, ignore.case = T), "slc6a20", "other"))))))
# check distribution of sequences across classes that are labelled as "other
table(redseq[redseq$gene == "other", ]$class, useNA = "ifany")
table(redseq[redseq$gene == "other" & redseq$class =="Actinopteri", ]$order, useNA = "ifany")

# remove sequences with gene = "other" and not from class Actinopteri
fishseq <- as.data.table(redseq[(!redseq$gene == "other") & (redseq$class == "Actinopteri") & (!is.na(redseq$class)),])

# determine the number of genes represented for each species
# so we find the unique combinations of species and gene
comb <- unique(fishseq[, c("species", "gene")])
# number of genes per species
count <- comb[, .N ,by = species]
# list of fish species with 4 or more genes represented (still have to refine since genes may not be 15, 17, 19, 20)
fourormore <- count[count$N >= 4,]$species
# reduce fishseq to entries with four or more genes represented
chanseq <- fishseq[fishseq$species %in% fourormore,]

# find species that do not have sequences for slc6a15, 17, 19, 20
changen <- c("slc6a15", "slc6a17", "slc6a19", "slc6a20")
# subset chanseq to show those genes only
temp <- chanseq[chanseq$gene %in% changen,]
# get counts per gene
counts <- as.data.table(table(temp$species, temp$gene))
colnames(counts) <- c("species", "gene", "N")
counts[order(counts$species, counts$gene)]
# species with 0 for any of the four genes
counts[counts$N == 0,]
# species with 0 sequence counts for relevant genes
inc <- unique(counts[counts$N == 0,]$species) 

# remove from chanseq entries with no sequence for any one of slc6a15, 17, 19, 20
chanseq <- chanseq[!(chanseq$species %in% inc),]
# check if species have more than one sequence per gene
chanseq[, (.N), by = .(species, gene)]

# get only one sequence for each gene for each species
chanseq <- chanseq[!duplicated(chanseq[,c("species", "gene")]),]
# chanseq above contains only species which have sequences for at least the 4 genes present represented in Chanos chanos: slc6a15, 17, 19, 20
# there is also only one sequence selected per gene
# to chanseq, we need to append sequences of the outgroup (Mammalia)

# now, we create a subset that represents all species from class Actinopteri
# like chanseq, there should only be one sequence selected per gene per species
# we start with fish seq, and just remove duplicates (based on species and gene column)

# check if species have more than one sequence per gene
fishseq[, (.N), by = .(species, gene)]

# remove duplicates based on species and gene column 
fishseq <- fishseq[!duplicated(fishseq[,c("species", "gene")]),]
# to fishseq above, we need to append outgroup sequences (Mammalia)

# outgroup
# remove sequences with gene = "other" and not from class Actinopteri
mamseq <- as.data.table(redseq[(!redseq$gene == "other") & (redseq$class == "Mammalia") & (!is.na(redseq$class)),])

# determine the number of genes represented for each species
# find the unique combinations of species and gene
comb <- unique(mamseq[, c("species", "gene")])
# number of genes per species
count <- comb[, .N ,by = species]
# all 6 genes must be represented so we find the species that satisfy that
six <- count[count$N == 6,]$species
# reduce mamseq to entries with all six genes represented
mamseq <- mamseq[mamseq$species %in% six,]

# find four species with the most number of sequences
most <- mamseq[, .N, by = species][order(N, decreasing = TRUE)]$species[1:4]

# reduce mamseq to just those four species
mamseq <- mamseq[mamseq$species %in% most,]
# get just one sequence per gene
# remove duplicates based on species and gene column 
mamseq <- mamseq[!duplicated(mamseq[,c("species", "gene")]),]

# rbind mamseq with chanseq and fishseq
# chanseq + mamseq
fish4genes <- rbind(chanseq, mamseq)
# saveRDS(fish4genes, "fish4genes.RDS")
# fish4genes <- readRDS("fish4genes.RDS")

# fishseq + mamseq
fishallspecies <- rbind(fishseq, mamseq)
# saveRDS(fishallspecies, "fishallspecies.RDS")
# fishallspecies <- readRDS("fishallspecies.RDS")

# JM wrote code for breaking sequences into pieces 
# script: "D:/Documents/Repositories/nutrigenomics/breakforblast.R"
# outputs are "seqbreaks4blast.RDS" and "seqbreaks4blastall.RDS" in RDS folder
# fish only - seqbreaks4blast.rds
# all species - seqbreaks4blastall.rds

# read in outputs
seqbreaks4blast <- readRDS("seqbreaks4blast.RDS") 
seqbreaks4blastall <- readRDS("seqbreaks4blastall.RDS")

# create fasta files from outputs
# but first combine id and breakNum
writeFasta(id = paste0(seqbreaks4blast$id, " breakNum_", seqbreaks4blast$breakNum), seq = seqbreaks4blast$breakSeq, file = "seqbreaks4blast.fas")
writeFasta(id = paste0(seqbreaks4blastall$id, " breakNum_", seqbreaks4blastall$breakNum), seq = seqbreaks4blastall$breakSeq, file = "seqbreaks4blastall.fas")




