# change directory to location of Chanos chanos genome
cd /mnt/d/Documents/Nutrigenomics/Blast/

# unzip genome file
# downloaded from NCBI
gunzip fChaCha1.pri.cur.20190805.fasta.gz

# load environment (where ncbi-blast is installed)
conda activate obiUse

# build a blast database from the genome and output as genome.fasta
makeblastdb -in fChaCha1.pri.cur.20190805.fasta -out genome.fasta -parse_seqids -dbtype nucl

# blast cut Slc6a sequences against milkfish genome
blastn -db genome.fasta -query seqbreaks4blast.fas -max_target_seqs 20 -outfmt '6 qaccver saccver stitle pident length mismatch gapopen qstart qend sstart send'    
