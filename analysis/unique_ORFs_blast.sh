#!/bin/sh
### Set the job name (for your reference)
#PBS -N unique_ORFs_blast
### Set the project name, your department code by default
#PBS -P menon.onetime
### Request email when job begins and ends
#PBS -m bea
### Specify email address to use for notification.
#PBS -M $USER@iitd.ac.in
####
#PBS -l select=1:ncpus=10
### Specify "wallclock time" required for this job, hhh:mm:ss
#PBS -l walltime=30:00:00


cd /scratch/bioschool/phd/blz227562/nm_nr/ORFs
module load compiler/R/4.1.2/intel2020
echo "Running for ATG"
Rscript get_unique_ORFs.R ATG 03052024_formatted_pep_ATG_bifunc_genes.csv bifunc

cd /scratch/bioschool/phd/blz227562/nm_nr/ORFs/unique_ORFs
export PATH=$PATH:/home/bioschool/phd/blz227562/ncbi-blast-2.14.0+/bin/
echo $PATH
which blastn
blastp -db hsa_hg38_prot -query 2024-05-05_fasta_uniqueORFs_bifunc_ATG.txt -out 05052024_blastp-fast_hybridGenes_uniqueORFs_ATG.csv -task blastp-fast  -mt_mode 1 -num_threads 10 -max_target_seqs 15 -outfmt "10 qseqid qtitle sseqid stitle pident length mismatch gapopen qlen qstart qend slen sstart send evalue bitscore qcovs qcovhsp qcovus"

echo "------------------------------------------------------------------"
echo "Running for CTG"
cd /scratch/bioschool/phd/blz227562/nm_nr/ORFs
Rscript get_unique_ORFs.R CTG 03052024_formatted_pep_CTG_bifunc_genes.csv bifunc

cd /scratch/bioschool/phd/blz227562/nm_nr/ORFs/unique_ORFs
blastp -db hsa_hg38_prot -query 2024-05-05_fasta_uniqueORFs_bifunc_CTG.txt -out 05052024_blastp-fast_hybridGenes_uniqueORFs_CTG.csv -task blastp-fast  -mt_mode 1 -num_threads 10 -max_target_seqs 15 -outfmt "10 qseqid qtitle sseqid stitle pident length mismatch gapopen qlen qstart qend slen sstart send evalue bitscore qcovs qcovhsp qcovus"

echo "------------------------------------------------------------------"
echo "Running for GTG"
cd /scratch/bioschool/phd/blz227562/nm_nr/ORFs
Rscript get_unique_ORFs.R GTG 03052024_formatted_pep_GTG_bifunc_genes.csv bifunc

cd /scratch/bioschool/phd/blz227562/nm_nr/ORFs/unique_ORFs
blastp -db hsa_hg38_prot -query 2024-05-05_fasta_uniqueORFs_bifunc_GTG.txt -out 05052024_blastp-fast_hybridGenes_uniqueORFs_GTG.csv -task blastp-fast  -mt_mode 1 -num_threads 10 -max_target_seqs 15 -outfmt "10 qseqid qtitle sseqid stitle pident length mismatch gapopen qlen qstart qend slen sstart send evalue bitscore qcovs qcovhsp qcovus"

echo "------------------------------------------------------------------"
echo "Running for TTG"
cd /scratch/bioschool/phd/blz227562/nm_nr/ORFs
Rscript get_unique_ORFs.R TTG 03052024_formatted_pep_TTG_bifunc_genes.csv bifunc

cd /scratch/bioschool/phd/blz227562/nm_nr/ORFs/unique_ORFs
blastp -db hsa_hg38_prot -query 2024-05-05_fasta_uniqueORFs_bifunc_TTG.txt -out 05052024_blastp-fast_hybridGenes_uniqueORFs_TTG.csv -task blastp-fast  -mt_mode 1 -num_threads 10 -max_target_seqs 15 -outfmt "10 qseqid qtitle sseqid stitle pident length mismatch gapopen qlen qstart qend slen sstart send evalue bitscore qcovs qcovhsp qcovus"