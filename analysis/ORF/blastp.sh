#!/bin/sh
### Set the job name (for your reference)
#PBS -N blastp
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

cd /scratch/bioschool/phd/blz227562/nm_nr/
echo "changed path"
export PATH=$PATH:/home/bioschool/phd/blz227562/ncbi-blast-2.14.0+/bin/
echo $PATH
which blastn

echo "------------------------------------------------------------------"
echo "Running for ATG"
cd /scratch/bioschool/phd/blz227562/nm_nr/
blastp -db hsa_hg38_prot -query ./ORFs/unique_ORFs/2024-05-05_fasta_uniqueORFs_bifunc_ATG.txt -out ./ORFs/unique_ORFs/15052024_blastp-fast_hybridGenes_uniqueORFs_ATG.tsv -task blastp-fast  -mt_mode 1 -num_threads 10 -max_target_seqs 15 -outfmt "6 qseqid qtitle sseqid stitle pident length mismatch gapopen qlen qstart qend slen sstart send evalue bitscore qcovs qcovhsp qcovus"
echo "blastp for ATG complete!"

echo "------------------------------------------------------------------"
echo "Running for CTG"
cd /scratch/bioschool/phd/blz227562/nm_nr/
blastp -db hsa_hg38_prot -query ./ORFs/unique_ORFs/2024-05-05_fasta_uniqueORFs_bifunc_CTG.txt -out ./ORFs/unique_ORFs/15052024_blastp-fast_hybridGenes_uniqueORFs_CTG.tsv -task blastp-fast  -mt_mode 1 -num_threads 10 -max_target_seqs 15 -outfmt "6 qseqid qtitle sseqid stitle pident length mismatch gapopen qlen qstart qend slen sstart send evalue bitscore qcovs qcovhsp qcovus"
echo "blastp for CTG complete!"

echo "------------------------------------------------------------------"
cd /scratch/bioschool/phd/blz227562/nm_nr/
blastp -db hsa_hg38_prot -query ./ORFs/unique_ORFs/2024-05-05_fasta_uniqueORFs_bifunc_GTG.txt -out ./ORFs/unique_ORFs/15052024_blastp-fast_hybridGenes_uniqueORFs_GTG.tsv -task blastp-fast  -mt_mode 1 -num_threads 10 -max_target_seqs 15 -outfmt "6 qseqid qtitle sseqid stitle pident length mismatch gapopen qlen qstart qend slen sstart send evalue bitscore qcovs qcovhsp qcovus"

echo "------------------------------------------------------------------"
cd /scratch/bioschool/phd/blz227562/nm_nr
blastp -db hsa_hg38_prot -query ./ORFs/unique_ORFs/2024-05-05_fasta_uniqueORFs_bifunc_TTG.txt -out ./ORFs/unique_ORFs/15052024_blastp-fast_hybridGenes_uniqueORFs_TTG.tsv -task blastp-fast  -mt_mode 1 -num_threads 10 -max_target_seqs 15 -outfmt "6 qseqid qtitle sseqid stitle pident length mismatch gapopen qlen qstart qend slen sstart send evalue bitscore qcovs qcovhsp qcovus"