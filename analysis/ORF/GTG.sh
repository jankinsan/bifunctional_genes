#!/bin/sh
### Set the job name (for your reference)
#PBS -N GTG_format
### Set the project name, your department code by default
#PBS -P menon.onetime
### Request email when job begins and ends
#PBS -m bea
### Specify email address to use for notification.
#PBS -M $USER@iitd.ac.in
####
#PBS -l select=1:ncpus=10
### Specify "wallclock time" required for this job, hhh:mm:ss
#PBS -l walltime=10:00:00

cd /scratch/bioschool/phd/blz227562/nm_nr/
conda activate nm_nr
python

import re
import pandas as pd
import numpy as np
import argparse

#trial basis only to show the output and whether it is worth doing
#already run on HPC
def format_orf(infile, outfile):
    #read the file with predicted ORFs
    file_pep=open(infile, 'r')
    #read file data
    lines= file_pep.readlines()
    print(len(lines))
    print("File read!")
    file_pep.close()#open gene symbols and accession document
    dat_gene=pd.read_csv("/scratch/bioschool/phd/blz227562/nm_nr/18122023_GRCh38.p14_gene_accessions.csv")
    #empty dataframe for information
    orf_info=pd.DataFrame(columns=['acc_orf', 'gene', 'orf_number', 'coordinates', 'type', 'length', 'frame', 'start', 'stop', 'orf_sequence'])
    #get the information out for sequences one by one
    for line in lines:
        line=re.split("\n", line)[0]
        if line.startswith('>'): #encountered a new ORF
            #get new gene information
            split=re.split(" ", line)
            split[1]=re.split("\\(", split[1])[0] #removing strand information
            split.insert(1, re.split("_", split[0])[2]) #adding orf number in a new place
            split[0]=re.split("_ORF", split[0])[0] #getting only the accession
            pos=np.where(dat_gene["accession"] == split[0])
            split.insert(1, dat_gene.loc[pos[0][0], "gene"]) #adding the gene symbol too
            split.append("")
            idx=len(orf_info) #setting an index
            orf_info.loc[idx]=split
            print('index:' + str(idx))
        else:
            orf_info.loc[idx, 'orf_sequence'] = orf_info.loc[idx, 'orf_sequence'] + line
    orf_info.to_csv(outfile)

print("Running for GTG..")
format_orf(infile="/scratch/bioschool/phd/blz227562/nm_nr/coding_transcripts/orfs/15052024_pep_FTG_coding_bifunc_genes.fa",
           outfile="/scratch/bioschool/phd/blz227562/nm_nr/coding_transcripts/orfs/15072024_formatted_allGenes_codingPEP_GTGorfs.faa")