import re
import pandas as pd
import numpy as np

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
    dat_gene=pd.read_csv("18122023_GRCh38.p14_gene_accessions.csv")
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


format_orf(infile="./ORFs/23122023_allpep_ATG-bifunc_genes.fa",
           outfile="./ORFs/23122023_formatted_allpep_ATG_genes.csv")

format_orf(infile="./ORFs/23122023_allpep_ATG-TTG-GTG-CTC-bifunc_genes.fa",
           outfile="./ORFs/23122023_formatted_allpep_ATG-TTG-GTG-CTC-bifunc_genes.csv")

