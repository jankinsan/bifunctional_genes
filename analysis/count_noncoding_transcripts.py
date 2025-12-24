#HPC COMMANDS BEFORE RUNNING R
#conda activate nm_nr
#srun --partition=compute --nodes=1 --ntasks-per-node=1 --time=12:00:00 --pty bash -i
#cd /lustre/sonam.dhamija/janki_nr/
#python

import re
import pandas as pd

#open file
file1=open("GCF_000001405.40_GRCh38.p14_rna.fna", 'r')
#read file data
Lines= file1.readlines()
print("File read!")

gene_list =[]
fasta_details_list=[]
#check for the '>' at the start of the seqeunce: indicates new sequence start
# here we need to quantify all accession and all transcripts related to a gene
for line in Lines:
    x=re.search("^>", line)
    if x is not None:
        split=re.split(", ", line)
        split_idx = [i for i in range(len(split)) if re.search("transcript variant", split[i]) is not None]
        if len(split_idx)==1:
            split=re.split("\(", split[split_idx[0]-1])
        else:
            split=re.split("\(", split[-2])
        split=re.split("\)", split[-1])
        print(split[0])
        fasta_details_list.append(line)
        gene_list.append(split[0]) #will allows us to get gene symbols for later use!

print("Gene List made!")
#find unique genes
genes_unique= pd.unique(gene_list)
print(genes_unique[0:100])
print(len(genes_unique))
#make a pandas data frame with gene symbols and number of transcripts for XR_, XM_, NM_, NP_ and total 
human_df= pd.DataFrame(0, columns= ["total", "NM", "NR", "XM", "XR"], index= genes_unique)
gene_acc= pd.DataFrame(0, columns= ["accession", "gene"], index=range(0, len(fasta_details_list)))
i=0
for line in fasta_details_list:
    print(line)
    split=re.split(",", line)
    split_idx = [i for i in range(len(split)) if re.search(" transcript variant", split[i]) is not None]
    if len(split_idx)==1:
        split=re.split("\(", split[split_idx[0]-1])
    else:
        split=re.split("\(", split[-2])
    split=re.split("\)", split[-1])
    gene_acc.loc[i, "accession"]=re.split(" ", line)[0]
    gene_acc.loc[i, "gene"]=split[0]
    if re.search(">NM_", line):
        human_df.loc[split[0], 'NM']+= 1
        human_df.loc[split[0], 'total']+= 1
    elif re.search(">NR_", line):
        human_df.loc[split[0], 'NR']+= 1
        human_df.loc[split[0], 'total']+= 1
    elif re.search(">XM_", line):
        human_df.loc[split[0], 'XM']+= 1
        human_df.loc[split[0], 'total']+= 1
    elif re.search(">XR_", line):
        human_df.loc[split[0], 'XR']+= 1
        human_df.loc[split[0], 'total']+= 1
    i=i+1
           
#printing
print(human_df.head())
print("Total number of lines in the data: ", len(Lines))
print("Number of trnascripts found (including all accessions:", len(fasta_details_list))
print("Unique genes found:", len(genes_unique))

#writing to file
human_df["NM_XM"]=human_df["NM"]+human_df["XM"]
human_df["NR_XR"]=human_df["NR"]+human_df["XR"]
human_df.to_csv("16122023_GRCh38.p14_genewise_counts.csv")
gene_acc.to_csv("18122023_GRCh38.p14_gene_accessions.csv")

human_df.to_csv("21032024_GRCh38.p14_genewise_counts.csv")
gene_acc.to_csv("21032024_GRCh38.p14_gene_accessions.csv")

file_fasta = open("21032024_hg38_fasta_details_list.txt", "w+")
file_fasta.writelines(fasta_details_list)
file_fasta.close()
file1.close()

file_fasta = open("hg38_fasta_details_list.txt", "w+")
file_fasta.writelines(fasta_details_list)
file_fasta.close()
file1.close()

#open file for T2T
file2=open("./16122023_t2t/rna.fna", 'r')
#read file data
Lines= file2.readlines()
print("File read!")

gene_list =[]
fasta_details_list=[]
#check for the '>' at the start of the seqeunce: indicates new sequence start
# here we need to quantify all accession and all transcripts related to a gene
for line in Lines:
    x=re.search("^>", line)
    if x is not None:
        split=re.split(",", line)
        split_idx = [i for i in range(len(split)) if re.search(" transcript variant", split[i]) is not None]
        if len(split_idx)==1:
            split=re.split("\(", split[split_idx[0]-1])
        else:
            split=re.split("\(", split[-2])
        split=re.split("\)", split[-1])
        print(split[0])
        fasta_details_list.append(line)
        gene_list.append(split[0]) #will allows us to get gene symbols for later use!

print("Gene List made!")
#find unique genes
genes_unique= pd.unique(gene_list)
print(genes_unique[0:100])
#make a pandas data frame with gene symbols and number of transcripts for XR_, XM_, NM_, NP_ and total 
human_df= pd.DataFrame(0, columns= ["total", "NM", "NR", "XM", "XR"], index= genes_unique)
gene_acc= pd.DataFrame(0, columns= ["accession", "gene"])
i=0
for line in fasta_details_list:
    print(line)
    split=re.split(", ", line)
    split_idx = [i for i in range(len(split)) if re.search("transcript variant", split[i]) is not None]
    if len(split_idx)==1:
        split=re.split("\(", split[split_idx[0]-1])
    else:
        split=re.split("\(", split[-2])
    split=re.split("\)", split[-1])
    gene_acc.loc[i, "accession"]=re.split(" ", line)[0]
    gene_acc.loc[i, "gene"]=split[0]
    if re.search(">NM_", line):
        human_df.loc[split[0], 'NM']+= 1
        human_df.loc[split[0], 'total']+= 1
    elif re.search(">NR_", line):
        human_df.loc[split[0], 'NR']+= 1
        human_df.loc[split[0], 'total']+= 1
    elif re.search(">XM_", line):
        human_df.loc[split[0], 'XM']+= 1
        human_df.loc[split[0], 'total']+= 1
    elif re.search(">XR_", line):
        human_df.loc[split[0], 'XR']+= 1
        human_df.loc[split[0], 'total']+= 1
    i=i+1
    
#printing
print(human_df.head())
print("Total number of lines in the data: ", len(Lines))
print("Number of trnascripts found (including all accessions:", len(fasta_details_list))
print("Unique genes found:", len(genes_unique))

#writing to file
human_df["NM_XM"]=human_df["NM"]+human_df["XM"]
human_df["NR_XR"]=human_df["NR"]+human_df["XR"]
human_df.to_csv("18122023_T2T-CHM13v2.0_genewise_counts.csv")
gene_acc.to_csv("18122023_T2T-CHM13v2.0_gene_accessions.csv")

file_fasta = open("T2T-CHM13v2.0_fasta_details_list.txt", "w+")
file_fasta.writelines(fasta_details_list)
file_fasta.close()
file2.close()


#which function
#link: https://alexmiller.phd/posts/python-pandas-which-function-indices-similar-to-R/#:~:text=If%20you're%20coming%20from,this%20as%20a%20good%20practice.
def which(self):
    try:
        self = list(iter(self))
    except TypeError as e:
        raise Exception("""'which' method can only be applied to iterables.
        {}""".format(str(e)))
    indices = [i for i, x in enumerate(self) if bool(x) == True]
    return(indices)

#### script to count transcripts for any species given the rna.fna file from NCBI!!!! 
##ON IIT HPC, after job resources have been allocated
# conda activate nm_nr
# cd /home/
# python
import re
import pandas as pd

###MAKING A FUNCTION FOR OTHER SPECIES!
def count_transcripts(rnaFilePath, outputFileName, outputFileAcc, outputFileFasta):
    file=open(rnaFilePath, 'r')
    Lines= file.readlines()
    print("File read!")
    gene_list =[]
    fasta_details_list=[]
    #check for the '>' at the start of the sequence: indicates new sequence start
    #here we need to quantify all accession and all transcripts related to a gene
    for line in Lines:
        x=re.search("^>", line)
        if x is not None:
            split=re.split(",", line)
            split_idx = [i for i in range(len(split)) if re.search(" transcript variant", split[i]) is not None]
            if len(split_idx)==1:
                split=re.split("\(", split[split_idx[0]-1])
            else:
                split=re.split("\(", split[-2])
            split=re.split("\)", split[-1])
            print(split[0])
            fasta_details_list.append(line)
            gene_list.append(split[0])
    print("Gene List made!")
    #find unique genes
    genes_unique= pd.unique(gene_list)
    print(genes_unique[0:100])
    #make a pandas data frame with gene symbols and number of transcripts for XR_, XM_, NM_, NP_ and total
    org_df= pd.DataFrame(0, columns= ["total", "NM", "NR", "XM", "XR"], index= genes_unique)
    gene_acc= pd.DataFrame(columns= ["accession", "gene"])
    i=0
    for line in fasta_details_list:
        print(line)
        split=re.split(",", line)
        split_idx = [i for i in range(len(split)) if re.search(" transcript variant", split[i]) is not None]
        i=i+1
        if len(split_idx)==1:
            split=re.split("\(", split[split_idx[0]-1])
        else:
            split=re.split("\(", split[-2])
        split=re.split("\)", split[-1])
        gene_acc.loc[i, "accession"]=re.split(" ", line)[0]
        gene_acc.loc[i, "gene"]=split[0]
        if re.search(">NM_", line):
            org_df.loc[split[0], 'NM']+= 1
            org_df.loc[split[0], 'total']+= 1
        elif re.search(">NR_", line):
            org_df.loc[split[0], 'NR']+= 1
            org_df.loc[split[0], 'total']+= 1
        elif re.search(">XM_", line):
            org_df.loc[split[0], 'XM']+= 1
            org_df.loc[split[0], 'total']+= 1
        elif re.search(">XR_", line):
            org_df.loc[split[0], 'XR']+= 1
            org_df.loc[split[0], 'total']+= 1
        i=i+1
    genes_unique= pd.unique(gene_list)
    print(org_df.head())
    print("Total number of lines in the data: ", len(Lines))
    print("Number of transcripts found (including all accessions):", len(fasta_details_list))
    print("Unique genes found:", len(genes_unique))
    #writing to file
    org_df["NM_XM"]=org_df["NM"]+org_df["XM"]
    org_df["NR_XR"]=org_df["NR"]+org_df["XR"]
    org_df.to_csv(outputFileName)
    gene_acc.to_csv(outputFileAcc)
    file_fasta = open(outputFileFasta, "w+")
    file_fasta.writelines(fasta_details_list)
    file_fasta.close()
    file.close()

count_transcripts("mouse_rna.fna", "05032024_mouse_transcripts_counts.csv", "05032024_mouse_transcripts_accessions.csv", "05032024_mouse_transcripts_headers.csv")
count_transcripts("zebrafish_rna.fna", "05032024_zebrafish_transcripts_counts.csv", "05032024_zebrafish_transcripts_accessions.csv", "05032024_zebrafish_transcripts_headers.csv")
count_transcripts("chimpanzee_rna.fna", "05032024_chimpanzee_transcripts_counts.csv", "05032024_chimpanzee_transcripts_accessions.csv", "05032024_chimpanzee_transcripts_headers.csv")

#18-02-2025 
#directory was changed to /home/species_dat/
count_transcripts("rat_ncbi_dataset/data/GCF_036323735.1/rna.fna", "18022025_rat_transcripts_counts.csv", "18022025_rat_transcripts_accessions.csv", "18022025_rat_transcripts_headers.csv")
count_transcripts("fruitfly_ncbi_dataset/data/GCF_000001215.4/rna.fna", "18022025_fruitfly_transcripts_counts.csv", "18022025_fruitfly_transcripts_accessions.csv", "18022025_fruitfly_transcripts_headers.csv")
count_transcripts("celegans_ncbi_dataset/data/GCF_000002985.6/rna.fna", "18022025_celegans_transcripts_counts.csv", "18022025_celegans_transcripts_accessions.csv", "18022025_celegans_transcripts_headers.csv")
count_transcripts("/home/bioschool/phd/blz227562/species_dat/chicken_ncbi_dataset/data/GCF_016699485.2/rna.fna", "18022025_chicken_transcripts_counts.csv", "18022025_chicken_transcripts_accessions.csv", "18022025_chicken_transcripts_headers.csv")
count_transcripts("/home/bioschool/phd/blz227562/species_dat/xenopus_ncbi_dataset/data/GCF_017654675.1/rna.fna", "18022025_xenopus_transcripts_counts.csv", "18022025_xenopus_transcripts_accessions.csv", "18022025_xenopus_transcripts_headers.csv")
count_transcripts("/home/bioschool/phd/blz227562/species_dat/yeast_ncbi_dataset/data/GCF_000146045.2/rna.fna", "18022025_yeast_transcripts_counts.csv", "18022025_yeast_transcripts_accessions.csv", "18022025_yeast_transcripts_headers.csv")

#22-03-2025
count_transcripts("/scratch/bioschool/phd/blz227562/bifunc/refseq/cattle_dataset/ncbi_dataset/data/GCF_002263795.3/rna.fna", "./counts/22032025_cattle_transcripts_counts.csv", "./counts/22032025_cattle_transcripts_accessions.csv", "./counts/22032025_cattle_transcripts_headers.csv")
count_transcripts("/scratch/bioschool/phd/blz227562/bifunc/refseq/arabidopsis_dataset/ncbi_dataset/data/GCF_000001735.4/rna.fna", "./counts/22032025_arabidopsis_transcripts_counts.csv", "./counts/22032025_arabidopsis_transcripts_accessions.csv", "./counts/22032025_arabidopsis_transcripts_headers.csv")
count_transcripts("/scratch/bioschool/phd/blz227562/bifunc/refseq/wheat_dataset/ncbi_dataset/data/GCF_018294505.1/rna.fna", "./counts/22032025_wheat_transcripts_counts.csv", "./counts/22032025_wheat_transcripts_accessions.csv", "./counts/22032025_wheat_transcripts_headers.csv")
count_transcripts("/scratch/bioschool/phd/blz227562/bifunc/refseq/rice_dataset/ncbi_dataset/data/GCF_034140825.1/rna.fna", "./counts/22032025_rice_transcripts_counts.csv", "./counts/22032025_rice_transcripts_accessions.csv", "./counts/22032025_rice_transcripts_headers.csv")
count_transcripts("/scratch/bioschool/phd/blz227562/bifunc/refseq/spombe_dataset/ncbi_dataset/data/GCF_000002945.2/rna.fna", "./counts/22032025_spombe_transcripts_counts.csv", "./counts/22032025_spombe_transcripts_accessions.csv", "./counts/22032025_spombe_transcripts_headers.csv")
count_transcripts("/scratch/bioschool/phd/blz227562/bifunc/refseq/rhizopus_dataset/ncbi_dataset/data/GCF_002708625.1/rna.fna", "./counts/22032025_rhizopus_transcripts_counts.csv", "./counts/22032025_rhizopus_transcripts_accessions.csv", "./counts/22032025_rhizopus_transcripts_headers.csv")
count_transcripts("/scratch/bioschool/phd/blz227562/bifunc/refseq/candida_dataset/ncbi_dataset/data/GCF_000182965.3/rna.fna", "./counts/22032025_candida_transcripts_counts.csv", "./counts/22032025_candida_transcripts_accessions.csv", "./counts/22032025_candida_transcripts_headers.csv")
count_transcripts("/scratch/bioschool/phd/blz227562/bifunc/refseq/aniger_dataset/ncbi_dataset/data/GCF_000002855.4/rna.fna", "./counts/22032025_aniger_transcripts_counts.csv", "./counts/22032025_aniger_transcripts_accessions.csv", "./counts/22032025_aniger_transcripts_headers.csv")
count_transcripts("/scratch/bioschool/phd/blz227562/bifunc/refseq/afumi_dataset/ncbi_dataset/data/GCF_000002655.1/rna.fna", "./counts/22032025_afumi_transcripts_counts.csv", "./counts/22032025_afumi_transcripts_accessions.csv", "./counts/22032025_afumi_transcripts_headers.csv")

######RUNNNNING SEPARATE SCRIPT FOR FRUITFLY!
##ON HPC, after job resources have been allocated
# conda activate nm_nr
# cd /home/
# python
import re
import pandas as pd

###MAKING A FUNCTION FOR OTHER SPECIES!
def count_transcripts(rnaFilePath, outputFileName, outputFileAcc, outputFileFasta):
    file=open(rnaFilePath, 'r')
    Lines= file.readlines()
    print("File read!")
    gene_list =[]
    fasta_details_list=[]
    #check for the '>' at the start of the sequence: indicates new sequence start
    #here we need to quantify all accession and all transcripts related to a gene
    for line in Lines:
        x=re.search("^>", line)
        if x is not None:
            split=re.split("\\(", line)
            split=re.split("\\)", split[1])
            print(split[0])
            fasta_details_list.append(line)
            gene_list.append(split[0])
    print("Gene List made!")
    
    #find unique genes
    genes_unique = pd.Series(gene_list).unique()
    print(genes_unique[0:100])
    #make a pandas data frame with gene symbols and number of transcripts for XR_, XM_, NM_, NP_ and total
    org_df= pd.DataFrame(0, columns= ["total", "NM", "NR", "XM", "XR"], index= genes_unique)
    gene_acc= pd.DataFrame(columns= ["accession", "gene"])
    
    i=0 #index start
    for line in fasta_details_list:
        print(line)
        split=re.split("\\(", line)
        split=re.split("\\)", split[1])
        gene_acc.loc[i, "accession"]=re.split(" ", line)[0]
        gene_acc.loc[i, "gene"]=split[0]
        if re.search(">NM_", line):
            org_df.loc[split[0], 'NM']+= 1
            org_df.loc[split[0], 'total']+= 1
        elif re.search(">NR_", line):
            org_df.loc[split[0], 'NR']+= 1
            org_df.loc[split[0], 'total']+= 1
        elif re.search(">XM_", line):
            org_df.loc[split[0], 'XM']+= 1
            org_df.loc[split[0], 'total']+= 1
        elif re.search(">XR_", line):
            org_df.loc[split[0], 'XR']+= 1
            org_df.loc[split[0], 'total']+= 1
        i=i+1
    genes_unique= pd.unique(gene_list)
    print(org_df.head())
    print("Total number of lines in the data: ", len(Lines))
    print("Number of transcripts found (including all accessions):", len(fasta_details_list))
    print("Unique genes found:", len(genes_unique))
    #writing to file
    org_df["NM_XM"]=org_df["NM"]+org_df["XM"]
    org_df["NR_XR"]=org_df["NR"]+org_df["XR"]
    org_df.to_csv(outputFileName)
    gene_acc.to_csv(outputFileAcc)
    file_fasta = open(outputFileFasta, "w+")
    file_fasta.writelines(fasta_details_list)
    file_fasta.close()
    file.close()
#20-04-2025
count_transcripts("fruitfly_ncbi_dataset/data/GCF_000001215.4/rna.fna", "20042025_fruitfly_transcripts_counts.csv", "20042025_fruitfly_transcripts_accessions.csv", "20042025_fruitfly_transcripts_headers.csv")


