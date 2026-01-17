import pandas as pd
import re
#changing the working directory 
#reading ORF information to match to later!
orfs_info = pd.read_csv("23122023_formatted_allpep_ATG-TTG-GTG-CTC-bifunc_genes.csv", header = 0)
#adding a new column where the sequence at kozak coordinates will be added
orfs_info.insert(len(orfs_info.columns), "kozak_seq",'')
#read NR_/XR_ sequences to get the sequences from -4 to +4 position!
nr_file = open("23122023_fasta_bifunc_genes_noncoding_transcripts.txt")
nr_lines = nr_file.readlines()
#which function, taken from an online source for ease of use here!
def which(self):
    try:
        self = list(iter(self))
    except TypeError as e:
        raise Exception("""'which' method can only be applied to iterables.
        {}""".format(str(e)))
    indices = [i for i, x in enumerate(self) if bool(x) == True]
    return(indices)
#saving bifunctional gene accessions and transcript lengths separately for other uses!
bifunc_genes_df = pd.DataFrame(columns= ["accession", "gene_symbol", "length"])
#getting the sequence at the kozak coordinates!
seq=""
row_idx=[]
for line in nr_lines:
    if line[0] == '>':
        if (len(row_idx)>0 and len(seq)>0): 
            for row in row_idx:
                coor = orfs_info.loc[row, 'coordinates']
                coord=(int(re.split("-", coor)[0][1:]))
                orfs_info.loc[row, 'kozak_seq']=seq[coord-4:coord+4]
        split=re.split(",", line)
        split_idx = [i for i in range(len(split)) if re.search(" transcript variant", split[i]) is not None]
        if len(split_idx)==1:
            split=re.split("\(", split[split_idx[0]-1])
        else:
            split=re.split("\(", split[-2])
        split=re.split("\)", split[-1])
        bifunc_genes_df.loc[len(bifunc_genes_df)]= [re.split(" ", line)[0], split[0], 0]
        print("Running for "+ split[0])
        print("Seq length"+ str(len(seq)))
        if len(bifunc_genes_df)>1:
            idx=len(bifunc_genes_df)-2
            bifunc_genes_df.loc[idx, "length"] = len(seq)
        seq="" #initializing it when a new sequence is encountered
        row_idx = which(orfs_info.acc_orf == re.split(" ", line)[0]) #getting indexes of rows which contain the accession-encoded ORFs
    else:
        seq=seq+ re.split("\n", line)[0]
#last one!
bifunc_genes_df.loc[len(bifunc_genes_df)-1, "length"] =len(seq)


#write to files
bifunc_genes_df.to_csv("26122023_bifunc-gene_noncoding-transcripts_length-acc.csv", index=False)
orfs_info.to_csv("26122023_bifunc_ORFs_allstart_kozak.csv", index=False)

