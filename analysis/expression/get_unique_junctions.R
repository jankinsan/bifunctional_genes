#Unique exon and unique exon-junction quantification for noncoding versus coding isoforms of bifunctional genes!
library(GenomicFeatures)
library(tidyverse)
library(txdbmaker)

#This was run-once on 28.04.2026, no need to re-run it!
#read gtf file and replace NCBI chromosome names with UCSC-style chromsome names to enable mapping!
setwd("D:/Janki/bifunctional genes/unique_exons")
chr_dat <- read.delim("E:/MSR BLY 10022025/THESIS WORK/NP NR Project/15122023/GCF_000001405.40_GRCh38.p14_assembly_report.txt", skip=63)
hg38 <- read.delim("C:/Users/janki/Downloads/GCF_000001405.40_GRCh38.p14_genomic.gtf", header=FALSE)
hg38[!is.na(match(hg38[,1], chr_dat$RefSeq.Accn)),1] <- chr_dat$UCSC.style.name[na.omit(match(hg38[,1], chr_dat$RefSeq.Accn))]
write.table(hg38, file="UCSC_GCF_000001405.40_GRCh38.p14_genomic.gtf", sep="\t", quote=F, row.names=F, col.names=F)

#remove all data for now!
rm(list=ls())

# load the updated gtf file
txdb <- makeTxDbFromGFF("UCSC_GCF_000001405.40_GRCh38.p14_genomic.gtf", format="gtf")

# Identify Coding vs Non-Coding Transcripts
all_tx <- transcripts(txdb, columns=c("tx_name"))
nc_tx_names <- all_tx$tx_name[c(grep("^NR_", all_tx$tx_name), grep("^XR_", all_tx$tx_name))]
coding_tx_names <- all_tx$tx_name[c(grep("^NM_", all_tx$tx_name), grep("^XM_", all_tx$tx_name))]

#read gene status!
hg38_genes <- read.csv("20062024__GRCh38.p14_genewise_status.csv.csv")
bifunctional_ids <-hg38_genes$X[hg38_genes$status=="hybrid"]
coding_ids <- hg38_genes$X[hg38_genes$status=="only_mRNA"]
noncoding_ids <- hg38_genes$X[hg38_genes$status=="only_ncRNA"]

# Categorize transcripts by the genetype!
all_tx_df <- transcripts(txdb, columns = c("tx_name", "gene_id")) %>% as.data.frame()

# Filter transcripts by coding, noncoding and bifunctional genes
bifunc_tx <- all_tx_df %>% filter(gene_id %in% bifunctional_ids)
coding_tx <- all_tx_df %>% filter (gene_id %in% coding_ids)
noncoding_tx <- all_tx_df %>% filter(gene_id %in% noncoding_ids)

# Group names by Coding (NM/XM) and Non-Coding (NR/XR)
coding_tx_bifunc <- bifunc_tx$tx_name[grepl("^[NX]M_", bifunc_tx$tx_name)]
nc_tx_bifunc <- bifunc_tx$tx_name[grepl("^[NX]R_", bifunc_tx$tx_name)]

# EXTRACT AND CATEGORIZE JUNCTIONS
all_introns <- intronsByTranscript(txdb, use.names = TRUE)

# Function to Extract unique coordinate strings with transcript names!
get_junc_strings <- function(tx_list) {
  juncs <- as.character(unlist(all_introns[names(all_introns) %in% tx_list]))
  names(juncs) <- names(as.character(unlist(all_introns[names(all_introns) %in% tx_list])))
  return(juncs)
}

#run the function for coding and noncoding transcripts
junc_coding_bifunc <- get_junc_strings(coding_tx_bifunc)
junc_nc_bifunc    <- get_junc_strings(nc_tx_bifunc)
junc_coding_other  <- get_junc_strings(coding_tx$tx_name)
junc_nc_other <- get_junc_strings(noncoding_tx$tx_name)

# make a table of coding and non-coding transcript frequencies
#how many times does each junction come up in coding/non-coding transcripts of bifunctional genes
bifunc_junc_coding_freq <- table(junc_coding_bifunc)
bifunc_junc_nc_freq <- table(junc_nc_bifunc)

# Categorize junctions so that they can be matched with expression data!
## JUNCTIONS IN BIFUNCTIONAL GENES
###1. JUNCTIONS SHARED BY CODING & NON-CODING TRANSCRIPTS
shared_juncs <- data.frame(junction = intersect(junc_coding_bifunc, junc_nc_bifunc))

#add more details such as frequencies!
shared_juncs$coding_freq <- table(junc_coding_bifunc)[match(shared_juncs$junction, names(bifunc_junc_coding_freq))]
shared_juncs$nc_freq <- table(junc_nc_bifunc)[match(shared_juncs$junction, names(bifunc_junc_nc_freq))]
shared_juncs$coding_transcript <- names(junc_coding_bifunc)[match(shared_juncs$junction, junc_coding_bifunc)]
shared_juncs$nc_transcript <- names(junc_nc_bifunc)[match(shared_juncs$junction, junc_nc_bifunc)]
shared_juncs$gene <- as.character(bifunc_tx$gene_id[match(shared_juncs$coding_transcript, bifunc_tx$tx_name)])

###2.Junctions only in coding transcripts of bifunctional genes
coding_only_juncs <- data.frame(junction=setdiff(junc_coding_bifunc, junc_nc_bifunc))

#adding details!
coding_only_juncs$coding_freq <- table(junc_coding_bifunc)[match(coding_only_juncs$junction, names(bifunc_junc_coding_freq))]
coding_only_juncs$nc_freq <- table(junc_nc_bifunc)[match(coding_only_juncs$junction, names(bifunc_junc_nc_freq))]
coding_only_juncs$coding_transcript <- names(junc_coding_bifunc)[match(coding_only_juncs$junction, junc_coding_bifunc)]
coding_only_juncs$nc_transcript <- names(junc_nc_bifunc)[match(coding_only_juncs$junction, junc_nc_bifunc)]
coding_only_juncs$gene <- as.character(bifunc_tx$gene_id[match(coding_only_juncs$coding_transcript, bifunc_tx$tx_name)])

###3. JUNCTIONS ONLY IN NON-CODING TRANSCRIPTS OF BIFUNCTIONAL GENES
nc_only_juncs <- data.frame(junction=setdiff(junc_nc_bifunc, junc_coding_bifunc))

#adding details!
nc_only_juncs$coding_freq <- table(junc_coding_bifunc)[match(nc_only_juncs$junction, names(bifunc_junc_coding_freq))]
nc_only_juncs$nc_freq <- table(junc_nc_bifunc)[match(nc_only_juncs$junction, names(bifunc_junc_nc_freq))]
nc_only_juncs$coding_transcript <- names(junc_coding_bifunc)[match(nc_only_juncs$junction, junc_coding_bifunc)]
nc_only_juncs$nc_transcript <- names(junc_nc_bifunc)[match(nc_only_juncs$junction, junc_nc_bifunc)]
nc_only_juncs$gene <- as.character(bifunc_tx$gene_id[match(nc_only_juncs$nc_transcript, bifunc_tx$tx_name)])

## JUNCTIONS IN ONLY CODING GENES
other_junc_coding_freq <- table(junc_coding_other)
other_junc_nc_freq <- table(junc_nc_other)

### 1. JUNCTIONS IN CODING GENES! 
coding_other_juncs <- data.frame(junction=unique(junc_coding_other))
coding_other_juncs$coding_freq <- table(junc_coding_other)[match(coding_other_juncs$junction, names(other_junc_coding_freq))]
coding_other_juncs$nc_freq <- table(junc_nc_other)[match(coding_other_juncs$junction, names(other_junc_nc_freq))]
coding_other_juncs$coding_transcript <- names(junc_coding_other)[match(coding_other_juncs$junction, junc_coding_other)]
coding_other_juncs$nc_transcript <- names(junc_nc_other)[match(coding_other_juncs$junction, junc_nc_other)]
coding_other_juncs$gene <- as.character(coding_tx$gene_id[match(coding_other_juncs$coding_transcript, coding_tx$tx_name)])

### 2. JUNCTIONS IN NON-CODING GENES! 
nc_other_juncs <- data.frame(junction=unique(junc_nc_other))
nc_other_juncs$coding_freq <- table(junc_coding_other)[match(nc_other_juncs$junction, names(other_junc_coding_freq))]
nc_other_juncs$nc_freq <- table(junc_nc_other)[match(nc_other_juncs$junction, names(other_junc_nc_freq))]
nc_other_juncs$coding_transcript <- names(junc_coding_other)[match(nc_other_juncs$junction, junc_coding_other)]
nc_other_juncs$nc_transcript <- names(junc_nc_other)[match(nc_other_juncs$junction, junc_nc_other)]
nc_other_juncs$gene <- as.character(noncoding_tx$gene_id[match(nc_other_juncs$nc_transcript, noncoding_tx$tx_name)])

#save files!
write.table(shared_juncs, file="bifunc_genes_shared_junctions.csv", sep=",")
write.table(coding_only_juncs, file="bifunc_genes_coding_junctions.csv", sep=",")
write.table(nc_only_juncs, file="bifunc_genes_noncoding_junctions.csv", sep=",")
write.table(coding_other_juncs, file="coding_genes_junctions.csv", sep=",")
write.table(nc_other_juncs, file="noncoding_genes_junctions.csv", sep=",")
